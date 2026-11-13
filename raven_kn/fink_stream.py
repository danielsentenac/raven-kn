"""Kilonova candidates from the Fink broker's Kafka livestream (ZTF or Rubin).

Each alert of the subscribed Fink topics becomes (or updates) a
:class:`~raven_kn.connection.KNCandidate` in the listener's store, through
:func:`raven_kn.listener.add_candidate`: it is then matched against stored GRB
notices, neutrino alerts and GW alerts, exactly like a candidate added by
hand. Optionally, :func:`raven_kn.search.assess` is re-run on every update
and the decision written to ``<store>/assess/``; a change to "trigger" is
logged, so the kilonova gate (light curve, spectrum or classifier score)
decides when a targeted search becomes worth it.

Credentials come from ``fink-client`` (``finkctl auth register -survey ...``,
stored in ``~/.finkclient/<survey>_credentials.yml``) after Fink sends them
(registration form linked from https://github.com/astrolabsoftware/fink-client).

Alert contents:

* ZTF: ``objectId``; ``candidate`` and ``prv_candidates`` with ``jd``,
  ``fid`` (1, 2, 3 = g, r, i), ``magpsf``/``sigmapsf`` for detections,
  ``diffmaglim`` for non-detections (``magpsf`` null), ``ra``/``dec``; a
  Fink kilonova score, if present in the alert, is used as ``kn_score``.
* Rubin: ``diaObject.diaObjectId``; ``diaSource``, ``prvDiaSources`` and
  ``prvDiaForcedSources`` with ``midpointMjdTai``, ``band``,
  ``psfFlux``/``psfFluxErr`` (nJy) and ``ra``/``dec``, converted as in
  :class:`~raven_kn.connection.FinkLSSTConnector` (S/N >= 5 detection, else
  a 5-sigma limit; TAI -> UTC).

Both parsers were checked on live alerts on 2026-11-13 (ZTF:
``fink_kn_candidates_ztf``, with ``rf_kn_vs_nonkn``; Rubin:
``fink_extragalactic_new_candidate_lsst``). Rubin alerts carry no kilonova
score (``clf`` holds supernova classifier scores only).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

import numpy as np

from .connection import FinkLSSTConnector, KNCandidate, PhotometryPoint

log = logging.getLogger("raven_kn.fink_stream")

# Fink's ZTF kilonova topics (doc.ztf.fink-broker.org, "Available topics", 2026-11):
# machine learning, crossmatch + property cuts, and rate-based (Andreoni et al. 2021).
ZTF_KN_TOPICS = ("fink_kn_candidates_ztf", "fink_early_kn_candidates_ztf",
                 "fink_rate_based_kn_candidates_ztf")
# Fink-added kilonova classifier scores that may be present in a ZTF alert
KN_SCORE_FIELDS = ("rf_kn_vs_nonkn",)
ZTF_BANDS = {1: "g", 2: "r", 3: "i"}


def _ztf_points(rows: Iterable[dict]) -> list[PhotometryPoint]:
    pts = []
    for r in rows or []:
        if r is None or r.get("jd") is None:
            continue
        mjd = float(r["jd"]) - 2400000.5
        band = ZTF_BANDS.get(int(r.get("fid", 0)), str(r.get("fid")))
        if r.get("magpsf") is not None and np.isfinite(r["magpsf"]):
            pts.append(PhotometryPoint(mjd, band, float(r["magpsf"]),
                                       float(r.get("sigmapsf") or 0.1), survey="ZTF"))
        elif r.get("diffmaglim") is not None:
            pts.append(PhotometryPoint(mjd, band, float(r["diffmaglim"]), upper_limit=True,
                                       survey="ZTF"))
    return pts


def ztf_alert_to_candidate(alert: dict) -> KNCandidate:
    c = alert["candidate"]
    kn = KNCandidate(name=str(alert["objectId"]), ra=float(c["ra"]), dec=float(c["dec"]),
                     pos_err_arcsec=0.5, source="Fink ZTF livestream")
    kn.add_photometry(_ztf_points([c] + list(alert.get("prv_candidates") or [])))
    for k in KN_SCORE_FIELDS:
        v = alert.get(k)
        if v is not None and np.isfinite(v):
            kn.kn_score = float(v)
            break
    return kn


def lsst_alert_to_candidate(alert: dict, snr_det: float = 5.0) -> KNCandidate:
    conv = FinkLSSTConnector(snr_det=snr_det)
    obj = alert.get("diaObject") or {}
    src = alert.get("diaSource") or {}
    oid = obj.get("diaObjectId", src.get("diaObjectId"))
    ra = obj.get("ra", src.get("ra"))
    dec = obj.get("dec", src.get("dec"))
    kn = KNCandidate(name=f"LSST{oid}", ra=float(ra), dec=float(dec), pos_err_arcsec=0.2,
                     source="Fink Rubin livestream")
    rows = [src] + list(alert.get("prvDiaSources") or []) + list(alert.get("prvDiaForcedSources") or [])
    pts, seen = [], set()
    for r in rows:
        if not r:
            continue
        p = conv._point({f"r:{k}": v for k, v in r.items()})
        if p is None:
            continue
        key = r.get("visit") or (round(p.mjd, 5), p.band)
        if key in seen:                          # a source and its forced measurement
            continue
        seen.add(key)
        pts.append(p)
    kn.add_photometry(pts)
    return kn


def alert_to_candidate(alert: dict, survey: str) -> KNCandidate:
    if survey == "ztf":
        return ztf_alert_to_candidate(alert)
    if survey == "lsst":
        return lsst_alert_to_candidate(alert)
    raise ValueError(f"unknown survey {survey!r}")


def merge_into_store(kn: KNCandidate, store) -> KNCandidate:
    """Merge a new alert's candidate into the stored one of the same name
    (new photometry added; the stored classification, redshift, GRB,
    neutrino and host information kept)."""
    for old in store.candidates():
        if old.name == kn.name:
            old.add_photometry(kn.photometry)
            if kn.kn_score is not None:
                old.kn_score = kn.kn_score
            return old
    return kn


def handle_alert(alert: dict, survey: str, store, assess_kw: Optional[dict] = None,
                 find_host: bool = False, **add_kw) -> dict:
    """One Fink alert: candidate created or updated, matched against stored
    GRB, neutrino and GW alerts, and (with ``assess_kw`` not None) assessed.
    Returns a summary (name, coincidences, decision, trigger change)."""
    from .listener import add_candidate
    new = not any(c.name == alert_to_candidate(alert, survey).name for c in store.candidates())
    kn = merge_into_store(alert_to_candidate(alert, survey), store)
    # the kilonova gate needs a distance (absolute magnitude): look the host up
    # once per candidate (recorded in its notes, also when nothing is found)
    if find_host and kn.redshift is None and "host lookup:" not in kn.notes:
        from .connection import find_host as _find_host
        try:
            res = _find_host(kn)
            msg = "found" if kn.redshift is not None else "no host redshift"
        except Exception as e:                   # network or catalogue failure: keep going
            log.warning("host lookup failed for %s: %s", kn.name, e)
            msg = f"failed ({type(e).__name__})"
        kn.notes = (kn.notes + f" host lookup: {msg}.").strip()
    coincs = add_candidate(kn, store, **add_kw)
    out = {"candidate": kn.name, "new": new, "n_photometry": len(kn.photometry),
           "coincidences": coincs}
    if assess_kw is not None:
        import json
        from .search import assess
        adir = Path(store.root) / "assess"
        prev_path = adir / f"{kn.name}_decision.json"
        prev = json.loads(prev_path.read_text()).get("trigger_targeted_search") \
            if prev_path.exists() else None
        d = assess(kn, out_dir=adir, **assess_kw)
        out["trigger"] = d["trigger_targeted_search"]
        out["trigger_changed"] = bool(d["trigger_targeted_search"]) != bool(prev)
        out["kilonova_evidence"] = d["kilonova_evidence"]["by"]
    return out


def stream(store, survey: str, topics: Iterable[str], credentials: Optional[dict] = None,
           assess_kw: Optional[dict] = None, find_host: bool = False,
           on_result: Optional[Callable[[dict], None]] = None, timeout_s: Optional[float] = None,
           max_alerts: Optional[int] = None, **add_kw) -> int:
    """Consume a Fink livestream forever (or ``max_alerts``). Credentials:
    ``fink_client.configuration.load_credentials(survey)`` unless given.
    Returns the number of alerts processed."""
    from fink_client.configuration import load_credentials
    from fink_client.consumer import AlertConsumer
    creds = credentials or load_credentials(survey=survey)
    config = {"username": creds["username"], "bootstrap.servers": creds["servers"],
              "group.id": creds["group_id"]}
    if creds.get("password"):
        config["password"] = creds["password"]
    timeout = timeout_s if timeout_s is not None else creds.get("maxtimeout", 10)
    n = 0
    with AlertConsumer(list(topics), config, survey=survey) as consumer:
        log.info("Fink %s livestream: %s", survey, ", ".join(topics))
        while max_alerts is None or n < max_alerts:
            try:
                topic, alert, _ = consumer.poll(timeout)
            except Exception:
                log.exception("Fink poll failed")
                continue
            if alert is None:
                continue
            n += 1
            try:
                res = handle_alert(alert, survey, store, assess_kw, find_host, **add_kw)
                res["topic"] = topic
                if res.get("trigger_changed") and res.get("trigger"):
                    log.info("search now worth it for %s (%s)", res["candidate"],
                             ", ".join(res.get("kilonova_evidence") or []))
                if on_result:
                    on_result(res)
            except Exception:                    # one bad alert must not stop the stream
                log.exception("failed on a %s alert from %s", survey, topic)
    return n
