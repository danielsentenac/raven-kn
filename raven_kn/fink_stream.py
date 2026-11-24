"""Optical transients from the Fink broker's Kafka livestream (ZTF or Rubin).

Each alert of the subscribed Fink topics becomes (or updates) a
:class:`~raven_kn.connection.OpticalTransient` in the listener's store, through
:func:`raven_kn.listener.add_optical`: it is then matched against stored GRB
notices, neutrino alerts and GW alerts, exactly like an optical transient added by
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

from .connection import FinkLSSTConnector, OpticalTransient, PhotometryPoint

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


def ztf_alert_to_optical(alert: dict) -> OpticalTransient:
    c = alert["candidate"]
    kn = OpticalTransient(name=str(alert["objectId"]), ra=float(c["ra"]), dec=float(c["dec"]),
                     pos_err_arcsec=0.5, source="Fink ZTF livestream")
    kn.add_photometry(_ztf_points([c] + list(alert.get("prv_candidates") or [])))
    for k in KN_SCORE_FIELDS:
        v = alert.get(k)
        if v is not None and np.isfinite(v):
            kn.kn_score = float(v)
            break
    return kn


def lsst_alert_to_optical(alert: dict, snr_det: float = 5.0) -> OpticalTransient:
    conv = FinkLSSTConnector(snr_det=snr_det)
    obj = alert.get("diaObject") or {}
    src = alert.get("diaSource") or {}
    oid = obj.get("diaObjectId", src.get("diaObjectId"))
    ra = obj.get("ra", src.get("ra"))
    dec = obj.get("dec", src.get("dec"))
    kn = OpticalTransient(name=f"LSST{oid}", ra=float(ra), dec=float(dec), pos_err_arcsec=0.2,
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


def alert_to_optical(alert: dict, survey: str) -> OpticalTransient:
    if survey == "ztf":
        return ztf_alert_to_optical(alert)
    if survey == "lsst":
        return lsst_alert_to_optical(alert)
    raise ValueError(f"unknown survey {survey!r}")


def utc_now_iso() -> str:
    from .journal import iso, utcnow
    return iso(utcnow())


def merge_into_store(kn: OpticalTransient, store) -> OpticalTransient:
    """Merge a new alert's optical transient into the stored one of the same name
    (new photometry added; the stored classification, redshift, GRB,
    neutrino and host information kept)."""
    for old in store.transients():
        if old.name == kn.name:
            old.add_photometry(kn.photometry)
            if kn.kn_score is not None:
                old.kn_score = kn.kn_score
            return old
    return kn


def handle_alert(alert: dict, survey: str, store, assess_kw: Optional[dict] = None,
                 find_host: bool = False, **add_kw) -> dict:
    """One Fink alert: optical transient created or updated, matched against stored
    GRB, neutrino and GW alerts, and (with ``assess_kw`` not None) assessed.
    Returns a summary (name, coincidences, decision, trigger change)."""
    from .listener import add_optical
    new = not any(c.name == alert_to_optical(alert, survey).name for c in store.transients())
    kn = merge_into_store(alert_to_optical(alert, survey), store)
    # the kilonova gate needs a distance (absolute magnitude): look the host up
    # once per optical transient (recorded in its notes, also when nothing is found)
    if find_host and kn.redshift is None and "host lookup:" not in kn.notes:
        from .connection import find_host as _find_host
        from .journal import journal_of
        try:
            res = _find_host(kn)
            msg = "found" if kn.redshift is not None else "no host redshift"
            journal_of(store).log("host_lookup", f"{kn.name}: host lookup, {msg}", optical=kn.name,
                                  host=(res or {}).get("host") if isinstance(res, dict) else None,
                                  redshift=kn.redshift.__dict__ if kn.redshift else None)
        except Exception as e:                   # network or catalogue failure: keep going
            log.warning("host lookup failed for %s: %s", kn.name, e)
            msg = f"failed ({type(e).__name__})"
            journal_of(store).log("lookup_failed", f"{kn.name}: host lookup failed "
                                  f"({type(e).__name__}: {e})"[:300], level="warning",
                                  optical=kn.name, lookup="host", error=repr(e)[:500])
        kn.notes = (kn.notes + f" host lookup: {msg}.").strip()
    coincs = add_optical(kn, store, **add_kw)
    out = {"optical": kn.name, "new": new, "n_photometry": len(kn.photometry),
           "coincidences": coincs}
    if assess_kw is not None:
        import json
        from .search import assess
        adir = Path(store.root) / "assess"
        prev_path = adir / f"{kn.name}_decision.json"
        from .journal import journal_of, log_decision
        prev_d = json.loads(prev_path.read_text()) if prev_path.exists() else None
        prev = prev_d.get("trigger_targeted_search") if prev_d else None
        d = assess(kn, out_dir=adir, **assess_kw)
        log_decision(journal_of(store), d, prev_d)
        from .notify import on_decision
        on_decision(store, d)                    # ready for a search: announced, once
        out["trigger"] = d["trigger_targeted_search"]
        out["trigger_changed"] = bool(d["trigger_targeted_search"]) != bool(prev)
        out["kilonova_evidence"] = d["kilonova_evidence"]["by"]
    return out


def stream(store, survey: str, topics: Iterable[str], credentials: Optional[dict] = None,
           assess_kw: Optional[dict] = None, find_host: bool = False,
           on_result: Optional[Callable[[dict], None]] = None, timeout_s: Optional[float] = None,
           max_alerts: Optional[int] = None, max_new_per_night: Optional[int] = None,
           count_only: bool = False, **add_kw) -> int:
    """Consume a Fink livestream forever (or ``max_alerts``). Credentials:
    ``fink_client.configuration.load_credentials(survey)`` unless given.
    Returns the number of alerts processed.

    ``max_new_per_night``: once that many optical transients were created in a
    UTC day (a night at Palomar or in Chile falls in one), the alerts of
    objects not yet in the store are left out until the next day, with one
    warning in the journal; the transients already stored are still updated.
    A selection that floods on a bad night cannot then fill the store.

    ``count_only``: nothing is written in the store; each alert is only
    reported (``on_result``: topic, object, number of points, whether the
    object was seen before in this run), to measure what a topic sends
    before putting it in service. Fink does not send those alerts again."""
    from fink_client.configuration import load_credentials
    from fink_client.consumer import AlertConsumer
    creds = credentials or load_credentials(survey=survey)
    config = {"username": creds["username"], "bootstrap.servers": creds["servers"],
              # fink-client 12 writes "groupid", earlier versions "group_id"
              "group.id": creds.get("groupid") or creds["group_id"]}
    if creds.get("password"):
        config["password"] = creds["password"]
    timeout = timeout_s if timeout_s is not None else creds.get("maxtimeout", 10)
    n = 0
    hb = None
    night, new_tonight, left_out, seen = None, 0, 0, set()
    if not count_only and getattr(store, "journal", None) is not None:
        from .journal import Heartbeat
        hb = Heartbeat(store.root, store.journal.source, store.journal, topics)
    with AlertConsumer(list(topics), config, survey=survey) as consumer:
        log.info("Fink %s livestream: %s", survey, ", ".join(topics))
        while max_alerts is None or n < max_alerts:
            if hb:
                hb.tick()
            try:
                topic, alert, _ = consumer.poll(timeout)
            except Exception as e:
                log.exception("Fink poll failed")
                if hb:
                    hb.kafka_error(f"poll failed: {e!r}")
                continue
            if alert is None:
                continue
            n += 1
            if hb:
                hb.message(topic)
            try:
                if count_only:
                    kn = alert_to_optical(alert, survey)
                    if on_result:
                        on_result({"optical": kn.name, "topic": topic, "first_in_run": kn.name not in seen,
                                   "n_photometry": len(kn.photometry), "kn_score": kn.kn_score,
                                   "received": utc_now_iso()})
                    seen.add(kn.name)
                    continue
                if max_new_per_night is not None:
                    today = utc_now_iso()[:10]
                    if today != night:
                        night, new_tonight, left_out = today, 0, 0
                    name = alert_to_optical(alert, survey).name
                    if not any(c.name == name for c in store.transients()):
                        if new_tonight >= max_new_per_night:
                            left_out += 1
                            if left_out == 1:
                                from .journal import journal_of
                                journal_of(store).log(
                                    "stream_capped", f"Fink {survey}: {max_new_per_night} new optical transients "
                                    f"on {today} UTC, the limit: the new objects of {', '.join(topics)} are left "
                                    f"out until tomorrow", level="warning", survey=survey, limit=max_new_per_night)
                                log.warning("limit of %d new transients reached for %s", max_new_per_night, today)
                            continue
                        new_tonight += 1
                res = handle_alert(alert, survey, store, assess_kw, find_host, **add_kw)
                res["topic"] = topic
                if res.get("trigger_changed") and res.get("trigger"):
                    log.info("search now worth it for %s (%s)", res["optical"],
                             ", ".join(res.get("kilonova_evidence") or []))
                if on_result:
                    on_result(res)
            except Exception as e:               # one bad alert must not stop the stream
                log.exception("failed on a %s alert from %s", survey, topic)
                if hb:
                    hb.failure(topic, e)
    return n
