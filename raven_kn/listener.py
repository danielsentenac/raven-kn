"""GCN Kafka listener: the real-time entry point of RAVEN-KN.

RAVEN listens to GraceDB/GCN and, for each new GW or external event, looks
for a counterpart in a fixed window. RAVEN-KN does the same with kilonova
candidates, whose merger-time windows are hours to days long:

* ``igwn.gwalert`` (LVK JSON alerts): for each kilonova candidate in the
  store whose t0 window contains the GW merger time, compute the joint FAR
  (temporal x sky x distance, :func:`raven_kn.search.coinc_far`) from the
  multi-order sky map carried in the alert;
* GRB notices in GCN's unified JSON schema (``trigger_time``, ``ra``,
  ``dec``, ``ra_dec_error``): attach the GRB to candidates at a compatible
  position whose first detection follows it, which narrows their window to
  RAVEN's [-5, +1] s;
* a kilonova candidate added later (from a broker, a circular, by hand) is
  matched against the GW alerts already received.

Handlers work on decoded JSON, so they can be tested without Kafka;
:func:`listen` wires them to ``gcn_kafka.Consumer`` (``pip install
gcn-kafka``, credentials from https://gcn.nasa.gov). Retractions remove the
superevent from the store. Results are written as JSON under the store
directory and returned; nothing is uploaded to GraceDB from here.
"""

from __future__ import annotations

import base64
import io
import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
import astropy.units as u

from .connection import GRBAssociation, KNCandidate, load_candidate, save_candidate
from .search import coinc_far
from .treatment.distance_inclination import posterior_for
from .treatment.t0 import estimate_window

log = logging.getLogger("raven_kn.listener")

GW_TOPIC = "igwn.gwalert"
# GCN unified-schema JSON GRB topics (others can be passed to listen())
GRB_TOPICS = ("gcn.notices.swift.bat.guano", "gcn.notices.einstein_probe.wxt.alert",
              "gcn.notices.svom.eclairs.alert")
MAX_FLAT_ORDER = 9          # flattened GW maps: nside 512


# --------------------------------------------------------------------------
# Multi-order (MOC) sky maps
# --------------------------------------------------------------------------

def read_moc_bytes(data: bytes) -> Table:
    return Table.read(io.BytesIO(data), format="fits")


def uniq_to_order_ipix(uniq: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    uniq = np.asarray(uniq, dtype=np.int64)
    order = (np.floor(np.log2(uniq) / 2) - 1).astype(np.int64)
    return order, uniq - 4 * (np.int64(1) << (2 * order))


def moc_lookup(moc: Table, ra: float, dec: float) -> dict[str, float]:
    """Row of the multi-order map containing (ra, dec), degrees."""
    order, ipix = uniq_to_order_ipix(moc["UNIQ"])
    for o in np.unique(order):
        pix = hp.ang2pix(1 << int(o), ra, dec, lonlat=True, nest=True)
        hit = np.flatnonzero((order == o) & (ipix == pix))
        if hit.size:
            row = moc[int(hit[0])]
            return {c: float(row[c]) for c in moc.colnames if c != "UNIQ"}
    raise ValueError("position not covered by the multi-order map")


def moc_to_flat(moc: Table, max_order: int = MAX_FLAT_ORDER) -> dict:
    """Flat NESTED map (the dict of :func:`raven_kn.skymap.read_gw_skymap`).

    Pixels finer than ``max_order`` are merged (probability summed, distance
    parameters probability-weighted), coarser ones are repeated.
    """
    order, ipix = uniq_to_order_ipix(moc["UNIQ"])
    target = int(min(order.max(), max_order))
    npix = hp.nside2npix(1 << target)
    dens = np.asarray(moc["PROBDENSITY"], dtype=float)
    area = 4 * np.pi / (12 * 4.0 ** order)
    prob_row = dens * area
    layers = [c for c in ("DISTMU", "DISTSIGMA", "DISTNORM") if c in moc.colnames]
    prob = np.zeros(npix)
    acc = {c: np.zeros(npix) for c in layers}
    wsum = np.zeros(npix)
    for o in np.unique(order):
        sel = order == o
        if o <= target:
            k = 4 ** (target - o)
            idx = (ipix[sel][:, None] * k + np.arange(k)[None, :]).ravel()
            prob[idx] = np.repeat(prob_row[sel] / k, k)
            for c in layers:
                acc[c][idx] = np.repeat(np.asarray(moc[c], float)[sel], k)
            wsum[idx] = 1.0
        else:
            parent = ipix[sel] >> (2 * (o - target))
            w = prob_row[sel]
            prob += np.bincount(parent, w, npix)
            ws = np.bincount(parent, w, npix)
            for c in layers:
                v = np.asarray(moc[c], float)[sel]
                good = np.isfinite(v)
                acc[c] += np.bincount(parent[good], (w * v)[good], npix)
            wsum += ws
    out = {"PROB": prob, "nested": True, "meta": dict(moc.meta)}
    fine = order.max() > target
    for c in layers:
        if fine:
            with np.errstate(invalid="ignore", divide="ignore"):
                merged = np.where(wsum > 0, acc[c] / np.where(wsum > 0, wsum, 1), np.inf)
            out[c] = merged
        else:
            out[c] = acc[c]
    return out


# --------------------------------------------------------------------------
# Alert parsing
# --------------------------------------------------------------------------

def _gps(iso: str) -> float:
    return float(Time(iso.replace("Z", ""), format="isot", scale="utc").gps)


def parse_gw_alert(msg: dict | bytes | str) -> dict[str, Any]:
    """Decode an ``igwn.gwalert`` JSON message (sky map kept as a Table)."""
    if isinstance(msg, (bytes, str)):
        msg = json.loads(msg)
    ev = msg.get("event") or {}
    out = {"superevent_id": msg["superevent_id"], "alert_type": msg["alert_type"],
           "time_created": msg.get("time_created"), "gpstime": None, "far": None,
           "significant": ev.get("significant"), "group": ev.get("group"),
           "pipeline": ev.get("pipeline"), "search": ev.get("search"),
           "instruments": ev.get("instruments"), "classification": ev.get("classification"),
           "properties": ev.get("properties"), "skymap": None}
    if ev.get("time"):
        out["gpstime"] = _gps(ev["time"])
    if ev.get("far") is not None:
        out["far"] = float(ev["far"])
    if ev.get("skymap"):
        sm = ev["skymap"]
        out["skymap"] = read_moc_bytes(base64.b64decode(sm) if isinstance(sm, str) else sm)
    return out


def parse_grb_notice(msg: dict | bytes | str) -> dict[str, Any]:
    """GCN unified-schema notice: time, position and its error (degrees)."""
    if isinstance(msg, (bytes, str)):
        msg = json.loads(msg)
    err = msg.get("ra_dec_error")
    return {"grb": GRBAssociation(gps=_gps(msg["trigger_time"]),
                                  instrument=msg.get("instrument", msg.get("mission", "")),
                                  graceid=msg.get("id") if isinstance(msg.get("id"), str) else None),
            "ra": msg.get("ra"), "dec": msg.get("dec"),
            "err_deg": float(err) if err is not None else None}


# --------------------------------------------------------------------------
# Store and handlers
# --------------------------------------------------------------------------

class CandidateStore:
    """Kilonova candidates, received GW alerts and coincidence results,
    kept as JSON (and FITS sky maps) under one directory."""

    def __init__(self, root: str | Path, gw_memory_days: float = 14.0):
        self.root = Path(root)
        for d in ("candidates", "gw", "coinc"):
            (self.root / d).mkdir(parents=True, exist_ok=True)
        self.gw_memory_s = gw_memory_days * 86400

    # kilonova candidates
    def candidates(self) -> list[KNCandidate]:
        return [load_candidate(p) for p in sorted((self.root / "candidates").glob("*.json"))]

    def save(self, kn: KNCandidate) -> None:
        save_candidate(kn, self.root / "candidates" / f"{kn.name}.json")

    # GW alerts (latest per superevent)
    def gw_alerts(self) -> list[dict]:
        return [json.loads(p.read_text()) for p in sorted((self.root / "gw").glob("*.json"))]

    def gw_skymap(self, superevent_id: str) -> Optional[Table]:
        p = self.root / "gw" / f"{superevent_id}.multiorder.fits"
        return Table.read(p) if p.exists() else None

    def save_gw(self, alert: dict) -> None:
        sid = alert["superevent_id"]
        meta = {k: v for k, v in alert.items() if k != "skymap"}
        if alert.get("skymap") is not None:
            alert["skymap"].write(self.root / "gw" / f"{sid}.multiorder.fits", overwrite=True)
        elif (self.root / "gw" / f"{sid}.json").exists():
            meta.setdefault("skymap_from", "previous alert")
        (self.root / "gw" / f"{sid}.json").write_text(json.dumps(meta, indent=1, default=str))

    def forget_gw(self, superevent_id: str) -> None:
        for p in (self.root / "gw").glob(f"{superevent_id}.*"):
            p.unlink()

    def save_coinc(self, res: dict) -> Path:
        p = self.root / "coinc" / f"{res['superevent_id']}__{res['external_event']}.json"
        p.write_text(json.dumps(res, indent=1, default=float))
        return p


def _evaluate(kn: KNCandidate, gw: dict, moc: Optional[Table], **coinc_kw) -> Optional[dict]:
    """Joint FAR if the GW time is inside the kilonova window, else None."""
    if gw.get("gpstime") is None or gw.get("far") is None:
        return None
    window = estimate_window(kn)
    if not window.gps_start <= gw["gpstime"] <= window.gps_end:
        return None
    flat = moc_to_flat(moc) if moc is not None else None
    post = None
    if flat is not None and (kn.redshift or kn.distance or kn.lc_fit_samples):
        post = posterior_for(kn)
    res = coinc_far(kn, {"superevent_id": gw["superevent_id"], "far": gw["far"]}, window,
                    gw_skymap=flat, post=post, **coinc_kw)
    res.update(alert_type=gw.get("alert_type"), gw_gpstime=gw["gpstime"],
               window=asdict(window), dt_from_first_detection_h=None)
    if kn.first_detection is not None:
        res["dt_from_first_detection_h"] = (kn.first_detection.gps - gw["gpstime"]) / 3600
    if moc is not None:
        try:
            res["gw_map_at_kn"] = moc_lookup(moc, kn.ra, kn.dec)
        except ValueError:
            pass
    return res


def handle_gw_alert(msg: dict | bytes | str, store: CandidateStore, **coinc_kw) -> list[dict]:
    alert = parse_gw_alert(msg)
    sid = alert["superevent_id"]
    if alert["alert_type"] == "RETRACTION":
        store.forget_gw(sid)
        log.info("%s retracted", sid)
        return []
    moc = alert["skymap"] if alert["skymap"] is not None else store.gw_skymap(sid)
    store.save_gw(alert)
    out = []
    for kn in store.candidates():
        res = _evaluate(kn, alert, moc, **coinc_kw)
        if res is not None:
            store.save_coinc(res)
            log.info("RAVEN-KN: %s in the window of %s, FAR %.3g Hz", sid, kn.name,
                     res.get("spatiotemporal_distance_coinc_far") or res.get("spatiotemporal_coinc_far") or res["temporal_coinc_far"])
            out.append(res)
    return out


def handle_grb_notice(msg: dict | bytes | str, store: CandidateStore,
                      max_lead_days: float = 3.0, min_err_deg: float = 0.0) -> list[str]:
    """Attach a GRB to candidates inside its error circle (or within
    ``min_err_deg``) whose first detection follows it by < ``max_lead_days``.
    Candidates that already have a GRB are left alone. Returns their names."""
    n = parse_grb_notice(msg)
    if n["ra"] is None or n["dec"] is None:
        return []
    radius = max(n["err_deg"] or 0.0, min_err_deg)
    grb_pos = SkyCoord(n["ra"] * u.deg, n["dec"] * u.deg)
    hits = []
    for kn in store.candidates():
        fd = kn.first_detection
        if kn.grb is not None or fd is None:
            continue
        lead = fd.gps - n["grb"].gps
        sep = grb_pos.separation(SkyCoord(kn.ra * u.deg, kn.dec * u.deg)).deg
        if 0 < lead < max_lead_days * 86400 and sep <= radius + kn.pos_err_arcsec / 3600:
            kn.grb = n["grb"]
            kn.notes = (kn.notes + f" GRB {n['grb'].graceid or ''} ({n['grb'].instrument}) "
                        f"attached by listener, sep {sep:.2f} deg.").strip()
            store.save(kn)
            hits.append(kn.name)
    return hits


def add_candidate(kn: KNCandidate, store: CandidateStore, now_gps: Optional[float] = None,
                  **coinc_kw) -> list[dict]:
    """Store a new or updated candidate and match it against GW alerts
    received in the last ``gw_memory_days``."""
    store.save(kn)
    now = now_gps if now_gps is not None else float(Time.now().gps)
    out = []
    for gw in store.gw_alerts():
        if gw.get("gpstime") is None or now - gw["gpstime"] > store.gw_memory_s:
            continue
        res = _evaluate(kn, gw, store.gw_skymap(gw["superevent_id"]), **coinc_kw)
        if res is not None:
            store.save_coinc(res)
            out.append(res)
    return out


# --------------------------------------------------------------------------
# Kafka loop
# --------------------------------------------------------------------------

def listen(store: CandidateStore, client_id: str, client_secret: str,
           grb_topics: Iterable[str] = GRB_TOPICS, domain: Optional[str] = None,
           on_result: Optional[Callable[[dict], None]] = None, **coinc_kw) -> None:
    """Consume GCN Kafka forever (``gcn_kafka.Consumer``)."""
    from gcn_kafka import Consumer
    kw = {"domain": domain} if domain else {}
    consumer = Consumer(client_id=client_id, client_secret=client_secret, **kw)
    topics = [GW_TOPIC, *grb_topics]
    consumer.subscribe(topics)
    log.info("listening to %s", ", ".join(topics))
    while True:
        for message in consumer.consume(timeout=1):
            if message.error():
                log.warning("kafka: %s", message.error())
                continue
            try:
                if message.topic() == GW_TOPIC:
                    for res in handle_gw_alert(message.value(), store, **coinc_kw):
                        if on_result:
                            on_result(res)
                else:
                    for name in handle_grb_notice(message.value(), store):
                        log.info("GRB attached to %s", name)
            except Exception:       # one bad message must not stop the listener
                log.exception("failed on a %s message", message.topic())
