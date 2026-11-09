"""GCN Kafka listener: the real-time entry point of RAVEN-KN.

RAVEN listens to GraceDB/GCN and, for each new GW or external event, looks
for a counterpart in a fixed window. RAVEN-KN does the same with kilonova
candidates, whose merger-time windows are hours to days long:

* ``igwn.gwalert`` (LVK JSON alerts): for each kilonova candidate in the
  store whose t0 window contains the GW merger time, compute the joint FAR
  (temporal x sky x distance, :func:`raven_kn.search.coinc_far`) from the
  multi-order sky map carried in the alert. When LVK's own RAVEN found a
  GW + GRB coincidence, the alert's ``external_coinc`` section (observatory,
  time difference, temporal and sky-position joint FARs, combined sky map)
  is kept with the superevent and copied into every RAVEN-KN result; the
  combined map is used for the sky factor only on request (``use_combined``);
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
              "gcn.notices.svom.eclairs.alert",
              # classic VOEvent (XML): Fermi-GBM flight / ground / final positions
              "gcn.classic.voevent.FERMI_GBM_FLT_POS", "gcn.classic.voevent.FERMI_GBM_GND_POS",
              "gcn.classic.voevent.FERMI_GBM_FIN_POS")
GBM_SYS_ERR_DEG = 3.7       # GBM systematic localisation error (Connaughton et al. 2015, core)
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
        out["skymap"] = _decode_moc(ev["skymap"])
    out["external_coinc"], out["combined_skymap"] = parse_external_coinc(msg.get("external_coinc"))
    return out


def _decode_moc(sm: bytes | str) -> Table:
    return read_moc_bytes(base64.b64decode(sm) if isinstance(sm, str) else sm)


# fields of the alert's ``external_coinc`` section (LVK alert schema)
EXTERNAL_COINC_KEYS = ("gcn_notice_id", "ivorn", "observatory", "search", "time_difference",
                       "time_coincidence_far", "time_sky_position_coincidence_far")


def parse_external_coinc(ec: Optional[dict]) -> tuple[Optional[dict], Optional[Table]]:
    """RAVEN coincidence block of an LVK alert -> (metadata, combined MOC map).

    The metadata keeps every field except the map and adds ``raven_joint_far``,
    the sky-position FAR if given, else the temporal one."""
    if not ec:
        return None, None
    meta = {k: v for k, v in ec.items() if k != "combined_skymap"}
    for k in ("time_difference", "time_coincidence_far", "time_sky_position_coincidence_far"):
        if meta.get(k) is not None:
            meta[k] = float(meta[k])
    meta["raven_joint_far"] = (meta.get("time_sky_position_coincidence_far")
                               if meta.get("time_sky_position_coincidence_far") is not None
                               else meta.get("time_coincidence_far"))
    combined = _decode_moc(ec["combined_skymap"]) if ec.get("combined_skymap") else None
    meta["has_combined_skymap"] = combined is not None
    return meta, combined


def parse_voevent(xml: bytes | str) -> dict[str, Any]:
    """Classic GCN VOEvent (Fermi-GBM, Swift-BAT...): ISOTime, C1/C2 and
    Error2Radius (degrees) from WhereWhen; TrigID from the parameters. For
    GBM the systematic error is added in quadrature."""
    import xml.etree.ElementTree as ET
    root = ET.fromstring(xml.encode() if isinstance(xml, str) else xml)

    def find(tag):
        return next((e for e in root.iter() if e.tag.split("}")[-1] == tag), None)
    params = {e.get("name"): e.get("value") for e in root.iter()
              if e.tag.split("}")[-1] == "Param"}
    ivorn = root.get("ivorn", "")
    t = find("ISOTime").text.strip()
    ra, dec, err = (find(k) for k in ("C1", "C2", "Error2Radius"))
    instrument = "Fermi-GBM" if "Fermi" in ivorn or "GBM" in ivorn else ivorn.split("/")[-1].split("#")[0]
    e = float(err.text) if err is not None else None
    if e is not None and instrument == "Fermi-GBM":
        e = float(np.hypot(e, GBM_SYS_ERR_DEG))
    return {"grb": GRBAssociation(gps=_gps(t), instrument=instrument,
                                  graceid=params.get("TrigID")),
            "ra": float(ra.text) if ra is not None else None,
            "dec": float(dec.text) if dec is not None else None,
            "err_deg": e, "role": root.get("role")}


def parse_grb_notice(msg: dict | bytes | str) -> dict[str, Any]:
    """GCN notice, unified JSON schema (time, position and error in
    degrees) or classic VOEvent XML (:func:`parse_voevent`)."""
    if isinstance(msg, (bytes, str)):
        text = msg.decode() if isinstance(msg, bytes) else msg
        if text.lstrip().startswith("<"):
            return parse_voevent(text)
        msg = json.loads(text)
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

    def gw_skymap(self, superevent_id: str, combined: bool = False) -> Optional[Table]:
        """GW-only map, or with ``combined`` the RAVEN GW x external map."""
        kind = "combined.multiorder" if combined else "multiorder"
        p = self.root / "gw" / f"{superevent_id}.{kind}.fits"
        return Table.read(p) if p.exists() else None

    def save_gw(self, alert: dict) -> None:
        """Latest alert per superevent. A RAVEN coincidence (``external_coinc``)
        seen in an earlier alert is kept when a later one omits it."""
        sid = alert["superevent_id"]
        jpath = self.root / "gw" / f"{sid}.json"
        prev = json.loads(jpath.read_text()) if jpath.exists() else {}
        meta = {k: v for k, v in alert.items() if k not in ("skymap", "combined_skymap")}
        if alert.get("skymap") is not None:
            alert["skymap"].write(self.root / "gw" / f"{sid}.multiorder.fits", overwrite=True)
        elif prev:
            meta.setdefault("skymap_from", "previous alert")
        if alert.get("combined_skymap") is not None:
            alert["combined_skymap"].write(self.root / "gw" / f"{sid}.combined.multiorder.fits",
                                           overwrite=True)
        if not meta.get("external_coinc") and prev.get("external_coinc"):
            meta["external_coinc"] = prev["external_coinc"]
            meta["external_coinc_from"] = "previous alert"
        jpath.write_text(json.dumps(meta, indent=1, default=str))

    def forget_gw(self, superevent_id: str) -> None:
        for p in (self.root / "gw").glob(f"{superevent_id}.*"):
            p.unlink()

    def save_coinc(self, res: dict) -> Path:
        p = self.root / "coinc" / f"{res['superevent_id']}__{res['external_event']}.json"
        p.write_text(json.dumps(res, indent=1, default=float))
        return p


def _evaluate(kn: KNCandidate, gw: dict, moc: Optional[Table], combined: Optional[Table] = None,
              use_combined: bool = False, **coinc_kw) -> Optional[dict]:
    """Joint FAR if the GW time is inside the kilonova window, else None.

    ``combined`` is RAVEN's GW x external-trigger map. It is looked up at the
    kilonova position and reported; it is used for the FAR only with
    ``use_combined`` (it already folds in the GRB localisation, so the
    sky factor is then no longer GW-only)."""
    if gw.get("gpstime") is None or gw.get("far") is None:
        return None
    # same window as raven_kn.search.assess: the 95th-percentile EM distance
    # decides which non-detections constrain the merger time
    post = posterior_for(kn) if (kn.redshift or kn.distance or kn.lc_fit_samples) else None
    d95 = float(np.quantile(post.d_mpc, 0.95)) if post is not None and post.constraints else None
    window = estimate_window(kn, d_mpc_upper=d95)
    if not window.gps_start <= gw["gpstime"] <= window.gps_end:
        return None
    far_map = combined if (use_combined and combined is not None) else moc
    flat = moc_to_flat(far_map) if far_map is not None else None
    res = coinc_far(kn, {"superevent_id": gw["superevent_id"], "far": gw["far"]}, window,
                    gw_skymap=flat, post=post, **coinc_kw)
    res.update(alert_type=gw.get("alert_type"), gw_gpstime=gw["gpstime"],
               window=asdict(window), dt_from_first_detection_h=None)
    if kn.first_detection is not None:
        res["dt_from_first_detection_h"] = (kn.first_detection.gps - gw["gpstime"]) / 3600
    for key, m in (("gw_map_at_kn", moc), ("combined_map_at_kn", combined)):
        if m is not None:
            try:
                res[key] = moc_lookup(m, kn.ra, kn.dec)
            except ValueError:
                pass
    res["skymap_used"] = "combined" if far_map is combined and combined is not None else "gw"
    res["raven_external_coinc"] = gw.get("external_coinc")
    return res


def handle_gw_alert(msg: dict | bytes | str, store: CandidateStore, use_combined: bool = False,
                    **coinc_kw) -> list[dict]:
    alert = parse_gw_alert(msg)
    sid = alert["superevent_id"]
    if alert["alert_type"] == "RETRACTION":
        store.forget_gw(sid)
        log.info("%s retracted", sid)
        return []
    moc = alert["skymap"] if alert["skymap"] is not None else store.gw_skymap(sid)
    store.save_gw(alert)
    ec = alert["external_coinc"]
    if ec:
        log.info("%s: RAVEN coincidence with %s %s, dt %s s, joint FAR %s Hz", sid,
                 ec.get("observatory"), ec.get("search"), ec.get("time_difference"),
                 ec.get("raven_joint_far"))
    gw = json.loads((store.root / "gw" / f"{sid}.json").read_text())   # with carried-over coinc
    combined = store.gw_skymap(sid, combined=True)
    out = []
    for kn in store.candidates():
        res = _evaluate(kn, gw, moc, combined, use_combined, **coinc_kw)
        if res is not None:
            store.save_coinc(res)
            log.info("RAVEN-KN: %s in the window of %s, FAR %.3g Hz", sid, kn.name,
                     res.get("spatiotemporal_distance_coinc_far") or res.get("spatiotemporal_coinc_far") or res["temporal_coinc_far"])
            out.append(res)
    return out


def handle_grb_notice(msg: dict | bytes | str, store: CandidateStore,
                      max_lead_days: float = 3.0, min_err_deg: float = 0.0,
                      err_scale: float = 2.0) -> list[str]:
    """Attach a GRB to candidates within ``err_scale`` x its error radius
    (or ``min_err_deg``) whose first detection follows it by
    < ``max_lead_days``. Notice radii are ~68 % (1-sigma-like) circles;
    x2 is ~99 % for a 2D Gaussian (GW170817: NGC 4993 lies 24 deg from the
    GBM ground position, radius 11.6 deg stat. + 3.7 deg sys.). Candidates
    that already have a GRB are left alone; test notices are ignored.
    Returns the names of the candidates updated."""
    n = parse_grb_notice(msg)
    if n.get("role") == "test" or n["ra"] is None or n["dec"] is None:
        return []
    radius = max(err_scale * (n["err_deg"] or 0.0), min_err_deg)
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
                  use_combined: bool = False, **coinc_kw) -> list[dict]:
    """Store a new or updated candidate and match it against GW alerts
    received in the last ``gw_memory_days``."""
    store.save(kn)
    now = now_gps if now_gps is not None else float(Time.now().gps)
    out = []
    for gw in store.gw_alerts():
        if gw.get("gpstime") is None or now - gw["gpstime"] > store.gw_memory_s:
            continue
        sid = gw["superevent_id"]
        res = _evaluate(kn, gw, store.gw_skymap(sid), store.gw_skymap(sid, combined=True),
                        use_combined, **coinc_kw)
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
