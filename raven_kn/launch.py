"""Targeted searches asked for from the logbook application.

The service does not run anything: on a request carrying the Virgo
targeted-search token it writes a **launch request**
(``<store>/launch/<id>.json``) and an entry in the journal
(:func:`request_search`). The **launcher** (``raven-kn launcher``, one
long-lived process on the Virgo site) takes the requests in turn
(:func:`process`):

1. writes the search request of the optical transient (window, position,
   priors, chirp-mass range), also when its decision is *held*: a person
   asked for the search;
2. looks at the data of the site (:class:`raven_kn.targeted.datasource.SiteSource`):
   a detector is used where its frames cover what the search reads (the
   window, the longest template before it, the padding) and its state vector
   says observing. With fewer than two detectors the request is **refused**,
   each detector with its reason;
3. runs ``raven-kn targeted`` in ``<store>/runs/<optical>/<id>/``, where the
   service reads the state, the journal and the result of the run.

``data``: ``site``, the online strain of the detectors; or ``replay``, the
O4 data replayed on the site, to test the chain. The replay has nothing to do
with the sky of the alert: the window keeps its length and is moved to the
latest stretch of replay data analysable in the most detectors, and the run
says so in its name and its note.

Status of a request: ``requested``, ``running``, ``done``, ``failed``,
``refused``.
"""

from __future__ import annotations

import json
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable, Optional

from .journal import Heartbeat, Journal, RunLog, _safe, _write_atomic, iso, journal_of, utcnow

log = logging.getLogger("raven_kn.launch")

DATA = ("site", "replay")
STATISTICS = ("coherent", "network")
ACTIVE = ("requested", "running")
F_LOW = 30.0
TAIL_S = 32.0                    # read after the window: 4 x the padding of the search
LEAD_FALLBACK_S = 110.0          # before it, when the template duration cannot be computed


def _dir(root: str | Path) -> Path:
    return Path(root) / "launch"


def launches(root: str | Path, optical: Optional[str] = None) -> list[dict]:
    """Launch requests of a store, oldest first (those of one optical transient)."""
    out = []
    for p in _dir(root).glob("*.json"):
        try:
            r = json.loads(p.read_text())
        except (OSError, ValueError):
            continue
        if optical is None or r.get("optical") == optical:
            out.append(r)
    return sorted(out, key=lambda r: r.get("requested") or "")


def _save(root: str | Path, rec: dict) -> dict:
    d = _dir(root)
    d.mkdir(parents=True, exist_ok=True)
    _write_atomic(d / f"{rec['id']}.json", json.dumps(rec, indent=1))
    return rec


def request_search(store: Any, optical: str, data: str = "site", statistic: str = "coherent",
                   by: str = "application", journal: Any = None) -> dict:
    """Record that a targeted search is asked for. Raises ValueError when
    the options are not among :data:`DATA` and :data:`STATISTICS`, when the
    optical transient is not in the store, or when a search of that
    transient is already waiting or running."""
    if data not in DATA or statistic not in STATISTICS:
        raise ValueError(f"data must be one of {DATA} and statistic one of {STATISTICS}")
    root = Path(store.root)
    from .listener import optical_files
    if not any(p.stem == optical for p in optical_files(root)):
        raise ValueError(f"no optical transient {optical} in the store")
    if any(r.get("status") in ACTIVE for r in launches(root, optical)):
        raise ValueError(f"a search of {optical} is already waiting or running")
    now = utcnow()
    rec = {"id": f"{_safe(optical)}_{data}_{statistic}_{now:%Y%m%dT%H%M%S}", "optical": optical,
           "data": data, "statistic": statistic, "requested": iso(now), "requested_by": by,
           "status": "requested", "detail": "", "test": data == "replay"}
    _save(root, rec)
    (journal or journal_of(store)).log("search_requested", f"{optical}: targeted search asked for from the {by} "
                          f"({'replay data, a test' if data == 'replay' else 'site data'}, "
                          f"{statistic} statistic)", optical=optical, launch=rec["id"], data=data,
                          statistic=statistic, requested_by=by)
    return rec


def _lead(mchirp_min: float) -> float:
    try:
        from .targeted.search import template_duration
        return float(template_duration(mchirp_min, F_LOW)) + TAIL_S
    except Exception:
        return LEAD_FALLBACK_S


def _run(cmd: list[str], log_path: Path, tick: Optional[Callable[[], None]] = None) -> int:
    with open(log_path, "ab") as f:
        proc = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        while proc.poll() is None:
            if tick:
                tick()
            time.sleep(2)
        return int(proc.returncode)


def process(store: Any, rec: dict, source: Any = None, now_gps: Optional[float] = None,
            runner: Callable[..., int] = _run, tick: Optional[Callable[[], None]] = None) -> dict:
    """Carry out one launch request (steps in the module documentation).
    ``source``, ``now_gps`` and ``runner`` are replaced in tests. Returns the
    request with its final status."""
    root = Path(store.root)
    jr = journal_of(store)
    name, label = rec["optical"], rec["id"]

    def end(status: str, detail: str, level: str = "info", **data) -> dict:
        rec.update(status=status, detail=detail, ended=iso(utcnow()))
        _save(root, rec)
        jr.log(f"search_{status}", f"{name}: targeted search {label} {status}"
               + (f": {detail}" if detail else ""), level=level, optical=name, launch=label, **data)
        return rec

    try:
        kn = next((t for t in store.transients() if t.name == name), None)
        if kn is None:
            return end("refused", "the optical transient is no longer in the store", "warning")
        run_dir = root / "runs" / _safe(name) / label
        run_dir.mkdir(parents=True, exist_ok=True)
        from .search import assess
        assess(kn, out_dir=run_dir, write_request=True)
        stem = name.replace(" ", "_")
        (run_dir / f"{stem}_skymap.fits").unlink(missing_ok=True)       # 100 MB, not read by the search
        req_path = run_dir / f"{stem}_targeted_search_request.json"
        req = json.loads(req_path.read_text())
        if source is None:
            from .targeted.datasource import SiteSource
            source = SiteSource.named("live" if rec["data"] == "site" else "replay")
        lead = _lead(float((req.get("mchirp_range") or [1.0])[0]))
        length = float(req["gps_end"]) - float(req["gps_start"])
        note = f"asked for from the {rec.get('requested_by', 'application')}"
        if rec["data"] == "replay":
            if now_gps is None:
                from astropy.time import Time
                now_gps = float(Time.now().gps)
            found = source.latest_window(length, lead, TAIL_S, now_gps)
            if found is None:
                return end("refused", "no replay data analysable in two detectors in the last 6 hours",
                           "warning")
            start, ifos = found
            note += (f"; TEST on replayed O4 data: the window of {length:.0f} s is moved from GPS "
                     f"{req['gps_start']:.0f} to {start:.0f}, the latest replay data")
            req.update(gps_start=start, gps_end=start + length, original_gps_start=req["gps_start"],
                       note=note)
        else:
            cov = source.coverage(source.ifos(), float(req["gps_start"]) - lead,
                                  float(req["gps_end"]) + TAIL_S)
            ifos = [i for i, c in cov.items() if c["ok"]]
            if len(ifos) < 2:
                why = "; ".join(f"{i}: {c['why']}" for i, c in cov.items() if not c["ok"])
                return end("refused", f"fewer than two detectors with analysable data ({why})",
                           "warning", coverage={i: c["why"] or "ok" for i, c in cov.items()})
            note += "; detectors without analysable data: " + (
                ", ".join(f"{i} ({c['why']})" for i, c in cov.items() if not c["ok"]) or "none")
        req_path.write_text(json.dumps(req, indent=2))
        RunLog(run_dir, name, label, visibility="lvk").submitted(note)
        rec.update(status="running", detail=note, started=iso(utcnow()), ifos=ifos, run=label,
                   gps_start=req["gps_start"], gps_end=req["gps_end"])
        _save(root, rec)
        cmd = [sys.executable, "-m", "raven_kn.cli", "targeted", "--request", str(req_path),
               "--source", rec["data"], "--ifos", ",".join(ifos), "--statistic", rec["statistic"],
               "--label", label, "--out", str(run_dir), "--visibility", "lvk"]
        code = runner(cmd, run_dir / f"{label}.log", tick)
        if code == 0:
            return end("done", note, ifos=ifos)
        return end("failed", f"the search stopped with exit code {code} (see {label}.log in its "
                   "run directory)", "error")
    except Exception as e:                           # the launcher goes on with the next request
        log.exception("launch %s failed", label)
        return end("failed", f"{type(e).__name__}: {e}"[:300], "error")


def launcher(store: Any, interval_s: float = 10.0, once: bool = False) -> None:
    """Carry out the launch requests of a store, one at a time, forever (or
    until none is waiting, with ``once``). Requests left ``running`` by a
    launcher that stopped are marked failed at the start."""
    root = Path(store.root)
    hb = Heartbeat(root, "launcher", Journal.for_store(root, "launcher"), ["launch requests"])
    for r in launches(root):
        if r.get("status") == "running":
            r.update(status="failed", detail="the launcher stopped while this search was running",
                     ended=iso(utcnow()))
            _save(root, r)
            journal_of(store).log("search_failed", f"{r['optical']}: targeted search {r['id']} failed: "
                                  "the launcher stopped while it was running", level="error",
                                  optical=r["optical"], launch=r["id"])
    while True:
        hb.tick()
        waiting = [r for r in launches(root) if r.get("status") == "requested"]
        if waiting:
            hb.message("launch requests")
            done = process(store, waiting[0], tick=hb.tick)
            log.info("launch %s: %s %s", done["id"], done["status"], done.get("detail", ""))
        elif once:
            return
        else:
            time.sleep(interval_s)
