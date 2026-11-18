"""Logbook of the pipeline: journal, listener heartbeats, job state, sky-map versions.

The store keeps only the latest state (latest alert per superevent, latest
decision per candidate, alerts deleted on retraction and after 14 days). The
journal keeps what happened: one JSON line per step, with the values of that
moment, so that it can be read later without the files it refers to.

* **Journal**: append-only ``<store>/journal/<source>.jsonl``, one file per
  writing process (``listen``, ``fink-ztf``, ``cli``...), so that no lock is
  needed, also on a network file system. A search or PE job writes
  ``<run>_journal.jsonl`` in its output directory, which is copied to
  ``<store>/runs/<candidate>/<directory>/``. :func:`read_events` merges all files;
  its cursor is the number of lines already read in each file, which is
  exact for append-only files whatever the order in which they arrive.
* **Heartbeat**: ``<store>/status/<source>.json``, rewritten every 30 s by a
  listener (:class:`Heartbeat`). At a start, the previous heartbeat gives the
  length of the outage, written in the journal.
* **Job state**: ``<run>_state.json`` in the output directory of a search or
  PE (:class:`RunLog`; the run is the ``--label`` of ``raven-kn targeted``): stage, sub-segment i of N, heartbeat, exit, error.
* **Sky-map versions**: ``<store>/gw_versions/<superevent>/``, one reduced
  map per alert carrying a map (:func:`save_sky_version`): the credible level
  of each pixel at nside 128, one byte per pixel, gzip. Never deleted.

Entry: ``{"time", "source", "kind", "level", "summary", "optical",
"superevent", "run", "visibility", "data"}``; ``optical`` is the name of the
optical transient the entry is about. Entries written before the rename of
the optical transient (key ``candidate``, kinds ``candidate_new`` and
``candidate_update``) are given with the new names by :func:`read_events`. ``level`` is ``info``,
``warning`` or ``error``; ``visibility`` is ``lvk`` for entries about
non-public data (hidden by ``raven-kn serve`` without the LVK token).
"""

from __future__ import annotations

import base64
import gzip
import json
import os
import socket
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

import numpy as np

SKY_ORDER = 7                    # reduced maps: nside 128, 196608 pixels
CRED_SCALE = 250                 # byte = ceil(credible level * 250); 90 % region: byte <= 225
HEARTBEAT_S = 30.0
RUN_HEARTBEAT_S = 60.0


def utcnow() -> datetime:
    return datetime.now(timezone.utc)


def iso(t: datetime) -> str:
    return t.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"


def parse_iso(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def _clean(o: Any) -> Any:
    """JSON-safe copy: numpy scalars and arrays, dataclasses, non-finite floats."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return f if np.isfinite(f) else None
    if hasattr(o, "__dataclass_fields__"):
        from dataclasses import asdict
        return _clean(asdict(o))
    if o is None or isinstance(o, (str, int)):
        return o
    return str(o)


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(name))


def _write_atomic(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, path)


class Journal:
    """Append-only journal of one writing process. ``clock`` can be replaced
    (a replay writes the time of the replayed event, not of the run)."""

    def __init__(self, path: str | Path, source: str,
                 clock: Callable[[], datetime] = utcnow, defaults: Optional[dict] = None):
        self.path = Path(path)
        self.source = source
        self.clock = clock
        self.defaults = defaults or {}

    @classmethod
    def for_store(cls, root: str | Path, source: str = "store", **kw) -> "Journal":
        return cls(Path(root) / "journal" / f"{_safe(source)}.jsonl", source, **kw)

    def log(self, kind: str, summary: str, level: str = "info", optical: Optional[str] = None,
            superevent: Optional[str] = None, **data) -> dict:
        entry = {"time": iso(self.clock()), "source": self.source, "kind": kind, "level": level,
                 "summary": summary, "optical": optical, "superevent": superevent,
                 **self.defaults, "data": _clean(data)}
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.path, "a") as f:        # one write per line
                f.write(json.dumps(entry, separators=(",", ":")) + "\n")
        except OSError:                            # the logbook must never stop the pipeline
            pass
        return entry


class NullJournal:
    source = "none"
    clock = staticmethod(utcnow)

    def log(self, *a, **kw) -> dict:
        return {}


def journal_of(store: Any) -> Any:
    return getattr(store, "journal", None) or NullJournal()


# --------------------------------------------------------------------------
# Reading: merged view with a cursor
# --------------------------------------------------------------------------

OLD_KINDS = {"candidate_new": "optical_new", "candidate_update": "optical_update"}


def journal_files(root: str | Path) -> dict[str, Path]:
    """Journal files of a store, keyed by a stable short name."""
    root = Path(root)
    out = {p.stem: p for p in sorted((root / "journal").glob("*.jsonl"))}
    for p in sorted((root / "runs").glob("*/*/*_journal.jsonl")):
        out[f"run/{p.parent.parent.name}/{p.parent.name}/{p.name[:-len('_journal.jsonl')]}"] = p
    return out


def encode_cursor(pos: dict[str, int]) -> str:
    return base64.urlsafe_b64encode(json.dumps(pos, separators=(",", ":"), sort_keys=True).encode()).decode()


def decode_cursor(cursor: Optional[str]) -> dict[str, int]:
    if not cursor:
        return {}
    try:
        d = json.loads(base64.urlsafe_b64decode(cursor.encode()))
        return {str(k): int(v) for k, v in d.items()}
    except Exception:
        return {}


def read_events(root: str | Path, cursor: Optional[str] = None, limit: int = 2000,
                allow_lvk: bool = True) -> tuple[list[dict], str, bool]:
    """Entries not yet seen by ``cursor``, oldest first, at most ``limit``.
    Returns (entries, new cursor, more). Each entry gets ``id``
    (``<file>:<line>``), stable for append-only files. A last line still
    being written (no newline) is left for the next call."""
    pos = decode_cursor(cursor)
    new: list[tuple[str, int, dict]] = []
    counts: dict[str, int] = {}
    for key, path in journal_files(root).items():
        try:
            text = path.read_text()
        except OSError:
            continue
        lines = text.split("\n")[:-1]               # drops a last line still being written
        counts[key] = len(lines)
        for n in range(pos.get(key, 0), len(lines)):
            try:
                e = json.loads(lines[n])
            except ValueError:
                e = {"time": "", "source": key, "kind": "journal_error", "level": "error",
                     "summary": "unreadable journal line", "data": {"line": lines[n][:200]}}
            e["id"] = f"{key}:{n}"
            if "candidate" in e and "optical" not in e:      # entries written before the rename
                e["optical"] = e.pop("candidate")
            e["kind"] = OLD_KINDS.get(e.get("kind"), e.get("kind"))
            new.append((key, n, e))
    new.sort(key=lambda x: (x[2].get("time") or "", x[0], x[1]))
    more = len(new) > limit
    out_pos = dict(pos)
    if more:
        # keep, per file, only the lines that follow each other from the cursor,
        # so that the new cursor (a line count per file) skips nothing
        head = new[:limit]
        for key in {k for k, _, _ in head}:
            run = pos.get(key, 0)
            for n in sorted(n for k, n, _ in head if k == key):
                if n != run:
                    break
                run += 1
            out_pos[key] = run
        new = [x for x in head if x[1] < out_pos.get(x[0], 0)]
        if not new:                                # cannot happen for time-ordered files
            new, more = head[:1], True
            out_pos[new[0][0]] = new[0][1] + 1
    else:
        out_pos.update(counts)
    events = [e for _, _, e in new if allow_lvk or e.get("visibility") != "lvk"]
    return events, encode_cursor(out_pos), more


# --------------------------------------------------------------------------
# Listener heartbeat
# --------------------------------------------------------------------------

class Heartbeat:
    """``<store>/status/<source>.json`` of a long-lived listener."""

    def __init__(self, root: str | Path, source: str, journal: Optional[Journal] = None,
                 topics: Iterable[str] = (), interval_s: float = HEARTBEAT_S):
        self.path = Path(root) / "status" / f"{_safe(source)}.json"
        self.source = source
        self.journal = journal
        self.interval_s = interval_s
        self._last_write = 0.0
        self._last_error_text = ""
        self._last_error_logged = 0.0
        prev = {}
        try:
            prev = json.loads(self.path.read_text())
        except (OSError, ValueError):
            pass
        now = utcnow()
        self.state = {"source": source, "host": socket.gethostname(), "pid": os.getpid(),
                      "started": iso(now), "heartbeat": iso(now), "topics": list(topics),
                      "last_message": dict(prev.get("last_message") or {}),
                      "messages": 0, "skipped": 0, "failed": 0, "kafka_errors": 0,
                      "last_error": None, "restarts": int(prev.get("restarts", -1)) + 1}
        down = None
        if prev.get("heartbeat"):
            try:
                down = (now - parse_iso(prev["heartbeat"])).total_seconds()
            except ValueError:
                pass
        if journal is not None:
            journal.log("listener_start",
                        f"{source} started on {self.state['host']}"
                        + (f" after {_duration(down)} without heartbeat" if down is not None else ""),
                        level="warning" if down is not None and down > 600 else "info",
                        host=self.state["host"], pid=self.state["pid"], topics=list(topics),
                        previous_heartbeat=prev.get("heartbeat"), down_s=down,
                        restarts=self.state["restarts"])
        self.write(force=True)

    def write(self, force: bool = False) -> None:
        t = time.monotonic()
        if not force and t - self._last_write < self.interval_s:
            return
        self._last_write = t
        self.state["heartbeat"] = iso(utcnow())
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            _write_atomic(self.path, json.dumps(self.state, indent=1))
        except OSError:
            pass

    def tick(self) -> None:
        self.write()

    def message(self, topic: str, skipped: bool = False) -> None:
        self.state["last_message"][topic] = iso(utcnow())
        self.state["messages"] += 1
        if skipped:
            self.state["skipped"] += 1
        self.write()

    def kafka_error(self, text: str) -> None:
        """Counted always; written in the journal when the text changes, or
        every 10 minutes while the same error repeats."""
        self.state["kafka_errors"] += 1
        self.state["last_error"] = {"time": iso(utcnow()), "text": text[:500]}
        t = time.monotonic()
        if self.journal is not None and (text != self._last_error_text
                                         or t - self._last_error_logged > 600):
            self._last_error_text, self._last_error_logged = text, t
            self.journal.log("kafka_error", f"{self.source}: {text[:200]}", level="warning",
                             count=self.state["kafka_errors"])
        self.write()

    def failure(self, topic: str, exc: BaseException) -> None:
        self.state["failed"] += 1
        self.state["last_error"] = {"time": iso(utcnow()), "text": f"{topic}: {exc!r}"[:500]}
        if self.journal is not None:
            self.journal.log("message_failed", f"{self.source}: a {topic} message could not be handled "
                             f"({type(exc).__name__}: {exc})"[:300], level="error", topic=topic,
                             traceback=traceback.format_exc()[-2000:])
        self.write(force=True)


def _duration(s: Optional[float]) -> str:
    if s is None:
        return "?"
    if s < 120:
        return f"{s:.0f} s"
    if s < 7200:
        return f"{s / 60:.0f} min"
    if s < 172800:
        return f"{s / 3600:.1f} h"
    return f"{s / 86400:.1f} d"


# --------------------------------------------------------------------------
# Search / PE job state
# --------------------------------------------------------------------------

class RunLog:
    """State and journal of one search or PE job, in its output directory.

    ``with RunLog(out, optical, run) as rl:`` marks the job running, keeps
    the heartbeat alive from a thread (so a job killed without a word is seen
    as stale), and records the exception if the block fails."""

    def __init__(self, out_dir: str | Path, optical: str, run: str, visibility: str = "lvk",
                 kind: str = "search", clock: Callable[[], datetime] = utcnow,
                 heartbeat_s: float = RUN_HEARTBEAT_S):
        self.dir = Path(out_dir)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.clock = clock
        self.heartbeat_s = heartbeat_s
        # files named after the run: several jobs may share one output directory
        self.state_path = self.dir / f"{_safe(run)}_state.json"
        self.journal = Journal(self.dir / f"{_safe(run)}_journal.jsonl", f"run:{run}", clock=clock,
                               defaults={"run": run, "visibility": visibility})
        self.optical = optical                  # name of the optical transient the run is for
        prev = {}
        try:
            prev = json.loads(self.state_path.read_text())
        except (OSError, ValueError):
            pass
        self.state = {"optical": optical, "run": run, "kind": kind, "visibility": visibility,
                      "stage": "starting", "submitted": prev.get("submitted"),
                      "started": None, "heartbeat": None, "ended": None,
                      "host": socket.gethostname(), "subsegments_done": 0,
                      "subsegments_total": None, "loudest": None, "p_value": None,
                      "pe": dict(prev.get("pe") or {}), "data_source": None,
                      "exit": None, "error": None, "note": prev.get("note")}
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._lock = threading.Lock()

    def _write(self) -> None:
        with self._lock:
            self.state["heartbeat"] = iso(utcnow())
            try:
                _write_atomic(self.state_path, json.dumps(_clean(self.state), indent=1))
            except OSError:
                pass

    def log(self, kind: str, summary: str, level: str = "info", **data) -> None:
        self.journal.log(kind, summary, level=level, optical=self.optical, **data)

    def submitted(self, note: str = "") -> None:
        self.state.update(stage="submitted", submitted=iso(self.clock()), note=note or None)
        self._write()
        self.log("run_submitted", f"{self.state['kind']} {self.state['run']} submitted"
                 + (f": {note}" if note else ""), note=note)

    def stage(self, stage: str, summary: str = "", **data) -> None:
        self.state["stage"] = stage
        self.state.update({k: v for k, v in data.items() if k in self.state})
        self._write()
        self.log("run_stage", summary or f"{self.state['run']}: {stage}", stage=stage, **data)

    def subsegment(self, done: int, total: int, loudest: Optional[dict] = None) -> None:
        self.state.update(stage="searching", subsegments_done=done, subsegments_total=total)
        if loudest:
            best = self.state["loudest"]
            if not best or loudest.get("rho_net", 0) > best.get("rho_net", 0):
                self.state["loudest"] = loudest
        self._write()
        self.log("run_progress", f"{self.state['run']}: sub-segment {done} of {total}"
                 + (f", loudest {loudest['rho_net']:.2f}" if loudest and "rho_net" in loudest else ""),
                 done=done, total=total, loudest=loudest)

    def pe(self, which: str, stage: str, summary: Optional[dict] = None,
           warnings: Optional[list[str]] = None) -> None:
        self.state["stage"] = "pe"
        self.state["pe"][which] = {"stage": stage, "time": iso(self.clock()),
                                   "warnings": warnings or []}
        self._write()
        if stage == "done":
            self.log("pe_result", f"{self.state['run']}: PE ({which} priors) finished"
                     + (f", {len(warnings)} warning(s): " + "; ".join(warnings) if warnings else ""),
                     level="warning" if warnings else "info", which=which, pe_summary=summary,
                     warnings=warnings or [])
        else:
            self.log("run_stage", f"{self.state['run']}: PE ({which} priors) {stage}",
                     stage=f"pe_{which}_{stage}")

    def __enter__(self) -> "RunLog":
        self.state.update(stage="starting", started=iso(self.clock()), exit=None, error=None, ended=None)
        self._write()
        self.log("run_start", f"{self.state['kind']} {self.state['run']} started on {self.state['host']}",
                 host=self.state["host"])

        def beat():
            while not self._stop.wait(self.heartbeat_s):
                self._write()
        self._thread = threading.Thread(target=beat, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, etype, exc, tb) -> bool:
        self._stop.set()
        self.state["ended"] = iso(self.clock())
        if etype is None or etype is SystemExit and not getattr(exc, "code", 0):
            self.state.update(stage="done", exit="ok")
            self._write()
            self.log("run_end", f"{self.state['run']} finished", p_value=self.state.get("p_value"),
                     loudest=self.state.get("loudest"))
        else:
            text = "".join(traceback.format_exception_only(etype, exc)).strip()
            self.state.update(stage="failed", exit="failed", error=text[-2000:])
            self._write()
            self.log("run_failed", f"{self.state['run']} failed: {text[-300:]}", level="error",
                     error=text[-2000:], traceback="".join(traceback.format_tb(tb))[-3000:])
        return False


def pe_warnings(summary: dict, search_snr: Optional[float] = None) -> list[str]:
    """Checks of a PE summary before it is trusted (see TODO.md, 2026-11-18):
    posterior on a prior edge, Bayes factor above what the SNR allows,
    maximum-likelihood SNR far from the search value."""
    out = []
    if summary.get("distance_at_prior_edge"):
        out.append(f"distance posterior on the {summary['distance_at_prior_edge']} prior bound")
    ml = (summary.get("maximum_likelihood") or {}).get("snr")
    lbf = summary.get("log_bayes_factor")
    ref = search_snr or ml
    if lbf is not None and ref and lbf > 0.5 * ref ** 2 * 1.05 + 5:
        out.append(f"log Bayes factor {lbf:.0f} above the limit of an SNR of {ref:.1f} "
                   f"({0.5 * ref ** 2:.0f})")
    if ml and search_snr and abs(ml - search_snr) > 0.2 * search_snr:
        out.append(f"maximum-likelihood SNR {ml:.1f} against {search_snr:.1f} in the search")
    return out


# --------------------------------------------------------------------------
# Sky-map versions
# --------------------------------------------------------------------------

def credible_levels(prob: np.ndarray) -> np.ndarray:
    """Credible level of each pixel: the probability enclosed by the pixels
    at least as probable."""
    p = np.asarray(prob, dtype=float)
    order = np.argsort(-p, kind="stable")
    cum = np.cumsum(p[order])
    if cum[-1] > 0:
        cum = cum / cum[-1]
    out = np.empty_like(cum)
    out[order] = cum
    return out


def reduce_map(moc) -> tuple[bytes, dict]:
    """Multi-order map -> (gzip of one byte per pixel at nside 128 NESTED,
    summary). Byte = ceil(credible level x 250), 255 where the map is empty."""
    import healpy as hp
    from .listener import moc_to_flat
    flat = moc_to_flat(moc)                                  # up to nside 512: areas
    cred = credible_levels(flat["PROB"])
    npix = len(cred)
    area = 4 * 180 ** 2 / np.pi / npix
    nside = hp.npix2nside(npix)
    imax = int(np.argmax(flat["PROB"]))
    ra, dec = hp.pix2ang(nside, imax, nest=True, lonlat=True)
    summary = {"area50_deg2": float(np.sum(cred <= 0.5) * area),
               "area90_deg2": float(np.sum(cred <= 0.9) * area),
               "max_ra": float(ra), "max_dec": float(dec), "nside": 1 << SKY_ORDER,
               "nested": True, "scale": CRED_SCALE}
    small = moc_to_flat(moc, max_order=SKY_ORDER)["PROB"]
    if len(small) < 12 * 4 ** SKY_ORDER:                     # coarse map: repeat pixels
        small = np.repeat(small / (12 * 4 ** SKY_ORDER // len(small)), 12 * 4 ** SKY_ORDER // len(small))
    c = credible_levels(small)
    b = np.minimum(np.ceil(c * CRED_SCALE), CRED_SCALE).astype(np.uint8)
    b[small <= 0] = 255
    return gzip.compress(b.tobytes(), 6), summary


def save_sky_version(root: str | Path, alert: dict, moc, combined=None) -> Optional[dict]:
    """Keep the reduced map(s) of one GW alert. Returns the version record
    (also appended to ``gw_versions/<superevent>/index.json``)."""
    if moc is None and combined is None:
        return None
    sid = alert["superevent_id"]
    d = Path(root) / "gw_versions" / _safe(sid)
    d.mkdir(parents=True, exist_ok=True)
    idx_path = d / "index.json"
    try:
        index = json.loads(idx_path.read_text())
    except (OSError, ValueError):
        index = []
    n = len(index) + 1
    rec = {"version": n, "alert_type": alert.get("alert_type"), "time_created": alert.get("time_created"),
           "far": alert.get("far"), "gpstime": alert.get("gpstime"), "received": iso(utcnow()),
           # LVK's own GW x GRB coincidence, kept with the version: the alert file is deleted later
           "external_coinc": alert.get("external_coinc"),
           "classification": alert.get("classification"), "properties": alert.get("properties"),
           "instruments": alert.get("instruments"), "maps": {}}
    for kind, m in (("gw", moc), ("combined", combined)):
        if m is None:
            continue
        data, summary = reduce_map(m)
        name = f"{n:02d}_{kind}.cred.gz"
        (d / name).write_bytes(data)
        rec["maps"][kind] = {"file": name, **summary}
    index.append(rec)
    _write_atomic(idx_path, json.dumps(_clean(index), indent=1))
    return rec


def sky_versions(root: str | Path, superevent_id: str) -> list[dict]:
    try:
        return json.loads((Path(root) / "gw_versions" / _safe(superevent_id) / "index.json").read_text())
    except (OSError, ValueError):
        return []


# --------------------------------------------------------------------------
# Decisions
# --------------------------------------------------------------------------

def decision_digest(d: Optional[dict]) -> Optional[tuple]:
    """What makes two decisions different for the logbook."""
    if not d:
        return None
    ev, det = d.get("kilonova_evidence") or {}, d.get("detectability") or {}
    return (bool(d.get("trigger_targeted_search")), bool(ev.get("confirmed")),
            tuple(ev.get("by") or ()), bool(det.get("worth_search")),
            (d.get("t0_window") or {}).get("method"))


def log_decision(journal: Any, decision: dict, previous: Optional[dict] = None,
                 always: bool = False) -> Optional[dict]:
    """Write the decision in the journal if it differs from the previous one
    (trigger, kilonova evidence, detectability, window method)."""
    if not always and decision_digest(decision) == decision_digest(previous):
        return None
    ev, det = decision.get("kilonova_evidence") or {}, decision.get("detectability") or {}
    trig = bool(decision.get("trigger_targeted_search"))
    was = None if previous is None else bool(previous.get("trigger_targeted_search"))
    name = decision.get("optical") or decision.get("candidate")      # "candidate": the key before the rename
    if trig:
        text = f"{name}: targeted search triggered"
    else:
        why = []
        if not ev.get("confirmed") and ev.get("required", True):
            why.append("not confirmed as a kilonova")
        if not det.get("worth_search"):
            why.append("not detectable")
        text = f"{name}: held ({', '.join(why) or 'no trigger'})"
    data = {k: decision.get(k) for k in ("t0_window", "distance_inclination", "search_priors",
                                         "detectability", "kilonova_evidence",
                                         "trigger_targeted_search", "files")}
    return journal.log("decision", text, optical=name, trigger=trig, trigger_changed=was is not None and was != trig,
                       first=previous is None, **{k: v for k, v in data.items() if k != "trigger_targeted_search"})
