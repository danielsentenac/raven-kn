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

With a **farm** (``raven-kn launcher --farm farmn1,farmn2,...`` or
``RAVEN_KN_FARM``), the search does not run on the launcher's machine but on
the least loaded node of the list, started over ``ssh`` (:func:`pick_node`).
The nodes must be reached without a password and see the same store, the same
Python and the same frames at the same paths. When no node qualifies the
search runs on the launcher's machine, and its note says so.

A window longer than one sub-segment of the search (1024 s; a window from
the light curve alone lasts hours to days) is cut in **pieces**, one
sub-segment each (:func:`pieces`). Each piece is checked for data on its own
and skipped, with the reason, where fewer than two detectors have any; the
others are searched as separate jobs, on as many nodes of the farm at once as
``max_nodes`` allows (one after the other on the launcher's machine without a
farm). A piece that fails is tried once more, on another node. At the end the
results of the pieces are merged in one result (their p-values combined with
the trials factor, :func:`raven_kn.targeted.search.merge_results`): the
application shows one run, its progress in pieces, and the pieces skipped or
failed.

Status of a request: ``requested``, ``running``, ``done``, ``failed``,
``refused``.
"""

from __future__ import annotations

import json
import logging
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Optional, Sequence

from .journal import Heartbeat, Journal, RunLog, _safe, _write_atomic, iso, journal_of, utcnow

log = logging.getLogger("raven_kn.launch")

DATA = ("site", "replay")
STATISTICS = ("coherent", "network")
ACTIVE = ("requested", "running")
F_LOW = 30.0
TAIL_S = 32.0                    # read after the window: 4 x the padding of the search
LEAD_FALLBACK_S = 110.0          # before it, when the template duration cannot be computed
NODE_MAX_LOAD = 1.0              # a farm node is taken when its 1-minute load is below this
NODE_MIN_MEM_GB = 10.0           # and it has this much memory available (6 processes of 1.5 GB)
PIECE_S = 1024.0                 # a piece of a long window: one sub-segment of the search
MIN_PIECE_S = 16.0               # a leftover shorter than this at the end of a long window is dropped
MAX_NODES = 10                   # farm nodes a search uses at once
NODE_WAIT_S = 1800.0             # pieces waiting this long with no node free and none running: failed
SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", "-o", "ServerAliveInterval=30"]


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


def farm_hosts(text: Optional[str]) -> list[str]:
    """Host names of a farm given as ``farmn1,farmn2`` (commas or spaces)."""
    return [h for h in (text or "").replace(",", " ").split() if h]


def _remote(script: str) -> str:
    """The remote command running ``script`` with bash, whatever the login shell."""
    return "/bin/bash -c " + shlex.quote(script)


def _probe(host: str, paths: Sequence[str]) -> Optional[dict]:
    """Cores, load and available memory of a node, or None when it is not
    reached without a password or does not see one of ``paths``."""
    script = ("; ".join(f"test -e {shlex.quote(str(p))} || exit 3" for p in paths) + "; " if paths else "") + (
        "echo $(nproc) $(cut -d' ' -f1 /proc/loadavg) $(awk '/MemAvailable/{print $2}' /proc/meminfo)")
    try:
        out = subprocess.run(SSH + [host, _remote(script)], capture_output=True, text=True, timeout=30)
        cores, load, mem_kb = out.stdout.split()[-3:]
        if out.returncode != 0:
            return None
        return {"host": host, "cores": int(cores), "load": float(load), "mem_gb": float(mem_kb) / 1048576.0}
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def pick_node(hosts: Sequence[str], paths: Sequence[str] = (),
              probe: Callable[..., Optional[dict]] = _probe) -> tuple[Optional[dict], str]:
    """The node of ``hosts`` to run a search on: among those reached, seeing
    ``paths``, with a load below :data:`NODE_MAX_LOAD` and
    :data:`NODE_MIN_MEM_GB` of memory available, the least loaded (then the
    one with the most memory). Returns (node or None, why none)."""
    if not hosts:
        return None, "no farm configured"
    with ThreadPoolExecutor(max_workers=min(20, len(hosts))) as ex:
        nodes = [n for n in ex.map(lambda h: probe(h, paths), hosts) if n]
    if not nodes:
        return None, f"none of the {len(hosts)} farm nodes reached"
    free = [n for n in nodes if n["load"] < NODE_MAX_LOAD and n["mem_gb"] >= NODE_MIN_MEM_GB]
    if not free:
        return None, f"the {len(nodes)} farm nodes reached are all busy"
    return min(free, key=lambda n: (n["load"], -n["mem_gb"])), ""


def _workers(cores: int) -> int:
    """Processes of a search on a machine of ``cores`` cores: two are left free."""
    return max(1, int(cores) - 2)


def pieces(gps_start: float, gps_end: float) -> list[tuple[float, float]]:
    """The pieces a window is searched in: itself if it is no longer than
    :data:`PIECE_S`, else the sub-segments ``raven-kn targeted`` would cut
    (a leftover shorter than :data:`MIN_PIECE_S` dropped, as it does)."""
    if gps_end - gps_start <= PIECE_S:
        return [(gps_start, gps_end)]
    out, a = [], gps_start
    while a < gps_end:
        b = min(a + PIECE_S, gps_end)
        if b - a >= MIN_PIECE_S:
            out.append((a, b))
        a = b
    return out


def _command(req_path: Path, data: str, ifos: Sequence[str], statistic: str, label: str, out: Path,
             cores: int, host: Optional[str] = None) -> list[str]:
    """The command of a search job, on this machine or on ``host`` by ssh."""
    cmd = [sys.executable, "-m", "raven_kn.cli", "targeted", "--request", str(req_path),
           "--source", data, "--ifos", ",".join(ifos), "--statistic", statistic,
           "--label", label, "--out", str(out), "--visibility", "lvk", "--workers", str(_workers(cores))]
    if host is None:
        return cmd
    libs, pp = os.environ.get("RAVEN_KN_LIBS", ""), os.environ.get("PYTHONPATH", "")
    script = (f"cd {shlex.quote(str(out))} && "
              + (f"export LD_LIBRARY_PATH={shlex.quote(libs)}:$LD_LIBRARY_PATH && " if libs else "")
              + (f"export PYTHONPATH={shlex.quote(pp)} && " if pp else "")     # the launcher's code
              + "exec " + " ".join(shlex.quote(c) for c in cmd))
    return SSH + [host, _remote(script)]


def _spawn(cmd: list[str], log_path: Path):
    """Start a job; returns an object with ``poll()`` (None while it runs)."""
    with open(log_path, "ab") as f:
        return subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)


def _run_pieces(rec: dict, root: Path, run_dir: Path, name: str, label: str, jobs: list[dict],
                skipped: list[dict], note: str, window: tuple[float, float], farm: Sequence[str],
                picker: Callable[..., tuple], spawn: Callable[..., Any],
                tick: Optional[Callable[[], None]], max_nodes: int, poll_s: float) -> tuple[str, str, dict]:
    """Search the pieces ``jobs`` of a window as separate jobs and merge their
    results (module documentation). Returns (status, detail, journal data)."""
    from .targeted.search import merge_results
    total = len(jobs)
    pending, running, done, failed = list(jobs), {}, [], []
    hosts: list[str] = []
    waited = 0.0

    def progress() -> None:
        rec.update(pieces_total=total + len(skipped), pieces_done=len(done), pieces_skipped=len(skipped),
                   pieces_failed=len(failed), pieces_running=sorted(running), hosts=hosts)
        _save(root, rec)

    with RunLog(run_dir, name, label, visibility="lvk") as rlog:
        rlog.state["data_source"] = f"site-{'replay' if rec['data'] == 'replay' else 'live'}"
        rlog.stage("searching", f"{label}: {total} pieces of {PIECE_S:.0f} s to search"
                   + (f", {len(skipped)} skipped for lack of data" if skipped else ""))
        rlog.state["subsegments_total"] = total
        while pending or running:
            if tick:
                tick()
            for slot, (proc, job) in list(running.items()):
                code = proc.poll()
                if code is None:
                    continue
                del running[slot]
                result = job["dir"] / f"{job['label']}_targeted_result.json"
                if code == 0 and result.exists():
                    done.append(job)
                    lo = json.loads(result.read_text()).get("loudest") or {}
                    rlog.subsegment(len(done), total, {"rho_net": lo.get("loudest_rho_net", 0.0),
                                                       "gps": lo.get("loudest_gps"),
                                                       "mchirp": lo.get("loudest_mchirp"),
                                                       "p_value": lo.get("p_value")})
                elif job["tries"] < 2:                     # once more, on another node
                    job.update(tries=job["tries"] + 1, avoid=slot)
                    pending.append(job)
                    rlog.log("piece_retry", f"{label}: piece {job['index']} stopped on {slot} "
                             f"(exit code {code}), tried again", level="warning", piece=job["index"],
                             host=slot, exit=code)
                else:
                    failed.append({**job, "exit": code, "host": slot})
                    rlog.log("piece_failed", f"{label}: piece {job['index']} failed twice "
                             f"(exit code {code} on {slot})", level="error", piece=job["index"],
                             host=slot, exit=code)
                progress()
            started = False
            while pending and len(running) < (max_nodes if farm else 1):
                job = pending[0]
                if farm:
                    node, _ = picker([h for h in farm if h not in running and h != job.get("avoid")],
                                     [str(job["dir"]), sys.executable])
                    if not node:
                        break                              # none free now: looked at again later
                    slot, cores = node["host"], node["cores"]
                else:
                    slot, cores = "the launcher's machine", os.cpu_count() or 1
                pending.pop(0)
                cmd = _command(job["req"], rec["data"], job["ifos"], rec["statistic"], job["label"],
                               job["dir"], cores, slot if farm else None)
                running[slot] = (spawn(cmd, job["dir"] / f"{job['label']}.log"), job)
                if slot not in hosts:
                    hosts.append(slot)
                started = True
                progress()
            if pending and not running and not started:    # nothing moves: no node is free
                waited += max(poll_s, 1.0)
                if waited > NODE_WAIT_S:
                    failed += [{**j, "exit": None, "host": None} for j in pending]
                    rlog.log("pieces_abandoned", f"{label}: {len(pending)} pieces not searched, no farm "
                             f"node free for {NODE_WAIT_S / 60:.0f} minutes", level="error")
                    pending = []
                    progress()
            else:
                waited = 0.0
            if pending or running:
                time.sleep(poll_s)
        if not done:
            raise RuntimeError(f"none of the {total} pieces could be searched")
        res = merge_results([j["dir"] / f"{j['label']}_targeted_result.json" for j in done], label,
                            run_dir, window[0], window[1],
                            skipped + [{"gps_start": j["gps_start"], "gps_end": j["gps_end"],
                                        "why": f"failed (exit code {j['exit']})"} for j in failed])
        rlog.state["p_value"] = res.p_value
        rlog.log("search_result", f"{label}: loudest {res.loudest.get('loudest_rho_net', 0):.2f}, p-value "
                 f"{res.p_value:.3g}, {len(done)} pieces", p_value=res.p_value, bank_size=res.bank_size,
                 loudest=res.loudest, result_file=f"{label}_targeted_result.json")
    detail = (f"{note}; {len(done)} of {total + len(skipped)} pieces searched on "
              f"{len(hosts)} machine{'s' if len(hosts) != 1 else ''}"
              + (f", {len(skipped)} skipped for lack of data" if skipped else "")
              + (f", {len(failed)} failed" if failed else ""))
    progress()
    return "done", detail, {"pieces_done": len(done), "pieces_skipped": len(skipped),
                            "pieces_failed": len(failed), "hosts": hosts}


def process(store: Any, rec: dict, source: Any = None, now_gps: Optional[float] = None,
            runner: Callable[..., int] = _run, tick: Optional[Callable[[], None]] = None,
            farm: Sequence[str] = (), picker: Callable[..., tuple] = pick_node,
            spawn: Callable[..., Any] = _spawn, max_nodes: int = MAX_NODES, poll_s: float = 5.0) -> dict:
    """Carry out one launch request (steps in the module documentation).
    ``source``, ``now_gps``, ``runner``, ``picker`` and ``spawn`` are replaced
    in tests. Returns the request with its final status."""
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
        run_dir = (root / "runs" / _safe(name) / label).resolve()     # the same path on a farm node
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
        elif len(pieces(float(req["gps_start"]), float(req["gps_end"]))) > 1:
            ifos = list(source.ifos())               # a long window: the data is checked piece by piece
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
        cut = pieces(float(req["gps_start"]), float(req["gps_end"]))
        if len(cut) > 1:
            jobs, skipped = [], []
            for n, (a, b) in enumerate(cut, 1):
                use = list(ifos)
                if rec["data"] != "replay":                # each piece on the detectors it has
                    cov = source.coverage(source.ifos(), a - lead, b + TAIL_S)
                    use = [i for i, c in cov.items() if c["ok"]]
                    if len(use) < 2:
                        skipped.append({"gps_start": a, "gps_end": b, "why": "; ".join(
                            f"{i}: {c['why']}" for i, c in cov.items() if not c["ok"])})
                        continue
                d = run_dir / "pieces" / f"p{n:04d}"
                d.mkdir(parents=True, exist_ok=True)
                (d / "request.json").write_text(json.dumps({**req, "gps_start": a, "gps_end": b}, indent=2))
                jobs.append({"index": n, "gps_start": a, "gps_end": b, "ifos": use, "dir": d,
                             "req": d / "request.json", "label": f"{label}_p{n:04d}", "tries": 1})
            if not jobs:
                return end("refused", f"none of the {len(cut)} pieces of the window has analysable "
                           "data in two detectors", "warning")
            note += f"; window of {(cut[-1][1] - cut[0][0]) / 3600:.1f} h searched in {len(cut)} pieces"
            RunLog(run_dir, name, label, visibility="lvk").submitted(note)
            rec.update(status="running", detail=note, started=iso(utcnow()), run=label,
                       ifos=sorted({i for j in jobs for i in j["ifos"]}),
                       gps_start=req["gps_start"], gps_end=req["gps_end"])
            _save(root, rec)
            status, detail, data = _run_pieces(rec, root, run_dir, name, label, jobs, skipped, note,
                                               (float(req["gps_start"]), float(req["gps_end"])), farm,
                                               picker, spawn, tick, max_nodes, poll_s)
            return end(status, detail, "warning" if data["pieces_failed"] else "info", **data)
        host, cores = None, os.cpu_count() or 1
        if farm:
            node, why = picker(farm, [str(run_dir), sys.executable])
            if node:
                host, cores = node["host"], node["cores"]
                note += f"; run on {host} ({node['cores']} cores)"
            else:
                note += f"; run on the launcher's machine ({why})"
        cmd = _command(req_path, rec["data"], ifos, rec["statistic"], label, run_dir, cores, host)
        RunLog(run_dir, name, label, visibility="lvk").submitted(note)
        rec.update(status="running", detail=note, started=iso(utcnow()), ifos=ifos, run=label,
                   host=host, gps_start=req["gps_start"], gps_end=req["gps_end"])
        _save(root, rec)
        code = runner(cmd, run_dir / f"{label}.log", tick)
        if code == 0:
            return end("done", note, ifos=ifos)
        return end("failed", f"the search stopped with exit code {code} (see {label}.log in its "
                   "run directory)", "error")
    except Exception as e:                           # the launcher goes on with the next request
        log.exception("launch %s failed", label)
        return end("failed", f"{type(e).__name__}: {e}"[:300], "error")


def launcher(store: Any, interval_s: float = 10.0, once: bool = False,
             farm: Sequence[str] = (), max_nodes: int = MAX_NODES) -> None:
    """Carry out the launch requests of a store, one at a time, forever (or
    until none is waiting, with ``once``), on a node of ``farm`` when given. Requests left ``running`` by a
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
            done = process(store, waiting[0], tick=hb.tick, farm=farm, max_nodes=max_nodes)
            log.info("launch %s: %s %s", done["id"], done["status"], done.get("detail", ""))
        elif once:
            return
        else:
            time.sleep(interval_s)
