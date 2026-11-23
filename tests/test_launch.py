"""Targeted searches asked for from the application: the request, the data
check of the launcher, the replay test mode, and the service's endpoint."""
import json
import threading
import urllib.error
import urllib.request

import numpy as np
import pytest

from raven_kn.connection import GRBAssociation, OpticalTransient
from raven_kn.journal import read_events
from raven_kn.launch import launches, process, request_search
from raven_kn.listener import AlertStore, add_optical
from raven_kn.serve import Store, make_server
from raven_kn.targeted.datasource import SiteSource

T0 = 1479043087.0


@pytest.fixture()
def store(tmp_path):
    s = AlertStore(tmp_path / "store", source="t")
    add_optical(OpticalTransient(name="KN1", ra=10.0, dec=20.0, source="GCN circular 1",
                                 grb=GRBAssociation(gps=T0, instrument="Fermi-GBM", graceid="1")), s,
                now_gps=T0 + 10)
    return s


def site(tmp_path, observing, t_files=(T0 - 2000, T0 + 2000)):
    """A site with one 4000 s file per detector; ``observing``: detector ->
    list of (start, end) where its state vector says observing."""
    layout = {}
    for ifo in observing:
        d = tmp_path / "frames" / ifo
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{ifo[0]}-Test-{int(t_files[0])}-{int(t_files[1] - t_files[0])}.gwf").write_bytes(b"")
        (d / f"{ifo[0]}-Test-{int(t_files[1])}-4000.gwf_NOT_YET_CLOSED").write_bytes(b"")
        layout[ifo] = {"dir": str(d), "channel": f"{ifo}:STRAIN", "state": f"{ifo}:STATE", "mask": 3}
    src = SiteSource.named("live", layout=layout)

    def state(path, channel, start, end):
        ifo = channel.split(":")[0]
        t = np.arange(start, end, 1.0)
        v = np.full(t.size, 496, dtype="int64")                       # h(t) not OK, not observing
        for a, b in observing[ifo]:
            v[(t >= a) & (t < b)] = 0b111                              # bits 0 and 1 set
        return float(start), 1.0, v
    src._state_samples = state
    return src


def test_site_source(tmp_path):
    src = site(tmp_path, {"H1": [(T0 - 500, T0 + 500)], "L1": [(T0 - 500, T0 + 20)], "V1": []})
    assert [(a, b) for a, b, _ in src.files("H1", T0 - 10, T0 + 10)] == [(T0 - 2000, T0 + 2000)]
    assert src.files("H1", T0 + 2100, T0 + 2200) == []                 # the file still being written
    assert src.analysable("H1", T0 - 600, T0 + 600) == [(T0 - 500, T0 + 500)]
    assert src.available("H1", T0 - 100, T0 + 32) and not src.available("L1", T0 - 100, T0 + 32)
    cov = src.coverage(["H1", "L1", "V1", "K1"], T0 - 100, T0 + 32)
    assert cov["H1"] == {"ok": True, "why": "", "analysable": [(T0 - 100, T0 + 32)]}
    assert cov["L1"]["why"] == "observing for 120 s of 132 s only"
    assert cov["V1"]["why"] == "not in observing mode (state vector)"
    assert cov["K1"]["why"] == "no frame file for that time"
    assert src.frame_urls("H1", T0, T0 + 1)[0].startswith("file://") and src.channel("H1") == "H1:STRAIN"
    # the latest stretch analysable in the most detectors: H1 and L1 until T0 + 20
    start, ifos = src.latest_window(6.0, 100.0, 32.0, now=T0 + 1000)
    assert ifos == ["H1", "L1"] and start == T0 + 20 - 32 - 6
    assert src.latest_window(6.0, 100.0, 32.0, now=T0 + 1000, min_ifos=3) is None


def test_request(store):
    with pytest.raises(ValueError):
        request_search(store, "KN1", data="elsewhere")
    with pytest.raises(ValueError):
        request_search(store, "unknown")
    rec = request_search(store, "KN1", "replay", "network")
    assert rec["status"] == "requested" and rec["test"] and rec["id"].startswith("KN1_replay_network_")
    assert launches(store.root, "KN1") == [rec] and launches(store.root, "other") == []
    with pytest.raises(ValueError, match="already waiting or running"):
        request_search(store, "KN1")
    assert [e["summary"] for e in read_events(store.root)[0] if e["kind"] == "search_requested"] == [
        "KN1: targeted search asked for from the application (replay data, a test, network statistic)"]


def test_refused_without_two_detectors(store, tmp_path):
    src = site(tmp_path, {"H1": [(T0 - 500, T0 + 500)], "L1": [], "V1": []})
    rec = process(store, request_search(store, "KN1"), source=src,
                  runner=lambda *a: pytest.fail("no search must run"))
    assert rec["status"] == "refused"
    assert rec["detail"] == ("fewer than two detectors with analysable data (L1: not in observing mode "
                             "(state vector); V1: not in observing mode (state vector))")
    e = [e for e in read_events(store.root)[0] if e["kind"] == "search_refused"]
    assert len(e) == 1 and e[0]["level"] == "warning"
    assert not list((store.root / "runs" / "KN1" / rec["id"]).glob("*_state.json"))   # no run was started
    request_search(store, "KN1")                                       # a refused one does not block another


def test_run_on_site_data(store, tmp_path):
    src = site(tmp_path, {"H1": [(T0 - 500, T0 + 500)], "L1": [(T0 - 500, T0 + 500)], "V1": []})
    seen = {}

    def runner(cmd, log_path, tick):
        seen["cmd"] = cmd
        return 0
    rec = process(store, request_search(store, "KN1", "site", "coherent"), source=src, runner=runner)
    assert rec["status"] == "done" and rec["ifos"] == ["H1", "L1"]
    cmd = seen["cmd"]
    assert cmd[1:4] == ["-m", "raven_kn.cli", "targeted"]
    opt = dict(zip(cmd[4::2], cmd[5::2]))
    assert (opt["--source"], opt["--ifos"], opt["--statistic"], opt["--label"], opt["--visibility"]) == (
        "site", "H1,L1", "coherent", rec["id"], "lvk")
    run_dir = store.root / "runs" / "KN1" / rec["id"]
    assert opt["--out"] == str(run_dir)
    req = json.loads(open(opt["--request"]).read())
    assert (req["gps_start"], req["gps_end"]) == (T0 - 5, T0 + 1)      # the GRB window, as it is
    assert not list(run_dir.glob("*_skymap.fits"))
    # the run is in the logbook from the start, with who asked for it
    state = json.loads((run_dir / f"{rec['id']}_state.json").read_text())
    assert state["stage"] == "submitted" and "V1 (not in observing mode" in state["note"]
    view = Store(store.root).optical_detail("KN1", True)
    assert [r["run"] for r in view["runs"]] == [rec["id"]]
    assert [l["status"] for l in view["launches"]] == ["done"]
    # a search that stops with an error
    bad = process(store, request_search(store, "KN1"), source=src, runner=lambda *a: 3)
    assert bad["status"] == "failed" and "exit code 3" in bad["detail"]


def test_run_on_a_farm_node(store, tmp_path, monkeypatch):
    from raven_kn.launch import farm_hosts, pick_node
    assert farm_hosts("farmn1, farmn2  farmn3") == ["farmn1", "farmn2", "farmn3"] and farm_hosts(None) == []
    nodes = {"a": {"host": "a", "cores": 8, "load": 0.4, "mem_gb": 13.0},
             "b": {"host": "b", "cores": 8, "load": 0.1, "mem_gb": 3.0},      # too little memory
             "c": {"host": "c", "cores": 8, "load": 0.2, "mem_gb": 11.0},
             "d": {"host": "d", "cores": 8, "load": 15.0, "mem_gb": 14.0},    # busy
             "e": None}                                                       # not reached
    probe = lambda h, paths: nodes[h]
    assert pick_node(list(nodes), probe=probe)[0]["host"] == "c"
    assert pick_node(["b", "d"], probe=probe) == (None, "the 2 farm nodes reached are all busy")
    assert pick_node(["e"], probe=probe) == (None, "none of the 1 farm nodes reached")

    src = site(tmp_path, {"H1": [(T0 - 500, T0 + 500)], "L1": [(T0 - 500, T0 + 500)]})
    monkeypatch.setenv("RAVEN_KN_LIBS", "/soft/lib")
    seen = {}
    rec = process(store, request_search(store, "KN1"), source=src, farm=["a", "c"],
                  picker=lambda hosts, paths: (seen.update(paths=paths) or nodes["c"], ""),
                  runner=lambda cmd, lp, tick: seen.update(cmd=cmd) or 0)
    run_dir = store.root / "runs" / "KN1" / rec["id"]
    assert rec["status"] == "done" and rec["host"] == "c" and "run on c (8 cores)" in rec["detail"]
    assert seen["paths"][0] == str(run_dir)                           # the node must see the run directory
    cmd = seen["cmd"]
    assert cmd[0] == "ssh" and cmd[-2] == "c" and cmd[-1].startswith("/bin/bash -c ")
    import shlex
    script = shlex.split(cmd[-1])[2]
    assert script.startswith(f"cd {run_dir} && export LD_LIBRARY_PATH=/soft/lib:$LD_LIBRARY_PATH && exec ")
    assert "-m raven_kn.cli targeted --request" in script and "--ifos H1,L1" in script
    # no node free: on the launcher's machine, and the note says why
    rec = process(store, request_search(store, "KN1"), source=src, farm=["d"],
                  picker=lambda hosts, paths: (None, "the 1 farm nodes reached are all busy"),
                  runner=lambda cmd, lp, tick: seen.update(cmd=cmd) or 0)
    assert rec["host"] is None and "launcher's machine (the 1 farm nodes reached are all busy)" in rec["detail"]
    assert seen["cmd"][1:4] == ["-m", "raven_kn.cli", "targeted"]


def test_long_window_in_pieces(store, tmp_path, monkeypatch):
    """A window longer than one piece: each piece checked for data and run as
    its own job on the farm, a failed one tried again elsewhere, one result."""
    import raven_kn.launch as L
    monkeypatch.setattr(L, "PIECE_S", 2.0)                             # the 6 s window: 3 pieces
    monkeypatch.setattr(L, "MIN_PIECE_S", 0.5)
    assert L.pieces(0.0, 5.2) == [(0.0, 2.0), (2.0, 4.0), (4.0, 5.2)] and L.pieces(0.0, 2.0) == [(0.0, 2.0)]
    assert L.pieces(0.0, 4.2) == [(0.0, 2.0), (2.0, 4.0)]              # the leftover dropped
    # H1 stops observing before what the last piece reads: that piece has one detector only
    src = site(tmp_path, {"H1": [(T0 - 500, T0 + L.TAIL_S)], "L1": [(T0 - 500, T0 + 500)]})
    calls = []

    class Job:
        def __init__(self, cmd, log_path):
            script = cmd[-1] if cmd[0] == "ssh" else " ".join(cmd)
            import shlex
            words = shlex.split(shlex.split(script)[2].split("exec ", 1)[1]) if cmd[0] == "ssh" else cmd
            opt = dict(zip(words[4::2], words[5::2]))
            self.host = cmd[-2] if cmd[0] == "ssh" else None
            calls.append((self.host, opt["--label"]))
            n = int(opt["--label"][-4:])
            self.code = 1 if (n == 1 and self.host == "a") else 0     # piece 1 fails on node a
            if self.code == 0:
                req = json.loads(open(opt["--request"]).read())
                sub = {"gps_start": req["gps_start"], "gps_end": req["gps_end"], "ifos": opt["--ifos"].split(","),
                       "loudest_rho_net": 5.0 + n, "loudest_gps": req["gps_start"] + 1, "loudest_mchirp": 1.2,
                       "loudest_rho_ifo": {"H1": 4.0, "L1": 4.0}, "background": [6.5] * 99, "p_value": 1.0 / n}
                (tmp_path / "x").write_text("")
                from pathlib import Path
                Path(opt["--out"], f"{opt['--label']}_targeted_result.json").write_text(json.dumps(
                    {"config": {"label": opt["--label"]}, "bank_size": 7, "subsegments": [sub], "p_value": 1.0 / n,
                     "loudest": {k: v for k, v in sub.items() if k != "background"}, "data": {"source": "site-live"}}))

        def poll(self):
            return self.code

    free = {"a": {"host": "a", "cores": 8, "load": 0.0, "mem_gb": 14.0},
            "b": {"host": "b", "cores": 8, "load": 0.0, "mem_gb": 14.0}}
    picker = lambda hosts, paths: ((free[hosts[0]], "") if hosts else (None, "none"))
    rec = process(store, request_search(store, "KN1"), source=src, farm=["a", "b"], picker=picker,
                  spawn=Job, poll_s=0.0)
    assert rec["status"] == "done", rec["detail"]
    assert (rec["pieces_done"], rec["pieces_skipped"], rec["pieces_failed"]) == (2, 1, 0)
    assert "2 of 3 pieces searched on 2 machines, 1 skipped for lack of data" in rec["detail"]
    assert calls[0] == ("a", f"{rec['id']}_p0001") and ("b", f"{rec['id']}_p0001") in calls   # tried again on b
    run_dir = store.root / "runs" / "KN1" / rec["id"]
    res = json.loads((run_dir / f"{rec['id']}_targeted_result.json").read_text())
    assert [s["gps_start"] for s in res["subsegments"]] == [T0 - 5, T0 - 3]
    assert res["loudest"]["loudest_rho_net"] == 7.0 and res["config"]["gps_end"] == T0 + 1
    assert res["p_value"] == pytest.approx(1 - (1 - 0.5) ** 2)         # trials factor over the pieces
    assert res["data"]["pieces"] == 2 and "H1" in res["data"]["skipped"][0]["why"]
    state = json.loads((run_dir / f"{rec['id']}_state.json").read_text())
    assert (state["stage"], state["subsegments_done"], state["subsegments_total"]) == ("done", 2, 2)
    view = Store(store.root).optical_detail("KN1", True)
    assert [r["run"] for r in view["runs"]] == [rec["id"]]            # one run, not one per piece
    # without a farm: the pieces one after the other on the launcher's machine
    calls.clear()
    rec = process(store, request_search(store, "KN1"), source=src, spawn=Job, poll_s=0.0)
    assert rec["status"] == "done" and [h for h, _ in calls] == [None, None]


def test_replay_moves_the_window(store, tmp_path):
    now = T0 + 90000
    src = site(tmp_path, {"H1": [(now - 3000, now - 100)], "L1": [(now - 3000, now - 100)],
                          "V1": [(now - 3000, now - 400)]}, t_files=(now - 4000, now))
    seen = {}
    rec = process(store, request_search(store, "KN1", "replay"), source=src, now_gps=now,
                  runner=lambda cmd, lp, tick: seen.update(cmd=cmd) or 0)
    assert rec["status"] == "done" and rec["ifos"] == ["H1", "L1", "V1"] and rec["test"]
    req = json.loads(open(dict(zip(seen["cmd"][4::2], seen["cmd"][5::2]))["--request"]).read())
    assert req["gps_end"] - req["gps_start"] == 6.0 and req["original_gps_start"] == T0 - 5
    assert req["gps_end"] == now - 400 - 32                            # the three detectors until there
    assert "TEST on replayed O4 data" in rec["detail"]
    # no replay data at all
    empty = site(tmp_path / "e", {"H1": [], "L1": []}, t_files=(now - 4000, now))
    r2 = process(store, request_search(store, "KN1", "replay"), source=empty, now_gps=now)
    assert r2["status"] == "refused" and "no replay data" in r2["detail"]


def test_endpoint(store):
    srv = make_server(store.root, "127.0.0.1", 0, token="pub", lvk_token="lvk", virgo_search_token="vts")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_address[1]}"

    def post(path, body, headers):
        req = urllib.request.Request(base + path, data=json.dumps(body).encode(), method="POST", headers=headers)
        try:
            with urllib.request.urlopen(req) as r:
                return r.status, json.loads(r.read())
        except urllib.error.HTTPError as e:
            return e.code, json.loads(e.read())
    try:
        ok = {"Authorization": "Bearer lvk", "X-Raven-KN-Virgo-Search": "vts"}
        assert post("/api/optical/KN1/search", {}, {})[0] == 401
        assert post("/api/optical/KN1/search", {}, {"Authorization": "Bearer lvk"})[0] == 403
        assert post("/api/optical/KN1/search", {}, {**ok, "X-Raven-KN-Virgo-Search": "no"})[0] == 403
        assert post("/api/optical/KN1/search", {"data": "rm -rf"}, ok)[0] == 400
        assert post("/api/optical/unknown/search", {}, ok)[0] == 400
        assert post("/api/optical/KN1/other", {}, ok)[0] == 404
        assert launches(store.root) == []
        code, out = post("/api/optical/KN1/search", {"data": "replay", "statistic": "network"}, ok)
        assert code == 200 and out["launch"]["status"] == "requested" and out["launch"]["data"] == "replay"
        assert post("/api/optical/KN1/search", {}, ok)[0] == 400       # one at a time
        assert [e["source"] for e in read_events(store.root)[0] if e["kind"] == "search_requested"] == ["serve"]
    finally:
        srv.shutdown()
        srv.server_close()
