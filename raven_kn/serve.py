"""Read-only HTTP service on a store, for the logbook application.

Read-only but for one request: ``POST /api/optical/<name>/search`` records
that a targeted search is asked for (:mod:`raven_kn.launch`); the service
itself runs nothing. Standard library only. Every request needs ``Authorization: Bearer <token>``;
without it the answer is a JSON 401, so that a client can tell this service
from a login page put in front of it. A second token (``lvk_token``) also
shows the runs and journal entries marked ``visibility: lvk``.

A third token, the *Virgo targeted-search token*, is not a reading token: it
is sent beside the reading token in the header ``X-Raven-KN-Virgo-Search`` and
will be required to start a targeted search on the Virgo site. ``/api/health``
says whether the service has one and whether the request carried it.

====================================  =========================================
``GET /api/health``                   listeners, store counts, age of the run sync
``GET /api/events?cursor=&limit=``    journal entries after the cursor, oldest first
``GET /api/optical``                  optical transients with their steps
``GET /api/optical/<name>``           optical transient, decision, coincidences, runs, steps
``GET /api/runs/<optical>/<run>``     state, results, SNR series, PE summaries and samples
``GET /api/gw``                       superevents and their sky-map versions
``GET /api/gw/<sid>/map/<n>/<kind>``  reduced map (gzip, one byte per pixel, nside 128 NESTED)
``GET /api/followups``                alerts grouped by the links the pipeline made, each named
                                      by its detection type; which groups are candidates and why
``GET /api/sky``                      positions and error radii for the sky view, and which
                                      alerts each optical transient is attached to
``GET /api/circulars``                GCN Circulars linked to the alerts of the store
``GET /api/circulars/<number>``       one of them, with its text
``GET /api/reported?keys=``           what the circulars linked to the alerts ``<type>:<id>,...``
                                      report, by category: prompt emission, afterglow,
                                      spectroscopy, photometry, links (not verified)
====================================  =========================================

Step status: ``pending``, ``running``, ``done``, ``held``, ``warning``,
``failed``, ``stale``.
"""

from __future__ import annotations

import hmac
import json
import logging
import re
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Optional
from urllib.parse import parse_qs, unquote, urlparse

from .circulars import summary as circular_summary
from .journal import _safe, iso, parse_iso, read_events, sky_versions, utcnow
from .listener import optical_files

log = logging.getLogger("raven_kn.serve")
API_VERSION = 1
LISTENER_ALIVE_S = 180.0         # heartbeat every 30 s
RUN_STALE_S = 600.0              # heartbeat every 60 s

_cache: dict[str, tuple[float, Any]] = {}


def _load(path: Path) -> Any:
    """JSON file, cached on its modification time; None if absent or unreadable."""
    try:
        m = path.stat().st_mtime
    except OSError:
        return None
    hit = _cache.get(str(path))
    if hit and hit[0] == m:
        return hit[1]
    try:
        v = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    _cache[str(path)] = (m, v)
    return v


def _age(t: Optional[str], now=None) -> Optional[float]:
    if not t:
        return None
    try:
        return ((now or utcnow()) - parse_iso(t)).total_seconds()
    except ValueError:
        return None


def _mjd_iso(mjd: float) -> str:
    from datetime import datetime, timedelta, timezone
    return iso(datetime(1858, 11, 17, tzinfo=timezone.utc) + timedelta(days=float(mjd)))


def _gps_iso(gps: Optional[float]) -> Optional[str]:
    if gps is None:
        return None
    from astropy.time import Time
    return Time(float(gps), format="gps").utc.isot + "Z"


class Store:
    """Read-only view of a store directory."""

    def __init__(self, root: str | Path):
        self.root = Path(root)

    # ---- runs -----------------------------------------------------------
    def sync_time(self) -> Optional[str]:
        """Time of the last copy of the runs from the cluster (``runs/.last_sync``,
        touched by the sync script); None when runs are written in place."""
        p = self.root / "runs" / ".last_sync"
        try:
            from datetime import datetime, timezone
            return iso(datetime.fromtimestamp(p.stat().st_mtime, timezone.utc))
        except OSError:
            return None

    def condor(self) -> dict[str, dict]:
        """``runs/condor.json`` (``condor_q -json``, written by the sync script),
        by the name of the job's working directory."""
        out = {}
        for j in _load(self.root / "runs" / "condor.json") or []:
            iwd = str(j.get("Iwd", "")).rstrip("/")
            if iwd:
                out[Path(iwd).name] = j
        return out

    def run_files(self, optical: Optional[str] = None) -> list[Path]:
        """State files of the runs (``<run>_state.json``; several runs may
        share a directory)."""
        pat = f"{_safe(optical)}/*/*_state.json" if optical else "*/*/*_state.json"
        return sorted((self.root / "runs").glob(pat))

    def run_state(self, f: Path, allow_lvk: bool) -> Optional[dict]:
        st = _load(f)
        if st is None:
            return None
        if st.get("visibility", "lvk") == "lvk" and not allow_lvk:
            return {"hidden": True}
        st = dict(st)
        st["run"] = f.name[:-len("_state.json")]
        d = f.parent
        sync = self.sync_time()
        ref = parse_iso(sync) if sync else utcnow()
        st["heartbeat_age_s"] = _age(st.get("heartbeat"), ref)
        # the queue is matched by working directory: only when the directory holds one run
        cq = self.condor().get(d.name) if len(list(d.glob("*_state.json"))) == 1 else None
        st["condor"] = ({"status": {1: "idle", 2: "running", 5: "held"}.get(cq.get("JobStatus"),
                                                                             str(cq.get("JobStatus"))),
                         "hold_reason": cq.get("HoldReason"), "cluster": cq.get("ClusterId")}
                        if cq else None)
        status, detail = "pending", ""
        stage = st.get("stage")
        if st.get("exit") == "failed" or stage == "failed":
            status, detail = "failed", (st.get("error") or "failed")[-300:]
        elif st.get("exit") == "ok":
            status = "done"
            if st.get("p_value") is not None:
                detail = f"p-value {st['p_value']:.3g}"
        elif st["condor"] and st["condor"]["status"] == "held":
            status, detail = "held", f"job held: {st['condor']['hold_reason'] or 'no reason given'}"
        elif stage == "submitted":
            status = "pending"
            detail = "idle in the queue" if st["condor"] and st["condor"]["status"] == "idle" \
                else "submitted, not started"
        elif st["heartbeat_age_s"] is not None and st["heartbeat_age_s"] > RUN_STALE_S:
            status = "stale"
            detail = f"no heartbeat for {st['heartbeat_age_s'] / 60:.0f} min at the last look"
        else:
            status = "running"
            if st.get("subsegments_total"):
                detail = f"sub-segment {st.get('subsegments_done', 0)} of {st['subsegments_total']}"
            elif stage == "pe":
                detail = "parameter estimation"
            else:
                detail = stage or ""
        st["status"], st["detail"] = status, detail
        return st

    def runs(self, optical: str, allow_lvk: bool) -> tuple[list[dict], int]:
        out, hidden = [], 0
        for f in self.run_files(optical):
            st = self.run_state(f, allow_lvk)
            if st is None:
                continue
            if st.get("hidden"):
                hidden += 1
            else:
                out.append(st)
        out.sort(key=lambda s: s.get("started") or s.get("submitted") or "")
        return out, hidden

    def run_detail(self, optical: str, run: str, allow_lvk: bool) -> Optional[dict]:
        f = next(iter(sorted((self.root / "runs" / _safe(optical)).glob(f"*/{_safe(run)}_state.json"))), None)
        st = self.run_state(f, allow_lvk) if f is not None else None
        if st is None or st.get("hidden"):
            return None
        d, label = f.parent, st["run"]
        snr = _load(d / f"{label}_snr_timeseries.json")
        if snr is None and len(list(d.glob("*_state.json"))) == 1:
            snr = _load(d / "snr_timeseries.json")
        out = {"state": st, "results": [], "snr": snr, "pe": []}
        r = _load(d / f"{label}_targeted_result.json")
        if r:
            out["results"].append({"file": f"{label}_targeted_result.json", **r})
        for p in sorted(d.glob(f"{label}_pe_*_pe_summary.json")):
            summ = _load(p)
            if summ:
                out["pe"].append({"file": p.name, "summary": summ,
                                  "samples": _load(d / p.name.replace("_pe_summary", "_pe_samples"))})
        return out

    # ---- optical transients ---------------------------------------------
    def optical_names(self) -> list[str]:
        return [p.stem for p in optical_files(self.root)]

    def optical(self, name: str) -> Optional[dict]:
        return _load(next((p for p in optical_files(self.root) if p.stem == name), self.root / "optical" / f"{name}.json"))

    def decision(self, name: str) -> Optional[dict]:
        return _load(self.root / "assess" / f"{name.replace(' ', '_')}_decision.json")

    def request(self, name: str) -> Optional[dict]:
        return _load(self.root / "assess" / f"{name.replace(' ', '_')}_targeted_search_request.json")

    def coincidences(self, name: str) -> list[dict]:
        out = []
        for p in sorted((self.root / "coinc").glob(f"*__{name}.json")):
            r = _load(p)
            if r:
                sid = r.get("superevent_id", p.stem.split("__")[0])
                r = dict(r)
                r["retracted"] = (self.root / "gw_versions" / _safe(sid) / "retracted").exists()
                r["alert_in_store"] = (self.root / "gw" / f"{sid}.json").exists()
                out.append(r)
        return out

    def steps(self, name: str, allow_lvk: bool, cand: Optional[dict] = None) -> list[dict]:
        c = cand or self.optical(name) or {}
        notes = c.get("notes") or ""
        steps = []

        def add(key, label, status, detail=""):
            steps.append({"key": key, "label": label, "status": status, "detail": detail})

        n = len(c.get("photometry") or [])
        add("optical", "Optical transient received", "done",
            f"{n} photometry point{'s' if n != 1 else ''}" + (f", from {c['source']}" if c.get("source") else ""))
        z = c.get("redshift")
        if z:
            add("host", "Host and redshift", "done",
                f"z = {z.get('z'):.4g} ({z.get('kind', '')}; {z.get('source', '')})")
        elif c.get("distance"):
            add("host", "Host and redshift", "done", "redshift-independent distance")
        elif re.search(r"host lookup: failed", notes):
            add("host", "Host and redshift", "failed", re.search(r"host lookup: (failed[^.]*)", notes).group(1))
        elif "host lookup: no host redshift" in notes:
            add("host", "Host and redshift", "held", "host lookup found no redshift")
        else:
            add("host", "Host and redshift", "pending", "no redshift yet")
        assoc = []
        if c.get("grb"):
            assoc.append(f"GRB {c['grb'].get('graceid') or ''} ({c['grb'].get('instrument', '')})")
        if c.get("neutrino"):
            assoc.append(f"neutrino {c['neutrino'].get('event_name') or ''}")
        detached = "detached" in notes or "retracted" in notes
        add("association", "GRB / neutrino", "done" if assoc else ("held" if detached else "pending"),
            ", ".join(assoc) if assoc else ("association removed (see the log)" if detached
                                           else "none attached"))
        d = self.decision(name)
        if d is None:
            add("decision", "Search decision", "pending", "not assessed yet")
        elif d.get("trigger_targeted_search"):
            ev = d.get("kilonova_evidence") or {}
            add("decision", "Search decision", "done", "search triggered"
                + (f" (kilonova: {', '.join(ev.get('by') or [])})" if ev.get("by") else ""))
        else:
            ev, det = d.get("kilonova_evidence") or {}, d.get("detectability") or {}
            why = []
            if not ev.get("confirmed") and ev.get("required", True):
                why.append("not confirmed as a kilonova: " + "; ".join(ev.get("reasons") or [])[:200])
            if not det.get("worth_search"):
                why.append("; ".join(det.get("reasons") or ["not detectable"])[:200])
            add("decision", "Search decision", "held", " | ".join(why) or "no trigger")
        co = self.coincidences(name)
        if co:
            best = min(co, key=lambda r: _far(r) if _far(r) is not None else 1e99)
            retr = [r["superevent_id"] for r in co if r.get("retracted")]
            add("coincidence", "GW alert coincidence", "warning" if retr else "done",
                f"{best['superevent_id']}: joint FAR {_far(best):.3g} Hz" if _far(best) is not None
                else best["superevent_id"])
            if retr:
                steps[-1]["detail"] += f"; retracted: {', '.join(retr)}"
        else:
            add("coincidence", "GW alert coincidence", "pending", "no GW alert in the window")
        req = self.request(name)
        triggered = bool(d and d.get("trigger_targeted_search"))
        if req:
            add("request", "Search request", "done",
                f"window {req['gps_end'] - req['gps_start']:.0f} s, chirp mass "
                f"{req['mchirp_range'][0]:.2f}-{req['mchirp_range'][1]:.2f}")
        else:
            add("request", "Search request", "pending", "written when the search is triggered")
        runs, hidden = self.runs(name, allow_lvk)
        searches = [r for r in runs if r.get("kind") != "pe"]
        hid = f"; {hidden} run(s) on LVK data not shown" if hidden else ""
        if searches:
            r = searches[-1]
            bad = [x["run"] for x in searches[:-1] if x["status"] in ("failed", "stale", "held")]
            add("search", "Targeted search", r["status"], f"{r['run']}: {r['detail']}"
                + (f"; earlier: {', '.join(bad)} did not finish" if bad else "") + hid)
        else:
            add("search", "Targeted search", "held" if triggered else "pending",
                ("triggered; no search launched yet" if triggered else "not triggered") + hid)
        pes = [(r, k, v) for r in runs for k, v in (r.get("pe") or {}).items()]
        if pes:
            warn = [f"{k}: {w}" for _, k, v in pes for w in v.get("warnings") or []]
            running = [k for r, k, v in pes if v.get("stage") != "done"]
            last = pes[-1][0]
            if last["status"] in ("failed", "stale", "held") and running:
                add("pe", "Parameter estimation", last["status"], f"{last['run']}: {last['detail']}")
            elif running:
                add("pe", "Parameter estimation", "running", f"{', '.join(running)} priors running")
            elif warn:
                add("pe", "Parameter estimation", "warning", "; ".join(warn)[:300])
            else:
                add("pe", "Parameter estimation", "done", ", ".join(k for _, k, _ in pes) + " priors")
        else:
            add("pe", "Parameter estimation", "pending", "not run")
        return steps

    def optical_summary(self, name: str, allow_lvk: bool) -> Optional[dict]:
        c = self.optical(name)
        if c is None:
            return None
        det = [p for p in c.get("photometry") or [] if not p.get("upper_limit")]
        first = min(det, key=lambda p: p["mjd"]) if det else None
        d = self.decision(name)
        steps = self.steps(name, allow_lvk, c)
        order = ["failed", "stale", "warning", "held", "running", "pending", "done"]
        worst = min((s["status"] for s in steps[1:]), key=order.index)
        return {"name": c.get("name", name), "ra": c.get("ra"), "dec": c.get("dec"),
                "source": c.get("source"), "n_photometry": len(c.get("photometry") or []),
                "first_detection": _mjd_iso(first["mjd"]) if first else None,
                "grb": c.get("grb"), "neutrino": c.get("neutrino"),
                "trigger": d.get("trigger_targeted_search") if d else None,
                "steps": steps, "attention": worst if worst in order[:4] else None}

    def optical_detail(self, name: str, allow_lvk: bool) -> Optional[dict]:
        c = self.optical(name)
        if c is None:
            return None
        c = dict(c)
        for k in ("lc_fit_samples", "t0_samples_gps"):           # large, not for display
            if isinstance(c.get(k), (list, dict)) and len(json.dumps(c[k])) > 20000:
                c[k] = None
        phot = [{**p, "time": _mjd_iso(p["mjd"])} for p in c.get("photometry") or []]
        c["photometry"] = sorted(phot, key=lambda p: p["mjd"])
        runs, hidden = self.runs(name, allow_lvk)
        d = self.decision(name)
        marks = []
        if c.get("grb"):
            marks.append({"kind": "grb", "label": f"GRB ({c['grb'].get('instrument', '')})",
                          "time": _gps_iso(c["grb"].get("gps"))})
        if c.get("neutrino"):
            marks.append({"kind": "nu", "label": "neutrino", "time": _gps_iso(c["neutrino"].get("gps"))})
        co = self.coincidences(name)
        for r in co:
            if r.get("gw_gpstime"):
                marks.append({"kind": "gw", "label": r["superevent_id"], "time": _gps_iso(r["gw_gpstime"])})
        win = (d or {}).get("t0_window")
        from .circulars import alert_keys, reported as reported_facts, reported_photometry
        from .connection import load_optical
        facts = None
        try:                                    # magnitudes reported in GCN Circulars, not verified
            path = next(p for p in optical_files(self.root) if p.stem == name)
            kn = load_optical(path)
            reported = reported_photometry(self.root, kn)
            facts = reported_facts(self.root, alert_keys(self.root, kn))
        except Exception:
            log.exception("%s: reported photometry not read", name)
            reported = []
        from .launch import launches
        return {"optical": c, "decision": d, "request": self.request(name), "coincidences": co,
                "reported_photometry": reported,
                # what the circulars about it and about its burst report, by category
                "reported": facts,
                # what the alert brokers have at its position (raven_kn.brokers), None: not asked yet
                "brokers": _load(self.root / "brokers" / f"{_safe(name)}.json"),
                # targeted searches asked for from the application, with what became of them
                "launches": launches(self.root, name),
                "runs": runs, "hidden_runs": hidden, "steps": self.steps(name, allow_lvk, c),
                "marks": marks,
                "window": ({"start": _gps_iso(win["gps_start"]), "end": _gps_iso(win["gps_end"]),
                            "method": win.get("method")} if win else None)}

    # ---- circulars ------------------------------------------------------
    def circulars(self, body: bool = False) -> list[dict]:
        """GCN Circulars linked to alerts of the store, oldest first, each with
        the quantities (``facts``) and the web addresses (``web_links``) read
        in its text (:mod:`raven_kn.circular_facts`)."""
        out = []
        for p in (self.root / "circulars").glob("*.json"):
            c = _load(p)
            if c:
                if "facts" not in c:             # stored before the quantities of the text were read
                    from .circular_facts import facts_in, links_in
                    c["facts"] = facts_in(c.get("subject") or "", c.get("body") or "", c.get("bands"))
                    c["web_links"] = links_in(c.get("body") or "")
                out.append(c if body else {k: v for k, v in c.items() if k != "body"})
        return sorted(out, key=lambda c: c.get("id") or 0)

    def reported(self, keys: str) -> dict:
        """What the circulars linked to the alerts ``keys`` report
        (``grb:<id>,optical:<name>,...``; :func:`raven_kn.circulars.reported`)."""
        from .circulars import reported
        pairs = [tuple(k.split(":", 1)) for k in keys.split(",") if ":" in k]
        return reported(self.root, pairs)

    def circular(self, number: str) -> Optional[dict]:
        return _load(self.root / "circulars" / f"{int(number)}.json") if number.isdigit() else None

    def circular_digests(self) -> dict[str, dict]:
        """``<type>:<id>`` of an alert -> digest of its circulars
        (:func:`raven_kn.circulars.summary`)."""
        by: dict[str, list] = {}
        for c in self.circulars():
            for ln in c.get("links") or []:
                by.setdefault(f"{ln['type']}:{ln['id']}", []).append(c)
        return {k: circular_summary(v) for k, v in by.items()}

    # ---- sky ------------------------------------------------------------
    def superevents(self) -> list[dict]:
        ids = {p.name for p in (self.root / "gw_versions").glob("*") if p.is_dir()}
        ids |= {p.stem for p in (self.root / "gw").glob("*.json")}
        grbs = {p.stem: (_load(p) or {}) for p in sorted((self.root / "grb").glob("*.json"))}
        out = []
        for sid in sorted(ids):
            vers = sky_versions(self.root, sid)
            alert = _load(self.root / "gw" / f"{sid}.json")
            # LVK's own GW x GRB coincidence: from the alert, else from the latest version that has one
            ec = (alert or {}).get("external_coinc") or next(
                (v["external_coinc"] for v in reversed(vers) if v.get("external_coinc")), None)
            gps = (alert or {}).get("gpstime") or next(
                (v["gpstime"] for v in reversed(vers) if v.get("gpstime")), None)
            m = match_external_grb(ec, gps, grbs)
            out.append({"superevent_id": sid, "alert": alert,
                        "retracted": (self.root / "gw_versions" / _safe(sid) / "retracted").exists(),
                        "versions": vers, "external_coinc": ec,
                        "grb_id": m[0] if m else None, "grb_match": m[1] if m else None})
        return out

    def sky(self) -> dict:
        cands = []
        for n in self.optical_names():
            c = self.optical(n)
            if c:
                g, v = c.get("grb") or {}, c.get("neutrino") or {}
                det = [p for p in c.get("photometry") or [] if not p.get("upper_limit")]
                cands.append({"name": c.get("name", n), "ra": c.get("ra"), "dec": c.get("dec"),
                              "pos_err_arcsec": c.get("pos_err_arcsec"),
                              "time": _mjd_iso(min(p["mjd"] for p in det)) if det else None,
                              # who found it: the survey of the first detection, else the stated source of the transient
                              "origin": (min(det, key=lambda p: p["mjd"]).get("survey") if det else None)
                              or c.get("source"),
                              "grb": g.get("graceid"), "neutrino": v.get("event_name"),
                              # the alerts this transient is attached to, by their id in this answer
                              "grb_id": _notice_key(g.get("graceid"), g) if g else None,
                              "nu_id": _notice_key(v.get("event_name"), v) if v else None,
                              "superevents": [r["superevent_id"] for r in self.coincidences(c.get("name", n))]})
        grb = []
        for p in sorted((self.root / "grb").glob("*.json")):
            r = _load(p) or {}
            g, x = r.get("grb"), r.get("xray")
            if g and r.get("ra") is not None:
                grb.append({"id": p.stem, "instrument": g.get("instrument"), "trigger": g.get("graceid"),
                            "time": _gps_iso(g.get("gps")), "ra": r["ra"], "dec": r["dec"],
                            "err_deg": r.get("err_deg"), "xray": None})
            if x and x.get("ra") is not None:
                e = {"instrument": x.get("instrument"), "ra": x["ra"], "dec": x["dec"],
                     "err_deg": x.get("err_deg")}
                if g and grb and grb[-1]["id"] == p.stem:
                    grb[-1]["xray"] = e
                else:
                    grb.append({"id": p.stem, "instrument": x.get("instrument"), "trigger": p.stem,
                                "time": None, "ra": x["ra"], "dec": x["dec"],
                                "err_deg": x.get("err_deg"), "xray": e, "burst_missing": True})
        nu = []
        for p in sorted((self.root / "nu").glob("*.json")):
            v = (_load(p) or {}).get("nu") or {}
            if v.get("ra") is not None:
                nu.append({"id": p.stem, "event": v.get("event_name"), "instrument": v.get("instrument"),
                           "time": _gps_iso(v.get("gps")), "ra": v["ra"], "dec": v["dec"],
                           "err_deg": v.get("err_deg"), "p_astro": v.get("p_astro")})
        gw = self.superevents()
        # what the GCN Circulars say about each alert (None: no circular linked to it)
        dig = self.circular_digests()
        for kind, items, key in (("optical", cands, "name"), ("grb", grb, "id"), ("nu", nu, "id"),
                                 ("gw", gw, "superevent_id")):
            for it in items:
                it["circulars"] = dig.get(f"{kind}:{it[key]}")
        return {"optical": cands, "grb": grb, "nu": nu, "gw": gw}

    # ---- follow-ups -----------------------------------------------------
    def followups(self, allow_lvk: bool) -> list[dict]:
        """Alerts grouped by the links the pipeline made, each alert named by
        its detection type (optical, GW, GRB, neutrino), and whether the
        group is a **candidate**.

        Links: a GRB or a neutrino attached to an optical transient; a
        superevent in the merger-time window of an optical transient; the
        GW x GRB coincidence carried by an LVK alert.

        A follow-up is a candidate when any of these holds:

        * two detection types are linked;
        * its optical transient is confirmed as a kilonova by the decision;
        * its GW alert, not retracted, is likely to involve a neutron star:
          P(BNS) + P(NSBH) >= 0.5, or HasRemnant >= 0.5;
        * a GCN Circular gives its GRB a redshift within the range of a GW
          search (:func:`raven_kn.circulars.distance_gate`).

        A GRB whose redshift puts it beyond that range is written in
        ``notes`` and sets ``too_far``. Such a follow-up is not a candidate
        unless it holds a GW alert: a burst at z = 3.7 linked to its own
        optical afterglow is two detection types, and no GW source.
        """
        sky = self.sky()
        nodes: dict[str, dict] = {}
        links: dict[str, set] = {}

        def link(a, b, why):
            if a in nodes and b in nodes:
                links.setdefault(a, set()).add(b)
                links.setdefault(b, set()).add(a)
                reasons.setdefault(frozenset((a, b)), why)
        reasons: dict[frozenset, str] = {}
        circ = self.circular_digests()
        for g in sky["grb"]:
            nodes[f"grb:{g['id']}"] = {"type": "grb", "id": g["id"], "time": g.get("time"),
                                       "instrument": g.get("instrument"),
                                       "name": f"{g.get('instrument') or ''} {g.get('trigger') or ''}".strip()}
        for n in sky["nu"]:
            nodes[f"nu:{n['id']}"] = {"type": "nu", "id": n["id"], "time": n.get("time"),
                                      "instrument": n.get("instrument"),
                                      "name": str(n.get("event") or n["id"])}
        for g in sky["gw"]:
            vers = g.get("versions") or []
            src = g.get("alert") or (vers[-1] if vers else {})
            nodes[f"gw:{g['superevent_id']}"] = {
                "type": "gw", "id": g["superevent_id"], "name": g["superevent_id"],
                "time": _gps_iso(src.get("gpstime")), "retracted": bool(g.get("retracted")),
                "instrument": " ".join(src.get("instruments") or []) or None,
                "classification": src.get("classification"), "properties": src.get("properties"),
                "far": src.get("far")}
        for c in sky["optical"]:
            nodes[f"optical:{c['name']}"] = {"type": "optical", "id": c["name"], "name": c["name"],
                                             "time": c.get("time"), "instrument": c.get("origin")}
        for c in sky["optical"]:
            o = f"optical:{c['name']}"
            if c.get("grb_id"):
                link(o, f"grb:{c['grb_id']}", "GRB attached to the optical transient")
            if c.get("nu_id"):
                link(o, f"nu:{c['nu_id']}", "neutrino attached to the optical transient")
            for sid in c.get("superevents") or []:
                link(o, f"gw:{sid}", "GW alert in the merger-time window of the optical transient")
        for g in sky["gw"]:
            if g.get("grb_id"):
                link(f"gw:{g['superevent_id']}", f"grb:{g['grb_id']}", "GW x GRB coincidence reported by LVK")
        out, seen = [], set()
        rank = {"optical": 0, "gw": 1, "grb": 2, "nu": 3}
        for start in nodes:
            if start in seen:
                continue
            comp, todo = set(), [start]
            while todo:
                x = todo.pop()
                if x in comp:
                    continue
                comp.add(x)
                todo.extend(links.get(x, ()))
            seen |= comp
            members = sorted(({"key": k, **nodes[k], "circulars": circ.get(k)} for k in comp),
                             key=lambda m: (m["time"] is None, m["time"] or "", rank[m["type"]]))
            why = sorted({r for pair, r in reasons.items() if pair <= comp})
            optical = next((m["id"] for m in members if m["type"] == "optical"), None)
            steps, attention = None, None
            if optical:
                d = self.decision(optical) or {}
                ev = d.get("kilonova_evidence") or {}
                if ev.get("confirmed"):
                    why.append("optical transient confirmed as a kilonova (" + ", ".join(ev.get("by") or []) + ")")
                summ = self.optical_summary(optical, allow_lvk) or {}
                steps, attention = summ.get("steps"), summ.get("attention")
            for m in members:
                if m["type"] != "gw" or m.get("retracted"):
                    continue
                cl, pr = m.get("classification") or {}, m.get("properties") or {}
                p_ns = float(cl.get("BNS") or 0) + float(cl.get("NSBH") or 0)
                if p_ns >= 0.5 or float(pr.get("HasRemnant") or 0) >= 0.5:
                    why.append(f"GW alert likely to involve a neutron star (BNS {float(cl.get('BNS') or 0):.2g}, "
                               f"NSBH {float(cl.get('NSBH') or 0):.2g}, HasRemnant {float(pr.get('HasRemnant') or 0):.2g})")
            notes = []
            for m in members:
                gate = (m.get("circulars") or {}).get("distance")
                if m["type"] != "grb" or not gate:
                    continue
                text = f"{m['circulars'].get('name') or 'GRB ' + m['name']}: {gate['text']}"
                (why if gate["within_reach"] else notes).append(text)
            anchor = min(members, key=lambda m: rank[m["type"]])
            # title: the GW name when the follow-up has a live superevent, else its optical, GRB or neutrino name
            lead = next((m for t in ("gw", "optical", "grb", "nu") for m in members
                         if m["type"] == t and not m.get("retracted")), members[0])
            label = {"gw": "GW", "optical": "Optical", "grb": "GRB", "nu": "Neutrino"}
            types = [t for t in ("gw", "grb", "nu", "optical") if any(m["type"] == t for m in members)]
            out.append({"id": anchor["key"], "name": anchor["name"],
                        "title": f"{label[lead['type']]} {lead['name']}", "origin": members[0]["type"],
                        "origin_name": members[0]["name"], "time": members[0]["time"], "types": types,
                        "members": members, "optical": optical, "reasons": why,
                        "candidate": bool(why) and not (notes and "gw" not in types),
                        "notes": notes, "too_far": bool(notes),
                        "steps": steps, "attention": attention})
        out.sort(key=lambda f: f["time"] or "", reverse=True)
        return out

    # ---- health ---------------------------------------------------------
    def health(self, allow_lvk: bool, virgo_search: Optional[dict] = None) -> dict:
        listeners = []
        for p in sorted((self.root / "status").glob("*.json")):
            s = _load(p)
            if s:
                s = dict(s)
                s["age_s"] = _age(s.get("heartbeat"))
                s["alive"] = s["age_s"] is not None and s["age_s"] < LISTENER_ALIVE_S
                listeners.append(s)
        sync = self.sync_time()
        count = lambda d, pat="*.json": sum(1 for _ in (self.root / d).glob(pat))   # noqa: E731
        return {"service": "raven-kn", "api": API_VERSION, "time": iso(utcnow()), "lvk": allow_lvk,
                # who may start a targeted search: on the Virgo site with its own token
                "targeted_search": {"virgo": virgo_search or {"configured": False, "granted": False}},
                "listeners": listeners,
                "runs_sync": {"time": sync, "age_s": _age(sync)} if sync else None,
                "store": {"optical": len(optical_files(self.root)), "gw": count("gw"), "grb": count("grb"),
                          "nu": count("nu"), "coinc": count("coinc"), "circulars": count("circulars"),
                          "runs": len(self.run_files())}}


def match_external_grb(ec: Optional[dict], gw_gps: Optional[float], grbs: dict[str, dict],
                       tolerance_s: float = 2.0) -> Optional[tuple[str, str]]:
    """Stored GRB notice of the GW x GRB coincidence carried by an LVK alert
    (``external_coinc``): (id in the store, how it was matched), or None.

    1. ``gcn_notice_id`` equal to the trigger number of a stored notice, of a
       compatible mission when the alert names its observatory;
    2. else a notice of that observatory whose time differs from the GW time
       by ``time_difference``, within ``tolerance_s`` (either sign: the
       convention was not checked on a real alert).
    ``grbs``: store id -> stored record (``grb``: instrument, gps, graceid)."""
    if not ec:
        return None

    def norm(s):
        return "".join(c for c in str(s or "").lower() if c.isalnum())
    obs = norm(ec.get("observatory"))

    def mission_ok(g):
        ins = norm(g.get("instrument"))
        return not obs or not ins or obs in ins or ins in obs
    nid = ec.get("gcn_notice_id")
    if nid not in (None, ""):
        for key, rec in grbs.items():
            g = rec.get("grb") or {}
            if g.get("graceid") is not None and str(g["graceid"]) == str(nid) and mission_ok(g):
                return key, "notice number"
    td = ec.get("time_difference")
    if td is None or gw_gps is None or not obs:
        return None
    best = None
    for key, rec in grbs.items():
        g = rec.get("grb") or {}
        if g.get("gps") is None or not norm(g.get("instrument")) or not mission_ok(g):
            continue
        off = abs(abs(float(g["gps"]) - float(gw_gps)) - abs(float(td)))
        if off <= tolerance_s and (best is None or off < best[1]):
            best = (key, off)
    return (best[0], "mission and time") if best else None


def _notice_key(ident: Optional[str], rec: dict) -> str:
    """File stem of a stored GRB notice or neutrino alert (as ``AlertStore``)."""
    return _safe(ident or f"{rec.get('instrument', '')}_{float(rec.get('gps', 0)):.2f}")


def _far(r: dict) -> Optional[float]:
    return (r.get("spatiotemporal_distance_coinc_far") or r.get("spatiotemporal_coinc_far")
            or r.get("temporal_coinc_far"))


# --------------------------------------------------------------------------

def make_handler(store: Store, token: Optional[str], lvk_token: Optional[str],
                 virgo_search_token: Optional[str] = None):
    class Handler(BaseHTTPRequestHandler):
        server_version = "raven-kn"
        protocol_version = "HTTP/1.1"

        def log_message(self, fmt, *args):
            log.info("%s %s", self.address_string(), fmt % args)

        def _send(self, code: int, body: bytes, ctype: str = "application/json") -> None:
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Raven-KN", str(API_VERSION))
            self.end_headers()
            self.wfile.write(body)

        def _json(self, obj: Any, code: int = 200) -> None:
            self._send(code, json.dumps(obj, default=str, separators=(",", ":")).encode())

        def _scope(self) -> Optional[bool]:
            """None: refused; False: public runs only; True: also LVK runs."""
            h = self.headers.get("Authorization", "")
            got = h[7:].strip() if h.lower().startswith("bearer ") else ""
            if not got:
                return None
            if lvk_token and hmac.compare_digest(got.encode(), lvk_token.encode()):
                return True
            if token and hmac.compare_digest(got.encode(), token.encode()):
                return False
            return None

        def _virgo_search(self) -> dict:
            """Whether the service has a Virgo targeted-search token, and
            whether this request carries it."""
            got = (self.headers.get("X-Raven-KN-Virgo-Search") or "").strip()
            ok = bool(virgo_search_token and got
                      and hmac.compare_digest(got.encode(), virgo_search_token.encode()))
            return {"configured": bool(virgo_search_token), "granted": ok}

        def do_POST(self):                                   # noqa: N802
            try:
                self._post()
            except BrokenPipeError:
                pass
            except Exception as e:
                log.exception("request failed: %s", self.path)
                self._json({"service": "raven-kn", "error": f"{type(e).__name__}: {e}"}, 500)

        def _post(self):
            """``POST /api/optical/<name>/search`` with ``{"data": site | replay,
            "statistic": coherent | network}``: ask for a targeted search. Needs
            a reading token and the Virgo targeted-search token."""
            u = urlparse(self.path)
            p = [unquote(x) for x in u.path.strip("/").split("/")]
            n = int(self.headers.get("Content-Length") or 0)
            body = self.rfile.read(min(n, 4096)) if n > 0 else b""
            if self._scope() is None:
                return self._json({"service": "raven-kn", "api": API_VERSION, "error": "unauthorized"}, 401)
            if len(p) != 4 or p[0] != "api" or p[1] != "optical" or p[3] != "search" \
                    or any(x in ("", ".", "..") or "/" in x or "\\" in x for x in p):
                return self._json({"service": "raven-kn", "error": "not found"}, 404)
            if not self._virgo_search()["granted"]:
                return self._json({"service": "raven-kn", "error": "the Virgo targeted-search token is "
                                   "required to start a search"}, 403)
            try:
                opts = json.loads(body or b"{}")
                if not isinstance(opts, dict) or n > 4096:
                    raise ValueError("the request must be a small JSON object")
                from .launch import request_search
                from .journal import Journal
                rec = request_search(store, p[2], str(opts.get("data", "site")),
                                     str(opts.get("statistic", "coherent")),
                                     journal=Journal.for_store(store.root, "serve"))
            except ValueError as e:
                return self._json({"service": "raven-kn", "error": str(e)}, 400)
            return self._json({"launch": rec})

        def do_GET(self):                                    # noqa: N802
            try:
                self._get()
            except BrokenPipeError:
                pass
            except Exception as e:                           # never a traceback page
                log.exception("request failed: %s", self.path)
                self._json({"service": "raven-kn", "error": f"{type(e).__name__}: {e}"}, 500)

        def _get(self):
            u = urlparse(self.path)
            parts = [unquote(x) for x in u.path.strip("/").split("/")]
            q = {k: v[0] for k, v in parse_qs(u.query).items()}
            lvk = self._scope()
            if lvk is None:
                return self._json({"service": "raven-kn", "api": API_VERSION,
                                   "error": "unauthorized"}, 401)
            if parts[:1] != ["api"]:
                return self._json({"service": "raven-kn", "error": "not found"}, 404)
            p = parts[1:]
            if any(x in ("", ".", "..") or "/" in x or "\\" in x for x in p):
                return self._json({"service": "raven-kn", "error": "not found"}, 404)
            if p == ["health"]:
                return self._json(store.health(lvk, self._virgo_search()))
            if p == ["events"]:
                limit = max(1, min(int(q.get("limit", 2000)), 5000))
                ev, cur, more = read_events(store.root, q.get("cursor"), limit, allow_lvk=lvk)
                return self._json({"events": ev, "cursor": cur, "more": more})
            if p == ["optical"]:
                out = [store.optical_summary(n, lvk) for n in store.optical_names()]
                return self._json({"optical": [c for c in out if c]})
            if len(p) == 2 and p[0] == "optical":
                d = store.optical_detail(p[1], lvk)
                return self._json(d) if d else self._json({"service": "raven-kn", "error": "not found"}, 404)
            if len(p) == 3 and p[0] == "runs":
                d = store.run_detail(p[1], p[2], lvk)
                return self._json(d) if d else self._json({"service": "raven-kn", "error": "not found"}, 404)
            if p == ["gw"]:
                return self._json({"superevents": store.superevents()})
            if len(p) == 5 and p[0] == "gw" and p[2] == "map" and p[3].isdigit() \
                    and p[4] in ("gw", "combined"):
                f = store.root / "gw_versions" / _safe(p[1]) / f"{int(p[3]):02d}_{p[4]}.cred.gz"
                if f.is_file():
                    return self._send(200, f.read_bytes(), "application/octet-stream")
                return self._json({"service": "raven-kn", "error": "not found"}, 404)
            if p == ["sky"]:
                return self._json(store.sky())
            if p == ["followups"]:
                return self._json({"followups": store.followups(lvk)})
            if p == ["circulars"]:
                return self._json({"circulars": store.circulars()})
            if p == ["reported"]:
                return self._json(store.reported(q.get("keys", "")))
            if len(p) == 2 and p[0] == "circulars":
                d = store.circular(p[1])
                return self._json(d) if d else self._json({"service": "raven-kn", "error": "not found"}, 404)
            return self._json({"service": "raven-kn", "error": "not found"}, 404)

    return Handler


def make_server(root: str | Path, host: str = "127.0.0.1", port: int = 8765,
                token: Optional[str] = None, lvk_token: Optional[str] = None,
                virgo_search_token: Optional[str] = None) -> ThreadingHTTPServer:
    if not token and not lvk_token:
        raise ValueError("an access token is required")
    if virgo_search_token and virgo_search_token in (token, lvk_token):
        raise ValueError("the Virgo targeted-search token must differ from the reading tokens")
    srv = ThreadingHTTPServer((host, port), make_handler(Store(root), token, lvk_token,
                                                         virgo_search_token))
    srv.daemon_threads = True
    return srv


def serve(root: str | Path, host: str = "127.0.0.1", port: int = 8765,
          token: Optional[str] = None, lvk_token: Optional[str] = None,
          virgo_search_token: Optional[str] = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    srv = make_server(root, host, port, token, lvk_token, virgo_search_token)
    log.info("serving %s on http://%s:%d (LVK runs: %s; Virgo targeted-search token: %s)", root, host,
             port, "with the LVK token" if lvk_token else "never shown",
             "set" if virgo_search_token else "not set")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()
