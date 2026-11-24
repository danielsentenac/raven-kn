"""Announce that an optical transient is ready for a targeted search.

Nothing in the pipeline starts a search: a person does, from the logbook
application. When the search decision of a transient becomes positive
(confirmed as a kilonova and within reach of a GW search,
:func:`raven_kn.search.assess`), this module tells them:

* an ``announcement`` entry in the journal (level ``warning``), which the
  application shows as a banner and at the top of its candidates;
* a mail, when the environment variable ``RAVEN_KN_MAIL_TO`` gives an address
  (several separated by commas), sent with the local ``sendmail``
  (``RAVEN_KN_SENDMAIL``, default ``/usr/sbin/sendmail``). The address is a
  setting of the machine: it is not written in the store or in the journal.

A transient is announced once when it becomes ready, and once more if the
decision is later withdrawn (new photometry that no longer looks like a
kilonova, a revised distance); a decision that stays as it is announces
nothing. What was announced last is kept in ``<store>/announce/<name>.json``.
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
from email.message import EmailMessage
from pathlib import Path
from typing import Any, Callable, Optional

from .journal import _safe, _write_atomic, iso, journal_of, utcnow

log = logging.getLogger("raven_kn.notify")


def _path(root: str | Path, name: str) -> Path:
    return Path(root) / "announce" / f"{_safe(name)}.json"


def announced(root: str | Path, name: str) -> Optional[dict]:
    """What was announced last about a transient: ``{"state": "ready" |
    "withdrawn", "time", "history": [...]}``, or None."""
    try:
        return json.loads(_path(root, name).read_text())
    except (OSError, ValueError):
        return None


def describe(decision: dict) -> list[str]:
    """Lines of an announcement: why it is ready (or no longer), where, when."""
    ev, det = decision.get("kilonova_evidence") or {}, decision.get("detectability") or {}
    di, win = decision.get("distance_inclination") or {}, decision.get("t0_window") or {}
    d = di.get("luminosity_distance_mpc") or {}
    lines = [f"Position: RA {decision.get('ra'):.5f}, Dec {decision.get('dec'):+.5f} deg"]
    if ev.get("confirmed"):
        lines.append("Kilonova evidence: " + ", ".join(ev.get("by") or ["confirmed"]))
    for r in (ev.get("reasons") or [])[:4]:
        lines.append(f"  - {r}")
    if d and di.get("distance_informed", True):
        lines.append(f"Distance: {d['median']:.0f} Mpc ({d['lower']:.0f} to {d['upper']:.0f})")
    if det:
        lines.append(f"Expected network SNR {det.get('snr_median', float('nan')):.1f}, "
                     f"P(detect) {det.get('p_detect', float('nan')):.2f}")
        for r in (det.get("reasons") or [])[:3]:
            lines.append(f"  - {r}")
    if win:
        lines.append(f"Merger-time window: {win.get('duration_h', 0) * 3600:.0f} s from {win.get('method')}, "
                     f"GPS {win.get('gps_start'):.0f} to {win.get('gps_end'):.0f}")
    return lines


def send_mail(subject: str, body: str, to: Optional[str] = None,
              run: Callable[..., Any] = subprocess.run) -> bool:
    """Send a mail with the local sendmail. False (and a warning in the log)
    when no address is set or the mail could not be handed over."""
    to = to if to is not None else os.environ.get("RAVEN_KN_MAIL_TO", "")
    rcpt = [a.strip() for a in to.split(",") if a.strip()]
    if not rcpt:
        return False
    msg = EmailMessage()
    msg["Subject"], msg["To"] = subject, ", ".join(rcpt)
    msg["From"] = os.environ.get("RAVEN_KN_MAIL_FROM") or f"raven-kn@{os.uname().nodename}"
    msg.set_content(body)
    try:
        r = run([os.environ.get("RAVEN_KN_SENDMAIL", "/usr/sbin/sendmail"), "-t", "-oi"],
                input=msg.as_bytes(), timeout=30, capture_output=True)
        if getattr(r, "returncode", 0) != 0:
            raise RuntimeError(f"sendmail exit {r.returncode}")
        return True
    except Exception as e:
        log.warning("announcement not mailed (%s: %s)", type(e).__name__, e)
        return False


def on_decision(store: Any, decision: dict, mail: Callable[..., bool] = send_mail) -> Optional[dict]:
    """Announce a decision if it changes what was announced: ready for the
    first time (or again after a withdrawal), or withdrawn after having been
    announced. Returns the journal entry written, else None. Never raises:
    an announcement that fails must not stop the listener."""
    try:
        name = decision.get("optical") or decision.get("candidate")
        root = Path(store.root)
        ready = bool(decision.get("trigger_targeted_search"))
        last = announced(root, name)
        if ready == (bool(last) and last["state"] == "ready"):
            return None
        if not ready and last is None:
            return None                              # never announced: nothing to withdraw
        state = "ready" if ready else "withdrawn"
        lines = describe(decision)
        if ready:
            summary = f"{name}: ready for a targeted search (confirmed as a kilonova, within reach)"
            head = (f"The optical transient {name} now passes the kilonova test and is within reach of a\n"
                    f"gravitational-wave search. No search has been started: start one from the logbook\n"
                    f"application (the transient, tab GW search) if you judge it worth it.")
        else:
            summary = f"{name}: no longer ready for a targeted search (announced {last['time']})"
            head = (f"The optical transient {name}, announced as ready for a targeted search on {last['time']},\n"
                    f"no longer passes: see the reasons below.")
        mailed = mail(f"[RAVEN-KN] {summary}", head + "\n\n" + "\n".join(lines) + "\n")
        now = iso(utcnow())
        rec = {"optical": name, "state": state, "time": now, "mailed": bool(mailed),
               "history": [*((last or {}).get("history") or []), {"state": state, "time": now}][-20:]}
        p = _path(root, name)
        p.parent.mkdir(parents=True, exist_ok=True)
        _write_atomic(p, json.dumps(rec, indent=1))
        return journal_of(store).log("announcement", summary, level="warning", optical=name,
                                     state=state, mailed=bool(mailed), details=lines)
    except Exception:
        log.exception("announcement failed")
        return None
