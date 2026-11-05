"""RAVEN-KN: kilonova-triggered coincidence and targeted-search trigger pipeline.

Built with the same layout as LVK RAVEN (``ligo.raven``) so that it can be
merged into it, or called from gwcelery, later:

* :mod:`raven_kn.connection` -- event classes and connectors (brokers, GraceDB,
  GCN, local files), the analogue of ``ligo.raven.gracedb_events``.
* :mod:`raven_kn.treatment` -- merger-time window, distance/inclination
  combination, GW detectability ("is a search worth it?") and joint FARs,
  the analogue of the FAR part of ``ligo.raven.search``.
* :mod:`raven_kn.skymap` -- kilonova 2D/3D sky maps, overlap integrals with GW
  sky maps (spatial and distance).
* :mod:`raven_kn.search` -- RAVEN-style orchestration: ``query``, ``search``,
  ``coinc_far``, ``calc_signif_gracedb`` plus ``assess`` (targeted-search
  trigger decision).
"""

__version__ = "0.1.0"
