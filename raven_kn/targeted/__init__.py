"""Targeted GW search at a kilonova's position and merger-time window.

:mod:`raven_kn.targeted.search` is a light, self-contained PyCBC-based
search (matched filter with a reduced bank, detector SNRs combined with the
fixed time delays of the known sky position, time-slide background, p-value),
following the approach shown by Jamsin et al. (arXiv:2610.28775). It is
meant for quick assessment and validation; the production LVK targeted
search is PyGRB (``pycbc_make_offline_grb_workflow``), with a coherent
statistic and signal-consistency tests.
"""
