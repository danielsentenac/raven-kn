"""Run NMMA's ``lightcurve-analysis`` with one upstream bug worked around.

Executed by the *NMMA environment's* Python (``python nmma_launcher.py
<lightcurve-analysis arguments>``), so it must not import raven_kn.

Bug (NMMA 1.0.1 / 1.0.2, ``nmma.em.utils.check_model_time_consistency``):
the check takes ``.min()`` of the detection times of *every* filter, and
crashes with "zero-size array to reduction operation minimum" when a filter
has only non-detections. Such filters matter here: the last image before
the merger is often in a filter that has no later detection. The
likelihood itself handles them correctly (upper limits only), so the
patch runs the check on the filters with detections and returns the data
unchanged.
"""
import sys

import numpy as np


def _patch():
    import nmma.em.utils as u
    original = u.check_model_time_consistency

    def check_model_time_consistency(light_curve_data, light_curve_model, priors,
                                     injection=None):
        times, mags, errs, trigger = light_curve_data
        keep = [k for k in times if (np.isfinite(mags[k]) & np.isfinite(errs[k])).any()]
        if injection is not None or len(keep) == len(times):
            return original(light_curve_data, light_curve_model, priors, injection)
        original(({k: times[k] for k in keep}, {k: mags[k] for k in keep},
                  {k: errs[k] for k in keep}, trigger), light_curve_model, priors, injection)
        return light_curve_data

    u.check_model_time_consistency = check_model_time_consistency


if __name__ == "__main__":
    _patch()
    from nmma.em.analysis import main
    sys.argv[0] = "lightcurve-analysis"
    sys.exit(main())
