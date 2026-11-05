"""Treatment module: merger-time window, distance/inclination, detectability, FAR."""

from .t0 import (MergerTimeWindow, estimate_window, window_from_bracket,  # noqa: F401
                 window_from_grb, window_from_samples)
from .distance_inclination import (DistanceInclinationPosterior, combine,  # noqa: F401
                                   constraints_from_candidate, posterior_for)
from .detectability import (DEFAULT_RANGES_MPC, Detectability,  # noqa: F401
                            assess_detectability, network_snr_samples)
from .far import (coinc_far_kn, combine_subsegments, targeted_far,  # noqa: F401
                  targeted_pvalue)
