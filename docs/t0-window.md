# Merger-time window

The window [t_l, t_h] is where the merger happened. It sets which GW alerts count as coincident and how
much strain data a targeted search must cover, so its width directly sets the trials factor and the joint
FAR. `raven_kn.treatment.t0` builds it from the best information available, in this order.

## 1. Prompt GRB: seconds

If a GRB is associated with the candidate, the window is RAVEN's [−5, +1] s around the GRB trigger.
Nothing else comes close.

## 2. Photometric bracket

Without a GRB, the simplest window runs from the **last constraining non-detection** to the **first
detection**:

- a kilonova is assumed at least as bright as M = −15 from 0.5 d after the merger;
- a non-detection constrains the merger only if it is deeper than that at the 95th-percentile distance;
- a non-detection in another band less than 0.25 d before the first detection is ignored (it reflects
  colour, not age);
- without a constraining non-detection, the early edge is 7 days before the first detection.

## 3. Light-curve fit with a model-systematic margin

A kilonova light-curve fit with [NMMA](https://github.com/nuclear-multimessenger-astronomy/nmma)
(`raven-kn nmma`) gives a posterior on the merger time. RAVEN-KN writes NMMA's data and prior files from the
candidate (non-detections included; distance and inclination priors from the EM measurements, clipped to
the model's training range) and reads the posterior back.

Fits on real kilonovae showed the statistical 90 % interval **misses the merger by about a day** when no
recent non-detection anchors it (see [the examples](examples/kilonova-fits.md)). The models do not follow
real light curves well enough. So the interval is widened by a margin:

$$
\text{margin} = f \times (\text{gap from the median } t_0 \text{ to the data that anchor it}),
\qquad f = 0.75 \text{ by default}
$$

with at least 1 h, and the result is intersected with the photometric bracket. f = 0.75 is the smallest
round value that recovers the merger in every real-data fit; on synthetic kilonovae it covers 97–100 % of
fits, depending on the survey ([injection calibration](examples/injections.md)).

Only the first 3 days of data are fitted by default (`fit_days`), because the models fit the early light
curve better.

## What makes the window narrow

A deep non-detection a few hours before the merger. On AT2017gfo, adding synthetic g and r limits 3 h
before the merger gives a 19 h window and removes the bias (−1.3 h). Survey depth alone does not help:
Rubin's 3-day revisit gives wider windows (83 h median) than ZTF's nightly one (35 h median),
and the combination gives 50 h.
