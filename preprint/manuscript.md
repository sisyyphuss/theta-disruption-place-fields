# Cumulative Theta-Rhythm Disruption Predicts a Reduction in Hippocampal Place-Field Count: An Exploratory Reanalysis of Public Medial Septal Cooling Data

**Saad Yusuf**  
*Department of Computer Science, Jamia Millia Islamia, New Delhi, India*  
Correspondence: `syed.saad7527@gmail.com` | ORCID: `0009-0000-2333-0078`  

**Funding:** This research received no specific funding.  
**Competing interests:** The author declares no competing interests.

---

## Abstract

The medial septum paces hippocampal theta oscillations, and transient cooling of the medial septum acutely slows theta and impairs spatial-alternation performance without disrupting the spatial firing correlates of individual place cells (Petersen & Buzsáki, 2020). That result was established within single recording sessions, comparing brief cooling epochs to immediately preceding baseline trials. Whether *repeated* exposure to this transient disruption, accumulated across the many-week span of such an experiment, leaves any trace on the hippocampal spatial code has not been tested. We reanalyzed the publicly archived dataset underlying that study (Dandiset 000059, DANDI Archive), reconstructing baseline place-field statistics for every well-isolated unit in every available recording session from four subjects (15 usable sessions of 50 attempted, spanning up to six weeks per animal), using an occupancy-normalized rate-map pipeline built on an openly available teaching notebook. We regressed three session-level place-field summary statistics — Skaggs spatial information, place-field count, and peak firing rate — against each animal's cumulative prior cooling exposure, using permutation-based inference to accommodate the small sample size, and cross-checked results against ordinary least squares to guard against sample-size artifacts. We found that cumulative cooling exposure is closely confounded with elapsed experimental time in this dataset's design (within-subject r = 0.90–0.98), a structural limitation that prevents fully attributing any observed trend to the cooling manipulation specifically. With that caveat stated explicitly, place-field count declined and peak firing rate increased with cumulative cooling exposure (permutation p = 0.022 and p = 0.044, respectively; parametric OLS in close agreement), while overall spatial information showed only a non-significant trend (p = 0.083). These results are preliminary given the sample size, but they are internally consistent, survive a documented artifact-correction step, and motivate a properly powered follow-up study using interleaved cooling and sham sessions to dissociate chronic theta disruption from the passage of time. All code and cached intermediate data are openly released.

---

## Introduction

Place cells in hippocampal area CA1 fire selectively at specific locations in an animal's environment, forming a population code for space (O'Keefe & Dostrovsky, 1971). The medial septum is the principal pacemaker of the hippocampal theta rhythm (6–10 Hz), and theta-phase timing is thought to coordinate the activity of hippocampal cell assemblies during navigation. Petersen and Buzsáki (2020) used reversible cooling of the medial septum in freely behaving rats performing a figure-8 spatial-alternation task to test this directly: cooling slowed and weakened theta, increased behavioral choice errors, and reduced distance-time compression of theta sequences — but left the spatial firing correlates of individual place cells intact, including a population-vector correlation between baseline and cooling epochs computed *within* each session. Their conclusion was that theta-phase coordination organizes hippocampal assemblies independently of, and is dispensable for, immediate single-trial spatial coding.

That finding characterizes the *acute* consequence of theta disruption. It leaves open a distinct question: does *repeated* transient disruption of theta timing, accumulated over the many sessions and weeks typical of such an experiment, leave any cumulative trace on the stability or quality of the hippocampal spatial map — even if each individual disruption is, by itself, undetectable at the single-session level? This is a natural extension of the original result's own logic: if theta-phase coordination is a mechanism (rather than a mere byproduct) of assembly organization, chronic interference with it could plausibly degrade map quality on a timescale invisible to a single 40-minute recording session.

We addressed this question by reanalyzing the publicly archived Neurodata Without Borders (NWB) dataset underlying Petersen and Buzsáki (2020), streamed directly from the DANDI Archive (Dandiset 000059). We computed standard place-field summary statistics for every usable session from every subject in the dataset, restricted to pre-cooling baseline trials, and tested whether these statistics covary with each animal's cumulative cooling exposure accumulated up to that point in the experiment.

## Methods

### Data source

We used Dandiset 000059 (version 0.230907.2101), the NWB-formatted dataset accompanying Petersen and Buzsáki (2020), comprising extracellular silicon-probe recordings from dorsal CA1 in five male Long-Evans rats (subjects MS10, MS12, MS13, MS21, MS22) performing a figure-8 spatial-alternation task, with reversible medial septal cooling applied within most recording sessions. Data were streamed on demand from the DANDI Archive (no local download) using `pynwb`, `remfile`, and the `dandi` Python client.

### Code provenance

The data-streaming pattern and the core occupancy-normalized place-field computation (spike-to-position alignment via `numpy.searchsorted`, two-dimensional spatial binning via `numpy.histogram2d`, and occupancy normalization) were adapted directly from an openly available teaching notebook series for NWB-based hippocampal data analysis (Denovellis & Joshi, *Neural Circuits for Navigation* course materials, NCBS; https://github.com/edeno/ncbs-neural-circuits-navigation), generalized from a single hard-coded session to iterate over every available session per subject. All subsequent code — cooling-epoch extraction, spatial-information and field-count computation, session-level aggregation, and statistical modeling — was written for this study and is not derived from that source (see *Code availability*).

### Unit and session selection

For each session, units were restricted to those with a curated quality label of `good` or `good2` (well-isolated single units), matching the original study's convention. A session was excluded if it lacked a `units` table, lacked a `trials` table, or lacked a `cooling state` trial-annotation column entirely (a subset of sessions in this dataset were not run with the cooling manipulation, or predate/postdate the annotated protocol). Within an included session, a unit was excluded if it fired fewer than 30 spikes during the pre-cooling baseline window.

### Baseline window and cumulative cooling dose

Each session's trials table includes a `cooling state` column with values `Pre-Cooling`, `Cooling on`, `Cooling off`, and `Post-Cooling`. We restricted all place-field computation to the contiguous span of `Pre-Cooling` trials at the start of each session, so that reported place-field statistics reflect baseline spatial coding uncontaminated by the acute effect of cooling itself. This session's own cooling dose was defined as the summed duration of its `Cooling on`/`Cooling off` trials. Each session's **cumulative prior cooling exposure** was defined as the running sum of cooling-trial duration across all of that animal's earlier included sessions (i.e., dose accumulated strictly *before* the session being scored), and **elapsed days** as calendar days since that animal's first included session.

### Place-field computation

For each unit, spikes were restricted to the pre-cooling baseline window and matched to the animal's position at spike time via nearest-preceding-sample lookup (`numpy.searchsorted`). Position and spike locations were binned into a 2 cm square grid spanning the recorded arena. A firing-rate map was computed as spike count divided by occupancy time in each bin. **Bins with less than 0.1 s of total occupancy were excluded** (treated as unvisited) from the rate map, occupancy weighting, and all downstream statistics; this threshold follows convention in the place-field literature (Muller & Kubie, 1987) and was necessary here because this dataset's ~120 Hz position sampling rate means a single-sample bin (~8 ms) coinciding with even one incidental spike produces a formally correct but spurious ~120 Hz "rate" — an artifact we identified during development (see *Data and code availability* for the affected/corrected commit) and corrected before final analysis.

Three session-level summary statistics were computed as the median across that session's included units:

- **Spatial information** (bits/spike), following Skaggs, McNaughton, Gothard, and Markus (1993): $I = \sum_i p_i \frac{r_i}{\bar{r}} \log_2\frac{r_i}{\bar{r}}$, where $p_i$ is occupancy probability and $r_i$ is firing rate in bin $i$.
- **Place-field count**: the rate map was smoothed with a Gaussian kernel (σ = 1.5 bins) and the number of spatially contiguous regions (4-connected) exceeding 20% of the peak rate and spanning at least 3 bins was counted.
- **Peak firing rate** (Hz): the maximum of the (unsmoothed, occupancy-masked) rate map.

### Statistical analysis

Session-level statistics from all included sessions across all four contributing subjects were pooled. We modeled each outcome as a linear function of a single time-related predictor (cumulative prior cooling exposure, or separately, elapsed days) plus a subject fixed effect (`outcome ~ predictor + C(subject)`), fit by ordinary least squares.

An initial model including *both* time-related predictors simultaneously produced results in stark disagreement between permutation and parametric inference (permutation p ≈ 0 against a parametric p = 0.21 for the same coefficient), traced to severe multicollinearity between the two predictors (variance inflation reflected in an OLS condition number of 2410) — cumulative cooling dose and elapsed days are correlated within subject at r = 0.90–0.98, because cooling was applied in nearly every session and therefore accumulates almost linearly with session order. We do not report that joint model as a result; we flag its instability here because a smaller apparent effect size that agrees between two independent inferential methods is more trustworthy than a larger one that does not. All results reported below use the single-predictor specification, for which permutation and parametric p-values converge closely, consistent with the joint model's instability being an artifact of collinearity rather than any single-predictor result being spurious.

Statistical significance was assessed via a permutation test rather than relying solely on OLS asymptotics, given the small sample (n = 15 sessions). For each of 5,000 permutations, the time-related predictor was shuffled *within* each subject's own sessions (preserving each animal's session count, dose range, and subject fixed effect structure while breaking the link between dose/time and the outcome), the model was refit, and the resulting coefficient recorded. The two-sided permutation p-value is the fraction of permuted coefficients with absolute value at least as large as the observed coefficient.

### Software

Python 3.13; `pynwb`, `dandi`, `remfile`, `h5py` for data access; `numpy`, `scipy` for computation; `pandas`, `statsmodels` for tabulation and regression; `matplotlib` for figures.

## Results

### Data yield

Of 50 sessions attempted across five subjects, 15 (30%) met inclusion criteria (Figure 1; Table 1). Subject MS10 contributed no usable sessions (all 9 attempted sessions lacked a `units` table in this NWB release) and was excluded from all downstream analysis. The remaining four subjects contributed between 2 and 6 usable sessions each.

**Table 1. Session yield by subject.**

| Subject | Sessions attempted | Sessions usable | Primary exclusion reason(s) |
|---|---|---|---|
| MS10 | 9 | 0 | No `units` table in any session |
| MS12 | 8 | 2 | No `units` table (5 sessions); no unit passed spike-count threshold (1 session) |
| MS13 | 5 | 2 | No `units` table (3 sessions) |
| MS21 | 15 | 5 | No `units` table (7); no `cooling state` column (1); no `trials` table (1); no unit passed threshold (1) |
| MS22 | 13 | 6 | No `units` table (7) |
| **Total** | **50** | **15** | |

**Table 2. Session-level dataset (n = 15).**

| Subject | Date | Elapsed days | Cumulative prior cooling (s) | Spatial info (bits/spike) | Field count | Peak rate (Hz) | Units used |
|---|---|---|---|---|---|---|---|
| MS12 | 2017-07-16 | 0 | 0.0 | 1.934 | 4.0 | 52.2 | 34 |
| MS12 | 2017-07-19 | 3 | 280.8 | 1.702 | 4.5 | 41.3 | 34 |
| MS13 | 2017-11-28 | 0 | 0.0 | 2.192 | 3.0 | 24.0 | 51 |
| MS13 | 2017-12-01 | 3 | 177.1 | 2.064 | 3.5 | 28.7 | 40 |
| MS21 | 2018-06-21 | 0 | 0.0 | 2.853 | 7.5 | 11.9 | 16 |
| MS21 | 2018-06-25 | 4 | 119.9 | 1.730 | 3.0 | 32.0 | 70 |
| MS21 | 2018-06-28 | 7 | 315.8 | 1.618 | 3.0 | 45.6 | 100 |
| MS21 | 2018-06-29 | 8 | 459.8 | 2.381 | 3.0 | 30.0 | 141 |
| MS21 | 2018-07-19 | 28 | 671.9 | 1.746 | 3.0 | 42.5 | 122 |
| MS22 | 2018-06-21 | 0 | 0.0 | 1.809 | 9.0 | 12.6 | 27 |
| MS22 | 2018-06-28 | 7 | 151.6 | 1.849 | 5.0 | 22.9 | 57 |
| MS22 | 2018-06-29 | 8 | 342.4 | 1.651 | 3.0 | 38.4 | 127 |
| MS22 | 2018-07-11 | 20 | 520.8 | 1.277 | 3.0 | 42.0 | 62 |
| MS22 | 2018-07-19 | 28 | 668.5 | 1.701 | 3.0 | 26.5 | 58 |
| MS22 | 2018-07-20 | 29 | 788.8 | 1.350 | 4.0 | 39.9 | 76 |

### Cumulative cooling exposure predicts place-field count and peak rate

**Table 3. Single-predictor regression results (outcome ~ predictor + subject fixed effect), n = 15.**

| Outcome | Predictor | Coefficient | SE | Permutation *p* (5,000 perms) | Parametric *p* (OLS) |
|---|---|---:|---:|---:|---:|
| Spatial information | Elapsed days | −0.0172 | 0.0095 | 0.114 | 0.102 |
| Spatial information | Cumulative cooling (s) | −0.00069 | 0.00037 | 0.083 | 0.089 |
| Field count | Elapsed days | −0.112 | 0.049 | 0.069 | 0.044 |
| **Field count** | **Cumulative cooling (s)** | **−0.00487** | **0.00176** | **0.022** | **0.020** |
| Peak rate (Hz) | Elapsed days | +0.627 | 0.285 | 0.069 | 0.052 |
| **Peak rate (Hz)** | **Cumulative cooling (s)** | **+0.0265** | **0.0105** | **0.044** | **0.031** |

Place-field count decreased, and peak firing rate increased, with cumulative prior cooling exposure (Figure 2), both reaching conventional significance under permutation testing and in close agreement with parametric OLS. Spatial information showed a trend in the same direction (decreasing with cumulative exposure) but did not reach significance (p = 0.083–0.089 across methods). The consistent direction of effect across all three outcome measures — fewer, sharper fields with preserved-to-declining overall information content — is a coherent pattern rather than three independent draws, though we note this coherence qualitatively rather than testing it as a joint hypothesis.

We note that field count in this dataset is bounded below at 3 by construction (the minimum contiguous-bin threshold for a detected field) and a majority of sessions with substantial cumulative cooling exposure sit at exactly this floor (Table 2); the effect is therefore better characterized as "high initial field counts collapse to a floor with cumulative exposure" than as a smooth continuous decline, a distinction relevant to interpreting effect size.

### The elapsed-time/cooling-dose confound

As described in *Methods*, cumulative cooling exposure and elapsed experimental days are highly correlated within subject in this dataset (MS21: r = 0.90; MS22: r = 0.977), because the cooling manipulation was applied in nearly every session. This is a structural property of the source experiment's design, not a limitation of sample size, and it cannot be resolved by analyzing additional sessions from the same design. We therefore cannot statistically attribute the observed trend to chronic theta disruption specifically, as opposed to any other process that changes monotonically with time in this experiment (e.g., electrode drift, task overtraining, or tissue response to the chronic implant). Both single-predictor models are reported in Table 3 to make this ambiguity transparent rather than to imply the two predictors are independently confirmed effects.

## Discussion

Petersen and Buzsáki (2020) established that a single, transient episode of medial septal cooling disrupts theta timing and behavioral choice accuracy without measurably altering the spatial firing correlates of hippocampal place cells within that session. Our reanalysis of the same public dataset asks a complementary question at a longer timescale: across the many sessions and weeks typical of such a chronic recording experiment, does the cumulative history of these transient disruptions leave a trace on baseline (pre-cooling) place-field properties? We find preliminary evidence that it may — place-field count declines and peak firing rate rises with cumulative cooling exposure — while overall spatial information content shows only a non-significant trend in the same direction. This pattern is consistent with a interpretation in which chronic theta disruption reshapes place fields toward fewer, sharper, higher-amplitude fields without a correspondingly large change in the information those fields carry about position, but we regard this as a hypothesis for future targeted testing rather than an established result.

The central limitation of this study is structural rather than statistical: cumulative cooling dose and elapsed experimental time are closely confounded in the source dataset's design, because cooling was applied in nearly every recording session. No reanalysis of this particular dataset, however large, can fully dissociate "effect of repeated theta disruption" from "effect of time since the experiment began." We report this explicitly rather than resolve it by predictor selection, and we consider identifying this confound to be a contribution of this work in its own right: it specifies exactly what a properly powered follow-up experiment would need to do differently — namely, interleave cooling and sham (no-cooling) sessions across time, so that dose and elapsed time vary independently and their effects can be separated.

A second limitation is sample size. Fifteen sessions across four subjects is small by the standards of a definitive claim, though we note the two significant results agree closely between an assumption-light permutation test and parametric OLS, which is reassuring evidence against those specific results being an artifact of small-sample estimation instability (in contrast to the joint two-predictor model, which we explicitly discarded for exactly this instability). We also note that the achievable sample size was constrained by the structure of the public data release itself: 35 of 50 attempted sessions (70%) lacked a usable spike-sorted `units` table, a `trials` table, or cooling-epoch annotations, and this ceiling — rather than any filtering choice we made — is the binding constraint on statistical power available from this dataset as currently archived.

Third, we caution that place-field count as computed here is bounded below by a minimum-field-size threshold, and most higher-cooling-exposure sessions in our sample sit at that floor; the effect is more accurately described as a collapse toward a floor than a smooth monotonic decline, which should temper claims about effect size or dose-response shape.

Finally, we highlight a methodological finding independent of the biological question: an early version of our peak-rate computation was confounded by this dataset's high (~120 Hz) position-sampling rate, which allowed single-sample occupancy bins coinciding with one incidental spike to register spurious ~120 Hz "peaks." We identified and corrected this via a minimum-occupancy exclusion criterion (Methods) before the results reported here. We mention this specifically because it is the kind of artifact that is easy to miss in any place-field pipeline built on high-sampling-rate position data, and because our corrected results differ non-trivially (and, in this case, more favorably) from the uncorrected ones — underscoring the value of occupancy-based quality control as standard practice.

## Limitations Summary

- Cumulative cooling exposure and elapsed time are confounded by experimental design; results cannot be attributed to cooling specifically.
- n = 15 sessions, 4 subjects; one subject (MS10) contributed no usable data.
- Field count is floor-bounded (minimum 3 contiguous bins), affecting interpretation of dose-response shape.
- Analysis is restricted to baseline (pre-cooling) trials within each session; we did not attempt cross-day matching of individual neurons (see *Code and data availability* for why this was assessed and ruled out).

## Data and Code Availability

Source data: DANDI Archive, Dandiset 000059, version 0.230907.2101 (https://dandiarchive.org/dandiset/000059), streamed directly with no local copy retained beyond per-session analysis caches. All analysis code, cached session-level statistics, and figure-generation scripts are available at https://github.com/sisyyphuss/theta-disruption-place-fields. Place-field-computation code is adapted from the openly licensed NCBS course notebooks (https://github.com/edeno/ncbs-neural-circuits-navigation); all cooling-dose extraction, statistical modeling, and permutation-testing code is original to this work.

## Acknowledgments

We thank Eric Denovellis for the openly available NCBS course materials on which the data-loading and place-field methodology in this work is built, and Peter Petersen and György Buzsáki for the original data collection and public archival of Dandiset 000059.

## References

Muller RU, Kubie JL (1987). The firing of hippocampal place cells predicts the future position of freely moving rats. *Journal of Neuroscience*, 7(7):1951–1968.

O'Keefe J, Dostrovsky J (1971). The hippocampus as a spatial map: preliminary evidence from unit activity in the freely-moving rat. *Brain Research*, 34(1):171–175.

Petersen PC, Buzsáki G (2020). Cooling of Medial Septum Reveals Theta Phase Lag Coordination of Hippocampal Cell Assemblies. *Neuron*, 107(4):731–744.e3. https://doi.org/10.1016/j.neuron.2020.05.023

Rübel O, Tritt A, Ly R, et al. (2022). The Neurodata Without Borders ecosystem for neurophysiological data science. *eLife*, 11:e78362.

Skaggs WE, McNaughton BL, Gothard KM, Markus EJ (1993). An information-theoretic approach to deciphering the hippocampal code. *Advances in Neural Information Processing Systems*, 5:1030–1037.

---

## Figure Legends

**Figure 1.** Session yield per subject. Of 50 attempted sessions across five subjects in Dandiset 000059, 15 met inclusion criteria (well-isolated units present, trials table with cooling-state annotation present, ≥30 baseline spikes for at least one unit). Subject MS10 contributed no usable sessions.

**Figure 2.** Session-level place-field statistics (rows: spatial information, field count, peak rate) plotted against elapsed days since each animal's first session (left column) and cumulative prior cooling exposure (right column). Points are individual sessions, colored by subject and connected in chronological order within subject. Field count and peak rate both show a significant relationship with cumulative cooling exposure (Table 3); spatial information shows a non-significant trend in the same direction.

*(Figures: `figures/fig1_session_yield.png`, `figures/fig2_trends.png`)*
