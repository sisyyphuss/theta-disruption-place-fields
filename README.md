# Cumulative Theta-Rhythm Disruption and Hippocampal Place-Field Plasticity

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Data: DANDI 000059](https://img.shields.io/badge/DANDI-000059-green.svg)](https://dandiarchive.org/dandiset/000059)

This repository provides the computational analysis pipeline, executed Jupyter notebook, and precomputed regression statistics for the study:

> **Cumulative Theta-Rhythm Disruption Predicts a Reduction in Hippocampal Place-Field Count: An Exploratory Reanalysis of Public Medial Septal Cooling Data**  
> **Saad Yusuf**  
> *Department of Computer Science, Jamia Millia Islamia, New Delhi, India*  
> Correspondence: `syed.saad7527@gmail.com` | ORCID: `0009-0000-2333-0078`  
> Manuscript preprint available in [`preprint/`](preprint/theta_disruption_place_fields_manuscript.pdf).

---

## Scientific Overview

The medial septum paces hippocampal theta oscillations (6–10 Hz). While transient, acute cooling of the medial septum slows theta and degrades spatial-alternation choices without disrupting the spatial firing fields of individual place cells within a single session (Petersen & Buzsáki, 2020, *Neuron*), whether *chronic, repeated* exposure to this disruption alters baseline place-field representation across multi-week timescales remains unexamined.

This project reanalyzes chronic dorsal CA1 extracellular recordings from rats performing a spatial working-memory task ([DANDI:000059](https://dandiarchive.org/dandiset/000059)). Reconstructing baseline (pre-cooling) rate maps and place-field statistics across 15 usable sessions spanning up to six weeks per animal, we evaluate whether cumulative prior cooling exposure predicts long-term place-field coding metrics.

### Primary Findings

1. **Place-Field Count Decline:** Place-field count per unit declines significantly as a function of cumulative prior cooling exposure (permutation $p = 0.022$, parametric OLS $p = 0.020$).
2. **Elevated Peak Firing Rates:** Peak firing rate increases significantly with cumulative cooling exposure (permutation $p = 0.044$, parametric OLS $p = 0.031$), consistent with fewer but sharper fields.
3. **Preserved Information Content:** Overall Skaggs spatial information shows a modest non-significant downward trend (permutation $p = 0.083$, parametric OLS $p = 0.089$).
4. **Structural Confound Identification:** In the original experimental protocol, cumulative cooling dose is strongly correlated with elapsed calendar time within animals ($r = 0.90\text{--}0.98$). We report this limitation explicitly to guide future interleaved experimental designs.

---

## Repository Structure

```
theta-disruption-place-fields/
├── src/                                  # Core computational analysis modules
│   ├── data_loading.py                   # DANDI Archive streaming and NWB session parsing
│   ├── cooling_dose.py                   # Cooling epoch extraction and cumulative dose tracking
│   ├── place_fields.py                   # Occupancy-normalized rate maps and spatial information
│   ├── pooled_regression.py              # Permutation testing and within-subject OLS models
│   ├── track_b_analysis.py               # Track B exploratory spatial analysis
│   ├── make_figures.py                   # Figure 1 (yield) and Figure 2 (trends) generation
│   ├── smoke_test.py                     # Pipeline verification test
│   └── run_pipeline.py                   # Master analysis pipeline runner and summary reporter
├── notebooks/                            # Executed Jupyter analysis notebooks
│   └── full_analysis.ipynb               # Fully executed end-to-end reproducible analysis notebook
├── results/                              # Figures and precomputed summary tables
│   ├── figures/                          # Publication figures (fig1_session_yield.png, fig2_trends.png)
│   ├── tables/                           # CSV session dataset and summary tables
│   └── cache/                            # Per-subject session CSVs and regression pickle (final_results.pkl)
├── preprint/                             # Final camera-ready manuscript
│   ├── theta_disruption_place_fields_manuscript.pdf # PDF manuscript
│   └── manuscript.md                     # Markdown manuscript draft
├── requirements.txt                      # Python library dependencies
├── LICENSE                               # MIT License
└── README.md
```

---

## Description of Analysis Modules (`src/`)

- **`data_loading.py`:** Streams `.nwb` sessions directly from Dandiset 000059 via `pynwb`/`remfile`, filters well-isolated units (`good`/`good2`), and extracts position tracking coordinates.
- **`cooling_dose.py`:** Parses trial intervals, isolates pre-cooling baseline periods, and tracks cumulative seconds of active medial-septal cooling.
- **`place_fields.py`:** Generates 2D occupancy-normalized rate maps ($2\text{ cm}$ bins, excluding bins with $<0.1\text{ s}$ occupancy to prevent high-frequency sampling artifacts), calculates Skaggs spatial information, and counts place fields.
- **`pooled_regression.py`:** Pools session-level metrics across subjects and fits within-subject OLS models with 5,000-iteration permutation null testing.
- **`make_figures.py`:** Renders publication-ready plots for session yield across subjects (Figure 1) and longitudinal trend lines across elapsed time and cooling dose (Figure 2).
- **`run_pipeline.py`:** Command-line driver that reproduces the session-yield statistics, regression models, and summary outputs.

---

## Data Source & Inclusion Criteria

- **Archive:** [DANDI Archive, Dandiset 000059](https://dandiarchive.org/dandiset/000059) (Petersen & Buzsáki, 2020).
- **Paradigm:** Figure-8 spatial alternation task with reversible medial septal cooling.
- **Included Cohort:** 15 sessions across 4 Long-Evans rats (`MS12`, `MS13`, `MS21`, `MS22`) that provided spike-sorted units, annotated cooling trial states, and $\ge 30$ baseline spikes.

---

## Reproduction & Usage

### 1. Environment Setup

```bash
git clone https://github.com/sisyyphuss/theta-disruption-place-fields.git
cd theta-disruption-place-fields
pip install -r requirements.txt
```

### 2. Running Pipelines

All precomputed session-level statistics and regression tables are stored in `results/tables/` and `results/cache/`. Execute the pipeline driver directly:

```bash
python src/run_pipeline.py
```

To recompute the 5,000-iteration permutation null models from scratch:

```bash
python src/run_pipeline.py --recompute
```

### 3. Interactive Notebook

Launch Jupyter Lab to inspect the complete reproducible walkthrough:

```bash
jupyter lab notebooks/full_analysis.ipynb
```

---

## Citation

If you use this pipeline or data tables in your research, please cite:

```bibtex
@article{yusuf2026theta,
  title={Cumulative Theta-Rhythm Disruption Predicts a Reduction in Hippocampal Place-Field Count: An Exploratory Reanalysis of Public Medial Septal Cooling Data},
  author={Yusuf, Saad},
  journal={bioRxiv},
  year={2026},
  note={Preprint}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
