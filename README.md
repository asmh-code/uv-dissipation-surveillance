# UV Dissipation and Mutational Signature Analysis for Multi-Modal Pandemic Surveillance

Companion code for:

> Harman, A. S. M. (2026). *Resilience, Entropy and Model Choice: Integrating
> Thermodynamic Dissipation Theory into Multi-Modal Pandemic Surveillance.*

## Overview

This repository contains the computational analyses reported in the paper:

1. **H3 Power Analysis** (`h3_power_analysis.py`) — Simulation-based power
   analysis for detecting UV-correlated G→T transversions (the 8-oxoguanine
   oxidative damage signature) in SARS-CoV-2 phylodynamic data.

2. **Virus-Specific UV Dissipation Analysis** (`virus_dissipation_analysis.py`)
   — Computation of UV absorption properties and H3 detectability scores for
   19 RNA viruses across seven families, under both direct UVB and
   photosensitised UVA/ROS photochemical models.

## Requirements

- Python 3.8+
- NumPy
- SciPy
- statsmodels
- matplotlib

Install dependencies:

```bash
pip install numpy scipy statsmodels matplotlib
```

## Usage

### H3 Power Analysis

```bash
python h3_power_analysis.py
```

Generates:
- Console output with full power grid across effect sizes (0–30%) and
  sequencing volumes (50–1000 sequences/week)
- `h3_power_results_v2.json` — machine-readable results
- Figures reproduced by running `h3_final_figure.py` (included)

The simulation uses:
- Published SARS-CoV-2 baseline mutational frequencies
- Seasonal UV variation (Northern Hemisphere sinusoidal + weather noise)
- Adversarial lineage turnover with seasonally correlated G→T baseline shifts
- Binomial GLM with overdispersion correction and lineage indicator controls
- 500–1000 simulations per parameter combination

Key result: A +5% UV-correlated increase in G→T transversions is detectable
with ≥80% power given ≥100 sequences/week over 104 weeks, even under
adversarial lineage confounding. False positive rates are well controlled
at 2–4% (nominal α = 0.05).

### Virus-Specific Dissipation Analysis

```bash
python virus_dissipation_analysis.py
```

Generates:
- Console output with dual-score rankings for 19 RNA viruses
- `virus_dissipation_figure_v2.pdf` and `.png` — publication figure
- Ranking comparison between direct UVB and photosensitised ROS models

The analysis uses:
- Published ribonucleotide molar absorption coefficients (Cavaluzzi & Borer, 2004)
- NCBI reference genome compositions
- Two detectability scores:
  - **Direct UVB:** ε₂₈₀/nt × G-fraction × genome length (laboratory context)
  - **Photosensitised ROS:** G-fraction × genome length (solar/environmental
    context, where photon capture is matrix-driven via organic photosensitisers)

Key result: Coronaviruses are the optimal phylodynamic targets under both
models (2.1× non-coronavirus mean under the spectrally appropriate ROS model),
while flaviviruses show the highest per-nucleotide target densities for
laboratory assays. The spectral correction strengthens the microlayer
hypothesis: under solar conditions, the organic photosensitiser content of
the microlayer is the rate-limiting variable for oxidative damage.

## Reproducibility

All results are exactly reproducible by running the scripts with the
dependencies listed above. No external datasets are required — the analyses
use published physical constants and reference genome compositions. Random
seeds are set for deterministic output.

## Citation

If you use this code, please cite:

```
Harman, A. S. M. (2026). Resilience, Entropy and Model Choice: Integrating
Thermodynamic Dissipation Theory into Multi-Modal Pandemic Surveillance.
[Journal, volume, pages — to be updated on publication]
```

## License

MIT License. See [LICENSE](LICENSE).

## Contact

Adem S. M. Harman — Independent Scholar
