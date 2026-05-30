# pQTL-MR Colocalization: GABARAP/GABARAPL1 locus and blood pressure

Analysis scripts for "GABARAP/GABARAPL1 locus pQTL colocalizes with blood pressure traits in Europeans" (Communications Biology, under review).

## Overview

Pathway-based cis-pQTL Mendelian randomization screen of 30 autophagy-related proteins, combined with multi-method Bayesian colocalization (coloc.abf, coloc.susie with 1000 Genomes LD, coloc.signals), positive control calibration, SuSiE fine-mapping, and cross-ancestry assessment.

## Requirements

R >= 4.0 with packages:
- `data.table`, `coloc` (>= 5.2.3), `susieR`, `ieugwasr`, `TwoSampleMR`
- `ggplot2`, `viridis`, `RColorBrewer` (figures only)

## Key Scripts

| Script | Description |
|--------|-------------|
| `coloc_positive_controls.R` | Positive control colocalization (6 proteins x 4 BP GWAS) |
| `coloc_susie_real_ld.R` | coloc.susie with real 1000G LD matrix |
| `coloc_susie_gabarap.R` | coloc.signals conditional analysis |
| `conditional_coloc_gabarap.R` | GABARAP vs GABARAPL1 conditional coloc |
| `generate_figures_nature.R` | Nature-quality figures |
| `generate_supp_tables.R` | Supplementary tables S1-S17 |
| `proteome_wide_mr_screen.R` | Proteome-wide MR screen (1124 proteins) |
| `coloc_followup_candidates.R` | Coloc follow-up for top MR candidates |

## Data Availability

- GWAS: IEU OpenGWAS (https://gwas.mrcieu.ac.uk)
- pQTL: deCODE genetics (https://download.decode.is)
- LD reference: 1000 Genomes Phase 3 EUR
- Supplementary tables: S1-S17 in `supplementary/`

## License

MIT

## Citation

Zhu Q. GABARAP/GABARAPL1 locus pQTL colocalizes with blood pressure traits in Europeans. *Commun Biol*. 2026. Code: https://github.com/tianchaizhuqi/pqtl-mr-gabarap-bp
