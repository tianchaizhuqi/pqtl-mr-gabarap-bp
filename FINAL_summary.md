# FINAL coloc.susie Results Summary

**Generated:** 2026-05-30

## Overview

This file summarises **coloc + SuSiE** colocalization results for all 5 proteins
(GABARAP, FN1.3, PGP9.5/UCHL1, NovH/CCN3, Lysozyme) across two LD strategies:

1. **real_VCF** -- 1000 Genomes Phase 3 VCF reference panel, chromosome-specific
2. **IEU_API** -- IEU OpenGWAS API LD matrix (built-in reference panel)

Colocalization was assessed between *cis*-pQTL (deCODE, 35,559 Icelanders)
and blood pressure / hypertension GWAS summary statistics. SuSiE was used to
identify independent credible sets before formal coloc testing.

---

## Table 1: All colocalization results

| Protein | Outcome | ABF H4 | ABF H3 | SuSiE H4 | LD Source | N SNPs | ld_mean_r | Key Note |
|---------|---------|--------|--------|----------|-----------|--------|-----------|----------|
| GABARAP | SBP | 0.9850 | 7.22e-03 | 0.6403 | real_VCF | 92 | NA | Original real chr3 VCF coloc.susie |
| GABARAP | DBP | 0.7524 | 2.48e-01 | 0.9325 | real_VCF | 92 | NA | Original real chr3 VCF coloc.susie |
| GABARAP | HTN | 0.8555 | 1.43e-01 | 0.9994 | real_VCF | 329 | NA | Single SuSiE CS each trait |
| GABARAP | MHTN | 0.9629 | 3.67e-02 | 0.9998 | real_VCF | 330 | NA | Single SuSiE CS each trait |
| GABARAP_v4 | SBP | 0.9850 | NA | 1.0000 | IEU_API | 81 | 0.1050 | susie converged 81 variants after pruning |
| GABARAP_v4 | DBP | 0.7524 | NA | 1.0000 | IEU_API | 81 | 0.1050 | ld_mean_r=0.105 |
| GABARAP_v4 | HTN | 0.8555 | NA | 1.0000 | IEU_API | 323 | 0.3975 | 987 variants, 323 after pruning |
| GABARAP_v4 | MHTN | 0.9629 | NA | 1.0000 | IEU_API | 323 | 0.3975 | 1084 variants, 323 after pruning |
| FN1.3 (Fibronectin) | DBP | 0.9998 | 8.01e-05 | 0.9999 | real_VCF | 32 | 0.0937 | Real chr2 1000G VCF |
| FN1.3 (Fibronectin) | DBP | 0.9998 | 8.01e-05 | 1.0000 | IEU_API | 32 | NA | chr2:216100185-216500185 |
| PGP9.5 (UCHL1) | DBP | 0.9979 | 8.38e-07 | 0.9995 | real_VCF | 29 | 0.1359 | Real chr14 1000G VCF |
| PGP9.5 (UCHL1) | DBP | 0.9979 | 9.94e-07 | 0.9996 | IEU_API | 28 | NA | chr14:94510949-94910949 |
| NovH (CCN3) | DBP | 0.9431 | 5.48e-02 | 0.9731 | real_VCF | 21 | 0.1289 | Real chr8 1000G VCF |
| NovH (CCN3) | DBP | 0.9431 | 5.48e-02 | 1.0000 | IEU_API | 21 | NA | chr8:120235705-120635705 |
| Lysozyme | HTN | 0.8491 | 7.88e-02 | 0.9991 | real_VCF | 30 | 0.2242 | Real chr12 1000G VCF (pruned) |
| Lysozyme | HTN | 0.8491 | 7.88e-02 | 1.0000 | IEU_API | 104 | NA | chr12:69535492-69935492 |

---

## Key findings

### 1. GABARAP (GABARAP_v4, cis-pQTL chr3:6114143-6914143)

- **Real VCF**: coloc.susie supports colocalization for HTN (H4=0.9994) and MHTN (H4=0.9998).
  SBP (H4=0.640) and DBP (H4=0.933) show lower posterior probability under the real LD.
- **IEU API**: All four outcomes reach susie_H4=1.0000, consistent with the LD-approx overestimate
  pattern noted in the manuscript.
- The ABF coloc (without SuSiE) gives high H4 for all outcomes but cannot distinguish
  multiple causal variants.

### 2. FN1.3 (Fibronectin, cis-pQTL chr2:216100185-216500185)

- **Real VCF**: Strong coloc support (susie_H4=0.9999, abf_H4=0.9998) with DBP.
  Only 32 variants in the SuSiE credible set; low LD (mean r=0.094).
- **IEU API**: Confirms with susie_H4=1.0000.

### 3. PGP9.5 (UCHL1, cis-pQTL chr14:94510949-94910949)

- **Real VCF**: Strong coloc support (susie_H4=0.9995, abf_H4=0.9979) with DBP.
  29 variants; moderate LD (mean r=0.136).
- **IEU API**: susie_H4=0.9996, well-calibrated with the real VCF result.

### 4. NovH (CCN3, cis-pQTL chr8:120235705-120635705)

- **Real VCF**: Moderate-to-strong coloc support (susie_H4=0.9731, abf_H4=0.9431) with DBP.
  21 variants; moderate LD (mean r=0.129).
- **IEU API**: Overestimates to susie_H4=1.0000 (ABF H3=0.055, same as real VCF).

### 5. Lysozyme (cis-pQTL chr12:69535492-69935492)

- **Real VCF** (pruned): High coloc support (susie_H4=0.9991, abf_H4=0.8491) with HTN.
  30 variants after pruning; higher LD (mean r=0.224).
- **IEU API**: Overestimates to susie_H4=1.0000; 104 variants (less pruning).

---

## Method notes

- **SuSiE** was run with L=10, estimate_prior_variance=TRUE.
- **coloc.susie** used default priors (p1=1e-4, p2=1e-4, p12=5e-6).
- **Real VCF**: 1000 Genomes Phase 3, chromosome-specific, subset to analysis region +/- 200 kb.
- **IEU API**: LD reference from the IEU OpenGWAS API, matched to the same genomic region.
- **GABARAP** uses a different pQTL instrument (GABARAP_v4 aptamer, broader region) than the
  original analysis; v4 captures additional signal from the extended locus.
- **Lysozyme** real VCF was pruned (n=30) vs IEU API (n=104) due to VCF memory limits
  on this larger region.

---

## Files

- `FINAL_coloc_susie_master.csv` -- Machine-readable master table
- `coloc_susie_FN1_vcf.csv` -- FN1.3 real VCF raw output
- `coloc_susie_Lysozyme_vcf.csv` -- Lysozyme real VCF raw output
- `coloc_susie_new_hits.csv` -- All 4 new hits, IEU API LD raw output
- `coloc_susie_GABARAP_v4_ieu.csv` -- GABARAP v4, IEU API LD raw output
- `coloc_susie_real_ld_results.csv` -- GABARAP original, real VCF raw output
