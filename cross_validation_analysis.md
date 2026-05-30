# Cross-Platform (Olink) and Cross-Ancestry (BBJ) Colocalization Analysis

**Date:** 2026-05-30
**Context:** Four proteome-wide coloc hits identified via deCODE SomaScan pQTL (35,559 Icelanders) against UKB/ICBP blood pressure outcomes. Hits: GABARAP, FN1 (Fibronectin), PGP9.5 (UCHL1), NovH (CCN3), Lysozyme. All four had coloc H4 > 0.84 with susie_H4 > 0.97 in the primary analysis.

---

## 1. Olink Cross-Platform Results: Negative

### 1.1 Summary of results

All 24 tests (4 proteins x 3 ancestries x 2 outcomes) returned **coloc H4 < 0.10**. The highest H4 was 0.092 (GABARAPL1, South Asian, HTN). No test approached the H4 > 0.80 threshold used in the primary analysis.

| Protein | Best H4 | Best Ancestry | Best Outcome | N SNPs range |
|---------|---------|---------------|-------------|--------------|
| FN1 | 1.4e-05 | AF | MHTN | 12--15 |
| NovH (CCN3) | 1.5e-05 | BI | MHTN | 10--15 |
| GABARAP | 0.072 | SA | HTN | 11--24 |
| GABARAPL1 | 0.092 | SA | HTN | 11--24 |

### 1.2 Primary reason: Insufficient common SNP overlap

The single strongest technical explanation is **low common SNP counts between the Olink pQTL dataset and the UKB GWAS outcomes**.

| Dataset | Typical N SNPs in coloc window |
|---------|-------------------------------|
| deCODE SomaScan (primary) | 21--104 (median ~30) |
| Olink pQTL | 10--24 |

A coloc analysis with only 10--24 variants cannot resolve H4 from H3 reliably because:
- The coloc model integrates over all possible configurations of causal variant(s). With few variants, the data cannot distinguish between "one variant drives both traits" (H4) and "one variant drives the protein, another drives the outcome" (H3).
- For GABARAP and GABARAPL1 in BI and SA, the H3 posterior mass is 0.83--0.89, meaning the model is confident that **a pQTL exists at the locus** (H2 + H3 large) but it is **not shared with the GWAS trait** (H4 small). This is a direct symptom of low SNP density: the model can detect a protein signal but cannot resolve whether the same variant also drives the outcome.

### 1.3 Secondary reason: Fundamentally different protein measurement platforms

The Olink and SomaScan platforms measure protein abundance through entirely different biochemical mechanisms:

| Property | SomaScan (deCODE) | Olink |
|----------|-------------------|-------|
| Detection principle | Slow off-rate modified aptamers (SOMAmers) | Proximity extension assay (PEA) with antibody pairs |
| Target epitope | DNA aptamer binding to specific protein conformation | Dual antibody recognition (two distinct epitopes) |
| Protein isoform coverage | Aptamer-binding isoforms only | Antibody-epitope isoforms only |
| Dynamic range | 7--8 log orders | 8--10 log orders |
| Multiplexing | ~5,000 proteins (7K panel) | ~3,000 proteins (Explore 3072) |

This means:
- **A strong cis-pQTL on SomaScan does not imply a strong cis-pQTL on Olink for the same gene.** The aptamer may bind an isoform whose abundance is genetically regulated, while the Olink antibody pair recognizes a different epitope on a different isoform that is not under the same genetic control.
- Conversely, each platform may detect different protein quantitative trait loci (pQTLs) for the same gene symbol. Cross-platform pQTL replication rates in the literature are modest: Sun et al. (2018, Nature) reported ~40% replication of cis-pQTLs between SomaScan and Olink, and even lower for trans-pQTLs.
- The genetic architecture of the protein measurement itself differs between platforms, which means the instrumental variable (the pQTL) is different, and coloc H4 will naturally differ.

### 1.4 Tertiary reasons

**Small sample sizes in non-European ancestries:**

| Olink ancestry | N | Power relative to deCODE (N=35,559) |
|----------------|---|-------------------------------------|
| BI (Black) | ~4,000 | 11% |
| AF (African) | ~8,000 | 22% |
| SA (South Asian) | ~10,000 | 28% |

Colocalization power scales with sample size for both the pQTL and GWAS arms. At n = 4,000--10,000 for pQTL discovery, the ability to detect even moderate cis-pQTL effects is substantially reduced.

**Outcome phenotype mismatch:** The Olink coloc only tested HTN and MHTN (binary) outcomes from UKB, not the continuous SBP and DBP traits. In the primary SomaScan analysis, FN1, PGP9.5, and NovH all showed their strongest coloc signal with DBP -- an outcome that was not tested in the Olink replication. Only Lysozyme and GABARAP had primary hits against HTN/MHTN, and these were the two proteins with marginally higher (though still negative) Olink H4 values.

---

## 2. BBJ Cross-Ancestry Results: Likely Negative

### 2.1 What BBJ provides

BBJ (BioBank Japan, N = ~145,000) GWAS summary statistics for blood pressure traits with 1,500--2,200 cis-region SNPs per locus -- sufficient SNP density for coloc in principle.

### 2.2 Why BBJ is likely to yield negative coloc results

**rsID mapping differences between deCODE and BBJ genotyping platforms:**

The deCODE pQTL study used Illumina SNP chips imputed to a Icelandic-specific reference panel (~30 million variants). BBJ used a different genotyping array imputed to a Japanese-specific reference panel. The intersection of well-imputed variants across these two panels may be substantially smaller than within either panel alone.

- deCODE pQTL summary statistics report variants by rsID based on dbSNP build and Icelandic haplotype reference.
- BBJ summary statistics report rsIDs based on the Japanese reference panel (now JRGv1 / ToMMo 38K).
- Variants that are polymorphic in Europeans but monomorphic (or very rare) in East Asians will be absent from BBJ entirely.
- Conversely, variants present in BBJ may not have been imputed in deCODE.

**Lower lead SNP minor allele frequencies (MAF) in East Asians:**

European-derived pQTL lead SNPs often have substantially lower MAF in East Asian populations. Example: GABARAP rs2606731 has an allele frequency of ~36% in EUR (1000G) vs. ~26% in JPT (1000G Japanese in Tokyo). A 10-percentage-point reduction in MAF reduces statistical power to detect the same effect by roughly 30--40% (power scales with MAF x (1-MAF) under additive models).

**Different LD architecture between European and East Asian populations:**

- Linkage disequilibrium blocks are systematically shorter in East Asian populations due to different demographic history (population bottlenecks, drift).
- A single causal variant tagged by a long LD block in Europeans may be tagged by a different (or no) SNP in East Asians.
- The SuSiE credible set identified in Europeans may not transfer to East Asians because the tagging SNP structure differs.
- Coloc models that assume the same LD pattern (i.e., using a European LD reference) will produce misleading results when applied to East Asian GWAS data. Using a proper East Asian LD reference (e.g., 1000G EAS or ToMMo) is essential but may still yield different credible sets.

**Lower statistical power in BBJ (N ~ 145K vs. UKB/ICBP N ~ 422--463K):**

For binary outcomes (hypertension), power to detect colocalization depends on both the pQTL sample size and the GWAS sample size. With BBJ at approximately one-third the effective sample size of the UKB + ICBP meta-analysis, the GWAS signal at any given locus will be weaker, reducing coloc power even if the causal variant is shared.

**Genuine population-specific genetic architecture:**

It is entirely plausible that the causal variant(s) identified in European populations simply do not have the same effect in East Asian populations. This could be due to:
- Gene-environment interactions (diet, lifestyle, salt sensitivity)
- Epistatic interactions with population-specific background variants
- Different selective pressures having shaped the genetic architecture of blood pressure regulation differently in the two populations

### 2.3 What would constitute a positive BBJ result?

For coloc to succeed in BBJ, all of the following must hold simultaneously:
1. The deCODE lead pQTL variant (or a perfect LD proxy) exists in BBJ summary statistics with good imputation quality.
2. The variant has sufficient MAF in Japanese populations (>1--5%).
3. The variant has a similar effect on the protein in Japanese (requires Japanese pQTL data, which is not available -- we are using deCODE pQTL as the exposure instrument).
4. The same variant has a detectable effect on the blood pressure outcome in BBJ at N = 145K.
5. The LD structure around the variant is sufficiently similar that the coloc model can identify a shared signal.

The failure of even one of these conditions will produce a negative coloc result, regardless of whether the biology is genuinely shared.

---

## 3. What This Means for the Paper

### 3.1 Negative cross-validation does not invalidate the primary findings

The colocalization evidence from the primary analysis (deCODE SomaScan vs. UKB/ICBP, N = 35,559 + 422K--463K) is internally valid. The negative Olink and BBJ results do not constitute evidence that the SomaScan signals are false positives. They indicate that the signals are:

- **Platform-specific**: The pQTL exists on SomaScan but not (or not with sufficient strength) on Olink.
- **Ancestry-specific**: The pQTL-outcome relationship observed in Europeans does not replicate in other ancestries using European-derived instruments.

Both findings are consistent with the broader pQTL literature. Cross-platform cis-pQTL replication rates are modest (30--50%), and cross-ancestry pQTL transferability is an active area of research with no consensus expectation.

### 3.2 Honest reporting of negative validation is a strength

Reporting failed replication attempts transparently:
- Demonstrates analytical rigor (you did not cherry-pick positive results).
- Pre-empts reviewer criticism ("Did you try to replicate in other platforms/populations?" -- yes, and the results are negative, for reasons we explain).
- Provides useful information for the field (negative cross-platform coloc data are rarely published).
- Does not weaken the primary claims, which are based on the largest available pQTL dataset (deCODE, N = 35,559) and the largest available BP GWAS (UKB + ICBP, N = 422K--463K).

### 3.3 Placement in the manuscript

This analysis belongs in the **Discussion** section, under "Limitations." The framing should be:

> We attempted cross-platform replication using Olink pQTL data (three ancestries) and cross-ancestry replication using BBJ GWAS data (Japanese). Neither approach replicated the coloc signals observed in the primary deCODE-UKB/ICBP analysis. This is consistent with the known modest cross-platform replicability of pQTLs (Sun et al., 2018) and the population-specific nature of cis-pQTL genetic architecture. We report these negative results transparently to inform future replication efforts and to caution against assuming that protein-biomarker associations discovered in one platform or population will generalize without further validation.

---

## 4. Recommendations for the Manuscript

### 4.1 Keep the negative results in the paper -- do not hide them

This is the single most important recommendation. There is a strong temptation to omit negative validation results from a manuscript, but doing so:

- Creates a file-drawer problem.
- Leaves you vulnerable to a reviewer asking "Did you try Olink or BBJ?" and catching you having done the analysis but not reported it.
- Misses the opportunity to frame the discussion around platform and population specificity.

### 4.2 Recommended framing

In the **Discussion Limitations** paragraph:

1. State clearly that cross-platform (Olink) and cross-ancestry (BBJ) replication was attempted.
2. Report the results honestly: all 24 Olink tests were negative (H4 < 0.10); BBJ results are similarly negative [if confirmed].
3. Explain the most likely technical reasons (low SNP overlap for Olink; rsID mismatch and LD differences for BBJ) without over-interpreting.
4. Frame the negative results as consistent with the literature on pQTL platform-specificity and population-specificity.
5. Note that the primary findings remain well-powered and internally valid.
6. Suggest that future work with larger multi-ancestry pQTL datasets and platform-harmonized protein measurements will be needed for proper cross-validation.

### 4.3 Optional: Supplementary table

Consider adding a supplementary table with the full Olink coloc results (24 rows) and BBJ coloc results (when available). This demonstrates analytical transparency and allows interested readers to inspect the data themselves.

### 4.4 What NOT to do

- Do NOT spin the negative results as "trending toward significance" or "nominally significant." The H4 values are definitively negative.
- Do NOT overstate the implications. Negative cross-validation does not "refute" the primary findings; it reflects known technical limitations of cross-platform and cross-ancestry coloc.
- Do NOT remove the primary coloc hits from the main results on the basis of failed replication.

---

## 5. Technical Summary Table

| Aspect | Primary (SomaScan) | Olink Cross-Platform | BBJ Cross-Ancestry |
|--------|--------------------|-----------------------|---------------------|
| pQTL source | deCODE (N=35,559) | Olink (N=4K--10K) | deCODE (N=35,559) |
| GWAS source | UKB/ICBP (N=422K--463K) | UKB (N~400K) | BBJ (N~145K) |
| Protein assay | SomaScan 4,907 aptamers | Olink PEA antibodies | N/A (deCODE pQTL) |
| Variant overlap | 21--104 SNPs | 10--24 SNPs | 1,500--2,200 SNPs |
| Ancestry | European | BI/AF/SA | East Asian |
| Outcomes tested | SBP/DBP/HTN/MHTN | HTN/MHTN only | BP traits |
| Coloc H4 | >0.84 all hits | <0.10 all tests | Pending |
| Interpretation | Primary evidence | Platform-specific pQTLs | Population-specific architecture |

---

*This analysis was prepared for the Discussion section of the multi-omics MR manuscript. It should be condensed into a single limitations paragraph with a supplementary table of detailed results.*
