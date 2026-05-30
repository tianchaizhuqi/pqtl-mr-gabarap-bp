# ===========================================================================
# Positive control coloc: NPPA, NPPB, AGT, REN, ACE2 × 4 BP GWAS
# ===========================================================================
library(data.table)
library(coloc)
library(ieugwasr)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ0aWFuY2hhaXpodXFpQGdtYWlsLmNvbSIsImlhdCI6MTc3OTc1NDY2MSwiZXhwIjoxNzgwOTY0MjYxfQ.O7wUvLZC46uFLdGOJXZO9_NqKYNRMz9CPNZ24ufdnsHJekgAUykzfcHa7dFITU9yjmlfNO_4DZoLbMJNGKC-stvFBmm6Di420rhgaD7wydKXTsULKJIgXd8NuPnXq0c9v-86aTQqgcYGIeBzAc_KeXYhRxxqevXzBzEm9IOJVwAZek96NM5L_Ez5t6tEm7f7dHArQxDM_VtHiLc8qpS3jXsrWeZqTKA6zsArC_25k60ZEYoUELQeaI2G5Mw0qf9bzMH5SSr_SB7yaFxTf48hc6HYdODmyCj3iULCifJs1Td8_3wMqApJuCkgxkeuImg42yNmZgZqtEfWgcir5NqFYw")

setwd("d:/桌面/测试/dkd_multiomics_mr")
pqtl_dir <- "D:/多组学MR数据/pqtl/final somascan smp"

# Positive control proteins
controls <- list(
  NPPA = list(file = file.path(pqtl_dir, "Proteomics_SMP_PC0_5443_62_NPPA_ANP_10032022.txt.gz"), seqid = 5443),
  NPPB = list(file = file.path(pqtl_dir, "Proteomics_SMP_PC0_16751_15_NPPB_BNP_10032022.txt.gz"), seqid = 16751),
  AGT  = list(file = file.path(pqtl_dir, "Proteomics_SMP_PC0_3484_60_AGT_Angiotensinogen_10032022.txt.gz"), seqid = 3484),
  REN  = list(file = file.path(pqtl_dir, "Proteomics_SMP_PC0_3396_54_REN_Renin_10032022.txt.gz"), seqid = 3396),
  ACE2 = list(file = file.path(pqtl_dir, "Proteomics_SMP_PC0_2805_6_ACE2_ACE2_10032022.txt.gz"), seqid = 2805),
  ACE  = list(file = file.path(pqtl_dir, "Proteomics_SMP_PC0_10714_7_ACE_ACE_10032022.txt.gz"), seqid = 10714)
)

# 4 BP GWAS outcomes
outcomes <- list(
  SBP  = list(id = "ebi-a-GCST90025981", N = 422713),
  DBP  = list(id = "ebi-a-GCST90025968", N = 422713),
  HTN  = list(id = "ukb-b-12493",       N = 463010),
  MHTN = list(id = "ukb-b-18167",       N = 426391)
)

results_all <- list()

for(gene_name in names(controls)) {
  ctrl <- controls[[gene_name]]
  cat(sprintf("\n========== %s ==========\n", gene_name))

  # Load pQTL data
  dat <- fread(ctrl$file)
  dat <- dat[!is.na(rsids) & rsids != "NA" & rsids != "." & rsids != ""]

  # Get lead cis-pQTL variant
  top_snp <- dat[which.min(Pval)]
  lead_rsid <- top_snp$rsids[1]
  lead_pos <- top_snp$Pos[1]

  cat(sprintf("  Lead SNP: %s at chr%s:%d, P=%.2e, F=%.1f\n",
    lead_rsid, top_snp$Chrom[1], lead_pos,
    top_snp$Pval[1], (top_snp$Beta[1]/top_snp$SE[1])^2))

  # For each BP GWAS outcome, run coloc
  for(oc_name in names(outcomes)) {
    oc <- outcomes[[oc_name]]
    cat(sprintf("  %s...", oc_name))

    # Query GWAS in cis-region (±200kb)
    cis_region <- sprintf("%s:%d-%d", gsub("chr", "", top_snp$Chrom[1]), lead_pos - 200000, lead_pos + 200000)
    gwas <- tryCatch(
      associations(variants = cis_region, id = oc$id),
      error = function(e) NULL
    )

    if(is.null(gwas)) { cat("API error\n"); next }
    gwas <- as.data.frame(gwas)

    # Match SNPs
    common <- intersect(dat$rsids, gwas$rsid)
    n_common <- length(common)
    cat(sprintf(" %d common,", n_common))

    if(n_common < 20) { cat(" too few\n"); next }

    # Build coloc datasets
    pqtl_idx <- match(common, dat$rsids)
    gwas_idx <- match(common, gwas$rsid)

    d1 <- list(pvalues = dat$Pval[pqtl_idx], N = dat$N[pqtl_idx[1]], MAF = dat$ImpMAF[pqtl_idx], type = "quant", snp = common)
    d2 <- list(pvalues = gwas$p[gwas_idx], N = oc$N, MAF = rep(0.3, n_common), type = "quant", snp = common)

    cres <- coloc.abf(dataset1 = d1, dataset2 = d2)
    s <- cres$summary
    pph4 <- as.numeric(s["PP.H4.abf"])
    pph3 <- as.numeric(s["PP.H3.abf"])

    verdict <- if(pph4 > 0.75) "*** PASSED ***" else if(pph4 > 0.5) "SUGGESTIVE" else "NEGATIVE"
    cat(sprintf(" H4=%.4f H3=%.4f %s\n", pph4, pph3, verdict))

    results_all[[length(results_all) + 1]] <- data.frame(
      gene = gene_name,
      snp = lead_rsid,
      outcome = oc_name,
      outcome_id = oc$id,
      n_common = n_common,
      pqtl_p = top_snp$Pval[1],
      PP.H0 = as.numeric(s["PP.H0.abf"]),
      PP.H1 = as.numeric(s["PP.H1.abf"]),
      PP.H2 = as.numeric(s["PP.H2.abf"]),
      PP.H3 = pph3,
      PP.H4 = pph4,
      stringsAsFactors = FALSE
    )
  }
}

# Save results
res <- do.call(rbind, results_all)
res <- res[order(-res$PP.H4), ]
write.csv(res, "results/tables/positive_control_coloc.csv", row.names = FALSE)
cat("\n========== SUMMARY ==========\n")
cat(sprintf("Total tests: %d\n", nrow(res)))
cat(sprintf("PASSED (H4>0.75): %d\n", sum(res$PP.H4 > 0.75)))
cat(sprintf("SUGGESTIVE (H4>0.5): %d\n", sum(res$PP.H4 > 0.5 & res$PP.H4 <= 0.75)))
print(res[, c("gene", "outcome", "n_common", "PP.H4")], row.names = FALSE)
cat("\n[Done]\n")
