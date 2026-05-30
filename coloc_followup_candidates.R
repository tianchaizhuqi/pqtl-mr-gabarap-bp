# Step 2: coloc follow-up for top proteome-wide MR hits
# Run after proteome_wide_mr_screen.R completes
# Selects top N proteins and downloads cis-region data for coloc

library(ieugwasr)
library(data.table)
library(coloc)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ0aWFuY2hhaXpodXFpQGdtYWlsLmNvbSIsImlhdCI6MTc3OTc1NDY2MSwiZXhwIjoxNzgwOTY0MjYxfQ.O7wUvLZC46uFLdGOJXZO9_NqKYNRMz9CPNZ24ufdnsHJekgAUykzfcHa7dFITU9yjmlfNO_4DZoLbMJNGKC-stvFBmm6Di420rhgaD7wydKXTsULKJIgXd8NuPnXq0c9v-86aTQqgcYGIeBzAc_KeXYhRxxqevXzBzEm9IOJVwAZek96NM5L_Ez5t6tEm7f7dHArQxDM_VtHiLc8qpS3jXsrWeZqTKA6zsArC_25k60ZEYoUELQeaI2G5Mw0qf9bzMH5SSr_SB7yaFxTf48hc6HYdODmyCj3iULCifJs1Td8_3wMqApJuCkgxkeuImg42yNmZgZqtEfWgcir5NqFYw")

setwd("d:/桌面/测试/dkd_multiomics_mr")

# ---- Load MR screen results ----
mr <- fread("results/tables/proteome_wide_mr_screen.csv")
cat(sprintf("Loaded %d MR results for %d unique proteins\n",
            nrow(mr), length(unique(mr$protein_id))))

# ---- Selection criteria for coloc follow-up ----
# 1. F > 10 (adequate instrument strength)
# 2. Nominal MR P < 0.05 for at least one BP outcome
# 3. Exclude proteins already analyzed (autophagy pathway)

mr[, significant := mr_pval < 0.05]
mr[, strong_instrument := lead_F > 10]

# Best MR result per protein
protein_summary <- mr[, .(
  min_mr_pval = min(mr_pval),
  max_F = max(lead_F),
  n_sig_outcomes = sum(significant),
  best_outcome = outcome[which.min(mr_pval)],
  lead_snp = lead_snp[1],
  protein_name = protein_name[1]
), by = protein_id]

# Filter for coloc candidates
candidates <- protein_summary[
  max_F > 10 & n_sig_outcomes >= 1
][order(min_mr_pval)]

cat(sprintf("\nColoc candidates (F>10, P<0.05): %d\n", nrow(candidates)))
print(head(candidates, 30))

# Save candidate list
write.csv(candidates, "results/tables/coloc_candidates.csv", row.names = FALSE)

# ---- Download cis-region data for top N candidates ----
top_n <- min(100, nrow(candidates))
cat(sprintf("\nWill download cis-region data for top %d candidates\n", top_n))

outcomes <- list(
  SBP  = list(id = "ebi-a-GCST90025981", N = 422713),
  DBP  = list(id = "ebi-a-GCST90025968", N = 422713),
  HTN  = list(id = "ukb-b-12493",       N = 463010),
  MHTN = list(id = "ukb-b-18167",       N = 426391)
)

coloc_results <- list()

for (i in 1:top_n) {
  pid <- candidates$protein_id[i]
  lead_snp <- candidates$lead_snp[i]
  prot_name <- candidates$protein_name[i]

  cat(sprintf("\n=== [%d/%d] %s (%s) ===\n", i, top_n, prot_name, pid))

  # Get lead SNP position from MR results
  snp_info <- mr[protein_id == pid][1]
  cis_region <- sprintf("%s:%d-%d",
    snp_info$lead_chr,
    snp_info$lead_pos - 200000,
    snp_info$lead_pos + 200000)

  # Get cis-pQTL data from OpenGWAS
  pqtl <- tryCatch(
    associations(variants = cis_region, id = pid),
    error = function(e) { cat(sprintf("  pQTL API error: %s\n", e$message)); return(NULL) }
  )
  if (is.null(pqtl) || nrow(pqtl) < 20) { cat("  Too few pQTL SNPs\n"); next }
  pqtl <- as.data.frame(pqtl)

  # Run coloc for each outcome
  for (oc_name in names(outcomes)) {
    oc <- outcomes[[oc_name]]

    gwas <- tryCatch(
      associations(variants = cis_region, id = oc$id),
      error = function(e) { cat(sprintf("  GWAS API error (%s): %s\n", oc_name, e$message)); return(NULL) }
    )
    if (is.null(gwas) || nrow(gwas) < 20) next
    gwas <- as.data.frame(gwas)

    # Match by rsID
    common <- intersect(pqtl$rsid, gwas$rsid)
    if (length(common) < 20) {
      cat(sprintf("  %s: %d common SNPs, skip\n", oc_name, length(common)))
      next
    }

    pqtl_idx <- match(common, pqtl$rsid)
    gwas_idx  <- match(common, gwas$rsid)

    # coloc.abf
    d1 <- list(pvalues = pqtl$p[pqtl_idx], N = 35275, MAF = rep(0.3, length(common)),
               type = "quant", snp = common)
    d2 <- list(pvalues = gwas$p[gwas_idx], N = oc$N, MAF = rep(0.3, length(common)),
               type = "quant", snp = common)

    cres <- tryCatch(
      coloc.abf(dataset1 = d1, dataset2 = d2),
      error = function(e) { cat(sprintf("  coloc error: %s\n", e$message)); return(NULL) }
    )
    if (is.null(cres)) next

    h4 <- as.numeric(cres$summary["PP.H4.abf"])
    h3 <- as.numeric(cres$summary["PP.H3.abf"])

    coloc_results[[length(coloc_results) + 1]] <- data.frame(
      protein_id   = pid,
      protein_name = prot_name,
      lead_snp     = lead_snp,
      outcome      = oc_name,
      n_common     = length(common),
      coloc_H4     = h4,
      coloc_H3     = h3,
      stringsAsFactors = FALSE
    )

    cat(sprintf("  %s: H4=%.4f H3=%.4f (%d SNPs)\n", oc_name, h4, h3, length(common)))
  }

  # Save intermediate
  if (length(coloc_results) > 0) {
    res_df <- do.call(rbind, coloc_results)
    write.csv(res_df, "results/tables/proteome_wide_coloc_results.csv", row.names = FALSE)
  }
}

# ---- Final summary ----
res_df <- do.call(rbind, coloc_results)
write.csv(res_df, "results/tables/proteome_wide_coloc_results.csv", row.names = FALSE)

cat(sprintf("\n========== FINAL SUMMARY ==========\n"))
cat(sprintf("Proteins with coloc: %d\n", length(unique(res_df$protein_id))))

# Proteins with H4 > 0.75
hits <- res_df[res_df$coloc_H4 > 0.75, ]
if (nrow(hits) > 0) {
  cat(sprintf("\nProteins with coloc H4 > 0.75: %d\n", nrow(hits)))
  print(hits[order(-hits$coloc_H4), ])
} else {
  cat("\nNo proteins with coloc H4 > 0.75 beyond GABARAP\n")
  # Show top coloc signals regardless
  cat("\nTop coloc signals (regardless of threshold):\n")
  top_coloc <- res_df[order(-res_df$coloc_H4), ]
  print(head(top_coloc, 20))
}

cat("\n[Done]\n")
