# Proteome-wide cis-pQTL MR screen for BP outcomes
# Uses 1,124 deCODE protein GWAS in OpenGWAS (prot-c-*)
# Step 1: For each protein, get lead cis-pQTL from OpenGWAS tophits
# Step 2: Run Wald ratio MR against 4 BP outcomes
# Step 3: Rank proteins by MR significance
# Step 4: Save results for coloc candidate selection

library(ieugwasr)
library(data.table)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ0aWFuY2hhaXpodXFpQGdtYWlsLmNvbSIsImlhdCI6MTc3OTc1NDY2MSwiZXhwIjoxNzgwOTY0MjYxfQ.O7wUvLZC46uFLdGOJXZO9_NqKYNRMz9CPNZ24ufdnsHJekgAUykzfcHa7dFITU9yjmlfNO_4DZoLbMJNGKC-stvFBmm6Di420rhgaD7wydKXTsULKJIgXd8NuPnXq0c9v-86aTQqgcYGIeBzAc_KeXYhRxxqevXzBzEm9IOJVwAZek96NM5L_Ez5t6tEm7f7dHArQxDM_VtHiLc8qpS3jXsrWeZqTKA6zsArC_25k60ZEYoUELQeaI2G5Mw0qf9bzMH5SSr_SB7yaFxTf48hc6HYdODmyCj3iULCifJs1Td8_3wMqApJuCkgxkeuImg42yNmZgZqtEfWgcir5NqFYw")

setwd("d:/桌面/测试/dkd_multiomics_mr")

# ---- Load deCODE protein IDs ----
# Read from the pre-downloaded gwasinfo JSON
library(jsonlite)
cat("Loading GWAS metadata...\n")
gwas_meta <- fromJSON("data/gwasinfo_full.json")
cat(sprintf("Total GWAS entries: %d\n", length(gwas_meta)))

# Extract prot-c- IDs (deCODE cis-pQTL)
decode_ids <- names(gwas_meta)[grepl("^prot-c-", names(gwas_meta))]
cat(sprintf("deCODE cis-pQTL proteins: %d\n", length(decode_ids)))

# Save for reference
write.csv(data.frame(id = decode_ids), "data/decode_prot_c_ids.csv", row.names = FALSE)

# ---- Outcomes ----
outcomes <- list(
  SBP  = list(id = "ebi-a-GCST90025981", N = 422713),
  DBP  = list(id = "ebi-a-GCST90025968", N = 422713),
  HTN  = list(id = "ukb-b-12493",       N = 463010),
  MHTN = list(id = "ukb-b-18167",       N = 426391)
)

# ---- Screen: get lead cis-pQTL for each protein, run MR ----
results <- list()
n_total <- length(decode_ids)
batch_size <- 50  # Process in batches to avoid API overload

cat(sprintf("\nProcessing %d proteins in batches of %d...\n", n_total, batch_size))

for (i in seq(1, n_total, by = batch_size)) {
  batch_end <- min(i + batch_size - 1, n_total)
  batch_ids <- decode_ids[i:batch_end]
  cat(sprintf("\n=== Batch %d-%d of %d ===\n", i, batch_end, n_total))

  for (pid in batch_ids) {
    # Get protein name from metadata
    trait_name <- gwas_meta[[pid]]$trait
    if (is.null(trait_name)) trait_name <- pid

    # Get top hits for this protein (handle proteins with no significant cis-pQTL)
    tops <- tryCatch(
      suppressWarnings(tophits(id = pid)),
      error = function(e) NULL
    )
    if (is.null(tops)) next
    tops <- as.data.frame(tops)
    if (nrow(tops) == 0) next

    # Take the SNP with smallest P-value as lead cis-pQTL
    lead <- tops[which.min(tops$p), ]

    # Basic QC
    if (is.na(lead$rsid) || lead$rsid == "." || lead$rsid == "") next

    # Run MR against each BP outcome
    for (oc_name in names(outcomes)) {
      oc <- outcomes[[oc_name]]

      # Extract outcome association for this SNP
      gwas_snp <- tryCatch(
        associations(variants = lead$rsid, id = oc$id),
        error = function(e) NULL
      )
      if (is.null(gwas_snp) || nrow(gwas_snp) == 0) next
      gwas_snp <- as.data.frame(gwas_snp)

      # Align effect alleles
      if (gwas_snp$ea[1] != lead$ea && gwas_snp$ea[1] == lead$nea) {
        gwas_snp$beta[1] <- -gwas_snp$beta[1]
      } else if (gwas_snp$ea[1] != lead$ea && gwas_snp$ea[1] != lead$nea) {
        next  # Alleles don't match, skip
      }

      # Wald ratio
      mr_beta <- gwas_snp$beta[1] / lead$beta
      mr_se   <- abs(gwas_snp$se[1] / lead$beta)
      mr_pval <- 2 * pnorm(abs(mr_beta / mr_se), lower.tail = FALSE)

      # Calculate F-statistic
      F_stat <- (lead$beta / lead$se)^2

      results[[length(results) + 1]] <- data.frame(
        protein_id    = pid,
        protein_name  = trait_name,
        lead_snp      = lead$rsid,
        lead_chr      = lead$chr,
        lead_pos      = lead$position,
        lead_ea       = lead$ea,
        lead_nea      = lead$nea,
        lead_beta     = lead$beta,
        lead_se       = lead$se,
        lead_pval     = lead$p,
        lead_F        = F_stat,
        outcome       = oc_name,
        outcome_N     = oc$N,
        mr_beta       = mr_beta,
        mr_se         = mr_se,
        mr_pval       = mr_pval,
        stringsAsFactors = FALSE
      )
    }

    # Progress
    if (length(results) %% 200 == 0) {
      cat(sprintf("  Progress: %d proteins done, %d MR results\n",
                  which(decode_ids == pid), length(results)))
    }
  }

  # Save intermediate results
  if (length(results) > 0) {
    res_df <- do.call(rbind, results)
    write.csv(res_df, "results/tables/proteome_wide_mr_screen.csv", row.names = FALSE)
  }
}

# ---- Final summary ----
res_df <- do.call(rbind, results)
write.csv(res_df, "results/tables/proteome_wide_mr_screen.csv", row.names = FALSE)

cat(sprintf("\n========== FINAL ==========\n"))
cat(sprintf("Proteins screened: %d\n", n_total))
cat(sprintf("Proteins with MR results: %d\n", length(unique(res_df$protein_id))))
cat(sprintf("Total MR tests: %d\n", nrow(res_df)))

# Top hits per outcome
for (oc_name in names(outcomes)) {
  oc_res <- res_df[res_df$outcome == oc_name, ]
  oc_res <- oc_res[order(oc_res$mr_pval), ]
  cat(sprintf("\n=== Top 10 for %s ===\n", oc_name))
  if (nrow(oc_res) > 0) {
    print(head(oc_res[, c("protein_name", "lead_snp", "lead_F", "mr_beta", "mr_pval")], 10))
  }
}

cat("\n[Done]\n")
