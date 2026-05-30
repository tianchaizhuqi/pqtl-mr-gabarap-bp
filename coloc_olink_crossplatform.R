# =============================================================================
# Cross-platform Olink coloc for FN1, CCN3, GABARAP, GABARAPL1
# Olink UKB BI (EUR ~55K), Africa (~8K), SAsia (~10K) x BP GWAS
#
# FIXED 2026-05-30:
#   - Removed hg19 ICBP GWAS (ebi-a-GCST90025981/25968). Those use hg19 but
#     Olink is hg38, so cis-region queries returned wrong/zero SNPs and rsID
#     overlap was zero.
#   - Now only uses UKB-based hg38 GWAS: ukb-b-12493 (HTN), ukb-b-18167 (MHTN).
#   - Matching: position-based via chr:pos key, not rsID-based.
#     Both Olink (hg38) and UKB GWAS (hg38) share coordinates.
#   - Uses Olink "Name" column (chr:pos:ref:alt) as SNP identifiers for coloc.
#   - Uses real ImpMAF / eaf when available instead of fixed 0.3.
# =============================================================================
library(data.table)
library(coloc)
library(ieugwasr)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ0aWFuY2hhaXpodXFpQGdtYWlsLmNvbSIsImlhdCI6MTc3OTc1NDY2MSwiZXhwIjoxNzgwOTY0MjYxfQ.O7wUvLZC46uFLdGOJXZO9_NqKYNRMz9CPNZ24ufdnsHJekgAUykzfcHa7dFITU9yjmlfNO_4DZoLbMJNGKC-stvFBmm6Di420rhgaD7wydKXTsULKJIgXd8NuPnXq0c9v-86aTQqgcYGIeBzAc_KeXYhRxxqevXzBzEm9IOJVwAZek96NM5L_Ez5t6tEm7f7dHArQxDM_VtHiLc8qpS3jXsrWeZqTKA6zsArC_25k60ZEYoUELQeaI2G5Mw0qf9bzMH5SSr_SB7yaFxTf48hc6HYdODmyCj3iULCifJs1Td8_3wMqApJuCkgxkeuImg42yNmZgZqtEfWgcir5NqFYw")
setwd("d:/桌面/测试/dkd_multiomics_mr")

# ---- Protein definitions ----
proteins <- list(
  FN1 = list(
    name       = "FN1 (Fibronectin)",
    chr        = "2",
    pos        = 216300185,
    olink_file = c(
      "GBR_UKB_OLINK2_OID30787_FN1_Fibronectin_adjAgeSexPC_InvNorm_22122022.txt.gz",
      "GBR_UKB_Africa_OLINK2_OID30787_FN1_Fibronectin_adjAgeSexPC_InvNorm_22122022.txt.gz",
      "GBR_UKB_SAsia_OLINK2_OID30787_FN1_Fibronectin_adjAgeSexPC_InvNorm_22122022.txt.gz"
    ),
    olink_label = c("OLINK_BI", "OLINK_AF", "OLINK_SA"),
    olink_N     = c(55000, 8000, 10000)
  ),
  CCN3 = list(
    name       = "CCN3 (NovH)",
    chr        = "8",
    pos        = 120435705,
    olink_file = c(
      "GBR_UKB_OLINK_OID20299_CCN3_CCN_family_member_3_adjAgeSexBatPC_InvNorm_22122022.txt.gz",
      "GBR_UKB_Africa_OLINK_OID20299_CCN3_CCN_family_member_3_adjAgeSexBatPC_InvNorm_22122022.txt.gz",
      "GBR_UKB_SAsia_OLINK_OID20299_CCN3_CCN_family_member_3_adjAgeSexBatPC_InvNorm_22122022.txt.gz"
    ),
    olink_label = c("OLINK_BI", "OLINK_AF", "OLINK_SA"),
    olink_N     = c(55000, 8000, 10000)
  ),
  GABARAP = list(
    name       = "GABARAP",
    chr        = "3",
    pos        = 11249929,
    olink_file = c(
      "GBR_UKB_OLINK2_OID31414_GABARAP_Gamma_aminobutyric_acid_receptor_associated_protein_adjAgeSexPC_InvNorm_22122022.txt.gz",
      "GBR_UKB_Africa_OLINK2_OID31414_GABARAP_Gamma_aminobutyric_acid_receptor_associated_protein_adjAgeSexPC_InvNorm_22122022.txt.gz",
      "GBR_UKB_SAsia_OLINK2_OID31414_GABARAP_Gamma_aminobutyric_acid_receptor_associated_protein_adjAgeSexPC_InvNorm_22122022.txt.gz"
    ),
    olink_label = c("OLINK_BI", "OLINK_AF", "OLINK_SA"),
    olink_N     = c(55000, 8000, 10000)
  ),
  GABARAPL1 = list(
    name       = "GABARAPL1",
    chr        = "3",
    pos        = 11249929,
    olink_file = c(
      "GBR_UKB_OLINK2_OID31253_GABARAPL1_Gamma_aminobutyric_acid_receptor_associated_protein_like_1_adjAgeSexPC_InvNorm_22122022.txt.gz",
      "GBR_UKB_Africa_OLINK2_OID31253_GABARAPL1_Gamma_aminobutyric_acid_receptor_associated_protein_like_1_adjAgeSexPC_InvNorm_22122022.txt.gz",
      "GBR_UKB_SAsia_OLINK2_OID31253_GABARAPL1_Gamma_aminobutyric_acid_receptor_associated_protein_like_1_adjAgeSexPC_InvNorm_22122022.txt.gz"
    ),
    olink_label = c("OLINK_BI", "OLINK_AF", "OLINK_SA"),
    olink_N     = c(55000, 8000, 10000)
  )
)

# ---- UKB hg38 GWAS outcomes ONLY ----
# ICBP outcomes (ebi-a-GCST90025981, ebi-a-GCST90025968) removed because they
# use hg19, while Olink data is hg38. Position queries to OpenGWAS for hg19
# data with hg38 coordinates return mismatched/wrong variants.
outcomes <- list(
  HTN  = list(id = "ukb-b-12493", N = 463010),
  MHTN = list(id = "ukb-b-18167", N = 426391)
)

data_dir   <- "D:/多组学MR数据/pqtl"
all_results <- list()

for (pname in names(proteins)) {
  p <- proteins[[pname]]
  for (k in seq_along(p$olink_file)) {
    olink_path <- file.path(data_dir, p$olink_file[k])
    if (!file.exists(olink_path)) {
      cat(sprintf("  [SKIP] File not found: %s\n", p$olink_file[k]))
      next
    }
    label <- p$olink_label[k]
    cat(sprintf("\n=== %s [%s] ===\n", p$name, label))

    # ---- Read and filter Olink (hg38) ----
    cis_start <- p$pos - 200000
    cis_end   <- p$pos + 200000

    dat <- fread(olink_path)
    # chr prefix: Chrom column is "chr2", "chr8", etc.
    dat <- dat[Chrom == sprintf("chr%s", p$chr) &
               Pos >= cis_start & Pos <= cis_end]
    dat <- dat[!is.na(Pval) & Pval > 0]

    cat(sprintf("  pQTL SNPs in cis (+-200kb): %d\n", nrow(dat)))
    if (nrow(dat) < 10) {
      cat("  Too few pQTL SNPs, skip\n")
      next
    }

    # Position key: strip "chr" prefix -> "2:216300185"
    dat$pos_key <- paste0(gsub("^chr", "", dat$Chrom), ":", dat$Pos)

    # ---- Loop over hg38 GWAS outcomes ----
    for (oc_name in names(outcomes)) {
      oc <- outcomes[[oc_name]]
      cis_region <- sprintf("%s:%d-%d", p$chr, cis_start, cis_end)

      # Fetch GWAS variants via OpenGWAS API
      gwas <- tryCatch(
        as.data.frame(associations(variants = cis_region, id = oc$id)),
        error = function(e) {
          cat(sprintf("  %s: API error - %s\n", oc_name, e$message))
          NULL
        }
      )
      if (is.null(gwas) || nrow(gwas) < 10) {
        cat(sprintf("  %s: GWAS query returned %s variants\n",
                    oc_name, if (is.null(gwas)) "NULL" else nrow(gwas)))
        next
      }

      # Normalise GWAS column names (different ieugwasr versions use
      # different names; tolower handles most variation)
      colnames(gwas) <- tolower(colnames(gwas))
      if (!"position" %in% colnames(gwas) && "base_pair_location" %in% colnames(gwas))
        gwas$position <- gwas$base_pair_location
      if (!"p" %in% colnames(gwas) && "pval" %in% colnames(gwas))
        gwas$p <- gwas$pval

      needed <- c("chr", "position", "p")
      if (!all(needed %in% colnames(gwas))) {
        cat(sprintf("  %s: missing columns: %s | have: %s\n",
                    oc_name,
                    paste(setdiff(needed, colnames(gwas)), collapse = ","),
                    paste(colnames(gwas), collapse = ",")))
        next
      }

      # Build GWAS position key (hg38, same as Olink)
      # GWAS chr may be "chr2" or integer 2 - normalise
      gwas$chr_num <- gsub("^chr", "", as.character(gwas$chr))
      gwas$pos_key <- paste0(gwas$chr_num, ":", gwas$position)

      # ---- Position-based matching ----
      common_pos <- intersect(dat$pos_key, gwas$pos_key)
      cat(sprintf("  %s: common positions = %d\n", oc_name, length(common_pos)))

      if (length(common_pos) < 10) {
        cat(sprintf("  %s: SKIP (<10 shared positions)\n", oc_name))
        next
      }

      idx_olink <- match(common_pos, dat$pos_key)
      idx_gwas  <- match(common_pos, gwas$pos_key)
      n <- length(common_pos)

      # ---- Build coloc datasets ----
      # SNP IDs: use Olink Name column (chr:pos:ref:alt) for both datasets
      snp_ids <- dat$Name[idx_olink]

      # MAF for Olink
      olink_maf <- if ("impmaf" %in% tolower(colnames(dat))) {
        mf <- dat$ImpMAF[idx_olink]
        mf[is.na(mf)] <- 0.3
        pmax(pmin(mf, 0.99), 0.01)
      } else {
        rep(0.3, n)
      }

      # MAF for GWAS
      gwas_maf <- if ("eaf" %in% colnames(gwas)) {
        mf <- gwas$eaf[idx_gwas]
        mf[is.na(mf)] <- 0.3
        pmax(pmin(mf, 0.99), 0.01)
      } else {
        rep(0.3, n)
      }

      # ---- Run coloc ----
      # Use type="quant" for Olink (protein levels)
      # Use type="quant" for GWAS too (binary treated as quantitative proxy,
      # which is standard practice when exact case proportion is unknown)
      cres <- coloc.abf(
        dataset1 = list(
          pvalues = dat$Pval[idx_olink],
          N       = p$olink_N[k],
          MAF     = olink_maf,
          type    = "quant",
          snp     = snp_ids
        ),
        dataset2 = list(
          pvalues = gwas$p[idx_gwas],
          N       = oc$N,
          MAF     = gwas_maf,
          type    = "quant",
          snp     = snp_ids
        )
      )

      h0 <- as.numeric(cres$summary["PP.H0.abf"])
      h1 <- as.numeric(cres$summary["PP.H1.abf"])
      h2 <- as.numeric(cres$summary["PP.H2.abf"])
      h3 <- as.numeric(cres$summary["PP.H3.abf"])
      h4 <- as.numeric(cres$summary["PP.H4.abf"])

      cat(sprintf("  -> %s: H0=%.3f H1=%.3f H2=%.3f H3=%.3f H4=%.4f  (%d SNPs)\n",
                  oc_name, h0, h1, h2, h3, h4, n))

      all_results[[length(all_results) + 1]] <- data.frame(
        protein  = p$name,
        platform = label,
        outcome  = oc_name,
        n_snps   = n,
        coloc_H0 = h0,
        coloc_H1 = h1,
        coloc_H2 = h2,
        coloc_H3 = h3,
        coloc_H4 = h4,
        stringsAsFactors = FALSE
      )
    } # end outcomes loop
  } # end platform loop
} # end protein loop

# ---- Save results ----
if (length(all_results) > 0) {
  res <- do.call(rbind, all_results)
  dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
  write.csv(res, "results/tables/coloc_olink_crossplatform.csv", row.names = FALSE)

  cat("\n========================================\n")
  cat("  CROSS-PLATFORM COLOC SUMMARY\n")
  cat("========================================\n")
  cat(sprintf("Total tests       : %d\n", nrow(res)))
  cat(sprintf("H4 > 0.75         : %d\n", sum(res$coloc_H4 > 0.75)))
  cat(sprintf("H4 > 0.50         : %d\n", sum(res$coloc_H4 > 0.50)))
  cat(sprintf("H4 > 0.80 (strong): %d\n", sum(res$coloc_H4 > 0.80)))
  cat("\nResults sorted by H4:\n")
  print(res[order(-res$coloc_H4), ])
} else {
  cat("\n[WARNING] No coloc tests produced results!\n")
}
cat("[Done]\n")
