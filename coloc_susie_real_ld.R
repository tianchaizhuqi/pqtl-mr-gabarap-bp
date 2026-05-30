# coloc.susie with real 1000G LD (from pre-extracted cis-region VCF)
# Strategy: read 117MB cis-region VCF → extract genotypes → compute real LD → run coloc.susie
library(data.table)
library(coloc)
library(susieR)
library(ieugwasr)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ0aWFuY2hhaXpodXFpQGdtYWlsLmNvbSIsImlhdCI6MTc3OTc1NDY2MSwiZXhwIjoxNzgwOTY0MjYxfQ.O7wUvLZC46uFLdGOJXZO9_NqKYNRMz9CPNZ24ufdnsHJekgAUykzfcHa7dFITU9yjmlfNO_4DZoLbMJNGKC-stvFBmm6Di420rhgaD7wydKXTsULKJIgXd8NuPnXq0c9v-86aTQqgcYGIeBzAc_KeXYhRxxqevXzBzEm9IOJVwAZek96NM5L_Ez5t6tEm7f7dHArQxDM_VtHiLc8qpS3jXsrWeZqTKA6zsArC_25k60ZEYoUELQeaI2G5Mw0qf9bzMH5SSr_SB7yaFxTf48hc6HYdODmyCj3iULCifJs1Td8_3wMqApJuCkgxkeuImg42yNmZgZqtEfWgcir5NqFYw")

setwd("d:/桌面/测试/dkd_multiomics_mr")

# ---- Step 1: Read cis-pQTL data ----
cat("Loading pQTL cis data...\n")
dat <- fread("data/cis_regions/GABARAP_cis.txt")
dat <- dat[!is.na(rsids) & rsids != "NA" & rsids != "." & rsids != ""]
dat <- dat[Chrom == "chr3" & Pos >= 11049929 & Pos <= 11449929]
cat(sprintf("pQTL SNPs in cis-region: %d\n", nrow(dat)))

# ---- Step 2: Read 1000G VCF cis-region ----
cat("Loading 1000G VCF genotypes...\n")
# Read VCF header to get sample count
vcf_head <- readLines("data/ld_reference/cis_region_1kg.vcf", n = 300)
header_line <- grep("^#CHROM", vcf_head, value = TRUE)
if (length(header_line) == 0) stop("Could not find #CHROM header line in VCF")
header_cols <- strsplit(header_line, "\t")[[1]]
n_samples <- length(header_cols) - 9
cat(sprintf("1000G samples: %d\n", n_samples))

# Read VCF data
vcf <- fread("data/ld_reference/cis_region_1kg.vcf", skip = "#CHROM", header = TRUE)
cat(sprintf("VCF variants in cis-region: %d\n", nrow(vcf)))

# Build position-based IDs for matching (1000G VCF ID column is all ".")
norm_chr <- function(x) gsub("^chr", "", x)
pqtl_pos_id <- paste(norm_chr(dat$Chrom), dat$Pos, sep = ":")
vcf_pos_id   <- paste(norm_chr(vcf[["#CHROM"]]), vcf$POS, sep = ":")
cat(sprintf("pQTL position IDs: %d unique\n", length(unique(pqtl_pos_id))))
cat(sprintf("VCF position IDs: %d unique\n", length(unique(vcf_pos_id))))

# ---- Step 3: Parse genotypes into dosage matrix ----
# Extract all sample columns (col 10 onward)
sample_cols <- 10:ncol(vcf)
geno_raw <- as.matrix(vcf[, ..sample_cols])

# Fast dosage parsing: vectorized over all genotypes at once
cat("Parsing genotypes to dosages...\n")
geno_vec <- as.vector(geno_raw)
gt <- substr(geno_vec, 1, 3)
dosage_vec <- as.integer((substr(gt, 1, 1) == "1") + (substr(gt, 3, 3) == "1"))
geno_dosage <- matrix(dosage_vec, nrow = nrow(geno_raw), ncol = ncol(geno_raw))
cat(sprintf("Genotype matrix: %d variants x %d samples\n", nrow(geno_dosage), ncol(geno_dosage)))

# ---- Step 4: Outcomes ----
outcomes <- list(
  SBP  = list(id = "ebi-a-GCST90025981", N = 422713),
  DBP  = list(id = "ebi-a-GCST90025968", N = 422713),
  HTN  = list(id = "ukb-b-12493",       N = 463010),
  MHTN = list(id = "ukb-b-18167",       N = 426391)
)

all_results <- list()

for (oc_name in names(outcomes)) {
  oc <- outcomes[[oc_name]]
  cat(sprintf("\n========== %s ==========\n", oc_name))

  # Query GWAS
  cis_region <- "3:11049929-11449929"
  gwas <- tryCatch(
    associations(variants = cis_region, id = oc$id),
    error = function(e) { cat(sprintf("  API error: %s\n", e$message)); return(NULL) }
  )
  if (is.null(gwas)) next
  gwas <- as.data.frame(gwas)
  cat(sprintf("  GWAS SNPs: %d\n", nrow(gwas)))

  # Match pQTL ↔ GWAS by rsID (build-independent)
  common_rsid <- intersect(dat$rsids, gwas$rsid)
  cat(sprintf("  pQTL-GWAS common rsIDs: %d\n", length(common_rsid)))
  if (length(common_rsid) < 20) { cat("  Too few, skip\n"); next }

  # Get GWAS positions (hg19) for matched SNPs, then match to VCF (also hg19)
  pqtl_idx_rsid <- match(common_rsid, dat$rsids)
  gwas_idx_rsid <- match(common_rsid, gwas$rsid)
  gwas_pos_id   <- paste(norm_chr(gwas$chr[gwas_idx_rsid]), gwas$position[gwas_idx_rsid], sep = ":")

  # Intersect with VCF positions
  common <- intersect(gwas_pos_id, vcf_pos_id)
  n_common <- length(common)
  cat(sprintf("  Common SNPs (pQTL+G WAS+VCF by position): %d\n", n_common))
  if (n_common < 20) { cat("  Too few, skip\n"); next }

  # Map back to each dataset
  gwas_pos_match <- match(common, gwas_pos_id)
  pqtl_idx <- pqtl_idx_rsid[gwas_pos_match]
  gwas_idx  <- gwas_idx_rsid[gwas_pos_match]
  vcf_idx   <- match(common, vcf_pos_id)

  # Use rsIDs as SNP labels for coloc
  snp_labels <- dat$rsids[pqtl_idx]

  # ---- coloc.abf (standard) ----
  d1_abf <- list(pvalues = dat$Pval[pqtl_idx], N = dat$N[pqtl_idx[1]],
                 MAF = dat$ImpMAF[pqtl_idx], type = "quant", snp = snp_labels)
  d2_abf <- list(pvalues = gwas$p[gwas_idx], N = oc$N,
                 MAF = rep(0.3, n_common), type = "quant", snp = snp_labels)
  cres_abf <- coloc.abf(dataset1 = d1_abf, dataset2 = d2_abf)
  abf_H4 <- as.numeric(cres_abf$summary["PP.H4.abf"])
  abf_H3 <- as.numeric(cres_abf$summary["PP.H3.abf"])
  cat(sprintf("  coloc.abf: H4=%.4f  H3=%.4f\n", abf_H4, abf_H3))

  # ---- Prune to top SNPs for coloc.susie (avoid >200 SNP convergence issues) ----
  susie_common <- common
  susie_snp_labels <- snp_labels
  susie_pqtl_idx <- pqtl_idx
  susie_gwas_idx <- gwas_idx
  susie_vcf_idx <- vcf_idx

  if (n_common > 200) {
    # Keep SNPs with strongest pQTL signal (P < 1e-4) or top 200
    p_snp <- dat$Pval[pqtl_idx]
    keep_susie <- if (sum(p_snp < 1e-4) >= 20) p_snp < 1e-4 else order(p_snp)[1:min(200, n_common)]
    susie_common <- common[keep_susie]
    susie_snp_labels <- snp_labels[keep_susie]
    susie_pqtl_idx <- pqtl_idx[keep_susie]
    susie_gwas_idx <- gwas_idx[keep_susie]
    susie_vcf_idx <- vcf_idx[keep_susie]
    cat(sprintf("  Pruned for coloc.susie: %d -> %d SNPs\n", n_common, length(susie_common)))
  }
  n_susie <- length(susie_common)
  if (n_susie < 20) { cat("  Too few after prune, skip susie\n"); next }

  # ---- Compute real LD from 1000G genotypes ----
  cat(sprintf("  Computing LD for %d SNPs...\n", n_susie))
  geno_sub <- geno_dosage[susie_vcf_idx, , drop = FALSE]
  # Remove monomorphic SNPs (no variation in 1000G)
  var_snp <- apply(geno_sub, 1, var, na.rm = TRUE) > 0
  if (sum(var_snp) < n_susie) {
    cat(sprintf("  Removed %d monomorphic SNPs\n", n_susie - sum(var_snp)))
    keep_idx <- which(var_snp)
    susie_common <- susie_common[keep_idx]
    susie_snp_labels <- susie_snp_labels[keep_idx]
    susie_pqtl_idx <- susie_pqtl_idx[keep_idx]
    susie_gwas_idx  <- susie_gwas_idx[keep_idx]
    geno_sub  <- geno_sub[keep_idx, , drop = FALSE]
    n_susie <- length(susie_common)
    if (n_susie < 20) next
  }

  R_ld <- cor(t(geno_sub), use = "pairwise.complete.obs")
  # Regularize for positive definiteness
  R_ld[is.na(R_ld)] <- 0
  R_ld <- R_ld + diag(0.001, n_susie)
  dimnames(R_ld) <- list(susie_snp_labels, susie_snp_labels)
  cat(sprintf("  LD matrix: %d x %d, mean |r| = %.4f\n",
              n_susie, n_susie, mean(abs(R_ld[upper.tri(R_ld)]))))

  # ---- coloc.susie with real LD ----
  d1_susie <- list(
    beta = dat$Beta[susie_pqtl_idx],
    varbeta = dat$SE[susie_pqtl_idx]^2,
    MAF = dat$ImpMAF[susie_pqtl_idx],
    N = dat$N[susie_pqtl_idx[1]],
    type = "quant",
    snp = susie_snp_labels,
    LD = R_ld
  )
  d2_susie <- list(
    beta = gwas$beta[susie_gwas_idx],
    varbeta = gwas$se[susie_gwas_idx]^2,
    MAF = rep(0.3, n_susie),
    N = oc$N,
    type = "quant",
    snp = susie_snp_labels,
    LD = R_ld
  )

  cat(sprintf("  Running coloc.susie (%d SNPs, real 1000G LD)...\n", n_common))
  cres_susie <- tryCatch(
    coloc.susie(dataset1 = d1_susie, dataset2 = d2_susie,
                p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    error = function(e) {
      cat(sprintf("  coloc.susie ERROR: %s\n", e$message))
      return(NULL)
    }
  )

  susie_H4 <- NA_real_; susie_H3 <- NA_real_
  susie_status <- "OK"
  if (!is.null(cres_susie) && !is.null(cres_susie$summary)) {
    s <- as.data.frame(cres_susie$summary)
    if (nrow(s) > 0) {
      # Take best H4 row (coloc.susie may return multiple CS pairs)
      best_idx <- which.max(s[["PP.H4.abf"]])
      susie_H4 <- as.numeric(s[best_idx, "PP.H4.abf"])
      susie_H3 <- as.numeric(s[best_idx, "PP.H3.abf"])
      if (length(susie_H4) == 0) susie_H4 <- NA_real_
      if (length(susie_H3) == 0) susie_H3 <- NA_real_
      cat(sprintf("  coloc.susie: H4=%.4f  H3=%.4f  (best of %d CS pairs)\n",
                  susie_H4, susie_H3, nrow(s)))
    } else {
      cat("  coloc.susie: summary empty\n")
      susie_status <- "SUMMARY_EMPTY"
    }
  } else {
    cat("  coloc.susie: failed or returned NULL\n")
    susie_status <- "SUSIE_FAILED"
  }

  all_results[[oc_name]] <- data.frame(
    outcome = oc_name,
    n_coloc_abf = n_common,
    n_coloc_susie = n_susie,
    coloc_abf_H4 = abf_H4, coloc_abf_H3 = abf_H3,
    coloc_susie_H4 = susie_H4, coloc_susie_H3 = susie_H3,
    status = susie_status,
    stringsAsFactors = FALSE
  )
}

# ---- Final summary ----
res <- do.call(rbind, all_results)
write.csv(res, "results/tables/coloc_susie_real_ld_results.csv", row.names = FALSE)

cat("\n========== FINAL SUMMARY ==========\n")
cat("LD source: 1000 Genomes Phase 3 EUR, chr3 cis-region genotypes\n")
cat(sprintf("pQTL: deCODE N=%d\n", dat$N[1]))
print(res)
cat("[Done]\n")
