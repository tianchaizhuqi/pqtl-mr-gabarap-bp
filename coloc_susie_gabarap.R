# coloc sensitivity: signals + conditional coloc for GABARAP x 4 BP GWAS
# Addresses reviewer C1 (single-causal-variant assumption violation)
library(data.table)
library(coloc)
library(ieugwasr)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ0aWFuY2hhaXpodXFpQGdtYWlsLmNvbSIsImlhdCI6MTc3OTc1NDY2MSwiZXhwIjoxNzgwOTY0MjYxfQ.O7wUvLZC46uFLdGOJXZO9_NqKYNRMz9CPNZ24ufdnsHJekgAUykzfcHa7dFITU9yjmlfNO_4DZoLbMJNGKC-stvFBmm6Di420rhgaD7wydKXTsULKJIgXd8NuPnXq0c9v-86aTQqgcYGIeBzAc_KeXYhRxxqevXzBzEm9IOJVwAZek96NM5L_Ez5t6tEm7f7dHArQxDM_VtHiLc8qpS3jXsrWeZqTKA6zsArC_25k60ZEYoUELQeaI2G5Mw0qf9bzMH5SSr_SB7yaFxTf48hc6HYdODmyCj3iULCifJs1Td8_3wMqApJuCkgxkeuImg42yNmZgZqtEfWgcir5NqFYw")

setwd("d:/桌面/测试/dkd_multiomics_mr")

# Load GABARAP cis-pQTL
dat <- fread("data/cis_regions/GABARAP_cis.txt")
dat <- dat[!is.na(rsids) & rsids != "NA" & rsids != "." & rsids != ""]
dat <- dat[Chrom == "chr3" & Pos >= 11049929 & Pos <= 11449929]
cat(sprintf("GABARAP cis-pQTL SNPs: %d\n", nrow(dat)))

outcomes <- list(
  SBP  = list(id="ebi-a-GCST90025981", N=422713),
  DBP  = list(id="ebi-a-GCST90025968", N=422713),
  HTN  = list(id="ukb-b-12493",       N=463010),
  MHTN = list(id="ukb-b-18167",       N=426391)
)

results_all <- list()

for(oc_name in names(outcomes)) {
  oc <- outcomes[[oc_name]]
  cat(sprintf("\n========== %s ==========\n", oc_name))

  cis_region <- "3:11049929-11449929"
  gwas <- tryCatch(associations(variants=cis_region, id=oc$id), error=function(e) NULL)
  if(is.null(gwas)) { cat("API error\n"); next }
  gwas <- as.data.frame(gwas)

  common <- intersect(dat$rsids, gwas$rsid)
  n_common <- length(common)
  cat(sprintf("Common SNPs: %d\n", n_common))
  if(n_common < 20) next

  pqtl_idx <- match(common, dat$rsids)
  gwas_idx <- match(common, gwas$rsid)

  d1 <- list(pvalues=dat$Pval[pqtl_idx], N=dat$N[pqtl_idx[1]],
             MAF=dat$ImpMAF[pqtl_idx], type="quant", snp=common)
  d2 <- list(pvalues=gwas$p[gwas_idx], N=oc$N,
             MAF=rep(0.3, n_common), type="quant", snp=common)

  # Standard coloc
  cres <- coloc.abf(dataset1=d1, dataset2=d2)
  s <- cres$summary
  pph4 <- as.numeric(s["PP.H4.abf"])
  pph3 <- as.numeric(s["PP.H3.abf"])
  cat(sprintf("coloc.abf: H4=%.4f H3=%.4f\n", pph4, pph3))

  # coloc.signals: test for multiple independent signals
  # This conditions on the top signal and tests for additional signals
  cat("coloc.signals (conditional analysis):\n")
  sig <- tryCatch(
    coloc.signals(dataset1=d1, dataset2=d2),
    error=function(e) { cat(sprintf("Error: %s\n", e$message)); return(NULL) }
  )

  if(!is.null(sig)) {
    n_signals <- nrow(sig)
    cat(sprintf("  Independent signals: %d\n", n_signals))
    for(i in 1:n_signals) {
      cat(sprintf("  Signal %d: SNP=%s, H4=%.4f\n", i, sig$snp[i], sig$PP.H4.abf[i]))
    }

    results_all[[length(results_all)+1]] <- data.frame(
      outcome=oc_name, n_common=n_common,
      coloc_abf_H4=pph4, coloc_abf_H3=pph3,
      n_independent_signals=n_signals,
      top_signal_snp=sig$snp[1],
      top_signal_H4=sig$PP.H4.abf[1],
      stringsAsFactors=FALSE
    )
  } else {
    results_all[[length(results_all)+1]] <- data.frame(
      outcome=oc_name, n_common=n_common,
      coloc_abf_H4=pph4, coloc_abf_H3=pph3,
      n_independent_signals=NA, top_signal_snp=NA, top_signal_H4=NA,
      stringsAsFactors=FALSE
    )
  }
}

res <- do.call(rbind, results_all)
write.csv(res, "results/tables/coloc_signals_sensitivity.csv", row.names=FALSE)
cat("\n========== SUMMARY ==========\n")
print(res)
cat("[Done]\n")
