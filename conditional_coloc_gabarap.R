# Conditional coloc: GABARAP vs GABARAPL1 gene assignment (Reviewer C2)
# ===========================================================================
library(data.table)
library(coloc)
library(ieugwasr)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ0aWFuY2hhaXpodXFpQGdtYWlsLmNvbSIsImlhdCI6MTc3OTc1NDY2MSwiZXhwIjoxNzgwOTY0MjYxfQ.O7wUvLZC46uFLdGOJXZO9_NqKYNRMz9CPNZ24ufdnsHJekgAUykzfcHa7dFITU9yjmlfNO_4DZoLbMJNGKC-stvFBmm6Di420rhgaD7wydKXTsULKJIgXd8NuPnXq0c9v-86aTQqgcYGIeBzAc_KeXYhRxxqevXzBzEm9IOJVwAZek96NM5L_Ez5t6tEm7f7dHArQxDM_VtHiLc8qpS3jXsrWeZqTKA6zsArC_25k60ZEYoUELQeaI2G5Mw0qf9bzMH5SSr_SB7yaFxTf48hc6HYdODmyCj3iULCifJs1Td8_3wMqApJuCkgxkeuImg42yNmZgZqtEfWgcir5NqFYw")

setwd("d:/桌面/测试/dkd_multiomics_mr")
smp_dir <- "D:/多组学MR数据/pqtl/final somascan smp"
decode_dir <- "D:/多组学MR数据/pqtl/decode"

# Load pQTL data
gabarap_file <- file.path(smp_dir, "Proteomics_SMP_PC0_5443_62_NPPA_ANP_10032022.txt.gz")
# Actually, GABARAP is NOT in the positive controls list... Let me check what file it is.
# GABARAP = SeqId not in the positive control list. Need to find it.

# Actually the GABARAP pQTL data is in data/cis_regions/
gabarap_cis <- fread("data/cis_regions/GABARAP_cis.txt")
gabarap_cis <- gabarap_cis[!is.na(rsids) & rsids != "NA" & rsids != "." & rsids != ""]
cat(sprintf("GABARAP cis SNPs loaded: %d\n", nrow(gabarap_cis)))

# Load GABARAPL1 full pQTL data from decode folder
cat("Loading GABARAPL1 data...\n")
gabarapl1 <- fread(file.path(decode_dir, "12661_44_GABARAPL1_GBRL1.txt.gz"))
gabarapl1 <- gabarapl1[!is.na(rsids) & rsids != "NA" & rsids != "." & rsids != ""]
cat(sprintf("GABARAPL1 SNPs loaded: %d\n", nrow(gabarapl1)))

# Find lead variants
gabarap_lead <- gabarap_cis[which.min(Pval)]
gabarapl1_lead <- gabarapl1[grep("chr3", Chrom)][which.min(Pval)]
cat(sprintf("GABARAP lead: %s at %s:%d, P=%.2e\n", gabarap_lead$rsids[1], gabarap_lead$Chrom[1], gabarap_lead$Pos[1], gabarap_lead$Pval[1]))
cat(sprintf("GABARAPL1 lead: %s at %s:%d, P=%.2e\n", gabarapl1_lead$rsids[1], gabarapl1_lead$Chrom[1], gabarapl1_lead$Pos[1], gabarapl1_lead$Pval[1]))

# Z-score correlation in cis-region
cis_snps <- gabarap_cis$rsids[gabarap_cis$Pos >= gabarap_lead$Pos[1] - 200000 & gabarap_cis$Pos <= gabarap_lead$Pos[1] + 200000]
gabarapl1_cis <- gabarapl1[rsids %in% cis_snps & Chrom == "chr3"]
common_gene <- intersect(gabarap_cis$rsids, gabarapl1_cis$rsids)
cat(sprintf("Common SNPs between GABARAP and GABARAPL1 in cis-region: %d\n", length(common_gene)))

if(length(common_gene) > 20) {
  gabarap_z <- gabarap_cis$Beta[match(common_gene, gabarap_cis$rsids)] / gabarap_cis$SE[match(common_gene, gabarap_cis$rsids)]
  gabarapl1_z <- gabarapl1_cis$Beta[match(common_gene, gabarapl1_cis$rsids)] / gabarapl1_cis$SE[match(common_gene, gabarapl1_cis$rsids)]
  z_cor <- cor(gabarap_z, gabarapl1_z, use="complete.obs")
  cat(sprintf("Z-score correlation: r=%.4f\n", z_cor))
}

# Conditional coloc for each BP GWAS
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

  cis_region <- sprintf("3:%d-%d", gabarap_lead$Pos[1] - 200000, gabarap_lead$Pos[1] + 200000)
  gwas <- tryCatch(associations(variants=cis_region, id=oc$id), error=function(e) NULL)
  if(is.null(gwas)) { cat("API error\n"); next }
  gwas <- as.data.frame(gwas)

  # GABARAP coloc
  common_ga <- intersect(gabarap_cis$rsids, gwas$rsid)
  if(length(common_ga) >= 20) {
    idx_p <- match(common_ga, gabarap_cis$rsids)
    idx_g <- match(common_ga, gwas$rsid)
    d1 <- list(pvalues=gabarap_cis$Pval[idx_p], N=gabarap_cis$N[idx_p[1]], MAF=gabarap_cis$ImpMAF[idx_p], type="quant", snp=common_ga)
    d2 <- list(pvalues=gwas$p[idx_g], N=oc$N, MAF=rep(0.3, length(common_ga)), type="quant", snp=common_ga)
    cres_ga <- coloc.abf(dataset1=d1, dataset2=d2)
    ga_H4 <- as.numeric(cres_ga$summary["PP.H4.abf"])
    ga_H3 <- as.numeric(cres_ga$summary["PP.H3.abf"])
    cat(sprintf("  GABARAP:    H4=%.4f H3=%.4f\n", ga_H4, ga_H3))
  } else { ga_H4 <- NA; ga_H3 <- NA }

  # GABARAPL1 coloc
  common_gl1 <- intersect(gabarapl1$rsids, gwas$rsid)
  if(length(common_gl1) >= 20) {
    idx_p <- match(common_gl1, gabarapl1$rsids)
    idx_g <- match(common_gl1, gwas$rsid)
    d1 <- list(pvalues=gabarapl1$Pval[idx_p], N=gabarapl1$N[idx_p[1]], MAF=gabarapl1$ImpMAF[idx_p], type="quant", snp=common_gl1)
    d2 <- list(pvalues=gwas$p[idx_g], N=oc$N, MAF=rep(0.3, length(common_gl1)), type="quant", snp=common_gl1)
    cres_gl1 <- coloc.abf(dataset1=d1, dataset2=d2)
    gl1_H4 <- as.numeric(cres_gl1$summary["PP.H4.abf"])
    gl1_H3 <- as.numeric(cres_gl1$summary["PP.H3.abf"])
    cat(sprintf("  GABARAPL1:  H4=%.4f H3=%.4f\n", gl1_H4, gl1_H3))
  } else { gl1_H4 <- NA; gl1_H3 <- NA }

  # Conditional coloc: condition GABARAPL1 on GABARAP lead variant
  # Using coloc.signals to check if GABARAPL1 signal is independent
  if(length(common_gl1) >= 20) {
    idx_p <- match(common_gl1, gabarapl1$rsids)
    idx_g <- match(common_gl1, gwas$rsid)
    d1_cond <- list(pvalues=gabarapl1$Pval[idx_p], N=gabarapl1$N[idx_p[1]],
                    MAF=gabarapl1$ImpMAF[idx_p], type="quant", snp=common_gl1)
    d2_cond <- list(pvalues=gwas$p[idx_g], N=oc$N, MAF=rep(0.3, length(common_gl1)),
                    type="quant", snp=common_gl1)
    sig <- tryCatch(
      coloc.signals(dataset1=d1_cond, dataset2=d2_cond),
      error=function(e) { cat(sprintf("  coloc.signals error: %s\n", e$message)); return(NULL) }
    )
    if(!is.null(sig)) {
      cat(sprintf("  coloc.signals (GABARAPL1): %d independent signal(s)\n", nrow(sig)))
      for(i in 1:nrow(sig)) {
        cat(sprintf("    Signal %d: %s, H4=%.4f\n", i, sig$snp[i], sig$PP.H4.abf[i]))
      }
    }
  }

  results_all[[length(results_all)+1]] <- data.frame(
    outcome=oc_name,
    gabarap_H4=ga_H4, gabarap_H3=ga_H3,
    gabarapl1_H4=gl1_H4, gabarapl1_H3=gl1_H3,
    stringsAsFactors=FALSE
  )
}

res <- do.call(rbind, results_all)
write.csv(res, "results/tables/conditional_coloc_gabarap.csv", row.names=FALSE)
cat("\n========== SUMMARY ==========\n")
print(res)
cat("[Done]\n")
