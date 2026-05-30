# coloc.susie FN1.3 (Fibronectin) — DBP, chr2:216100185-216500185
library(data.table); library(coloc); library(susieR); library(ieugwasr)
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ0aWFuY2hhaXpodXFpQGdtYWlsLmNvbSIsImlhdCI6MTc3OTc1NDY2MSwiZXhwIjoxNzgwOTY0MjYxfQ.O7wUvLZC46uFLdGOJXZO9_NqKYNRMz9CPNZ24ufdnsHJekgAUykzfcHa7dFITU9yjmlfNO_4DZoLbMJNGKC-stvFBmm6Di420rhgaD7wydKXTsULKJIgXd8NuPnXq0c9v-86aTQqgcYGIeBzAc_KeXYhRxxqevXzBzEm9IOJVwAZek96NM5L_Ez5t6tEm7f7dHArQxDM_VtHiLc8qpS3jXsrWeZqTKA6zsArC_25k60ZEYoUELQeaI2G5Mw0qf9bzMH5SSr_SB7yaFxTf48hc6HYdODmyCj3iULCifJs1Td8_3wMqApJuCkgxkeuImg42yNmZgZqtEfWgcir5NqFYw")
setwd("d:/桌面/测试/dkd_multiomics_mr")

chr <- "2"; cis_start <- 216100185; cis_end <- 216500185
outcome_id <- "ebi-a-GCST90025968"; outcome_N <- 422713; oc_name <- "DBP"
decode_paths <- c("D:/多组学MR数据/pqtl/decode/Proteomics_SMP_PC0_3434_34_FN1_FN1_3_10032022.txt.gz",
                  "data/pqtl/decode/Proteomics_SMP_PC0_3434_34_FN1_FN1_3_10032022.txt.gz")
vcf_base <- "D:/多组学MR数据/ld_reference"
vcf_chr  <- file.path(vcf_base, "ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

cat("=== FN1.3 (Fibronectin) ===\n")

# Read deCODE data
decode_path <- decode_paths[file.exists(decode_paths)][1]
if (is.na(decode_path)) stop("File not found")
dat <- fread(decode_path)
dat <- dat[!is.na(rsids) & rsids != "NA" & rsids != "." & rsids != ""]
dat <- dat[Chrom == sprintf("chr%s", chr) & Pos >= cis_start & Pos <= cis_end]
cat(sprintf("pQTL SNPs: %d\n", nrow(dat)))

# Extract VCF
out_vcf <- file.path(vcf_base, sprintf("cis_chr%s_%d_%d.vcf", chr, cis_start, cis_end))
if (!file.exists(out_vcf)) {
  cmd <- sprintf('gzip -cd "%s" | awk \'/^#/ || ($1=="%s" && $2>=%d && $2<=%d)\' > "%s"',
                 vcf_chr, chr, cis_start, cis_end, out_vcf)
  system(cmd)
}
vcf_head <- readLines(out_vcf, n = 300)
hl <- grep("^#CHROM", vcf_head, value = TRUE)
n_samples <- length(strsplit(hl, "\t")[[1]]) - 9
vcf <- fread(out_vcf, skip = "#CHROM", header = TRUE)
cat(sprintf("VCF: %d variants, %d samples\n", nrow(vcf), n_samples))

# GWAS
gwas <- as.data.frame(associations(variants = sprintf("%s:%d-%d", chr, cis_start, cis_end), id = outcome_id))
cat(sprintf("GWAS: %d SNPs\n", nrow(gwas)))

# Match
norm_chr <- function(x) gsub("^chr", "", x)
common_rsid <- intersect(dat$rsids, gwas$rsid)
cat(sprintf("Common rsID: %d\n", length(common_rsid)))
pi <- match(common_rsid, dat$rsids); gi <- match(common_rsid, gwas$rsid)
gw_pos <- paste(norm_chr(gwas$chr[gi]), gwas$position[gi], sep=":")
vc_pos <- paste(norm_chr(vcf[["#CHROM"]]), vcf$POS, sep=":")
common <- intersect(gw_pos, vc_pos)
cat(sprintf("Common (3-way): %d\n", length(common)))
gm <- match(common, gw_pos); vm <- match(common, vc_pos)
pqtl_idx <- pi[gm]; gwas_idx <- gi[gm]; vcf_idx <- vm
snp_labels <- dat$rsids[pqtl_idx]; n <- length(common)

# coloc.abf
cres <- coloc.abf(
  dataset1 = list(pvalues = dat$Pval[pqtl_idx], N = dat$N[1], MAF = dat$ImpMAF[pqtl_idx], type = "quant", snp = snp_labels),
  dataset2 = list(pvalues = gwas$p[gwas_idx], N = outcome_N, MAF = rep(0.3, n), type = "quant", snp = snp_labels))
cat(sprintf("coloc.abf: H4=%.4f H3=%.4f\n", cres$summary["PP.H4.abf"], cres$summary["PP.H3.abf"]))

# Real LD
susie_labels <- snp_labels; susie_pi <- pqtl_idx; susie_gi <- gwas_idx; susie_vi <- vcf_idx
if (n > 200) {
  p_snp <- dat$Pval[pqtl_idx]
  keep <- if (sum(p_snp < 1e-4) >= 20) p_snp < 1e-4 else order(p_snp)[1:200]
  susie_labels <- snp_labels[keep]; susie_pi <- pqtl_idx[keep]
  susie_gi <- gwas_idx[keep]; susie_vi <- vcf_idx[keep]
  cat(sprintf("Pruned: %d -> %d\n", n, length(susie_labels)))
}
ns <- length(susie_labels)
sample_cols <- 10:ncol(vcf)
geno_raw <- as.matrix(vcf[susie_vi, ..sample_cols])
gt <- substr(as.vector(geno_raw), 1, 3)
dosage <- matrix(as.integer((substr(gt,1,1)=="1") + (substr(gt,3,3)=="1")), nrow=ns, ncol=n_samples)
var_snp <- apply(dosage, 1, var, na.rm=TRUE) > 0
if (sum(var_snp) < ns) {
  dosage <- dosage[var_snp,,drop=FALSE]; susie_labels <- susie_labels[var_snp]
  susie_pi <- susie_pi[var_snp]; susie_gi <- susie_gi[var_snp]; ns <- sum(var_snp)
}
R <- cor(t(dosage), use="pairwise.complete.obs"); R[is.na(R)] <- 0; R <- R + diag(0.001, ns)
dimnames(R) <- list(susie_labels, susie_labels)
cat(sprintf("Real LD: %d x %d, |r|_mean=%.4f\n", ns, ns, mean(abs(R[upper.tri(R)]))))

# coloc.susie
d1 <- list(beta=dat$Beta[susie_pi], varbeta=dat$SE[susie_pi]^2, MAF=dat$ImpMAF[susie_pi],
           N=dat$N[1], type="quant", snp=susie_labels, LD=R)
d2 <- list(beta=gwas$beta[susie_gi], varbeta=gwas$se[susie_gi]^2, MAF=rep(0.3,ns),
           N=outcome_N, type="quant", snp=susie_labels, LD=R)
cs <- coloc.susie(dataset1=d1, dataset2=d2, p1=1e-4, p2=1e-4, p12=1e-5)
s <- as.data.frame(cs$summary)
best <- which.max(s[["PP.H4.abf"]])
cat(sprintf("coloc.susie: H4=%.4f H3=%.4f (%d pairs)\n", s[best,"PP.H4.abf"], s[best,"PP.H3.abf"], nrow(s)))

write.csv(data.frame(
  protein="FN1.3 (Fibronectin)", outcome=oc_name, n_abf=n, n_susie=ns,
  abf_H4=as.numeric(cres$summary["PP.H4.abf"]), abf_H3=as.numeric(cres$summary["PP.H3.abf"]),
  susie_H4=as.numeric(s[best,"PP.H4.abf"]), susie_H3=as.numeric(s[best,"PP.H3.abf"]),
  ld_mean_r=mean(abs(R[upper.tri(R)])), stringsAsFactors=FALSE),
  "results/tables/coloc_susie_FN1_vcf.csv", row.names=FALSE)
cat("[Done]\n")
