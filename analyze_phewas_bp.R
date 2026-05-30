library(data.table)
phewas <- readRDS('results/phewas_scan.rds')

# All BP-related trait patterns
bp_patterns <- 'blood pressure|hypertension|systolic|diastolic|SBP|DBP|HTN|BP'
bp_idx <- grep(bp_patterns, phewas$outcome_trait, ignore.case = TRUE)
bp <- phewas[bp_idx, ]

cat(sprintf('BP-related: %d unique outcomes, %d total rows\n',
            length(unique(bp$outcome_id)), nrow(bp)))
cat(sprintf('Unique genes with BP association: %d\n', length(unique(bp$gene))))

# Best P per gene
gene_best <- aggregate(outcome_pval ~ gene, data = bp, FUN = min)
gene_best <- gene_best[order(gene_best$outcome_pval), ]

cat('\n=== Top 50 proteins by best BP P-value ===\n')
print(gene_best[1:min(50, nrow(gene_best)), ], row.names = FALSE)

# Save full results
write.csv(gene_best, 'results/tables/phewas_bp_gene_ranking.csv', row.names = FALSE)
cat(sprintf('\nSaved %d genes to results/tables/phewas_bp_gene_ranking.csv\n', nrow(gene_best)))

# How many have P < 0.05? P < 0.01? P < 1e-4?
cat(sprintf('P < 0.05: %d genes\n', sum(gene_best$outcome_pval < 0.05)))
cat(sprintf('P < 0.01: %d genes\n', sum(gene_best$outcome_pval < 0.01)))
cat(sprintf('P < 1e-4: %d genes\n', sum(gene_best$outcome_pval < 1e-4)))
cat(sprintf('Bonferroni (P < 1.96e-6): %d genes\n', sum(gene_best$outcome_pval < 1.96e-6)))
