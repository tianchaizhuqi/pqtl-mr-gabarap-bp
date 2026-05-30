# Nature-quality figures for Communications Biology
# Requires: ggplot2, viridis, ggsci, RColorBrewer

library(ggplot2)
library(viridis)      # perceptually uniform, colorblind-safe
library(ggsci)        # Nature Publishing Group palettes
library(RColorBrewer) # diverging/sequential palettes

setwd("d:/桌面/测试/dkd_multiomics_mr")

# ---- Figure 4: Positive Control Calibration Heatmap ----
# Color scheme: viridis "rocket" or scico "vik" for diverging H4 values
# Best for heatmap: use a diverging palette centered at H4=0.5
# RColorBrewer "RdBu" reversed: low=blue, mid=white, high=red

pos_data <- data.frame(
  Gene = factor(c('NPPA','NPPA','NPPA','NPPA',
                   'AGT','AGT','AGT','AGT',
                   'ACE2','ACE2','ACE2','ACE2',
                   'REN','REN','REN','REN',
                   'ACE','ACE','ACE','ACE',
                   'NPPB','NPPB','NPPB','NPPB'),
                levels = rev(c('NPPA','NPPB','AGT','REN','ACE','ACE2'))),
  Outcome = rep(c('SBP','DBP','HTN','MHTN'), 6),
  H4 = c(0.936,0.088,0.002,0.009,
         0.001,0.001,0.679,0.027,
         0.116,0.534,0.006,0.004,
         0.001,0.001,0.004,0.002,
         0.002,0.006,0.003,0.002,
         0.001,0.001,0.001,0.001)
)
pos_data$Label <- sprintf('%.3f', pos_data$H4)
pos_data$Font <- ifelse(pos_data$H4 > 0.5, 'bold', 'plain')

# ---- Palette A: RColorBrewer RdBu (diverging, classic) ----
pdf('submission/figures/Figure4_PositiveControls_A.pdf', width=8, height=5)
ggplot(pos_data, aes(x=Outcome, y=Gene, fill=H4)) +
  geom_tile(color='white', linewidth=0.5) +
  geom_text(aes(label=Label, fontface=Font), size=3.8) +
  scale_fill_distiller(palette='RdBu', direction=1, limits=c(0,1),
                       name=expression(PP.H[4])) +
  labs(title='Figure 4. Positive control colocalization calibration',
       subtitle=expression('6 proteins × 4 GWAS; bold = ' * PP.H[4] > 0.5 * ' (suggestive)'),
       x='', y='') +
  theme_minimal(base_size=14) +
  theme(panel.grid=element_blank(),
        plot.title=element_text(face='bold'),
        legend.position='right')
dev.off()

# ---- Palette B: viridis magma (warm, modern, colorblind-safe) ----
pdf('submission/figures/Figure4_PositiveControls_B.pdf', width=8, height=5)
ggplot(pos_data, aes(x=Outcome, y=Gene, fill=H4)) +
  geom_tile(color='white', linewidth=0.5) +
  geom_text(aes(label=Label, fontface=Font), size=3.8) +
  scale_fill_viridis(option='magma', limits=c(0,1), direction=1,
                     name=expression(PP.H[4])) +
  labs(title='Figure 4. Positive control colocalization calibration (viridis magma)',
       subtitle=expression('6 proteins × 4 GWAS; bold = ' * PP.H[4] > 0.5),
       x='', y='') +
  theme_minimal(base_size=14) +
  theme(panel.grid=element_blank(),
        plot.title=element_text(face='bold'),
        legend.position='right')
dev.off()

# ---- Palette C: scico vik (diverging, perceptually uniform) ----
# Requires: install.packages('scico')
message("Palette C (scico::vik) requires scico package. Install with install.packages('scico')")
# library(scico)
# ggplot(...) + scale_fill_scico(palette='vik', limits=c(0,1), ...)

cat("Generated Figure 4 variants: A (RColorBrewer RdBu), B (viridis magma)\n")
cat("Both saved to submission/figures/\n")
cat("\nFor Nature submission, recommend Palette A (RdBu) — classic, familiar to reviewers.\n")
cat("Palette B (magma) is more modern and completely colorblind-safe.\n")
