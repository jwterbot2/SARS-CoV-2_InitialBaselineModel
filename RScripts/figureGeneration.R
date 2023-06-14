#Generates figures of range and means of the number of SNPs in each model for all given combinations of mu, K, and initial bottleneck size.

mu0_k0_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,0,0, ymax = 100);
mu1_k0_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,0,0, ymax = 1000);
mu2_k0_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,0,0, ymax = 10000);
mu0_k1_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,1,0, ymax = 100);
mu1_k1_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,1,0, ymax = 1000);
mu2_k1_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,1,0, ymax = 10000);
mu0_k2_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,2,0, ymax = 100);
mu1_k2_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,2,0, ymax = 1000);
mu2_k2_b0 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,2,0, ymax = 10000);

mu0_k0_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,0,1, ymax = 1500);
mu1_k0_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,0,1, ymax = 2000);
mu2_k0_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,0,1, ymax = 10000);
mu0_k1_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,1,1, ymax = 1500);
mu1_k1_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,1,1, ymax = 2000);
mu2_k1_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,1,1, ymax = 10000);
mu0_k2_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,2,1, ymax = 1500);
mu1_k2_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,2,1, ymax = 2000);
mu2_k2_b1 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,2,1, ymax = 10000);

mu0_k0_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,0,2, ymax = 3000);
mu1_k0_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,0,2, ymax = 3000);
mu2_k0_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,0,2, ymax = 12000);
mu0_k1_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,1,2, ymax = 3000);
mu1_k1_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,1,2, ymax = 3000);
mu2_k1_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,1,2, ymax = 12000);
mu0_k2_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 0,2,2, ymax = 3000);
mu1_k2_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 1,2,2, ymax = 3000);
mu2_k2_b2 = getRequiredPlot_v1(setNames("1000", "_partial_1000"), setNames("30000", "Full"), "filtSNPs", 2,2,2, ymax = 12000);

bottleneck0ArrFigure = grid.arrange(mu0_k0_b0, mu0_k1_b0, mu0_k2_b0, mu1_k0_b0, mu1_k1_b0, mu1_k2_b0, mu2_k0_b0, mu2_k1_b0, mu2_k2_b0, nrow=3, ncol=3)
ggsave("Figure3.tiff", bottleneck0ArrFigure, width = 12, height = 10, units = "in", dpi=1200)
bottleneck0mu0ArrFigure = grid.arrange(mu0_k0_b0, mu0_k1_b0, mu0_k2_b0, nrow=3, ncol=1)
ggsave("Figure3alt.tiff", bottleneck0mu0ArrFigure, width = 11, height = 7.5, units = "in", dpi=600)

bottleneck0ArrFigure = grid.arrange(mu0_k0_b1, mu0_k1_b1, mu0_k2_b1, mu1_k0_b1, mu1_k1_b1, mu1_k2_b1, mu2_k0_b1, mu2_k1_b1, mu2_k2_b1, nrow=3, ncol=3)
ggsave("suppFigure1.tiff", bottleneck0ArrFigure, width = 12, height = 10, units = "in", dpi=1200)

bottleneck0ArrFigure = grid.arrange(mu0_k0_b2, mu0_k1_b2, mu0_k2_b2, mu1_k0_b2, mu1_k1_b2, mu1_k2_b2, mu2_k0_b2, mu2_k1_b2, mu2_k2_b2, nrow=3, ncol=3)
ggsave("suppFigure2.tiff", bottleneck0ArrFigure, width = 12, height = 10, units = "in", dpi=1200)
