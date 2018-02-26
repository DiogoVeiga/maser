#
# system.file("extdata", file.path("MATS_output", "SE.MATS.ReadsOnTargetAndJunctionCounts.txt"),
#             package = "maser")
#
# path <- system.file("extdata", file.path("MATS_output"),
#                     package = "maser")
# mats <- importEvents(path, c("Hypoxia 0h", "Hypoxia 24h"))
#
# mats.filt <- filterByCoverage(mats, avg_reads = 5)
# dotplot(mats.filt, fdr = 0.05, deltaPSI = 0.1, type = "SE")
# volcano(mats.filt, fdr = 0.05, deltaPSI = 0.1, type = "SE")
#
# splicingDistribution(mats.filt, fdr = 0.05, deltaPSI = 0.1)
#
# PSI_levels(mats.filt, type = "MXE")
#
# pca(mats.filt, type = "SE")
#
# mats.top <- topEvents(mats.filt, fdr = 0.05, deltaPSI = 0.1 )
# dotplot(mats.top, fdr = 0.05, deltaPSI = 0.1, type = "SE")
# volcano(mats.top, fdr = 0.05, deltaPSI = 0.1, type = "SE")
#
# pca(mats.top, type = "SE")
#
