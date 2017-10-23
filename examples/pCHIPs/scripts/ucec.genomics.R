#!/usr/bin/env	Rscript

source("lib/gistic.R")

gisticObj <- parse.gistic("gdac.broadinstitute.org_UCEC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt", 
	"gdac.broadinstitute.org_UCEC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/amp_genes.conf_99.txt",
	"gdac.broadinstitute.org_UCEC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/del_genes.conf_99.txt")

# sample groups
serious.samples <- as.character(read.table("ucec/ucec-serious.txt", header=F)[,1])
better.samples <- as.character(read.table("ucec/better.surv.samples.txt", header=F)[,1])


# write out CNV events
cnv.events <- proportions(gisticObj, serious.samples, better.samples)
write.table(sig.dels.HUGO, file="ucec/DIFF_A.cnv.dels.txt", quote=F, row.names=F)
write.table(sig.amps.HUGO, file="ucec/DIFF_A.cnv.amps.txt", quote=F, row.names=F)

source("lib/chasm.R")

mutations <- read.table("ucec/mutations.tab", header=T, sep='\t', row.names=1, check.names=F)
colnames(mutations) <- unlist(lapply(strsplit(colnames(data), ''), function(x) paste(x[1:15], collapse='')))


sig.muts <- proportions(mutations, serious.samples, better.samples, correct=F)
write.table(sig.muts, file="ucec/DIFF_A.muts.txt", quote=F, row.names=F)

