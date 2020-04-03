args<-commandArgs(TRUE)

name <- args[1]

title <- paste(name," Engineering Signatures",sep = "")

library(rtracklayer)

library(karyoploteR)

pdf(file="eng_sig_figure.pdf")

custom.genome <- toGRanges("genome.bed")

gff.data <- import("blastn.gff")

pp <- getDefaultPlotParams(plot.type = 1)

pp$leftmargin <- 0.15

pp$topmargin <- 2000

kp <- plotKaryotype(genome = custom.genome, plot.type = 1, plot.params = pp, main = paste0(title))

kpPlotMarkers(kp, data = gff.data,labels=gff.data$ID,ignore.chromosome.ends = T,line.color = "blue",cex=0.7)

dev.off()
