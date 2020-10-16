args<-commandArgs(TRUE)

name <- args[1]

title <- paste(name," annotated genome",sep = "")

library(chromoMap)

library(htmltools)

map <- chromoMap('genome.bed', 'cmap.txt', data_based_color_map = T, data_type = "categorical",legend = T, lg_x = 100, lg_y = 300, left_margin = 55, title = title)

save_html(map, 'genome_annotations.html', background = "white", libdir = "lib")

while (!is.null(dev.list()))  dev.off()
