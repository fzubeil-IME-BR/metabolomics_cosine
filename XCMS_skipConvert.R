#Script for generating Heatmaps based on metabolomic cosine similarity
# Author: Florian Zubeil, 2017

calc_cosine_path <- "calc_cosine.exe" 	#Set path to calc_cosine.exe here
#msconvert_path <- "msconvert.exe" 		#Set path to msconvert here. Skipped in this version
filter_pattern <- "**"
num_cores <- 12

#Load libraries
library(xcms)
library(BiocParallel)
library("d3heatmap")
library(htmlwidgets)
library(scales)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


#batch convert *.d files to mzXML using Proteowizard msconvert.exe
#files <- list.files(path=".", pattern=filter_pattern, full.names=T, recursive=F)
#system.time(lapply(files, function(file) {
#	msconvert <- msconvert_path
#	system2(msconvert, c(file, "--mzXML", "-o mzXML\\", "--zlib", "--filter \"peakPicking vendor msLevel=1-1\"", "--filter \"threshold absolute 10000 most-intense\""))
#}))

#perform peak picking etc. using xcms
setwd("mzXML")
xset <- xcmsSet(method="centWave",ppm=15,peakwidth=c(5,20),nSlaves=num_cores)
setwd("..")
pdf("AA_retcor.pdf")
xset1 <- retcor(xset,method="obiwarp",plottype=c("deviation"))
dev.off()
xset2 <- group(xset1, bw = 2, minfrac = 0.5, mzwid = 0.015)
xset3 <- fillPeaks(xset2)
pdf("AA_QC.pdf",paper="a4r")
par(mfrow=c(2,3))
plotQC(xset3)
dev.off()
pt <- peakTable(xset3, filebase="peakList")
save.image()

#run calc_cosine
system2(command = calc_cosine_path, args = c("peakList.tsv","0"))

#build heatmap
data = read.csv("peakListCosine.csv", header = TRUE, row.names = 1, sep=";")
matrix <- as.matrix(data)
dim(matrix)
map <- d3heatmap(matrix, symm = TRUE,  distfun=function(x) as.dist((1-x)/2), colors=col_numeric("RdYlBu", domain = c(0,1)))
saveWidget(map, "heatmap.html", selfcontained=FALSE)
print(map)
