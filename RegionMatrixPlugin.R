library("derfinder")
library("derfinderData")
library("GenomicRanges")
library("knitr")

#####################################################################################
# REGION MATRIX
input <- function(inputfile) {
#pheno <- subset(brainspanPheno, structure_acronym == "AMY")
#p <- pheno[, -which(colnames(pheno) %in% c(
#    "structure_acronym",
#    "structure_name", "file"
#))]
#rownames(p) <- NULL
	print(inputfile)
files <<- rawFiles(inputfile,
    samplepatt = "bw", fileterm = NULL
)
}

run <- function() {}

output <- function(outputfile) {
names(files) <<- gsub(".bw", "", names(files))
fullCov <- fullCoverage(
    files = files, chrs = "chr21",
    totalMapped = rep(1, length(files)), targetSize = 1
)
filteredCov <- lapply(fullCov, filterData, cutoff = 2)
round(c(
    fullCov = object.size(fullCov),
    filteredCov = object.size(filteredCov)
) / 1024^2, 1)
regionMat <- regionMatrix(fullCov, cutoff = 30, L = 76)
write.csv(regionMat$chr21$regions, paste(outputfile, "regions", "csv", sep="."))
write.csv(regionMat$chr21$coverageMatrix, paste(outputfile, "coverage", "csv", sep="."))
saveRDS(regionMat, paste(outputfile, "rds", sep="."))
saveRDS(regionMat$chr21$regions, paste(outputfile, "chr21", "rds", sep="."))
saveRDS(fullCov, paste(outputfile, "fullCov", "rds", sep="."))
saveRDS(filteredCov, paste(outputfile, "filteredCov", "rds", sep="."))
}
#####################################################################################

