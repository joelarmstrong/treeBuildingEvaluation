#!/usr/bin/env Rscript
# Plot a combined coalescence/coverage plot for multiple different runs to "combined.pdf".
# Usage: plotCombinedStats.R <coalescence results, comma-separated> <coverage results, comma-separated> [(optional) names, comma-separated]
require(gridExtra)
require(stringr)
require(plyr)
source("plotCoverageStats.R")
source("plotCoalescenceStats.R")

# parse args
args <- commandArgs(TRUE)
coalescences <- str_split(args[1], ",")
coverages <- str_split(args[2], ",")
stopifnot(length(coalescences) == length(coverages))
if (length(args) > 2) {
    names <- str_split(args[3], ",")
    stopifnot(length(coalescences) == length(names))
} else {
    names <- sapply(coalescences, basename)
}

# get the coalescence data frames
coalescenceDfs <- list()
tmp <- data.frame(path=unlist(coalescences), name=names, stringsAsFactors=FALSE)
for (i in 1:nrow(tmp)) {
    row <- tmp[i,]
    coalescenceDfs[[i]] <- getCoalescenceResults(row$path, row$name)
}
coalescenceDf <- rbind.fill(coalescenceDfs)

# get the coverage data frames
coverageDfs <- list()
tmp <- data.frame(path=unlist(coverages), name=names, stringsAsFactors=FALSE)
print(tmp)
for (i in 1:nrow(tmp)) {
    row <- tmp[i,]
    print(row)
    coverageDfs[[i]] <- getCoverageDataFrame(row$path, row$name)
}
print(coverageDfs)
coverageDf <- rbind.fill(coverageDfs)

print(coverageDf)
coalescencePlot <- getCoalescencePlot(coalescenceDf, "human")
coveragePlot <- getCoveragePlot(coverageDf)

pdf("combined.pdf")
print(grid.arrange(coveragePlot, coalescencePlot))
dev.off()
