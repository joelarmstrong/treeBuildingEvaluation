#!/usr/bin/env Rscript
# Plot a combined coalescence/coverage plot for multiple different runs to "combined.pdf".
# Usage: plotCombinedStats.R <coalescence results, comma-separated> <coverage results, comma-separated> [(optional) names, comma-separated]
require(gridExtra)
require(stringr)
require(plyr)
require(ggplot2)
require(XML)
require(reshape2)

# Data frame: genome1 genome2 identical early late identicalFraction earlyFraction lateFraction name
# Aggregate results are indicated by genome2="aggregate"
# Total aggregate results are genome1="aggregate", genome2="aggregate"
getCoalescenceResults <- function(xmlPath, name) {
    xmlFile <- xmlInternalTreeParse(xmlPath)
    # Genome-to-genome results
    genomeResults <- data.frame(t(sapply(xpathApply(xmlFile, "//genomeCoalescenceTest/coalescenceResults", xmlAttrs), unlist)))

    # Aggregate results per-genome
    genomes <- unlist(xpathApply(xmlFile, "//genomeCoalescenceTest", xmlGetAttr, "genome"))
    genomeAggregateResults <- data.frame(t(sapply(xpathApply(xmlFile, "//genomeCoalescenceTest/aggregateCoalescenceResults", xmlAttrs), unlist)))
    genomeAggregateResults$genome1 <- genomes
    genomeAggregateResults$genome2 <- rep("aggregate", length(genomeAggregateResults$genome1))

    # Total aggregate results
    aggregateResults <- data.frame(t(sapply(xpathApply(xmlFile, "/coalescenceTest/aggregateCoalescenceResults", xmlAttrs), unlist)))
    aggregateResults$genome1 <- rep("aggregate", length(aggregateResults$identical))
    aggregateResults$genome2 <- rep("aggregate", length(aggregateResults$identical))

    combined <- rbind(genomeResults, genomeAggregateResults, aggregateResults)
    combined$name <- rep(name, length(combined$genome1))
    # reorder name factors so they will appear in the order given
    combined$name <- factor(combined$name, as.character(combined$name))
    combined <- transform(combined, identical=as.numeric(as.character(identical)), late=as.numeric(as.character(late)), early=as.numeric(as.character(early)),
              identicalFraction=as.numeric(as.character(identicalFraction)), lateFraction=as.numeric(as.character(lateFraction)), earlyFraction=as.numeric(as.character(earlyFraction)))
    return(combined)
}

# get a base ggplot2 plot comparing multiple coalescence results (in a
# single data frame) against the same reference.
getCoalescencePlot <- function(coalescenceResults, refGenome) {
    meltedCoalResults <- melt(coalescenceResults, id.vars=c("genome1", "genome2", "name"), measure.vars=c("identicalFraction", "lateFraction", "earlyFraction"))
    return(ggplot(subset(meltedCoalResults, genome1 == refGenome & genome2 != "aggregate"), aes(y=value, x=name, fill=variable)) + geom_bar(stat="identity") + facet_grid(~ genome2) + theme_classic())
}

# Get data frame with columns: genome, coverageFraction,
# multiCoverageFraction, sitesMapping1Times, sitesMapping2Times, ...
getCoverageDataFrame <- function(path, name) {
    df <- read.table(path, sep=",", header=T)
    numRefSites <- max(df$sitesMapping1Times)
    df$coverageFraction <- df$sitesMapping1Times / numRefSites
    df$multiCoverageFraction <- df$sitesMapping2Times / numRefSites
    df$name <- rep(name, length(df$Genome))
    # reorder name factors so they will appear in the order given
    df$name <- factor(df$name, as.character(df$name))
    return(df)
}

# get a base ggplot2 plot comparing multiple coverage results (in a
# single data frame).
getCoveragePlot <- function(coverageResults) {
    return(ggplot(coverageResults, aes(x=name, y=coverageFraction)) + geom_bar(stat="identity") + facet_grid( ~ Genome)  + theme_classic())
}

# parse args
args <- commandArgs(TRUE)
coalescences <- str_split(args[1], ",")
coverages <- str_split(args[2], ",")
stopifnot(length(coalescences) == length(coverages))
if (length(args) > 2) {
    names <- unlist(str_split(args[3], ","))
} else {
    names <- unlist(coalescences)
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
coalescencePlot <- getCoalescencePlot(coalescenceDf, "human") + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + theme(legend.position="top")
coveragePlot <- getCoveragePlot(coverageDf) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

pdf("combined.pdf", width=10, height=5)
print(grid.arrange(coveragePlot, coalescencePlot, ncol=2))
dev.off()
