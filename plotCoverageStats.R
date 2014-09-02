# Parse and plot stats from halCoverage
require(plyr)
require(ggplot2)

# Get data frame with columns: genome, coverageFraction,
# multiCoverageFraction, sitesMapping1Times, sitesMapping2Times, ...
getCoverageDataFrame <- function(path, name) {
    df <- read.table(path, sep=",", header=T)
    numRefSites <- max(df$sitesMapping1Times)
    df$coverageFraction <- df$sitesMapping1Times / numRefSites
    df$multiCoverageFraction <- df$sitesMapping2Times / numRefSites
    df$name <- rep(name, length(df$Genome))
    return(df)
}

# get a base ggplot2 plot comparing multiple coverage results (in a
# single data frame).
getCoveragePlot <- function(coverageResults) {
    return(ggplot(coverageResults, aes(x=name, y=coverageFraction, fill=name)) + geom_bar(stat="identity") + facet_grid( ~ Genome)  + theme_classic())
}

