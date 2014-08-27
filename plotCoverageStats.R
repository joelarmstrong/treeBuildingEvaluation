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

df <- rbind.fill(getCoverageDataFrame("development_humanCoverage", "development"), getCoverageDataFrame("bestRecon_10trees_3out_nucLikelihood_1.0breakpoint_humanCoverage", "Breakpoint1.0"), getCoverageDataFrame("bestRecon_10trees_3out_nucLikelihood_5.0breakpoint_humanCoverage", "Breakpoint5.0"), getCoverageDataFrame("bestRecon_10trees_3out_nucLikelihood_noBreakpoint_humanCoverage", "Breakpoint0.0"))
pdf("testCoverage.pdf")
coveragePlot <- ggplot(df, aes(x=name, y=multiCoverageFraction, fill=name)) + geom_bar(stat="identity") + facet_grid(~ Genome) + theme_classic()
print(coveragePlot)
dev.off()
