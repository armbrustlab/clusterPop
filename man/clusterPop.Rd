\name{clusterPop}
\alias{clusterPop}
\title{Cluster Phytoplankton populations from the SeaFlow repository}
\description{
The package analyzes a continuous series of OPP files collected by SeaFlow. The package performs classification of plankton populations based on k-means clustering algorithm and customizable parameters used to define the pre-gating regions. The centroides of each population obtained from the last file is used to clustering the populations for the current file.
}

\usage{
clusterPop(filename,  save.path=getCruisePath(filename), pop.baseline = c("noise","beads", "cocco", "crypto","synecho"), numc = c(2,2),  pop.kmean = c("ultra","nano", "pico","pico"), pop.def.path, noise = NULL, error = 2, diff = 3, min.opp = 800, min.beads = 150, min.synecho=200, lim.largeP = NULL, lim.cocco = 15000, lim.elong = 20000, diff.lim.beads = 0, diff.lim.synecho =0, diff.lim.noise =0, do.plot=TRUE)

}
\arguments{
\item{filename}{list of OPP files}
\item{save.path}{Path to the directory where you wish to output data. Cruise folder by default}
\item{numc}{Number of Phytoplankton smaller (first) and larger (second) than Synecho (excluded Beads, Cryptophytes, Coccolithophores and Synecho). Default is 2 and 2, respectively}
\item{pop.kmean}{Name of the populations defined in 'numc'. Default is c("ultra","nano", "pico","pico"). }
\item{pop.def.tab}{Path to the file that defines how to gate & cluster the populations defined in pop.baseline.}
\item{noise}{Set up the noise threshold for phytoplankton cells. Only cells with chlorophyll value higher than the noise will be clustered}
\item{error}{Fold-difference tolerance between the centroides of two consecutive files. Default is 2}
\item{diff}{Fold-difference tolerance between the size of populations of two consecutive files. Default is 3}
\item{min.opp}{Minimum number of cells for clustering. Default is 800}
\item{min.beads}{Minimum number of beads for reference. Default is 150}
\item{min.synecho}{Minimum number of Synechococcus for reference, if beads number not sufficient. Default is 200}
\item{lim.largeP}{Minimal size for large particles with low red fluorescence}
\item{lim.cocco}{Minimal red fluorescence for Coccolithophores}
\item{lim.elong}{Minimal size for elongated particles (small size particles with high fluorescence). }
\item{diff.lim.beads}{correction for size reference from beads. Default is 0}
\item{diff.lim.synecho}{correction for size reference from Synechococcus. Default is 0}
\item{diff.lim.noise}{correction for noise reference from synecho. Default is 0}
\item{do.plot}{whether saving the cytograms showing the output of the classification. Default is TRUE}

}
\value{
1) A consensus vector file composed of one single column indicating population assignment of each particle from theFCM file.
2) Cytograms showing the phytoplankton populations
3) A list a flagged files where significant differences where observed between two consecutive files 
}

\examples{

example.cruise.name <- 'seaflow_cruise'
seaflow.path <- system.file("extdata", example.cruise.name, package="flowPhyto")
pop.path <- system.file("extdata", example.cruise.name, 'pop.def.tab',package="flowPhyto")
save.path <- '.'

clusterPop(path.cruise=seaflow.path, save.path=save.path, pop.def.path=pop.path, min.beads=0)

}