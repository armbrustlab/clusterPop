\name{clusterPop}
\alias{clusterPop}
\title{Rea-time clustering of phytoplankton populations from the SeaFlow instrument}
\description{
The package analyzes a continuous series of OPP files collected by SeaFlow. The package performs classification of plankton populations based on k-means clustering algorithm and customizable parameters used to define the pre-gating regions. 
}

\usage{
clusterPop(opp.filelist, save.path = getCruisePath(opp.filelist), pop.def.path, initialize=F,  noise = c(0,0), n.pop = c(3,3), tol=0.2,h0=0.5, h=1.0, min.opp = 800, min.beads = 200, min.synecho=200, lim.largeP = NULL,  lim.beads = 0, lim.synecho =0, lim.noise=c(10000,4000), save.plot=F)

}
\arguments{
\item{opp.filelist}{list of OPP files}
\item{save.path}{Path to the directory where you wish to output data. Cruise folder by default}
\item{pop.def.tab}{Path to the file that defines how to gate & cluster the populations defined in pop.baseline.}
\item{initialize}{Cluster the first of the OPP files to help set up each parameters.}
\item{noise}{Noise value for fsc_small and chl_small, respectively. Only particles above noise will be clustered into phytoplankton populations. Default is c(0,0)}
\item{n.pop}{Number of  Smaller and Larger phytoplankton than 1 micron Beads or Synecho if no beads. Default is c(3, 3)}
\item{tol}{Peaks merging tolerance level makes the density function of all modes similar to each other. See flowPeaks package for details. Default is 0.2}
\item{h0}{Scaling factor of the common variance matrix makes the density function smoother. See flowPeaks package for details. Default is 0.5}
\item{h}{Scaling factor of the individual variance matrix makes neighboring peaks merging easier. See flowPeaks package for details. Default is 1.0}
\item{min.opp}{Minimum number of cells for clustering. Default is 800}
\item{min.beads}{Minimum number of beads for reference. Default is 150}
\item{min.synecho}{Minimum number of Synechococcus for reference, if beads number not sufficient. Default is 200}
\item{lim.largeP}{Minimal size for large particles with low red fluorescence}
\item{diff.lim.beads}{correction for size reference from beads. Default is 0}
\item{lim.synecho}{correction for size reference from Synechococcus. Default is 0}
\item{lim.noise}{Chl_small median and fsc_small width minima of Small phytoplankton. If below, population converted to 'noise'. Default is c(10000,4000)}
\item{save.plot}{whether saving the cytograms showing the output of the classification. Default is FALSE}

}
\value{
1) A vector file composed of one single column indicating population assignment of each particle.
2) Cytograms showing the phytoplankton populations
}

\examples{

example.cruise.name <- 'seaflow_cruise'
seaflow.path <- system.file("extdata", example.cruise.name, package="flowPhyto")
pop.path <- system.file("extdata", example.cruise.name, 'pop.def.tab',package="flowPhyto")
save.path <- '.'

#clusterPop(path.cruise=seaflow.path, save.path=save.path, pop.def.path=pop.path, min.beads=0)

}