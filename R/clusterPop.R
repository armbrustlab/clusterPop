.Prompt <- function(){
			answer <- readline("What do you want to do ? \n	1: Delete these settings\n	2: Keep these settings and stop \n	3: Apply these settings to the rest of the files\n Answer =")
			if(answer == 1) x <- 0
			if(answer == 2) x <- 1
			if(answer == 3) x <- 2
			return(x)
			}

.nextplot <- function(){
			answer <- readline("Next Plot (y / n)?")
			if(answer == 'y') x <- 0
			if(answer == 'n') x <- 1
			return(x)
			}



clusterPop <- function(opp.filelist, save.path = getCruisePath(opp.filelist), numc = c(2,1),  pop.kmean = c("ultra","nano", "pico","pico"), pop.def.path, initialize=F, noise = c(0,0), error = 10000, min.opp = 800, min.beads = 200, min.synecho=200, lim.largeP = NULL, lim.cocco = 15000, lim.pe = 20000, diff.lim.beads = 0, diff.lim.synecho =0, lim.noise=100000, do.plot=TRUE){
	

#library(cluster)	
require(flowPhyto)
require(flowMeans)
require(flowPeaks)
require(caroline)
	
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))
	
# read pop.def.tab
pop.def <- readPopDef(pop.def.path) 
	
# list of phytoplankton populations
phyto <- c("noise","beads", "cocco", "crypto","synecho", pop.kmean) 
	
	numc1 <- numc[1]
	numc2 <- numc[2]
	
	iter.max <- 300 # kmeans calculation
	nstart <- 50 # random points to start flowMeans
	breaks <- 25 # histograms
	do <- 2

	
prev.file <- opp.filelist[1]

#################################################
### check that there is no cluster settings already ###
#################################################

if(initialize){

###################
###################
#### prev FILE ###
###################	
###################
	
	print(paste("Clustering the first file in the list:", prev.file))
	
	opp <- readSeaflow(prev.file) 
	opp$pop <- 0

	lim <- 0.5*2^16

	###############################
	### clustering pop.baseline ###
	###############################
	
	#Cluster Beads
	x <- subset(opp, pop==0)
	beads <- subset(x, pe > 65500) ## PE-saturated particles assigned as 'beads
	opp[row.names(beads),'pop'] <- 2
	x <- subset(opp, pop==0)
	xvar <- pop.def["beads", "xvar"]
	yvar <- pop.def["beads", "yvar"]
	beads <- subset(x, x[,yvar] > 0.5*x[,xvar] + pop.def["beads", "lim"] & x[,xvar] > pop.def["beads", "xmin"] & x[,yvar] > pop.def["beads", "ymin"] & x[,xvar] < pop.def["beads", "xmax"] & x[,yvar] < pop.def["beads", "ymax"])
	opp[row.names(beads),'pop'] <- 2

	if(nrow(beads) > min.beads){
		lim <- quantile(beads$fsc_small,0.1) -diff.lim.beads
		print(paste("Light Scatter cut-off between particles 'smaller' and 'larger' than Beads:", lim))
		}
	
	#Cluster Synecho
	x <- subset(opp, pop==0)
	yvar <- pop.def["synecho", "yvar"]
	xvar <- pop.def["synecho", "xvar"]
	synecho <- subset(x, x[,yvar] > x[,xvar] - pop.def["synecho", "lim"] & x[,xvar] > pop.def["synecho", "xmin"] & x[,yvar] > pop.def["synecho", "ymin"] & x[,xvar] < pop.def["synecho", "xmax"] & x[,yvar] < pop.def["synecho", "ymax"])
	opp[row.names(synecho), 'pop'] <- 5

	if(nrow(beads) < min.beads & nrow(synecho) > min.synecho){
		lim <- median(synecho$chl_small) -diff.lim.synecho
		print(paste("CHL cut-off between particles 'smaller' and 'larger' than Synecho:", lim))

			}
	
	#Cluster Cryptophytes
	x <- subset(opp, pop==0)
	yvar <- pop.def["crypto", "yvar"]
	xvar <- pop.def["crypto", "xvar"]
	crypto <- subset(x, x[,yvar] > x[,xvar] + pop.def["crypto", "lim"] & x[,xvar] > pop.def["crypto", "xmin"] & x[,yvar] > pop.def["crypto", "ymin"])
	opp[row.names(crypto), 'pop'] <- 4


	#Cluster Large Particles
	if(!is.null(lim.largeP)){
		x <- subset(opp, pop==0)
		hetero <- subset(x, x[,"fsc_small"] > 1.5*x[,"chl_small"] + lim.largeP)
		opp[row.names(hetero), 'pop'] <- 1
		}


	#Cluster Coccolithophores
	
	x <- subset(opp, pop==0)
	yvar <- pop.def["cocco", "yvar"] 
	xvar <- pop.def["cocco", "xvar"]
	cocco <- subset(x, x[,yvar] > x[,xvar] + pop.def["cocco", "lim"] & x[,xvar] > pop.def["cocco", "xmin"] & x[,yvar] > pop.def["cocco", "ymin"]& 	x[,"chl_small"] > lim.cocco)
	opp[row.names(cocco), 'pop'] <- 3

	
	#Cluster Noise
	if(!is.null(noise)){
		x <- subset(opp, pop==0)
		bg <- subset(x, chl_small < noise[1] | fsc_small < noise[2])
		opp[row.names(bg),'pop'] <- 1
		}
	

	###############################################
	####   K-means pe-negative phytoplankton   ####
	###############################################

	### cluster phytoplankton into 'large' and 'small' cell size populations
	x <- subset(opp, pop == 0 & chl_small > 50000 & fsc_small > 55000) # Chl-saturating particles and largest cells assigned as 'nano'
	opp[row.names(x),'pop'] <- 5 + numc1

	x <- subset(opp, pop == 0)
	x <- x[,c(5,9)]

	km <- flowPeaks(x, tol=0.1, h0=0.5, h=1.5)

	opp[row.names(x),'pop'] <- km$peaks.cluster + 100
	
		for(i in 1:((max(km$peaks.cluster)))){
			df <- subset(opp, pop == i+100)
			print(paste("Light Scatter mediam for cluster",i,": ", median(df$chl_small)))
				if(median(df$chl_small) > lim) opp[row.names(df),'pop'] <- "y"
				else opp[row.names(df),'pop'] <- "z"	
			}


	### cluster cells larger than 1 um Beads or Synecho (if no beads)
	y <- subset(opp, pop == "y")
	y <- y[,c(5,5,5,9)]

	prev.km.big <- try(flowMeans(y, NumC=numc1,MaxN=numc1+numc2+2, nstart=nstart, Standardize=F, Update='None'))
	if(class(prev.km.big) == "try-error") prev.km.big <- try(flowMeans(y, NumC=numc1,MaxN=numc1+numc2+1, nstart=nstart, Standardize=F, Update='None'))
	if(class(prev.km.big) == "try-error") prev.km.big <- try(flowMeans(y, NumC=numc1,MaxN=numc1+numc2, nstart=nstart, Standardize=F, Update='None'))
	
	opp[row.names(y),'pop'] <- prev.km.big@Label + 100
	
	fsc.sort <- names(sort(by(opp[as.numeric(opp$pop) > 100,"fsc_small"], opp[as.numeric(opp$pop) > 100,"pop"], median)))

		for(i in 1:numc1){
			df <- subset(opp, pop == as.numeric(fsc.sort[i]))
			opp[row.names(df),'pop'] <- i+5
			}
	


	### cluster cells smaller than 1 um Beads or Synecho (if no beads)
	z <- subset(opp, pop == "z")
	z <- z[,c(5,9)]

	if(median(z$fsc_small) < 5000 & median(z$chl_small) < 10000){
		print("Electrical noise detected!")
		prev.km.small <- flowPeaks(z, tol=0.3, h0=0.1)
		}else prev.km.small <- flowPeaks(z, tol=0.1, h0=0.5, h=1.5)

	opp[row.names(z),'pop'] <- prev.km.small$peaks.cluster + 100
	
		for(i in 1:(max(prev.km.small$peaks.cluster))){
			df <- subset(opp, pop == (i+100))
			quant.chl <- median(df$chl_small)
			print(paste("25% Red Fluo=", quant.chl, ", OPP number=", nrow(df)))
			if(quant.chl < lim.noise){
							opp[row.names(df),'pop'] <- 1
							print(paste("pop", i+numc1+5, "converts to", 1))
							# opp[row.names(subset(df, chl_small > lim2)), 'pop'] <- numc1 + i + 5

					}else{
							opp[row.names(df),'pop'] <- numc1 + 1 + 5
							}
		}

	####################
	### PLOT CLUSTER ###
	####################
		
	hist1 <- hist(x$fsc_small, breaks=breaks, plot=FALSE)
	hist2 <- hist(x$chl_small, breaks=breaks, plot=FALSE)
	hist3 <- hist(y$fsc_small, breaks=breaks, plot=FALSE)
	hist4 <- hist(y$chl_small, breaks=breaks, plot=FALSE)
	hist5 <- hist(z$fsc_small, breaks=breaks, plot=FALSE)
	hist6 <- hist(z$chl_small, breaks=breaks, plot=FALSE)
		
	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),4,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)	

	par(mar=c(6,6,1,1))
	plot(km,drawboundary=T,drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)
	mtext("Cut-off Large and Small Phytoplankton", 3,cex=0.7)
	abline(h=lim,col='red',lwd=4)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
	
	mtext(paste("file: ",flowPhyto:::.getYearDay(prev.file),'/',basename(prev.file), sep=""), side=3, line=-4, outer=T,cex=1.2)

	par(mar=c(6,6,1,1))
	plot(y[,"fsc_small"],y[,"chl_small"], col=prev.km.big@Label,cex=0.5,xlab='fsc_small', ylab='chl_small')
	mtext("Large Phytoplankton", 3, cex=0.7)
	par(mar=c(0,6,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
		
	par(mar=c(6,6,1,1))
	plot(prev.km.small,drawboundary=T,drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)
	mtext("Small Phytoplankton", 3, cex=0.7)
	par(mar=c(0,6,1,1))
	barplot(hist5$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist6$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(def.par)


	go <- .nextplot()	
	
	if(go == 1) do <- 0
		

	if(go == 0){

	######################
	### SAVE SETTINGS  ###
	######################

	for(i in 1:(length(phyto))){
		p <- subset(opp, pop == i)
		opp[row.names(p),'pop'] <- phyto[i]
		}
	
	all.settings <- aggregate(opp[,c("fsc_small", "chl_small")], FUN=median, by=list(opp$pop))		
	settings <- all.settings[all.settings$Group.1 == "nano" | all.settings$Group.1 == "ultra" | all.settings$Group.1 == "synecho",]
	write.delim(settings, paste(save.path,"prev.settings",sep=""))



	##################
	### PLOT PHYTO ###
	##################

	
	
		
	hist1 <- hist(opp$fsc_small, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist2 <- hist(opp$chl_small, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist3 <- hist(opp$pe, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist4 <- hist(opp$fsc_perp, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)

	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),4,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'chl_small', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
	
	mtext(paste("file: ",flowPhyto:::.getYearDay(prev.file),'/',basename(prev.file), sep=""), side=3, line=-4, outer=T,cex=1.2)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'chl_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'fsc_perp', pop.def=pop.def, add.legend=TRUE)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(def.par)





	#######################
	### Keep settings ? ###
	#######################
	do <- .Prompt()	
	
	if(do == 1) stop

	if(do == 0){
		system(paste("rm ",save.path,"prev.settings",sep=""))
		stop
		} 
	}
	
}

################################################
################################################
##### ALL CRUISES FILES based on prev FILE ####
################################################
################################################
if(initialize == F){ 
	lim <- 0.5*2^16
	}
	
if(do == 2){ 


outlier.table <- try(read.csv(paste(save.path,"c.outliers",sep="")),silent=T)
	
	if(class(outlier.table) == "try-error"){
		outlier.table <- data.frame(day=NA,file=NA)
		write.csv(outlier.table, file=paste(save.path,"c.outliers", sep=""), row.names=FALSE, quote=FALSE)
		}


for (file in opp.filelist){
	opp <- readSeaflow(file)
	opp$pop <- 0

	day <- flowPhyto:::.getYearDay(file)
	
	## Check number OPP before Clustering 
	
	if(nrow(opp[opp[,"chl_small"]>20000,]) < min.opp){
		print(paste("not enough OPP for clustering in:",file))
		outlier <- data.frame(day=day, file=getFileNumber(file))
		outlier.table <- rbind(outlier.table, outlier)
		write.csv(outlier.table, file=paste(save.path,"c.outliers", sep=""), row.names=FALSE, quote=FALSE)
	
		png(paste(save.path, day,"/",basename(file),".",getFileNumber(file),".outlier.gif", sep=""),width=9, height=12, unit='in', res=100)
	
			par(mfrow=c(1,2), pty='s')
				plot(opp[,"fsc_small"], opp[,"chl_small"], pch=1, cex=0.7, xlab="fsc_small", ylab="chl_small",ylim=c(0,2^16),xlim=c(0,2^16),col= densCols(opp[,"fsc_small"], opp[,"chl_small"], colramp=.rainbow.cols))
				plot(opp[,"fsc_small"], opp[,"pe"], pch=1, cex=0.7, xlab="fsc_small", ylab="pe",ylim=c(0,2^16),xlim=c(0,2^16),col= densCols(opp[,"fsc_small"], opp[,"pe"], colramp=.rainbow.cols))
				mtext(paste("Nb cells <", min.opp,": FLAGGED !"), side=3, outer=T, line=-3,cex=3,col='red')		
				mtext(paste("file: ",flowPhyto:::.getYearDay(file),'/',basename(file), sep=""), side=3, line=-4, outer=T,cex=1.2)

		dev.off()
		
		
		next
		}
					
		

	###############################
	### clustering pop.baseline ###
	###############################
	print(paste("clustering", file))

	#Cluster Beads
	x <- subset(opp, pop==0)
	beads <- subset(x, pe > 65500) ## PE-saturated particles assigned as 'beads
	opp[row.names(beads),'pop'] <- 2
	xvar <- pop.def["beads", "xvar"]
	yvar <- pop.def["beads", "yvar"]
	beads <- subset(x, x[,yvar] > 0.5*x[,xvar] + pop.def["beads", "lim"] & x[,xvar] > pop.def["beads", "xmin"] & x[,yvar] > pop.def["beads", "ymin"] & x[,xvar] < pop.def["beads", "xmax"] & x[,yvar] < pop.def["beads", "ymax"])
	opp[row.names(beads),'pop'] <- 2

		if(nrow(beads) > min.beads){
			prev.beads <- beads
			lim <- quantile(beads$fsc_small,0.1) -diff.lim.beads
			}
			else print("not enough beads")
	

	#Cluster Synecho
	x <- subset(opp, pop==0)
	yvar <- pop.def["synecho", "yvar"]
	xvar <- pop.def["synecho", "xvar"]
	synecho <- subset(x, x[,yvar] > x[,xvar] - pop.def["synecho", "lim"] & x[,xvar] > pop.def["synecho", "xmin"] & x[,yvar] > pop.def["synecho", "ymin"] & x[,xvar] < pop.def["synecho", "xmax"] & x[,yvar] < pop.def["synecho", "ymax"])
	opp[row.names(synecho), 'pop'] <- 5

		if(nrow(synecho) > min.synecho){
			lim.synecho <- median(synecho$fsc_small)-diff.lim.synecho
				if(nrow(beads) < min.beads) lim <- lim.synecho
			}
		if(nrow(beads) < min.beads & nrow(synecho) < min.synecho){
			print("too few Synecho")
			}

	#Cluster Cryptophytes
	x <- subset(opp, pop==0)
	yvar <- pop.def["crypto", "yvar"]
	xvar <- pop.def["crypto", "xvar"]
	crypto <- subset(x, x[,yvar] > x[,xvar] + pop.def["crypto", "lim"] & x[,xvar] > pop.def["crypto", "xmin"] & x[,yvar] > pop.def["crypto", "ymin"])
	opp[row.names(crypto), 'pop'] <- 4
	
	
	#Cluster "heterotrophs"
	if(!is.null(lim.largeP)){
		x <- subset(opp, pop==0)
		hetero <- subset(x, x[,"fsc_small"] > 1.5*x[,"chl_small"] + lim.largeP)
		opp[row.names(hetero), 'pop'] <- 1
		}

	#Cluster Coccolithophores
	x <- subset(opp, pop==0)
	yvar <- pop.def["cocco", "yvar"] 
	xvar <- pop.def["cocco", "xvar"]
	cocco <- subset(x, x[,yvar] > x[,xvar] + pop.def["cocco", "lim"] & x[,xvar] > pop.def["cocco", "xmin"] & x[,yvar] > pop.def["cocco", "ymin"]& 	x[,"chl_small"] > lim.cocco)
	opp[row.names(cocco), 'pop'] <- 3


	#Cluster Noise
	if(!is.null(noise)){
		x <- subset(opp, pop==0)
		bg <- subset(x, chl_small < noise[1] | fsc_small < noise[2])
		opp[row.names(bg),'pop'] <- 1
		}

			
	########################
	#####   K-means    #####
	########################
	
	#### cluster phytoplankton into a number of populations specified in NUMC
	x <- subset(opp, pop == 0 & chl_small > 50000 & fsc_small > 55000) # Chl-saturating particles and largest cells assigned as 'nano'
	opp[row.names(x),'pop'] <- 5 + numc1

	x <- subset(opp, pop == 0)
	x <- x[,c(5,9)]

	km <- flowPeaks(x, tol=0.1, h0=0.5, h=1.5)

	opp[row.names(x),'pop'] <- km$peaks.cluster + 100
	
		for(i in 1:((max(km$peaks.cluster)))){
			df <- subset(opp, pop == i+100)
			#print(paste("Light Scatter mediam for cluster",i,": ", median(df$chl_small)))
				if(median(df$chl_small) > lim) opp[row.names(df),'pop'] <- "y"
				else opp[row.names(df),'pop'] <- "z"	
			}


	#### cluster cells larger than 1 um Beads or Synecho (if no beads)
	y <- subset(opp, pop == "y")
	y <- y[,c(5,5,5,9)]
	
	prev.km.big <- try(flowMeans(y, NumC=numc1,MaxN=numc1+numc2+2, nstart=nstart, Standardize=F, Update='None'),silent=F)
	if(class(prev.km.big) == "try-error") prev.km.big <- try(flowMeans(y, NumC=numc1,MaxN=numc1+numc2+1, nstart=nstart, Standardize=F, Update='None'))
	if(class(prev.km.big) == "try-error") prev.km.big <- try(flowMeans(y, NumC=numc1,MaxN=numc1+numc2, nstart=nstart, Standardize=F, Update='None'))
	
			if(class(prev.km.big) == "try-error"){
				print("flowMeans error, FLAGGED FILE")
				outlier <- data.frame(day=day, file=getFileNumber(file))
				outlier.table <- rbind(outlier.table, outlier)
				write.csv(outlier.table, file=paste(save.path,"c.outliers", sep=""), row.names=FALSE, quote=FALSE)
				
			png(paste(save.path, day,"/",basename(file),".",getFileNumber(file),".outlier.gif", sep=""),width=9, height=12, unit='in', res=100)

				hist1 <- hist(x$fsc_small, breaks=breaks, plot=FALSE)
				hist2 <- hist(x$chl_small, breaks=breaks, plot=FALSE)
				hist3 <- hist(y$fsc_small, breaks=breaks, plot=FALSE)
				hist4 <- hist(y$chl_small, breaks=breaks, plot=FALSE)
		
				def.par <- par(no.readonly = TRUE) # save default, for resetting...
				nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),4,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)	

				par(mar=c(6,6,1,1))
				plot(km,drawboundary=T,drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)
				mtext("Cut-off Large and Small Phytoplankton", 3,cex=0.7)
				mtext(paste(file), side=3, line=-2, outer=T)
				abline(h=lim,col='red',lwd=4)
				par(mar=c(0,6,1,1))
				barplot(hist1$counts, axes=FALSE, space=0, col=NA)
				par(mar=c(6,0,1,1))
				barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
	
				mtext(paste("file: ",flowPhyto:::.getYearDay(file),'/',basename(file), sep=""), side=3, line=-4, outer=T,cex=1.2)

				par(mar=c(6,6,1,1))
				plot(y[,"fsc_small"],y[,"chl_small"], col=prev.km.big@Label,cex=0.5,xlab='fsc_small', ylab='chl_small')
				mtext("Large Phytoplankton", 3, cex=0.7)
				par(mar=c(0,6,1,1))
				barplot(hist3$counts, axes=FALSE, space=0, col=NA)
				par(mar=c(6,0,1,1))
				barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
		
				par(def.par)
				
			dev.off()
				
				next
				}
	opp[row.names(y),'pop'] <- prev.km.big@Label + 100
							
	fsc.sort <- names(sort(by(opp[as.numeric(opp$pop) > 100,"fsc_small"], opp[as.numeric(opp$pop) > 100,"pop"], median)))

		for(i in 1:numc1){
			df <- subset(opp, pop == as.numeric(fsc.sort[i]))
			opp[row.names(df),'pop'] <- i+5
			}
	
	
	#### cluster cells smaller than 1 um Beads or Synecho (if no beads)
	z <- subset(opp, pop == "z")
	z <- z[,c(5,9)]

	if(nrow(z) < 3){
		print(paste("found only", nrow(z), "small cells, FLAGGED FILE"))
			outlier <- data.frame(day=day, file=getFileNumber(file))
			outlier.table <- rbind(outlier.table, outlier)
			write.csv(outlier.table, file=paste(save.path,"c.outliers", sep=""), row.names=FALSE, quote=FALSE)
			
			png(paste(save.path, day,"/",basename(file),".",getFileNumber(file),".cluster.gif", sep=""),width=9, height=12, unit='in', res=100)
		
				hist1 <- hist(x$fsc_small, breaks=breaks, plot=FALSE)
				hist2 <- hist(x$chl_small, breaks=breaks, plot=FALSE)
				hist3 <- hist(y$fsc_small, breaks=breaks, plot=FALSE)
				hist4 <- hist(y$chl_small, breaks=breaks, plot=FALSE)
			
				def.par <- par(no.readonly = TRUE) # save default, for resetting...
				nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),4,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)	

				par(mar=c(6,6,1,1))
				plot(km,drawboundary=T,drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)
				mtext("Cut-off Large and Small Phytoplankton", 3,cex=0.7)
				mtext(paste(file), side=3, line=-2, outer=T)
				abline(h=lim,col='red',lwd=4)
				par(mar=c(0,6,1,1))
				barplot(hist1$counts, axes=FALSE, space=0, col=NA)
				par(mar=c(6,0,1,1))
				barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
	
				par(mar=c(6,6,1,1))
				plot(y[,"fsc_small"],y[,"chl_small"], col=prev.km.big@Label,cex=0.5,xlab='fsc_small', ylab='chl_small')
				mtext("Large Phytoplankton", 3, cex=0.7)
				par(mar=c(0,6,1,1))
				barplot(hist3$counts, axes=FALSE, space=0, col=NA)
				par(mar=c(6,0,1,1))
				barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

				mtext(paste("file: ",flowPhyto:::.getYearDay(file),'/',basename(file), sep=""), side=3, line=-4, outer=T,cex=1.2)
			
			dev.off()
			
				}
	
	if(nrow(z) > 3) {
		
		if(median(z$fsc_small) < 5000 & median(z$chl_small) < 10000){
		
		print("Electrical noise detected!")
		prev.km.small <- try(flowPeaks(z, tol=0.3, h0=0.1))
		
		}else{ 
		
		prev.km.small <- try(flowPeaks(z))
		
		}

	    if(class(prev.km.small) == "try-error"){
				print("flowPeaks error, FLAGGED FILE")
				outlier <- data.frame(day=day, file=getFileNumber(file))
				outlier.table <- rbind(outlier.table, outlier)
				write.csv(outlier.table, file=paste(save.path,"c.outliers", sep=""), row.names=FALSE, quote=FALSE)
				
			png(paste(save.path, day,"/",basename(file),".",getFileNumber(file),".cluster.gif", sep=""),width=9, height=12, unit='in', res=100)
		
				hist1 <- hist(x$fsc_small, breaks=breaks, plot=FALSE)
				hist2 <- hist(x$chl_small, breaks=breaks, plot=FALSE)
				hist5 <- hist(z$fsc_small, breaks=breaks, plot=FALSE)
				hist6 <- hist(z$chl_small, breaks=breaks, plot=FALSE)
		
				def.par <- par(no.readonly = TRUE) # save default, for resetting...
				nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),4,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)	

				par(mar=c(6,6,1,1))
				plot(km,drawboundary=T,drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)
				mtext("Cut-off Large and Small Phytoplankton", 3,cex=0.7)
				mtext(paste(file), side=3, line=-2, outer=T)
				abline(v=lim,col='red',lwd=4)
				par(mar=c(0,6,1,1))
				barplot(hist1$counts, axes=FALSE, space=0, col=NA)
				par(mar=c(6,0,1,1))
				barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
	
				mtext(paste("file: ",flowPhyto:::.getYearDay(file),'/',basename(file), sep=""), side=3, line=-4, outer=T,cex=1.2)

				par(mar=c(6,6,1,1))
				plot(prev.km.small,drawboundary=T,drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)
				mtext("Small Phytoplankton", 3, cex=0.7)
				par(mar=c(0,6,1,1))
				barplot(hist5$counts, axes=FALSE, space=0, col=NA)
				par(mar=c(6,0,1,1))
				barplot(hist6$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

			dev.off()
			
		par(def.par)

				}

		if(class(prev.km.small) != "try-error"){

			opp[row.names(z),'pop'] <- prev.km.small$peaks.cluster + 100
	
			for(i in 1:max((prev.km.small$peaks.cluster))){
				df <- subset(opp, pop == (i+100))
				quant.chl <- median(df$chl_small)
				if(quant.chl < lim.noise){
							opp[row.names(df),'pop'] <- 1

					}else{
							opp[row.names(df),'pop'] <- numc1 + 1 + 5
							}
				}	
			}
		}
					
	#################
	### class vct ###
	#################

	for(i in 1:(length(phyto))){
		p <- subset(opp, pop == i)
		opp[row.names(p),'pop'] <- phyto[i]
			}
		
	write.table(opp$pop, paste(save.path, day,"/",basename(file),".",getFileNumber(file),'-class.vct',sep=""), row.names=FALSE, col.names='pop', quote=FALSE)
	
	#################
	### Consensus ###
	#################
	
	opp$support <- 3
	caroline:::write.delim(opp[,c("pop","support")],paste(save.path, day,"/",basename(file),'.consensus.vct',sep=""))
	
	##############
	### Census ###
	##############
	census <- census(opp$pop, pop.def=pop.def)
	write.delim(census,paste(save.path, day,"/",basename(file),'.census.tab',sep=""))
	
	
	###############				
	###	flagged ###
	###############
	prev.settings <- read.delim(paste(save.path,"prev.settings",sep=""))

	all.settings <- aggregate(opp[,c("fsc_small", "chl_small")], FUN=median, by=list(opp$pop))		
	settings <- all.settings[all.settings$Group.1 == "nano" | all.settings$Group.1 == "ultra" | all.settings$Group.1 == "synecho",]
	
	if(any(abs(settings[,"fsc_small"]-prev.settings[,"fsc_small"])) > error | 
		any(abs(settings[,"chl_small"]-prev.settings[,"chl_small"])) > error) {
			flag <- 1
			print("FLAGGED FILE")
			outlier <- data.frame(day=day, file=getFileNumber(file))
			outlier.table <- rbind(outlier.table, outlier)
			write.csv(outlier.table, file=paste(save.path,"c.outliers", sep=""), row.names=FALSE, quote=FALSE)

		}else{
			flag <- 0
			write.delim(settings, paste(save.path,"prev.settings",sep=""))
		}


	############
	### PLOT ###
	############
	if(do.plot==TRUE){

	png(paste(save.path, day,"/",basename(file),".",getFileNumber(file),".class.gif", sep=""),width=9, height=12, unit='in', res=100)
	if(class(opp) != "try-error"){

	hist1 <- hist(opp$fsc_small, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist2 <- hist(opp$chl_small, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist3 <- hist(opp$pe, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist4 <- hist(opp$fsc_perp, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)


	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),4,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'chl_small', pop.def=pop.def)
			if(flag == 1) mtext("flagged !", side=3, outer=T, line=-3,cex=3,col='red')
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	mtext(paste("file: ",flowPhyto:::.getYearDay(file),'/',basename(file), sep=""), side=3, line=-4, outer=T,cex=1.2)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'chl_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'fsc_perp', pop.def=pop.def, add.legend=TRUE)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(def.par)

	dev.off()
		}
	
	}
	
	
	####################
	### EMPTY MEMORy ###	
	####################	
	rm(opp, noise, synecho, crypto, beads, cocco, hetero, prev.km, prev.km.big, prev.km.small, p, df, x, y, z, hist1, hist2, hist3, hist4)	
	}

	print("Combining census.tab ...")
	caroline:::write.delim(combineCensusFiles(save.path),paste(save.path,'census.tab',sep=""))
	
	
	}		

}