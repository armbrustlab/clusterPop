.Prompt <- function(){
			answer <- readline("Do you want to apply these settings (y or n)? ")
			if(substr(answer,1,1) == "n" | substr(answer,1,1) == "N"| substr(answer,1,1) == "NO" | substr(answer,1,1) == "no") x <- 0
			if(substr(answer,1,1) == "y" | substr(answer,1,1) == "Y"| substr(answer,1,1) == "YES" | substr(answer,1,1) == "yes") x <- 1
			return(x)
			}





clusterPop <- function(filename, save.path = getCruisePath(filename), pop.baseline = c("noise","beads", "cocco", "crypto","synecho"), numc = c(2,2),  pop.kmean = c("ultra","nano", "pico","pico"), pop.def.path, noise = NULL, error = 2, diff = 3, min.opp = 800, min.beads = 150, min.synecho=200, lim.largeP = NULL, lim.cocco = 15000, lim.elong = 20000, diff.lim.beads = 0, diff.lim.synecho =0, diff.lim.noise=0, do.plot=TRUE){
	

	
require(flowPhyto, quietly=T)
require(flowPeaks, quietly=T)
require(cluster, quietly=T)
	
# read pop.def.tab
pop.def <- readPopDef(pop.def.path) 
	
# list of phytoplankton populations
phyto <- c(pop.baseline, pop.kmean) 
	

	numc1 <- numc[1]
	numc2 <- numc[2]
	
	iter.max <- 300 # kmeans calculation
	breaks <- 25 # histograms

	
prev.file <- filename[1]
lim <- lim2 <- 0

##############################
### LOAD results last file ###
##############################
prev.km <- try(load(paste(save.path,"prev.km.RData",sep="")), silent=T)
prev.km.big <- try(load(paste(save.path,"prev.km.big.RData",sep="")), silent=T)
prev.km.small <- try(load(paste(save.path,"prev.km.small.RData",sep="")), silent=T)
prev.synecho <- try(load(paste(save.path,"prev.synecho.RData",sep="")), silent=T)
do <- 1	

if(class(prev.km) == "try-error"){

###################
###################
#### prev FILE ###
###################	
###################
	
	print(paste("Clustering the first file in the list:", filename[1]))
	
	opp <- readSeaflow(prev.file) 
	opp$pop <- 0

	
	###############################
	### clustering pop.baseline ###
	###############################
	
	#Cluster Beads
	x <- subset(opp, pop==0)
	xvar <- pop.def["beads", "xvar"]
	yvar <- pop.def["beads", "yvar"]
	beads <- subset(x, x[,yvar] > 0.5*x[,xvar] + pop.def["beads", "lim"] & x[,xvar] > pop.def["beads", "xmin"] & x[,yvar] > pop.def["beads", "ymin"] & x[,xvar] < pop.def["beads", "xmax"] & x[,yvar] < pop.def["beads", "ymax"])
	opp[row.names(beads),'pop'] <- 2

	if(nrow(beads) > min.beads){
		lim <- median(beads$chl_small) -diff.lim.beads
		print(paste("Red fluorescence cut-off between particles 'smaller' and 'larger' than beads:", lim))
		}
	
	#Cluster "heterotrophs"
	if(!is.null(lim.largeP)){
		x <- subset(opp, pop==0)
		hetero <- subset(x, x[,"fsc_small"] > 1.5*x[,"chl_small"] + lim.largeP)
		opp[row.names(hetero), 'pop'] <- 1
		}

	#Cluster Synecho
	x <- subset(opp, pop==0)
	yvar <- pop.def["synecho", "yvar"]
	xvar <- pop.def["synecho", "xvar"]
	synecho <- subset(x, x[,yvar] > x[,xvar] - pop.def["synecho", "lim"] & x[,xvar] > pop.def["synecho", "xmin"] & x[,yvar] > pop.def["synecho", "ymin"] & x[,xvar] < pop.def["synecho", "xmax"] & x[,yvar] < pop.def["synecho", "ymax"])
	prev.synecho <- synecho
	save(prev.synecho, file=paste(save.path,"prev.synecho.RData",sep=""))
	
	opp[row.names(synecho), 'pop'] <- 5

	if(nrow(beads) < min.beads & nrow(synecho) > min.synecho){
		lim <- median(synecho$chl_small) -diff.lim.synecho
		print(paste("Red fluorescence cut-off between particles 'smaller' and 'larger' than synecho:", lim))

			}
	
	lim2 <- median(synecho$chl_small) - diff.lim.noise

	#Cluster Noise
	if(!is.null(noise)){
		x <- subset(opp, pop==0)
		bg <- subset(x, chl_small < noise)
		opp[row.names(bg),'pop'] <- 1
		}
	
	#Cluster Coccolithophores
	
	if(length(which(pop.baseline == "cocco"))>0){
		x <- subset(opp, pop==0)
		yvar <- pop.def["cocco", "yvar"] 
		xvar <- pop.def["cocco", "xvar"]
		cocco <- subset(x, x[,yvar] > x[,xvar] + pop.def["cocco", "lim"] & x[,xvar] > pop.def["cocco", "xmin"] & x[,yvar] > pop.def["cocco", "ymin"] & x[,"chl_small"] > lim.cocco & x[,"pe"] < lim.cocco)
		opp[row.names(cocco), 'pop'] <- 3
		}
	
	#Cluster Cryptophytes
	x <- subset(opp, pop==0)
	yvar <- pop.def["crypto", "yvar"]
	xvar <- pop.def["crypto", "xvar"]
	crypto <- subset(x, x[,yvar] > x[,xvar] + pop.def["crypto", "lim"] & x[,xvar] > pop.def["crypto", "xmin"] & x[,yvar] > pop.def["crypto", "ymin"])
	opp[row.names(crypto), 'pop'] <- 4

	#cluster Elongated cells
	if(length(which(pop.baseline == "elong"))>0){

	x <- subset(opp, pop==0 & chl_small > lim & fsc_small > pop.def["elong", "xmin"] & chl_small > pop.def["elong", "ymin"])
	x1 <- quantile(x[,"fsc_small"], 0.1) ; x2 <- quantile(x[,"fsc_small"], 0.9) 
	y1 <- quantile(x[,"chl_small"], 0.1) ; y2 <- quantile(x[,"chl_small"], 0.9) 
	z <- (y2-y1)/(x2-x1)
	b <- median(x[,"chl_small"]) - z*median(x[,"fsc_small"])
	elongated <- subset(x, x[,"chl_small"] > z*x[,"fsc_small"] + b + lim.elongated)
	
	fp <- flowPeaks(elongated[,c("fsc_small","chl_small")])
	fpc <- assign.flowPeaks(fp, fp$x, tol=0.2)
	elongated$pop <- fpc

	opp[row.names(elongated[which(elongated$pop == 1),]), 'pop'] <- 6
			
		}		


	#################################################
	####   K-means   pe-negative phytoplankton   ####
	#################################################

	
	# cluster phytoplankton into a number of populations specified in NUMC

	x <- subset(opp, pop == 0)
	x <- x[,c(5,8:10)]
	km <- kmeans(x, numc1+numc2,iter.max)#, nstart=numc1+numc2)
	prev.km <- kmeans(x, km$centers,iter.max)
	save(prev.km,file= paste(save.path,"prev.km.RData",sep=""))
	opp[row.names(x),'pop'] <- prev.km$cluster + length(pop.baseline)
	opp.b <- opp
	
		for(i in 1:(numc1+numc2)){
			df <- subset(opp.b, pop == i+length(pop.baseline))
			print(paste("RED fluo mediam for cluster",i+length(pop.baseline),": ", median(df$chl_small)))
				if(median(df$chl_small) > lim) opp[row.names(df),'pop'] <- "y"
				else opp[row.names(df),'pop'] <- "z"	
			}


	# cluster cells larger than 1 um Beads or Synecho (if no beads)
	y <- subset(opp, pop == "y")
	y <- y[,c(5,8:10)]
	if(numc1 > 1 & nrow(y) > 2){
	prev.km.big <- kmeans(y, numc1,iter.max)#, nstart=numc1)
	save(prev.km.big,file= paste(save.path,"prev.km.big.RData",sep=""))
	opp[row.names(y),'pop'] <- prev.km.big$cluster + length(pop.baseline)
	opp.b <- opp
		
		for(i in 1:numc1){
			fsc <- names(sort(prev.km.big$centers[,"fsc_small"]))[i]
			df <- subset(opp.b, pop == (as.numeric(fsc)+length(pop.baseline)))
			opp[row.names(df),'pop'] <- i+ length(pop.baseline)
			}
	}	
	
	# cluster cells smaller than 1 um Beads or Synecho (if no beads)
	z <- subset(opp, pop == "z")
	z <- z[,c(5,8:10)]
	if(numc2 > 1 & nrow(z) > 2){
	prev.km.small <- kmeans(z, numc2,iter.max)#, nstart=numc2)
	save(prev.km.small,file= paste(save.path,"prev.km.small.RData",sep=""))
	opp[row.names(z),'pop'] <- prev.km.small$cluster + numc1 + length(pop.baseline)
	
		for(i in 1:(numc2)){
			df <- subset(opp, pop == (i+numc1+length(pop.baseline)))
			quant.fsc <- quantile(df$fsc_small, probs=0.9) - quantile(df$fsc_small, probs=0.1)
			quant.chl <- quantile(df$chl_small, probs=0.25)
			#print(paste(quant.fsc, "> 20000 and", quant.chl, "<", lim2))
			if(quant.chl < lim2){
							opp[row.names(df),'pop'] <- 1
							print(paste("pop", i+numc1+length(pop.baseline), "converts to", 1))
							# opp[row.names(subset(df, chl_small > lim2)), 'pop'] <- numc1 + i + length(pop.baseline)

					}else{
							opp[row.names(df),'pop'] <- numc1 + i + length(pop.baseline)
							}
			}

	}	
		
	##################
	### PLOT PHYTO ###
	##################
	

	for(i in 1:(length(phyto))){
		p <- subset(opp, pop == i)
		opp[row.names(p),'pop'] <- phyto[i]
			}

	hist1 <- hist(opp$fsc_small, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist2 <- hist(opp$chl_small, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist3 <- hist(opp$pe, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist4 <- hist(opp$chl_big, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist5 <- hist(opp$fsc_perp, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)

	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),6,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)
	
	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'chl_small', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'chl_big', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'chl_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'fsc_perp', pop.def=pop.def, add.legend=TRUE)
	mtext(paste(prev.file), side=1, line=5, at=-30000)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist5$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(def.par)





	#########################
	### keep prev file ? ###
	#########################
	do <- .Prompt()	


	if(do == 0){
		
		system(paste("rm ",save.path,"*.RData",sep=""))

		stop
	} 
	
	
	



}

################################################
################################################
##### ALL CRUISES FILES based on prev FILE ####
################################################
################################################

	




if(do == 1){ 
	km.big <- km.small <- km <- NULL
	
	outlier.table <- data.frame(day=NA,file=NA)
	write.csv(outlier.table, file=paste(save.path,"c.outliers", sep=""), row.names=FALSE, quote=FALSE)

for (file in filename) {
	
	flag <- 0
	prev.km.big$center <- prev.km.small$center <- prev.km.big$size <- prev.km.small$size <- prev.km.big <- prev.km.small <- NULL
	opp <- readSeaflow(file)
	opp$pop <- 0

	day <- flowPhyto:::.getYearDay(file)
	
	
	## Load previous outlier table
	outlier.table <- read.csv(paste(save.path,"c.outliers",sep=""))
	
	## Load previous cluster output
	
	try(load(paste(save.path,"prev.km.RData",sep="")), silent=T)
	try(load(paste(save.path,"prev.km.big.RData",sep="")), silent=T)
	try(load(paste(save.path,"prev.km.small.RData",sep="")), silent=T)
	try(load(paste(save.path,"prev.synecho.RData",sep="")), silent=T)
	lim2 <- median(prev.synecho$chl_small) - diff.lim.noise
	lim <- lim.synecho <- median(prev.synecho$chl_small)-diff.lim.synecho
			
	## Check number OPP before Clustering 
	
	if(nrow(opp) < min.opp){
		print(paste("not enough OPP for clustering in:",file))
		outlier <- paste(day, ",", getFileNumber(file), sep="")
		outlier.table <- rbind(outlier.table, outlier)
		next
		}
					
	else{
			

	###############################
	### clustering pop.baseline ###
	###############################
	print(paste("clustering", file))

	#Cluster Beads
	x <- subset(opp, pop==0)
	xvar <- pop.def["beads", "xvar"]
	yvar <- pop.def["beads", "yvar"]
	beads <- subset(x, x[,yvar] > 0.5*x[,xvar] + pop.def["beads", "lim"] & x[,xvar] > pop.def["beads", "xmin"] & x[,yvar] > pop.def["beads", "ymin"] & x[,xvar] < pop.def["beads", "xmax"] & x[,yvar] < pop.def["beads", "ymax"])
	opp[row.names(beads),'pop'] <- 2

		if(nrow(beads) > min.beads){
			prev.beads <- beads
			lim.beads <- median(prev.beads$chl_small)-diff.lim.beads
			lim <- lim.beads
			}else print("not enough beads")
	
	#Cluster "heterotrophs"
	if(!is.null(lim.largeP)){
		x <- subset(opp, pop==0)
		hetero <- subset(x, x[,"fsc_small"] > 1.5*x[,"chl_small"] + lim.largeP)
		opp[row.names(hetero), 'pop'] <- 1
		}

	#Cluster Synecho
	x <- subset(opp, pop==0)
	yvar <- pop.def["synecho", "yvar"]
	xvar <- pop.def["synecho", "xvar"]
	synecho <- subset(x, x[,yvar] > x[,xvar] - pop.def["synecho", "lim"] & x[,xvar] > pop.def["synecho", "xmin"] & x[,yvar] > pop.def["synecho", "ymin"] & x[,xvar] < pop.def["synecho", "xmax"] & x[,yvar] < pop.def["synecho", "ymax"])
	opp[row.names(synecho), 'pop'] <- 5

		if(nrow(synecho) > min.synecho){
			prev.synecho <- synecho
			save(prev.synecho, file=paste(save.path,"prev.synecho.RData",sep=""))
			lim.synecho <- median(prev.synecho$chl_small)-diff.lim.synecho
			lim2 <- median(prev.synecho$chl_small)-diff.lim.noise
				if(nrow(beads) < min.beads) lim <- lim.synecho
			}
		if(nrow(beads) < min.beads & nrow(synecho) < min.synecho){
			print("too few Synecho")
			}

	#Cluster Noise
	if(!is.null(noise)){
		x <- subset(opp, pop==0)
		bg <- subset(x, chl_small < noise)
		opp[row.names(bg),'pop'] <- 1
		}
	
	#Cluster Coccolithophores
	if(length(which(pop.baseline == "cocco"))>0){
		x <- subset(opp, pop==0)
		yvar <- pop.def["cocco", "yvar"] 
		xvar <- pop.def["cocco", "xvar"]
		cocco <- subset(x, x[,yvar] > x[,xvar] + pop.def["cocco", "lim"] & x[,xvar] > pop.def["cocco", "xmin"] & x[,yvar] > pop.def["cocco", "ymin"] & x[,"chl_small"] > lim.cocco & x[,"pe"] < lim.cocco)
		opp[row.names(cocco), 'pop'] <- 3
		}
	
	#Cluster Cryptophytes
	x <- subset(opp, pop==0)
	yvar <- pop.def["crypto", "yvar"]
	xvar <- pop.def["crypto", "xvar"]
	crypto <- subset(x, x[,yvar] > x[,xvar] + pop.def["crypto", "lim"] & x[,xvar] > pop.def["crypto", "xmin"] & x[,yvar] > pop.def["crypto", "ymin"])
	opp[row.names(crypto), 'pop'] <- 4

	#cluster Elongated cells
	if(length(which(pop.baseline == "elong"))>0){
	yvar <- pop.def["elong", "yvar"]
	xvar <- pop.def["elong", "xvar"]
	x <- subset(opp, pop==0 & chl_small > lim & fsc_small > pop.def["elong", "xmin"] & chl_small > pop.def["elong", "ymin"])
	x1 <- quantile(x[,"fsc_small"], 0.1) ; x2 <- quantile(x[,"fsc_small"], 0.9) 
	y1 <- quantile(x[,"chl_small"], 0.1) ; y2 <- quantile(x[,"chl_small"], 0.9) 
	z <- (y2-y1)/(x2-x1)
	b <- median(x[,"chl_small"]) - z*median(x[,"fsc_small"])
	elongated <- subset(x, x[,"chl_small"] > z*x[,"fsc_small"] + b + lim.elongated)
	
	fp <- flowPeaks(elongated[,c("fsc_small","chl_small")])
	fpc <- assign.flowPeaks(fp, fp$x, tol=0.2)
	elongated$pop <- fpc

	opp[row.names(elongated[which(elongated$pop == 1),]), 'pop'] <- 6
			
		}		
			
	########################
	#####   K-means    #####
	########################
	
	# cluster phytoplankton into a number of  populations specify in NUMC
	x <- subset(opp, pop == 0) 
	x <- x[,c(5,8:10)]
	km <- try(kmeans(x, prev.km$centers,iter.max))
	if(class(km) == "try-error"){
		km <- kmeans(x, numc1+numc2, iter.max); km <- kmeans(x,km$medoids,iter.max)
		print("reset km, flagged file")
		outlier <- data.frame(day, file=getFileNumber(file))
		outlier.table <- rbind(outlier.table, outlier)
		flag <- 1
		}
		opp[row.names(x),'pop'] <- km$cluster + length(pop.baseline)
	opp.b <- opp

		for(i in 1:(numc1+numc2)){
			chl <- names(sort(prev.km$centers[,"chl_small"]))[i]
			df <- subset(opp.b, pop == (as.numeric(chl)+length(pop.baseline)))
				if(median(df$chl_small) > lim) opp[row.names(df),'pop'] <- "y"
				else opp[row.names(df),'pop'] <- "z"
					}
				
	# cluster cells larger than Synecho OR 1 um Beads
	y <- subset(opp, pop == "y")
	y <- y[,c(5,8:10)]
	if(nrow(y) > 30 & numc1 > 1){
	km.big <- try(kmeans(y, prev.km.big$centers,iter.max))
	if(class(km.big) == "try-error"){
		km.big <- kmeans(y,numc1 ,iter.max); km.big <- kmeans(y, km.big$centers, iter.max)
		print("reset km.big, flagged file")
		outlier <- data.frame(day, file=getFileNumber(file))
		outlier.table <- rbind(outlier.table, outlier)
		flag <- 1
			}
	opp[row.names(y),'pop'] <- km.big$cluster + length(pop.baseline)
	opp.b <- opp

		for(i in 1:numc1){
			fsc <- names(sort(km.big$centers[,"fsc_small"]))[i]
			df <- subset(opp.b, pop == (as.numeric(fsc)+length(pop.baseline)))
			opp[row.names(df),'pop'] <- i + length(pop.baseline)
			}
		}else{	
			print("not enough cells for clustering ultra-nanoplankton")
			# flag <- 1
			outlier <- data.frame(day, file=getFileNumber(file))
			outlier.table <- rbind(outlier.table, outlier)
				} 
	
	# cluster cells smaller than Synecho OR 1 um Beads
	z <- subset(opp, pop == "z")
	z <- z[,c(5,8:10)]
	if(nrow(z) > 30 & numc2 > 1){
		km.small <- try(kmeans(z, prev.km.small$centers,iter.max))
			if(class(km.small) == "try-error"){
				km.small <- kmeans(z,numc2,iter.max);km.small <- kmeans(z, km.small$centers, iter.max)
				print("reset km.small, flagged file")
				outlier <- data.frame(day, file=getFileNumber(file))
				outlier.table <- rbind(outlier.table, outlier)
				flag <- 1
					}
		opp[row.names(z),'pop'] <- km.small$cluster + numc1 + length(pop.baseline)
	
			for(i in 1:(numc2)){
				df <- subset(opp, pop == (i+numc1+length(pop.baseline)))
				quant.fsc <- quantile(df$fsc_small, probs=0.9) - quantile(df$fsc_small, probs=0.1)
				quant.chl <- quantile(df$chl_small, probs=0.25)
				if(quant.chl < lim2){
								opp[row.names(df),'pop'] <- 1
								# print(paste("pop", i+numc1+length(pop.baseline), "converts to", 1))
								# opp[row.names(subset(df, chl_small > lim2)), 'pop'] <- numc1 + i + length(pop.baseline)

						}else{
								opp[row.names(df),'pop'] <- numc1 + i + length(pop.baseline)
								#print(paste("pop", i+numc1+length(pop.baseline), "converts to", numc1 + i + length(pop.baseline)))
						}
				}
		
			}else{
				print("not enough cells for clustering picoplankton")
				# flag <- 1
				outlier <- data.frame(day, file=getFileNumber(file))
				outlier.table <- rbind(outlier.table, outlier)
				}
			
			

	# ##########################				
	# ###	K-means validation ###
	# ##########################
	if(flag == 0 & any(prev.km.big$centers[,"chl_small"]/km.big$centers[,"chl_small"] > 1 * error) |any(prev.km.big$centers[,"chl_small"]/km.big$centers[,"chl_small"] < 1 / error) |any(prev.km.big$centers[,"fsc_small"]/km.big$centers[,"fsc_small"] > 1 * error) |any(prev.km.big$centers[,"fsc_small"]/km.big$centers[,"fsc_small"] < 1 / error) |any(prev.km.big$size/km.big$size > diff) |any(prev.km.big$size/km.big$size < 1/diff) |any(prev.km.small$centers[,"chl_small"]/km.small$centers[,"chl_small"] > 1 * error) |any(prev.km.small$centers[,"chl_small"]/km.small$centers[,"chl_small"] < 1 / error) |any(prev.km.small$centers[,"fsc_small"]/km.small$centers[,"fsc_small"] > 1 * error) |any(prev.km.small$centers[,"fsc_small"]/km.small$centers[,"fsc_small"] < 1 / error) |any(prev.km.small$size/km.small$size > diff) |any(prev.km.small$size/km.small$size < 1/diff)){
0				
			print(paste("file",basename(file),"flagged"))

			flag <- 1
			
			outlier <- data.frame(day, file=getFileNumber(file))
			outlier.table <- rbind(outlier.table, outlier)
						
			}else{
			
				# save kmeans output
				
				prev.km.big <- km.big 
				prev.km.small <- km.small 
				prev.km <- km 
				save(prev.km,file= paste(save.path,"prev.km.RData",sep=""))
				save(prev.km.big,file= paste(save.path,"prev.km.big.RData",sep=""))
				save(prev.km.small,file= paste(save.path,"prev.km.small.RData",sep=""))
				save(prev.synecho,file= paste(save.path,"prev.synecho.RData",sep=""))

			}		


	write.csv(outlier.table, file=paste(save.path,"c.outliers", sep=""), row.names=FALSE, quote=FALSE)
		
				
	#################
	### class vct ###
	#################

	for(i in 1:(length(phyto))){
		p <- subset(opp, pop == i)
		opp[row.names(p),'pop'] <- phyto[i]
			}

	write.table(opp$pop, paste(save.path, day,"/",basename(file),".",getFileNumber(file),'-class.vct',sep=""), row.names=FALSE, col.names='pop', quote=FALSE)
	
	
	
	
	############
	### PLOT ###
	############
	if(do.plot==TRUE){

	png(paste(save.path, day,"/",basename(file),".",getFileNumber(file),".class.gif", sep=""),width=9, height=12, unit='in', res=100)
	if(class(opp) != "try-error"){

	hist1 <- hist(opp$fsc_small, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist2 <- hist(opp$chl_small, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist3 <- hist(opp$pe, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist4 <- hist(opp$chl_big, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)
	hist5 <- hist(opp$fsc_perp, breaks=seq(0,2^16, by=2^16/breaks), plot=FALSE)

	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),6,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)
	
	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'chl_small', pop.def=pop.def)
		
		if(flag == 1) mtext("flagged !", side=1, line=4,cex=3,col='red', at=100000)
		
		if(getFileNumber(file)-getFileNumber(prev.file) > 1) mtext("discontinued file !", side=1, line=7,cex=3,col='red', at=100000)

	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'chl_big', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'chl_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'fsc_perp', pop.def=pop.def, add.legend=TRUE)
	mtext(paste(file), side=1, line=5, at=-30000)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist5$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(def.par)
	
	dev.off()
		}
	
	}
	
	
	###########
	
	prev.file <- file
	
	
			}
		}
	
	
	
	
	}
}