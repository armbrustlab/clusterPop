.nextplot <- function(){
			answer <- readline("Continue ? (y/n)")
			if(answer == 'y') x <- 1
			if(answer == 'n') x <- 0
			return(x)
			}

clusterPopBaseline <- function(opp, pop.def, noise){
	
	#Cluster Beads
	x <- subset(opp, pop=='x')
	beads <- subset(x, pe > 65500) ## PE-saturated particles assigned as 'beads
	opp[row.names(beads),'pop'] <- 'beads'
	
	x <- subset(opp, pop=='x')
	xvar <- pop.def["beads", "xvar"]
	yvar <- pop.def["beads", "yvar"]
	beads <- subset(x, x[,yvar] > 0.5*x[,xvar] + pop.def["beads", "lim"] & x[,xvar] > pop.def["beads", "xmin"] & x[,yvar] > pop.def["beads", "ymin"] & x[,xvar] < pop.def["beads", "xmax"] & x[,yvar] < pop.def["beads", "ymax"])
	opp[row.names(beads),'pop'] <- 'beads'

	#Cluster Synecho
	x <- subset(opp, pop=='x')
	yvar <- pop.def["synecho", "yvar"]
	xvar <- pop.def["synecho", "xvar"]
	synecho <- subset(x, x[,yvar] > x[,xvar] - pop.def["synecho", "lim"] & x[,xvar] > pop.def["synecho", "xmin"] & x[,yvar] > pop.def["synecho", "ymin"] & x[,xvar] < pop.def["synecho", "xmax"] & x[,yvar] < pop.def["synecho", "ymax"])
	opp[row.names(synecho), 'pop'] <- 'synecho'
	
	#Cluster Cryptophytes
	x <- subset(opp, pop=='x')
	yvar <- pop.def["crypto", "yvar"]
	xvar <- pop.def["crypto", "xvar"]
	crypto <- subset(x, x[,yvar] > x[,xvar] + pop.def["crypto", "lim"] & x[,xvar] > pop.def["crypto", "xmin"] & x[,yvar] > pop.def["crypto", "ymin"])
	opp[row.names(crypto), 'pop'] <- 'crypto'

	#Cluster Coccolithophores
	x <- subset(opp, pop=='x')
	yvar <- pop.def["cocco", "yvar"] 
	xvar <- pop.def["cocco", "xvar"]
	cocco <- subset(x, x[,yvar] > x[,xvar] + pop.def["cocco", "lim"] & x[,xvar] > pop.def["cocco", "xmin"] & x[,yvar] > pop.def["cocco", "ymin"]& 	x[,"chl_small"] > quantile(synecho$chl_small,0.9))
	opp[row.names(cocco), 'pop'] <- 'cocco'

	
	#Cluster Noise
	if(!is.null(noise)){
		x <- subset(opp, pop=='x')
		bg <- subset(x, chl_small < noise[1] | fsc_small < noise[2])
		opp[row.names(bg),'pop'] <- 'noise'
		}
	
	
	return(opp)
		
}

clusterPENegative <- function(opp, n.pop, lim, lim.debris, tol, h0, h, test=F){
	
	##################################################################
	### Separate phytoplankton into 'large' and 'small' populations ###
	##################################################################
	x <- subset(opp, pop == 'x')
	x <- x[,c(5,9)]

	km <- flowPeaks(x, tol=tol, h0=h0, h=h)
	
	opp[row.names(x),'pop'] <- km$peaks.cluster 
	
	p1 <- median(opp[row.names(x),"chl_small"] / opp[row.names(x),"fsc_small"]) * 1.25
	
	if(p1 < 0.8) p1 <- 1 
		for(i in 1:((max(km$peaks.cluster)))){
			df <- subset(opp, pop == i)
				if(median(df$chl_small/df$fsc_small) > p1) opp[row.names(df),'pop'] <- 'pennates'	
				if(median(df$chl_small/df$fsc_small) < p1 & median(df$chl_small) > lim) opp[row.names(df),'pop'] <- 'y'
				if(median(df$chl_small/df$fsc_small) < p1 & median(df$chl_small) < lim) opp[row.names(df),'pop'] <- 'z'
			}


	#####################################################################
	### cluster cells larger than 1 um Beads or Synecho (if no beads) ###
	#####################################################################
	y <- subset(opp, pop == 'y' & chl_small > 50000 & fsc_small > 50000) # Chl-saturating particles and largest cells assigned as 'nano'
	opp[row.names(y),'pop'] <- 'nano'

	y <- subset(opp, pop == "y")
	y <- y[,c(5,9)]

	km.big <- try(flowMeans(y, NumC=n.pop[1],MaxN=n.pop[1]*2, Standardize=F, Update='None'))
	if(class(km.big) == 'try-error')	km.big <- flowMeans(y, NumC=n.pop[1],MaxN=n.pop[1]+2, Standardize=F, Update='None')
	if(class(km.big) == 'try-error')	km.big <- flowMeans(y, NumC=n.pop[1],MaxN=n.pop[1]+1, Standardize=F, Update='None')

	opp[row.names(y),'pop'] <- km.big@Label 
	
	fsc.med <- by(opp[as.numeric(opp$pop) > 0,"fsc_small"], opp[as.numeric(opp$pop) > 0,"pop"], median)
	fsc.sort <- names(sort(fsc.med))
	
	if(n.pop[1]  == 2){
				df <- subset(opp, pop == as.numeric(fsc.sort[2])); opp[row.names(df),'pop'] <- 'nano'
				df <- subset(opp, pop == as.numeric(fsc.sort[1])); opp[row.names(df),'pop'] <- 'ultra'
				}
	
	if(n.pop[1]  == 3){
			
			chl.med <- by(opp[as.numeric(opp$pop) > 0,"chl_small"], opp[as.numeric(opp$pop) > 0,"pop"], median)
			diato <- chl.med/fsc.med
			p2 <- 	median(opp[as.numeric(opp$pop) > 0,"chl_small"] / opp[as.numeric(opp$pop) > 0,"fsc_small"], na.rm=T) * 1.25
				for( i in 1:3){	
							df <- subset(opp, pop == as.numeric(names(diato[i])))
							if (diato[i] > p2) opp[row.names(df),'pop'] <- 'pennates'
							if (diato[i] < p2) opp[row.names(df),'pop'] <- 'y'
									}
																
			w <- subset(opp, pop == "y")
			w <- w[,c(5,9)]

			km.big2 <- try(flowMeans(exp(w/2^14), NumC=2,MaxN=2*2, Standardize=F, Update='None'))
			if(class(km.big2) == 'try-error')	km.big <- flowMeans(exp(w/2^14), NumC=2,MaxN=n.pop[1]+2, Standardize=F, Update='None')
			if(class(km.big2) == 'try-error')	km.big <- flowMeans(exp(w/2^14), NumC=2,MaxN=n.pop[1]+1, Standardize=F, Update='None')
			
			opp[row.names(w),'pop'] <- km.big2@Label 
	
			fsc.med <- by(opp[as.numeric(opp$pop) > 0,"fsc_small"], opp[as.numeric(opp$pop) > 0,"pop"], median)
			fsc.sort <- names(sort(fsc.med))
	
			df <- subset(opp, pop == as.numeric(fsc.sort[2])); opp[row.names(df),'pop'] <- 'nano'
			df <- subset(opp, pop == as.numeric(fsc.sort[1])); opp[row.names(df),'pop'] <- 'ultra'
			
		}
				
	######################################################################
	### cluster cells smaller than 1 um Beads or Synecho (if no beads) ###
	######################################################################
	
	z <- subset(opp, pop == "z")
	z <- z[,c(5,9,9)]

	km.small <- flowMeans(z, NumC=n.pop[2],MaxN=n.pop[2]*2, Standardize=F, Update='None')
	
	opp[row.names(z),'pop'] <- km.small@Label 
	
		for(i in 1:as.numeric(n.pop[2])){
			df <- subset(opp, pop == i)
			if(nrow(df) > 2){
				d.quant <- peaks(density(df[,"chl_small"],from=0, to=2^16, n=2^16))
				width.chl <- round(d.quant[which(d.quant$x == max(d.quant$x)),"w"],3)
				quant.chl <- round(d.quant[which(d.quant$x == max(d.quant$x)),"x"],3)
				if(quant.chl < lim.debris[1]){
							opp[row.names(df),'pop'] <- 'noise'
							print(paste("LOW CHL,",quant.chl, "pop",i, "converted to Noise"))	
							}
				if(quant.chl > lim.debris[1] & width.chl < lim.debris[2]){	
							opp[row.names(df),'pop'] <- 'noise'
							print(paste("Narrow CHL width,",width.chl, "pop",i, "converted to Noise"))	
							}
				if(quant.chl > lim.debris[1] & width.chl > lim.debris[2]){
							opp[row.names(df),'pop'] <- 'pico'
					}
			}else{opp[row.names(df),'pop'] <- 'noise'
					print(paste("Only 1 cell found, pop",i, "converted to Noise"))	
					}
		}
	
	
	if(test == T){
		####################
		### PLOT CLUSTER ###
		####################
			breaks <- 50
				
			hist1 <- hist(x$fsc_small, breaks=breaks, plot=FALSE)
			hist2 <- hist(x$chl_small, breaks=breaks, plot=FALSE)
			hist3 <- hist(y$fsc_small, breaks=breaks, plot=FALSE)
			hist4 <- hist(y$chl_small, breaks=breaks, plot=FALSE)				
			hist5 <- hist(w$fsc_small, breaks=breaks, plot=FALSE)
			hist6 <- hist(w$chl_small, breaks=breaks, plot=FALSE)
			hist7 <- hist(y$fsc_small, breaks=breaks, plot=FALSE)
			hist8 <- hist(y$chl_small, breaks=breaks, plot=FALSE)				
		
			def.par <- par(no.readonly = TRUE) # save default, for resetting...
			nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),4,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)	

			par(mar=c(6,6,1,1))
			plot(km,drawboundary=T,drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)
			mtext("Cut-off Large and Small Phytoplankton", 3,cex=0.7)
			abline(h=lim, col='red',lwd=4)
			par(mar=c(0,6,1,1))
			barplot(hist1$counts, axes=FALSE, space=0, col=NA)
			par(mar=c(6,0,1,1))
			barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
	
			par(mar=c(6,6,1,1))
			plot(y[,"fsc_small"],y[,"chl_small"], col=km.big@Label,cex=0.5,xlab='fsc_small', ylab='chl_small', pch=16)
			mtext("Large Phytoplankton", 3, cex=0.7)
			par(mar=c(0,6,1,1))
			barplot(hist3$counts, axes=FALSE, space=0, col=NA)
			par(mar=c(6,0,1,1))
			barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
		
			par(mar=c(6,6,1,1))
			plot(w[,"fsc_small"],w[,"chl_small"], col=km.big2@Label,cex=0.5,xlab='fsc_small', ylab='chl_small', pch=16)
			mtext("Large Phytoplankton", 3, cex=0.7)
			par(mar=c(0,6,1,1))
			barplot(hist5$counts, axes=FALSE, space=0, col=NA)
			par(mar=c(6,0,1,1))
			barplot(hist6$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
		
			par(mar=c(6,6,1,1))
			plot(z[,"fsc_small"],z[,"chl_small"], col=km.small@Label,cex=0.5,xlab='fsc_small', ylab='chl_small', pch=16)
			mtext("Small Phytoplankton", 3, cex=0.7)
			par(mar=c(0,6,1,1))
			barplot(hist7$counts, axes=FALSE, space=0, col=NA)
			par(mar=c(6,0,1,1))
			barplot(hist8$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

			par(def.par)
	
		}
	
		
		
	return(opp)
}

clusterPop <- function(opp.filelist, save.path = getCruisePath(opp.filelist), pop.def.path, test=F,  noise = c(0,0), n.pop = c(3,3), tol=0.1,h0=0.5, h=1.0, min.opp = 800, min.beads = 200, min.synecho=200, lim.beads = 0, lim.synecho =0, lim.debris=c(10000,4000), save.plot=F){
	

#library(cluster)	
require(flowPhyto)
require(flowMeans)
require(flowPeaks)
require(caroline)
require(IDPmisc)
	
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))
	
# read pop.def.tab
pop.def <- readPopDef(pop.def.path) 
	
	
prev.file <- opp.filelist[1]

#################################################
### check that there is no cluster settings already ###
#################################################

if(test){

###################
###################
#### prev FILE ###
###################	
###################
	
	print(paste("Initialization", prev.file))
	
	opp <- readSeaflow(prev.file) 
	opp$pop <- "x"

	lim <- 0.25*2^16 

	###############################
	### clustering pop.baseline ###
	###############################
	
 	opp <- clusterPopBaseline(opp, pop.def, noise)
	
	beads <- subset(opp, pop == "beads")
	synecho <- subset(opp, pop == "synecho")

	if(nrow(beads) > min.beads){
		lim <- quantile(beads$chl_small,0.1) + lim.beads
	print(paste("cut-off between particles 'smaller' and 'larger' than Beads:", lim))
		}

	if(nrow(beads) < min.beads & nrow(synecho) > min.synecho){
		lim <- median(synecho$chl_small) + lim.synecho
	print(paste("cut-off between particles 'smaller' and 'larger' than Synecho:", lim))
		}	

	###############################################
	#### clustering pe-negative phytoplankton  ####
	###############################################

	opp <- clusterPENegative(opp, n.pop, lim, lim.debris, tol, h0, h, test=test)

	go <- .nextplot()	
	
	if(go == 1){

	######################
	### SAVE SETTINGS  ###
	######################
	
	all.settings <- aggregate(opp[,c("fsc_small", "chl_small")], FUN=median, by=list(opp$pop))		
	settings <- all.settings[all.settings$Group.1 == "nano" | all.settings$Group.1 == "ultra" | all.settings$Group.1 == "synecho",]
	write.delim(settings, paste(save.path,"prev.settings",sep=""))



	##################
	### PLOT PHYTO ###
	##################
	breaks <- 50
		
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

	}
	
}

################################################
################################################
##### ALL CRUISES FILES based on prev FILE ####
################################################
################################################
if(test == F){ 
	lim <- 0.5*2^16

for (file in opp.filelist){
	opp <- readSeaflow(file)
	opp$pop <- 'x'

	day <- flowPhyto:::.getYearDay(file)
	
	##########################################
	### Check number OPP before Clustering ### 
	##########################################
	
	if(nrow(opp[opp[,"chl_small"]>20000,]) < min.opp){
		print(paste("not enough OPP for clustering in:",file))
		opp$pop <- 'noise'
		write.table(opp$pop, paste(save.path, day,"/",basename(file),".",getFileNumber(file),'-class.vct',sep=""), row.names=FALSE, col.names='pop', quote=FALSE)
			
		png(paste(save.path, day,"/",basename(file),".",getFileNumber(file),".outlier.class.gif", sep=""),width=9, height=12, unit='in', res=100)
	
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

 	opp <- clusterPopBaseline(opp, pop.def, noise)

	beads <- subset(opp, pop == "beads")
	synecho <- subset(opp, pop == "synecho")

	if(nrow(beads) > min.beads){
		lim <- quantile(beads$chl_small,0.1) + lim.beads
	print(paste("cut-off between particles 'smaller' and 'larger' than Beads:", lim))
		}

	if(nrow(beads) < min.beads & nrow(synecho) > min.synecho){
		lim <- median(synecho$chl_small) + lim.synecho
	print(paste("cut-off between particles 'smaller' and 'larger' than Synecho:", lim))
		}
				
	###############################################
	#### clustering pe-negative phytoplankton  ####
	###############################################
	
	opp1 <- clusterPENegative(opp, n.pop, lim, lim.debris, tol, h0, h, test=F)	
	if(class(opp1) == "try-error"){
				print("Error in clustering PE negative cells, FLAGGED FILE")
				opp$pop <- 'noise'
				write.table(opp$pop, paste(save.path, day,"/",basename(file),".",getFileNumber(file),'-class.vct',sep=""), row.names=FALSE, col.names='pop', quote=FALSE)

				next
				}
	
		
					
	#################
	### class vct ###
	#################
	opp <- opp1
	write.table(opp$pop, paste(save.path, day,"/",basename(file),".",getFileNumber(file),'-class.vct',sep=""), row.names=FALSE, col.names='pop', quote=FALSE)
	
	
	
	############
	### PLOT ###
	############

	if(save.plot==TRUE){

	png(paste(save.path, day,"/",basename(file),".",getFileNumber(file),".class.gif", sep=""),width=9, height=12, unit='in', res=100)
	if(class(opp) != "try-error"){
	
	breaks <- 50

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
	### EMPTY MEMORY ###	
	####################	

	rm(opp, noise, synecho, crypto, beads, cocco, hetero, p, df, opp1, hist1, hist2, hist3, hist4)	
	}

	print("Combining census.tab ...")
	caroline:::write.delim(combineCensusFiles(save.path),paste(save.path,'census.tab',sep=""))
	
	
	}		

}