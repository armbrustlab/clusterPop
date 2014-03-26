.nextplot <- function(){
			answer <- readline("Continue ? (y/n)")
			if(answer == 'y') x <- 1
			if(answer == 'n') x <- 0
			return(x)
			}


.clusterPopBaseline <- function(opp, pop.def, lim.largeP,...){
	
	#Cluster Beads
	x <- subset(opp, pop=='x')
	beads <- subset(x, pe == max(opp$pe, na.rm=T)) ## PE-saturated particles assigned as 'beads
	opp[row.names(beads),'pop'] <- 'beads'
	
	x <- subset(opp, pop=='x')
	xvar <- pop.def["beads", "xvar"]
	yvar <- pop.def["beads", "yvar"]
	beads <- subset(x, x[,yvar] > 0.5*x[,xvar] + pop.def["beads", "lim"] & x[,xvar] > pop.def["beads", "xmin"] & x[,yvar] > pop.def["beads", "ymin"] & x[,xvar] < pop.def["beads", "xmax"] & x[,yvar] < pop.def["beads", "ymax"])
	opp[row.names(beads),'pop'] <- 'beads'

	#Cluster Synechococcus
	x <- subset(opp, pop=='x')
	yvar <- pop.def["synecho", "yvar"]
	xvar <- pop.def["synecho", "xvar"]
	synecho <- subset(x, x[,yvar] > x[,xvar] - pop.def["synecho", "lim"] & x[,xvar] > pop.def["synecho", "xmin"] & x[,yvar] > pop.def["synecho", "ymin"] & x[,xvar] < pop.def["synecho", "xmax"] & x[,yvar] < pop.def["synecho", "ymax"])
	opp[row.names(synecho), 'pop'] <- 'synecho'
	
 	#Cluster Coccolithophores
	x <- subset(opp, pop=='x')
	yvar <- pop.def["cocco", "yvar"] 
	xvar <- pop.def["cocco", "xvar"]
	cocco <- subset(x, x[,yvar] > x[,xvar] + pop.def["cocco", "lim"] & x[,xvar] > pop.def["cocco", "xmin"] & x[,yvar] > pop.def["cocco", "ymin"])
	opp[row.names(cocco), 'pop'] <- 'cocco'
	
	#Cluster Cryptophytes
	x <- subset(opp, pop=='x')
	yvar <- pop.def["crypto", "yvar"]
	xvar <- pop.def["crypto", "xvar"]
	crypto <- subset(x, x[,yvar] > x[,xvar] + pop.def["crypto", "lim"] & x[,xvar] > pop.def["crypto", "xmin"] & x[,yvar] > pop.def["crypto", "ymin"])
	opp[row.names(crypto), 'pop'] <- 'crypto'


	#Cluster Large Particles
	if(!is.null(lim.largeP)){
		x <- subset(opp, pop=='x')
		hetero <- subset(x, x[,"fsc_small"] > 1.5*x[,"chl_small"] + lim.largeP)
		opp[row.names(hetero), 'pop'] <- 'noise'
		}


	
# # 	#Cluster Noise
	# if(!is.null(noise)){
		# x <- subset(opp, pop=='x')
		# bg <- subset(x, chl_small < noise[1] | fsc_small < noise[2])
		# opp[row.names(bg),'pop'] <- 'noise'
		# }
	
	
	return(opp)
		
}

.clusterPENegative <- function(opp, lim, lim.debris, tol, h0, h, test=F,...){
	
	##################################################################
	### Separate phytoplankton into 'large' and 'small' populations ###
	##################################################################
	x <- subset(opp, pop == 'x')
	x <- x[,c(5,9)]
	
	km <- suppressMessages(flowPeaks(x, tol=tol, h0=h0, h=h))
	
	opp[row.names(x),'pop'] <- km$peaks.cluster 
	
	p1 <- median(opp[row.names(x),"chl_small"] / opp[row.names(x),"fsc_small"]) * 1.25
	if(p1 < 1.1) p1 <- 1.1 
			
		for(i in 1:((max(km$peaks.cluster)))){
			df <- subset(opp, pop == i)
				if(median(df$chl_small) < lim) opp[row.names(df),'pop'] <- 'z'
				if(median(df$chl_small/df$fsc_small) > p1 & median(df$chl_small) > lim) opp[row.names(df),'pop'] <- 'pennates'
				if(median(df$chl_small/df$fsc_small) < p1 & median(df$chl_small) > lim) opp[row.names(df),'pop'] <- 'y'
			}


	#####################################################################
	### cluster cells larger than 1 um Beads or Synecho (if no beads) ###
	#####################################################################
	y <- subset(opp[which(opp$pop == 'y' | opp$pop == 'pennates'),], chl_small > 50000 & fsc_small > 50000) # High Chl- and Fsc-particles assigned as 'nano'
	opp[row.names(y),'pop'] <- 'nano'

	y <- subset(opp, pop == "y")
	y <- y[,c(5,5,5,9)]
	
	if(nrow(y) > 10){
		km.big <- suppressMessages(flowMeans(y, NumC=3,Standardize=F, Update='Mahalanobis'))
	
		opp[row.names(y),'pop'] <- km.big@Label
	
		fsc.med <- by(opp[rownames(y),"fsc_small"], opp[rownames(y),"pop"], median); fsc.sort <- names(sort(fsc.med))
		chl.med <- by(opp[rownames(y),"chl_small"], opp[rownames(y),"pop"], median); diato <- chl.med/fsc.med
		p2 <- median(opp[rownames(y),"chl_small"] / opp[rownames(y),"fsc_small"], na.rm=T) * 1.25; if(p2 < 1.1) p2 <- 1.1
	
				for( i in 1:3){	
							df <- subset(opp, pop == as.numeric(names(diato[i])))
							if (diato[i] > p2)opp[row.names(df),'pop'] <- 'pennates'
							if (diato[i] < p2) opp[row.names(df),'pop'] <- 'y'
									}
	
	w <- subset(opp, pop == "y")
	w <- w[,c(5,5,5,9)]

	if(nrow(w) > 10){
		km.big2 <- flowMeans(w, NumC=2, Standardize=F, Update='Mahalanobis')
			
		opp[row.names(w),'pop'] <- km.big2@Label
	
		fsc.med <- by(opp[rownames(w),"fsc_small"], opp[rownames(w),"pop"], median)
		fsc.sort <- names(sort(fsc.med, decreasing =T))
			
		df <- subset(opp, pop == as.numeric(fsc.sort[1]))
		opp[row.names(df),'pop'] <- 'nano'
			
		df <- subset(opp, pop == as.numeric(fsc.sort[2]))
		opp[row.names(df),'pop'] <- 'ultra'
		}
	}
		
	######################################################################
	### cluster cells smaller than 1 um Beads or Synecho (if no beads) ###
	######################################################################
	
	z <- subset(opp, pop == "z")
	z <- z[,c(5,9)]
	
	if(nrow(z) > 10){

	km.small <- suppressMessages(flowPeaks(z, tol=tol*1.5, h0=h0, h=h))
	
	opp[row.names(z),'pop'] <- km.small$peaks.cluster
		
		for(i in unique(km.small$peaks.cluster)){
			df <- subset(opp, pop == i)
			med.chl <- median(df$chl_small)
			width.chl <- diff(c(quantile(df$chl_small, 0.05),quantile(df$chl_small, 0.95)))
			print(paste("pop",i, ": MED,",round(med.chl) , "& WIDTH", round(width.chl)), quote=F)	
				if(all(c(med.chl  > lim.debris[1], width.chl > lim.debris[2]))){
								opp[row.names(df),'pop'] <- 'pico'
							}else{
								opp[row.names(df),'pop'] <- 'noise'
							    print(paste("... converted to Noise"), quote=F)	
								}
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
			hist7 <- hist(z$fsc_small, breaks=breaks, plot=FALSE)
			hist8 <- hist(z$chl_small, breaks=breaks, plot=FALSE)				
		
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
			mtext("Looking for remaining Pennates in Large Phytoplankton", 3, cex=0.7)
			par(mar=c(0,6,1,1))
			barplot(hist3$counts, axes=FALSE, space=0, col=NA)
			par(mar=c(6,0,1,1))
			barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

			par(mar=c(6,6,1,1))
			plot(w[,"fsc_small"],w[,"chl_small"], col=km.big2@Label,cex=0.5,xlab='fsc_small', ylab='chl_small', pch=16)
			mtext("Clustering Large Phytoplankton in 2 populations", 3, cex=0.7)
			par(mar=c(0,6,1,1))
			barplot(hist5$counts, axes=FALSE, space=0, col=NA)
			par(mar=c(6,0,1,1))
			barplot(hist6$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)
		
			par(mar=c(6,6,1,1))
			plot(km.small,drawboundary=T,drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)
			mtext("Isolating Small Phytoplankton from debris", 3, cex=0.7)
			par(mar=c(0,6,1,1))
			barplot(hist7$counts, axes=FALSE, space=0, col=NA)
			par(mar=c(6,0,1,1))
			barplot(hist8$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

			par(def.par)
	
		}
	
		
		
	return(opp)
}

clusterPop <- function(opp.filelist, save.path = getCruisePath(opp.filelist), pop.def.path, test=F, lim.largeP = 0.5*2^16, tol=0.1, h0=0.5, h=1.0, min.opp = 800, min.beads = 200, min.synecho=200, lim.beads = 0, lim.synecho =0, lim.debris=c(10000,4000), save.plot=F,...){
	

require(flowPhyto)
require(flowMeans)
require(flowPeaks)
	
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
	
	message(paste("Initialization", prev.file))
	
	opp <- readSeaflow(prev.file, transform=F) 
	opp$pop <- "x"

	lim <- 0.25*2^16 

	###############################
	### clustering pop.baseline ###
	###############################
	
 	opp <- .clusterPopBaseline(opp, pop.def, lim.largeP)
	
	beads <- subset(opp, pop == "beads")
	synecho <- subset(opp, pop == "synecho")

	if(nrow(beads) > min.beads){
		lim <- quantile(beads$chl_small,0.1) + lim.beads
	print(paste("cut-off between particles 'smaller' and 'larger' than Beads:", lim), quote=F)
		}

	if(nrow(beads) < min.beads & nrow(synecho) > min.synecho){
		lim <- median(synecho$chl_small) + lim.synecho
	print(paste("cut-off between particles 'smaller' and 'larger' than Synecho:", lim), quote=F)
		}	

	###############################################
	#### clustering pe-negative phytoplankton  ####
	###############################################

	opp <- .clusterPENegative(opp, lim, lim.debris, tol, h0, h, test=test)

	go <- .nextplot()	
	
	if(go == 1){

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
	
	message(paste("clustering", file))

	opp <- readSeaflow(file)
	opp$pop <- 'x'

	day <- flowPhyto:::.getYearDay(file)
	
	##########################################
	### Check number OPP before Clustering ### 
	##########################################
	
	if(nrow(opp[opp[,"chl_small"]>20000,]) < min.opp){
		print(paste("not enough OPP for clustering in:",file), quote=F)
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

 	opp <- .clusterPopBaseline(opp, pop.def, lim.largeP)

	beads <- subset(opp, pop == "beads")
	synecho <- subset(opp, pop == "synecho")

	if(nrow(beads) > min.beads){
		lim <- quantile(beads$chl_small,0.1) + lim.beads
	print(paste("cut-off between particles 'smaller' and 'larger' than Beads:", lim), quote=F)
		}

	if(nrow(beads) < min.beads & nrow(synecho) > min.synecho){
		lim <- median(synecho$chl_small) + lim.synecho
	print(paste("cut-off between particles 'smaller' and 'larger' than Synecho:", lim), quote=F)
		}
				
	###############################################
	#### clustering pe-negative phytoplankton  ####
	###############################################
	
	opp1 <- try(.clusterPENegative(opp, lim, lim.debris, tol, h0, h, test=F))	
	if(class(opp1) == "try-error"){
				print("Error in clustering PE negative cells, FLAGGED FILE", quote=F)
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
	
	
	####################
	### EMPTY MEMORY ###	
	####################	

	rm(opp, opp1, synecho, beads, hist1, hist2, hist3, hist4)	
	}

	}		

}