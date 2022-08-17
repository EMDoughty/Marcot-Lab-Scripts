
plotRateShiftsBM <- function(this.rez, this.tree, this.alpha=0.5, n.increments=10) {
	# max.rates.BM <- max(which(!sapply(lapply(this.rez, function(x) x$BM), is.null))) 
	this.break.dates <- c(max(this.tree$node.date), this.rez$BM$break.dates, min(this.tree$node.date))
	inv.break.dates <- max(this.tree$node.date) - this.break.dates
	x.left <- par()$usr[1]
	x.right <- par()$usr[2]
	heat.increment <- max(this.rez$BM$sigma) / n.increments
	heat.values <- (as.vector(round(this.rez$BM$sigma/heat.increment, 0))) + 1
	this.colors <- c(adjustcolor("white", alpha.f=this.alpha), rev(heat.colors(n=n.increments, alpha=this.alpha)))
	for (i in seq_len(length(inv.break.dates)-1)) {
		rect(xleft = x.left, ybottom = inv.break.dates[i], xright=x.right, ytop = inv.break.dates[i+1], col=this.colors[heat.values[i]], border=NA)
		text(x=x.left, y=mean(inv.break.dates[c(i, i+1)]), labels=round(this.rez$BM$sigma[[i]], 3), pos = 4)
		text(x=x.left, y=inv.break.dates[i], labels=paste0(round(this.break.dates[i], 2), "Ma"), pos = 4, cex=0.5, col="firebrick4")
	}
}

plotRateShiftsOU <- function(this.rez, this.tree, this.alpha=0.5, n.increments=10) {
	# max.rates.BM <- max(which(!sapply(lapply(this.rez, function(x) x$BM), is.null))) 
	this.break.dates <- c(max(this.tree$node.date), this.rez$OU$break.dates, min(this.tree$node.date))
	inv.break.dates <- max(this.tree$node.date) - this.break.dates
	x.left <- par()$usr[1]
	x.right <- par()$usr[2]
	heat.increment <- max(this.rez$OU$theta) / n.increments
	heat.values <- (as.vector(round(this.rez$OU$theta/heat.increment, 0))) + 1
	this.colors <- c(adjustcolor("white", alpha.f=this.alpha), rev(heat.colors(n=n.increments, alpha=this.alpha)))
	for (i in seq_len(length(inv.break.dates)-1)+1) {
		rect(xleft = x.left, ybottom = inv.break.dates[i-1], xright=x.right, ytop = inv.break.dates[i], col=this.colors[heat.values[i]], border=NA)
		text(x=x.left, y=mean(inv.break.dates[c(i, i-1)]), labels=round(this.rez$OU$theta[[i]], 3), pos = 4)
		text(x=x.left, y=inv.break.dates[i], labels=paste0(round(this.break.dates[i], 2), "Ma"), pos = 4, cex=0.5, col="firebrick4")
	}
}

evoModelHists <- function(evoResults, intervals, runHistOn = c("AllRates"), model=c("BMOU"), plotPercent = FALSE, ylim = NULL, relativeFreq = NULL) {
	
	quartz(width=12, height = 6)
	par(mfrow=c(1,1), mar=c(4,3,1.5,2.5), mgp=c(2, 1,0), oma=c(0,0,0,4), xaxs="i")
	# optList_tax_allReps <- optList_tax_allReps_heuristic
	# optList_tax_allReps <- optList_tax_allReps_full
	
	#if(class(intervals) == "list"){ #will compile the intervals from a list if running multiple trees so that no bins are being left out due to younger node dates
	#	int.Top <- vector()
	#	int.Base <- vector()
	#	for(kk in seq(1, length(intervals))){
	#		int.Top <- append(int.Top, intervals[[kk]]$ageTop)
	#		int.Base <- append(int.Base, intervals[[kk]]$ageBase)
	#	}
	#	intervals <- unique(cbind(ageTop = int.Top, ageBase = int.Base))
	#}
	
	#may need to set the max date of the interval vector to be the smallest node date in the tree.list so it can be applied across all trees	
	
	int.Vec <- unique(append(intervals[,"ageTop"], intervals[,"ageBase"]))
	
	if(runHistOn == "AllRates"){
	if(model == "BM" | model == "BMOU"){
	#BM model
	BMdat.List <- unlist(sapply(evoResults, function(x) sapply(x, function (y) y$BM$break.dates)))
	names(BMdat.List) <- NULL
	
	BM.hist <- hist(BMdat.List,breaks=sort(unique(unlist(intervals))), plot=FALSE)
	if(is.null(ylim)) {
		if(plotPercent == FALSE) ylimits <- c(0, max(BM.hist$counts)); ylabel <- "Number of Replicates with Rate Shift" #must imbed this after BM.hsitcounts is generated
		if(plotPercent == TRUE) ylimits <- c(0,1); ylabel <- "Proportion of Replicates with a Rate Shift (%)"
	} 
	if(!is.null(ylim)){
		ylimits <- ylim
		if(plotPercent == FALSE) ylabel <- "Number of Replicates with Rate Shift" #must imbed this after BM.hsitcounts is generated
		if(plotPercent == TRUE)  ylabel <- "Proportion of Replicates with a Rate Shift (%)"
	} 
	plot(BM.hist, col=NA, border=NA, labels=FALSE, freq = !plotPercent, cex=0.3, xaxp =c(55,5,5), xlim=rev(range(intervals)), 
			 ylim=ylimits, xaxt = "n", ylab = ylabel, xlab="Time (Ma)", main=NULL)#, main="BM Model: Number of Replicates with a Rate Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(BM.hist, col="orchid4", border="orchid1", labels=TRUE, freq=!plotPercent, cex=0.3, #xaxp =c(55,5,5), 
			 xlim=c(55, 0), ylim=ylimits, add=TRUE)
	}
	
	if(model == "OU" | model == "BMOU"){
	#OU model
	OUdat.List <- unlist(sapply(evoResults, function(x) sapply(x, function (y) y$OU$break.dates)))
	names(OUdat.List) <- NULL
	
	print(OUdat.List)
	
	OU.hist <- hist(OUdat.List,breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(BM.hist, col=NA, border=NA, labels=FALSE, freq=!plotPercent, cex=0.3, xaxp =c(55,5,5), xlim=rev(range(intervals)), 
			 ylim=ylimits, main="OU Model: Number of Replicates with a Rate Shift", xlab="Time (Ma)", 
			 fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #xaxp =c(55,5,10),
	overlayCzTimescale(do.subepochs=TRUE)
	axis(side = 1, at=rev(seq(4,56,2)),tcl=-0.5, labels = FALSE, 
			 fg="gray75", bg="gray75", col.axis="gray75")
	axis(side = 1, at= rev(seq(4,56,10)),tcl=-1, labels = TRUE,
			 fg="gray75", bg="gray75", col.axis="gray75")

	plot(OU.hist, col="firebrick4", border="firebrick1", labels=TRUE, freq=!plotPercent, cex=0.3, xaxp =c(55,5,5), 
			 xlim=c(55, 0), ylim=ylimits, add=TRUE)
	}}
	
	if(runHistOn == "BestRates"){
		if(model == "BM" | model == "BMOU"){
			#BM model
			BMdat.List <- unlist(sapply(evoResults, function (y) y$BM$break.dates))
			names(BMdat.List) <- NULL
			
			BM.hist <- hist(BMdat.List,breaks=sort(unique(unlist(intervals))), plot=FALSE)
			if(is.null(ylim)) {
				if(plotPercent == FALSE) ylimits <- c(0, max(BM.hist$counts)); ylabel <- "Number of Replicates with Rate Shift" #must imbed this after BM.hsitcounts is generated
				if(plotPercent == TRUE) ylimits <- c(0,1); ylabel <- "Proportion of Replicates with a Rate Shift (%)"
			} 
			if(!is.null(ylim)){
				ylimits <- ylim
				if(plotPercent == FALSE) ylabel <- "Number of Replicates with Rate Shift" #must imbed this after BM.hsitcounts is generated
				if(plotPercent == TRUE)  ylabel <- "Proportion of Replicates with a Rate Shift (%)"
			} 
			if(relativeFreq == TRUE) BM.hist$density <- BM.hist$counts/sum(BM.hist$counts)
			
			plot(BM.hist, col=NA, border=NA, labels=FALSE, freq=!plotPercent, cex=0.3, xlim=c(56,4), #xaxp =c(55,5,5), #rev(range(intervals)), 
					 ylim=ylimits, xaxt = "n", ylab = ylabel, xlab="Time (Ma)", main=NULL,
					 fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #xaxp =c(55,5,10),
			overlayCzTimescale(do.subepochs=TRUE)
			axis(side = 1, at=rev(seq(4,56,2)),tcl=-0.5, labels = FALSE, 
					 fg="gray75", bg="gray75", col.axis="gray75")
			axis(side = 1, at= rev(seq(4,56,10)),tcl=-1, labels = TRUE,
					 fg="gray75", bg="gray75", col.axis="gray75")
			plot(BM.hist, col="orchid4", border="orchid1", labels=FALSE, freq=!plotPercent, cex=0.3, xaxp =c(55,5,5), 
					 xlim=c(55, 0), ylim=ylimits, add=TRUE)
		}
		
		if(model == "OU" | model == "BMOU"){
			#OU model
			OUdat.List <- unlist(sapply(evoResults, function (y) y$OU$break.dates))
			names(OUdat.List) <- NULL
			
			OU.hist <- hist(OUdat.List,breaks=sort(unique(unlist(intervals))), plot=FALSE)
			plot(BM.hist, col=NA, border=NA, labels=FALSE, freq=!plotPercent, cex=0.3, xaxp =c(55,5,5), xlim=rev(range(intervals)), 
					 ylim=ylimits, main="OU Model: Number of Replicates with a Rate Shift", xlab="Time (Ma)")
			overlayCzTimescale(do.subepochs=TRUE)
			plot(OU.hist, col="firebrick4", border="firebrick1", labels=TRUE, freq=!plotPercent, cex=0.3, xaxp =c(55,5,5), 
					 xlim=c(55, 0), ylim=ylimits, add=TRUE)
		}
	}
	return()
}

evoModelRates <- function(evoResults, runOnRates = "AllRates"){
	
	if(runOnRates == "AllRates"){
		#BM model
		BMdat.List <- unlist(sapply(evoResults, function(x) sapply(x, function (y) length(y$BM$sigma))))
		names(BMdat.List) <- NULL
		cat("Number of Breaks in BM Models", "\n")
		cat(names(table(BMdat.List)),"\n")
		cat(table(BMdat.List),"\n")
	
		cat("\n")
	
		#OU model
		OUdat.List <- unlist(sapply(evoResults, function(x) sapply(x, function (y) length(y$OU$theta))))
		names(OUdat.List) <- NULL
		cat("Number of Breaks in OU Models", "\n")
		cat(names(table(OUdat.List)),"\n")
		cat(table(OUdat.List),"\n")
	}

	if(runOnRates == "BestRates"){
		#BM model
		BMdat.List <- unlist(sapply(evoResults, function (y) length(y$BM$break.dates)))
		names(BMdat.List) <- NULL
		cat("Number of Breaks in BM Models", "\n")
		cat(names(table(BMdat.List)),"\n")
		cat(table(BMdat.List),"\n")
		hist(BMdat.List, breaks = c(seq(0,max(BMdat.List)+1,1)), right = FALSE) # how to handle 1 rate models? #seems to break when there are no breaks wthin given column
		cat("median # breaks", median(BMdat.List))
		cat("\n")
		
		#OU model
		OUdat.List <- unlist(sapply(evoResults, function (y) length(y$OU$breaks)))
		names(OUdat.List) <- NULL
		cat("Number of Breaks in OU Models", "\n")
		cat(names(table(OUdat.List)),"\n")
		cat(table(OUdat.List),"\n")
	}
	
}

plotRatesBM <- function(this.rez, intervals, PlotXmax = NULL, PlotXmin = NULL, num.Trees = 1) {
	quartz(width=12, height = 6)
	par(mfrow=c(1,1),  mar=c(4,3,1.5,2.5), mgp=c(2, 1,0), oma=c(0,0,0,4),xaxs="i")
	
	if(num.Trees == 1){
		#this.break.dates <- c(max(this.tree$node.date), this.rez$BM$break.dates, min(this.tree$node.date))
		this.break.dates <- c(56, this.rez$BM$break.dates, 2)
		
		breakMat <- matrix(nrow= length(this.break.dates)-1, ncol=3)
		for(jj in seq_len(length(this.break.dates)-1)){
			breakMat[jj,] <- c(this.break.dates[[jj]], this.break.dates[[jj+1]], this.rez$BM$sigma[jj])
		}
		colnames(breakMat) <- c("FO","LO","Sigma")
		
		if(is.null(PlotXmax)) {
			xlimMax <- max(breakMat[,"FO"])
		} else {xlimMax <- PlotXmax}
		if(is.null(PlotXmin)) {
			xlimMin <- min(breakMat[,"LO"])
		} else {xlimMin <- PlotXmin}
	
		plot(breakMat[,"FO"], breakMat[,"Sigma"], xlim=c(xlimMax, xlimMin), 
				 ylim = c(0,max(breakMat[,"Sigma"])+0.01), type="n", ylab="Sigma", 
				 xlab="Time (Ma)", cex.axis=0.5, cex.lab=0.5, xaxt = "n", col=alphaColor("black",0.7),
				 fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #xaxp =c(55,5,10),
		overlayCzTimescale(do.subepochs=TRUE)
		axis(side = 1, at=rev(seq(4,56,2)),tcl=-0.5, labels = FALSE, 
				 fg="gray75", bg="gray75", col.axis="gray75")
		axis(side = 1, at= rev(seq(4,56,10)),tcl=-1, labels = TRUE,
				 fg="gray75", bg="gray75", col.axis="gray75")
		for (i in seq_len(nrow(breakMat))) {
			if (is.finite(breakMat[,"FO"][i]) & is.finite(breakMat[,"LO"][i]) & breakMat[,"FO"][i] != breakMat[,"LO"][i]) lines(x=breakMat[i,c("FO","LO")], y=c(breakMat[,"Sigma"][i], breakMat[,"Sigma"][i]), lwd=2, pch=21, col=alphaColor("gray0", 0.5)) #alphaColor(orderColors[i], 0.5)
		}
	} else {
		for(kk in seq_len(length(this.rez))){
			#this.break.dates <- c(max(this.tree[[kk]]$node.date), this.rez[[kk]]$BM$break.dates, min(this.tree[[kk]]$node.date))
			#this.break.dates <- c(max(intervals[[kk]]), this.rez[[kk]]$BM$break.dates, min(intervals[[kk]]))
			this.break.dates <- c(56, this.rez[[kk]]$BM$break.dates, 2)
			
			breakMat <- matrix(nrow= length(this.break.dates)-1, ncol=3)
			for(jj in seq_len(length(this.break.dates)-1)){
				breakMat[jj,] <- c(this.break.dates[[jj]], this.break.dates[[jj+1]], this.rez[[kk]]$BM$sigma[jj])
			}
			colnames(breakMat) <- c("FO","LO","Sigma")
			
			if(kk == 1){
			par(mar=c(3,4, 2.5,0.5))
			if(is.null(PlotXmax)) {
				xlimMax <- max(breakMat[,"FO"])
			} else {xlimMax <- PlotXmax}
			if(is.null(PlotXmin)) {
				xlimMin <- min(breakMat[,"LO"])
			} else {xlimMin <- PlotXmin}
			
			plot(breakMat[,"FO"], breakMat[,"Sigma"], xlim=c(xlimMax, xlimMin), 
					 #ylim = c(0,max(unlist(sapply(this.rez, function (y) y$BM$sigma)))+0.01), 
					 ylim = c(0,0.15), 
					 type="n", ylab="Sigma",  xlab="Time (Ma)", xaxt = "n", cex.axis=1, cex.lab=1, 
					 col=alphaColor("black",0.7), fg="gray75", bg="gray75", 
					 col.axis="gray75", col.lab="gray75") #xaxp =c(55,5,10),
			overlayCzTimescale(do.subepochs=TRUE)
			axis(side = 1, at=rev(seq(4,56,2)),tcl=-0.5, labels = FALSE, 
					 fg="gray75", bg="gray75", col.axis="gray75")
			axis(side = 1, at= rev(seq(4,56,10)),tcl=-1, labels = TRUE,
						fg="gray75", bg="gray75", col.axis="gray75")
			}
			for (i in seq_len(nrow(breakMat))) {
				if (is.finite(breakMat[,"FO"][i]) & is.finite(breakMat[,"LO"][i]) & breakMat[,"FO"][i] != breakMat[,"LO"][i]) lines(x=breakMat[i,c("FO","LO")], y=c(breakMat[,"Sigma"][i], breakMat[,"Sigma"][i]), lwd=2, pch=21, col=alphaColor("gray0", 0.25)) #alphaColor(orderColors[i], 0.5)
			}
		}
	}
	
	Int.biggest <- intervals
	#get max size of intervals across all trees
	#for(kk in seq(1, length(intervals),1)){
	#	if(length(intervals[[kk]]$ageTop) > length(Int.biggest$ageTop)) Int.biggest <- intervals[[kk]]
	#}
	
	median.sigmaMat <- matrix(nrow = length(this.rez), ncol = length(Int.biggest$ageTop))
	#add in line for median sigma through time
	#maybe incorporate shoulder plot
	for(kk in seq(1, length(this.rez), 1)){
		median.sigmaOneRep <- getSigmaVecForOneRep(this.rez = this.rez[[kk]],intervals = intervals)
		names(median.sigmaOneRep) <- rowMeans(intervals)
		
		if(length(median.sigmaOneRep) < length(Int.biggest$ageTop)) median.sigmaOneRep <- append(median.sigmaOneRep, rep(NA, length(Int.biggest$ageTop)-length(median.sigmaOneRep)))
	
		median.sigmaMat[kk,] <- median.sigmaOneRep
	}
	#add each entry to matrix if value then value else 0
	#add # of 0 to make interval equal to full size or should I use NULL
	median.sigma <- apply(median.sigmaMat, 2, median)
	#median.sigma <- apply(median.sigmaMat, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975),na.rm=TRUE)
	lines(rev(rowMeans(Int.biggest)), median.sigma, lwd=6, col="orange")
	lines(rev(rowMeans(Int.biggest)), median.sigma, lwd=4, col="red")
	#polygon(c( rev(rowMeans(Int.biggest)), rowMeans(Int.biggest)), c(median.sigma[2,], rev(median.sigma[4,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
	#polygon(c( rev(rowMeans(Int.biggest)), rowMeans(Int.biggest)), c(median.sigma[1,], rev(median.sigma[5,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
	
	
	#dev.off()
	return()
}

getSigmaVecForOneRep <- function(this.rez, intervals) {
	this <- this.rez[["BM"]]
	this.dates <- rev(intervals[,1]) #how does this work? intervals not an argument
	date.vec <- sort(c(max(intervals), this$break.dates, 0), decreasing=TRUE)
	this.sigma <- array(dim=length(this.dates))
	for (i in seq(from=1, to=length(date.vec)-1)) {
		this.sigma[this.dates > date.vec[i+1] & this.dates <= date.vec[i]] <- this$sigma[i]
	}
	this.sigma
}

plotSigmasForOneRep <- function(this.rez, intervals) {
	this <- this.rez[[length(this.rez)]]["BM"]
	# this.dates <- sort(c(max(intervals), this$break.dates + 0.5, this$break.dates - 0.5, 0), decreasing=TRUE)
	# lines(this.dates, 	unlist(lapply(this$sigma, rep, times=2)), col=adjustcolor("black", alpha.f=0.25))

	this.dates <- sort(c(max(intervals), this$break.dates, 0), decreasing=TRUE)
	for (this.interval in seq_len(length(this.dates)-1)) {
		lines(c(this.dates[this.interval], this.dates[this.interval+1]), y=rep(this$sigma[this.interval], 2), col=adjustcolor("black", alpha.f=0.25))
	}	
}

plotAllRezListSigmas <- function(rez.list, intervals) {
	source("https://dl.dropbox.com/s/dozeb8o2pxu4sbj/CzTimescale.R")
	plot(x=0, xlim=c(55, 0), ylim=c(0,0.15))
	overlayCzTimescale(do.subepochs=TRUE)
	sapply(rez.list, plotSigmasForOneRep, intervals=intervals)
	median.sigma <- apply(sapply(rez.list, getSigmaVecForOneRep, intervals=intervals), 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
	lines(rev(rowMeans(intervals)), median.sigma[3,], lwd=2, col="firebrick4")
}

getSigmaRateRiseDecline <- function(rez.list){
	
#	unlist(rez.list[[1]]$BM$sigma)
#	rez.list[[1]]$BM$break.dates
	
	increaseRateBreak <- vector()
	declineRateBreak <- vector()
	
	sigmaMat <- matrix(nrow = length(rez.list$BM$break.dates),ncol = 2)
	#sigmaMat[1,] <- cbind(SigmaBase=rez.list[[1]]$BM$sigma[1], SigmaTop=rez.list[[1]]$BM$sigma[2])
	
	for(ii in seq(1, length(rez.list$BM$break.dates),1)) {
	sigmaMat[ii,] <- cbind(rez.list$BM$sigma[ii], rez.list$BM$sigma[ii+1])
	}
	
	sigmaMat <- cbind(sigmaMat,rez.list$BM$break.dates)
	rownames(sigmaMat) <- rez.list$BM$break.dates
	colnames(sigmaMat) <- c("SigmaBase", "SigmaTop","BreakDate")
	
	for(jj in seq(1, length(rownames(sigmaMat)),1))	{
		if(sigmaMat[jj,"SigmaBase"] < sigmaMat[jj,"SigmaTop"]) increaseRateBreak <- append(increaseRateBreak, sigmaMat[jj,3])
		if(sigmaMat[jj,"SigmaBase"] > sigmaMat[jj,"SigmaTop"]) declineRateBreak <- append(declineRateBreak, sigmaMat[jj,3])
	}
	
	if(!is.numeric(increaseRateBreak)) increaseRateBreak <- NULL
	if(!is.numeric(declineRateBreak)) declineRateBreak <- NULL
	
	breakList <- list(increaseBreak = increaseRateBreak, declineBreak = declineRateBreak)
	
	return(breakList)
}

plotBreaksRiseDecline <- function(RiseDeclineList, intervals, plotPercent = FALSE){
	
	#plot increase break frequencies
	BMdat.ListPos <- unlist(sapply(RiseDeclineList, function (y) y$increaseBreak))
	#names(BMdat.ListPos) <- NULL
	BM.histPos <- hist(BMdat.ListPos,breaks=sort(unique(unlist(intervals))), plot=FALSE)
	#plot declining break frequencies
	BMdat.ListNeg <- unlist(sapply(RiseDeclineList, function (y) y$declineBreak))
	#names(BMdat.ListPos) <- NULL
	BM.histNeg <- hist(BMdat.ListNeg,breaks=sort(unique(unlist(intervals))), plot=FALSE)

	totalCountPos <- sum(BM.histPos$counts)
	totalCountNeg <- sum(BM.histNeg$counts)
	
	BM.histNeg$counts = - BM.histNeg$counts
	hmax = max(BM.histPos$counts) +30
	hmin = min(BM.histNeg$counts) -30
	X = c(BM.histPos$breaks, BM.histNeg$breaks)
	xmax = 56#max(X)
	xmin = 4 #min(X)
	
	quartz(width=12, height = 6)
	#par(mfrow=c(1,1), mgp=c(2, 1,0), mar=c(4,4,0.5,4))
	par(mfrow=c(1,1), mgp=c(2, 1,0), mar=c(4,3,1.5,4), oma=c(0,0,0,0), xaxs="i")
	
	if(!plotPercent == TRUE){
		plot(1, ylim=c(hmin, hmax), col="orchid4", xlim=c(xmax, xmin), xlab = "Time (Ma)", 
				 ylab = "Frequency of Replicates with a Rate Shift", yaxt = "n",xaxt= "n",
				 fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75")
		axis(side = 1, at=rev(seq(4,56,2)),tcl=-0.5, labels = FALSE, 
				 fg="gray75", bg="gray75", col.axis="gray75")
		axis(side = 1, at= rev(seq(4,56,10)),tcl=-1, labels = TRUE,
				 fg="gray75", bg="gray75", col.axis="gray75")
		axis(side = 2, at= c(seq(0,hmax,100)), labels = TRUE, 
				 fg="gray75", bg="gray75", col.axis="gray75")
		axis(side = 2, at= c(-1*seq(0,hmax,100)), labels = abs(c(-1*seq(0,hmax,100))), 
				 fg="gray75", bg="gray75", col.axis="gray75")
		#plot(BM.histPos, ylim=c(hmin, hmax), col="orchid4", xlim=c(xmax, xmin), labels = TRUE)
		overlayCzTimescale(do.subepochs=TRUE)
		
		lines(BM.histPos, col="orchid4", labels = FALSE)
		lines(BM.histNeg, col="blue", labels = FALSE )#ifelse(BM.histNeg$counts < 0, BM.histNegLabels, BM.histNeg$counts))
		#vector for x/y pos for negative bars
		labelPos <- cbind(BM.histNeg$mids,BM.histNeg$counts-5)
	#	text(labelPos, labels = ifelse(-1*BM.histNeg$counts > 0, -1*BM.histNeg$counts, ""))
	}
	if(plotPercent == TRUE){
		plot(1, ylim=c(-0.10, 0.10), col="orchid4", xlim=c(xmax, xmin), xlab = "Time (Ma)", 
				 ylab = "Proportion of Replicates with a Rate Shift (%)", yaxt="n", xaxt= "n",
				 fg="gray75", bg="gray75", 
				 col.axis="gray75", col.lab="gray75")
		axis(side = 1, at=rev(seq(4,56,2)),tcl=-0.5, labels = FALSE, 
				 fg="gray75", bg="gray75", col.axis="gray75")
		axis(side = 1, at= rev(seq(4,56,10)),tcl=-1, labels = TRUE,
				 fg="gray75", bg="gray75", col.axis="gray75")
		axis(side = 2, at= c(seq(0,0.2,0.05)), labels = TRUE, 
				 fg="gray75", bg="gray75", col.axis="gray75")
		axis(side = 2, at= c(-1*seq(0,0.2,0.05)), labels = abs(c(-1*seq(0,0.2,0.05))), 
				 fg="gray75", bg="gray75", col.axis="gray75")
		#plot(BM.histPos, ylim=c(hmin, hmax), col="orchid4", xlim=c(xmax, xmin), labels = TRUE)
		overlayCzTimescale(do.subepochs=TRUE)
		
	BM.histPos$density <- BM.histPos$counts/(totalCountPos + totalCountNeg)
	BM.histNeg$density <- BM.histNeg$counts/(totalCountNeg + totalCountPos)
	#	BM.histNeg$density <- -1*BM.histNeg$density
		plot(BM.histPos,freq = FALSE, add = TRUE, col="orchid4")
		plot(BM.histNeg,freq = FALSE, add = TRUE, col="blue")
	}
#	dev.off()
	return()
}
