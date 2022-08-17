#function needs to search between intervals and pulls out changing taxa
#read in objects

setwd("~/Dropbox/ungulate_RA/RCode/Orignation_Extinction_Results")
#setwd("~/Dropbox/ungulate_RA/RCode/Orignation_Extinction_Results")
#install.packages("ggplot2")
library(ggplot2)

load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=TRUE##------ Mon Jan 14 17:48:35 2019 ------##.Rdata")

#check jon ks code
	#ksMatrix()
	#pairwiseKSTestsSubsequent

getUngulateOnly <- function(repIntSp, bigList, shortFam, thisMat)
{
	#Need to remove non-target taxonomic designations and/or species (i.e. family and genera names, species not in target clades)
	#testList <- lapply(repIntSp[[1]], function(y) sapply(y, function(x) x[x %in% bigList[as.character(bigList$family) %in% shortFam,]$accepted_name], simplify="array"))
	
	###Crop dataset to onnly be ungulates
	repCulled_Int <- list()
	repCulled_All <- list()
	ungOrder <- c("Artiodactyla", "Perissodactyla")
	#repIntSp[[1]][[1]][repIntSp[[1]][[1]] %in% bigList$accepted_name[as.character(bigList$family) %in% shortFam]]
	for(xx in seq_len(length(repIntSp))){
		for(yy in seq_len(length(repIntSp[[xx]]))){
			#intCulled <- repIntSp[[xx]][[yy]][repIntSp[[xx]][[yy]] %in% bigList$accepted_name[as.character(bigList$family) %in% shortFam & bigList$order == c("Artiodactyla", "Perissodactyla")]]
			intCulled <- repIntSp[[xx]][[yy]][repIntSp[[xx]][[yy]] %in% bigList$accepted_name[as.character(bigList$family) %in% shortFam]]
			
			#remove species lacking a body mass
			intCulled <- intCulled[intCulled %in% thisMat$species]
			#remove non-ungulate species
			intCulled <- intCulled[intCulled %in% bigList$accepted_name[as.character(bigList$order) %in% ungOrder]]
			
			repCulled_Int[[yy]] <- intCulled
		}
		names(repCulled_Int) <- names(repIntSp[[xx]])
		repCulled_All[[xx]] <- repCulled_Int
		#	print(xx)
	}
	return(repCulled_All)
}

#change this to handle a single instance of the repIntSp, 
#to get full plot make sure to run using an apply
findOrigExt_Breaks <- function(repIntSpSingle, optList, measureMat, bigList, shortFam, breaks = NULL, useBreaks = TRUE, startDate=56,endDate=3,plot.check = FALSE,plot.type = "lines")
{
		if(is.null(breaks)){
			print("Get breaks\n")
			#locate shifts (get optimal breaks)
			optBreaks <- 9999 
			for(ii in seq_len(length(optList))){
				#breaks <- optList[[length(optList)]]$optBreaks
				if(optBreaks > optList[[ii]]$AICc) {
					optBreaks <- optList[[ii]]$AICc
					#print(optList[[ii]]$AICc)
					breaks <- optList[[ii]]$optBreaks
				}
			}
		}
		names(repIntSpSingle) <- gsub(pattern = " Ma", replacement = "", x = names(repIntSpSingle))
		interval.names <- vector()
		taxList <- list()
		
		print("Bin taxa \n")
		#bin taxa within each shift
		for(nn in seq(1,length(breaks) + 1,1)){
			# print(nn)
			#locate intervals
			if(nn == 1){
				#startDate <- 56
				topDate <- startDate
				bottomDate <- breaks[1]
				#print("First break)
			} 
			if(nn == length(breaks) + 1){
				topDate <- breaks[nn-1]
				bottomDate <- endDate
				#print("Last break")
			} 
			if(nn > 1 & nn < length(breaks) + 1){
				topDate <- breaks[nn-1]
				bottomDate<- breaks[nn]
				#print(nn)  
			}
			
			print("Extract Species\n")
			#extract species
			tax.vec <- vector()
			
			for(xx in seq_len(length(repIntSpSingle))){
			#	for(kk in as.numeric(sub("* Ma", "", names(repIntSpSingle[[xx]])))){
					if(xx < topDate & xx > bottomDate){
						tax.vec <- append(tax.vec, repIntSpSingle[[xx]]) #[[kk]])
						#print(paste(nn,paste(kk, paste(topDate, bottomDate))))
						#print(tax.vec)
					}
					taxList[[nn]] <- unique(tax.vec)
					#run.count <- run.count + 1
			#	}
				interval.names <- append(interval.names, paste(paste(topDate,"_",sep=""),bottomDate,sep=""))
			}
			#print(taxList)
		}
		names(taxList) <- unique(interval.names)
		#need to find way to impliment for each interval/regime

		#optList <-optList_bm_median
		# "all" #c("all", "artio", "perisso")
		#"bm_median" #c("bm_median", , )
		#repIntSp
		
		plot.row <- ceiling((length(taxList)-1)/2)
		OrigExtList <- list()
		thisMat <- measureMat
		
		print("Generate Image File\n")
		#file.png <- paste(paste(file.name,paste("_",paste(timestamp(),".png",sep=""),sep=""),sep=""),sep="")
		#if(plot.check == TRUE){
		#	png(filename= file.png)}
		#par(mfrow=c(plot.row,2), mar=c(2,4,1,0.5))
		
		for(mm in seq(1,length(taxList)-1,1)){
		
			spec.list <- list()
			#find taxa that went extinct
			tax.extinct <- thisMat[taxList[[mm]][!taxList[[mm]] %in% taxList[[mm+1]]],]
			tax.extinct <- tax.extinct[complete.cases(tax.extinct[,"bodyMass"]),]
			spec.list[[1]] <- tax.extinct[,c("bodyMass")]
			
			#find taxa that originate
			tax.originate <- thisMat[taxList[[mm+1]][!taxList[[mm]] %in% taxList[[mm+1]]],]
			tax.originate <- tax.originate[complete.cases(tax.originate[,"bodyMass"]),]
			spec.list[[2]] <-tax.originate[,c("bodyMass")]
			
			names(spec.list) <- c("extinction", "origination")
			
			OrigExtList[[mm]] <- spec.list
		}
		#dev.off()
		
		names(OrigExtList) <- breaks

return(OrigExtList)
}

plotOrigExt <- function(OrigExtList, plot.type = "lines")# plot.OrigExt = c("Orig", "Ext", "Both"))
{
	for(yy in names(OrigExtList[[1]]))
	{
		for(xx in seq(1, length(OrigExtList),1))
		{
			if(xx == 1) 
			{
				#get initial plot window for extinction
				if(plot.type == "points"){
					plot(OrigExtList[[xx]][[yy]]$extinction, col="blue", pch=19, 
							 ylab = "Body Mass", main = breaks[mm])
				}
				if(plot.type == "lines"){
					plot(density(OrigExtList[[xx]][[yy]]$extinction), col="blue", 
							 xlab = "Body Mass", main = yy, ylim=c(0,1))
				}
			}
			if(xx > 1)
			{
				if(plot.type == "points"){
					points(OrigExtList[[xx]][[yy]]$extinction, col="blue", pch=19)}
				if(plot.type == "lines"){
					lines(density(OrigExtList[[xx]][[yy]]$extinction), col="blue")}
			}
		}
	
		for(xx in seq(1, length(OrigExtList),1))
		{
			if(plot.type == "points"){
				points(OrigExtList[[xx]][[yy]]$origination, col="red", pch=19)}
			if(plot.type == "lines"){
				lines(density(OrigExtList[[xx]][[yy]]$origination), col="red")}
		}
	}
}

findOrigExt_Interval <- function(repIntSpSingle, measureMat, bigList, shortFam, intervals = NULL, startDate=56,endDate=3)
{
	
	names(repIntSpSingle) <- gsub(pattern = " Ma", replacement = "", x = names(repIntSpSingle))
	interval.names <- vector()
	taxList <- list()
	
	print("Bin taxa \n")
	#bin taxa within each shift
	for(nn in rev(seq(1,nrow(intervals),1))){
		# print(nn)
		#locate intervals
		if(nn == nrow(intervals)){
			#startDate <- 56
			topDate <- startDate
			bottomDate <- intervals[nn,1]
			#print("First break)
		} 
		if(nn == 1){
			topDate <- intervals[nn+1,1]
			bottomDate <- 1
			#print("Last break")
		} 
		if(nn > endDate & nn < nrow(intervals)){
			topDate <- intervals[nn+1,1]
			bottomDate <- intervals[nn,1] # which one will be the -1/+1
			#print(nn)  
		}
		
		print("Extract Species\n")
		#extract species
		tax.vec <- vector()
		
		for(xx in seq_len(length(repIntSpSingle))){
			#	for(kk in as.numeric(sub("* Ma", "", names(repIntSpSingle[[xx]])))){
			if(xx <= topDate & xx >= bottomDate){
				tax.vec <- append(tax.vec, repIntSpSingle[[xx]]) #[[kk]])
				#print(paste(nn,paste(kk, paste(topDate, bottomDate))))
				#print(tax.vec)
			}
			taxList[[nn]] <- unique(tax.vec)
			#run.count <- run.count + 1
			#	}
			interval.names <- append(interval.names, paste(paste(topDate,"_",sep=""),bottomDate,sep=""))
		}
		#print(taxList)
	}
	names(taxList) <- rev(unique(interval.names))
	#need to find way to impliment for each interval/regime
	
	#optList <-optList_bm_median
	# "all" #c("all", "artio", "perisso")
	#"bm_median" #c("bm_median", , )
	#repIntSp
	
#	plot.row <- ceiling((length(taxList)-1)/2)
	OrigExtList <- list()

	print("Generate Image File\n")
	#file.png <- paste(paste(file.name,paste("_",paste(timestamp(),".png",sep=""),sep=""),sep=""),sep="")
	#if(plot.check == TRUE){
	#	png(filename= file.png)}
	#par(mfrow=c(plot.row,2), mar=c(2,4,1,0.5))
	
	for(mm in seq(endDate,length(taxList),1)){
		
		#need to think on how to fill the output file
		## origination should fill all but first bin whereas extinction fills all but last bin
		##to get at the fact that they are offset
		
		spec.list <- list()
		#find taxa that went extinct
		#tax.extinct <- thisMat[taxList[[mm]][!taxList[[mm-1]] %in% taxList[[mm]]],]
		#tax.extinct <- tax.extinct[complete.cases(tax.extinct[,"bodyMass"]),]
		#spec.list[[1]] <- tax.extinct[,c("bodyMass")]
		
			#find taxa that originate
			tax.originate <- thisMat[taxList[[mm-1]][!taxList[[mm]] %in% taxList[[mm-1]]],]
			tax.originate <- tax.originate[complete.cases(tax.originate[,"bodyMass"]), ]#"bodyMass"]
			#spec.list <-tax.originate[,c("bodyMass")]
			
			#names(spec.list) <- c("origination")
			
			OrigExtList[[mm]] <- tax.originate
		
	}
	#dev.off()
	
	names(OrigExtList) <- rev(intervals[nrow(intervals):1,2]-.5)
	
	return(OrigExtList)
}

###########################################################################
###########################################################################
###########################################################################

getBm4List <- function(repIntSp, thisMat)
{
	repBM_Int <- list()
	repBM_All <- list()
	if(length(repIntSp) <= 60){
		for(xx in seq_len(length(repIntSp))){
				intBM <- thisMat[thisMat$species %in% repIntSp[[xx]],"bodyMass"]
				#thisMat[thisMat %in% repIntSp[[mm]],"bodyMass"]
			#	if(intBM == numeric(0))
				repBM_Int[[xx]] <- intBM
				repBM_All <- repBM_Int
		}
	}
	if(length(repIntSp) > 60){
		for(xx in seq_len(length(repIntSp))){
			for(yy in seq_len(length(repIntSp[[xx]]))){
				
				intBM <- thisMat[thisMat$species %in% repIntSp[[xx]][[yy]],"bodyMass"]
				#thisMat[thisMat %in% repIntSp[[mm]],"bodyMass"]
				repBM_Int[[yy]] <- intBM
			}
			names(repBM_Int) <- names(repIntSp[[xx]])
			repBM_All[[xx]] <- repBM_Int
			#	print(xx)
		}
	}
	return(repBM_All)
}

getKS4IntervalsOrig <- function(repOrigBmUng, repBmUng,result.out = c("SigInt", "AllInt"), out.putType = c("Intervals","Breaks"))
{
	sig.IntervalsKS <- matrix(nrow = 2, ncol = length(repOrigBmUng))
	all.IntervalsKS <- matrix(nrow = 2, ncol = length(repOrigBmUng))
	if(out.putType == "Intervals"){
		for(pp in seq(1, length(repOrigBmUng),1))
		{
			print(pp)
			ks.result <- ks.test(repBmUng[[pp+1]],repOrigBmUng[[pp]])
			inputs <- as.matrix(c(as.numeric(ks.result$p.value), as.numeric(ks.result$statistic)))
			colnames(inputs) <- paste(pp+1, paste("to ", pp, sep=""),sep=" ")
			
			all.IntervalsKS[,pp] <-inputs
			
			if(ks.result$p.value < 0.05) sig.IntervalsKS[,pp] <- inputs
		}
	}
	if(out.putType == "Breaks"){
		for(pp in seq(1, length(repOrigBmUng),1))
		{
			ks.result <- ks.test(repBmUng[[pp+1]],repOrigBmUng[[pp]])
			inputs <- as.matrix(c(as.numeric(ks.result$p.value), as.numeric(ks.result$statistic)))
			colnames(inputs) <- paste(pp+1, paste("to ", pp, sep=""),sep=" ")
			
			all.IntervalsKS[,pp] <-inputs
			
			if(ks.result$p.value < 0.05) sig.IntervalsKS[,pp] <- inputs
		}
	}
	rownames(sig.IntervalsKS) <- rownames(all.IntervalsKS) <- c("P-value", "Statstic")
	
	if( result.out == "SigInt")	return(sig.IntervalsKS)
	if( result.out == "AllInt")	return(all.IntervalsKS)
}

getKS4IntervalsExt <- function(repExtBmUng, repBmUng,result.out = c("SigInt", "AllInt"), out.putType = c("Intervals","Breaks"))
{
	sig.IntervalsKS <- matrix(nrow = 2, ncol = length(repExtBmUng))
	all.IntervalsKS <- matrix(nrow = 2, ncol = length(repExtBmUng))
	if(out.putType == "Intervals"){
		for(pp in seq(1, length(repExtBmUng),1))
		{
			ks.result <- ks.test(repBmUng[[pp]],repExtBmUng[[pp]])
			inputs <- as.matrix(c(as.numeric(ks.result$p.value), as.numeric(ks.result$statistic)))
			colnames(inputs) <- paste(pp+1, paste("to ", pp, sep=""),sep=" ")
			
			all.IntervalsKS[,pp] <-inputs
			
			if(ks.result$p.value < 0.05) sig.IntervalsKS[,pp] <- inputs
		}
	}
	if(out.putType == "Breaks"){
		for(pp in seq(1, length(repExtBmUng),1))
		{
			ks.result <- ks.test(repExtBmUng[[pp]], repBmUng[[pp]])
			inputs <- as.matrix(c(as.numeric(ks.result$p.value), as.numeric(ks.result$statistic)))
			colnames(inputs) <- paste(pp+1, paste("to ", pp, sep=""),sep=" ")
			
			all.IntervalsKS[,pp] <-inputs
			
			if(ks.result$p.value < 0.05) sig.IntervalsKS[,pp] <- inputs
		}
	}
	rownames(sig.IntervalsKS) <- rownames(all.IntervalsKS) <- c("P-value", "Statstic")
	
	if( result.out == "SigInt")	return(sig.IntervalsKS)
	if( result.out == "AllInt")	return(all.IntervalsKS)
}

getBreakTaxList <- function(repIntSpSingle, breaks, startDate = 56, endDate = 3)
{
	taxList <- list()
	
	print("Bin taxa \n")
	#bin taxa within each shift
	for(nn in seq(1,length(breaks) + 1,1)){
		# print(nn)
		#locate intervals
		if(nn == 1){
			#startDate <- 56
			topDate <- startDate
			bottomDate <- breaks[1]
			#print("First break)
		} 
		if(nn == length(breaks) + 1){
			topDate <- breaks[nn-1]
			bottomDate <- endDate
			#print("Last break")
		} 
		if(nn > 1 & nn < length(breaks) + 1){
			topDate <- breaks[nn-1]
			bottomDate<- breaks[nn]
			#print(nn)  
		}
		
		print("Extract Species\n")
		#extract species
		tax.vec <- vector()
		
		for(xx in seq_len(length(repIntSpSingle))){
			#	for(kk in as.numeric(sub("* Ma", "", names(repIntSpSingle[[xx]])))){
			if(xx < topDate & xx > bottomDate){
				tax.vec <- append(tax.vec, repIntSpSingle[[xx]]) #[[kk]])
				#print(paste(nn,paste(kk, paste(topDate, bottomDate))))
				#print(tax.vec)
			}
			taxList[[nn]] <- unique(tax.vec)
			#run.count <- run.count + 1
			#	}
			interval.names <- append(interval.names, paste(paste(topDate,"_",sep=""),bottomDate,sep=""))
		}
		#print(taxList)
	}
	#names(taxList) <- unique(interval.names)
	return(taxList)
}

check4num0 <- function(repIntSp, repNum, out.type = c("RepsList","num0List"))
{
	count <- 1
	num0 <- matrix(ncol = 2)
	int.numeric.is0 <- list()
	for(gg in seq(1,length(repIntSp),1))
	{
		if(length(repIntSp[[gg]]) == 0)
		{
			num0 <- c(repNum,gg)
			if(is.numeric(num0)) 
			{
				int.numeric.is0[[count]] <- num0
				count <- count + 1
				rm(num0)
			}
			repIntSp[[gg]] <- 0
		}
	}	

	if(out.type =="RepList") return(repIntSp)
	if(out.type =="num0List") return(int.numeric.is0)
}

getListNum0Ints <- function(numList)
{
	numVec <- unlist(numList)
	if(is.null(numVec)) 
	{
		numMat <- matrix(c(0,0), nrow=1,ncol = 2)
		return(numMat)
	}
	count <- 1
	seqVec<- seq(1, length(numVec)-1,2)
	numMat <- matrix(nrow=length(seqVec),ncol = 2)
	for(yy in seq(1, length(numVec)-1,2))
	{
		numMat[count,] <- c(numVec[yy], numVec[yy+1])
		count <- count +1
	}
	return(numMat)
}

checkRangeThrough1MaBins <- function(repIntTest, ints)
{
	#function to find if rangethrough is working
	
	species <- list()
	for(xx in ints) species[[xx-(min(ints)-1)]] <- repIntTest[[xx]]
	uni.sp <- unique(unlist(species))
	
	checkMat <- matrix(nrow = length(uni.sp), ncol = length(ints))
	rownames(checkMat) <- uni.sp
	colnames(checkMat) <- ints+0.5
	
	# find way to parse through each species and mark if they are in a given interval
	for(xx in seq(1, nrow(checkMat),1))
	{
		for(yy in seq(1, ncol(checkMat),1))
		{
			#for(zz in repIntTest[[yy+49]])
			##{
			test.true <- which(rownames(checkMat)[xx] == repIntTest[[yy+(min(ints)-1)]])
			if(length(test.true) == 0)
			{
				checkMat[xx,yy] <- 0
			}
			if(length(test.true) >= 1) checkMat[xx,yy] <- 1
			#}
		}
	}
	return(checkMat)
}

rangeThrough_alt <- function(IntMat)
{
	#IntMat <- repIntCheck
	#IntMat.rangeThrough <- IntMat
	for(xx in seq_len(nrow(IntMat)))
	{
		int.index <- which(IntMat[xx,] == 1)
		int.new <- min(int.index):max(int.index)
		for (yy in int.new)
		{
			IntMat[xx,yy] <- 1
		}
	}
	
	repIntRangeThrough <- list()
	
	for(gg in seq_len(ncol(IntMat)))
	{
		specVec <- vector()
		for(hh in seq_len(nrow(IntMat)))
		{
			if(IntMat[hh,gg] == 1) specVec <- append(specVec, rownames(IntMat)[hh])
		}
		repIntRangeThrough[[gg]] <- specVec
	}
	names(repIntRangeThrough) <- colnames(IntMat)
	return(repIntRangeThrough)
}

checkRangeThrough <- function(repIntTest, ints)
{
	#function to find if rangethrough is working
	int.label<- vector()
	species <- list()
	for(xx in seq(1, length(ints)-1,1)) 
	{
		species[[xx]] <- repIntTest[[xx]] #species[[xx-(min(ints)-1)]] 
		int.label[xx] <- (ints[xx]+ints[xx+1])/2
	}
	uni.sp <- unique(unlist(species))
	
	checkMat <- matrix(nrow = length(uni.sp), ncol = length(ints)-1)
	rownames(checkMat) <- uni.sp
	colnames(checkMat) <- int.label
	
	# find way to parse through each species and mark if they are in a given interval
	for(xx in seq(1, nrow(checkMat),1))
	{
		for(yy in seq(1, ncol(checkMat),1))
		{
			#for(zz in repIntTest[[yy+49]])
			##{
			test.true <- which(rownames(checkMat)[xx] == repIntTest[[yy]]) #repIntTest[[yy+(min(ints)-1)]])
			if(length(test.true) == 0)
			{
				checkMat[xx,yy] <- 0
			}
			if(length(test.true) >= 1) checkMat[xx,yy] <- 1
			#}
		}
	}
	return(checkMat)
}
