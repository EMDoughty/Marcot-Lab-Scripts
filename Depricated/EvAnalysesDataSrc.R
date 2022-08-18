
############################################################################################################################################

#function is meant to be ran prior to analysis to bring all measurement datasets together into a single entity
getSingleSpeciesMatrix <- function() {
	#compile and lable dental measurments for specimens
	specimenMat <- getSpecimenMatFromMeasurements(filename="https://dl.dropbox.com/s/423x0zn3mpxuwc7/specimens.csv")
	specimenMat <- merge(specimenMat, getBirlenbachBlastoSpecimens(filename="https://dl.dropbox.com/s/943dq4lb1kd1h1e/blastoBirlenbach20140207.csv", info.file="https://dl.dropbox.com/s/nhrqzrclwtr5c8e/info2.csv"), all=TRUE)
	specimenMat <- merge(specimenMat, getLiteratureSpecimenMat(filename="https://dl.dropbox.com/s/qef8ts9a73j5ukb/literature.csv"), all=TRUE)
	
	specimenMat[sapply(specimenMat, is.nan)] <- NA
	specimenMat$species <- getCurrentTaxa(tax.vec = specimenMat$species)
	
	upLabels<-c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W") #"P2_L","P2_W",
	loLabels <- casefold(upLabels)
	
	### node that specimens are aggregated by their medians, so as to minimize the effect of outlier measurements
	thisMat <- aggregate(specimenMat[,c(upLabels, loLabels)], by=list(species=specimenMat$species), median, na.rm=TRUE)
	
	thisMat[,sapply(thisMat, is.numeric)] <- thisMat[,sapply(thisMat, is.numeric)] / 10  #converts mm measurements to cm for compatibility with Janis regressions
	thisMat <- transform(thisMat, p4_a=p4_l*p4_w, m1_a=m1_l*m1_w, m2_a=m2_l*m2_w, m3_a=m3_l*m3_w, M2_A=M2_L*M2_W)
	thisMat[sapply(thisMat, is.nan)] <- NA
	
	#thisMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = thisMat$species)
	rownames(thisMat) <- thisMat$species
	
	return (thisMat)
}

############################################################################################################################################

#Compile and format matrix of all measurements from multiple sources
appendRegressionCategories <- function(thisMat, regMat) {
	
	uniqTax <- lapply(c("Artiodactyla", "Perissodactyla"), FUN=getTaxonomyForOneBaseTaxon)
	uniqTax <- rbind(uniqTax[[1]], uniqTax[[2]])
	uniqTax$taxon_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$taxon_name)
	
	#delete unique taxa that lack family and genus
	thisMat$family <- uniqTax$family[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
	# thisMat$family <- sapply(thisMat$family, as.character)
	thisMat$family[thisMat$family == ""] <- NA
	
	thisMat$genus <- uniqTax$genus[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
	# thisMat$genus <- sapply(thisMat$genus, as.character)
	thisMat$genus[thisMat$genus == ""] <- NA
	
	#### 
	#append regression catagories to each species
	####
	
	family.names <- uniqTax$family[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
	reg.vec <- regMat$reg[!is.na(regMat$family)][match(family.names, regMat$family[!is.na(regMat$family)])]  # this is the regression "labels" for the species from measure.mat in the correct order, based on family name
	
	genus.names <- uniqTax$genus[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
	reg.vec[is.na(reg.vec)] <- regMat$reg[match(genus.names, regMat$genus)][is.na(reg.vec)]   #this is the regression "labels" for the species from measure.mat in the correct order, based on genus name; it appears that having regMat$genus[!is.na(regMat$genus)] will cause the index to improperly assign regressions
	
	thisMat$reg.vec <- reg.vec

	#check for taxa that are not recieving a regression
	#missingReg <- measure.mat[is.na(measure.mat $reg.vec),]
	#missingReg <- missingReg[!grepl("sp.",missingReg$species),]
	#missingReg <- missingReg[!grepl("indet.",missingReg$species),]
	#missingReg <- missingReg[!grepl("cf.",missingReg$species),]
	#write.csv(missingReg, '/Users/evandoughty/Dropbox/ungulate_RA/RCode/JonCode/2017_2_27_missingReg.csv')
	
	# nrow(thisMat) #803 species remain after final removal 
	
	# rownames(thisMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisMat))
	
	return(thisMat)
}

############################################################################################################################################

approxBodyMass <- function(thisMat) {
	#Approximate body mass
	thisMat[,"bodyMass"] <- getBodyMassVectorFromThisMatAllMeasures(thisMat, linked.files=TRUE)
	thisMat$bodyMass <- fillMissingBodyMasses(thisMat)	# this fills taxa missing their body mass with the average body mass of its cogeners
	thisMat[!sapply(thisMat, is.finite)] <- NA
	#rownames(thisMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisMat))
	thisMat$species <- rownames(thisMat)
	return(thisMat)
}

############################################################################################################################################

#Function for checking for, isolating, and appending species entries that were not present in the initial regression catagorizing matrix.
#Entries are mearly located and appended. The reg.vec column in the regMat file will retain all NA values for the taxa missing those entries
# unless the user changes or removes them via manual or automatic means.

checkMissingReg<- function(measReg, uniqTax) {
	missingReg <- measure.matReg[is.na(measure.matReg$reg.vec),]
	#merge/match with occs fiel to get order, genus, and family columns
	occurrences.order_name <- uniqTax$order[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.order_name)
	
	occurrences.family_name <- uniqTax$family[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.family_name)
	
	occurrences.genus_name <- uniqTax$genus[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.genus_name)
	
	missingReg <- missingReg[,c("occurrences.order_name","occurrences.family_name","occurrences.genus_name","species","reg.vec")]
	head(missingReg)
	colnames(missingReg)[colnames(missingReg) == "species"] <- "taxon"
	colnames(missingReg)[colnames(missingReg) == "reg.vec"] <- "reg"
	
	regMat <- rbind(regMat, missingReg)
	
	return(regMat) 
}

#############################################################################################################
#Get frequency of combinations for breaks ecology analysis
BreakComboFreq <- function(optList = NULL)
{
	optFrame <- matrix(nrow = 1000, ncol = 15)
	for(ii in seq(1, length(optList),1))
	{
		for(mm in seq(1,length(optList[[ii]][[length(optList[[ii]])]]$optBreaks),1))
		{
			optFrame[[ii,mm]] <- optList[[ii]][[length(optList[[ii]])]]$optBreaks[[mm]]
		}
	}
	return(count(optFrame))
}

#regimeHist is meant to parse through the repIntSp list object to print the distribution of traits
#within the regimes demarcated by breaks or designated by optList (if breaks are not provided).
#netFreq and regime Freq sets whether to use proportions or actually counts for the net gain/loss and regime 
#distribution histograms, respectively.
regimeHist <- function(repIntSp = NULL, breaks = NULL, optList, thisMat, netFreq = TRUE, regimeFreq=FALSE,
											 netPlotType = "absolute", plot.together = FALSE)
{
#function to generate histograms of species within each regime
#get list of intervals that comprise a regime
repIntSp_regimes <- repIntSp
if(is.null(breaks)) breaks <- optList_bm_allReps[[1]][[10]]$optBreaks

regimeBM.list <- list()
regimeSp.List <- list()

for(ii in seq(1, length(repIntSp_regimes),1))
{
	names(repIntSp_regimes[[ii]]) <- str_remove(names(repIntSp_regimes[[ii]]), " Ma")
}

	regimeSp <- unique(unlist(repIntSp_regimes[[1]][which(as.double(names(repIntSp_regimes[[ii]])) > breaks[1])]))
	regimeSp <- regimeSp[regimeSp %in% thisMat$species]
	regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
	
	regimeSp.List[[1]] <- regimeSp
	regimeBM.list[[1]] <- regimeBM
	
	rm(regimeSp, regimeBM)
	
for(mm in seq(2,length(breaks),1))
{
#get regimes for remaining sections
			
			#seq(min(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),
			#max(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),1))
		regimeSp <- unique(unlist(repIntSp_regimes[[1]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] 
																													& as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])]))
		regimeSp <- regimeSp[regimeSp %in% thisMat$species]
		regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
	
		regimeSp.List[[mm]] <- regimeSp
		regimeBM.list[[mm]] <- regimeBM
	
		rm(regimeSp, regimeBM)
}

	regimeSp <- unique(unlist(repIntSp_regimes[[1]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[length(breaks)])]))
	regimeSp <- regimeSp[regimeSp %in% thisMat$species]
	regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
	
	regimeSp.List[[length(breaks)+1]] <- regimeSp
	regimeBM.list[[length(breaks)+1]] <- regimeBM
	
	rm(regimeSp, regimeBM)
	
	regimeNames <- vector()
	regimeNames[1] <- paste(">",max(breaks,sep=""))
	for(ii in seq(2,length(regimeBM.list)-1,1))
	{
		regimeNames[ii] <- paste(paste(breaks[ii-1]," to ",sep=""),breaks[ii],sep="")
	}
	regimeNames[length(regimeBM.list)] <- paste(breaks[length(regimeBM.list)-1],">",sep="")
	names(regimeBM.list) <- regimeNames 
	bmBreaks <- c(0, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, 4.0)
	breakCol <- rainbow(length(bmBreaks))
	
	hist.list <- list()
	quartz()
	par(mfrow=c(2,length(regimeBM.list)/2+1))
	for(jj in seq(1, length(regimeBM.list),1))
	{
		hist.list[[jj]] <- hist(regimeBM.list[[jj]], main=names(regimeBM.list[[jj]]), breaks = bmBreaks, 
				 las=1, col = breakCol, freq = regimeFreq, xlab="log Body Mass (kg)", ylim = c(0,1))
	}
	
	regimeNetChange <- hist.list
	regimeNetChange[[length(regimeNetChange)]] <- NULL
	names(regimeNetChange) <- breaks
	
	for(tt in seq(2,length(hist.list),1))
	{
		regimeNetChange[[tt-1]]$counts <- hist.list[[tt]]$counts - hist.list[[tt-1]]$counts
		tt <- tt+1
	}
	
	if(netPlotType == "pos/neg")
	{
		#net change between regimes
		quartz()
		par(mfrow=c(2, length(hist.list)/2))
		for(hh in seq(1, length(regimeNetChange),1))
		{
		plot(regimeNetChange[[hh]],main=names(regimeNetChange[hh]), freq = netFreq,
				 las=1, col = breakCol, xlab="log Body Mass (kg)")  #####for some reason is identical to regular hist
		}
	}
	if(netPlotType == "absolute")
		{
		#absolute net change between regimes
		quartz()
		par(mfrow=c(2, length(hist.list)/2))
		for(dd in seq(1, length(regimeNetChange),1)) regimeNetChange[[dd]]$counts <- abs(regimeNetChange[[dd]]$counts)
		for(hh in seq(1, length(regimeNetChange),1))
		{
			plot(regimeNetChange[[hh]],main=names(regimeNetChange[hh]), freq = netFreq,
					 las=1, col = breakCol, xlab="log Body Mass (kg)")  #####for some reason is identical to regular hist
		}
	}
	
	####Get a hsitogram for each entry of repIntSp to avoid missing unique species in each regime
	####and then take median across all of the histograms.	
	
	return()
}

#regimeHist_countbox is meant to parse through the taxcube or countcube objects (matrices) to print the distribution 
#of traits within the regimes demarcated by breaks or designated by optList (if breaks are not provided).
#netFreq and regime Freq sets whether to use proportions or actually counts for the net gain/loss and regime 
#distribution histograms, respectively.
regimeHist_countBox <- function(countBox = NULL, breaks = NULL, optList, thisMat, netFreq = TRUE, regimeFreq=FALSE,
																netPlotType = "absolute", plot.together = FALSE, grayscale = TRUE, plot.axes = FALSE)
{
	if(class(countBox) == "matrix"){
		if(is.null(breaks)) breaks <- optList_bm_allReps[[1]][[length(optList_bm_allReps)]]$optBreaks
		
		regimeBM.list <- list()
		regimeSp.List <- list()
		
		colnames(countBox) <- str_remove(colnames(countBox), " Ma")
		
		regimeBM <- countBox[,which(as.double(colnames(countBox)) > breaks[1])]
		regimeBM.list[[1]] <- apply(regimeBM, c(1), sum)
		rm(regimeBM)
		
		for(mm in seq(2,length(breaks),1))
		{
			#get regimes for remaining sections
			
			#seq(min(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),
			#max(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),1))
			regimeBM <- countBox[,which(as.double(colnames(countBox)) < breaks[mm-1] 
																	& as.double(colnames(countBox)) > breaks[mm])]
			regimeBM.list[[mm]] <- apply(regimeBM, c(1), sum)
			
			rm(regimeBM)
		}  
		
		regimeBM <- countBox[,which(as.double(colnames(countBox)) < breaks[length(breaks)])]
		regimeBM.list[[length(breaks)+1]] <- apply(regimeBM, c(1),sum)
		
		rm(regimeBM)
		
		regimeNames <- vector()
		regimeNames[1] <- paste(">",max(breaks,sep=""))
		
		for(ii in seq(2,length(regimeBM.list)-1,1))
		{
			regimeNames[ii] <- paste(paste(breaks[ii-1]," to ",sep=""),breaks[ii],sep="")
		}
		regimeNames[length(regimeBM.list)] <- paste(breaks[length(regimeBM.list)-1],">",sep="")
		names(regimeBM.list) <- regimeNames 
		bmBreaks <- c(0, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, 4.0)
		
		btwnBreaks <- vector()
		for(ii in seq(1, length(bmBreaks)-1,1)) btwnBreaks[ii] <- (bmBreaks[ii]+bmBreaks[ii+1])/2 
		
		if(grayscale == TRUE) breakCol <- gray.colors(length(bmBreaks), start=0, end=1)
		else breakCol <- rainbow(length(bmBreaks))
		
		#make vector that conbtains the btwBreak values at qunatity listed by countBox
		for(hh in seq(1, length(regimeBM.list),1))
		{
			
			regimeBM.list[[hh]] <- c(rep(btwnBreaks[1],regimeBM.list[[hh]][1]),rep(btwnBreaks[2],regimeBM.list[[hh]][2]),
															 rep(btwnBreaks[3],regimeBM.list[[hh]][3]),rep(btwnBreaks[4],regimeBM.list[[hh]][4]),
															 rep(btwnBreaks[5],regimeBM.list[[hh]][5]),rep(btwnBreaks[6],regimeBM.list[[hh]][6]))
		}
		
		hist.list <- list()
		for(jj in seq(1, length(regimeBM.list),1))
		{
			
			hist.list[[jj]] <- hist(regimeBM.list[[jj]], breaks = bmBreaks, plot = FALSE)
			hist.list[[jj]]$density <- hist.list[[jj]]$counts/sum(hist.list[[jj]]$counts)
		}
		
		net.plot <- hist.list
		net.plot[[length(regimeNetChange)]] <- NULL
		names(net.plot) <- breaks
		
		for(tt in seq(2,length(hist.list),1))
		{
			net.plot[[tt-1]]$counts <- hist.list[[tt]]$counts - hist.list[[tt-1]]$counts
			#tt <- tt+1
		}
		
		max(sapply(net.plot, function(x) max(x$counts)))
		min(sapply(net.plot, function(x) min(x$counts)))
		#get limits for plot dimensions and scale
		hist.plot_ylimMax <- max(sapply(hist.list, function(x) max(x$density)))
		net.plot_ylimMax <- max(sapply(net.plot, function(x) max(x$counts)))
		net.plot_ylimMin <- min(sapply(net.plot, function(x) min(x$counts)))
		
		if(plot.axes == TRUE) marginPlot <- c(2,3,0,0)
		else marginPlot <- c(0,0.5,0,0)
		
		quartz(width = 11, height = 4)
		par(mfrow = c(1, length(hist.list)),mar=marginPlot)
		axesCheck <- FALSE
		for(ii in seq(1, length(hist.list),1))
		{
			if(plot.axes == TRUE)
			{
				if(ii == 1) axesCheck <- TRUE 
				else axesCheck <- FALSE
			}
			plot(hist.list[[ii]], col = breakCol, axes = axesCheck, ylim = c(0,1), main = NULL, freq = FALSE)
		}
		quartz(width = 11, height = 4)
		par(mfrow = c(1, length(net.plot)+1),mar=marginPlot) #have +1 to plot lenght so plots are same size as those for hist.plot
		for(ii in seq(1, length(net.plot),1))
		{
			if(plot.axes == TRUE)
			{
				if(ii == 1) axesCheck <- TRUE 
				else axesCheck <- FALSE
			}
			plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = axesCheck, ylim = c(net.plot_ylimMin,net.plot_ylimMax), 
					 main = NULL)
		}
	}
	
	return()
}
#############################################################################################################################
regimeHist_HistMedian<- function(repIntSp = NULL, breaks = NULL, optList, thisMat, netFreq = TRUE, regimeFreq=FALSE,
																 netPlotType = "absolute", plot.together = FALSE)
{
	repIntSp_regimes <- repIntSp
	if(is.null(breaks)) print("Error: Input a vector of break dates")
	
	regimeBM.list <- list()
	regimeSp.List <- list()
	all.hist <- list()
	
	for(ii in seq(1, length(repIntSp_regimes),1))
	{
		names(repIntSp_regimes[[ii]]) <- str_remove(names(repIntSp_regimes[[ii]]), " Ma")
	}
	
	for (ii in seq(1, length(repIntSp_regimes),1)) {
		regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) > breaks[1])]))
		regimeSp <- regimeSp[regimeSp %in% thisMat$species]
		regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
		
		regimeSp.List[[1]] <- regimeSp
		regimeBM.list[[1]] <- regimeBM
		
		rm(regimeSp, regimeBM)
		
		for(mm in seq(2,length(breaks),1))
		{
			#get regimes for remaining sections
			
			#seq(min(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),
			#max(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),1))
			regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] 
																														 & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])]))
			regimeSp <- regimeSp[regimeSp %in% thisMat$species]
			regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
			
			regimeSp.List[[mm]] <- regimeSp
			regimeBM.list[[mm]] <- regimeBM
			
			rm(regimeSp, regimeBM)
		}
		
		regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[length(breaks)])]))
		regimeSp <- regimeSp[regimeSp %in% thisMat$species]
		regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
		
		regimeSp.List[[length(breaks)+1]] <- regimeSp
		regimeBM.list[[length(breaks)+1]] <- regimeBM
		
		rm(regimeSp, regimeBM)
		
		regimeNames <- vector()
		regimeNames[1] <- paste(">",max(breaks,sep=""))
		for(nn in seq(2,length(regimeBM.list)-1,1))
		{
			regimeNames[nn] <- paste(paste(breaks[nn-1]," to ",sep=""),breaks[nn],sep="")
		}
		regimeNames[length(regimeBM.list)] <- paste(breaks[length(regimeBM.list)-1],">",sep="")
		names(regimeBM.list) <- regimeNames 
		bmBreaks <- c(0, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, 4.0)
		breakCol <- rainbow(length(bmBreaks))
		
		hist.list <- list()
		#quartz()
		#par(mfrow=c(2,length(regimeBM.list)/2+1))
		for(jj in seq(1, length(regimeBM.list),1))
		{
			hist.list[[jj]] <- hist(regimeBM.list[[jj]], main=names(regimeBM.list[[jj]]), breaks = bmBreaks, 
															las=1, col = breakCol, freq = regimeFreq, xlab="log Body Mass (kg)", plot = FALSE)
		}
		all.hist[[ii]] <- hist.list 
	}
	
	#get median across all hist counts and test to get median density
	median.density <- list()
	median.counts <- list()
	
	for(rr in seq(1, length(all.hist[[1]]),1)) #to go through each regime
	{
		regimeDensityMat <- matrix(nrow = length(all.hist),ncol= length(bmBreaks)-1)
		regimeCountMat <- matrix(nrow = length(all.hist),ncol= length(bmBreaks)-1)
		
		for(gg in seq(1, length(all.hist),1))
		{
			regimeDensityMat[gg,] <- all.hist[[gg]][[rr]]$density 
			regimeCountMat[gg,] <- all.hist[[gg]][[rr]]$counts
		}
		
		median.density[[rr]] <- apply(regimeDensityMat,2,median)
		median.counts[[rr]] <- apply(regimeCountMat,2,median)
	}
	hist.plot <- hist.list
	
	#overwite with median density and counts
	quartz(width = 11, height = 4)
	par(mfrow = c(1, length(hist.plot)),mar=c(0,0.5,0,0))
	for(pp in seq(1, length(hist.plot),1))
	{
		hist.plot[[pp]]$density <- median.density[[pp]]
		hist.plot[[pp]]$counts <- median.counts[[pp]]
		#plot(hist.plot[[pp]], col = breakCol, axes = FALSE)
	}
	
	#get net change between regimes
	net.plot <- hist.plot; net.plot[[length(hist.plot)]] <- NULL
	for(ii in seq(1, length(net.plot),1))
	{
		net.plot[[ii]]$counts <- hist.plot[[ii+1]]$counts-hist.plot[[ii]]$counts
		#plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = FALSE)
	}
	max(sapply(net.plot, function(x) max(x$counts)))
	min(sapply(net.plot, function(x) min(x$counts)))
	#get limits for plot dimensions and scale
	hist.plot_ylimMax <- max(sapply(hist.plot, function(x) max(x$density)))
	net.plot_ylimMax <- max(sapply(net.plot, function(x) max(x$counts)))
	net.plot_ylimMin <- min(sapply(net.plot, function(x) min(x$counts)))
	
	#plot both the hist and net changes
	quartz(width = 11, height = 4)
	par(mfrow = c(1, length(hist.plot)),mar=c(0,0.5,0,0))
	for(ii in seq(1, length(hist.plot),1))
	{
		plot(hist.plot[[ii]], col = breakCol, axes = TRUE, ylim = c(0,hist.plot_ylimMax), main = NULL)
	}
	quartz(width = 11, height = 4)
	par(mfrow = c(1, length(net.plot)),mar=c(0,0.5,0,0))
	for(ii in seq(1, length(net.plot),1))
	{
		plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = TRUE, ylim = c(net.plot_ylimMin,net.plot_ylimMax), 
				 main = NULL)
	}
	return() 
}

################################################################################################################
###Determine if range-through is working properly and then fill in the gaps if it isnt 
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
################################################################################################################
