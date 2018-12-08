# # setwd("C:/Users/Evan/Dropbox/ungulate_RA/RCode")

# #all ungulates
#load("~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata")

# #artio only
# load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/artio/handleyResult##------ Wed Nov 15 16:57:31 2017 ------##_artiodactyla.Rdata")

# #perisso only
# load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/perisso/handleyResult##------ Thu Nov 16 10:27:58 2017 ------##.Rdata")
startTime <- Sys.time()
#sources for Jon Marcot's code and specimen measurements 
source("https://dl.dropbox.com/s/8jy9de5owxj72p7/strat.R")
source("https://dl.dropbox.com/s/253p4avcvb66795/occFns.R")g
source("https://dl.dropbox.com/s/9gdafsqss2b586x/phy_dateTree.R")
source("https://dl.dropbox.com/s/9tdawj35qf502jj/amandaSrc.R")
source("https://dl.dropbox.com/s/rlof7juwr2q4y77/blasto_Birlenbach.R")
source("https://dl.dropbox.com/s/pd5nmg1it13noc2/sampling.R") 
source("https://dl.dropbox.com/s/643op7ye4s49w8p/utils_marcot.R")
source("https://dl.dropbox.com/s/dozeb8o2pxu4sbj/CzTimescale.R") 
# source("C:/Users/Evan/Dropbox/ungulate_RA/RCode/isotopes.R")

if(Sys.info()["sysname"] == "Darwin"){
	source("~/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R", chdir = TRUE) #call cource file for functions
	source("~/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R", chdir = TRUE) #call cource file for functions
} else if(Sys.info()["sysname"] == "Windows"){
	source('C:/Users/Blaire/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R', chdir = TRUE) #call cource file for functions
}
#######
####### 
####### 

#occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Artiodactyla,Perissodactyla&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
#occs <- occs[occs$cc %in% c("US", "CA", "MX"), ]
#occs <- occs[!occs$order %in% c("Desmostylia", "Perissodactyla"), ]

occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]

# occs <- occs[occs$cc %in% c("US", "CA", "MX"), ]

#occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)

ranges <- getTaxonRangesFromOccs(occs=occs, random=TRUE)
rownames(ranges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(ranges))

int_length <- 1
intervals <- makeIntervals(1, 55.5, int_length)
intList <- listifyMatrixByRow(intervals)

####################################################################################################################################
print("Building measurement matrix...")

thisMat <- getSingleSpeciesMatrix()
thisMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = thisMat$species)

#regCatMat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv", na.strings=c("", "NA"), stringsAsFactors = TRUE)
regCatMat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv", na.strings=c("", "NA"), stringsAsFactors = TRUE)
thisMat <- appendMissingPaleoDBSpecies(thisMat, tax.vec=rownames(ranges))		# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files

print("Building body mass estimates...")
# thisMat <- setMeasureReg(occs= occs, regMat= regCatMat)
thisMat <- appendRegressionCategories(thisMat, regCatMat)
# nrow(thisMat)
# thisMat$species <- rownames(thisMat)
thisMat[!is.na(thisMat$family) & thisMat$family=="Entelodontidae", c("P2_L", "P3_L", "p2_w", "m2_w", "m3_w")] <- NA	### entelodont tooth widths were generating >5 ton body masses, so dropped here.
thisMat$bodyMass <- getBodyMassVectorFromThisMatAllMeasures(thisMat = thisMat, linked.files=TRUE)
thisMat$bodyMass <- fillMissingBodyMasses(thisMat)			# this fills taxa missing their body mass with the average body mass of its cogeners
# thisMat[!sapply(thisMat, is.finite)] <- NA
thisMat <- thisMat[is.finite(thisMat$bodyMass),]
# 2.2*(10^(thisMat$bodyMass[!is.na(thisMat$family) & thisMat$family=="Entelodontidae"]))
rownames(thisMat) <- thisMat$species
# length(rownames(thisMat))
####################################################################################################################################

# focal.order <- "Artiodactyla"
# focal.order <- "Perissodactyla"
focal.order <- c("Artiodactyla", "Perissodactyla")
bigList <- unique(occs[occs$accepted_rank =="species", c("order","family", "accepted_name")])
bigList <- bigList[order(bigList$order, bigList$family, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFam <- sort(unique(bigList$family[bigList$order %in% focal.order]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)

# matrix(thisMat$species[!thisMat$species %in% bigList$accepted_name[bigList$order %in% focal.order]], ncol=1)
thisMat <- thisMat[thisMat$species %in% bigList$accepted_name[bigList$order %in% focal.order], ]

####################################################################################################################################


########
#try getting rid of rows that have NA values in reg.vec column
########

	# richness <- matrix(0, nrow=nrow(intervals), ncol=4)
	# for (intv in seq_len((nrow(intervals)))) {
		# thisInt <- rownames(ranges)[ranges[,"FO"] >= intervals$ageTop[intv] & ranges[, "LO"] < intervals$ageBase[intv]]
		# richness[intv,1] <- length(thisInt)
		# richness[intv,2] <- sum(thisInt %in% rownames(thisMat))
		# richness[intv,3] <- sum(thisInt %in% rownames(thisMat)[is.finite(thisMat$bodyMass)])
		# richness[intv,4] <- sum(thisInt %in% rownames(thisMat)[is.finite(thisMat$PC3)])	
		# print(thisInt[!thisInt %in% rownames(thisMat)])	
	# }
	
	# par(mar=c(3,4,0.5,0.5), cex=0.66)
	# plot(rowMeans(intervals), richness[,3]/richness[,1], type="n", xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)), ylim=c(0,1), main="", xlab="Time (Ma)", ylab="% sampled")
	# overlayCzTimescale()
	# lines(rowMeans(intervals), richness[,2]/richness[,1], lwd=1.5, type="o", pch=15, col="firebrick")
	# lines(rowMeans(intervals), richness[,3]/richness[,1], lwd=1, type="o", pch=21, col="green4", bg="green1")
	# # lines(rowMeans(intervals), richness[,4]/richness[,1], lwd=1.5, type="o", pch=17, col="dodgerblue")

#runID <- 
do.parallel <- TRUE
	if (do.parallel) require(parallel)
reps <- 10
do.subsample <- TRUE
quota <- 0.4
do.disparity <- FALSE
bootstrapSpecimens <- FALSE
bootstrapSpecies <- FALSE
bootstrapSpeciesWithinIntervals <- FALSE
plotHist <- FALSE
do.heuristic = FALSE
	extra.intvs = 0
do.rangethrough <- TRUE

if (bootstrapSpecies) holderMat <- thisMat

if (plotHist) {
	quartz("Guild Histograms")
	par(mfrow=c((nrow(intervals)), 3), mar=c(0,0,0.75,0), cex.axis=0.5, cex.main=0.75)
}

repIntSp <- list()

######
# get species within intervals
#####

for (rep in seq_len(reps)) {
	cat("Beginning Rep", rep, "of", reps, "...\r")
	##################################################We need to update this sbootstrap section
	if (bootstrapSpecimens) {
		thisMat <- specimenMat[sample.int(nrow(specimenMat), size=nrow(specimenMat), replace=TRUE),]
		thisMat <- aggregate(thisMat, by=list(species=specimenMat$species), mean, na.rm=TRUE)
		# thisMat <- thisMat[,apply(!sapply(thisMat, is.na), 2, any)]
		rownames(thisMat) <- thisMat$species
		thisMat[sapply(thisMat, is.nan)] <- NA
		# thisMat<-cbind(thisMat, cbind(FO=vector(length=nrow(thisMat), mode="numeric"), LO=vector(length=nrow(thisMat), mode="numeric")))
		thisMat[,"reg"] <- as.character(famList$reg[match(thisMat$species,famList$taxon)])
		thisMat[,"bodyMass"] <- makeBodyMasses(thisMat, regList, best.only=TRUE)
		thisMat[,"PC2"] <- pcVec[match(thisMat$species, names(pcVec))]
		thisMat[,"PC3"] <- pcaLo$x[match(thisMat$species, rownames(pcaLo$x)),3]
	}
	if (bootstrapSpecies) thisMat <- holderMat[sample.int(n=nrow(thisMat), size=nrow(thisMat), replace=TRUE),]

	col.dates <- getCollectionAgesFromOccs(occs=occs[, c("collection_no", "max_ma", "min_ma")], random=TRUE)
	occDates <- col.dates$collection_age[match(occs$collection_no, col.dates$collection_no)]
	intOccs <- apply(intervals, 1, function(thisIntv) occs$occurrence_no[occDates > thisIntv[1] & occDates <= thisIntv[2]])
	# intTaxa <- sapply(intOccs, function(x) unique(occs$accepted_name[occs$occurrence_no %in% x]))
	# x <- intOccs
	intSp <- sapply(intOccs, function(x) match(sort(unique(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$accepted_rank =="species" & occs$occurrence_no %in% x]))), rownames(thisMat)))
	
	#which(occs$occurrence_no %in% x == TRUE) # none are being returned as TRUE
	
	if (do.subsample) { 
		nOccs <- sapply(intOccs, length)
	 	nTaxa <- sapply(intSp, length)
	 	nTaxa <- 0
		quota <- max(c(max(nTaxa), min(nOccs)))
		cat("Subsampling quota set to", quota, "occurrences")

		# intSp <- sapply(intOccs, function(x) match(sort(unique(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$occurrence_no %in% sample(x=x, size=quota)]))), rownames(thisMat)))
		intSp <- sapply(intOccs, function(x) sort(unique(as.character(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$occurrence_no %in% sample(x=x, size=quota)])))))
	}	else intSp <- sapply(intOccs, function(x) sort(unique(as.character(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$occurrence_no %in% x])))))

	########
	#range through method
	########
	
	if(do.rangethrough == TRUE) intSp <- makeRangeThroughOneRep(intSp)
	countInt <- intSp
	repIntSp[[rep]] <- intSp 
}

print("Completed getting species with intervals")
if(Sys.info()["sysname"] == "Darwin"){
	save(repIntSp, file=paste0("~/Dropbox/ungulate_RA/EcologyResults/RepIntSp_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
	#load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
} else if(Sys.info()["sysname"] == "Windows"){
	save(repIntSp, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/RepIntSp_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
	# load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
}


####################################################################################################################################
### Handley analysis of taxonomic distributions
print("Beginning median taxonomic Handley analysis...")

# bigList <- bigList[bigList$order %in% focal.order,]
# shortFam <- sort(unique(bigList$family))

taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
med.n <- median(sapply(repIntSp, function(x) length(unique(unlist(sapply(x, function(y) y))))))
optList_tax_median <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	

print("Beginning taxonomic Handley analysis for all reps...")
optList_tax_allReps <- list()
for (this.rep in seq_len(reps)) {
	taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
	this.n <- length(unique(unlist(sapply(repIntSp [[this.rep]], function(x) x))))
	optList_tax_allReps[[this.rep]] <- doHandleyTest(taxCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	
	if(this.rep %% 100 == 0) cat("Handley Rep:",this.rep, "\n")
}

####################################################################################################################################
if(Sys.info()["sysname"] == "Darwin"){
	save(repIntSp, optList_tax_median, optList_tax_allReps, file=paste0("~/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
	#load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
} else if(Sys.info()["sysname"] == "Windows"){
	save(repIntSp, optList_tax_median, optList_tax_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
	# load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
}

####################################################################################################################################
### Handley analysis of body mass distributions
print("Beginning median body mass Handley analysis...")

bmBreaks <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, Inf) #Janis 2000  max(thisMat$bodyMass, na.rm=TRUE)
# bmBreaks <- c(-Inf, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, Inf) #Badgely and Fox 2000
# bmBreaks <- hist(thisMat$bodyMass, plot=FALSE)$breaks

countCube <- sapply(repIntSp, function(y) sapply(y, function(x) hist(thisMat$bodyMass[match(x, rownames(thisMat))], breaks=bmBreaks, plot=FALSE)$counts), simplify = "array")
countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

optList_bm_median <- doHandleyTest(countBox[2,,], n=nrow(thisMat), do.heuristic=do.heuristic, extra.intvs=extra.intvs)
	# quartz(width=3.3, height=9.8)
	# par(mfrow=c(ncol(countBox[2,,])/2,2), mar=c(0,3,0.5, 0.5), mfg=c(2,1))
	# for (i in seq(from=1, to=ncol(countBox[2,,]), by=1)) barplot(countBox[2,,][,i], width=c(0.68897, 0.68897, 0.778151, 0.522879, 1.18786),space=0, cex.axis=0.5, ylim=c(0,30))
	
	# quartz()
	# thisTab <- table(unlist(sapply (optList_bm_allReps, function(x) x[[(length(x) - 1)]]$optBreaks)))
	# thisTab <- array(thisTab[match(seq_len(nrow(intervals)), names(thisTab))], dimnames=list(rownames(intervals)))
	# barplot(rev(thisTab)/reps, cex.names=0.5, ylim=c(0,1))
	# abline(h=c(0.95, 0.75, 0.5), lty=3, col="gray50")

print("Beginning body mass Handley analysis for all reps...")
optList_bm_allReps <- list()
for (this.rep in seq_len(reps)) {
	this.n <- length(unique(unlist(sapply(repIntSp [[this.rep]], function(x) x))))
	optList_bm_allReps[[this.rep]] <- doHandleyTest(countCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)
	if(this.rep %% 100 == 0) cat("Handley Rep:",this.rep, "\n")
}

####################################################################################################################################
if(Sys.info()["sysname"] == "Darwin"){
	save(repIntSp, optList_bm_median, optList_bm_allReps, file=paste0("~/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
	#load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
	} else if(Sys.info()["sysname"] == "Windows"){
		save(repIntSp, optList_bm_median, optList_bm_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
		# load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
	}

#save(repIntSp, optList_tax_median, optList_tax_allReps, optList_bm_median, optList_bm_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/handleyResult", timestamp(),".Rdata"))


# load('~/Dropbox/ungulate_RA/EcologyResults/handleyResult##------ Sun May 13 01:32:21 2018 ------##.Rdata')
####################################################################################################################################

	### Handley-block histogram series
	this.rep <- 1
	old.shift <- nrow(intervals)
	shift.ints <- optList_bm_median[[length(optList_bm_median) - 1]]$optBreaks
	# shift.ints <- rev(which(intervals$ageBase %in% optList_bm_allReps[[this.rep]][[length(optList_bm_allReps[[this.rep]]) - 1]]$optBreaks))
	par(mfrow=c(1,length(shift.ints)+1), mar=c(4,2,0.5,0.5), col.axis="gray50", col.lab="gray50", fg="gray50")
	hist.breaks <- c(min(thisMat$bodyMass, na.rm=TRUE),bmBreaks[2:6],max(thisMat$bodyMass, na.rm=TRUE))
	for (this.shift in c(shift.ints, 0)) {
		hist(thisMat$bodyMass[unique(unlist(repIntSp[[this.rep]][seq(from=old.shift, to=(this.shift+1))]))], freq=FALSE, breaks=hist.breaks, col=rainbow(length(hist.breaks)), main="", xlab="log Body Mass", ylab="", ylim=c(0,1))
	}

	quartz()
	par(mfrow=c(2,1), mar=c(5, 4, 4, 1))
	hist(sapply(optList_tax_allReps, function(x) length(x) - 2), breaks=seq(-0.5, 10, 1.0), col="orchid4", main="Number of taxonomic shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))
	hist(sapply(optList_bm_allReps, function(x) length(x) - 2), breaks=seq(-0.5, 10.5, 1.0), col="firebrick4", main="Number of body mass shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))
	
quants <- apply(sapply(repIntSp, function(y) sapply(y, function(x) quantile(thisMat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
####################################################################################################################################
	### number of Replicates with a shift in that interval
	quartz()
	par(mfrow=c(2,1), mar=c(4, 4, 1, 1))
	# optList_tax_allReps <- optList_tax_allReps_heuristic
	# optList_tax_allReps <- optList_tax_allReps_full
	breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_tax_allReps, function(x) x[[length(x) - 1]]$optBreaks))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=rev(range(intervals)), ylim=c(0,reps), main="Number of Replicates with a Taxonomic Distribution Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(breakHist, col="orchid4", border="orchid1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=c(55, 0), ylim=c(0,reps), add=TRUE)

	# par(mfrow=c(2,1), mar=c(3.5, 3.5, 1, 1))
	# optList_bm_allReps <- optList_bm_allReps_heuristic
	# optList_bm_allReps <- optList_bm_allReps_full
	### number of Replicates with a shift in that interval
	breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_bm_allReps, function(x) unique(x[[length(x) - 1]]$optBreaks)))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=rev(range(intervals)), ylim=c(0,reps), main="Number of Replicates with a Body Mass Distribution Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(breakHist, col="firebrick4", border="firebrick1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=c(55, 0), ylim=c(0,reps), add=TRUE)


####################################################################################################################################


####################################################################################################################################
######
# make three-panel figure
#####
	quartz(width=6.89)
		par(mfrow=c(4,1), mar=c(0,4,0.5,0.5), mgp=c(2, 1,0))
		
	### isotope panel
		if(Sys.info()["sysname"] == "Darwin"){
				source("~/Dropbox/ungulate_RA/RCode/isotopes.R")
			} else if(Sys.info()["sysname"] == "Windows"){
				source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
			}
		
	#source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
	# source("~/Dropbox/code/R/common_src/isotopes.R")
		
	optList_topes <- doTopesRateAnalysis(intervals)
	plotTopesRateAnalysis(optList_topes, intervals, x.axis=FALSE) #
	box(lwd=1)
	getAlroyStatistics(intervals)
			
	### taxonomy panel
		par(mar=c(0,4, 2.5,0.5))
		
		# taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% rownames(thisMat)[x]], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
		taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
		dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
		
		# prop <- t(apply(taxCube, c(1,2), median, na.rm=TRUE))
		prop <- t(apply(taxCube, c(1,2), mean, na.rm=TRUE))
		colnames(prop)[colnames(prop)==""] <- "indeterminate"
		# dimnames(prop) <- list(rownames(intervals), shortFam)
		source("https://dl.dropbox.com/s/iy0tu983xesbig2/taxonomicEv.R")
		plotStackedRichness(this.box=prop, intervals=intervals, do.log=FALSE, overlay.labels=TRUE, numbers.only=TRUE, legend=FALSE, xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)))
		#med.n <- median(length(unique(unlist(sapply(repIntSp[[this.rep]], function(x) rownames(thisMat)[x]))))) #what is this.rep set to during this function?  variable is used in for loop in Handley
		# med.n <- median(sapply(repIntSp, function(x) length(unique(unlist(sapply(x, function(y) rownames(thisMat)[y]))))))
		# optList_tax <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
		abline(v=sort(c(intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col="darkorchid4")
		text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_tax_median[[length(optList_tax_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="darkorchid4")
		text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0), cex=0.5, col="darkorchid4")
		box(lwd=1)
	
	### body mass panel
		thisRanges <- getTaxonRangesFromOccs(occs=occs, random=FALSE)
		rownames(thisRanges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisRanges))
		thisMat[,c("FO","LO")] <- thisRanges[match(rownames(thisMat), rownames(thisRanges)),]

		par(mar=c(0,4,2.5,0.5))
		# quartz(width=12, height=6)
		plot(thisMat$FO, thisMat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(55,5,10), xlab="Time (Ma)", cex.axis=1.5, cex.lab=1.5)
		# plot(thisMat$FO, thisMat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(50,0,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, col="gray75", fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #alter xaxpto change x-axis values
		# rect(-10e6, -10e6, 10e6, 10e6, col="white")
		overlayCzTimescale(do.subepochs=TRUE)
		
		famColors <- rainbow(length(shortFam))
		colorList <- famColors[match(bigList$family[as.character(bigList$accepted_name) %in% rownames(thisMat)], shortFam)]
		colorList[is.na(colorList)] <- "gray25"
		
		orderColors <- array(NA, dim=nrow(thisMat))
		# orderColors[bigList$order[match(rownames(thisMat), bigList$accepted_name)]=="Perissodactyla"] <- "dodgerblue4"
		# orderColors[bigList$order[match(rownames(thisMat), bigList$accepted_name)] =="Artiodactyla"] <- "deeppink4"

		for (i in seq_len(nrow(thisMat))) {
			# lines(x=c(this["FO"], x["LO"]), y=c(x["bodyMass"], x["bodyMass"]), lwd=3, pch=21, col=famColors[match(bigList[match(rownames(thisMat), bigList[,1]),2], shortFam)])
			# lines(x=c(thisMat$FO[i], thisMat$LO[i]), y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor(colorList[i], 0.75))
			# lines(x=c(thisRanges[match(rownames(thisMat)[i], rownames(thisRanges)),"FO"], thisRanges[match(rownames(thisMat)[i], rownames(thisRanges)),"LO"]), y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor("gray0", 0.75)) #
			if (is.finite(thisMat$FO[i]) & is.finite(thisMat$LO[i]) & thisMat$FO[i] != thisMat$LO[i]) lines(x=thisMat[i,c("FO","LO")], y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.75, pch=21, col=alphaColor("gray0", 0.5)) #alphaColor(orderColors[i], 0.5)
		}
		points(thisMat[complete.cases(thisMat[ ,c("FO","LO")]) & thisMat$FO == thisMat$LO, c("FO","bodyMass")], pch=21, col=alphaColor("gray0", 0.5), cex=0.25) #this line is not generating the proper output for the final graph due to c("FO","bodyMass") causing a  "undefined columns selected" error
		
		# optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
		# optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on median
		abline(v=sort(c(intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col="firebrick4")
		text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_bm_median[[length(optList_bm_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="firebrick4")
		text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0),cex=0.5, col="firebrick4")
		
		quants <- apply(sapply(repIntSp, function(y) sapply(y, function(x) quantile(thisMat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
		polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[1,], rev(quants[5,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
		polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[2,], rev(quants[4,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
		lines(rowMeans(intervals), quants[3,], col=alphaColor("goldenrod1", 0.5), lwd=5)
		lines(rowMeans(intervals), quants[3,], col=alphaColor("darkorange4", 1.0), lwd=3)
		points(rowMeans(intervals), quants[3,], col=alphaColor("darkorange1", 0.5), cex=0.5)
		box(lwd=1)

		endTime <- Sys.time()
		RunTime2Ma <- endTime - startTime

# write function to lot distribtion of body mass within each break interval
#Jons code for histograms below: run tests to see how they change iteration to iteration#make into function 
#	this.rep <- 1
 #   old.shift <- nrow(intervals)
  #  shift.ints <- rev(which(intervals$ageBase %in% optList[[this.rep]][[length(optList[[this.rep]])]]$optBreaks))
  #  par(mfrow=c(1,length(shift.ints)+1), mar=c(4,2,0.5,0.5))
  #  for (this.shift in c(shift.ints, 0)) 
  #  {
  #  	hist(thisMat$bodyMass[unique(unlist(repIntSp[[this.rep]][seq(from=old.shift, to=(this.shift+1))]))], freq=TRUE, breaks=c(min(thisMat$bodyMass, na.rm=TRUE),bmBreaks[2:6],max(thisMat$bodyMass, na.rm=TRUE)), col=rainbow(n=length(shift.ints)+1), main="", xlab="log Body Mass", ylab="", ylim=c(0,1))
  #  }


###
#Compiled Histogram of Proportional shifts in body mass
###

# # install.packages("ggplot2")
# require(ggplot2)
# install.packages("Rcmdr")
	# require(Rcmdr)
# quartz(width=6.89)
	# par(mfrow=c(2,1), mar=c(4,4, 1,0.5), mgp=c(2, 1,0))
		# #plot(thisMat$FO, thisMat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="% Species Diversity", xaxp =c(55,5,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, fg="black", bg="black", col.axis="black", col.lab="black")
		# #rect(-10e6, -10e6, 10e6, 10e6, col="white")
		# #overlayCzTimescale(do.subepochs=TRUE)
		
		# #assign species bin from Janis 2000
		# #need a vector like countCube but made with raw data prior to distributions
		# #get list of taxa throughout all intervals by running single run of species interval and rangethrough
		# countInt
				
		# #load("/Users/evandoughty/Dropbox/ungulate_RA/RCode/EcologyAnalysisResults/countInt.RData")
		# #countInt <- countInt_bin1myr
		# #countInt <- countInt_bin2myr		
			# y <- countInt
		# #countDist <- sapply(countInt, function(y) sapply(y, function(x) hist(this.column[x], breaks=breaks, plot=FALSE)$counts), simplify = "array")
		# countDist <- sapply(countInt, function(x) hist(this.column[x], breaks=breaks, plot=FALSE)$counts)
		# #dimnames(countDist)
		
		# r <- matrix()
		# r <- countDist
		# colnames(r) <- intervals$ageTop
		# a <- r[,ncol(r):1]
		# r <- a
				
		# #Need to Remove Intervals from Handley Bins and into actual intervals
		# cols <- palette(rainbow(6))
		# distBar <- barplot(r, ylim=c(0,100), xlab="Time(Ma)", ylab = "Number of Species", col =cols, space = 0)

	# par(mar=c(6,4,2.5,0.5))	
		# summedDist <- colPercents(r)
		# sumBar <- rev(barplot(summedDist[-c(7:9),], col = cols, xlab="Time(Ma)", ylab = "% Species Diversity"))
		
	
# ####
# #NMDS
# ###
# install.packages("vegan")
# require(vegan)
# bmBox <- sapply(repIntSp, function(thisRep) t(sapply(thisRep, function(taxa) hist(thisMat$bodyMass[taxa], freq=TRUE, xlim=c(min(bmBreaks), max(bmBreaks)), ylim=c(0,22), breaks=bmBreaks, col="firebrick1", ylab="", plot=plotHist)$counts)), simplify="array") #, main=paste(intervals[intv, "ageTop"], "-", intervals[intv, "ageBase"], "Ma (n = ", length(intSp[[intv]]), ")")
# par(mfrow=c(2,1), mar=c(3,3,0.5,0.5), mgp=c(1.25, 0.5, 0))
# tax <- plotNMDS(thisBox=t(apply(taxCube, c(1,2), mean, na.rm=TRUE)), intervals, polygon.ints=optList_tax[[(length(optList_tax)-1)]]$optBreaks, scaler=5, title="taxonomy") #, filename="~/Desktop/NMDS_bodyMass.pdf"
# bm <- plotNMDS(thisBox=bmBox, intervals, polygon.ints=oneOpt[[(length(oneOpt)-1)]]$optBreaks, scaler=5, title="Body Mass") #, filename="~/Desktop/NMDS_bodyMass.pdf"
