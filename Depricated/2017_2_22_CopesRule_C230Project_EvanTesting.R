#source("/Users/evandoughty/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_MainProg.R")  # for copy/paste torun in terminal

require(phytools)
require(mvMORPH)
require(stringr)
require(parallel)

#sources for Jon Marcot's code and specimen measurements 
source("https://dl.dropbox.com/s/8jy9de5owxj72p7/strat.R")
source("https://dl.dropbox.com/s/253p4avcvb66795/occFns.R")
source("https://dl.dropbox.com/s/9gdafsqss2b586x/phy_dateTree.R")
source("https://dl.dropbox.com/s/9tdawj35qf502jj/amandaSrc.R")
source("https://dl.dropbox.com/s/rlof7juwr2q4y77/blasto_Birlenbach.R")
source("https://dl.dropbox.com/s/643op7ye4s49w8p/utils_marcot.R")

if(Sys.info()["sysname"] == "Darwin"){
################MAC
	source('~/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call cource file for functions
  source("~/Dropbox/ungulate_RA/RCode/EvAnalysesTreeSrc.R")
	source("~/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R")
	source("~/Dropbox/ungulate_RA/RCode/EvAnalysesPlotSrc.R")
	tree.backbone <- read.nexus("~/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
	clade.definitions <- read.csv("~/Dropbox/ungulate_RA/2017_3_20_Clade_species_test.csv", stringsAsFactors = FALSE)
	wildcard.positions <- read.csv("~/Dropbox/ungulate_RA/2017_4_17_MCRA_Codes.csv", stringsAsFactors = FALSE)
	regressCat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")
	setwd('~/Dropbox/ungulate_RA/RCode/Results')
} else if(Sys.info()["sysname"] == "Windows"){
#################PC   #had to swap my name with Blaires in pathnames
	source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call source for use on Evan's PC
	source("C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysesTreeSrc.R")
	source("C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R")
	#setwd("/Users/Evan/Dropbox/ungulate_RA/RCode")
	tree.backbone <- read.nexus("/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
	clade.definitions <- read.csv("/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/2017_3_20_Clade_species.csv", stringsAsFactors = FALSE)
	wildcard.positions <- read.csv("/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/2017_4_17_MCRA_CodesAbsolute.csv", stringsAsFactors = FALSE)
	regressCat <- read.csv("/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")
	setwd("/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/Results")
} else print("Mac or Windows operating systems are not detected")

tree.backbone$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree.backbone$tip.label, fixed=TRUE)))
tree.list.wildcards <- list()
for (this.wildcard in seq_len(nrow(wildcard.positions))) {
	tree.list.wildcards[[this.wildcard]] <- read.nexus(file = wildcard.positions$Filename[this.wildcard])
	tree.list.wildcards[[this.wildcard]]$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree.list.wildcards[[this.wildcard]]$tip.label, fixed=TRUE)))
	#cat(this.wildcard, ":", length(tree.list.wildcards[[this.wildcard]]$tip.label),"\n")
}
class(tree.list.wildcards) <- "multiPhylo"

reps <- 10
load.files <- FALSE
save.files <- FALSE
if (!load.files) {
	
	occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
	occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
	occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
	#occs <- occs[!occs$order %in% c("Desmostylia", "Perissodactyla"), ]
	occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores
	
	
	print("Completed Occs")
	#####################################################################################################################################################
	
	measureMat <- getSingleSpeciesMatrix()
	print("getSingleSpMat completed")
	rownames(measureMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = measureMat$species)
	measureMat <- appendRegressionCategories(thisMat= measureMat, regMat= regressCat)
	print("regression append completed")
	#head(measureMat)
	measureMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = measureMat$species)
	#Approxiate body masses
	measureMat <-  approxBodyMass(thisMat = measureMat)
	print("bodymass approximated")
	#####################################################################################################################################################
	
	print("Start tree generation")
	tree.list <- replicate(n=reps, expr=tree_resolution_v3(tree.backbone = tree.backbone, clade.definitions = clade.definitions, wildcard.positions = wildcard.positions, tree.list.wildcards=tree.list.wildcards, occs = occs, thisMat = measureMat), simplify=FALSE)
	if (reps>1) class(tree.list) <- "multiPhylo"
	
	check.wildcards <- FALSE
	if(check.wildcards){
		par(mfrow=c(5,2))
		for(ii in 1:length(tree.list)) {
			paintWildcardClades(tree = tree.list[[ii]])
			}
	}
	
	print("tree.list completed")

 	measureMatCrop <- matrix()
 	dat.vec <- data.frame()
 	if (reps == 1) { measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[1]], thisMat = measureMat)
 									tree.list[[1]] <- comp_TreeSpecMat_outTree(this.tree = tree.list[[1]], thisMat = measureMat)
 	} else { 
 		for (ii in seq(1, length(tree.list))) { 
 			measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[ii]], thisMat = measureMat) 
 			tree.list[[ii]] <- comp_TreeSpecMat_outTree(this.tree = tree.list[[ii]], thisMat = measureMat)
 			}
 		}
	dat.vec <- as.data.frame(data.matrix(frame = measureMatCrop[,"bodyMass"]))
	dimnames(dat.vec) <- list(rownames(measureMatCrop), "body.mass")
	names(dat.vec) <- "body.mass"
	
	print("Dat vector completed")
}

if(save.files==TRUE){
	save(dat.vec, tree.list, file=paste(paste("~/Dropbox/ungulate_RA/RCode/Results/TreesLists/NAUngulateTrees_and_Data",timestamp(),sep=""),".Rdata",sep=""))
}

#attach $node.date to each tree and recaclucate root.length
if (load.files==TRUE) {
	#tree.list<- read.tree("DoughtyEvanSec1ALab2_NAUngulateTree_test.tre")
	load("~/Desktop/c230LAB/Project/DoughtyEvanSec1ALab2_NAUngulateTree_test_2018_3_8.Rdata")
}

modelType <- "temporal" # c("temporal", "clade")
#####################################################################

####################################################################


analysis.settings <- list(do.BMsimple = TRUE, do.BMM = TRUE, 
													do.OUsimple = TRUE, do.OUM = TRUE, 
													do.heuristic = TRUE, adjust.date.after=10, 
													adjust.date.increment = 0.1, this.constraint = FALSE, 
													do.parallel = FALSE) #heurisitc does not currently work.  will break in current iteration due to intbreaks being removed from output Evan Doughty 7/11/2018

analysis.settingsClade = list(do.BMsimple = TRUE, do.BMM = TRUE, 
															do.OUsimple = TRUE, do.OUM = TRUE, 
															do.heuristic = FALSE,this.constraint = FALSE,
															do.parallel = FALSE)

MaxTime <- max(tree.list[[1]]$node.date) # setting this value to 80+ will cause error when running non-hueristic model for the earliest break (i.e. near the base of tree)
MinTime <- 5
int_length <-2
#intervals <- makeIntervals(MinTime, MaxTime, int_length)

###################################For SVP 2018 Abstract
	intervals <- matrix(nrow=1, ncol=2); colnames(intervals) <- c("ageTop", "ageBase")
	rownames(intervals) <- "20MaBreak" #c(paste("42to",MaxTime,"Ma", sep=""), "34to42Ma","25to34Ma","16to25Ma","5to16Ma")
	intervals <- as.data.frame(intervals)

	intervals[1,] <- c(5,20)
		
#	intervals[1,] <- c(MinTime,16)	#MMCO 5 to 16
#	intervals[2,] <- c(16,25)	#Late Olig Warming 16 to 25
#	intervals[3,] <- c(25,34)	# E/O 25 to 34
#	intervals[4,] <- c(34, 42)	# Mid-Eocene CLiamtic Optimum 34 to 42
#	intervals[5,] <- c(42, 55)	# Early Eocene Cliamtic Optimum 42 to 52
	
#	intervals <- intervals[-5,]

intervals <- intervals.Radiant(breakDates = c(seq(5,55,3)),MinTime = MinTime,RootAge = MaxTime, int.Length = 0, int.Dist = 0)
	
#this.tree <- tree.list[[1]]

#tree.painted <- make.era.map(this.tree, limits = c(0,sort(MaxTime - intervals$ageBase)))
#plotSimmap(tree.painted, fsize = 0.2)

###################################

print("Intervals completed")

treeOptList_Hueristic <- list()
treeOptList_NonHueristic <- list()

#save(tree.list, dat.vec, occs,measureMat, file = "/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModels_treelist2.RData")
# load(file = "/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModels_treelist.RData")

cat("Model settings: \n","Number of Trees",reps,"\n Max Age of Intervals:", MaxTime, "\n Min Age of Intervals:", MinTime, "\n Interval Length:", int_length)
cat("analysis.settings:","\n do.BMsimple:", analysis.settings$do.BMsimple, "\n do.BMM:", analysis.settings$do.BMM, "\n do.OUsimple:", analysis.settings$do.OUsimple, "\n do.OUM:", analysis.settings$do.OUM, "\n do.heuristic", analysis.settings$do.heuristic, "\n adjust.date.after:", analysis.settings$adjust.date.after, "\n adjust.date.increment :", analysis.settings$adjust.date.increment, "\n this.constraint:",analysis.settings$this.constraint, "\n do.parallel:", analysis.settings$do.parallel)


	# print("Start evo modeling")
	# modelTime <- system.time(model.finalResults <- getTreeOptList(tree.list = tree.list, dat = dat.vec, intervals = intervals, analysis.settings = analysis.settings))
	# print("Models Completed")

filename_save_temp <- paste("~/Dropbox/ungulate_RA/RCode/Results/Temporal_NAUngulateTree_test_Hueristics=", 
														analysis.settings$do.heuristic,
														"_do.parallel=", analysis.settings$do.parallel,
														"_do.BMsimple=", analysis.settings$do.BMsimple,
														"_do.BMM=", analysis.settings$do.BMM, 
														"_do.OUsimple=", analysis.settings$do.OUsimple,
														"_do.OUMp=", analysis.settings$do.OUM,
														timestamp(),".RData", sep="")

count <- 1
resultsTemporalShifts <- list()
for (this.tree in tree.list) {
	#print(this.tree)
	if(modelType == "temporal") {

	print(count)
	
	#AIC ~ 1300
	resultsTemporalShifts[[count]] <- testRateShiftsMVMorphIntervals(tree=this.tree, dat = dat.vec, intervals = intervals, analysis.settings = analysis.settings)
	print("Test Rate Shifts Done")
	
	count <- count + 1
	
	save(resultsTemporalShifts, tree.list, file= filename_save_temp)
	}
	
	
	if(modelType == "clade") {
	resultClade <- orderEvoModel(tree = this.tree, dat = dat.vec, occs = occs, clade = "Artiodactyla", state.clade = "2", analysis.settings = analysis.settingsClade)
	filename_save_clade <- paste("~/Dropbox/ungulate_RA/RCode/Results/Clade_NAUngulateTree_test_clade=", 
															 clade_select, "_", timestamp(),".RData", sep="")
	
	save(resultClade, this.tree, tree.list, filename_save_clade)
	}
}

##########################################################################################################
#Post-Processing
##########################################################################################################
optTree <- list()
optBMOU <- list()

for(hh in seq(1,length(resultsTemporalShifts),1)){
cat(paste("Tree", paste(hh, "\n", sep=" "), sep=""))
print(getOptModelStatisticsFromOptList(resultsTemporalShifts[[hh]]))

#get opt for both bm and ou
optBMOU[[hh]] <- getOptModels(resultsTemporalShifts[[hh]], model = "BMOU")
}
	
#best model for each tree
evoModelRates(optBMOU, runOnRates = "BestRates")
evoModelHists(optBMOU, intervals, runHistOn = "BestRates", model = "BMOU")

#plot sigma for BM
this.rez <- optBMOU[[2]]
plot(this.tree, direction = "u", cex=0.05)
plotRateShiftsBM(this.rez, this.tree=tree.list[[1]], this.alpha=0.5, n.increments=10)

plotRatesBM(this.rez, this.tree=tree.list, num.Trees = length(tree.list))

plot(this.tree, direction = "u", cex=0.05)
plotRateShiftsOU(this.rez, this.tree=tree.list[[1]], this.alpha=0.5, n.increments=5)

#All entries in the opt.list
evoModelRates(resultsTemporalShifts, runOnRates = "AllRates")
evoModelHists(resultsTemporalShifts, intervals, runHistOn = "AllRates", model = "BMOU")

#plot sigma for BM
sigma.total <-  sapply(evoResults, function(y) y$BM$sigma)
names.vec <- sapply(evoResults[[1]]$BM$break.dates, function(y) y)

names(sigma.total) <- names.vec
boxplot(sigma.total)

#plot paramters for OU

#########################################################################################
#########################################################################################
#########################################################################################

#load("/Users/emdoughty/Desktop/c230LAB/Project/C230Project_NAUngulateTree_test_Hueristics= FALSE _do.parallel= TRUE _do.BMsimple= TRUE _do.BMM= TRUE _do.OUsimple= TRUE _do.OUMp= TRUE ##------ Tue Feb 27 12:17:37 2018 ------## .RData")

#need to pull and open files from folder
load("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/EvoAnalysis_Optlist_nrates_2##------ Fri Jul 13 14:27:06 2018 ------##.RData")
resultsTemporalShifts_2rate <- resultsTemporalShifts[[2]]
breakList2rate <- breakList
rez.list2rate <- rez.list
break.dates2rate <- sort(as.numeric(intervals[resultsTemporalShifts[[2]]$bm$BreakList.index, "ageBase"]), decreasing=TRUE)

intbreaks <- breakList2rate[resultsTemporalShifts[[2]]$bm$BreakList.index][[1]]
break.dates2rate <- sort(as.numeric(intervals[intbreaks, "ageBase"]), decreasing=TRUE)

this.limits <- max(this.tree$node.date) - c(max(this.tree$node.date), break.dates2rate)
tree.painted <- make.era.map(this.tree, limits = this.limits)
plotSimmap(tree.painted, fsize=0.01)

this.param.BM <- list(constraint = analysis.settings$this.constraint, smean = TRUE, trend = FALSE)    #default values to begin analysis

# AIC ~ 450
result <- mvBM(tree=tree.painted, data=dat.vec, model= "BMM",param=this.param.BM,
							 # method="rpf", control=list(maxit=1e8), optimization="subplex",
							 diagnostic=FALSE, echo=FALSE); result$AICc
# AIC ~ 450
master.BM <- doOneMultipleModel(this.tree = tree.painted, dat = dat.vec, intBreaks = NULL, break.dates = break.dates2rate, this.model = "BMM", analysis.settings = analysis.settings)
master.BM$AICc

# AIC ~450
topIntv <- min(which(intervals[,"ageBase"] > min(this.tree$node.date))) # node.date is non longer an attribute on the tree
baseIntv <- max(which(intervals[,"ageTop"] < max(this.tree$node.date[seq_along(this.tree$tip.label)])))
this.intervals <- intervals[topIntv:baseIntv,]
testBM1 <- getLkMVMorphOneIntervalSet(newBreak = breakList2rate[[3]], oldIntBreaks=NULL, this.tree = this.tree, dat = dat.vec, this.intervals, analysis.settings = analysis.settings)

# AIC ~1325
intervals <- matrix(nrow=2, ncol=2); colnames(intervals) <- c("ageTop", "ageBase")
rownames(intervals) <- c("20MaBreak","20to55") #c(paste("42to",MaxTime,"Ma", sep=""), "34to42Ma","25to34Ma","16to25Ma","5to16Ma")
intervals <- as.data.frame(intervals)

intervals[1,] <- c(5,20)
intervals[2,] <- c(20,MaxTime)

resAll <- testRateShiftsMVMorphIntervals(tree=this.tree, 
																				 dat = dat.vec, 
																				 intervals = intervals, 
																				 analysis.settings = analysis.settings)

###################################################################################################################################
#Check Outputs
###################################################################################################################################

filenames.OptTime <- file.info(list.files("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/Results_5MaBins/CompleteRun/", pattern="*.RData", full.names=TRUE))
filenames.IndivTime   <- file.info(list.files("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/Results_5MaBins/IndivRates/", pattern="*.RData", full.names=TRUE))

#filenames.Trees <- list.files()

#get files into order
filenames.OptList <- rownames(filenames.OptTime[with(filenames.OptTime, order(as.POSIXct(mtime))), ])
filenames.Indiv <- rownames(filenames.IndivTime[with(filenames.IndivTime, order(as.POSIXct(mtime))), ])

filenames.IndivIndex <- file.info(filenames.Indiv)

#need to create list that holds filenames for each tree
file.list <- list()

for(kk in seq(1, length(filenames.OptList),1)) {
	if(kk < 2) {
		file.tree <- filenames.Indiv[filenames.IndivIndex$mtime < filenames.OptTime$mtime[kk]] 
		names(file.tree) <- paste("Rates",seq(1, length(file.tree),1),sep="")
		file.list[[kk]] <- as.data.frame(file.tree)
		names(file.list[[kk]]) <- paste("CompletedFile", filenames.OptTime$mtime[1], sep="")
	} else {
		file.tree <- filenames.Indiv[filenames.IndivIndex$mtime > filenames.OptTime$mtime[kk-1] & filenames.IndivIndex$mtime < filenames.OptTime$mtime[kk]] 
		names(file.tree) <- paste("Rates",seq(1, length(file.tree),1),sep="")
		file.list[[kk]] <- as.data.frame(file.tree)
		names(file.list[[kk]]) <- paste("CompletedFile", filenames.OptTime$mtime[1], sep="")
	}
}

optTrees <- list()
optBM <- list()
optOU <- list()
optTreeList <- list()
optTreeModel <- list()
modelsType <- vector()
ratesCountAll <- vector()
ratesCountBM <- vector()
ratesCountOU <- vector()
thetaCount <- vector()
OUsigma <- vector()
OUalpha <- vector()
BMoptdatesAll <- vector()
break.datesAll <- vector()

OptMaster <- list()

for(ii in seq(1, length(file.list),1)){
	load(as.character(filenames.OptList[[ii]]))
	load(as.character(file.list[[ii]][1,])) #to get this.tree
	resultsTemporalShifts
	
	print(ii)

	param.list <- list(constraint = analysis.settings$this.constraint)
	OptResultBM <- mvBM(tree=this.tree, data=dat.vec, model="BM1", param=param.list, echo=FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
	OptResultBM$nrates <- 1
	OptResultBM$optBreaks <- NULL
	OptResultBM$model <- "BM"
	
	this.param <- list(sigma = NULL, alpha = NULL, vcv = "fixedRoot", decomp = c("cholesky", "spherical", "eigen", "qr", "diagonal", "upper", "lower")) 
	OptResultOU <- mvOU(tree=this.tree, data=dat.vec, model="OU1", param=param.list, echo=FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
	OptResultOU$nrates <- 1
	OptResultOU$optBreaks <- NULL
	OptResultOU$model <- "OU"
	
	OptMaster[[1]] <- list(bm=OptResultBM, ou=OptResultOU)
	
	for(jj in seq(2, length(file.list[[ii]][[1]]),1)){
	load(as.character(file.list[[ii]][jj-1,]))
	
	#rerun optimal model for rates
#	for(jj in seq(2, length(resultsTemporalShifts),1)){
		cat(jj,"\n")
		if(is.null(resultsTemporalShifts[[jj]]$bm)){
			OptResultBM <- NULL
		} else {Optbreaks <- resultsTemporalShifts[[jj]]$bm$BreakList.index
		
		intbreaks <- breakList[resultsTemporalShifts[[jj]]$bm$BreakList.index][[1]]
		break.dates <- sort(as.numeric(intervals[intbreaks, "ageBase"]), decreasing=TRUE)
		
		this.limits <- max(this.tree$node.date) - c(max(this.tree$node.date), break.dates)
		tree.painted <- make.era.map(this.tree, limits = this.limits)
		plotSimmap(tree.painted, fsize=0.01)
		
		this.param.BM <- list(constraint = analysis.settings$this.constraint, smean = TRUE, trend = FALSE)    #default values to begin analysis
		
		OptResultBM <- mvBM(tree=tree.painted, data=dat.vec, model= "BMM", param=this.param.BM,
									 # method="rpf", control=list(maxit=1e8), optimization="subplex",
									 diagnostic=FALSE, echo=FALSE)
		OptResultBM$nrates <- jj
		OptResultBM$optBreaks <- break.dates
		OptResultBM$model <- "BM"
		}
		
		if(is.null(resultsTemporalShifts[[jj]]$ou)){
			OptResultOU <- NULL
		} else {Optbreaks <- resultsTemporalShifts[[jj]]$ou$BreakList.index
		
		intbreaks <- breakList[resultsTemporalShifts[[jj]]$bm$BreakList.index][[1]]
		break.dates <- sort(as.numeric(intervals[intbreaks, "ageBase"]), decreasing=TRUE)
		
		this.limits <- max(this.tree$node.date) - c(max(this.tree$node.date), break.dates)
		tree.painted <- make.era.map(this.tree, limits = this.limits)
		plotSimmap(tree.painted, fsize=0.01)
		
		this.param.OU <- list(sigma = NULL, alpha = NULL, vcv = "fixedRoot", root=TRUE, decomp = c("cholesky", "spherical", "eigen", "qr", "diagonal", "upper", "lower")) 
		
		OptResultOU <- mvOU(tree=tree.painted, data=dat.vec, model= "OUM",
												param=this.param.OU, # method="rpf",
												control=list(maxit=1e8), 
												# optimization="subplex",
												diagnostic=FALSE, echo=FALSE)
		OptResultOU$nrates <- jj
		OptResultOU$optBreaks <- break.dates
		OptResultOU$model <- "OU"
		}
		
		OptMaster[[jj]] <- list(bm=OptResultBM, ou=OptResultOU)
	}
	
	optBM <- getOptModels(opt.list = OptMaster, model = "BM")
	optOU <- getOptModels(opt.list = OptMaster, model = "OU")	
	
	optTreeList[[ii]] <- list(bm=optBM, ou=optOU)
	
	if(optTreeList[[ii]]$bm$AICc < optTreeList[[ii]]$ou$AICc) optTreeModel[[ii]] <- optTreeList[[ii]]$bm
	
	if(optTreeList[[ii]]$ou$AICc < optTreeList[[ii]]$bm$AICc) optTreeModel[[ii]] <- optTreeList[[ii]]$ou
}

for(nn in seq(1, length(optTreeModel),1)){
	modelsType <- append(modelsType, optTreeModel[[nn]]$model)
	ratesCountAll <- append(ratesCountAll, optTreeModel[[nn]]$nrates)
	break.datesAll <- append(break.datesAll, optTreeModel[[nn]]$optBreaks)
}

#	modelsType <- append(modelsType, optTrees[[ii]][[1]]$model) # three missing opt models are OU but don't have model parameter added (check generating function)
#	ratesCountAll <- append(ratesCountAll, optTrees[[ii]][[1]]$nrates)
#	ratesCountBM <- append(ratesCountBM, optBM[[ii]]$nrates)
#	BMoptdatesAll <- append(BMoptdatesAll, optBM[[ii]]$opt.dates)
		
#	if(optTrees[[ii]][[1]]$nrates == 1) {
#		thetaCount <- append(thetaCount, optTrees[[ii]][[1]]$theta)
#		OUalpha <- append(OUalpha, optTrees[[ii]][[1]]$alpha)
#		OUsigma <- append(OUsigma, optTrees[[ii]][[1]]$sigma)}


multiRateOU <- which(ratesCount > 1)
length(multiRateOU)
optTrees[[multiRateOU[1]]][[1]]$opt.dates
optTrees[[multiRateOU[2]]][[1]]$opt.dates
optTrees[[multiRateOU[3]]][[1]]$opt.dates

quartz(6.18)
par(mfrow=c(1,2), oma = c(1,1,1,1))
hist(ratesCountAll,breaks = max(ratesCountAll)+1, xlim = c(1, max(ratesCountAll)))
hist(ratesCountBM, breaks = max(ratesCountBM)+1, xlim = c(1, max(ratesCountBM)))

quartz(6.18)
par(mfrow=c(1,2), oma = c(1,1,1,1))
hist(BMoptdatesAll, breaks = 70, xlim=c(1, 70))

#	if(optTrees[[ii]][[1]]$model == "OU") {
#		thetaCount <- append(thetaCount, optTrees[[ii]][[1]]$theta)
#		OUalpha <- append(ratesCount, optTrees[[ii]][[1]]$alpha)
#		OUsigma <- append(ratesCount, optTrees[[ii]][[1]]$sigma)
#	}
#	if(optTrees[[ii]][[1]]$model == "BM"){
#		BMsigma <- append(ratesCount, optTrees[[ii]][[1]]$sigma)
#	}

########################################################################################################################################
###Notes
########################################################################################################################################

#breaklist spliting idea
##master file 
#another fiel to read where we pull 10 or so to run
#set another column to workign status and have skip if true
#could be workaround for crashes or things that interupt 

#August 6-10th Jon be at teachign thing so cant do field


#Tara smiley
##Cohone pass work
##isotopes and paleo environment
##Oregan State under Rebecca Terry



##POTENTIAL REVIEWERS FOR MASTER PUB
#Ed Davis
#Julie Meachen
#Nick Famoso
#Richard Hulbert
#Larisa DeSantis
