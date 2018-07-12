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
	tree.backbone <- read.nexus("~/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
	clade.definitions <- read.csv("~/Dropbox/ungulate_RA/2017_3_20_Clade_species_test.csv", stringsAsFactors = FALSE)
	wildcard.positions <- read.csv("~/Dropbox/ungulate_RA/2017_4_17_MCRA_Codes.csv", stringsAsFactors = FALSE)
	regressCat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")
	setwd('~/Dropbox/ungulate_RA/RCode/Results')
} else if(Sys.info()["sysname"] == "Windows"){
#################PC   #had to swap my name with Blaires in pathnames
	source('C:/Users/Blaire/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call source for use on Evan's PC
	#setwd("/Users/Evan/Dropbox/ungulate_RA/RCode")
	tree.backbone <- read.nexus("/Users/Blaire/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
	clade.definitions <- read.csv("/Users/Blaire/Dropbox/ungulate_RA/2017_3_20_Clade_species.csv", stringsAsFactors = FALSE)
	wildcard.positions <- read.csv("/Users/Blaire/Dropbox/ungulate_RA/2017_4_17_MCRA_CodesAbsolute.csv", stringsAsFactors = FALSE)
	regressCat <- read.csv("/Users/Blaire/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")
	setwd("/Users/Blaire/Dropbox/ungulate_RA/RCode/Results")
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
save.files <- TRUE
if (!load.files) {
	
	occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
	occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
	occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
	#occs <- occs[!occs$order %in% c("Desmostylia", "Perissodactyla"), ]
	occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores
	
	
	print("Completed Occs")
	#####################################################################################################################################################
	
	#Compile Data for dental measurements, regression equation,
	# measureMat <- setMeasureReg(occs= occs, regMat= regressCat)
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
	
	#place function to request number of iterations
	#iter <- as.numeric(iterations()) #select number of iterations to run
	
	# tree.list <- list()
	#tree.list <- replicate(n = reps,expr = tree_resolution_v2(backbone_tree = tree.backbone,Species_file = clade.definitions,MRCA_file = wildcard.positions, occur = occs, thisMat = measureMat))
	# tree.list
	#quartz(width=14, height= 10)
	#	par(mfrow = c(1,reps))
	# jk <- 1
	# while(jk <= reps) {
		# print("jk")
		# print(jk)
		# resolvTime <- system.time(tree.current <- tree_resolution_v3(tree.backbone = tree.backbone, clade.definitions = clade.definitions, wildcard.positions = wildcard.positions, occs = occs, thisMat = measureMat))
		# tree.list[[jk]] <- tree.current
		# jk <- jk + 1
	# }
	print("Start tree generation")
	tree.list <- replicate(n=reps, expr=tree_resolution_v3(tree.backbone = tree.backbone, clade.definitions = clade.definitions, wildcard.positions = wildcard.positions, tree.list.wildcards=tree.list.wildcards, occs = occs, thisMat = measureMat), simplify=FALSE)
	if (reps>1) class(tree.list) <- "multiPhylo"
	
	check.wildcards <- FALSE
	if(check.wildcards){
		par(mfrow=c(5,2))
		#paintWildcardClades(tree.list[[ii]])
		for(ii in 1:length(tree.list)) {
			#print(length(tree.list[[ii]]$tip.label))
			paintWildcardClades(tree.list[[ii]])
			}
	}
	
	print("tree.list completed")

 	measureMatCrop <- matrix()
 	dat.vec <- data.frame()
 	#if(reps == 1) { measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[1]], thisMat = measureMat)
 	if (reps == 1) { measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[1]], thisMat = measureMat)
 									tree.list[[1]] <- comp_TreeSpecMat_outTree(this.tree = tree.list[[1]], thisMat = measureMat)
 	} else { 
 		for (ii in seq(1, length(tree.list))) { 
 			measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[ii]], thisMat = measureMat) 
 			tree.list[[ii]] <- comp_TreeSpecMat_outTree(this.tree = tree.list[[ii]], thisMat = measureMat)
 			}
 		}
	# measureMat[tree.list[[1]]$tip.label,] # drop taxa from the matrix not in the tree, also reorders matrix to match order in tree
	dat.vec <- as.data.frame(data.matrix(frame = measureMatCrop[,"bodyMass"]))
	dimnames(dat.vec) <- list(rownames(measureMatCrop), "body.mass")
	names(dat.vec) <- "body.mass"
	
	print("Dat vector completed")
}

if(save.files==TRUE){
	#write.tree(tree.list[[1]],file=paste("~/Dropbox/ungulate_RA/RCode/Results/NAUngualteTrees",timestamp,".tre",sep=""))
	#write.table(dat.vec, file="~/Dropbox/ungulate_RA/RCode/Results/DoughtyEvanSec1ALab2_NAUngulateTree_2018_3_8.txt", quote = TRUE, col.names = FALSE, row.names = TRUE)
	#write.csv(tree.dat,file="DoughtyEvanSec1ALab2_NAUngulateTree_test.csv")
	save(dat.vec, tree.list, file=paste(paste("~/Dropbox/ungulate_RA/RCode/Results/NAUngulateTrees_and_Data",timestamp(),sep=""),".Rdata",sep=""))
}

#attach $node.date to each tree and recaclucate root.length
if (load.files==TRUE) {
	#tree.list<- read.tree("DoughtyEvanSec1ALab2_NAUngulateTree_test.tre")
	load("~/Desktop/c230LAB/Project/DoughtyEvanSec1ALab2_NAUngulateTree_test_2018_3_8.Rdata")
}

modelType <- "temporal" # c("temporal", "clade")
#####################################################################
runSubclade <- FALSE
#Set the code below to subsample tree for specific clades
clade_select <- list(runSubclade=TRUE, clade_level = "order", cladename = "Perissodactyla")
if(runSubclade == TRUE) {
	clade_sp <- unique(occs[which(occs[clade_select$clade_level] == clade_select$cladename),c("accepted_name", "accepted_rank")])
	clade_sp <- clade_sp[which(clade_sp$accepted_rank == "species"),]
	clade_sp <- clade_sp[clade_sp$accepted_name %in% tree.list[[1]]$tip.label,]
	
	tree_clade <- extract.clade(phy=tree.list[[1]],node = getMRCA(phy = tree.list[[1]], clade_sp$accepted_name))
	tree.list[[1]] <- tree_clade
	#make sure that dat.vec is truncated before running evo model
}

####################################################################


analysis.settings <- list(do.BMsimple = TRUE, do.BMM = TRUE, 
													do.OUsimple = TRUE, do.OUM = TRUE, 
													do.heuristic = FALSE, adjust.date.after=10, 
													adjust.date.increment = 0.1, this.constraint = FALSE, 
													do.parallel = TRUE) #heurisitc does not currently work.  will break in current iteration due to intbreaks being removed from output Evan Doughty 7/11/2018

analysis.settingsClade = list(do.BMsimple = TRUE, do.BMM = TRUE, 
															do.OUsimple = TRUE, do.OUM = TRUE, 
															do.heuristic = FALSE,this.constraint = FALSE,
															do.parallel = FALSE)

MaxTime <- max(tree.list[[1]]$node.date) # setting this value to 80+ will cause error when running non-hueristic model for the earliest break (i.e. near the base of tree)
MinTime <- 5
int_length <-2
#intervals <- makeIntervals(MinTime, MaxTime, int_length)

###################################For SVP 2018 Abstract
#	intervals <- matrix(nrow=5, ncol=2); colnames(intervals) <- c("ageTop", "ageBase")
#	rownames(intervals) <- c(paste("42to",MaxTime,"Ma", sep=""), "34to42Ma","25to34Ma","16to25Ma","5to16Ma")
#	intervals <- as.data.frame(intervals)

#	intervals[1,] <- c(MinTime,16)	#MMCO 5 to 16
#	intervals[2,] <- c(16,25)	#Late Olig Warming 16 to 25
#	intervals[3,] <- c(25,34)	# E/O 25 to 34
#	intervals[4,] <- c(34, 42)	# Mid-Eocene CLiamtic Optimum 34 to 42
#	intervals[5,] <- c(42, 55)	# Early Eocene Cliamtic Optimum 42 to 52
	
#	intervals <- intervals[-5,]

intervals <- intervals.Radiant(breakDates = c(seq(5,55,3)),MinTime = MinTime,RootAge = MaxTime, int.Length = 0, int.Dist = 0)
	
this.tree <- tree.list[[1]]

tree.painted <- make.era.map(this.tree, limits = c(0,sort(MaxTime - intervals$ageBase)))
#col.ints <- c("black", "red", "green3", "blue","cyan")  ;names(col.ints) <- rownames(intervals)
#col.ints <- rainbow(nrow(intervals)) ;names(col.ints) <- rownames(intervals)
plotSimmap(tree.painted, fsize = 0.2)

resultsBM1 <- mvBM(tree = tree.painted, data = dat.vec, model = "BM1")
resultsBMM <- mvBM(tree = tree.painted, data = array(dat.vec[this.tree$tip.label,1], dimnames=list(this.tree$tip.label)), model = "BMM")
resultsBMM_brownie <- brownie.lite(tree = tree.painted, x = array(dat.vec[this.tree$tip.label,1], dimnames=list(this.tree$tip.label)))
resultsOU1 <- mvOU(tree = tree.painted, data = dat.vec, model = "OU1", param=list(root=TRUE))
resultsOU1_geiger <- fitContinuous(phy=tree.painted, dat=array(dat.vec[this.tree$tip.label,1], dimnames=list(this.tree$tip.label)), model="OU")
resultsOUM <- mvOU(tree = tree.painted, data = array(dat.vec[this.tree$tip.label,1], dimnames=list(this.tree$tip.label)), model = "OUM", param=list(root=TRUE))

eventResults <- list(resultsBM1, resultsBMM, resultsOU1, resultsOUM)

if(Sys.info()["sysname"] == "Darwin"){
save(eventResults, file = paste("~/Dropbox/ungulate_RA/RCode/Results/SVP_Abstract_", paste("#intervals=",nrow(intervals), sep=""),timestamp(), ".Rdata", sep=""))
} else if(Sys.info()["sysname"] == "Windows"){
	save(eventResults, file = paste("C:/Users/Blaire//Dropbox/ungulate_RA/RCode/Results/SVP_Abstract_", paste("#intervals=",nrow(intervals), sep=""),timestamp(), ".Rdata", sep=""))
}
###################################


print("Intervals completed")
###
#Need to impliment way to find and remove intervals that do not contain taxa
#Could retool code from ecology projject to determine how many taxa are in the specificed intervals and remove those that lack any taxa from the breaklist
###


treeOptList_Hueristic <- list()
treeOptList_NonHueristic <- list()

#save(tree.list, dat.vec, occs,measureMat, file = "/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModels_treelist2.RData")
# load(file = "/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModels_treelist.RData")

#tree.list <- tree.list[[1]]
#print(class(tree.list[[1]]))

cat("Model settings: \n","Number of Trees",reps,"\n Max Age of Intervals:", MaxTime, "\n Min Age of Intervals:", MinTime, "\n Interval Length:", int_length)
cat("analysis.settings:","\n do.BMsimple:", analysis.settings$do.BMsimple, "\n do.BMM:", analysis.settings$do.BMM, "\n do.OUsimple:", analysis.settings$do.OUsimple, "\n do.OUM:", analysis.settings$do.OUM, "\n do.heuristic", analysis.settings$do.heuristic, "\n adjust.date.after:", analysis.settings$adjust.date.after, "\n adjust.date.increment :", analysis.settings$adjust.date.increment, "\n this.constraint:",analysis.settings$this.constraint, "\n do.parallel:", analysis.settings$do.parallel)


	# print("Start evo modeling")
	# modelTime <- system.time(model.finalResults <- getTreeOptList(tree.list = tree.list, dat = dat.vec, intervals = intervals, analysis.settings = analysis.settings))
	# print("Models Completed")

for (this.tree in tree.list) {
	#print(this.tree)
	if(modelType == "temporal") {
	#All temporal shift model
	#tree.paintedTempAll <- paintTemporalRange(tree.list[[1]],maxAge = MaxTime, minAge = MinTime, int.length = int_length)
	#plotSimmap(tree.paintedTempAll, fsize = 0.05)
	
	#mvBM(tree=tree.paintedTempAll, data=dat.vec, model="BMM", optimization = "subplex", control=list(maxit =1e7))
	#mvOU(tree=tree.paintedTempAll, data=dat.vec, model="OUM", optimization = "subplex", control=list(maxit =1e7))
	
	#resultsTemporalAll <- 
	#mean(dat.vec)
	#Optimal Rate Shift
	resultsTemporalShifts <- list()
	resultsTemporalShifts <- testRateShiftsMVMorphIntervals(tree=this.tree, dat = dat.vec, intervals = intervals, analysis.settings = analysis.settings)
	print("Test Rate Shifts Done")
	filename_save_temp <- paste("~/Dropbox/ungulate_RA/RCode/Results/Temporal_NAUngulateTree_test_Hueristics=", 
															analysis.settings$do.heuristic,
															"_do.parallel=", analysis.settings$do.parallel,
															"_do.BMsimple=", analysis.settings$do.BMsimple,
															"_do.BMM=", analysis.settings$do.BMM, 
															"_do.OUsimple=", analysis.settings$do.OUsimple,
															"_do.OUMp=", analysis.settings$do.OUM,
															timestamp(),".RData", sep="")
	
	save(resultsTemporalShifts, file= filename_save_temp)
	}
	
	
	
	if(modelType == "clade") {
	resultClade <- orderEvoModel(tree = this.tree, dat = dat.vec, occs = occs, clade = "Artiodactyla", state.clade = "2", analysis.settings = analysis.settingsClade)
	filename_save_clade <- paste("~/Dropbox/ungulate_RA/RCode/Results/Clade_NAUngulateTree_test_clade=", 
															 clade_select, "_", timestamp(),".RData", sep="")
	
	save(resultClade, this.tree, tree.list, filename_save_clade)
	}
	#need to save or output list that contains all iterations of models#may need to input check for intervals to be included in filename
	# filename_save <- paste("/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModelResultsList_Hueristics=",analysis.settings$do.heuristic,sep="")
	#filename_save <- paste("~/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModelResultsList_Hueristics=", 
}

#load("/Users/emdoughty/Desktop/c230LAB/Project/C230Project_NAUngulateTree_test_Hueristics= FALSE _do.parallel= TRUE _do.BMsimple= TRUE _do.BMM= TRUE _do.OUsimple= TRUE _do.OUMp= TRUE ##------ Tue Feb 27 12:17:37 2018 ------## .RData")

#need to pull and open files from folder
filenames <- list.files("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/Results_Debugging", pattern="*.RData", full.names=TRUE)
optTrees <- list()
optBM <- list()
modelsType <- vector()
ratesCountAll <- vector()
ratesCountBM <- vector()
thetaCount <- vector()
OUsigma <- vector()
OUalpha <- vector()
BMoptdatesAll <- vector()
for(ii in seq(1, length(filenames),1)){
	print(ii)
	load(filenames[[ii]])
	resultsTemporalShifts
	optTrees[[ii]] <- getOptModels(opt.list = resultsTemporalShifts)
#	optBM[[ii]] <- getOptBM(opt.list = resultsTemporalShifts)
	
#	modelsType <- append(modelsType, optTrees[[ii]][[1]]$model) # three missing opt models are OU but don't have model parameter added (check generating function)
#	ratesCountAll <- append(ratesCountAll, optTrees[[ii]][[1]]$nrates)
#	ratesCountBM <- append(ratesCountBM, optBM[[ii]]$nrates)
#	BMoptdatesAll <- append(BMoptdatesAll, optBM[[ii]]$opt.dates)
		
#	if(optTrees[[ii]][[1]]$nrates == 1) {
#		thetaCount <- append(thetaCount, optTrees[[ii]][[1]]$theta)
#		OUalpha <- append(OUalpha, optTrees[[ii]][[1]]$alpha)
#		OUsigma <- append(OUsigma, optTrees[[ii]][[1]]$sigma)}
}

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




