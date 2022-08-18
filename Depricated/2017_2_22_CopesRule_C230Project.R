#source("/Users/evandoughty/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_MainProg.R")  # for copy/paste torun in terminal
if(Sys.info()["sysname"] == "Darwin"){
	source("~/Dropbox/ungulate_RA/RCode/EvAnalysisMainSetup.R")
} else if(Sys.info()["sysname"] == "Windows"){
	source("C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysisMainSetup.R")
}
sourceCall()

tree.list.wildcards <- GetWildcardTrees(wildcard.positions)

tree.backbone$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree.backbone$tip.label, fixed=TRUE)))
tree.list.wildcards <- list()


reps <- 1000
min.bl <- 1.0
load.files <- FALSE
# save.files <- TRUE
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
	tree.list <- replicate(n=reps, expr=tree_resolution_v3(tree.backbone = tree.backbone, clade.definitions = clade.definitions, wildcard.positions = wildcard.positions, tree.list.wildcards=tree.list.wildcards, occs = occs, thisMat = measureMat, min.bl=min.bl), simplify=FALSE)
	if (reps>1) class(tree.list) <- "multiPhylo"
	
	# check.wildcards <- FALSE
	# if(check.wildcards){
	# 	par(mfrow=c(5,2))
	# 	for(ii in 1:length(tree.list)) {
	# 		paintWildcardClades(tree.list[[ii]])
	# 		}
	# }
	# 
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



# if(save.files==TRUE){
	#save(dat.vec, tree.list, file=paste(paste("~/Dropbox/ungulate_RA/RCode/NA_Ungualte_Trees/NAUngulateTrees_and_Data",timestamp(),sep=""),".Rdata",sep=""))
# }
# 
# #attach $node.date to each tree and recaclucate root.length
# if (load.files==TRUE) {
# 	#tree.list<- read.tree("DoughtyEvanSec1ALab2_NAUngulateTree_test.tre")
# 	load("~/Desktop/c230LAB/Project/DoughtyEvanSec1ALab2_NAUngulateTree_test_2018_3_8.Rdata")
# }

modelType <- "temporal" # c("temporal", "clade")
#####################################################################

####################################################################


analysis.settings <- list(do.BMsimple = TRUE, do.BMM = TRUE, 
													do.OUsimple = FALSE, do.OUM = FALSE, 
													do.heuristic = TRUE, 
													trim.short.terminals = TRUE, min.terminal=0.5, 
													adjust.date.after=10, adjust.date.increment = 0.1, 
													this.constraint = FALSE, do.parallel = TRUE)

analysis.settingsClade = list(do.BMsimple = TRUE, do.BMM = TRUE, 
															do.OUsimple = TRUE, do.OUM = TRUE, 
															do.heuristic = FALSE,this.constraint = FALSE,
															do.parallel = FALSE)

# MaxTime <- max(tree.list[[1]]$node.date) # setting this value to 80+ will cause error when running non-hueristic model for the earliest break (i.e. near the base of tree)
MaxTime <- 57 # setting this value to 80+ will cause error when running non-hueristic model for the earliest break (i.e. near the base of tree)
MinTime <- 1
int_length <- 2
intervals <- makeIntervals(MinTime, MaxTime, int_length)
print("Intervals completed")

treeOptList_Hueristic <- list()
treeOptList_NonHueristic <- list()

cat("Model settings: \n","Number of Trees",reps,"\n Max Age of Intervals:", MaxTime, "\n Min Age of Intervals:", MinTime, "\n Interval Length:", int_length)
cat("analysis.settings:","\n do.BMsimple:", analysis.settings$do.BMsimple, "\n do.BMM:", analysis.settings$do.BMM, "\n do.OUsimple:", analysis.settings$do.OUsimple, "\n do.OUM:", analysis.settings$do.OUM, "\n do.heuristic", analysis.settings$do.heuristic, "\n adjust.date.after:", analysis.settings$adjust.date.after, "\n adjust.date.increment :", analysis.settings$adjust.date.increment, "\n this.constraint:",analysis.settings$this.constraint, "\n do.parallel:", analysis.settings$do.parallel)

#filename_save_temp <- paste("~/Dropbox/ungulate_RA/RCode/Results/BM_NAUngulateTree_test_Intervals=", 
filename_save_temp <- paste("~/Dropbox/ungulate_RA/RCode/Results/TestDate/BM_NAUngulateTree_test_Intervals=",
														int_length,"Ma_",
														"Hueristics=", analysis.settings$do.heuristic,
														"_do.parallel=", analysis.settings$do.parallel,sep="")
														#"_do.BMsimple=", analysis.settings$do.BMsimple,
														#"_do.BMM=", analysis.settings$do.BMM, 
														#"_do.OUsimple=", analysis.settings$do.OUsimple,
														#"_do.OUMp=", analysis.settings$do.OUM, sep="")

###############################################
#####Test Tree Dating 
tree.list <- tree_TestDate(tree.backbone = tree.backbone, clade.definitions = clade.definitions,
													 wildcard.positions = wildcard.positions, 
													 tree.list.wildcards=tree.list.wildcards, occs = occs, thisMat = measureMat, 
													 reps = 10, min.bl=min.bl)
for(ii in seq(1, length(tree.list),1))  print(max(tree.list[[ii]]$node.date))

####Test Tree Structure
###work in progress
###############################################

#set up code to time runs
startTime <- Sys.time()
for(this.tree in tree.list){

#set intervals 
	# MaxTime <- max(this.tree$node.date) # setting this value to 80+ will cause error when running non-hueristic model for the earliest break (i.e. near the base of tree)
	# MinTime <- 5
	# int_length <- 2
	# intervals <- makeIntervals(MinTime, MaxTime, int_length)
	# print("Intervals completed")

	filename_save_temp <- paste("~/Dropbox/ungulate_RA/RCode/Results/TestDate/BM_NAUngulateTree_test_Intervals=",
															int_length,"Ma_",
															"Hueristics=", analysis.settings$do.heuristic,
															"_do.parallel=", analysis.settings$do.parallel,sep="")
	#"_do.BMsimple=", analysis.settings$do.BMsimple,
	#"_do.BMM=", analysis.settings$do.BMM, 
	#"_do.OUsimple=", analysis.settings$do.OUsimple,
	#"_do.OUMp=", analysis.settings$do.OUM, sep="")
	
#run MvMorph 
resultsTemporalShifts <- testRateShiftsMVMorphIntervals(tree=this.tree, dat = dat.vec, intervals = intervals, analysis.settings = analysis.settings)
print("Test Rate Shifts Done")

#save each run and tree independently (give both the same timestamp)
time <- timestamp()
filename_temp <- paste(filename_save_temp, time,".Rdata", sep="")

#save(dat.vec, this.tree, intervals, file=paste(paste("~/Dropbox/ungulate_RA/RCode/NA_Ungulate_Trees/NAUngulateTrees_and_Data",time,sep=""),".Rdata",sep=""))
save(dat.vec, this.tree, intervals, file=paste(paste("~/Dropbox/ungulate_RA/RCode/Results/TestDate/NAUngulateTrees_and_Data",time,sep=""),".Rdata",sep=""))

save(resultsTemporalShifts, file= filename_temp)
}
endTime <- Sys.time()
RunTime2Ma <- endTime - startTime


#Time check
RunTime2Ma
#RunTime1Ma #closed 5:33-5:43 and 7:41:10 to 8:46

#search output files
#filenames.2Ma <- file.info(list.files("~/Google Drive/EvAnalysisResults20181017_dep/", pattern = "*.Rdata", full.names=TRUE))
#filenames.Trees <- file.info(list.files("~/Google Drive/EvAnalysisResults20181017_dep_Trees/", pattern = "*.Rdata", full.names=TRUE))

filenames.2Ma <- file.info(list.files("~/Dropbox/ungulate_RA/RCode/Results/TestDate/", pattern = "*.Rdata", full.names=TRUE))
#filenames.2Ma <- file.info(list.files("/Users/emdoughty/Google Drive/EvAnalysisResults20181023/", pattern = "*.Rdata", full.names=TRUE))
filenames.Trees <- file.info(list.files("~/Dropbox/ungulate_RA/RCode/NA_Ungulate_Trees/", pattern = "*.Rdata", full.names=TRUE))
#16
filenames.2Ma <- rownames(filenames.2Ma[with(filenames.2Ma, order(as.POSIXct(mtime))), ])
filenames.Trees <- rownames(filenames.Trees[with(filenames.Trees, order(as.POSIXct(mtime))), ])

#cbind(filenames.2Ma$mtime,filenames.Trees$mtime, filenames.2Ma$mtime - filenames.Trees$mtime)

#sub2Ma <- sapply(filenames.2Ma$mtime, substring, 1,16)
#subTree <- sapply(filenames.Trees$mtime, substring, 1,16)

Optbreaks <-vector()
OptSigma <- list()
optBM <- list()
Analysis.trees <- list()
Int.List <- list()

for(ii in seq(1, length(filenames.2Ma),1)){ #rerun and save data file for sigma and separate for break.dates
	load(as.character(filenames.2Ma[[ii]]))
#	load(as.character(filenames.Trees[[ii]]))
#	Analysis.trees[[ii]] <- this.tree
#	Int.List[[ii]] <- intervals
#	load(as.character(file.list[[ii]][1,])) #to get this.tree
#	resultsTemporalShifts
	
	bmBest <- getOptModels(opt.list = resultsTemporalShifts, model = "BM")
	optBM[[ii]] <-list(BM=bmBest[c("sigma","break.dates")], OU=NULL)
	
#	Optbreaks <- append(Optbreaks, optBM[[ii]]$break.dates)

	print(ii)
} #vector memory runs out at ~648 

#save sigma and breaks as separate files
#filenamesOptJon <- paste("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/OptBreaks_FinalData", timestamp(),".Rdata", sep="")
#filenamesOptJon <- paste("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/OptSigma_FinalData"
					#							 , timestamp(),".Rdata", sep="")

#save(optbreaks, file = filenamesOptJon)
#save(OptSigma, file = filenamesOptJon)

evoModelRates(evoResults = optBM, runOnRates = "BestRates")
evoModelHists(optBM, intervals, runHistOn = "BestRates", 
							model = "BM", plotPercent = TRUE, ylim=c(0, 0.15), relativeFreq = TRUE)

evoModelHists(optBM, intervals, runHistOn = "BestRates", 
							model = "BM", plotPercent = FALSE)

#plotRatesBM(this.rez=optBM, this.tree= Analysis.trees, intervals = Int.List, num.Trees = length(optBM))
plotRatesBM(this.rez=optBM, intervals = intervals, num.Trees = length(optBM), 
						PlotXmax = 56, PlotXmin = 4)

RiseDeclineList <- lapply(optBM, getSigmaRateRiseDecline)
plotBreaksRiseDecline(RiseDeclineList, intervals)
plotBreaksRiseDecline(RiseDeclineList, intervals, plotPercent = TRUE)


##############################################################################
#load("/Users/emdoughty/Desktop/c230LAB/Project/C230Project_NAUngulateTree_test_Hueristics= FALSE _do.parallel= TRUE _do.BMsimple= TRUE _do.BMM= TRUE _do.OUsimple= TRUE _do.OUMp= TRUE ##------ Tue Feb 27 12:17:37 2018 ------## .RData")
# 
# #need to pull and open files from folder
# load("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/EvoAnalysis_Optlist_nrates_2##------ Thu Jul 12 13:56:25 2018 ------##.RData")
# this.rez[[2]]
# breakList2rate <- breakList
# rez.list2rate <- rez.list
# break.dates2rate <- sort(as.numeric(intervals[this.rez[[2]]$bm$BreakList.index, "ageBase"]), decreasing=TRUE)
# 
# intbreaks <- breakList2rate[this.rez[[2]]$bm$BreakList.index][[1]]
# break.dates2rate <- sort(as.numeric(intervals[intbreaks, "ageBase"]), decreasing=TRUE)
# 
# this.limits <- max(this.tree$node.date) - c(max(this.tree$node.date), break.dates2rate)
# tree.painted <- make.era.map(this.tree, limits = this.limits)
# plotSimmap(tree.painted, fsize=0.01)
# 
# this.param.BM <- list(constraint = analysis.settings$this.constraint, smean = TRUE, trend = FALSE)    #default values to begin analysis
# 
# # AIC ~ 450
# result <- mvBM(tree=tree.painted, data=dat.vec, model= "BMM",param=this.param.BM,
# 							 # method="rpf", control=list(maxit=1e8), optimization="subplex",
# 							 diagnostic=FALSE, echo=FALSE)
# # AIC ~ 450
# master.BM <- doOneMultipleModel(this.tree = tree.painted, dat = dat.vec, intBreaks = NULL, break.dates = break.dates2rate, this.model = "BMM", analysis.settings = analysis.settings)
# 
# 
# load("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/EvoAnalysis_Optlist_nrates_3##------ Thu Jul 12 13:58:53 2018 ------##.RData")
# this.rez[[3]]
# breakList3rate <- breakList
# rez.list3rate <- rez.list
# 
# intbreaks <- breakList3rate[this.rez[[3]]$bm$BreakList.index][[1]]
# break.dates3rate <- sort(as.numeric(intervals[intbreaks, "ageBase"]), decreasing=TRUE); break.dates3rate
# 
# load("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/EvoAnalysis_Optlist_nrates_4##------ Thu Jul 12 14:08:42 2018 ------##.RData")
# this.rez[[4]]
# breakList4rate <- breakList
# rez.list4rate <- rez.list
# 
# intbreaks <- breakList4rate[this.rez[[4]]$bm$BreakList.index][[1]]
# break.dates4rate <- sort(as.numeric(intervals[intbreaks, "ageBase"]), decreasing=TRUE)
# #break.dates4rate <- sort(as.numeric(intervals[this.rez[[4]]$bm$BreakList.index, "ageBase"]), decreasing=TRUE)
# 
# 
# 
# filenames <- list.files("/Users/emdoughty/Dropbox/ungulate_RA/RCode/Results/Results_Debugging", pattern="*.RData", full.names=TRUE)
# optTrees <- list()
# optBM <- list()
# modelsType <- vector()
# ratesCountAll <- vector()
# ratesCountBM <- vector()
# thetaCount <- vector()
# OUsigma <- vector()
# OUalpha <- vector()
# BMoptdatesAll <- vector()
# for(ii in seq(1, length(filenames),1)){
# 	print(ii)
# 	load(filenames[[ii]])
# 	this.rez
# 	optTrees[[ii]] <- getOptModels(opt.list = this.rez)
# #	optBM[[ii]] <- getOptBM(opt.list = this.rez)
# 	
# #	modelsType <- append(modelsType, optTrees[[ii]][[1]]$model) # three missing opt models are OU but don't have model parameter added (check generating function)
# #	ratesCountAll <- append(ratesCountAll, optTrees[[ii]][[1]]$nrates)
# #	ratesCountBM <- append(ratesCountBM, optBM[[ii]]$nrates)
# #	BMoptdatesAll <- append(BMoptdatesAll, optBM[[ii]]$opt.dates)
# 		
# #	if(optTrees[[ii]][[1]]$nrates == 1) {
# #		thetaCount <- append(thetaCount, optTrees[[ii]][[1]]$theta)
# #		OUalpha <- append(OUalpha, optTrees[[ii]][[1]]$alpha)
# #		OUsigma <- append(OUsigma, optTrees[[ii]][[1]]$sigma)}
# }
# 
# multiRateOU <- which(ratesCount > 1)
# length(multiRateOU)
# optTrees[[multiRateOU[1]]][[1]]$opt.dates
# optTrees[[multiRateOU[2]]][[1]]$opt.dates
# optTrees[[multiRateOU[3]]][[1]]$opt.dates
# 
# quartz(6.18)
# par(mfrow=c(1,2), oma = c(1,1,1,1))
# hist(ratesCountAll,breaks = max(ratesCountAll)+1, xlim = c(1, max(ratesCountAll)))
# hist(ratesCountBM, breaks = max(ratesCountBM)+1, xlim = c(1, max(ratesCountBM)))
# 
# quartz(6.18)
# par(mfrow=c(1,2), oma = c(1,1,1,1))
# hist(BMoptdatesAll, breaks = 70, xlim=c(1, 70))
# 
# #	if(optTrees[[ii]][[1]]$model == "OU") {
# #		thetaCount <- append(thetaCount, optTrees[[ii]][[1]]$theta)
# #		OUalpha <- append(ratesCount, optTrees[[ii]][[1]]$alpha)
# #		OUsigma <- append(ratesCount, optTrees[[ii]][[1]]$sigma)
# #	}
# #	if(optTrees[[ii]][[1]]$model == "BM"){
# #		BMsigma <- append(ratesCount, optTrees[[ii]][[1]]$sigma)
# #	}
# 
# ########################################################################################################################################
# ###Notes
# ########################################################################################################################################
# 
# #breaklist spliting idea
# ##master file 
# #another fiel to read where we pull 10 or so to run
# #set another column to workign status and have skip if true
# #could be workaround for crashes or things that interupt 
# 
# #August 6-10th Jon be at teachign thing so cant do field
# 
# 
# #Tara smiley
# ##Cohone pass work
# ##isotopes and paleo environment
# ##Oregan State under Rebecca Terry
# 
