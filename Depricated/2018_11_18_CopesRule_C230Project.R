#source("/Users/evandoughty/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_MainProg.R")  # for copy/paste torun in terminal
if(Sys.info()["sysname"] == "Darwin"){
	source("~/Dropbox/ungulate_RA/RCode/EvAnalysisMainSetup.R")
} else if(Sys.info()["sysname"] == "Windows"){
	source("C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysisMainSetup.R")
}
sourceCall()

treeFileList <- getTreeFiles()
	tree.backbone <- treeFileList$tree.backbone
	clade.definitions <- treeFileList$clade.definitions
	wildcard.positions <- treeFileList$wildcard.positions
	regressCat <- treeFileList$regressCat
	tree.list.wildcards <- treeFileList$tree.list.wildcards
		
#############################################################################################
###Arguments
reps <- 5
min.bl <- 1.0
load.files <- FALSE
MaxTime <- 57
MinTime <- 1
int_length <- 2
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
treeAnalysis <- "Regular" #("Regular, testDate, testTipLength) Regular is to run full analysis, testDate allow testing of dating method on single tree, tstTipLength will modify tips that are near the boundary within a specified value)
#############################################################################################
if (!load.files) {
	occs <- getOccs() 
	measureMat <- getMeasureMat(occs = occs, regMat = regressCat)
	#####################################################################################################################################################
	#Tree Generation/Compilation
	print("Start tree generation")
	####Regular Run    ##########Need to add toggles to swap between testing and regular tree generation
	if(treeAnalysis == "Regular"){
		tree.list <- replicate(n=reps, expr=tree_resolution_v3(tree.backbone = tree.backbone, 
																													 clade.definitions = clade.definitions, 
																													 wildcard.positions = wildcard.positions, 
																													 tree.list.wildcards=tree.list.wildcards, 
																													 occs = occs, thisMat = measureMat, 
																													 min.bl=min.bl), simplify=FALSE)
	}
	###############################################
	#####Test Tree Dating 
	if(treeAnalysis == "testDate"){
		tree.list <- tree_TestDate(tree.backbone = tree.backbone, clade.definitions = clade.definitions,
															 wildcard.positions = wildcard.positions, 
															 tree.list.wildcards=tree.list.wildcards, occs = occs, thisMat = measureMat, 
															 reps = 10, min.bl=min.bl)
		for(ii in seq(1, length(tree.list),1))  print(max(tree.list[[ii]]$node.date))
	}
	###############################################
	####Test Tree Structure
	###work in progress
	if(treeAnalysis == "testTipLength"){ print("Not yet implimented")}
	###############################################
	if (reps>1) class(tree.list) <- "multiPhylo"
	checkWildcards(tree.list=tree.list,check.wildcards = FALSE)
	print("tree.list completed")
	
	dat.vec <- getDat.vec(tree.list, measureMat)
}
####################################################################
###Get intervals
#setting this value to 80+ will cause error when running non-hueristic model for the earliest break (i.e. near the base of tree)
intervals <- makeIntervals(startDate = MinTime, endDate = MaxTime, intervalLength = int_length)
print("Intervals completed")
####################################################################
###Display Model Settings
cat("Model settings: \n","Number of Trees",reps,"\n Max Age of Intervals:", MaxTime, "\n Min Age of Intervals:", MinTime, "\n Interval Length:", int_length)
cat("analysis.settings:","\n do.BMsimple:", analysis.settings$do.BMsimple, "\n do.BMM:", analysis.settings$do.BMM, "\n do.OUsimple:", analysis.settings$do.OUsimple, "\n do.OUM:", analysis.settings$do.OUM, "\n do.heuristic", analysis.settings$do.heuristic, "\n adjust.date.after:", analysis.settings$adjust.date.after, "\n adjust.date.increment :", analysis.settings$adjust.date.increment, "\n this.constraint:",analysis.settings$this.constraint, "\n do.parallel:", analysis.settings$do.parallel)
####################################################################
###Run Analysis
#Set up code to time runs
startTime <- Sys.time()
for(this.tree in tree.list){
	#run MvMorph 
	resultsTemporalShifts <- testRateShiftsMVMorphIntervals(tree=this.tree, dat = dat.vec, intervals = intervals, analysis.settings = analysis.settings)
	print("Test Rate Shifts Done")
	
	#save each run and tree independently (give both the same timestamp)
	filename_save_temp <- paste("~/Dropbox/ungulate_RA/RCode/Results/TestDate/Data/BM_NAUngulateTree_test_Intervals=",
															int_length,"Ma_",
															"Hueristics=", analysis.settings$do.heuristic,
															"_do.parallel=", analysis.settings$do.parallel,sep="")
	time <- timestamp()
	filename_temp <- paste(filename_save_temp, time,".Rdata", sep="")
	
	#save(dat.vec, this.tree, intervals, file=paste(paste("~/Dropbox/ungulate_RA/RCode/NA_Ungulate_Trees/NAUngulateTrees_and_Data",time,sep=""),".Rdata",sep=""))
	save(dat.vec, this.tree, intervals, file=paste(paste("~/Dropbox/ungulate_RA/RCode/Results/TestDate/Trees/NAUngulateTrees_and_Data",time,sep=""),".Rdata",sep=""))
	save(resultsTemporalShifts, file= filename_temp)
}
endTime <- Sys.time()
RunTime <- endTime - startTime
RunTime     
####################################################################
###Search output files
#filenames.2Ma <- file.info(list.files("~/Google Drive/EvAnalysisResults20181017_dep/", pattern = "*.Rdata", full.names=TRUE))
#filenames.Trees <- file.info(list.files("~/Google Drive/EvAnalysisResults20181017_dep_Trees/", pattern = "*.Rdata", full.names=TRUE))

filenames.2Ma <- file.info(list.files("~/Dropbox/ungulate_RA/RCode/Results/TestDate/Data/", pattern = "*.Rdata", full.names=TRUE))
#filenames.Trees <- file.info(list.files("~/Dropbox/ungulate_RA/RCode/NA_Ungulate_Trees/", pattern = "*.Rdata", full.names=TRUE))

filenames.2Ma <- rownames(filenames.2Ma[with(filenames.2Ma, order(as.POSIXct(mtime))), ])
#filenames.Trees <- rownames(filenames.Trees[with(filenames.Trees, order(as.POSIXct(mtime))), ])

optBM <- getOptModelValues(fileList = filenames.2Ma, ExtractAttributes = c("sigma","break.dates"))

evoModelRates(evoResults = optBM, runOnRates = "BestRates")
evoModelHists(optBM, intervals, runHistOn = "BestRates", 
							model = "BM", plotPercent = TRUE, ylim=c(0, 0.15), relativeFreq = TRUE)
evoModelHists(optBM, intervals, runHistOn = "BestRates", 
							model = "BM", plotPercent = FALSE)
plotRatesBM(this.rez=optBM, intervals = intervals, num.Trees = length(optBM), 
						PlotXmax = 56, PlotXmin = 4)
RiseDeclineList <- lapply(optBM, getSigmaRateRiseDecline)
plotBreaksRiseDecline(RiseDeclineList, intervals)
plotBreaksRiseDecline(RiseDeclineList, intervals, plotPercent = TRUE)
####################################################################
