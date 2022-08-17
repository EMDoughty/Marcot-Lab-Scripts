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

source('~/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call cource file for functions

tree.backbone <- read.nexus("~/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
clade.definitions <- read.csv("~/Dropbox/ungulate_RA/2017_3_20_Clade_species.csv", stringsAsFactors = FALSE)
wildcard.positions <- read.csv("~/Dropbox/ungulate_RA/2017_4_17_MCRA_Codes.csv", stringsAsFactors = FALSE)
regressCat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")

#PC
#source('C:/Users/Evan/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call source for use on Evan's PC
#setwd("/Users/Evan/Dropbox/ungulate_RA/RCode")

#tree.backbone <- read.nexus("/Users/Evan/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
#clade.definitions <- read.csv("/Users/Evan/Dropbox/ungulate_RA/2017_3_20_Clade_species.csv", stringsAsFactors = FALSE)
#wildcard.positions <- read.csv("/Users/Evan/Dropbox/ungulate_RA/2017_4_17_MCRA_Codes.csv", stringsAsFactors = FALSE)
#regressCat <- read.csv("/Users/Evan/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")

reps <- 1

tree.backbone$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree.backbone$tip.label, fixed=TRUE)))
tree.list.wildcards <- list()
for (this.wildcard in seq_len(nrow(wildcard.positions))) {
	tree.list.wildcards[[this.wildcard]] <- read.nexus(file = wildcard.positions$Filename[this.wildcard])
	tree.list.wildcards[[this.wildcard]]$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree.list.wildcards[[this.wildcard]]$tip.label, fixed=TRUE)))
}
class(tree.list.wildcards) <- "multiPhylo"

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
	
	tree.list <- replicate(n=reps, expr=tree_resolution_v3(tree.backbone = tree.backbone, clade.definitions = clade.definitions, wildcard.positions = wildcard.positions, tree.list.wildcards=tree.list.wildcards, occs = occs, thisMat = measureMat), simplify=FALSE)
	if(reps>1){class(tree.list) <- "multiPhylo"}

	print("tree.list completed")
	###
	#plot trees for check
	###
	#quartz(width=14, height= 10)
	#	par(mfrow = c(1,reps))
	#jk <- 1
	#while(jk < reps +1) {
	#	this.map <-	paintWildcardClades_Dated(tree.list[[jk]])
	#	jk <- jk+1
	#}
	# plot(tree.list[[1]])
		
 	measureMatCrop <- matrix()
 	dat.vec <- data.frame()
 	#if(reps == 1) { measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[1]], thisMat = measureMat)
 	if(reps == 1) { measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list, thisMat = measureMat)
 	} else { measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[1]], thisMat = measureMat) }
	# measureMat[tree.list[[1]]$tip.label,] # drop taxa from the matrix not in the tree, also reorders matrix to match order in tree
	dat.vec <- as.data.frame(data.matrix(frame = measureMatCrop[,"bodyMass"]))
	dimnames(dat.vec) <- list(rownames(measureMatCrop), "body.mass")
	names(dat.vec) <- "body.mass"
	
	print("Dat vector completed")

}

if (save.files==TRUE){
	tree.dat() <- list()
	write.tree(tree.list,file="DoughtyEvanSec1ALab2_NAUngulateTree_test.tre")
	for(ii in seq(1,reps,1)){
		tree.dat[[ii]] <-tree.list[[ii]]$node.date
		#	tree.dat[[2]]<-tree.list[[1]]$root.length
		#write.csv(tree.dat,file="DoughtyEvanSec1ALab2_NAUngulateTree_test.csv")
		save(tree.dat, tree.list, dat.vec, file="DoughtyEvanSec1ALab2_NAUngulateTree_test.Rdata")}}

#attach $node.date to each tree and recaclucate root.length
if(load.files==TRUE){
	#tree.list<- read.tree("DoughtyEvanSec1ALab2_NAUngulateTree_test.tre")
	load("~/Desktop/c230LAB/Project/DoughtyEvanSec1ALab2_NAUngulateTree_test.Rdata")
	for(ii in seq(1,length(tree.list),1)){
		tree.list[[ii]]$node.date <- tree.dat[[ii]]}}

# reps <- 10
analysis.settings <- list(do.BMsimple = TRUE, do.BMM = TRUE, do.OUsimple = TRUE, do.OUM = TRUE, do.heuristic = FALSE, adjust.date.after=10, adjust.date.increment = 0.1, this.constraint = FALSE, do.parallel = TRUE)

# do.BMsimple = TRUE
# do.BMM = TRUE
# do.OUsimple = TRUE
# do.OUM = TRUE
# do.heuristic = FALSE
# adjust.date.after=10
# adjust.date.increment = 0.1
# this.constraint = FALSE
# do.parallel = FALSE

MaxTime <- 65 # setting this value to 80+ will cause error when running non-hueristic model for the earliest break (i.e. near the base of tree)
MinTime <- 1
int_length <- 2

###
#Need to impliment way to find and remove intervals that do not contain taxa
#Could retool code from ecology projject to determine how many taxa are in the specificed intervals and remove those that lack any taxa from the breaklist
###

intervals <- makeIntervals(MinTime, MaxTime, int_length)
print("Intervals completed")

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
#	print(this.tree)
	results <- testRateShiftsMVMorphIntervals(tree=this.tree, dat = dat.vec, intervals = intervals, analysis.settings = analysis.settings)
	print("Test Rate Shifts Done")
	#need to save or output list that contains all iterations of models#may need to input check for intervals to be included in filename
	# filename_save <- paste("/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModelResultsList_Hueristics=",analysis.settings$do.heuristic,sep="")
	#filename_save <- paste("~/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModelResultsList_Hueristics=", 
	filename_save <- paste("~/Desktop/c230LAB/Project/C230Project_NAUngulateTree_test_Hueristics=", 
		analysis.settings$do.heuristic,
		"_do.parallel=", analysis.settings$do.parallel,
		"_do.BMsimple=", analysis.settings$do.BMsimple,
		"_do.BMM=", analysis.settings$do.BMM, 
		"_do.OUsimple=", analysis.settings$do.OUsimple,
		"_do.OUMp=", analysis.settings$do.OUM,
		timestamp(),".RData")

	save(results, file= filename_save)
}

#sapply(model.finalResults, function(x) x$optList$bm)

# # filename_save <- paste("/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModels_do.Bmsimple=",analysis.settings$do.BMsimple)
	# filename_save <- paste(filename_save, "_do.BMM=", sep="") 
	# filename_save<- paste(filename_save, analysis.settings$do.BMM, sep="")
	# filename_save <- paste(filename_save, "_do.OUsimple=", sep = "")
	# filename_save <- paste(filename_save, analysis.settings$do.OUsimple, sep ="")
	# filename_save <- paste(filename_save, "_do.OUMp=", sep="")
	# filename_save <- paste(filename_save, analysis.settings$do.OUM, sep="")
	# filename_save <- paste(filename_save,"TemoralInt", sep="")
	# filename_save <- paste(filename_save, "_TimeIntLength=", sep="")
	# filename_save <- paste(filename_save,int_length, sep="")	
	# filename_save <- paste(filename_save, "_MaxTime=", sep="")
	# filename_save <- paste(filename_save,MaxTime, sep="")	
	# filename_save <- paste(filename_save, "_MinTime=", sep="")
	# filename_save <- paste(filename_save,MinTime, sep="")	
	# filename_save <- paste(filename_save, "_rep=", sep="")
	# filename_save <- paste(filename_save,reps, sep="")
	# filename_save <- paste(filename_save, "_do.parallel=", sep="")
	# filename_save <- paste(filename_save, analysis.settings$do.parallel, sep="")
	# filename_save <- paste(filename_save, ".RData", sep="")
	
	# save(model.finalResults, modelTime,tree.list, file = filename_save)
	# print("Data Saved; Program Finished")
	
	# load(filename_save)
