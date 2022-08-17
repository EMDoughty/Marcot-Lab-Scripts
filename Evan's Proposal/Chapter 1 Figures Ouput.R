##make the neccessary objects

source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/occFns.R")
source("~/Dropbox/code/R/common_src/sampling.R") 
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R")
source("~/Dropbox/Code/R/common_src/taxonomicEv.R")

source('~/Dropbox/code/R/dentalMeasurements/src/src_dentalDataFns.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_ecologyAnalysisFns.R', chdir = TRUE)


source('~/Dropbox/Code/Recreate Dental Plots Source 2021_6_2.R')

require(plyr)
require(abind)
require(corrplot)

#run.taxon <-  "carnivores" #"ungulates"
this.rank <- "species" #"genus" "species"
interval.type <- "bins" #"nalma" "bins

bmBreaks_herb <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
bmBreaks_pred <- c(-Inf, 0, 0.845098, 1.322219, 2, Inf) 

#do.pred.diet <- TRUE
#diet.type <- c("hypercarnivore") #, "mesocarnivore", "hypocarnivore")

save.pathname <- "~/Dropbox/Proposal/Chapter 1 Figures/"

#if you want to load a repIntOccs or repIntTaxa from file put the pathname as this object.  otherwise keep as NUll to make a new repIntOccs and repIntTaxa using the settings below
  ####Full runs

if(interval.type %in% "bins")  repIntLoad <- "/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/repIntMaster__this.rank=species_timebin=2Mabins_start=64Ma_end=0Ma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Sun Mar  6 19:30:39 2022 ------##.Rdata"

if(interval.type %in% "nalma") repIntLoad <- "/Users/emdoughty/Dropbox/Code/R/Results/NALMA/repIntMaster__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 00:51:19 2022 ------##.Rdata"

load(repIntLoad)

###################################################################################################
occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=100&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", 
                                 "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", 
                                 "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Panotariidae","Phocidae", "Physeteridae", 
                                 "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

if(interval.type == "bins")
{
  int_length <- 2
  intervals <- makeIntervals(0, 64, int_length) # for some reason it adds an additional interval when doing odd midpoints (i.e. 0 to 64 Ma have a 64-66Ma interval while 1-64 Ma will not)....its late i may be crazy
  intList <- listifyMatrixByRow(intervals)
  
  save.path.bins <- paste0(int_length,"Ma",interval.type)
}
#######################

if(interval.type == "nalma")
{
  nalma.mark <- read.csv("/Users/emdoughty/Dropbox/Proposal/NOW_intervals_edit.csv")
  nalma.mark <- nalma.mark[,1:3]
  #nalma.add  <- rbind(c("Aquilian", 84, 70),c("Lancian", 70, 66),c("Puercan", 66, 64.81)); colnames(nalma.add) <- colnames(nalma.mark)
  # nalma.mark <- rbind(nalma.mark, nalma.add)
  #  nalma.mark[,2] <- as.numeric(nalma.mark[,2])
  # nalma.mark[,3] <- as.numeric(nalma.mark[,3])
  # nalma.mark <- nalma.mark[order(as.numeric(nalma.mark$Max_age), decreasing = TRUE),]
  rownames(nalma.mark) <- nalma.mark$CHRON
  colnames(nalma.mark) <- c("NALMA_Subdivision", "ageBase","ageTop")
  nalma.mark <- nalma.mark[,-1]
  nalma.mark <- nalma.mark[,c(2,1)]
  intervals <- nalma.mark[nrow(nalma.mark):1,]
  
  save.path.bins <- paste0(interval.type)
}

intervals.34Ma <- intervals[intervals$ageBase <= 34,]

repIntTaxa_master <- repIntTaxa
repIntTaxa.34Ma <- lapply(repIntTaxa, function(this.rep) {this.rep[rownames(intervals.34Ma)]})

intervals.24Ma <- intervals[intervals$ageBase <= 24,]
repIntTaxa.24Ma <- lapply(repIntTaxa, function(this.rep) {this.rep[rownames(intervals.24Ma)]})

measure.mat <- getMeasureMatWithBodyMasses()

####################################################################################################################################
#### reduces matrix to just the focal order(s)
####################################################################################################################################

archaic.ung <- read.csv("/Users/emdoughty/Dropbox/Code/ArchaicUngulate_UploadFile_2021_4_29.csv")

archaic.Mat.tot <- getMeasureMatCondylarths(data.raw = archaic.ung, occs = occs, 
                                            col.order = colnames(measure.mat), 
                                            all.bm = FALSE,
                                            regression = "ArchaicNonselenodonts")
#test <- jon_archicBodyMassEstim()
#test.bm <- 10^test$bodyMass
#test.bm <- test.bm/1000
#test.bm <- log10(test.bm)

measure.mat <- rbind(measure.mat, archaic.Mat.tot)

focal.order <- c("Artiodactyla", "Perissodactyla")
focal.family <- unique(occs[occs$order %in% focal.order,]$family)
add.family <- c("Arctocyonidae", "Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae")
focal.family <- c(as.character(focal.family), add.family)
focal.family <- focal.family[!focal.family %in% ""]
focal.family <- focal.family[order(focal.family)]

bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & (occs$order %in% focal.order | occs$family %in% add.family)), c("order","family", "genus", "accepted_name")])
#bigList.cond <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% add.family), c("order","family", "genus", "accepted_name")])
# bigList <- rbind(bigList,bigList.cond)

bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
shortFam <- sort(unique(bigList$family[bigList$family %in% focal.family]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name[bigList$family %in% shortFam], ]

#this.rank <- "genus"
if (this.rank=="genus") measure.mat<- makeOneGenusMatFromSpecimenMat(measure.mat)

#########################################################################################################################
pred.data <- read.csv("~/Dropbox/Proposal/predator_data_final.csv")
colnames(pred.data) <- c("family", "taxon",	"max_ma","min_ma",	"m1L",	"rbl",	"bodyMass",	"Citation")
pred.data$taxon <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.data$taxon)
rownames(pred.data) <- pred.data$taxon

pred.data[,c("bodyMass")] <- log10(pred.data[,c("bodyMass")])
pred.data <- pred.data[is.finite(pred.data$bodyMass),]

focal.orderPred <- c("Carnivora", "Creodonta","Hyaenodonta")
#occsPred <- occs[occs$order %in% focal.orderPred,]
focal.familyPred <- unique(occs[occs$order %in% focal.orderPred,]$family)
focal.familyPred <- c(as.character(focal.familyPred), "Viverravidae")

focal.familyPred <- focal.familyPred[!focal.familyPred%in% ""]
focal.familyPred <- focal.familyPred[order(focal.familyPred)]

bigListPred <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% focal.familyPred), c("order","family", "genus", "accepted_name")])
bigListPred <- bigListPred[order(bigListPred$order, bigListPred$family, bigListPred$genus, bigListPred$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]

bigListPred$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigListPred$accepted_name)
bigListPred <- bigListPred[bigListPred$accepted_name %in% pred.data$taxon,]
shortFamPred <- sort(unique(bigListPred$family[bigListPred$family %in% focal.familyPred]))	

#pred.data <- getDateTaxa(measure.mat = pred.data, occs=occsPred, this.rank = "species")
pred.data <- pred.data[complete.cases(pred.data),]

pred.diet <- read.csv("/Users/emdoughty/Dropbox/Proposal/Diet Data PPP CJ.csv")
pred.diet$MASTER_LIST <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.diet$MASTER_LIST)

pred.data.master  <- pred.data

pred.data <- pred.data[pred.data$taxon %in% pred.diet$MASTER_LIST,]
pred.data$Diet <- pred.diet[pred.diet$MASTER_LIST %in% pred.data$taxon, "PPP_diet_2"]

pred.data.hyper <- pred.data[pred.data$Diet %in% "hypercarnivore",]
pred.data.meso <- pred.data[pred.data$Diet %in% "mesocarnivore",]
pred.data.hypo <- pred.data[pred.data$Diet %in% "hypocarnivore",]
pred.data.hyper.meso <- pred.data[pred.data$Diet %in% c("hypercarnivore","mesocarnivore"),]
pred.data.hyper.hypo <- pred.data[pred.data$Diet %in% c("hypercarnivore","hypocarnivore"),]
pred.data.meso.hypo <- pred.data[pred.data$Diet %in% c("mesocarnivore","hypocarnivore"),]

test <- measure.mat[measure.mat$taxon %in% occs$accepted_name[occs$max_ma <= 24 & occs$min_ma <= 24],]
test <- pred.data[pred.data$taxon %in% occs$accepted_name[occs$max_ma <= 24 & occs$min_ma <= 24],]

######################################################################################################################

rem.col <- c("family","genus","reg.vec")
measure.mat <- measure.mat[,!colnames(measure.mat) %in% rem.col]

countCube_herb <- sapply(repIntTaxa, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(measure.mat[,"bodyMass"][match(this.intv, measure.mat$taxon)], 
         breaks= bmBreaks_herb, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countCube <- countCube[,,1]

sizecateg <- c("<5 kg", 
               "5-25 kg",
               "25-150kg", 
               "150-500 kg", 
               ">500 kg")

dimnames(countCube_herb) <- list(sizecateg, rownames(intervals), NULL)


pred.data <- pred.data[,!colnames(pred.data) %in% rem.col]

countCube_pred <- sapply(repIntTaxa, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(pred.data[,"bodyMass"][match(this.intv, pred.data$taxon)], 
         breaks= bmBreaks_pred, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countCube <- countCube[,,1]

sizecateg <- c("<1kg",
               "1-7kg",
               "7-21kg", 
               "21-100kg",
               ">100kg")

dimnames(countCube_pred) <- list(sizecateg, rownames(intervals), NULL)


pred.data.hyper <- pred.data.hyper[,!colnames(pred.data.hyper) %in% rem.col]

countCube_pred_hyper <- sapply(repIntTaxa.34Ma, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(pred.data.hyper[,"bodyMass"][match(this.intv, pred.data.hyper$taxon)], 
         breaks= bmBreaks_pred, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countCube <- countCube[,,1]

dimnames(countCube_pred_hyper) <- list(sizecateg, rownames(intervals.34Ma), NULL)


pred.data.meso <- pred.data.meso[,!colnames(pred.data.meso) %in% rem.col]

countCube_pred_meso <- sapply(repIntTaxa.34Ma, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(pred.data.meso[,"bodyMass"][match(this.intv, pred.data.meso$taxon)], 
         breaks= bmBreaks_pred, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countCube <- countCube[,,1]

dimnames(countCube_pred_meso) <- list(sizecateg, rownames(intervals.34Ma), NULL)


pred.data.hypo <- pred.data.hypo[,!colnames(pred.data.hypo) %in% rem.col]

countCube_pred_hypo <- sapply(repIntTaxa.34Ma, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(pred.data.hypo[,"bodyMass"][match(this.intv, pred.data.hypo$taxon)], 
         breaks= bmBreaks_pred, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countCube <- countCube[,,1]

dimnames(countCube_pred_hypo) <- list(sizecateg, rownames(intervals.34Ma), NULL)

pred.data.hyper.meso <- pred.data.hyper.meso[,!colnames(pred.data.hyper.meso) %in% rem.col]

countCube_pred_hyper.meso <- sapply(repIntTaxa.34Ma, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(pred.data.hyper.meso[,"bodyMass"][match(this.intv, pred.data.hyper.meso$taxon)], 
         breaks= bmBreaks_pred, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countCube <- countCube[,,1]

dimnames(countCube_pred_hyper.meso) <- list(sizecateg, rownames(intervals.34Ma), NULL)

pred.data.hyper.hypo<- pred.data.hyper.hypo[,!colnames(pred.data.hyper.hypo) %in% rem.col]

countCube_pred_hyper.hypo <- sapply(repIntTaxa.34Ma, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(pred.data.hyper.hypo[,"bodyMass"][match(this.intv, pred.data.hyper.hypo$taxon)], 
         breaks= bmBreaks_pred, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countCube <- countCube[,,1]

dimnames(countCube_pred_hyper.hypo) <- list(sizecateg, rownames(intervals.34Ma), NULL)

pred.data.meso.hypo<- pred.data.meso.hypo[,!colnames(pred.data.meso.hypo) %in% rem.col]

countCube_pred_meso.hypo <- sapply(repIntTaxa.34Ma, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(pred.data.meso.hypo[,"bodyMass"][match(this.intv, pred.data.meso.hypo$taxon)], 
         breaks= bmBreaks_pred, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countCube <- countCube[,,1]

dimnames(countCube_pred_meso.hypo) <- list(sizecateg, rownames(intervals.34Ma), NULL)
#####################################################################################################################
if(interval.type %in% "bins") load("/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/BM_handleyResult_ungulates_this.rank=species_timebin=2Mabins_start=66_end=0_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Thu Mar 10 19:47:14 2022 ------##.Rdata")
if(interval.type %in% "nalma") load("/Users/emdoughty/Dropbox/Code/R/Results/NALMA/BM_handleyResult_ungulates_this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 03:36:19 2022 ------##.Rdata")
#countCube_ung <- countCube
optList_bm_median_ung <- optList_bm_median
optList_bm_allReps_ung <- optList_bm_allReps

bm_quants_ung <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[as.character(x),"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)

if(interval.type %in% "bins") load("/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/BM_handleyResult_carnivores_allCateg_this.rank=species_timebin=2Mabins_start=66_end=0_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Thu Mar 10 21:44:05 2022 ------##.Rdata")
if(interval.type %in% "nalma") load("/Users/emdoughty/Dropbox/Code/R/Results/NALMA/BM_handleyResult_carnivores_allCateg_this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 15:58:31 2022 ------##.Rdata")
optList_bm_median_pred <- optList_bm_median
optList_bm_allReps_pred <- optList_bm_allReps

bm_quants_pred <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(pred.data[as.character(x),"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)

if(interval.type %in% "bins") load("/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/BM_handleyResult_carnivores_allCateg_hypercarnivores_this.rank=species_timebin=2Mabins_start=34_end=0_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Fri Mar 11 12:00:42 2022 ------##.Rdata")
if(interval.type %in% "nalma") load("/Users/emdoughty/Dropbox/Code/R/Results/NALMA/BM_handleyResult_carnivores_allCateg_hyper__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 17:06:41 2022 ------##.Rdata")
optList_bm_median_pred_hyper <- optList_bm_median
optList_bm_allReps_pred_hyper <- optList_bm_allReps

bm_quants_pred_hyper <- apply(sapply(repIntTaxa.34Ma, function(y) sapply(y, function(x) quantile(pred.data.hyper[as.character(x),"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)

if(interval.type %in% "bins") load("/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/BM_handleyResult_carnivores_allCateg_mesocarnivore__this.rank=species_timebin=2Mabins_start=34_end=0_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Fri Mar 11 14:17:49 2022 ------##.Rdata")
if(interval.type %in% "nalma") load("/Users/emdoughty/Dropbox/Code/R/Results/NALMA/BM_handleyResult_carnivores_allCateg_meso__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 19:47:43 2022 ------##.Rdata")
optList_bm_median_pred_meso <- optList_bm_median
optList_bm_allReps_pred_meso <- optList_bm_allReps

bm_quants_pred_meso <- apply(sapply(repIntTaxa.34Ma, function(y) sapply(y, function(x) quantile(pred.data.meso[as.character(x),"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)

if(interval.type %in% "bins") load("/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/BM_handleyResult_carnivores_allCateg_hypocarnivore__this.rank=species_timebin=2Mabins_start=34_end=0_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Fri Mar 11 17:18:18 2022 ------##.Rdata")
if(interval.type %in% "nalma") load("/Users/emdoughty/Dropbox/Code/R/Results/NALMA/BM_handleyResult_carnivores_allCateg_hypo__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Wed Mar 16 12:29:15 2022 ------##.Rdata")
optList_bm_median_pred_hypo <- optList_bm_median
optList_bm_allReps_pred_hypo <- optList_bm_allReps

bm_quants_pred_hypo <- apply(sapply(repIntTaxa.34Ma, function(y) sapply(y, function(x) quantile(pred.data.hypo[as.character(x),"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)

if(interval.type %in% "bins") load("/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/BM_handleyResult_carnivores_allCateg_hyper_meso__this.rank=species_timebin=2Mabins_start=66_end=0_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Mon Mar 21 18:02:46 2022 ------##.Rdata")
if(interval.type %in% "nalma") load("/Users/emdoughty/Dropbox/Code/R/Results/NALMA/BM_handleyResult_carnivores_allCateg_hyper_meso__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 21:53:56 2022 ------##.Rdata")
optList_bm_median_pred_hyper.meso <- optList_bm_median
optList_bm_allReps_pred_hypo.meso <- optList_bm_allReps

bm_quants_pred_hyper.meso <- apply(sapply(repIntTaxa.34Ma, function(y) sapply(y, function(x) quantile(pred.data.hyper.meso[as.character(x),"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)

if(interval.type %in% "bins") load("/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/BM_handleyResult_carnivores_allCateg_hyper_hypo__this.rank=species_timebin=2Mabins_start=66_end=0_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Mon Mar 21 16:24:01 2022 ------##.Rdata")
if(interval.type %in% "nalma") load("/Users/emdoughty/Dropbox/Code/R/Results/NALMA/BM_handleyResult_carnivores_allCateg_hyper_hypo__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 22:55:31 2022 ------##.Rdata")
optList_bm_median_pred_hyper.hypo <- optList_bm_median
optList_bm_allReps_pred_hyper.hypo <- optList_bm_allReps

bm_quants_pred_hyper.hypo <- apply(sapply(repIntTaxa.34Ma, function(y) sapply(y, function(x) quantile(pred.data.hyper.hypo[as.character(x),"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)

if(interval.type %in% "bins") load("/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/BM_handleyResult_carnivores_allCateg_meso_hypo__this.rank=species_timebin=2Mabins_start=66_end=0_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Wed Mar 16 16:40:06 2022 ------##.Rdata")
if(interval.type %in% "nalma") load("/Users/emdoughty/Dropbox/Code/R/Results/NALMA/BM_handleyResult_carnivores_allCateg_meso_hypo__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 23:50:54 2022 ------##.Rdata")
optList_bm_median_pred_meso.hypo <- optList_bm_median
optList_bm_allReps_pred_meso.hypo <- optList_bm_allReps

bm_quants_pred_meso.hypo <- apply(sapply(repIntTaxa.34Ma, function(y) sapply(y, function(x) quantile(pred.data.meso.hypo[as.character(x),"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)


###############################################################

countBox_herb <- apply(countCube_herb, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

countBox_pred <- apply(countCube_pred, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

countBox_pred_hyper <- apply(countCube_pred_hyper, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

countBox_pred_meso <- apply(countCube_pred_meso, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

countBox_pred_hypo <- apply(countCube_pred_hypo, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

countBox_pred_hyper.meso <- apply(countCube_pred_hyper.meso, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

countBox_pred_hyper.hypo <- apply(countCube_pred_hyper.hypo, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

countBox_pred_meso.hypo <- apply(countCube_pred_meso.hypo, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 
##############################################################
prop_herb <- t(apply(countCube_herb, c(1,2), median, na.rm=TRUE))
colnames(prop_herb)[colnames(prop_herb)==""] <- "indeterminate"

prop_pred <- t(apply(countCube_pred, c(1,2), median, na.rm=TRUE))
colnames(prop_pred)[colnames(prop_pred)==""] <- "indeterminate"

prop_pred_7kg <- prop_pred
prop_pred_7kg[,c(1,2)] <- 0

prop_pred_hyper <- t(apply(countCube_pred_hyper, c(1,2), median, na.rm=TRUE))
colnames(prop_pred_hyper)[colnames(prop_pred_hyper)==""] <- "indeterminate"

prop_pred_meso <- t(apply(countCube_pred_meso, c(1,2), median, na.rm=TRUE))
colnames(prop_pred_meso)[colnames(prop_pred_meso)==""] <- "indeterminate"

prop_pred_hypo <- t(apply(countCube_pred_hypo, c(1,2), median, na.rm=TRUE))
colnames(prop_pred_hypo)[colnames(prop_pred_hypo)==""] <- "indeterminate"

prop_pred_hyper.meso <- t(apply(countCube_pred_hyper.meso, c(1,2), median, na.rm=TRUE))
colnames(prop_pred_hyper.meso)[colnames(prop_pred_hyper.meso)==""] <- "indeterminate"

prop_pred_hyper.hypo <- t(apply(countCube_pred_hyper.hypo, c(1,2), median, na.rm=TRUE))
colnames(prop_pred_hyper.hypo)[colnames(prop_pred_hyper.hypo)==""] <- "indeterminate"

prop_pred_meso.hypo <- t(apply(countCube_pred_meso.hypo, c(1,2), median, na.rm=TRUE))
colnames(prop_pred_meso.hypo)[colnames(prop_pred_meso.hypo)==""] <- "indeterminate"

###############################################################



###############################################################
#quartz(width = 8, height =8)
if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Species Richness Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
     width = 8, height = 9, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Species Richness Plots Draft NALMA 2022_3_18.png",
      width = 8, height = 9, units = "in", res = 200)
}
  
par(mfrow = c(4,1), mar=c(2,4,0,2), mgp=c(2, 1,0), oma = c(2,2,2,0))
plotStackedRichness(this.box=prop_herb, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    xlim = c(65, min(intervals, na.rm=TRUE)), ylim = c(0, max(rowSums(prop_herb))+20), xaxp = c(65, 0, 13), yaxp = c(0, 80, 4),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Large Herbivores", side = 2, line = 3.5, cex = 1.25, at = -10)

plotStackedRichness(this.box=prop_herb/rowSums(prop_herb), intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    xlim = c(65, min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

plotStackedRichness(this.box=prop_pred, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    xlim = c(65,  min(intervals, na.rm=TRUE)), ylim = c(0, max(rowSums(prop_pred)+10)), xaxp = c(65, 0, 13), yaxp = c(0, 30, 6),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Predators", side = 2, line = 3.5, cex = 1.25, at = -10)

plotStackedRichness(this.box=prop_pred/rowSums(prop_pred), intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    xlim = c(65,  min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)
mtext(text = "Time (Ma)", side = 1, line = 2, cex = 0.75)
dev.off()


##### >7 kg
if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Species Richness Plots pred>7kg Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 8, height = 5, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Species Richness Plots pred>7kg  Draft NALMA 2022_3_18.png",
      width = 8, height = 5, units = "in", res = 200)
}

par(mfrow = c(2,1), mar=c(2,4,0,2), mgp=c(2, 1,0), oma = c(2,2,2,0))
plotStackedRichness(this.box=prop_herb, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    xlim = c(65, min(intervals, na.rm=TRUE)), ylim = c(0, max(rowSums(prop_herb))+20), xaxp = c(65, 0, 13), yaxp = c(0, 80, 4),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Large Herbivores", side = 2, line = 3.5, cex = 1.25, at = -10)

plotStackedRichness(this.box=prop_pred_7kg, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    xlim = c(65,  min(intervals, na.rm=TRUE)), ylim = c(0, max(rowSums(prop_pred)+10)), xaxp = c(65, 0, 13), yaxp = c(0, 30, 6),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Predators", side = 2, line = 3.5, cex = 1.25, at = -10)

dev.off()


# PREDATOR DIETARY CATEGORIES
if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Species Richness Plots for Pred Diet Draft 2 Ma Intervals 64 to 0 Ma 2022_4_8.png",
      width = 15, height = 9, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Species Richness Plots for Pred Diet Draft NALMA 2022_4_8.png",
      width = 15, height = 9, units = "in", res = 200)
}

par(mfrow = c(3,3), mar=c(1,4,1.6,0), mgp=c(2.5, 1,0), oma = c(4,2,2,2))
plotStackedRichness(this.box=prop_pred_hyper, intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 35), xaxp = c(65, 0, 13), yaxp = c(0, 30, 6), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Hypercarnivores", side = 2, line = 3.75, cex = 1.25)
mtext(text = "Median Species Richness", side = 3, line = 2, cex = 1.25)

plotStackedRichness(this.box=prop_pred_hyper/rowSums(prop_pred_hyper), intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Proportional Species Richness", side = 3, line = 2, cex = 1.25)
mtext(text = "Within Diet Category", side = 3, line = 0.70, cex = 1.25)


plotStackedRichness(this.box=prop_pred_hyper/rowSums(prop_pred[rownames(prop_pred) %in% rownames(intervals.34Ma),]), intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE, 
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Proportional Species Richness", side = 3, line = 2, cex = 1.25)
mtext(text ="Relative to All Predators ", side = 3, line = 0.70, cex = 1.25)

####
plotStackedRichness(this.box=prop_pred_meso, intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 35), xaxp = c(65, 0, 13), yaxp = c(0, 30, 6), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Mesocarnivores", side = 2, line = 3.75, cex = 1.25)

plotStackedRichness(this.box=prop_pred_meso/rowSums(prop_pred_meso), intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

plotStackedRichness(this.box=prop_pred_meso/rowSums(prop_pred[rownames(prop_pred) %in% rownames(intervals.34Ma),]), intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

####
plotStackedRichness(this.box=prop_pred_hypo, intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 35), xaxp = c(65, 0, 13), yaxp = c(0, 30, 6), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)

mtext(text = "Hypocarnivores", side = 2, line = 3.75, cex = 1.25)
mtext(text = "Time (Ma)", side = 1, line = 3, cex = 1)

plotStackedRichness(this.box=prop_pred_hypo/rowSums(prop_pred_hypo), intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)
mtext(text = "Time (Ma)", side = 1, line = 3, cex = 1)

plotStackedRichness(this.box=prop_pred_hypo/rowSums(prop_pred[rownames(prop_pred) %in% rownames(intervals.34Ma),]), intervals=intervals.34Ma, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    xlim = c(34,  min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5), cex.axis = 1.5, cex.lab = 1.5,
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)
mtext(text = "Time (Ma)", side = 1, line = 3, cex = 1)

dev.off()

###############################################################
if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Large Herbivore Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_4_25.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley nalysis Large Herbivore Plots Draft NALMA 2022_4_25.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1)) #mfrow = c(8,1),
shoulderPlot(measure.mat = measure.mat, plot.y = "bodyMass", intervals = intervals, occs = occs, this.rank = "species", 
             bigList = bigList, shortFam = shortFam, repIntTaxa = repIntTaxa_ung, quants = bm_quants_ung,
             optList_bm_median = optList_bm_median_ung,  
             specOcc.col = "gray0", specOcc.alpha = 0.5,
             plot.breaks = TRUE,  break.text = FALSE, break.col = "blue", #manual.breaks = c(62, 50, 46, 40, 34, 22, 4),
             ylab = "log bodymass (kg)", xlab = "Time (Ma)", xaxp = c(70,0,14), ylim = c(-0.5,4.5), yaxp = c(0,4,4),
             cex.axis = 1, cex.lab = 1, 
             overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.90, text.offset = 0.05,
             do.quants = TRUE,
             col.axis = "black", col.lab = "black", poly.col = "deepskyblue4", median.col = c("darkturquoise", "deepskyblue4", "deepskyblue1"))

dev.off()

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Predator Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_4_25.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley nalysis Predator Plots Draft NALMA 2022_4_25.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1)) #mfrow = c(8,1),
shoulderPlot(measure.mat = pred.data, plot.y = "bodyMass", intervals = intervals, occs = occs, this.rank = "species", 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = repIntTaxa_pred, quants = bm_quants_pred,
             optList_bm_median = optList_bm_median_pred,  #ylim = c(-0.5,  4),
             specOcc.col = "gray0", specOcc.alpha = 0.5,
             plot.breaks = TRUE,  break.text = FALSE, break.col = "red", #manual.breaks = c(57, 41, 25),
             ylab = "log bodymass (kg)", xlab = "Time (Ma)", xaxp = c(70,0,14), ylim = c(-1.5,4.5), yaxp = c(-1,3,4),
             cex.axis = 1, cex.lab = 1, 
             overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.90, text.offset = 0.05,
             do.quants = TRUE,
             col.axis = "black", col.lab = "black", poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))

dev.off()


if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Hypercarnivores Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_4_25.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Hypercarnivores Predator Plots Draft NALMA 2022_4_25.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1)) #mfrow = c(8,1),
shoulderPlot(measure.mat = pred.data.hyper, plot.y = "bodyMass", intervals = intervals.34Ma, occs = occs, this.rank = "species", 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = repIntTaxa, quants = bm_quants_pred_hyper,
             optList_bm_median = optList_bm_median_pred_hyper,  xlim = c(65, 0),
             specOcc.col = "gray0", specOcc.alpha = 0.5,
             plot.breaks = TRUE,  break.text = FALSE, break.col = "firebrick4", #manual.breaks = c(57, 41, 25),
             ylab = "log bodymass (kg)", xlab = "", xaxp = c(70,0,14), ylim = c(-1.5,4.5), yaxp = c(-1,3,4),
             cex.axis = 1, cex.lab = 1, 
             overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.90, text.offset = 0.05,
             do.quants = TRUE,
             col.axis = "black", col.lab = "black", poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))
mtext(text = "Hypercarnivores", side = 2, line = 3, cex = 1.25)
mtext(text = "Time (Ma)", side = 1, line = 1.75, cex = 1, at = 15)

dev.off()



if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Mesocarnivores Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Mesocarnivores Predator Plots Draft NALMA 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1)) #mfrow = c(8,1),
shoulderPlot(measure.mat = pred.data.meso, plot.y = "bodyMass", intervals = intervals.34Ma, occs = occs, this.rank = "species", 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = repIntTaxa, quants = bm_quants_pred_meso,
             optList_bm_median = optList_bm_median_pred_meso,  xlim = c(65, 0),
             specOcc.col = "gray0", specOcc.alpha = 0.5,
             plot.breaks = TRUE,  break.text = FALSE, break.col = "firebrick4", #manual.breaks = c(57, 41, 25),
             ylab = "log bodymass (kg)", xlab = "", xaxp = c(70,0,14), ylim = c(-1.5,4.5), yaxp = c(-1,3,4),
             cex.axis = 1, cex.lab = 1, 
             overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.90, text.offset = 0.05,
             do.quants = TRUE,
             col.axis = "black", col.lab = "black", poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))
mtext(text = "Mesocarnivores", side = 2, line = 3, cex = 1.25)
mtext(text = "Time (Ma)", side = 1, line = 1.75, cex = 1, at = 15)


dev.off()


if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Hypocarnivores Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Hypocarnivores Predator Plots Draft NALMA 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1)) #mfrow = c(8,1),
shoulderPlot(measure.mat = pred.data.hypo, plot.y = "bodyMass", intervals = intervals.34Ma, occs = occs, this.rank = "species", 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = repIntTaxa, quants = bm_quants_pred_hypo,
             optList_bm_median = optList_bm_median_pred_hypo, xlim = c(65, 0),
             specOcc.col = "gray0", specOcc.alpha = 0.5,
             plot.breaks = TRUE,  break.text = FALSE, break.col = "firebrick4", #manual.breaks = c(57, 41, 25),
             ylab = "log bodymass (kg)", xlab = "", xaxp = c(70,0,14), ylim = c(-1.5,4.5), yaxp = c(-1,3,4),
             cex.axis = 1, cex.lab = 1, 
             overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.90, text.offset = 0.05,
             do.quants = TRUE,
             col.axis = "black", col.lab = "black", poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))
mtext(text = "Hypocarnivores", side = 2, line = 3, cex = 1.25)
mtext(text = "Time (Ma)", side = 1, line = 1.75, cex = 1, at = 15)

dev.off()


if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Hyper_Meso Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Hyper_Meso Predator Plots Draft NALMA 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1)) #mfrow = c(8,1),
shoulderPlot(measure.mat = pred.data.hyper.meso, plot.y = "bodyMass", intervals = intervals.34Ma, occs = occs, this.rank = "species", 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = repIntTaxa, quants = bm_quants_pred_hyper.meso,
             optList_bm_median = optList_bm_median_pred_hyper.meso, xlim = c(65, 0),
             specOcc.col = "gray0", specOcc.alpha = 0.5,
             plot.breaks = TRUE,  break.text = FALSE, break.col = "firebrick4", #manual.breaks = c(57, 41, 25),
             ylab = "log bodymass (kg)", xlab = "", xaxp = c(70,0,14), ylim = c(-1.5,4.5), yaxp = c(-1,3,4),
             cex.axis = 1, cex.lab = 1, 
             overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.90, text.offset = 0.05,
             do.quants = TRUE,
             col.axis = "black", col.lab = "black", poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))
mtext(text = "Hyper- + Meso-", side = 2, line = 3, cex = 1.25)
mtext(text = "Time (Ma)", side = 1, line = 1.75, cex = 1, at = 15)

dev.off()


if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Hyper_Hypo Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Hyper_Hypo Predator Plots Draft NALMA 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1)) #mfrow = c(8,1),
shoulderPlot(measure.mat = pred.data.hyper.hypo, plot.y = "bodyMass", intervals = intervals.34Ma, occs = occs, this.rank = "species", 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = repIntTaxa, quants = bm_quants_pred_hyper.hypo,
             optList_bm_median = optList_bm_median_pred_hyper.hypo, xlim = c(65, 0),
             specOcc.col = "gray0", specOcc.alpha = 0.5,
             plot.breaks = TRUE,  break.text = FALSE, break.col = "firebrick4", #manual.breaks = c(57, 41, 25),
             ylab = "log bodymass (kg)", xlab = "", xaxp = c(70,0,14), ylim = c(-1.5,4.5), yaxp = c(-1,3,4),
             cex.axis = 1, cex.lab = 1, 
             overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.90, text.offset = 0.05,
             do.quants = TRUE,
             col.axis = "black", col.lab = "black", poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))
mtext(text = "Hyper- + Hypo-", side = 2, line = 3, cex = 1.25)
mtext(text = "Time (Ma)", side = 1, line = 1.75, cex = 1, at = 15)

dev.off()


if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Meso_Hypo Plots Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Handley Analysis Meso_Hypo Predator Plots Draft NALMA 2022_3_18.png",
      width = 8, height = 2.75, units = "in", res = 200)
}

par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1)) #mfrow = c(8,1),
shoulderPlot(measure.mat = pred.data.meso.hypo, plot.y = "bodyMass", intervals = intervals.34Ma, occs = occs, this.rank = "species", 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = repIntTaxa, quants = bm_quants_pred_meso.hypo,
             optList_bm_median = optList_bm_median_pred_meso.hypo, xlim = c(65, 0),
             specOcc.col = "gray0", specOcc.alpha = 0.5,
             plot.breaks = TRUE,  break.text = FALSE, break.col = "firebrick4", #manual.breaks = c(57, 41, 25),
             ylab = "log bodymass (kg)", xlab = "", xaxp = c(70,0,14), ylim = c(-1.5,4.5), yaxp = c(-1,3,4),
             cex.axis = 1, cex.lab = 1, 
             overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.90, text.offset = 0.05,
             do.quants = TRUE,
             col.axis = "black", col.lab = "black", poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))
mtext(text = "Meso- + Hypo-", side = 2, line = 3, cex = 1.25)
mtext(text = "Time (Ma)", side = 1, line = 1.75, cex = 1, at = 15)

dev.off()
###########################################################################################################################

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Large Herbivore Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Large Herbivores Draft NALMA 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

regimeHist_countBox(countBox = countBox_herb[2,,], 
                    breaks = intervals[optList_bm_median_ung[[length(optList_bm_median_ung)-1]]$optBreaks,2], 
                    bmbreaks <- bmBreaks_herb,
                    intervals = intervals,
                    optList = optList_bm_median_ung, 
                    regime.func = median,
                    plot.type = c("proportion"), 
                    ylim = c(0, 1),
                    grayscale = FALSE,
                    mtext.cex = 0.75, mtext.line = 1, cex.names = 1, cex.axis = 1,
                    mfrow = c(1, 8), mar = c(5.5,3,1,0), omi = c(0,0,0.25,0.25))
dev.off()

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist All predator Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist All predator Draft NALMA 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

regimeHist_countBox(countBox = countBox_pred[2,,], 
                    breaks = intervals[optList_bm_median_pred[[length(optList_bm_median_pred)-1]]$optBreaks,2], 
                    bmbreaks <- bmBreaks_pred,
                    intervals = intervals,
                    optList = optList_bm_median_pred, 
                    regime.func = median,
                    plot.type = c("proportion"), 
                    ylim = c(0, 1),
                    grayscale = FALSE,
                    mtext.cex = 0.75, mtext.line = 1, cex.names = 1, cex.axis = 1,
                    mfrow = c(1, 8), mar = c(5.5,3,1,0), omi = c(0,0,0.25,0.25))

dev.off()

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Hypercarnivores Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Hypercarnivores Draft NALMA 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

regimeHist_countBox(countBox = countBox_pred_hyper[2,,], 
                    breaks = intervals[optList_bm_median_pred_hyper[[length(optList_bm_median_pred_hyper)-1]]$optBreaks,2], 
                    bmbreaks <- bmBreaks_pred,
                    intervals = intervals.34Ma,
                    optList = optList_bm_median_pred_hyper, 
                    regime.func = median,
                    plot.type = c("proportion"), 
                    ylim = c(0, 1),
                    grayscale = FALSE,
                    mtext.cex = 0.75, mtext.line = 1, cex.names = 1, cex.axis = 1,
                    mfrow = c(1, 8), mar = c(5.5,3,1,0), omi = c(0,0,0.25,0.25))

dev.off()

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Mesocarnivores Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Mesocarnivores Draft NALMA 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

regimeHist_countBox(countBox = countBox_pred_meso[2,,], 
                    breaks = intervals[optList_bm_median_pred_meso[[length(optList_bm_median_pred_meso)-1]]$optBreaks,2], 
                    bmbreaks <- bmBreaks_pred,
                    intervals = intervals.34Ma,
                    optList = optList_bm_median_pred_meso, 
                    regime.func = median,
                    plot.type = c("proportion"), 
                    ylim = c(0, 1),
                    grayscale = FALSE,
                    mtext.cex = 0.75, mtext.line = 1, cex.names = 1, cex.axis = 1,
                    mfrow = c(1, 8), mar = c(5.5,3,1,0), omi = c(0,0,0.25,0.25))

dev.off()

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Hypocarnivores Draft 2 Ma Intervals 64 to 0 Ma 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Hypocarnivores Draft NALMA 2022_3_18.png",
      width = 9, height = 3, units = "in", res = 200)
}

regimeHist_countBox(countBox = countBox_pred_hypo[2,,], 
                    breaks = intervals[optList_bm_median_pred_hypo[[length(optList_bm_median_pred_hypo)-1]]$optBreaks,2], 
                    bmbreaks <- bmBreaks_pred,
                    intervals = intervals.34Ma,
                    optList = optList_bm_median_pred_hypo, 
                    regime.func = median,
                    plot.type = c("proportion"), 
                    ylim = c(0, 1),
                    grayscale = FALSE,
                    mtext.cex = 0.75, mtext.line = 1, cex.names = 1, cex.axis = 1,
                    mfrow = c(1, 8), mar = c(5.5,3,1,0), omi = c(0,0,0.25,0.25))

dev.off()

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Hyper_Meso Draft 2 Ma Intervals 64 to 0 Ma 2022_3_21.png",
      width = 9, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Hyper_Meso Draft NALMA 2022_3_21.png",
      width = 9, height = 3, units = "in", res = 200)
}

regimeHist_countBox(countBox = countBox_pred_hyper.meso[2,,], 
                    breaks = intervals[optList_bm_median_pred_hyper.meso[[length(optList_bm_median_pred_hyper.meso)-1]]$optBreaks,2], 
                    bmbreaks <- bmBreaks_pred,
                    intervals = intervals.34Ma,
                    optList = optList_bm_median_pred_hyper.meso, 
                    regime.func = median,
                    plot.type = c("proportion"), 
                    ylim = c(0, 1),
                    grayscale = FALSE,
                    mtext.cex = 0.75, mtext.line = 1, cex.names = 1, cex.axis = 1,
                    mfrow = c(1, 8), mar = c(5.5,3,1,0), omi = c(0,0,0.25,0.25))

dev.off()


if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Hyper_Hypo Draft 2 Ma Intervals 64 to 0 Ma 2022_3_21.png",
      width = 9, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Hyper_Hypo Draft NALMA 2022_3_21.png",
      width = 9, height = 3, units = "in", res = 200)
}

regimeHist_countBox(countBox = countBox_pred_hyper.hypo[2,,], 
                    breaks = intervals[optList_bm_median_pred_hyper.hypo[[length(optList_bm_median_pred_hyper.hypo)-1]]$optBreaks,2], 
                    bmbreaks <- bmBreaks_pred,
                    intervals = intervals.34Ma,
                    optList = optList_bm_median_pred_hyper.hypo, 
                    regime.func = median,
                    plot.type = c("proportion"), 
                    ylim = c(0, 1),
                    grayscale = FALSE,
                    mtext.cex = 0.75, mtext.line = 1, cex.names = 1, cex.axis = 1,
                    mfrow = c(1, 8), mar = c(5.5,3,1,0), omi = c(0,0,0.25,0.25))

dev.off()

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Meso_Hypo Draft 2 Ma Intervals 64 to 0 Ma 2022_3_21.png",
      width = 9, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Regime BMdist Meso_Hypo Draft NALMA 2022_3_21.png",
      width = 9, height = 3, units = "in", res = 200)
}

regimeHist_countBox(countBox = countBox_pred_meso.hypo[2,,], 
                    breaks = intervals[optList_bm_median_pred_meso.hypo[[length(optList_bm_median_pred_meso.hypo)-1]]$optBreaks,2], 
                    bmbreaks <- bmBreaks_pred,
                    intervals = intervals.34Ma,
                    optList = optList_bm_median_pred_meso.hypo, 
                    regime.func = median,
                    plot.type = c("proportion"), 
                    ylim = c(0, 1),
                    grayscale = FALSE,
                    mtext.cex = 0.75, mtext.line = 1, cex.names = 1, cex.axis = 1,
                    mfrow = c(1, 8), mar = c(5.5,3,1,0), omi = c(0,0,0.25,0.25))

dev.off()

#########################################################
# Correlations
##All predator
if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices All predator 2 Ma Intervals 64 to 0 Ma 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices All predator Draft NALMA 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

### raw
 par(mfrow = c(1,2))
  ungulates <- prop_herb
  predators <- prop_pred
  
  all.val <- cbind(predators, ungulates)
  
  corr.results.both <- cor(predators, ungulates, method = "spearman")
  cor.p <- cor.mtest(all.val)$p
  
  cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
  cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]
  
  corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                     cl.align.text = 'l', 
                     #addCoef.col = 'black',
                     addgrid.col = 'black',
                     method = "color",
                     na.label = "-",
                     sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                     pch.cex = 1,
                     tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
 # text(p1$corrPos$x, p1$corrPos$y+0.1, round(p1$corr, 2))
  #text(c(rep(1,1)), c(5:1), unique(p1$corrPos$yName)) 
 
#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.55, adj = 1.5)
 mtext("All Predator", side = 2, adj = 0.25, line = 3, cex = 1.5)
  
  
  pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred), intervals = intervals)
  
  UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb), intervals = intervals)
  
  #UngualteBMGroupDiv_FirstDiff
  
  #pred.prey.DivFirstDiff
  
  #should make list of plots of each catagpry through time to make sure data is stationary
  
  ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
  predators <- t(as.matrix(pred.prey.DivFirstDiff))
  
  all.val <- cbind(predators, ungulates)
  
  corr.results.both <- cor(predators, ungulates, method = "spearman")
  cor.p <- cor.mtest(all.val)$p
  
  cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
  cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]
  
  corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                     cl.align.text = 'l', 
                     #addCoef.col = 'black',
                     addgrid.col = 'black',
                     method = "color",
                     na.label = "-",
                     sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                     pch.cex = 1,
                     tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
  
#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.5)
 # mtext("All Predator", side = 2, adj = 0.6, line = 0, cex = 1.5)
dev.off()

#######################################
## Hypercarnivores

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Hypercarnivore 2 Ma Intervals 64 to 0 Ma 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Hypercarnivore Draft NALMA 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

### raw
par(mfrow = c(1,2))
prop_herb.34Ma <- prop_herb[rownames(prop_herb) %in% rownames(intervals.34Ma),]
ungulates <- prop_herb.34Ma
predators <- prop_pred_hyper

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black',
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
# text(p1$corrPos$x, p1$corrPos$y+0.1, round(p1$corr, 2))
#text(c(rep(1,1)), c(5:1), unique(p1$corrPos$yName)) 

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.55, adj = 1.5)
mtext("Hypercarnivores", side = 2, adj = 0.1, line = 3, cex = 1.5)

# First Differences
pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred_hyper), intervals = intervals.34Ma)

UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb.34Ma), intervals = intervals.34Ma)

#UngualteBMGroupDiv_FirstDiff

#pred.prey.DivFirstDiff

#should make list of plots of each catagpry through time to make sure data is stationary

ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
predators <- t(as.matrix(pred.prey.DivFirstDiff))

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black',
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.5)
#  mtext("Predator", side = 2, adj = 0.6, line = 0, cex = 1.5)
dev.off()


#######################################
## Mesocarnivores

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Mesocarnivore 2 Ma Intervals 64 to 0 Ma 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Mesocarnivore Draft NALMA 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

### raw
par(mfrow = c(1,2))
prop_herb.34Ma <- prop_herb[rownames(prop_herb) %in% rownames(intervals.34Ma),]
ungulates <- prop_herb.34Ma
predators <- prop_pred_meso

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black',
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
# text(p1$corrPos$x, p1$corrPos$y+0.1, round(p1$corr, 2))
#text(c(rep(1,1)), c(5:1), unique(p1$corrPos$yName)) 

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.55, adj = 1.5)
mtext("Mesocarnivores", side = 2, adj = 0.1, line = 3, cex = 1.5)

# First Differences
pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred_meso), intervals = intervals.34Ma)

UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb.34Ma), intervals = intervals.34Ma)

#UngualteBMGroupDiv_FirstDiff

#pred.prey.DivFirstDiff

#should make list of plots of each catagpry through time to make sure data is stationary

ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
predators <- t(as.matrix(pred.prey.DivFirstDiff))

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black',
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.5)
#  mtext("Predator", side = 2, adj = 0.6, line = 0, cex = 1.5)
dev.off()


#######################################
## Hypocarnivores

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Hypocarnivore 2 Ma Intervals 64 to 0 Ma 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Hypocarnivore Draft NALMA 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

### raw
par(mfrow = c(1,2))
prop_herb.34Ma <- prop_herb[rownames(prop_herb) %in% rownames(intervals.34Ma),]
ungulates <- prop_herb.34Ma
predators <- prop_pred_hypo

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black',
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
# text(p1$corrPos$x, p1$corrPos$y+0.1, round(p1$corr, 2))
#text(c(rep(1,1)), c(5:1), unique(p1$corrPos$yName)) 

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.55, adj = 1.5)
mtext("Hypocarnivores", side = 2, adj = 0.1, line = 3, cex = 1.5)

# First Differences
pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred_hypo), intervals = intervals.34Ma)

UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb.34Ma), intervals = intervals.34Ma)

#UngualteBMGroupDiv_FirstDiff

#pred.prey.DivFirstDiff

#should make list of plots of each catagpry through time to make sure data is stationary

ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
predators <- t(as.matrix(pred.prey.DivFirstDiff))

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black', 
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.5)
#  mtext("Predator", side = 2, adj = 0.6, line = 0, cex = 1.5)
dev.off()

#######################################
## Hypercarnivores + Mesocarnivores

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Hyper and Mesocarnivores 2 Ma Intervals 64 to 0 Ma 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Hyper and Mesocarnivores Draft NALMA 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

### raw
par(mfrow = c(1,2))
prop_herb.34Ma <- prop_herb[rownames(prop_herb) %in% rownames(intervals.34Ma),]
ungulates <- prop_herb.34Ma
predators <- prop_pred_hyper.meso

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black',
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
# text(p1$corrPos$x, p1$corrPos$y+0.1, round(p1$corr, 2))
#text(c(rep(1,1)), c(5:1), unique(p1$corrPos$yName)) 

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.55, adj = 1.5)
mtext("Hyper- + Meso-", side = 2, adj = 0.1, line = 3, cex = 1.5)

# First Differences
pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred_hyper.meso), intervals = intervals.34Ma)

UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb.34Ma), intervals = intervals.34Ma)

#UngualteBMGroupDiv_FirstDiff

#pred.prey.DivFirstDiff

#should make list of plots of each catagpry through time to make sure data is stationary

ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
predators <- t(as.matrix(pred.prey.DivFirstDiff))

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black', 
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.5)
#  mtext("Predator", side = 2, adj = 0.6, line = 0, cex = 1.5)
dev.off()


#######################################
## Hypercarnivores + Hypocarnivores

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Hyper and Hypocarnivores 2 Ma Intervals 64 to 0 Ma 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Hyper and Hypocarnivores Draft NALMA 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

### raw
par(mfrow = c(1,2))
prop_herb.34Ma <- prop_herb[rownames(prop_herb) %in% rownames(intervals.34Ma),]
ungulates <- prop_herb.34Ma
predators <- prop_pred_hyper.hypo

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black',
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
# text(p1$corrPos$x, p1$corrPos$y+0.1, round(p1$corr, 2))
#text(c(rep(1,1)), c(5:1), unique(p1$corrPos$yName)) 

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.55, adj = 1.5)
mtext("Hyper- + Hypo-", side = 2, adj = 0.1, line = 3, cex = 1.5)

# First Differences
pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred_hyper.hypo), intervals = intervals.34Ma)

UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb.34Ma), intervals = intervals.34Ma)

#UngualteBMGroupDiv_FirstDiff

#pred.prey.DivFirstDiff

#should make list of plots of each catagpry through time to make sure data is stationary

ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
predators <- t(as.matrix(pred.prey.DivFirstDiff))

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black', 
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.5)
#  mtext("Predator", side = 2, adj = 0.6, line = 0, cex = 1.5)
dev.off()


#######################################
## Mesocarnivores + Hypocarnivores

if(interval.type == "bins")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Meso and Hypocarnivores 2 Ma Intervals 64 to 0 Ma 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png("/Users/emdoughty/Dropbox/Proposal/Chapter 1 Figures/Correlation Matrices Meso and Hypocarnivores Draft NALMA 2022_4_21.png",
      width = 7, height = 3, units = "in", res = 200)
}

### raw
par(mfrow = c(1,2))
prop_herb.34Ma <- prop_herb[rownames(prop_herb) %in% rownames(intervals.34Ma),]
ungulates <- prop_herb.34Ma
predators <- prop_pred_meso.hypo

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black',
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
# text(p1$corrPos$x, p1$corrPos$y+0.1, round(p1$corr, 2))
#text(c(rep(1,1)), c(5:1), unique(p1$corrPos$yName)) 

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.55, adj = 1.5)
mtext("Meso- + Hypo-", side = 2, adj = 0.1, line = 3, cex = 1.5)

# First Differences
pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred_meso.hypo), intervals = intervals.34Ma)

UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb.34Ma), intervals = intervals.34Ma)

#UngualteBMGroupDiv_FirstDiff

#pred.prey.DivFirstDiff

#should make list of plots of each catagpry through time to make sure data is stationary

ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
predators <- t(as.matrix(pred.prey.DivFirstDiff))

all.val <- cbind(predators, ungulates)

corr.results.both <- cor(predators, ungulates, method = "spearman")
cor.p <- cor.mtest(all.val)$p

cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]

corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                   cl.align.text = 'l', 
                   #addCoef.col = 'black', 
                   addgrid.col = 'black',
                   method = "color",
                   na.label = "-",
                   sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                   pch.cex = 1,
                   tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1

#  mtext("Large Herbivore", side = 3, line = -3, cex = 1.5)
#  mtext("Predator", side = 2, adj = 0.6, line = 0, cex = 1.5)
dev.off()

