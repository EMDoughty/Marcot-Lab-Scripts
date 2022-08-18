
#sources for Jon Marcot's code and specimen measurements 
source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/occFns.R")
source("~/Dropbox/code/R/common_src/sampling.R") 
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R") 

source('~/Dropbox/code/R/dentalMeasurements/src/src_dentalDataFns.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_ecologyAnalysisFns.R', chdir = TRUE)

source('~/Dropbox/Code/R/dentalMeasurements/src/src_evanproposal.R')

####################################################################################################################################

occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
#occs <- read.csv("/Users/emdoughty/Dropbox/Code/Occs_2021_9_5.csv")
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Allodesminae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", 
                                 "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", 
                                 "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", 
                                 "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
occs <- occs[!occs$genus %in% c("Enaliarctos", "Pteronarctos", "Kolponomos", "Pacificotaria", "Pinnarctidion", "Pteronarctos"), ]
occs <- occs[!occs$accepted_name %in% c("Archaeoceti", "Pinnipedia", "Imagotariinae"), ]
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

#remove duplicate taxa occurrences within same collection; this happens when synonyms 
if(this.rank %in% names(occs)) { occs <- occs[!(occs[,this.rank] != "" & duplicated(occs[,c("collection_no", this.rank)])),]
} else if (this.rank == "species") {
  occs <- occs[!(occs$accepted_rank=="species" & duplicated(occs[,c("collection_no", "accepted_name")])),]
}

################################################  
#Settings
run.taxon <- "ungulates" #"ungulates"
this.rank <- "genus" #"genus"
interval.type <- "bins" #"nalma"
do.parallel <- TRUE
if (do.parallel) require(parallel)
save.files <- FALSE

herbivore.size.cat <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)

predator.size.cat  <- c(-Inf, 0, 0.845098, 1.322219, 2, Inf) #includes <1kg, adn 1-7 kg categories
#predator.size.cat  <- c(-Inf,0.845098,1.322219,2,Inf) #<7 kg as single category
pred.cat.name <- "allCateg" #"sub7kg"

#the primary location of where the output files will go.  Make sure this only includes the file destination as the rest of the filename is concatenated in its respective section.
save.pathname <- "~/Dropbox/Code/R/Results/"

#toggle so one can fire and forget without having to go back and forth hunting for each section individually
analysis.toggle <- c("bmHandley") #"taxHandley",

#if you want to load a repIntOccs or repIntTaxa from file put the pathname as this object.  otherwise keep as NUll to make a new repIntOccs and repIntTaxa using the settings below
repIntLoad <-
  #"/Users/emdoughty/Dropbox/Code/R/Results/repIntMaster__this.rank=species_timebin=2Mabins_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Fri Mar  4 21:58:08 2022 ------##.Rdata"
NULL 
#"/Users/emdoughty/Dropbox/Code/R/Results/NALMA/repIntMaster__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 00:51:19 2022 ------##.Rdata"
#NULL

################################################

######################################################################################################################################################
#Ungualtes
##################
if(run.taxon == "ungulates")
{
  measure.mat <- getMeasureMatWithBodyMasses()
  
  ####################################################################################################################################
  #### reduces matrix to just the focal order(s)
  ####################################################################################################################################
  
  archaic.ung <- read.csv("/Users/emdoughty/Dropbox/Code/ArchaicUngulate_UploadFile_2021_4_29.csv")
  
  archaic.Mat.tot <- getMeasureMatCondylarths(data.raw = archaic.ung, occs = occs, 
                                              col.order = colnames(measure.mat), 
                                              all.bm = FALSE,
                                              regression = "ArchaicNonselenodonts")

  measure.mat <- rbind(measure.mat, archaic.Mat.tot)
  
  
  
  if(this.rank=="genus") 
  {
    measure.mat <- makeOneGenusMatFromSpecimenMat(measure.mat) # need to reassign family and reg.vec fields
  
    measure.mat <- measure.mat[,!colnames(measure.mat) %in% c("genus", "reg.vec")] #remove genus and reg.vec fields
    
    #add an order column
    for(xx in unique(occs$order))
    {
      measure.mat$order[measure.mat$taxon %in% unique(occs$genus[occs$order %in% xx])] <- xx
    }
    
    #family
    for(xx in unique(occs$family))
    {
      measure.mat$family[measure.mat$taxon %in% unique(occs$genus[occs$family %in% xx])] <- xx
    }
  }
  
  focal.order <- c("Artiodactyla", "Perissodactyla", 
                   "Proboscidea", 
                   "Dinocerata", 
                   "Tillodontia")
  focal.family <- unique(occs[occs$order %in% focal.order,]$family)
  #search through those without order
  add.family <- c("Arctocyonidae", "Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae", #Condylarths
                  "Conoryctidae", "Stylinodontidae") #Taenodonts
  focal.family <- c(as.character(focal.family), add.family)
  focal.family <- focal.family[!focal.family %in% ""]
  focal.family <- focal.family[order(focal.family)]
  
 # bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
  bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & (occs$order %in% focal.order | occs$family %in% add.family)), c("order","family", "genus", "accepted_name")])
 # bigList.cond <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% add.family), c("order","family", "genus", "accepted_name")])
  #bigList <- rbind(bigList,bigList.cond)
  
  bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
  # bigList[order(bigList$family, bigList$accepted_name),]
  shortFam <- sort(unique(bigList$family[bigList$family %in% focal.family]))	
  
  bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
  measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name[bigList$family %in% shortFam], ]
}
#################

#################
#Predators
################
if(run.taxon == "carnivores")
{
  pred.data <- read.csv("/Users/emdoughty/Dropbox/Code/R/DentalMeasurements/dat/predator_data_final.csv")
  colnames(pred.data) <- c("family", "taxon",	"max_ma","min_ma",	"m1L",	"rbl", "bodyMass", 	"Citation") #"BM_all_carnivoran","BM_extant_reg")
  pred.data$taxon <- getCurrentTaxa(tax.vec = pred.data$taxon)
  pred.data$taxon <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.data$taxon)
  rownames(pred.data) <- pred.data$taxon 
  pred.data$genus <- sapply(pred.data$taxon, function(x) { strsplit(x, "_")[[1]][1]}) #need a genus field to make genus only runs work
  
  pred.data[,c("bodyMass")] <- log10(pred.data[,c("bodyMass")])
  pred.data <- pred.data[is.finite(pred.data$bodyMass),]
  
  if (this.rank=="genus") pred.data <- makeOneGenusMatFromSpecimenMat(pred.data)
  
  focal.order <- c("Carnivora", "Creodonta","Hyaenodonta")
  #occsPred <- occs[occs$order %in% focal.orderPred,]
  focal.family <- unique(occs[occs$order %in% focal.order,]$family)
  add.family <- c("Viverravidae", 
                  "Mesonychidae")
  
  focal.family <- c(as.character(focal.family), add.family)
  focal.family <- focal.family[!focal.family %in% ""]
  focal.family <- focal.family[order(focal.family)]
  
  bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
  bigList.cond <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% add.family), c("order","family", "genus", "accepted_name")])
  bigList <- rbind(bigList,bigList.cond)
  
  bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
  # bigList[order(bigList$family, bigList$accepted_name),]
  shortFam <- sort(unique(bigList$family[bigList$family %in% focal.family]))	
  
  bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
  measure.mat <- pred.data[pred.data$taxon %in% bigList$accepted_name[bigList$family %in% shortFam], ]
  measure.mat <- measure.mat[is.finite(measure.mat$bodyMass),]
}

###############################################################################################################################################

if(interval.type == "bins")
{
  int_length <- 2
  intervals <- makeIntervals(0, 64, int_length)
  intList <- listifyMatrixByRow(intervals)
  
  save.path.bins <- paste0(int_length,"Ma",interval.type)
}
#######################

if(interval.type == "nalma")
{
  nalma.mark <- read.csv("/Users/emdoughty/Dropbox/Proposal/Proposal/NOW_intervals_edit.csv")
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
  intervals <- nalma.mark
  
  save.path.bins <- paste0(interval.type)
}

##############################################################################################################################################

if(is.null(repIntLoad))
{
  do.parallel <- TRUE
  if (do.parallel) require(parallel)
  reps <- 100
  do.subsample <- TRUE
  quota <- 0.4
  do.disparity <- FALSE
  bootstrapSpecimens <- FALSE
  bootstrapSpecies <- FALSE
  bootstrapSpeciesWithinIntervals <- FALSE
  plotHist <- FALSE
  do.heuristic <- TRUE
  extra.intvs <- 0
  do.rangethrough <- TRUE
  save.files <- TRUE
  
  if (bootstrapSpecies) holderMat <- measure.mat
  
  if (plotHist) {
    quartz("Guild Histograms")
    par(mfrow=c((nrow(intervals)), 3), mar=c(0,0,0.75,0), cex.axis=0.5, cex.main=0.75)
  }
  
  repIntOccs <- list()
  
  ################################################################################################################################################
  #### get species within intervals
  ################################################################################################################################################
  
  for (rep in seq_len(reps)) {
    cat("Beginning Rep", rep, "of", reps, "...\r")
    ##################################################We need to update this sbootstrap section
    if (bootstrapSpecimens) {
      measure.mat <- specimenMat[sample.int(nrow(specimenMat), size=nrow(specimenMat), replace=TRUE),]
      measure.mat <- aggregate(measure.mat, by=list(taxon=specimenMat$taxon), mean, na.rm=TRUE)
      # measure.mat <- measure.mat[,apply(!sapply(measure.mat, is.na), 2, any)]
      rownames(measure.mat) <- measure.mat$taxon
      measure.mat[sapply(measure.mat, is.nan)] <- NA
      # measure.mat<-cbind(measure.mat, cbind(FO=vector(length=nrow(measure.mat), mode="numeric"), LO=vector(length=nrow(measure.mat), mode="numeric")))
      measure.mat[,"reg"] <- as.character(famList$reg[match(measure.mat$taxon,famList$taxon)])
      measure.mat[,"bodyMass"] <- makeBodyMasses(measure.mat, regList, best.only=TRUE)
      # measure.mat[,"PC2"] <- pcVec[match(measure.mat$taxon, names(pcVec))]
      # measure.mat[,"PC3"] <- pcaLo$x[match(measure.mat$taxon, rownames(pcaLo$x)),3]
    }
    if (bootstrapSpecies) measure.mat <- holderMat[sample.int(n=nrow(measure.mat), size=nrow(measure.mat), replace=TRUE),]
    
    col.dates <- getCollectionAgesFromOccs(occs=occs[, c("collection_no", "max_ma", "min_ma")], random=TRUE)
    occDates <- col.dates$collection_age[match(occs$collection_no, col.dates$collection_no)]
    intOccs <- apply(intervals, 1, function(thisIntv) occs$occurrence_no[occDates > thisIntv[1] & occDates <= thisIntv[2]]) # greater than ageTop, less than or equal to ageBase
    # intTaxa <- sapply(intOccs, function(x) unique(occs$accepted_name[occs$occurrence_no %in% x]))
    # x <- intOccs
    # intSp <- sapply(intOccs, function(x) match(sort(unique(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$accepted_rank  == "species" & occs$occurrence_no %in% x]))), measure.mat$taxon))
    # intSp <- sapply(intOccs, function(x) match(sort(unique(occs$accepted_name[occs$accepted_rank  %in% c("genus", "species") & occs$occurrence_no %in% x]))), measure.mat$taxon))
    # intSp <- sapply(intOccs, function(x) sort(unique(as.character(occs$accepted_name[occs$occurrence_no %in% x]))))
    
    #which(occs$occurrence_no %in% x == TRUE) # none are being returned as TRUE
    
    if (do.subsample) { 
      nOccs <- sapply(intOccs, length)
      # nTaxa <- sapply(intSp, length)			### if you want to set the quota no lower than the maximum number of SIB taxa; intSp is required for this to work, so has to be done above
      nTaxa <- 0									### set to zero to simply set the quota to the minimum number of occurrences
      quota <- max(c(max(nTaxa), min(nOccs)))		### quota is either the maximum number of observed taxa, or the minimum number of occurrences
      cat("Subsampling quota set to", quota, "occurrences")
      
      intOccs <- lapply(X=intOccs, FUN=sample, size=quota)
    }
    
    repIntOccs[[rep]] <- intOccs 
  }
  
  repIntTaxa <- getRepIntTaxaFromRepIntOccs(repIntOccs, this.rank=this.rank, do.rangethrough=do.rangethrough)
  
  print("Completed getting taxa with intervals")
  
  ###################################################################################################################################
  if(save.files)
  {
    if(Sys.info()["sysname"] == "Darwin"){
      if(run.taxon == "carnivores") run.taxon <- paste0(run.taxon, "_",pred.cat.name)
      save(repIntTaxa, repIntOccs, intervals, reps, do.subsample, quota,do.disparity, 
           bootstrapSpecimens,bootstrapSpecies,bootstrapSpeciesWithinIntervals ,
           plotHist,do.heuristic,extra.intvs,do.rangethrough,
           file=paste0(save.pathname,"repIntMaster_",
                       "_this.rank=", this.rank,
                       "_timebin=", save.path.bins,
                       "_SampleStandardized=", do.subsample, 
                       "_Reps=", reps, gsub("-","_",Sys.info()["nodename"]),
                       timestamp(),".Rdata"))
      #load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    } else if(Sys.info()["sysname"] == "Windows"){
      save(repIntTaxa, repIntOccs, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/repIntTaxa_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
      # load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    }
  }
}

if(!is.null(repIntLoad)) load(repIntLoad)

########################################################################################################################################################################

#make a fuction to designate specific unranked clades (e.g. Taeniodonta) so that these can be incorporated more easily into bigList
#taenodonts seem to be only group to have this issue other at least have order or family designations that can be called upon

data.coverage(temp,
              bigList.check,
              focal.archaic,
              "family")

data.coverage(collected.dat = measure.mat,  
              occs.list = occs, 
              taxon.vec = c("Artiodactyla", "Perissodactyla"), 
              column.taxon = "order",
              this.rank = "species")
data.coverage(collected.dat = measure.mat,  
              occs.list = occs, 
              taxon.vec = add.family, 
              column.taxon = "family",
              this.rank = "species")
nrow(probo.mat)


###########################################################################################################################
#export clade species list

clade.column = "order"
#focal.clade = c("Arctocyonidae", "Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae")
focal.clade = "Proboscidea"
output.name = paste("Condylarthra",".csv", sep = "")

exp.occs <- occs[occs[,clade.column] %in% focal.clade & occs$accepted_rank %in% "species",c("accepted_name",	"order", "family", "genus",	"max_ma",	"min_ma")]
exp.occs <- unique(exp.occs)

test1 <- aggregate(x = exp.occs$max_ma,
                   by = list(exp.occs$accepted_name, exp.occs$order, exp.occs$family, exp.occs$genus),
                   FUN = max)
test2 <- aggregate(x = exp.occs$min_ma,
                   by = list(exp.occs$accepted_name, exp.occs$order, exp.occs$family, exp.occs$genus),
                   FUN = min)
test.comb <- cbind(test1, test2[,"x"]); colnames(test.comb) <- c("accepted_name",	"order", "family", "genus",	"max_ma",	"min_ma")

test.comb <- test.comb[order(test.comb$family, test.comb$genus),]

#test.comb$bodyMass <- match(test.comb$accepted_name, )

write.csv(test.comb, file = paste0("/Users/emdoughty/Dropbox/Code/R/Results/Clade species lists/", output.name))


########################################################################################################################

data.mat <- measure.mat
data.mat$BMcat <-data.mat$bodyMass

Janis2000 <- read.csv("/Users/emdoughty/Dropbox/Papers/Datasets/Janis 2000 Appendix.csv")

temp <- unique(Janis2000[, c(5,7)])
temp.save <- temp
temp$OldGenus <- temp$Genus
temp$Genus <- getCurrentTaxa(tax.vec = temp$Genus)

#make Janis2000 data into same format as repIntSp

write.csv(temp, "/Users/emdoughty/Dropbox/Papers/Datasets/temp.csv")

data.mat$Janis2000_1 <- NA

data.mat$Janis2000_2 <- NA #archaeotherium has species that are in 4 and 2

for(xx in seq(1, length(bmBreaks)-1, 1)){
  data.mat$BMcat[data.mat$bodyMass > bmBreaks[xx] & data.mat$bodyMass < bmBreaks[xx+1]] <- xx
}  

for(xx in unique(Janis2000$Genus))
{
  print(paste0(xx, temp[temp$Genus %in% xx, 2], sep=" "))
  data.mat$Janis2000_1[data.mat$taxon %in% xx] <- temp[temp$Genus %in% xx, 2][1]
  if(length(temp[temp$Genus %in% xx, 2]) > 1) data.mat$Janis2000_2[data.mat$taxon %in% xx] <- temp[temp$Genus %in% xx, 2][2] #archaeotherium has species that are in 4 and 2
}

save.data.mat <- data.mat

paleogene.data.mat <- data.mat[data.mat$taxon %in% unique(paleogene.occs$genus),]

write.csv(paleogene.data.mat[,c("bodyMass", "BMcat","Janis2000_1", "Janis2000_2")], "/Users/emdoughty/Dropbox/Papers/Datasets/CompareMeasure.mat&Janis2000.csv")

##########################################################################################################################

#merge with Janis 2000 values
###fill in missing taxa using Janis 2000 values but keep our calls for those genera that we have.
#step 1 make column for size category calls

#step 2 merge janis 20000 with measure.mat

#step3 apply categories so that I can recreate countCube and 

#will need to alter how countCube is generated
#countCube_Alternaive
#or append Janis calls on  measure.mat....make new column for categories and alter code to make countcube to accept that

measure.mat$SizeCat <- measure.mat$bodyMass
bmBreaks <- herbivore.size.cat

for(xx in seq(1, length(bmBreaks)-1, 1)){
  measure.mat$SizeCat[measure.mat$bodyMass > bmBreaks[xx] & measure.mat$bodyMass < bmBreaks[xx+1]] <- xx
}  

#match with Janis 2000 but create new rows for janis uniques
Janis2000 <- read.csv("/Users/emdoughty/Dropbox/Papers/Datasets/Janis 2000 Appendix.csv")

temp <- unique(Janis2000[, c(5,7)])
temp.save <- temp
temp$OldGenus <- temp$Genus
temp$Genus <- getCurrentTaxa(tax.vec = temp$Genus)
temp <- unique(temp)

Janis.2add <- temp[!temp$Genus %in% measure.mat$taxon,] #get genera that are unique to Janis 2000, mainly those genera from the other archaic we haven't sampled (as of 7/13/2022)

Janis.2add <- unique(Janis.2add[,1:2]); colnames(Janis.2add) <- c("taxon", "SizeCat")
Janis.mat <- matrix(nrow=nrow(Janis.2add), ncol=ncol(measure.mat)-2); colnames(Janis.mat) <- colnames(measure.mat[,!colnames(measure.mat) %in% colnames(Janis.2add)]) #get Janis.2add to have same columns as measure.mat
Janis.mat <- cbind(Janis.mat, Janis.2add)
Janis.mat <- Janis.mat[,colnames(measure.mat)] #get in proper order

probo.mat <- as.data.frame(matrix(nrow=length(unique(occs$accepted_name[occs$order %in% "Proboscidea" & occs$accepted_rank %in% this.rank])), ncol=ncol(measure.mat))); colnames(probo.mat) <- colnames(measure.mat)
probo.mat$taxon <- unique(occs$accepted_name[occs$order %in% "Proboscidea" & occs$accepted_rank %in% this.rank]); probo.mat <- probo.mat[!probo.mat$taxon %in% "",]
probo.mat$SizeCat <- 5

measure.mat <- rbind(measure.mat, Janis.mat, probo.mat)

##########################################################################################################################

#make alternative countcube that works with the SizeCat column rather than using hist()
#get a count of 
repIntTest <- repIntTaxa
this.rep <- repIntTest 
this.intv <- this.rep[[1]]

countCube_herb <- sapply(repIntTaxa, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(measure.mat[,"SizeCat"][match(this.intv, measure.mat$taxon)], 
         breaks= c(-Inf, seq(1.5, 4.5,1), Inf), plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

sizecateg <- c("<5 kg", 
               "5-25 kg",
               "25-150kg", 
               "150-500 kg", 
               ">500 kg")

dimnames(countCube_herb) <- list(sizecateg, rownames(intervals), NULL)

prop_herb <- t(apply(countCube_herb, c(1,2), median, na.rm=TRUE))
colnames(prop_herb)[colnames(prop_herb)==""] <- "indeterminate"

plotStackedRichness(this.box=prop_herb, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    col.axis = "black", col.lab = "black", xaxt = 'n', yaxt = 'n',
                    cex.axis = 1.5, cex.lab = 1, las = 1,
                    ylab = "", #"Species Richness (Number of Subtaxa)",
                    xlim = c(65, min(intervals, na.rm=TRUE)), ylim = c(0, max(rowSums(prop_herb))+20), xaxp = c(65, 0, 13), yaxp = c(0, 120, 6),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)
axis(1,col.ticks="black", col.axis = "black", cex.axis =2, labels = TRUE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 1.5)
axis(2,col.ticks="black", col.axis = "black", cex.axis =2, labels = TRUE, at = c(seq(0,100,by = 20)), lwd = 2, tck = -0.05, las = 1, hadj = 1.5)

##########################################################################################################################

temp.mat <- pred.data

meson.mat <- read.csv("/Users/emdoughty/Dropbox/Code/R/Results/Mesonychidae_BodySize.csv")
meson.mat <- meson.mat[meson.mat$family %in% "Mesonychidae", c("accepted_name", "order", "family", "genus", "max_ma", "min_ma", "accepted_rank", "kg")]
if(this.rank == "genus") meson.mat <- makeOneGenusMatFromSpecimenMat(meson.mat)
meson.mat <- meson.mat[!meson.mat$taxon %in% "",]

bmBreaks_pred <- predator.size.cat
meson.mat$SizeCat <- meson.mat$kg

for(xx in seq(1, length(bmBreaks_pred)-1, 1)){
  meson.mat$SizeCat[log10(meson.mat$kg) > bmBreaks_pred[xx] & log10(meson.mat$kg) < bmBreaks_pred[xx+1]] <- xx
}  

attributes(pred.data)



