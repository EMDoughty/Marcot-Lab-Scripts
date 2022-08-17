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
this.rank <- "genus" #"genus" "species"
interval.type <- "bins" #"nalma" "bins
add.Janis2000 <- FALSE
add.probo <- TRUE
add.mesonychid <- TRUE
date.save <- paste0("Janis=", add.Janis2000,"_Probo=", add.probo,"_Mesonichid=",add.mesonychid,"_2022_7_15")

bmBreaks_herb <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
#bmBreaks_herb <- c(-Inf, 1.057, 1.54, 2.0265, 2.5785, Inf ) #k = 5
#bmBreaks_herb <- c(-Inf, 0.844, 1.269, 1.6245, 2.05665, 2.5785, Inf ) #k = 6

bmBreaks_pred <- c(-Inf, 0, 0.845098, 1.322219, 2, Inf) #PPP categories
#bmBreaks_pred <- c(-Inf, 0.8565, Inf)#k = 2
#bmBreaks_pred <- c(-Inf, -0.3235, 0.2505, 0.77, 1.538, Inf) #k = 5
#bmBreaks_pred <- c(-Inf, -0.3235, 0.2505, 0.7235, 1.3205, 2.0895, Inf)#k = 6

#do.pred.diet <- TRUE
#diet.type <- c("hypercarnivore") #, "mesocarnivore", "hypocarnivore")

save.pathname <- "~/Dropbox/Proposal/Proposal/Chapter 1 Figures/"

#if you want to load a repIntOccs or repIntTaxa from file put the pathname as this object.  otherwise keep as NUll to make a new repIntOccs and repIntTaxa using the settings below
  ####Full runs

#MASTER use this for normal runs
if(interval.type %in% "bins"){
  if(this.rank %in% "genus") repIntLoad <- "/Users/emdoughty/Dropbox/Code/R/Results/repIntMaster__this.rank=genus_timebin=2Mabins_SampleStandardized=TRUE_Reps=10000wifi_131_179_8_9.host.ucla.edu##------ Mon Jul 25 17:42:51 2022 ------##.Rdata"
  if(this.rank %in% "species") repIntLoad <- "/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/repIntMaster__this.rank=species_timebin=2Mabins_start=64Ma_end=0Ma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Sun Mar  6 19:30:39 2022 ------##.Rdata"
}

#kmeans
#if(interval.type %in% "bins")  repIntLoad <- "/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/JenksBreaks/repIntMaster__this.rank=species_timebin=2Mabins_start=66_end=0_SampleStandardized=TRUE_Reps=1000Jonathans_MBP.lan##------ Wed Apr 13 23:48:30 2022 ------##.Rdata"


#if(interval.type %in% "nalma") repIntLoad <- "/Users/emdoughty/Dropbox/Code/R/Results/repIntMaster__this.rank=genus_timebin=nalma_SampleStandardized=FALSE_Reps=1000wifi_131_179_8_30.host.ucla.edu##------ Wed Jul 13 15:40:28 2022 ------##.Rdata"

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

####################################################################################################################################

measure.mat <- getMeasureMatWithBodyMasses()
 
archaic.ung <- read.csv("/Users/emdoughty/Dropbox/Code/ArchaicUngulate_UploadFile_2021_4_29.csv")
archaic.Mat.tot <- getMeasureMatCondylarths(data.raw = archaic.ung, occs = occs, 
                                            col.order = colnames(measure.mat), 
                                            all.bm = FALSE,
                                            regression = "ArchaicNonselenodonts")
#archaic.Mat.tot$taxon <- getCurrentTaxa(tax.vec = archaic.Mat.tot$taxon) already happens in getMeasureMatCondylarths

measure.mat <- rbind(measure.mat, archaic.Mat.tot)

#bmBreaks <- herbivore.size.cat

if(this.rank=="genus") 
{
  measure.mat <- makeOneGenusMatFromSpecimenMat(measure.mat) # need to reassign family and reg.vec fields
  
  measure.mat <- measure.mat[,!colnames(measure.mat) %in% c("genus", "reg.vec")] #remove genus and reg.vec fields
}

  #add an order column
  for(xx in unique(occs$order))
  {
    measure.mat$order[measure.mat$taxon %in% unique(occs$accepted_name[occs$order %in% xx & occs$accepted_rank %in% this.rank])] <- xx
  }
  
  #family
  for(xx in unique(occs$family))
  {
    measure.mat$family[measure.mat$taxon %in% unique(occs$accepted_name[occs$family %in% xx & occs$accepted_rank %in% this.rank])] <- xx
  }

measure.mat$SizeCat <- measure.mat$bodyMass

for(xx in seq(1, length(bmBreaks_herb)-1, 1)){
  measure.mat$SizeCat[measure.mat$bodyMass > bmBreaks_herb[xx] & measure.mat$bodyMass < bmBreaks_herb[xx+1]] <- xx
} 

if(add.Janis2000){
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
  
  measure.mat <- rbind(measure.mat, Janis.mat)
  measure.mat <- measure.mat[order(measure.mat$taxon),]
}
  
if(add.probo){
  probo.mat <- as.data.frame(matrix(nrow=length(unique(occs$accepted_name[occs$order %in% "Proboscidea" & occs$accepted_rank %in% this.rank])), ncol=ncol(measure.mat))); colnames(probo.mat) <- colnames(measure.mat)
  probo.mat$taxon <- unique(occs$accepted_name[occs$order %in% "Proboscidea" & occs$accepted_rank %in% this.rank]); probo.mat <- probo.mat[!probo.mat$taxon %in% "",]
  probo.mat$SizeCat <- 5
  
  measure.mat <- rbind(measure.mat, probo.mat)
  measure.mat <- measure.mat[order(measure.mat$taxon),]
}

focal.order <- c("Artiodactyla", "Perissodactyla", 
                 "Proboscidea", 
                 "Dinocerata", 
                 "Tillodontia")
focal.family <- unique(occs[occs$order %in% focal.order,]$family)
#search through those without order
add.family <- c("Arctocyonidae", "Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae") #, #, #Condylarths
              #  "Conoryctidae", "Stylinodontidae") #Taenodonts
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
#measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name[bigList$family %in% shortFam], ]

#########################################################################################################################
pred.data <- read.csv("~/Dropbox/Proposal/Proposal/predator_data_final.csv")
colnames(pred.data) <- c("family", "taxon",	"max_ma","min_ma",	"m1L",	"rbl", "bodyMass", 	"Citation") #"BM_all_carnivoran","BM_extant_reg")


#add mesonychidae
if(add.mesonychid){
  meson.mat <- read.csv("/Users/emdoughty/Dropbox/Code/R/Results/Mesonychidae_BodySize.csv")
  meson.mat <- meson.mat[meson.mat$family %in% "Mesonychidae", c("family", "accepted_name", "max_ma", "min_ma", "kg")]
  meson.mat <- meson.mat[!meson.mat$accepted_name %in% "",]
  meson.mat$m1L <-  meson.mat$rbl <- meson.mat$Citation <- NA
  colnames(meson.mat) <- c("family", "taxon",	"max_ma","min_ma",	"bodyMass","m1L",	"rbl",	"Citation")
  
  meson.mat <- meson.mat[, c("family", "taxon", "max_ma", "min_ma", "m1L", "rbl", "bodyMass", "Citation")]
  
  pred.data <- rbind(pred.data, meson.mat)
}

pred.data$taxon <- getCurrentTaxa(tax.vec = pred.data$taxon)
pred.data$taxon <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.data$taxon)

if(this.rank=="genus") 
{
  pred.data$genus <- unlist(lapply(strsplit(as.character(pred.data$taxon),"_"), function(x) x[1]))
  pred.data <- makeOneGenusMatFromSpecimenMat(pred.data) # need to reassign family and reg.vec fields
  pred.data <- pred.data[,!colnames(pred.data) %in% c("genus", "reg.vec")] #remove genus and reg.vec fields
  
  #add an order column
  for(xx in unique(occs$order))
  {
    pred.data$order[pred.data$taxon %in% unique(occs$genus[occs$order %in% xx])] <- xx
  }
  
  #family
  for(xx in unique(occs$family))
  {
    pred.data$family[pred.data$taxon %in% unique(occs$genus[occs$family %in% xx])] <- xx
  }
}

pred.data[,c("bodyMass")] <- log10(pred.data[,c("bodyMass")])
pred.data <- pred.data[is.finite(pred.data$bodyMass),]

focal.orderPred <- c("Carnivora", "Creodonta","Hyaenodonta", "Acreodi")
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
#pred.data <- pred.data[complete.cases(pred.data),]

#pred.diet <- read.csv("/Users/emdoughty/Dropbox/Proposal/Proposal/Diet Data PPP CJ.csv")
#pred.diet$MASTER_LIST <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.diet$MASTER_LIST)

#pred.data.master  <- pred.data

#pred.data <- pred.data[pred.data$taxon %in% pred.diet$MASTER_LIST,]
#pred.data$Diet <- pred.diet[pred.diet$MASTER_LIST %in% pred.data$taxon, "PPP_diet_2"]


######################################################################################################################
#check # taxa in prey and herbs
length(unique(measure.mat$taxon[measure.mat$taxon %in% occs$accepted_name[occs$accepted_rank %in% this.rank]]))

length(unique(pred.data$taxon[pred.data$taxon %in% occs$accepted_name[occs$accepted_rank %in% this.rank]]))
######################################################################################################################

#rem.col <- c("order", "family","genus","reg.vec")
#measure.mat <- measure.mat[,!colnames(measure.mat) %in% rem.col]

if(add.Janis2000 | add.probo){
  countCube_herb <- sapply(repIntTaxa, function(this.rep) {
    sapply(this.rep, function(this.intv, this.rep) {
      hist(measure.mat[,"SizeCat"][match(this.intv, measure.mat$taxon)], 
           breaks= c(-Inf, seq(1.5, 4.5,1), Inf), plot=FALSE)$counts
    }, this.rep=this.rep)
  }, simplify = "array")
  
} else {
  countCube_herb <- sapply(repIntTaxa, function(this.rep) {
    sapply(this.rep, function(this.intv, this.rep) {
      hist(measure.mat[,"bodyMass"][match(this.intv, measure.mat$taxon)], 
           breaks= bmBreaks_herb, plot=FALSE)$counts
    }, this.rep=this.rep)
  }, simplify = "array")
}
#countCube <- countCube[,,1]

sizecateg <- c("<5 kg", 
               "5-25 kg",
               "25-150kg", 
               "150-500 kg", 
               ">500 kg")

dimnames(countCube_herb) <- list(sizecateg, rownames(intervals), NULL)


#pred.data <- pred.data[,!colnames(pred.data) %in% rem.col]

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
               ">100kg" )

#sizecateg <- c("<7", ">7")
              

dimnames(countCube_pred) <- list(sizecateg, rownames(intervals), NULL)

###############################################################

countBox_herb <- apply(countCube_herb, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

countBox_pred <- apply(countCube_pred, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

##############################################################
prop_herb <- t(apply(countCube_herb, c(1,2), median, na.rm=TRUE))
colnames(prop_herb)[colnames(prop_herb)==""] <- "indeterminate"

prop_pred <- t(apply(countCube_pred, c(1,2), median, na.rm=TRUE))
colnames(prop_pred)[colnames(prop_pred)==""] <- "indeterminate"

###############################################################



###############################################################
#quartz(width = 8, height =8)
if(interval.type == "bins")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Proposal/Chapter 1 Figures/", this.rank, " Richness Plots Draft 2 Ma Intervals 64 to 0 Ma ", date.save,".png"),
     width = 8, height = 9, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Proposal/Chapter 1 Figures/", this.rank, " Richness Plots Draft NALMA ", date.save,".png"),
      width = 8, height = 9, units = "in", res = 200)
}
  
par(mfrow = c(4,1), mar=c(2,4,0,2), mgp=c(2, 1,0), oma = c(2,2,2,0), bg=NA)
plotStackedRichness(this.box=prop_herb, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    col.axis = "black", col.lab = "black", xaxt = 'n', yaxt = 'n',
                    cex.axis = 1.5, cex.lab = 1, las = 1,
                    ylab = "", #"Species Richness (Number of Subtaxa)",
                    xlim = c(65, min(intervals, na.rm=TRUE)), ylim = c(0, max(rowSums(prop_herb))+20), xaxp = c(65, 0, 13), yaxp = c(0, 120, 6),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)
axis(1,col.ticks="black", col.axis = "black", cex.axis =2, labels = FALSE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 1.5)
axis(2,col.ticks="black", col.axis = "black", cex.axis =2, labels = TRUE, at = c(seq(0,100,by = 20)), lwd = 2, tck = -0.05, las = 1, hadj = 1.5)

#mtext(text = "Large Herbivores", side = 2, line = 3.5, cex = 1.25, at = -10)

plotStackedRichness(this.box=prop_herb/rowSums(prop_herb), intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    col.axis = "black", col.lab = "black", xaxt = 'n', yaxt = 'n',
                    cex.axis = 1.5, cex.lab = 1, las = 1,
                    ylab = "",# "Proportion of Subtaxa",
                    xlim = c(65, min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)
axis(1,col.ticks="black", col.axis = "black", cex.axis =2, labels = FALSE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 1.5)
axis(2,col.ticks="black", col.axis = "black", cex.axis =2, labels = TRUE, at = c(seq(0,1,by = 0.2)), lwd = 2, tck = -0.05, las = 1, hadj = 1.5)

plotStackedRichness(this.box=prop_pred, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, 
                    col.axis = "black", col.lab = "black", xaxt = 'n',yaxt = 'n',
                    cex.axis = 1.5, cex.lab = 1, las = 1,
                    ylab = "",# "Species Richness (Number of Subtaxa)",
                    xlim = c(65,  min(intervals, na.rm=TRUE)), ylim = c(0, max(rowSums(prop_pred)+10)), xaxp = c(65, 0, 13), yaxp = c(0, 45, 9),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)
axis(1,col.ticks="black", col.axis = "black", cex.axis =2, labels = FALSE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 1.5)
axis(2,col.ticks="black", col.axis = "black", cex.axis =2, labels = TRUE, at = c(seq(0,40,by = 5)), lwd = 2, tck = -0.05, las = 1,  hadj = 1.5)

#mtext(text = "Predators", side = 2, line = 3.5, cex = 1.25, at = -10)

plotStackedRichness(this.box=prop_pred/rowSums(prop_pred), intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, 
                    numbers.only=FALSE, add.legend=FALSE, prop.ylab = TRUE,
                    col.axis = "black", col.lab = "black", xaxt = 'n', yaxt = 'n',
                    cex.axis = 1.5, cex.lab = 1, las = 1,
                    ylab = "",# "Proportion of Subtaxa",
                    xlim = c(65,  min(intervals, na.rm=TRUE)), ylim = c(0, 1.2), xaxp = c(65, 0, 13), yaxp = c(0, 1, 5),
                    overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = TRUE, thisAlpha.text = 0.75, borderCol = "black", invertTime = FALSE, scale.cex = 1.25, scale.headers = 0.90, text.offset = 0.05)
axis(1,col.ticks="black", col.axis = "black", cex.axis =2, labels = TRUE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 0.75)
axis(2,col.ticks="black", col.axis = "black", cex.axis =2, labels = TRUE, at = c(seq(0,1,by = 0.2)), lwd = 2, tck = -0.05, las = 1,hadj = 1.5)
#mtext(text = "Time (Ma)", side = 1, line = 2, cex = 0.75, las = 1)
dev.off()

#########################################################
# Correlations
##All predator
if(interval.type == "bins")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Proposal/Chapter 1 Figures/Correlation Matrices All predator", this.rank, " 2 Ma Intervals 64 to 0 Ma ", date.save,".png"),
      width = 7, height = 3, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Proposal/Chapter 1 Figures/Correlation Matrices All predator", this.rank, " Draft NALMA ", date.save,".png"),
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

  






