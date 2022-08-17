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
date.save <- "2022_6_3"

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

######################################################################################################################


#bring over the setup for other objects
###########################################################################################################################

pred_prey_ratio <- rowSums(prop_pred)/rowSums(prop_herb)
pred_prey_ratio.hyper <- rowSums(prop_pred_hyper)/rowSums(prop_herb)[rownames(prop_herb) %in% rownames(intervals.34Ma)]
pred_prey_ratio.meso <- rowSums(prop_pred_meso)/rowSums(prop_herb)[rownames(prop_herb) %in% rownames(intervals.34Ma)]
pred_prey_ratio.hypo <- rowSums(prop_pred_hypo)/rowSums(prop_herb)[rownames(prop_herb) %in% rownames(intervals.34Ma)]
pred_prey_ratio.hyper.meso <- rowSums(prop_pred_hyper.meso)/rowSums(prop_herb)[rownames(prop_herb) %in% rownames(intervals.34Ma)]
pred_prey_ratio.hyper.hypo <- rowSums(prop_pred_hyper.hypo)/rowSums(prop_herb)[rownames(prop_herb) %in% rownames(intervals.34Ma)]
pred_prey_ratio.meso.hypo <- rowSums(prop_pred_hyper.hypo)/rowSums(prop_herb)[rownames(prop_herb) %in% rownames(intervals.34Ma)]

continental_richness <- rowSums(prop_pred) + rowSums(prop_herb)

if(interval.type == "bins")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/PredPrey Ratio Cenozoic Ma Intervals 64 to 0 Ma ", date.save,".png"),
      width = 8, height = 5, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/PredPrey Ratio Cenozoic NALMA ", date.save,".png"),
      width = 8, height = 5, units = "in", res = 200)
}

do.subepochs=TRUE
overlay.labels = FALSE
overlay.color=TRUE
thisAlpha.intervals=0.33
thisAlpha.text = 0.75
borderCol="black"
invertTime=FALSE
scale.cex=0.75
scale.headers = 0.95
text.offset = 0.025

par(mfrow = c(2,1), mar=c(2,4,0,2), mgp=c(2, 1,0), oma = c(2,2,2,0), bg=NA)
#par(mfrow = c(1,1), mar=c(3,3,0.5,0), mgp=c(2, 0.75 ,0), oma = c(0,1,0,1), bg = NA)
plot(1, 1, type="n", 
     xlim= c(65,0), xaxp = c(70,0,14), xlab = "", yaxt="n", xaxt = "n",
     ylim=c(0,1.2), yaxp = c(0,1.2,6), ylab = "")
overlayCzTimescale(do.subepochs= do.subepochs, color = overlay.color, thisAlpha.text = thisAlpha.text, thisAlpha.intervals = thisAlpha.intervals, borderCol = borderCol, invertTime = invertTime, scale.cex = scale.cex, scale.headers = scale.headers, text.offset = text.offset)

lines(rowMeans(intervals), pred_prey_ratio, lwd= 3)
axis(1,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 0.75)
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,1.2,by = 0.2)), lwd = 2, tck = -0.05, las = 1,hadj = 1.5)
#

dev.off()



if(interval.type == "bins")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/PredPrey Ratio Draft 2 Ma Intervals 64 to 0 Ma ", date.save,".png"),
      width = 2.5, height = 5.5, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/PredPrey Ratio Draft NALMA ", date.save,".png"),
      width = 2.5, height = 5.5, units = "in", res = 200)
}

check.diets <- pred_prey_ratio.hyper + pred_prey_ratio.meso + pred_prey_ratio.hypo

par(mar = c(0.5, 2, 0.5, 1), mfrow = c(4,1), oma = c(2,2,1,0), bg=NA)
plot(continental_richness, pred_prey_ratio, 
     cex = 1, cex.axis = 1, xaxt='n', xlim = c(20, 135), 
     ylim = c(0,1.1),
     xaxt = 'n', yaxt = 'n',
     xlab = "Continental Richness", ylab= "Predator Richness/Prey Richness", type = "n", las = 1,
     col.lab = "white")
text(continental_richness, pred_prey_ratio, labels = names(pred_prey_ratio), cex = 0.75)
#mtext(text = "All Predators", side = 3, line = 2.5, cex = 0.75, col = "white")

#if(interval.type %in% "bins") mtext(text = "2Ma-Interval Analysis", side = 3, line = 0.5, cex = 0.75)
#if(interval.type %in% "nalma") mtext(text = "sub-NALMA Analysis", side = 3, line = 0.5, cex = 0.75)

axis(side = 1, at = seq(0, 140, 20), labels = FALSE, tick = TRUE,
     col.ticks="white", col.axis = "white")
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,1.2,by = 0.2)), lwd = 2, tck = -0.05, las = 1,hadj = 1.5)


plot(continental_richness[rownames(prop_herb) %in% rownames(intervals.34Ma)], pred_prey_ratio.hyper, 
     cex = 1, cex.axis = 1, xaxt='n', xlim = c(20, 135), 
     ylim = c(0,1.1), yaxp = c(2, 0, 10),
     xaxt = 'n', yaxt = 'n',
     xlab = "Continental Richness", ylab= "Predator Richness/Prey Richness", type = "n", las = 1,
     col.lab = "white")
text(continental_richness[rownames(prop_herb) %in% rownames(intervals.34Ma)], 
     pred_prey_ratio.hyper, labels = names(pred_prey_ratio.hyper), cex = 0.75)
#mtext(text = "Hypercarnivores", side = 3, line = 2.5, cex = 0.75)
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, tick = TRUE,
     col.ticks="white", col.axis = "white")
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,1.2,by = 0.2)), lwd = 2, tck = -0.05, las = 1,hadj = 1.5)


plot(continental_richness[rownames(prop_herb) %in% rownames(intervals.34Ma)], pred_prey_ratio.meso, 
     cex = 1, cex.axis = 1, xaxt='n', xlim = c(20, 135), 
     ylim = c(0,1.1),
     xaxt = 'n', yaxt = 'n',
     xlab = "Continental Richness", ylab= "Predator Richness/Prey Richness", type = "n", las =1,
     col.lab = "white")
text(continental_richness[rownames(prop_herb) %in% rownames(intervals.34Ma)], 
     pred_prey_ratio.meso, labels = names(pred_prey_ratio.meso), cex = 0.75)
#mtext(text = "Mesocarnivores", side = 2, line = 2.5, cex = 0.75)
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, tick = TRUE,
     col.ticks="white", col.axis = "white")
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,1.2,by = 0.2)), lwd = 2, tck = -0.05, las = 1,hadj = 1.5)


plot(continental_richness[rownames(prop_herb) %in% rownames(intervals.34Ma)], pred_prey_ratio.hypo, 
     cex = 1, cex.axis = 1, xlim = c(20, 135), 
     ylim = c(0,1.1),
     xaxt = 'n', yaxt = 'n',
     xlab = "Continental Richness", ylab= "Predator Richness/Prey Richness", type = "n", las =1,
     col.lab = "white")
text(continental_richness[rownames(prop_herb) %in% rownames(intervals.34Ma)], 
     pred_prey_ratio.hypo, labels = names(pred_prey_ratio.hypo), cex = 0.75)
#mtext(text = "Hypocarnivores", side = 2, line = 2.5, cex = 0.75)
axis(side = 1, at = seq(0, 140, 20), labels = TRUE, tick = TRUE,
     col.ticks="white", col.axis = "white")
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,1.2,by = 0.2)), lwd = 2, tck = -0.05, las = 1,hadj = 1.5)


dev.off()

###################################

#Plots for Jon's code
subsamp <- FALSE

#load("~/Dropbox/Proposal/Chapter 2 Figures/OccsBox_2022_4_22.Rdata")
load("~/Dropbox/Proposal/Chapter 2 Figures/OccsBox_do.subample=FALSE_2022_4_22.Rdata")
# occ.box <- data.frame(occ.box, prop_ung=round(occ.box[,"n_occs_ung"]/occ.box[,"n_occs"], digits=3), prop_carn=round(occ.box[,"n_occs_carn"]/occ.box[,"n_occs"], digits=3))


if(interval.type == "bins")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/Occupancy Subsample=",subsamp," 2 Ma Intervals 64 to 0 Ma ", date.save,".png"),
      width = 8, height = 3.5, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/Occupancy Subsampled=",subsamp," NALMA ", date.save,".png"),
      width = 8, height = 3.5, units = "in", res = 200)
}

#quartz(height = 5.5, width = 8.5)
par(mfrow=c(1,1), mar=c(4,4,1.6,4), mgp=c(2, 1,0), bg=NA) #, height = 2.15, width = 8.21)

plot(apply(intervals, 1, mean), 
     occ.box[2,,"n_occs_ung"], 
    # xlim=c(65, 0), ylim=c(0,500),  #ylim=c(0,750), 
     xlim=c(65, 0), ylim=c(0,10000),
     xaxp = c(70,0,14), yaxp = c(0, 500, 5), #yaxp = c(0, 4000, 8), #
     type="n", yaxt = 'n', xaxt = 'n', ylab= "", xlab="") #, 
    # xlab="Time (Ma", ylab="Number of occurrences with taxon")
overlayCzTimescale(do.subepochs=TRUE)
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_occs_ung"], rev(occ.box[3,,"n_occs_ung"])), col=adjustcolor("dodgerblue1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_occs_carn"], rev(occ.box[3,,"n_occs_carn"])), col=adjustcolor("firebrick1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_occs_rod"], rev(occ.box[3,,"n_occs_rod"])), col=adjustcolor("darkolivegreen1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_occs"], rev(occ.box[3,,"n_occs"])), col=adjustcolor("grey25", alpha=0.2))
points(apply(intervals, 1, mean), occ.box[2,,"n_occs_rod"], pch=24, col="darkolivegreen4", type="o", bg="darkolivegreen3")
points(apply(intervals, 1, mean), occ.box[2,,"n_occs_carn"], pch=22, col="firebrick4", type="o", bg="firebrick1")
points(apply(intervals, 1, mean), occ.box[2,,"n_occs_ung"], pch=21, type="o", col="dodgerblue4", bg="dodgerblue1")
points(apply(intervals, 1, mean), occ.box[2,,"n_occs"], pch=21, type="o", col="black", bg="grey25")

axis(1,col.ticks="white", col.axis = "white", cex.axis =1, labels = FALSE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 0.75)
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,10000,by = 2000)), lwd = 2, tck = -0.05, las = 1, hadj = 1.5)

dev.off()


if(interval.type == "bins")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/Proportional Occupancy Subsampled=",subsamp," 2 Ma Intervals 64 to 0 Ma ", date.save,".png"),
      width = 8, height = 3.5, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/Proportional Occupancy Subsampled=",subsamp," NALMA ", date.save,".png"),
      width = 8, height = 3.5, units = "in", res = 200)
}

#quartz(height = 5.5, width = 8.5)
par(mfrow=c(1,1), mar=c(4,4,1.6,4), mgp=c(2, 1,0), bg=NA) #, height = 2.15, width = 8.21)

plot(apply(intervals, 1, mean), 
     occ.box[2,,"prop_occ_ung"], 
     xlim=c(65, 0), ylim=c(0, 0.75), 
     xaxp = c(70,0,14), yaxp = c(0, 1, 10),
     type="n",  yaxt = 'n', xaxt = 'n', ylab= "", xlab="")#, 
    # xlab="Time (Ma", ylab="Proportion of occurrences within guild")
overlayCzTimescale(do.subepochs=TRUE)
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_occ_ung"], rev(occ.box[3,,"prop_occ_ung"])), col=adjustcolor("dodgerblue1", alpha=0.2), border="dodgerblue1")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_occ_carn"], rev(occ.box[3,,"prop_occ_carn"])), col=adjustcolor("firebrick1", alpha=0.2), border="firebrick1")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_occ_other"], rev(occ.box[3,,"prop_occ_other"])), col=adjustcolor("darkgoldenrod1", alpha=0.2), border="darkgoldenrod4")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_occ_rod"], rev(occ.box[3,,"prop_occ_rod"])), col=adjustcolor("darkolivegreen3", alpha=0.2), border="darkolivegreen4")
points(apply(intervals, 1, mean), occ.box[2,,"prop_occ_rod"], pch=24, col="darkolivegreen4", type="o", bg="darkolivegreen3")
points(apply(intervals, 1, mean), occ.box[2,,"prop_occ_carn"], pch=22, col="firebrick4", type="o", bg="firebrick1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_occ_ung"], pch=21, type="o", col="dodgerblue4", bg="dodgerblue1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_occ_other"], pch=21, type="o", col="darkgoldenrod4", bg="darkgoldenrod1")

axis(1,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 0.75)
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,1.2,by = 0.1)), lwd = 2, tck = -0.05, las = 1,hadj = 1.5)

dev.off()



if(interval.type == "bins")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/Proportional Collections Subsampled=",subsamp," 2 Ma Intervals 64 to 0 Ma ", date.save,".png"),
      width = 8, height = 3.5, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/Proportional Collections Subsampled=",subsamp," NALMA ", date.save,".png"),
      width = 8, height = 3.5, units = "in", res = 200)
}

par(mfrow=c(1,1), mar=c(4,4,1.6,4), mgp=c(2, 1,0), bg=NA) #, height = 2.15, width = 8.21)

plot(apply(intervals, 1, mean), 
     occ.box[2,,"prop_col_ung"], 
     xlim=c(65, 0), ylim=c(0, 1), 
     type="n",  yaxt = 'n', xaxt = 'n', ylab= "", xlab="")#, 
    # xlab="Time (Ma", ylab="Proportion of collections with taxon")
overlayCzTimescale(do.subepochs=TRUE)
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_col_ung"], rev(occ.box[3,,"prop_col_ung"])), col=adjustcolor("dodgerblue1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_col_carn"], rev(occ.box[3,,"prop_col_carn"])), col=adjustcolor("firebrick1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_col_rod"], rev(occ.box[3,,"prop_col_rod"])), col=adjustcolor("darkolivegreen1", alpha=0.2))
points(apply(intervals, 1, mean), occ.box[2,,"prop_col_rod"], pch=24, col="darkolivegreen4", type="o", bg="darkolivegreen1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_col_carn"], pch=22, col="firebrick4", type="o", bg="firebrick1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_col_ung"], pch=21, type="o", col="dodgerblue4", bg="dodgerblue1")

axis(1,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 0.75)
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,1.2,by = 0.2)), lwd = 2, tck = -0.05, las = 1,hadj = 1.5)

dev.off()


if(interval.type == "bins")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/Collection Subsampled=",subsamp," 2 Ma Intervals 64 to 0 Ma ", date.save,".png"),
      width = 8, height = 3.5, units = "in", res = 200)
}

if(interval.type == "nalma")
{
  png(paste0("/Users/emdoughty/Dropbox/Proposal/Chapter 2 Figures/Collection Subsampled=",subsamp," NALMA ", date.save,".png"),
      width = 8, height = 3.5, units = "in", res = 200)
}

par(mfrow=c(1,1), mar=c(4,4,1.6,4), mgp=c(2, 1,0), bg=NA) #, height = 2.15, width = 8.21)

plot(apply(intervals, 1, mean), 
     occ.box[2,,"n_cols_ung"], 
     xlim=c(65, 0), ylim=c(0,2000), 
     type="n",  yaxt = 'n', xaxt = 'n', ylab= "", xlab="")#, 
   #  xlab="Time (Ma", ylab="Number of collections with taxon")
overlayCzTimescale(do.subepochs=TRUE, borderCol = "black")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_cols_ung"], rev(occ.box[3,,"n_cols_ung"])), col=adjustcolor("dodgerblue1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_cols_carn"], rev(occ.box[3,,"n_cols_carn"])), col=adjustcolor("firebrick1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_cols_rod"], rev(occ.box[3,,"n_cols_rod"])), col=adjustcolor("darkolivegreen1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_cols"], rev(occ.box[3,,"n_cols"])), col=adjustcolor("grey25", alpha=0.2))
points(apply(intervals, 1, mean), occ.box[2,,"n_cols_rod"], pch=24, col="darkolivegreen4", type="o", bg="darkolivegreen1")
points(apply(intervals, 1, mean), occ.box[2,,"n_cols_carn"], pch=22, col="firebrick4", type="o", bg="firebrick1")
points(apply(intervals, 1, mean), occ.box[2,,"n_cols_ung"], pch=21, type="o", col="dodgerblue4", bg="dodgerblue1")
points(apply(intervals, 1, mean), occ.box[2,,"n_cols"], pch=21, type="o", col="black", bg="grey25")

axis(1,col.ticks="white", col.axis = "white", cex.axis =1, labels = FALSE, at = c(seq(0,65,by = 5)), lwd = 2, tck = -0.05, las = 1, padj = 0.75)
axis(2,col.ticks="white", col.axis = "white", cex.axis =1, labels = TRUE, at = c(seq(0,2000,by = 400)), lwd = 2, tck = -0.05, las = 1, hadj = 1.5)


dev.off()

###########################################################################################################################################
#this for extracting csv from occs









