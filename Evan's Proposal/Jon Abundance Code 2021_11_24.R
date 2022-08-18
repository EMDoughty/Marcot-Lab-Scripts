source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/occFns.R")
source("~/Dropbox/code/R/common_src/sampling.R") 
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R") 

source('~/Dropbox/code/R/dentalMeasurements/src/src_dentalDataFns.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_ecologyAnalysisFns.R', chdir = TRUE)

source('~/Dropbox/Code/Recreate Dental Plots Source 2021_6_2.R')

####################################################################################################################################

occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
#occs <- read.csv("/Users/emdoughty/Dropbox/Code/Occs_2021_9_5.csv")
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", 
                                 "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", 
                                 "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", 
                                 "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

interval.type <- "bins"
load.occ.box = FALSE
do.subsample = FALSE

if(interval.type %in% "bins") { 
  repIntLoad <- 
  #"repIntMaster__this.rank=species_timebin=2Mabins_start=64Ma_end=0Ma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Sun Mar  6 19/30/39 2022 ------##.Rdata"
  "/Users/emdoughty/Dropbox/Code/R/Results/2Ma Interval/repIntMaster__this.rank=species_timebin=2Mabins_start=66_end=0_SampleStandardized=FALSE_Reps=10000Jonathans_MBP.lan##------ Thu Apr 14 20:19:28 2022 ------##.Rdata"
}

if(interval.type %in% "nalma") repIntLoad <- "/Users/emdoughty/Dropbox/Code/R/Results/NALMA/repIntMaster__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Tue Mar 15 00:51:19 2022 ------##.Rdata"

load(repIntLoad)

if(interval.type == "bins")
{
  int_length <- 2
  intervals <- makeIntervals(0, 64, int_length) # for some reason it adds an additional interval when doing odd midpoints (i.e. 0 to 64 Ma have a 64-66Ma interval while 1-64 Ma will not)....its late i may be crazy
  intList <- listifyMatrixByRow(intervals)
  
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

################################################################################################################################################

# sapply(repIntOccs, function(x) sapply(x, length))

#########
#set Condylarthra to have that as order designation
##could use occs to compile occurrences ofr order then append the family level for condylarthra and then reduce it to uniques.  That should work.....
occs.backup <- occs

#occs[occs$order %in% focal.order | occs$family %in% condylarth.families, ]
condylarth.families <- c("Arctocyonidae", "Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae")
focal.order <- c("Artiodactyla", "Perissodactyla")
#########
if(!load.occ.box) {
  splat <- function (this.rep) {
   # occ.box <- cbind(n_occs_ung=sapply(this.rep, function(this.intv) sum(occs[occs$occurrence_no %in% this.intv,]$order %in% focal.order)), 
    occ.box <- cbind(n_occs_ung=sapply(this.rep, function(this.intv) sum(occs[occs$occurrence_no %in% this.intv,]$order %in% focal.order | occs[occs$occurrence_no %in% this.intv,]$family %in% condylarth.families)), 
                     n_occs_carn=sapply(this.rep, function(this.intv) sum(occs[occs$occurrence_no %in% this.intv,]$order %in% c("Carnivora", "Creodonta"))), 
                     n_occs_rod=sapply(this.rep, function(this.intv) sum(occs[occs$occurrence_no %in% this.intv,]$order %in% c("Rodentia", "Lagomorpha"))), 
                     n_occs_other=sapply(this.rep, function(this.intv) sum(!occs[occs$occurrence_no %in% this.intv,]$order %in% c("Carnivora", "Creodonta","Rodentia", "Lagomorpha", focal.order) | !occs[occs$occurrence_no %in% this.intv,]$family %in% condylarth.families)), 
                     n_occs=sapply(this.rep, length),
                     n_cols_ung=sapply(this.rep, function(this.intv) length(sort(unique(occs$collection_no[occs$occurrence_no %in% this.intv & occs$order %in% focal.order])))),
                     n_cols_carn=sapply(this.rep, function(this.intv) length(sort(unique(occs$collection_no[occs$occurrence_no %in% this.intv & occs$order %in% c("Carnivora", "Creodonta")])))),
                     n_cols_rod=sapply(this.rep, function(this.intv) length(sort(unique(occs$collection_no[occs$occurrence_no %in% this.intv & occs$order %in% c("Rodentia", "Lagomorpha")])))),
                     n_cols_other=sapply(this.rep, function(this.intv) length(sort(unique(!occs$collection_no[occs$occurrence_no %in% this.intv & occs$order %in% c("Carnivora", "Creodonta","Rodentia", "Lagomorpha", focal.order) | !occs$occurrence_no %in% this.intv & occs$family %in% condylarth.families])))), 
                     n_cols=sapply(this.rep, function(this.intv) length(sort(unique(occs$collection_no[occs$occurrence_no %in% this.intv])))))
    occ.box <- cbind(occ.box, prop_occ_ung=round(occ.box[,"n_occs_ung"]/occ.box[,"n_occs"], digits=3), prop_occ_carn=round(occ.box[,"n_occs_carn"]/occ.box[,"n_occs"], digits=3), prop_occ_rod=round(occ.box[,"n_occs_rod"]/occ.box[,"n_occs"], digits=3), prop_occ_other=round(occ.box[,"n_occs_other"]/occ.box[,"n_occs"], digits=3))
    occ.box <- cbind(occ.box, prop_col_ung=round(occ.box[,"n_cols_ung"]/occ.box[,"n_cols"], digits=3), prop_col_carn=round(occ.box[,"n_cols_carn"]/occ.box[,"n_cols"], digits=3), prop_col_rod=round(occ.box[,"n_cols_rod"]/occ.box[,"n_cols"], digits=3), prop_col_other=round(occ.box[,"n_cols_other"]/occ.box[,"n_cols"], digits=3))
  }
  occ.box <- apply(sapply(repIntOccs, splat, simplify="array"), c(1,2), quantile, probs=c(0.025, 0.50, 0.975))
  
  save(occs, occ.box, file = paste0("~/Dropbox/Proposal/Chapter2/OccsBox_do.subample=",do.subsample,"_2022_4_22.Rdata"))
}

load("~/Dropbox/Proposal/Chapter2/OccsBox_2022_4_22.Rdata")
# occ.box <- data.frame(occ.box, prop_ung=round(occ.box[,"n_occs_ung"]/occ.box[,"n_occs"], digits=3), prop_carn=round(occ.box[,"n_occs_carn"]/occ.box[,"n_occs"], digits=3))

#quartz(height = 5.5, width = 8.5)
par(mfrow=c(2,1), mar=c(4,4,1.6,4), mgp=c(2, 1,0)) #, height = 2.15, width = 8.21)

plot(apply(intervals, 1, mean), 
     occ.box[2,,"n_occs_ung"], 
     xlim=c(65, 0), ylim=c(0,4000),  #ylim=c(0,750), 
     xaxp = c(70,0,14), yaxp = c(0, 4000, 8), #yaxp = c(0, 750, 15),
     type="n", 
     xlab="Time (Ma", ylab="Number of occurrences with taxon")
overlayCzTimescale(do.subepochs=TRUE)
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_occs_ung"], rev(occ.box[3,,"n_occs_ung"])), col=adjustcolor("dodgerblue1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_occs_carn"], rev(occ.box[3,,"n_occs_carn"])), col=adjustcolor("firebrick1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_occs_rod"], rev(occ.box[3,,"n_occs_rod"])), col=adjustcolor("darkolivegreen1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_occs"], rev(occ.box[3,,"n_occs"])), col=adjustcolor("grey25", alpha=0.2))
points(apply(intervals, 1, mean), occ.box[2,,"n_occs_rod"], pch=24, col="darkolivegreen4", type="o", bg="darkolivegreen1")
points(apply(intervals, 1, mean), occ.box[2,,"n_occs_carn"], pch=22, col="firebrick4", type="o", bg="firebrick1")
points(apply(intervals, 1, mean), occ.box[2,,"n_occs_ung"], pch=21, type="o", col="dodgerblue4", bg="dodgerblue1")
points(apply(intervals, 1, mean), occ.box[2,,"n_occs"], pch=21, type="o", col="black", bg="grey25")

plot(apply(intervals, 1, mean), 
     occ.box[2,,"prop_occ_ung"], 
     xlim=c(65, 0), ylim=c(0, 0.75), 
     xaxp = c(70,0,14), yaxp = c(0, 1, 10),
     type="n", 
     xlab="Time (Ma", ylab="Proportion of occurrences within guild")
overlayCzTimescale(do.subepochs=TRUE)
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_occ_ung"], rev(occ.box[3,,"prop_occ_ung"])), col=adjustcolor("dodgerblue1", alpha=0.2), border="dodgerblue1")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_occ_carn"], rev(occ.box[3,,"prop_occ_carn"])), col=adjustcolor("firebrick1", alpha=0.2), border="firebrick1")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_occ_other"], rev(occ.box[3,,"prop_occ_other"])), col=adjustcolor("darkgoldenrod1", alpha=0.2), border="darkgoldenrod4")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_occ_rod"], rev(occ.box[3,,"prop_occ_rod"])), col=adjustcolor("darkolivegreen1", alpha=0.2), border="darkolivegreen1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_occ_rod"], pch=24, col="darkolivegreen4", type="o", bg="darkolivegreen1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_occ_carn"], pch=22, col="firebrick4", type="o", bg="firebrick1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_occ_ung"], pch=21, type="o", col="dodgerblue4", bg="dodgerblue1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_occ_other"], pch=21, type="o", col="darkgoldenrod4", bg="darkgoldenrod1")

plot(apply(intervals, 1, mean), occ.box[2,,"prop_col_ung"], xlim=c(65, 0), ylim=c(0, 1), type="n", xlab="Time (Ma", ylab="Proportion of collections with taxon")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_col_ung"], rev(occ.box[3,,"prop_col_ung"])), col=adjustcolor("dodgerblue1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_col_carn"], rev(occ.box[3,,"prop_col_carn"])), col=adjustcolor("firebrick1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"prop_col_rod"], rev(occ.box[3,,"prop_col_rod"])), col=adjustcolor("darkolivegreen1", alpha=0.2))
points(apply(intervals, 1, mean), occ.box[2,,"prop_col_rod"], pch=24, col="darkolivegreen4", type="o", bg="darkolivegreen1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_col_carn"], pch=22, col="firebrick4", type="o", bg="firebrick1")
points(apply(intervals, 1, mean), occ.box[2,,"prop_col_ung"], pch=21, type="o", col="dodgerblue4", bg="dodgerblue1")

plot(apply(intervals, 1, mean), occ.box[2,,"n_cols_ung"], xlim=c(65, 0), ylim=c(0,65), type="n", xlab="Time (Ma", ylab="Number of collections with taxon")
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_cols_ung"], rev(occ.box[3,,"n_cols_ung"])), col=adjustcolor("dodgerblue1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_cols_carn"], rev(occ.box[3,,"n_cols_carn"])), col=adjustcolor("firebrick1", alpha=0.2))
polygon(x=c(apply(intervals, 1, mean), rev(apply(intervals, 1, mean))), c(occ.box[1,,"n_cols_rod"], rev(occ.box[3,,"n_cols_rod"])), col=adjustcolor("darkolivegreen1", alpha=0.2))
points(apply(intervals, 1, mean), occ.box[2,,"n_cols_rod"], pch=24, col="darkolivegreen4", type="o", bg="darkolivegreen1")
points(apply(intervals, 1, mean), occ.box[2,,"n_cols_carn"], pch=22, col="firebrick4", type="o", bg="firebrick1")
points(apply(intervals, 1, mean), occ.box[2,,"n_cols_ung"], pch=21, type="o", col="dodgerblue4", bg="dodgerblue1")

########################################################################################################################
run.taxon <-  "ungulates"

herbivore.size.cat <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
predator.size.cat  <- c(-Inf, 0, 0.845098, 1.322219, 2, Inf)

if(run.taxon == "ungulates") bmBreaks <- herbivore.size.cat
if(run.taxon == "carnivores") bmBreaks <- predator.size.cat


#make a new repIntTaxa that does not collapse things down to a single species entry

getAbundTaxaFromOneRepIntOccs <- function(this.repIntOccs, this.rank="species", do.rangethrough=TRUE) {
  if (this.rank=="species") {
    intTaxa <- lapply(this.repIntOccs, function(x) sort(as.character(occs$accepted_name[occs$occurrence_no %in% x & occs$accepted_rank=="species"])))
  } else if (this.rank=="genus") {
    intTaxa <- lapply(this.repIntOccs, function(x) sort(as.character(occs$genus[occs$occurrence_no %in% x & occs$accepted_rank %in% c("genus", "species")])))
  }
  
  if (do.rangethrough) intTaxa <- makeRangeThroughOneRep(intTaxa)
  intTaxa
}

getRepIntTaxaAbundFromRepIntOccs <- function(repIntOccs, this.rank="species", do.rangethrough=TRUE) {
                                        lapply(repIntOccs, getAbundTaxaFromOneRepIntOccs, this.rank=this.rank, 
                                               do.rangethrough= do.rangethrough)
                                }

repIntTaxa <- getRepIntTaxaAbundFromRepIntOccs(repIntOccs, this.rank=this.rank, do.rangethrough=do.rangethrough)

countCube <- sapply(repIntTaxa, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(measure.mat$bodyMass[match(this.intv, measure.mat$taxon)], breaks=bmBreaks, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

sizecateg <- c("rabbit (<5 kg)", "water deer (5-25 kg)",#"dog (5-25 kg)", could also do beaver for 5-25kg
               "antelope (25-150kg)", "horse (150-500 kg)", 
               "rhino (500-1000 kg)")

dimnames(countCube) <- list(sizecateg, rownames(intervals), NULL)


#countCube <- countCube[,,1]
#prop <- t(countCube)
herb_prop <- t(apply(countCube, c(1,2), mean, na.rm=TRUE))

colnames(prop)[colnames(herb_prop)==""] <- "indeterminate"
# dimnames(prop) <- list(rownames(intervals), shortFam)
plotStackedRichness(this.box=herb_prop, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, overlay.labels=FALSE, numbers.only=FALSE, legend=FALSE, xlim = c(65, #xlim=c(max(intervals, na.rm=TRUE)
                                                                                                                                                           min(intervals, na.rm=TRUE)))

countCube <- sapply(repIntTaxa, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(measure.mat$bodyMass[match(this.intv, measure.mat$taxon)], breaks=bmBreaks, plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

#countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

sizecateg <- c("rabbit (<5 kg)", "water deer (5-25 kg)",#"dog (5-25 kg)", could also do beaver for 5-25kg
               "antelope (25-150kg)", "horse (150-500 kg)", 
               "rhino (500-1000 kg)")

dimnames(countCube) <- list(sizecateg, rownames(intervals), NULL)


#countCube <- countCube[,,1]
#prop <- t(countCube)
prop <- t(apply(countCube, c(1,2), median, na.rm=TRUE))
colnames(prop)[colnames(prop)==""] <- "indeterminate"

par(mfrow=c(2,1), mar=c(4,4,1.6,4), mgp=c(2, 1,0)) #, height = 2.15, width = 8.21)
plotStackedRichness(this.box=prop, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, overlay.labels=FALSE, numbers.only=TRUE, legend=FALSE, xlim = c(65, #xlim=c(max(intervals, na.rm=TRUE)
                                                                                                                                                            min(intervals, na.rm=TRUE)))


pred_prop <- t(pred.preyDivCompare[2:6,])
colnames(prop)[colnames(prop)==""] <- "indeterminate"

#dimnames(prop) <- list(rownames(intervals), shortFam)
plotStackedRichness(this.box=pred_prop, intervals=intervals, reorder.taxa = FALSE, do.log=FALSE, overlay.labels=FALSE, numbers.only=TRUE, legend=FALSE, xlim = c(65, #xlim=c(max(intervals, na.rm=TRUE)
                                                                                                                                                                 min(intervals, na.rm=TRUE)))
