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

#######################################################################################################################################################
#PBDB generic ungulate diversity
ung.diversity <- read.csv("/Users/emdoughty/Downloads/Condylarth Papers/Archaic Ungualte vs Modern ungualte Generic Diversity.csv")

#plot(NULL, xlim=c(94,0), ylim = c(0,max(ung.diversity[,8:16], na.rm = TRUE)),
 #    xlab = "Time (Ma)", ylab = "Generic Diversity")

plot(NULL, xlim=c(94,0), ylim = c(0,80),
     xlab = "Time (Ma)", ylab = "Generic Diversity",
     xaxp = c(95,0,19), yaxp = c(0,80,16))
lines(ung.diversity$max_ma, ung.diversity$Terrestrial.Artio, col = "darkgreen")
lines(ung.diversity$max_ma, ung.diversity$Perissodactyla, col = "blue")
lines(ung.diversity$max_ma, ung.diversity$Acreodi, col = "red")
lines(ung.diversity$max_ma, ung.diversity$Ambylopoda, col = "orange")
lines(ung.diversity$max_ma, ung.diversity$Bulbulodentata, col = "purple")
lines(ung.diversity$max_ma, ung.diversity$Condylartha, col = "cyan")
lines(ung.diversity$max_ma, ung.diversity$Taeniodonta, col = "magenta")

abline(v=56)

legend("topleft", legend = c("Terrestrial Artiodactyla", "Perissodactyla", "Acreodi","Ambylopoda", "Bulbulodentata",
                             "Condylarthra", "Taeniodonta"),
       col = c("darkgreen", "blue", "red", "orange","purple","cyan","magenta"), lty = 1)


###########################################################################################################################

#get plot of species diversity form PBDB
occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=100&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", 
                                 "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", 
                                 "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", 
                                 "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

this.rank <- "species"
if (this.rank=="genus") measure.mat <- makeOneGenusMatFromSpecimenMat(measure.mat)

####################################################################################################################################
#### reduces matrix to just the focal order(s)
####################################################################################################################################

# focal.order <- "Artiodactyla"
# focal.order <- "Perissodactyla"
focal.order <- c("Artiodactyla", "Perissodactyla")
focal.family <- unique(occs[occs$order %in% focal.order,]$family)
focal.family <- c(as.character(focal.family),"Arctocyonidae", "Hyopsodontidae","Periptychidae","Phenacodontidae")
focal.family <- focal.family[!focal.family %in% ""]
focal.family <- focal.family[order(focal.family)]

bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% focal.family), c("order","family", "genus", "accepted_name")])
bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFam <- sort(unique(bigList$family[bigList$family %in% focal.family]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)

####################################################################################################################################
this.rank <- "species"
if (this.rank=="genus") measure.mat <- makeOneGenusMatFromSpecimenMat(measure.mat)

intervals <- getNALMAIntervals()

nalma.mark <- read.csv("/Users/emdoughty/Downloads/PPP_summary.csv")
nalma.mark <- nalma.mark[,1:3]
nalma.add  <- rbind(c("Aquilian", 84, 70),c("Lancian", 70, 66),c("Puercan", 66, 64.81)); colnames(nalma.add) <- colnames(nalma.mark)
nalma.mark <- rbind(nalma.mark, nalma.add)
nalma.mark[,2] <- as.numeric(nalma.mark[,2])
nalma.mark[,3] <- as.numeric(nalma.mark[,3])
nalma.mark <- nalma.mark[order(as.numeric(nalma.mark$Max_age), decreasing = TRUE),]
rownames(nalma.mark) <- nalma.mark$NALMA_Subdivision
colnames(nalma.mark) <- c("NALMA_Subdivision", "ageBase","ageTop")
nalma.mark <- nalma.mark[,-1]
nalma.mark <- nalma.mark[,c(2,1)]

intervals <- nalma.mark

rep.species <- getRepIntData(#bigList, 
  #col.nam = NULL,
  occs=occs,
  intervals = intervals,
  int_length = 2,
  int_min = 1,
  int_max = 95,
  do.parallel = FALSE, 
  reps = 1,
  this.rank = "species",
  do.rangethrough = TRUE)

#intervals <- makeIntervals(1,95,2)
midpoint <- array(data = rowMeans(intervals), dim = c(1,nrow(intervals),1))

#do thing order by order so not getting contaminants due to the wierd naming/data outputs of PBDB
#EX
#1)get family diversity of artios
#2) Set NO_FAMILY_SPECIFIED to be artio only

#Artiodactyla
artio.bigList <- getBigList(focal.order = c("Artiodactyla"))
artio.shortFam <- sort(unique(artio.bigList$family))	
artio.taxCube <- sapply(rep.species$RepsTaxa, function(y) sapply(y, function(x) tabulate(match(artio.bigList$family[as.character(artio.bigList$accepted_name) %in% x], artio.shortFam), nbins=length(artio.shortFam)), simplify="array"), simplify="array")
dimnames(artio.taxCube) <- list(artio.shortFam, rownames(intervals), NULL)

artio.total <- colSums(artio.taxCube); colnames(artio.total) <- "Artiodactyla"
artio.total.mat <- t(artio.total)

artio.total.array <- array(data = colMeans(artio.total.mat), dim = c(1,nrow(intervals),1))

artio.all <- abind(artio.taxCube, Artio.Sum = artio.total.array, midpoint=midpoint, along = 1)


#Perissodactyla
perisso.bigList <- getBigList(focal.order = c("Perissodactyla"))
perisso.shortFam <- sort(unique(perisso.bigList$family))	
perisso.taxCube <- sapply(rep.species$RepsTaxa, function(y) sapply(y, function(x) tabulate(match(perisso.bigList$family[as.character(perisso.bigList$accepted_name) %in% x], perisso.shortFam), nbins=length(perisso.shortFam)), simplify="array"), simplify="array")
dimnames(perisso.taxCube) <- list(perisso.shortFam, rownames(intervals), NULL)

perisso.total <- colSums(perisso.taxCube); colnames(perisso.total) <- "Perissodactyla"
perisso.total.mat <- t(perisso.total)

perisso.total.array <- array(data = colMeans(perisso.total.mat), dim = c(1,nrow(intervals),1))

perisso.all <- abind(perisso.taxCube, Perisso.Sum = perisso.total.array, midpoint=midpoint, along = 1)

#Acreodi
acreo.bigList <- getBigList(focal.order = c("Acreodi"))
acreo.shortFam <- sort(unique(acreo.bigList$family))	
acreo.taxCube <- sapply(rep.species$RepsTaxa, function(y) sapply(y, function(x) tabulate(match(acreo.bigList$family[as.character(acreo.bigList$accepted_name) %in% x], acreo.shortFam), nbins=length(acreo.shortFam)), simplify="array"), simplify="array")
dimnames(acreo.taxCube) <- list(acreo.shortFam, rownames(intervals), NULL)

acreo.total <- colSums(acreo.taxCube); colnames(acreo.total) <- "Acreodi"
acreo.total.mat <- t(acreo.total)

acreo.total.array <- array(data = colMeans(acreo.total.mat), dim = c(1,nrow(intervals),1))

acreo.all <- abind(acreo.taxCube, acreo.Sum = acreo.total.array, midpoint=midpoint, along = 1)

#Dinocerata
dino.bigList <- getBigList(focal.order = c("Dinocerata"))
dino.shortFam <- sort(unique(dino.bigList$family))	
dino.taxCube <- sapply(rep.species$RepsTaxa, function(y) sapply(y, function(x) tabulate(match(dino.bigList$family[as.character(dino.bigList$accepted_name) %in% x], dino.shortFam), nbins=length(dino.shortFam)), simplify="array"), simplify="array")
dimnames(dino.taxCube) <- list(dino.shortFam, rownames(intervals), NULL)

dino.total <- colSums(dino.taxCube); colnames(dino.total) <- "Acreodi"
dino.total.mat <- t(dino.total)

dino.total.array <- array(data = colMeans(dino.total.mat), dim = c(1,nrow(intervals),1))

dino.all <- abind(dino.taxCube, dino.Sum = dino.total.array, midpoint=midpoint, along = 1)

#"Condylarthra"
##Condylarthra
cond.bigList <- getBigList(focal.order = c("Condylarthra"))
cond.shortFam <- sort(unique(cond.bigList$family))	
cond.taxCube <- sapply(rep.species$RepsTaxa, function(y) sapply(y, function(x) tabulate(match(cond.bigList$family[as.character(cond.bigList$accepted_name) %in% x], cond.shortFam), nbins=length(cond.shortFam)), simplify="array"), simplify="array")
dimnames(cond.taxCube) <- list(cond.shortFam, rownames(intervals), NULL)

##"Procreodi"
proc.bigList <- getBigList(focal.order = c("Procreodi"))
proc.shortFam <- sort(unique(proc.bigList$family))	
proc.taxCube <- sapply(rep.species$RepsTaxa, function(y) sapply(y, function(x) tabulate(match(proc.bigList$family[as.character(proc.bigList$accepted_name) %in% x], proc.shortFam), nbins=length(proc.shortFam)), simplify="array"), simplify="array")
colnames(proc.taxCube) <- proc.shortFam
#dim(proc.taxCube) <- c(dim(proc.taxCube),1)
#colnames(proc.taxCube) <- rownames(intervals)
#rownames(proc.taxCube[1,48,1]) <- shortFam

##Phenacodontidae
phen.bigList <-  getBigList(focal.order = "NO_ORDER_SPECIFIED")
phen.shortFam <- "Phenacodontidae"
phen.taxCube <- sapply(rep.species$RepsTaxa, function(y) sapply(y, function(x) tabulate(match(phen.bigList$family[as.character(phen.bigList$accepted_name) %in% x], phen.shortFam), nbins=length(phen.shortFam)), simplify="array"), simplify="array")
colnames(phen.taxCube) <- phen.shortFam
#phen.taxCube <- t(phen.taxCube)

##Periptychidae
peri.bigList <-  getBigList(focal.order = "NO_ORDER_SPECIFIED")
peri.shortFam <- "Periptychidae"
peri.taxCube <- sapply(rep.species$RepsTaxa, function(y) sapply(y, function(x) tabulate(match(peri.bigList$family[as.character(peri.bigList$accepted_name) %in% x], peri.shortFam), nbins=length(peri.shortFam)), simplify="array"), simplify="array")
colnames(peri.taxCube) <- peri.shortFam
#peri.taxCube <- t(peri.taxCube)

##Combined
condy.combined <- cbind(proc.taxCube, phen.taxCube, peri.taxCube)
condy.combined <- t(condy.combined)
test <- array(condy.combined, dim = c(3,ncol(cond.taxCube),1))
dimnames(test)[[1]] <- rownames(condy.combined)
dimnames(test)[[2]] <- colnames(condy.combined)

condyl.all <- abind(test, cond.taxCube, along = 1)

condyl.total <- colSums(condyl.all); colnames(condyl.total) <- "Condylarthra"
condyl.total.mat <- t(condyl.total)

condyl.total.array <- array(data = condyl.total.mat, dim = c(1,nrow(intervals),1))

condyl.all <- abind(condyl.all, Condy.Sum = condyl.total.array, midpoint=midpoint, along = 1)


plot(NA, xlim = c(75,0), ylim = c(0,120), xaxp = c(0,95,19), xlab = "Time (Ma)", ylab = 'Species Richness',
     main = "Species Richness Through Cenozoic (Sub-NALMA)")
overlayCzTimescale(do.subepochs= TRUE)
lines(midpoint, artio.all["Artio.Sum",,], col = "darkgreen")
lines(midpoint, perisso.all["Perisso.Sum",,], col = "blue")
lines(midpoint, acreo.all["acreo.Sum",,], col = "red")
lines(midpoint, dino.all["dino.Sum",,], col = "orange")
#lines(midpoint, condyl.all["Condy.Sum",,], col = "cyan")
lines(midpoint, condyl.all["Arctocyonidae",,], col = "purple")
lines(midpoint, condyl.all["Phenacodontidae",,], col = "magenta")
lines(midpoint, condyl.all["Periptychidae",,], col = "black")
lines(midpoint, condyl.all["Hyopsodontidae",,], col = "green")

abline(v=56)

legend("topleft", legend = c("Terrestrial Artiodactyla", "Perissodactyla", "Acreodi","Ambylopoda","Condylarthra",
                             "    Arctocyonidae", "    Phenacodontidae", "    Periptychidae", "    Hyopsodontidae"),
       col = c("darkgreen", "blue", "red", "orange","white","purple", "magenta","black","green"), lty = 1)


ung.all <- rbind(artio.all["Artio.Sum",,], 
                 perisso.all["Perisso.Sum",,],
                 acreo.all["acreo.Sum",,],
                 dino.all["dino.Sum",,],
                 condyl.all["Condy.Sum",,])
rownames(ung.all) <- c("Artio.Sum", "Perisso.Sum", "Acreo.Sum", "Dino.Sum", "Condy.Sum")

ung.all <- rbind(ung.all, Total = colSums(ung.all))

###############################################################################################################
# Look at ungulate diversity for times where we see multiple predators over 21 kg
pred.data <- read.csv("~/Dropbox/Code/R/Carnivore Paleobio Seminar/Predator_data_all.csv")
colnames(pred.data) <- c("family", "taxon",	"max_ma","min_ma",	"m1L",	"rbl",	"BM_all_carnivoran",	"BM_extant_reg")
pred.data$taxon <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.data$taxon)
rownames(pred.data) <- pred.data$taxon

pred.data[,c("BM_all_carnivoran",	"BM_extant_reg")] <- log10(pred.data[,c("BM_all_carnivoran",	"BM_extant_reg")])

focal.orderPred <- c("Carnivora", "Creodonta","Hyaenodonta")
occsPred <- occs[occs$order %in% focal.orderPred,]
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

pred.data <- getDateTaxa(measure.mat = pred.data, occs=occsPred, this.rank = "species")
pred.data <- pred.data[complete.cases(pred.data),]

intervals
test.pred <- matrix(nrow = nrow(pred.data), ncol = nrow(intervals))
colnames(test.pred) <- rownames(intervals)
rownames(test.pred) <- rownames(pred.data)

#sort into bins using a AND(base date) or AND(top date) style filter
#may need to do rangethrough if using small bins.  i.e. 2 Ma bins may fall between FO and LO and not be filtered

#maybe could have set to 
###replace start and end dates with the interval name, 
###take unique from those columns
###index to fill those columns with a 1 for that species
###do rangethrough if needed

for(xx in seq(1, nrow(test.pred),1))
{
  FO <- pred.data[rownames(test.pred)[xx],"max_ma"]
  LO <- pred.data[rownames(test.pred)[xx],"min_ma"]
  
  for(yy in seq(1, ncol(test.pred),1))
  {
    
    int.top <- intervals[rownames(intervals) %in% colnames(test.pred)[yy],]$ageTop
    int.base <- intervals[rownames(intervals) %in% colnames(test.pred)[yy],]$ageBase
    if(FO < int.base & FO >= int.top) test.pred[xx,yy] <- 1
    if(LO < int.base & LO >= int.top) test.pred[xx,yy] <- 1
  }
  #do rangethrough
  col.index <- which(test.pred[xx,] >= 1)
  if(length(col.index) > 1) test.pred[xx,seq(col.index[1], col.index[2])] <- 1
}
#write.csv(test.pred, file = "/Users/emdoughty/Dropbox/Code/pred.diversitycheck.csv")

carniv.Diversity <- colSums(test.pred, na.rm = TRUE)

#lines(midpoint,carniv.Diversity)

pred.21up <- pred.data[pred.data$BM_all_carnivoran > log10(21),]
predDiv.21up <- colSums(test.pred[rownames(test.pred) %in% rownames(pred.21up),], na.rm = TRUE)

pred.100up <- pred.data[pred.data$BM_all_carnivoran > log10(100),]
predDiv.100up <- colSums(test.pred[rownames(test.pred) %in% rownames(pred.100up),], na.rm = TRUE)

pred.200up <- pred.data[pred.data$BM_all_carnivoran > log10(200),]
predDiv.200up <- colSums(test.pred[rownames(test.pred) %in% rownames(pred.200up),], na.rm = TRUE)

pred.preyDivCompare <- rbind(AllPredDiversity = carniv.Diversity, 
                             ">21kg" = predDiv.21up,
                             ">100kg" = predDiv.100up,
                             ">200kg" = predDiv.200up, 
                             UngulateDiv = ung.all["Total",])

# What are the distributions of ungulate bodymass during these times?
####what is the diversity of the diffrerent size catagories of prey?

measure.mat <- getMeasureMatWithBodyMasses()
measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name[bigList$family %in% focal.family], ]

#make sure condylarthra and archics are appended on

countCube <- sapply(rep.species$RepsTaxa, function(this.rep) {
  sapply(this.rep, function(this.intv, this.rep) {
    hist(measure.mat[,"bodyMass"][match(this.intv, measure.mat$taxon)], 
         breaks=c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, Inf), plot=FALSE)$counts
  }, this.rep=this.rep)
}, simplify = "array")

rownames(countCube) <- c("rabbit (<5 kg)","dog (5-25 kg)",
                         "antelope (25-150kg)", "horse (150-500 kg)", 
                         "rhino (500-1000 kg)","megaherbivore (>1000 kg)")

#######################################################################################
#plot total diversity curves
plot(NA, xlim = c(80,0), ylim = c(0,210), xaxp = c(80,0,16), yaxp = c(0,220,11),
     xlab = "Time (Ma)",ylab = "Species Diversity")
overlayCzTimescale(do.subepochs= TRUE)
lines(midpoint, pred.preyDivCompare["UngulateDiv",], col = "darkgreen")
lines(midpoint, pred.preyDivCompare["AllPredDiversity",], col = "red")
lines(midpoint, pred.preyDivCompare[">21kg",], col = "cyan")
lines(midpoint, pred.preyDivCompare[">100kg",], col = "lightblue")
lines(midpoint, pred.preyDivCompare[">200kg",], col = "darkblue")

legend("topleft", legend = c("All Ungulate", "All Predators", "Pred. >21kg","Pred. >100kg","Pred. >200kg"),
       col = c("darkgreen", "red", "cyan","lightblue","darkblue"), lty = 1)

#plot diversity of ungualte size catagories vs carnivoran diversity
plot(NA, xlim = c(80,0), ylim = c(0,200), xaxp = c(80,0,16), yaxp = c(0,200,8),
     xlab = "Time (Ma)",ylab = "Species Diversity")
overlayCzTimescale(do.subepochs= TRUE)
lines(midpoint, countCube["rabbit (<5 kg)",,], col = "darkgreen")
lines(midpoint, countCube["dog (5-25 kg)",,], col = "darkgreen")
lines(midpoint, countCube["antelope (25-150kg)",,], col = "darkgreen")
lines(midpoint, countCube["horse (150-500 kg)",,], col = "darkgreen")
lines(midpoint, countCube["rhino (500-1000 kg)",,], col = "darkgreen")
lines(midpoint, countCube["megaherbivore (>1000 kg)",,], col = "darkgreen")

lines(midpoint, pred.preyDivCompare["AllPredDiversity",], col = "red")
lines(midpoint, pred.preyDivCompare[">21kg",], col = "cyan")
lines(midpoint, pred.preyDivCompare[">100kg",], col = "lightblue")
lines(midpoint, pred.preyDivCompare[">200kg",], col = "darkblue")

####################################################################################
#plot carnivore vs ungulate diversity
plot(NA, xlim = c(0,200), ylim = c(0,100), xaxp = c(0,200,20), yaxp = c(0,100,10),
     xlab = "Ungulate Species Diversity",ylab = "Predator Species Diversity",
     main = "Ungulate vs Predator Diversity")
#lines(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare["AllPredDiversity",], lty = 1)
points(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare["AllPredDiversity",])
text(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare["AllPredDiversity",],
     labels = colnames(pred.preyDivCompare), cex = 0.5, adj = 1, pos = 4)

#plot carnivore vs ungulate diversity
plot(NA, xlim = c(0,210), ylim = c(0,40), xaxp = c(0,200,20), yaxp = c(0,40,8),
     xlab = "Ungulate Species Diversity",ylab = "Predator Species Diversity > 21kg",
     main = "Ungulate vs >21kg Predator Diversity")
#lines(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">21kg",], lty = 1)
points(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">21kg",])
text(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">21kg",],
     labels = colnames(pred.preyDivCompare), cex = 0.5, adj = 1, pos = 4)

plot(NA, xlim = c(0,210), ylim = c(0,15), xaxp = c(0,200,20), yaxp = c(0,40,8),
     xlab = "Ungulate Species Diversity",ylab = "Predator Species Diversity > 100kg",
     main = "Ungulate vs >100kg Predator Diversity")
#lines(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">100kg",], lty = 1)
points(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">100kg",])
text(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">100kg",],
     labels = colnames(pred.preyDivCompare), cex = 0.5, adj = 1, pos = 4)

plot(NA, xlim = c(0,210), ylim = c(0,5), xaxp = c(0,200,20), yaxp = c(0,40,40),
     xlab = "Ungulate Species Diversity",ylab = "Predator Species Diversity > 200kg",
     main = "Ungulate vs >200kg Predator Diversity")
#lines(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">200kg",], lty = 1)
points(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">200kg",])
text(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare[">200kg",],
     labels = colnames(pred.preyDivCompare), cex = 0.5, adj = 1, pos = 4)

############################################################################
#Carnivore vs Ungualte diversity for dif size catagories of Ungulate
#make some kind of grid pattern for automation of plots
##i.e. 6x4 layout of plots
###Cols the carnivore breeakdown, rows of plots are the ungulate catagories

colfunc <- colorRampPalette(colors = c("red","yellow","springgreen","royalblue"))
col.nalma <- colfunc(ncol(pred.preyDivCompare))

do.legend = TRUE

quartz(h=20,w=11)
par(mfrow = c(7,4), oma = c(1,4,1,1))

for(xx in seq(1, nrow(countCube),1))
{
  for(yy in seq(1, nrow(pred.preyDivCompare)-1,1))
  {
    plot(NA, 
         xlim = c(0,round_any(max(countCube[xx,,]),10,f = ceiling)), 
         ylim = c(0,round_any(max(pred.preyDivCompare[yy,]),10,f = ceiling)),
         xaxp = c(0,200,40), yaxp = c(0,100,20),
         xlab = "", ylab = "", cex.lab = 0.75)
    title(ylab = "Predator Diversity", line = 2)
    title(xlab = "Ungulate Diversity", line = 2)
    if(yy == 1) mtext(rownames(countCube)[xx], side = 2, line = 4, cex = 1.5)
    if(xx == 1) mtext(rownames(pred.preyDivCompare)[yy], side = 3, line = 2, cex = 1.5)
    #f(do.legend & xx == 1 & yy == 1) legend("topright", legend = c(colnames(pred.preyDivCompare)), pch = 16, col = col.nalma, cex = 0.33)

    #lines(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare["AllPredDiversity",], lty = 1)
    text(countCube[xx,,], pred.preyDivCompare[yy,],
         labels = colnames(pred.preyDivCompare), 
         cex = 0.25, adj = 1, pos = 4, 
         col = ifelse(countCube[xx,,] < 1,alphaColor("black",0.35),'black'))
    points(countCube[xx,,], pred.preyDivCompare[yy,], pch=16, col = ifelse(countCube[xx,,] < 1 | pred.preyDivCompare[yy,] < 1,alphaColor(col.nalma,0.35),col.nalma))
     }
}
plot(NA, 
     xlim = c(0,round_any(max(countCube[xx,,]),10,f = ceiling)), 
     ylim = c(0,round_any(max(pred.preyDivCompare[yy,]),10,f = ceiling)))
legend("topleft", legend = c(colnames(pred.preyDivCompare)), pch = 16, col = col.nalma, cex = 0.40, ncol = 3)

#FOR PRESENTATION
quartz(h=16,w=16)
par(mfrow = c(5,6), oma = c(1,4,1,1), mar = c(3,3,3,3))

for(xx in seq(1, nrow(pred.preyDivCompare)-1,1))
{
  for(yy in  seq(1, nrow(countCube),1))
  {
    plot(NA, 
         xlim = c(0,round_any(max(countCube[yy,,]),10,f = ceiling)), 
         ylim = c(0,round_any(max(pred.preyDivCompare[xx,]),10,f = ceiling)),
         xaxp = c(0,200,40), yaxp = c(0,100,20),
         xlab = "", ylab = "", cex.lab = 0.75)
    title(ylab = "Predator Diversity", line = 2)
    title(xlab = "Ungulate Diversity", line = 2)
    if(xx == 1) mtext(rownames(countCube)[yy], side = 3, line = 2, cex = 1.5)
    if(yy == 1) mtext(rownames(pred.preyDivCompare)[xx], side = 2, line = 3.5, cex = 1.5)
    #f(do.legend & xx == 1 & yy == 1) legend("topright", legend = c(colnames(pred.preyDivCompare)), pch = 16, col = col.nalma, cex = 0.33)
    
    #lines(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare["AllPredDiversity",], lty = 1)
    text(countCube[yy,,], pred.preyDivCompare[xx,],
         labels = colnames(pred.preyDivCompare), 
         cex = 0.25, adj = 1, pos = 4, 
         col = ifelse(countCube[yy,,] < 1,alphaColor("black",0.35),'black'))
    points(countCube[yy,,], pred.preyDivCompare[xx,], pch=16, col = ifelse(countCube[yy,,] < 1 | pred.preyDivCompare[xx,] < 1,alphaColor(col.nalma,0.35),col.nalma))
  }
}
plot(NA, 
     xlim = c(0,round_any(max(countCube[yy,,]),10,f = ceiling)), 
     ylim = c(0,round_any(max(pred.preyDivCompare[xx,]),10,f = ceiling)))
legend("topleft", legend = c(colnames(pred.preyDivCompare)), pch = 16, col = col.nalma, cex = 0.40, ncol = 3)


##############################################################################################################################
#First Differences for Diversity
pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = pred.preyDivCompare, intervals = intervals)

UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = countCube[,,1], intervals = intervals)


quartz(h=20,w=11)
par(mfrow = c(8,4), oma = c(1,4,1,1))

for(xx in seq(1, nrow(UngualteBMGroupDiv_FirstDiff),1))
{
  for(yy in seq(1, nrow(pred.prey.DivFirstDiff)-1,1))
  {
    plot(NA, 
         xlim = c(round_any(min(UngualteBMGroupDiv_FirstDiff[xx,]),10,f = floor),
                  round_any(max(UngualteBMGroupDiv_FirstDiff[xx,]),10,f = ceiling)), 
         ylim = c(round_any(min(pred.prey.DivFirstDiff[yy,]),10,f = floor),
                  round_any(max(pred.prey.DivFirstDiff[yy,]),10,f = ceiling)),
         xaxp = c(-200,200,80), yaxp = c(-100,100,40),
         xlab = "", ylab = "", cex.lab = 0.75)
    title(ylab = "Predator Diversity 1st Diff", line = 2)
    title(xlab = "Ungulate Diversity 1st Diff", line = 2)
    if(yy == 1) mtext(rownames(UngualteBMGroupDiv_FirstDiff)[xx], side = 2, line = 4, cex = 1.5)
    if(xx == 1) mtext(rownames(pred.prey.DivFirstDiff)[yy], side = 3, line = 2, cex = 1.5)
    #f(do.legend & xx == 1 & yy == 1) legend("topright", legend = c(colnames(pred.preyDivCompare)), pch = 16, col = col.nalma, cex = 0.33)
    
    #lines(pred.preyDivCompare["UngulateDiv",], pred.preyDivCompare["AllPredDiversity",], lty = 1)
    # text(UngualteBMGroupDiv_FirstDiff[xx,], pred.prey.DivFirstDiff[yy,],
    #       labels = colnames(pred.prey.DivFirstDiff), 
    #      cex = 0.25, adj = 1, pos = 4, 
    #     col = 'black')
    points(UngualteBMGroupDiv_FirstDiff[xx,], pred.prey.DivFirstDiff[yy,], pch = 16, 
           col = col.nalma)
    points(UngualteBMGroupDiv_FirstDiff[xx,], pred.prey.DivFirstDiff[yy,], pch = 1, lwd = 0.25,
           col = "black")
    
    best.fit <- lm(pred.prey.DivFirstDiff[yy,]~UngualteBMGroupDiv_FirstDiff[xx,])
    abline(best.fit)
    mtext(text = paste("y =", summary(best.fit)$coefficients[2], "*x +",summary(best.fit)$coefficients[1], 
                       sep= " "),
          side = 3, line = -0.75, cex = 0.25, adj = 0.1)
    mtext(text = paste("R^2=",summary(best.fit)$r.squared,
                       sep = " "),
          side = 3, line = -1.25, cex = 0.25, adj = 0.1)
    mtext(text = paste("adj R^2=", summary(best.fit)$adj.r.squared,
                       sep = " "),
          side = 3, line = -1.75, cex = 0.25, adj = 0.1)
    mtext(text = paste("P-value=", summary(best.fit)$coefficients[8], sep = " "),
          side = 3, line = -2.25, cex = 0.25, adj = 0.1)
    
  }
}

for(yy in seq(1, nrow(pred.prey.DivFirstDiff)-1,1))
{
  plot(NA, 
       xlim = c(round_any(min(pred.prey.DivFirstDiff["UngulateDiv",]),10,f = floor),
                round_any(max(pred.prey.DivFirstDiff["UngulateDiv",]),10,f = ceiling)), 
       ylim = c(round_any(min(pred.prey.DivFirstDiff[yy,]),10,f = floor),
                round_any(max(pred.prey.DivFirstDiff[yy,]),10,f = ceiling)),
       xaxp = c(-200,200,80), yaxp = c(-100,100,40),
       xlab = "", ylab = "", cex.lab = 0.75)
  title(ylab = "Predator Diversity 1st Diff", line = 2)
  title(xlab = "Ungulate Diversity 1st Diff", line = 2)
  if(yy == 1) mtext("All Ungulate", side = 2, line = 4, cex = 1.5)
  
  #text(pred.prey.DivFirstDiff["UngulateDiv",], pred.prey.DivFirstDiff[yy,],
  #     labels = colnames(pred.prey.DivFirstDiff), 
  #     cex = 0.25, adj = 1, pos = 4, 
  #     col = 'black')
  points(pred.prey.DivFirstDiff["UngulateDiv",], pred.prey.DivFirstDiff[yy,], pch = 16, col = col.nalma)
  points(pred.prey.DivFirstDiff["UngulateDiv",], pred.prey.DivFirstDiff[yy,], pch = 1, lwd = 0.25, col = "black")
  
  best.fit <- lm(pred.prey.DivFirstDiff[yy,]~pred.prey.DivFirstDiff["UngulateDiv",])
  abline(best.fit)
  mtext(text = paste("y =", summary(best.fit)$coefficients[2], "*x +",summary(best.fit)$coefficients[1], 
                     sep= " "),
        side = 3, line = -0.75, cex = 0.25, adj = 0.1)
  mtext(text = paste("R^2=",summary(best.fit)$r.squared,
                     sep = " "),
        side = 3, line = -1.25, cex = 0.25, adj = 0.1)
  mtext(text = paste("adj R^2=", summary(best.fit)$adj.r.squared,
                     sep = " "),
        side = 3, line = -1.75, cex = 0.25, adj = 0.1)
  mtext(text = paste("P-value=", summary(best.fit)$coefficients[8], sep = " "),
        side = 3, line = -2.25, cex = 0.25, adj = 0.1)
}

plot(NA, 
     xlim = c(0,round_any(max(UngualteBMGroupDiv_FirstDiff[xx,]),10,f = ceiling)), 
     ylim = c(0,round_any(max(pred.prey.DivFirstDiff[yy,]),10,f = ceiling)))
legend("topleft", legend = c(colnames(pred.prey.DivFirstDiff)), pch = 16, col = col.nalma, cex = 0.8, ncol = 2)


#############################################################################################################################

#What taxon group comprise these catagories for Carnivore and Ungulate, both in general and for size catagories

#make funciton to plot NAlma subsets like Jons geo tiemscale overlay
#make function to overlay habiat shifts onto plots (e.g. opening habitat, grassland appearance, C3->c4 transition)

#shoulder plot grpahs for GSA 2021 ppt
#ungulates


load("/Users/emdoughty/Dropbox/Code/R/Results/BM_handleyResult_ungulates_SampleStandardized=TRUE_Reps=5000Jonathans_MacBook_Pro.local##------ Tue Sep  7 17:23:16 2021 ------##.Rdata")
repIntTaxa_ung <- repIntTaxa
repIntOccs_ung <- repIntOccs
countCube_ung <- countCube
optList_bm_median_ung <- optList_bm_median
optList_bm_allReps_ung <- optList_bm_allReps
bm_quants_ung <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)

load("/Users/emdoughty/Dropbox/Code/R/Results/Taxon_handleyResult_ungulates_SampleStandardized=TRUE_Reps=5000Jonathans_MacBook_Pro.local##------ Tue Sep  7 15:56:16 2021 ------##.Rdata")
taxCube_ung <- taxCube
optList_tax_median_ung <- optList_tax_median
optList_tax_allReps_ung <- optList_tax_allReps

focal.order <- c("Artiodactyla", "Perissodactyla")
focal.family <- unique(occs[occs$order %in% focal.order,]$family)
#search through those without order
add.family <- c("Arctocyonidae", "Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae")
focal.family <- c(as.character(focal.family), add.family)
focal.family <- focal.family[!focal.family %in% ""]
focal.family <- focal.family[order(focal.family)]

bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
bigList.cond <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% add.family), c("order","family", "genus", "accepted_name")])
bigList <- rbind(bigList,bigList.cond)

bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFam_ung <- sort(unique(bigList$family[bigList$family %in% focal.family]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
bigList_ung <- bigList

# Number of shifts per rep
quartz()
par(mfrow=c(2,1), mar=c(5, 4, 4, 1))
hist(sapply(optList_tax_allReps, function(x) length(x) - 2), breaks=seq(-0.5, 10, 1.0), col="orchid4", main="Number of taxonomic shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))
hist(sapply(optList_bm_allReps, function(x) length(x) - 2), breaks=seq(-0.5, 11.5, 1.0), col="firebrick4", main="Number of body mass shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))

quants <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
####################################################################################################################################
### number of Replicates with a shift in that interval
quartz()
par(mfrow=c(2,1), mar=c(4, 4, 1, 1))
# optList_tax_allReps <- optList_tax_allReps_heuristic
# optList_tax_allReps <- optList_tax_allReps_full
breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_tax_allReps, function(x) x[[length(x) - 1]]$optBreaks))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xlim=rev(range(intervals)), xaxp =c(65,5,10), ylim=c(0,1000), main="Number of Replicates with a Taxonomic Distribution Shift", xlab="Time (Ma)")
overlayCzTimescale(do.subepochs=TRUE)
plot(breakHist, col="orchid4", border="orchid1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(65,5,10),  xlim=c(55, 0), ylim=c(0,1000), add=TRUE)

# par(mfrow=c(2,1), mar=c(3.5, 3.5, 1, 1))
# optList_bm_allReps <- optList_bm_allReps_heuristic
# optList_bm_allReps <- optList_bm_allReps_full
### number of Replicates with a shift in that interval
breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_bm_allReps, function(x) unique(x[[length(x) - 1]]$optBreaks)))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,10), xlim=rev(range(intervals)), ylim=c(0,1000), main="Number of Replicates with a Body Mass Distribution Shift", xlab="Time (Ma)")
overlayCzTimescale(do.subepochs=TRUE)
plot(breakHist, col="firebrick4", border="firebrick1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,10), xlim=c(55, 0), ylim=c(0,1000), add=TRUE)

#########################
#carnivores
load("/Users/emdoughty/Dropbox/Code/R/Results/BM_handleyResult_carnivores_SampleStandardized=FALSE_Reps=1000Armageddon_2.local##------ Wed Sep  8 04:43:09 2021 ------##.Rdata")
repIntTaxa_pred <- repIntTaxa
repIntOccs_pred <- repIntOccs
countCube_pred <- countCube
optList_bm_median_pred <- optList_bm_median
optList_bm_allReps_pred <- optList_bm_allReps
bm_quants_pred <- apply(sapply(repIntTaxa_pred, function(y) sapply(y, function(x) quantile(pred.data[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)


load("/Users/emdoughty/Dropbox/Code/R/Results/Taxon_handleyResult_carnivores_SampleStandardized=FALSE_Reps=1000Armageddon_2.local##------ Wed Sep  8 03:30:45 2021 ------##.Rdata")
taxCube_pred <- taxCube
optList_tax_median_pred <- optList_tax_median
optList_tax_allReps_pred <- optList_tax_allReps

# Number of shifts per rep
quartz()
par(mfrow=c(2,1), mar=c(5, 4, 4, 1))
hist(sapply(optList_tax_allReps, function(x) length(x) - 2), breaks=seq(-0.5, 10, 1.0), col="orchid4", main="Number of taxonomic shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))
hist(sapply(optList_bm_allReps, function(x) length(x) - 2), breaks=seq(-0.5, 11.5, 1.0), col="firebrick4", main="Number of body mass shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))

quants <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
####################################################################################################################################
### number of Replicates with a shift in that interval
quartz()
par(mfrow=c(2,1), mar=c(4, 4, 1, 1))
# optList_tax_allReps <- optList_tax_allReps_heuristic
# optList_tax_allReps <- optList_tax_allReps_full
breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_tax_allReps, function(x) x[[length(x) - 1]]$optBreaks))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xlim=rev(range(intervals)), xaxp =c(65,5,10), ylim=c(0,1000), main="Number of Replicates with a Taxonomic Distribution Shift", xlab="Time (Ma)")
overlayCzTimescale(do.subepochs=TRUE)
plot(breakHist, col="orchid4", border="orchid1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(65,5,10),  xlim=c(55, 0), ylim=c(0,1000), add=TRUE)

# par(mfrow=c(2,1), mar=c(3.5, 3.5, 1, 1))
# optList_bm_allReps <- optList_bm_allReps_heuristic
# optList_bm_allReps <- optList_bm_allReps_full
### number of Replicates with a shift in that interval
breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_bm_allReps, function(x) unique(x[[length(x) - 1]]$optBreaks)))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,10), xlim=rev(range(intervals)), ylim=c(0,1000), main="Number of Replicates with a Body Mass Distribution Shift", xlab="Time (Ma)")
overlayCzTimescale(do.subepochs=TRUE)
plot(breakHist, col="firebrick4", border="firebrick1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,10), xlim=c(55, 0), ylim=c(0,1000), add=TRUE)
  
focal.order <- c("Carnivora", "Creodonta","Hyaenodonta")
#occsPred <- occs[occs$order %in% focal.orderPred,]
focal.family <- unique(occs[occs$order %in% focal.order,]$family)
add.family <- c("Viverravidae")

focal.family <- c(as.character(focal.family), add.family)
focal.family <- focal.family[!focal.family %in% ""]
focal.family <- focal.family[order(focal.family)]

bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
bigList.cond <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% add.family), c("order","family", "genus", "accepted_name")])
bigList <- rbind(bigList,bigList.cond)

bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFam_pred <- sort(unique(bigList$family[bigList$family %in% focal.family]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
bigList_pred <- bigList

#bm shifts ungulates c(57, 51, 39, 21, 5)
par(mfrow=c(2,1), mar=c(3.9,4,0.5,0.5))#, mgp=c(2, 1,0))
shoulderPlot(measure.mat = measure.mat, plot.y = "bodyMass", occs = occs, this.rank = "species", 
             bigList = bigList_ung, shortFam = shortFam_ung, repIntTaxa = repIntTaxa_ung, quants = bm_quants_ung,
             intervals = makeIntervals(1,64,2), optList_bm_median = optList_bm_median_ung, plot.breaks = TRUE,
             ylab = "log bodymass (kg) (ungulate)", xlab = "Time (Ma)", xlim = c(max(intervals), min(intervals)), xaxp = c(70,0,14),                  
             cex.axis = 1, cex.lab = 1, 
             do.subepochs = TRUE,  do.quants = TRUE, manual.breaks = c(57, 51, 39,33, 21, 5))

shoulderPlot(measure.mat = pred.data, plot.y = "BM_all_carnivoran", intervals = intervals, occs = occs, 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = PredrepIntAll$RepsTaxa, quants = bm_quantsPred_allcarniv,
             optList_bm_median = optListBMAllCarniv$OptListTraitMedian, plot.breaks = TRUE, this.rank = "species",
             ylab = "log Bodymass (kg) (all_carnivoran)", xlab = "Time (Ma)", xaxp = c(75,0,15), cex.axis = 1, cex.lab = 1, 
             do.subepochs = TRUE,  do.quants = TRUE, manual.breaks = c(57,41,25,11))
#bm shifts at c(57,41,23) subsampled
###57 due to taxon turnover (all 1000 reps of taxon finds shift here)
###41 also relates to taxon turnover but this also poorly sampled Duchesnean
###25
  

