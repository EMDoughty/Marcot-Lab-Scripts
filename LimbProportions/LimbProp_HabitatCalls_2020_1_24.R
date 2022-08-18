mcmahon.master <- read.csv("/Users/emdoughty/Dropbox/CalcanealGearRatio/Limbs/mcMahon1975.csv", stringsAsFactors = FALSE)
mcmahon.master$taxon

mcmahon.master$taxon <- gsub(" ", "_", mcmahon.master$current)

test <- mcmahon.master

for(xx in seq(7, 16,1)) mcmahon.master[,xx] <- as.numeric(mcmahon.master[,xx])

this.dat <- data.frame(aggregate(mcmahon.master[sapply(mcmahon.master, class)=="numeric"], by=list(mcmahon.master$taxon), FUN=mean, na.rm=TRUE), row.names=1)
# this.dat <- mcmahon.master
#rownames(this.dat) <- paste(this.dat$taxon,rownames(this.dat), sep = "_"); taxa <- this.dat$taxon; this.dat <- this.dat[,c(7:16)]
 

taxComplete <- which(complete.cases(this.dat)); this.dat <- this.dat[complete.cases(this.dat),]
this.dat[,1:10] <- log(this.dat[,1:10])

pca <- prcomp(this.dat, scale=TRUE)	

############################################################################################################################################
#### append habitat to transformed data matrix
############################################################################################################################################
#Mendoza and Palmqvist 2007/2008
habitat <- read.csv("/Users/emdoughty/Dropbox/CalcanealGearRatio/ExtantUngulate_BM_Eco_Mendoza_etal2007.csv", stringsAsFactors = FALSE)
habitat$taxon <- paste(habitat$Genus, habitat$Species, sep = "_")

#Barr 2014/2017 (bovids) and Carran 20xx (cervids) habitat calls
dat.master <- read.csv("/Users/emdoughty/Dropbox/CalcanealGearRatio/Photos/10914_2018_9446_MOESM1_ESM.csv", stringsAsFactors=FALSE)
barr_carran <- data.frame(aggregate(dat.master[sapply(dat.master, class)=="numeric"], by=list(dat.master$taxon), FUN=mean, na.rm=TRUE), row.names=1)
#barr_carran <- barr_carran[complete.cases(barr_carran),]
barr_carran <- data.frame(barr_carran, taxon=rownames(barr_carran), habitat=sapply(rownames(barr_carran), function(x) unique(dat.master$habitat[dat.master$taxon==x])), Family=sapply(rownames(barr_carran), function(x) unique(dat.master$Family[dat.master$taxon == x])))

#Kovarovic 2004 thesis habitat calls
hab.Kova <- read.csv("/Users/emdoughty/Dropbox/Code/R/LimbProportions/Kovarovic_Thesis_HabCalls_2004.csv")
hab.Kova$Species <- gsub(" ","_", hab.Kova$Species)
hab.Kova$Habitat <- gsub(" ","_", hab.Kova$Habitat)
hab.Kova$Habitat <- gsub("-","_", hab.Kova$Habitat)

#Shellhorn 2015
hab.Shell <- read.csv("/Users/emdoughty/Dropbox/Code/R/LimbProportions/Shellhorn_HabCalls_2015.csv")
hab.Shell$Species <- gsub(" ","_", hab.Shell$Species)
hab.Shell$Habitat <- gsub(" ","_", hab.Shell$Habitat)
hab.Shell$Habitat <- gsub("-","_", hab.Shell$Habitat)

# this.dat$taxon <- taxa[taxComplete]

this.miss <- this.dat[!rownames(this.dat) %in% habitat$taxon,]
this.dat <- this.dat[rownames(this.dat) %in% habitat$taxon,]
this.dat <- data.frame(this.dat, taxon=rownames(this.dat), 
                       habitatMP07=sapply(rownames(this.dat), function(x) unique(habitat$Habitat[habitat$taxon==x])), 
                       Family=sapply(rownames(this.dat), function(x) unique(habitat$Family[habitat$taxon == x])))
this.miss$taxon <- rownames(this.miss)
this.miss$habitatMP07 <- "none"
this.miss$Family <- c("Bovidae", "Bovidae", "Bovidae", "Bovidae", "Bovidae", 
                      "Cervidae", 
                      "Equidae", "Equidae", "Equidae", "Equidae", "Equidae", "Equidae",
                      "Bovidae", "Bovidae", "Bovidae", "Bovidae", "Bovidae", "Bovidae", 
                      "Bovidae", "Bovidae", "Bovidae", "Bovidae","Bovidae", "Bovidae", "Bovidae")

this.dat <- rbind(this.dat, this.miss)

this.miss <- this.dat[!rownames(this.dat) %in% barr_carran$taxon,]
this.dat <- this.dat[rownames(this.dat) %in% barr_carran$taxon,]
this.dat <- data.frame(this.dat, 
                       habBarrCarran=sapply(rownames(this.dat), function(x) unique(barr_carran$habitat[barr_carran$taxon==x])))
this.miss$habBarrCarran <- "none"
this.dat <- rbind(this.dat, this.miss)

this.miss <- this.dat[!rownames(this.dat) %in% hab.Kova$Species,]
this.dat <- this.dat[rownames(this.dat) %in% hab.Kova$Species,]
this.dat <- data.frame(this.dat, 
                       habKov04=sapply(rownames(this.dat), function(x) unique(hab.Kova$Habitat[hab.Kova$Species==x])))
this.miss$habKov04 <- "none"
this.dat <- rbind(this.dat, this.miss)

this.miss <- this.dat[!rownames(this.dat) %in% hab.Shell$Species,]
this.dat <- this.dat[rownames(this.dat) %in% hab.Shell$Species,]
this.dat <- data.frame(this.dat, 
                       habShell15=sapply(rownames(this.dat), function(x) unique(hab.Shell$Habitat[hab.Shell$Species==x])))
this.miss$habShell15 <- "none"
this.dat <- rbind(this.dat, this.miss)

#this.dat[rownames(this.dat) %in% habitat$taxon,]

this.colors <- c(CH="sienna", Forest="sienna", HeavyCover="darkgreen", MH="lightgreen", 
                 LightCover="lightgreen",  OH="khaki", Open="khaki", Grassland="khaki", 
                 Wooded_bushed_Grassland="lightgreen", Light_Woodland_Bushland="green", 
                 Heavy_Woodland_Bushland="darkgreen", Montane_Light_Cover="orange", 
                 Mountainous="orange", Montane_Heavy_Cover="red", "none"="black")

this.colorsNA <- c(CH="sienna", Forest = "sienna", HeavyCover = "darkgreen", MH="lightgreen", 
                   LightCover="lightgreen",  OH="khaki", Open="khaki", Grassland="khaki", 
                   Wooded_bushed_Grassland="lightgreen", Light_Woodland_Bushland="green", 
                   Heavy_Woodland_Bushland="darkgreen", Montane_Light_Cover="orange", 
                   Mountainous="orange", Montane_Heavy_Cover="red")

shape.colors <- rainbow(length(unique(this.dat$taxon)))

fam.colors <- rainbow(length(unique(this.dat$Family)))
names(fam.colors) <- sort(unique(this.dat$Family))

this.x <- 2
this.y <- 3

habitatMP07_NA <- gsub("none", NA, this.dat$habitatMP07)
habBarrCarran_NA <- gsub("none", NA, this.dat$habBarrCarran)
this.dat_07 <- this.dat[!this.dat$habitatMP07 == "none",]
this.dat_BC <- this.dat[!this.dat$habBarrCarran == "none",]

# par(mfrow=c(1,2))
quartz(h=11, w=8)
par(mfrow = c(4,2))

#plot 1
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Mendoza $ Palmqvist 2007")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habitatMP07)], 
       bg=this.colors[as.character(this.dat$habitatMP07)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=this.dat$taxon, pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habitatMP07))])

#2
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Mendoza $ Palmqvist 2007")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habitatMP07)], 
       bg=this.colorsNA[as.character(this.dat$habitatMP07)])

#3
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Barr/Carran")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colors[as.character(this.dat$habBarrCarran)], 
       bg=this.colors[as.character(this.dat$habBarrCarran)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=this.dat$taxon, pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habBarrCarran))])

#4
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Barr/Carran")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habBarrCarran)], 
       bg=this.colorsNA[as.character(this.dat$habBarrCarran)])

#5
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Kovarovic 2004")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colors[as.character(this.dat$habKov04)], 
       bg=this.colors[as.character(this.dat$habKov04)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=rownames(this.dat), pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habKov04))])

#6
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Kovarovic 2004")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habKov04)], 
       bg=this.colorsNA[as.character(this.dat$habKov04)])

#7
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Shellhorn 2015")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colors[as.character(this.dat$habShell15)], 
       bg=this.colors[as.character(this.dat$habShell15)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=rownames(this.dat), pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habShell15))])

#8
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Shellhorn 2015")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habShell15)], 
       bg=this.colorsNA[as.character(this.dat$habShell15)])

################################################################################################################################
#Specimen PCA
################################################################################################################################
mcmahon.master <- read.csv("/Users/emdoughty/Dropbox/CalcanealGearRatio/Limbs/mcMahon1975.csv", stringsAsFactors = FALSE)
mcmahon.master

mcmahon.master$taxon <- gsub(" ", "_", mcmahon.master$current)

test <- mcmahon.master

for(xx in seq(7, 16,1)) mcmahon.master[,xx] <- as.numeric(mcmahon.master[,xx])

#this.dat <- data.frame(aggregate(mcmahon.master[sapply(mcmahon.master, class)=="numeric"], by=list(mcmahon.master$taxon), FUN=mean, na.rm=TRUE), row.names=1)
this.dat <- mcmahon.master
rownames(this.dat) <- paste(this.dat$taxon,rownames(this.dat), sep = "_"); taxa <- this.dat$taxon; this.dat <- this.dat[,c(7:16)]

taxComplete <- which(complete.cases(this.dat)); this.dat <- this.dat[complete.cases(this.dat),]
this.dat[,1:10] <- log(this.dat[,1:10])

pca <- prcomp(this.dat, scale=TRUE)	

############################################################################################################################################
#### append habitat to transformed data matrix
############################################################################################################################################
#Mendoza and Palmqvist 2007/2008
habitat <- read.csv("/Users/emdoughty/Dropbox/CalcanealGearRatio/ExtantUngulate_BM_Eco_Mendoza_etal2007.csv", stringsAsFactors = FALSE)
habitat$taxon <- paste(habitat$Genus, habitat$Species, sep = "_")

#Barr 2014/2017 (bovids) and Carran 20xx (cervids) habitat calls
dat.master <- read.csv("/Users/emdoughty/Dropbox/CalcanealGearRatio/Photos/10914_2018_9446_MOESM1_ESM.csv", stringsAsFactors=FALSE)
barr_carran <- data.frame(aggregate(dat.master[sapply(dat.master, class)=="numeric"], by=list(dat.master$taxon), FUN=mean, na.rm=TRUE), row.names=1)
barr_carran <- barr_carran[complete.cases(barr_carran),]
barr_carran <- data.frame(barr_carran, taxon=rownames(barr_carran), habitat=sapply(rownames(barr_carran), function(x) unique(dat.master$habitat[dat.master$taxon==x])), Family=sapply(rownames(barr_carran), function(x) unique(dat.master$Family[dat.master$taxon == x])))

#Kovarovic 2004 thesis habitat calls
hab.Kova <- read.csv("/Users/emdoughty/Dropbox/Code/R/LimbProportions/Kovarovic_Thesis_HabCalls_2004.csv")
hab.Kova$Species <- gsub(" ","_", hab.Kova$Species)
hab.Kova$Habitat <- gsub(" ","_", hab.Kova$Habitat)
hab.Kova$Habitat <- gsub("-","_", hab.Kova$Habitat)

#Shellhorn 2015
hab.Shell <- read.csv("/Users/emdoughty/Dropbox/Code/R/LimbProportions/Shellhorn_HabCalls_2015.csv")
hab.Shell$Species <- gsub(" ","_", hab.Shell$Species)
hab.Shell$Habitat <- gsub(" ","_", hab.Shell$Habitat)
hab.Shell$Habitat <- gsub("-","_", hab.Shell$Habitat)

#get taxon names only so code below doesn't break
this.dat$taxon <- taxa[taxComplete]

this.miss <- this.dat[!this.dat$taxon %in% habitat$taxon,]
this.dat <- this.dat[this.dat$taxon %in% habitat$taxon,]
this.dat <- data.frame(this.dat, #taxon=this.dat$taxon, 
                       habitatMP07=sapply(this.dat$taxon, function(x) unique(habitat$Habitat[habitat$taxon==x]))) 
                       #Family=sapply(this.dat), function(x) unique(habitat$Family[habitat$taxon == x])))
this.miss$habitatMP07 <- "none"
#this.miss$Family <- c("Bovidae", "Bovidae", "Bovidae", "Bovidae", "Bovidae", 
#                      "Cervidae", 
#                      "Equidae", "Equidae", "Equidae", "Equidae", "Equidae", "Equidae",
#                      "Bovidae", "Bovidae", "Bovidae", "Bovidae", "Bovidae", "Bovidae", 
#                      "Bovidae", "Bovidae", "Bovidae", "Bovidae","Bovidae", "Bovidae", "Bovidae")

this.dat <- rbind(this.dat, this.miss)

this.miss <- this.dat[!this.dat$taxon %in% barr_carran$taxon,]
this.dat <- this.dat[this.dat$taxon %in% barr_carran$taxon,]
this.dat <- data.frame(this.dat, 
                       habBarrCarran=sapply(this.dat$taxon, function(x) unique(barr_carran$habitat[barr_carran$taxon==x])))
this.miss$habBarrCarran <- "none"
this.dat <- rbind(this.dat, this.miss)

this.miss <- this.dat[!this.dat$taxon %in% hab.Kova$Species,]
this.dat <- this.dat[this.dat$taxon %in% hab.Kova$Species,]
this.dat <- data.frame(this.dat, 
                       habKov04=sapply(this.dat$taxon, function(x) unique(hab.Kova$Habitat[hab.Kova$Species==x])))
this.miss$habKov04 <- "none"
this.dat <- rbind(this.dat, this.miss)

this.miss <- this.dat[!this.dat$taxon %in% hab.Shell$Species,]
this.dat <- this.dat[this.dat$taxon %in% hab.Shell$Species,]
this.dat <- data.frame(this.dat, 
                       habShell15=sapply(this.dat$taxon, function(x) unique(hab.Shell$Habitat[hab.Shell$Species==x])))
this.miss$habShell15 <- "none"
this.dat <- rbind(this.dat, this.miss)

#this.dat[rownames(this.dat) %in% habitat$taxon,]

this.colors <- c(CH="sienna", Forest="sienna", HeavyCover="darkgreen", MH="lightgreen", 
                 LightCover="lightgreen",  OH="khaki", Open="khaki", Grassland="khaki", 
                 Wooded_bushed_Grassland="lightgreen", Light_Woodland_Bushland="green", 
                 Heavy_Woodland_Bushland="darkgreen", Montane_Light_Cover="orange", 
                 Mountainous="orange", Montane_Heavy_Cover="red", "none"="black")

this.colorsNA <- c(CH="sienna", Forest = "sienna", HeavyCover = "darkgreen", MH="lightgreen", 
                   LightCover="lightgreen",  OH="khaki", Open="khaki", Grassland="khaki", 
                   Wooded_bushed_Grassland="lightgreen", Light_Woodland_Bushland="green", 
                   Heavy_Woodland_Bushland="darkgreen", Montane_Light_Cover="orange", 
                   Mountainous="orange", Montane_Heavy_Cover="red")

#fam.colors <- rainbow(length(unique(this.dat$Family)))
#names(fam.colors) <- sort(unique(this.dat$Family))

this.x <- 2
this.y <- 3

habitatMP07_NA <- gsub("none", NA, this.dat$habitatMP07)
habBarrCarran_NA <- gsub("none", NA, this.dat$habBarrCarran)
this.dat_07 <- this.dat[!this.dat$habitatMP07 == "none",]
this.dat_BC <- this.dat[!this.dat$habBarrCarran == "none",]

#need to make sure that the data in this.dat in same order as PCA
this.dat <- this.dat[rownames(pca$x),]

# par(mfrow=c(1,2))
quartz(h=11, w=8)
par(mfrow = c(4,2))

#plot 1
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Mendoza $ Palmqvist 2007")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colors[as.character(this.dat$habitatMP07)], 
       bg=this.colors[as.character(this.dat$habitatMP07)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=rownames(this.dat), pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habitatMP07))])

#2
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Mendoza $ Palmqvist 2007")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habitatMP07)], 
       bg=this.colorsNA[as.character(this.dat$habitatMP07)])

#text(pca$x[ , this.x], pca$x[ , this.y], labels=this.dat$taxon, pos=4, cex=0.5, 
#    col=this.colorsNA[as.character(factor(this.dat$habitatMP07))])

#3
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Barr/Carran")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colors[as.character(this.dat$habBarrCarran)], 
       bg=this.colors[as.character(this.dat$habBarrCarran)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=rownames(this.dat), pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habBarrCarran))])

#4
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Barr/Carran")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habBarrCarran)], 
       bg=this.colorsNA[as.character(this.dat$habBarrCarran)])

#5
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Kovarovic 2004")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colors[as.character(this.dat$habKov04)], 
       bg=this.colors[as.character(this.dat$habKov04)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=rownames(this.dat), pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habKov04))])

#6
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Kovarovic 2004")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habKov04)], 
       bg=this.colorsNA[as.character(this.dat$habKov04)])

#7
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Shellhorn 2015")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colors[as.character(this.dat$habShell15)], 
       bg=this.colors[as.character(this.dat$habShell15)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=rownames(this.dat), pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habShell15))])

#8
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Shellhorn 2015")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colorsNA[as.character(this.dat$habShell15)], 
       bg=this.colorsNA[as.character(this.dat$habShell15)])

