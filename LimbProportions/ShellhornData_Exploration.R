require(MASS)

schell.dat <- read.csv("/Users/emdoughty/Dropbox/Code/R/LimbProportions/Dataset-Schellhorn_2009-Schellhorn_Pfretzschner_2015.csv", stringsAsFactors = FALSE)
schell.dat$species <- gsub(" ", "_", schell.dat$species)

schell2015 <- c("city", "family", "species", "coll..no.", "hfl", "hdt", " hdl", "rfl", "rpt", "rpl",
                "mc7fl", "mc7pt", "mc7pl",
                "mc7dt", "mc7dl", "mc3fl", "mc3pt", "mc3pl", "mc3dt", " mc3dl", "ffl",
                "tfl", "tpt", "mt7fl", "mt7pt", "mt7pl", "mt7dt", "mt7dl", "mt3fl",
                "mt3pt", " mt3pl", "mt3dt", "mt3dl")
schell.dat <- schell.dat[,colnames(schell.dat) %in% schell2015]

#Shellhorn 2015
hab.Shell <- read.csv("/Users/emdoughty/Dropbox/Code/R/LimbProportions/Shellhorn_HabCalls_2015.csv")
hab.Shell$Species <- gsub(" ","_", hab.Shell$Species)
hab.Shell$Habitat <- gsub(" ","_", hab.Shell$Habitat)
hab.Shell$Habitat <- gsub("-","_", hab.Shell$Habitat)

hab.MP07 <- read.csv("/Users/emdoughty/Dropbox/CalcanealGearRatio/ExtantUngulate_BM_Eco_Mendoza_etal2007.csv", stringsAsFactors = FALSE)
hab.MP07$taxon <- paste(hab.MP07$Genus,hab.MP07$Species, sep="_") 

for(xx in seq(5, ncol(shell.dat),1)) shell.dat[,xx] <- as.numeric(shell.dat[,xx])

#this.dat <- data.frame(aggregate(shell.dat[sapply(shell.dat, class)=="numeric"], by=list(shell.dat$species), FUN=mean, na.rm=TRUE), row.names=1)
this.dat <- schell.dat
this.dat <- this.dat[complete.cases(this.dat),]
#this.dat <- this.dat[this.dat$species %in% hab.Shell$Species,]

dat.string <- this.dat[,c(1:4)]
this.dat <- this.dat[,c(5:ncol(this.dat))]

this.dat <- log(this.dat)

pca <- prcomp(this.dat, scale=TRUE)	

this.dat <- cbind(dat.string, this.dat)

this.dat <- data.frame(this.dat, #taxon=this.dat$taxon, 
                       habitatShell=sapply(this.dat$species, function(x) unique(hab.Shell$Habitat[hab.Shell$Species==x]))) 

this.miss <- this.dat[!this.dat$species %in% hab.Shell$Species,]
this.dat <- this.dat[this.dat$species  %in% hab.Shell$Species,]

this.miss[!this.miss$species %in% hab.MP07$taxon,]

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

this.x <- 2
this.y <- 3

#quartz(h=11, w=8)
quartz()
#par(mfrow = c(1,2))

#plot 1
plot(pca$x[ , this.x], pca$x[ , this.y], type="n", 
     xlab=paste0("PC", this.x, " (", round(100*pca$sdev[this.x]^2/sum(pca$sdev^2), 2), "%)"), 
     ylab=paste0("PC", this.y, " (", round(100*pca$sdev[this.y]^2/sum(pca$sdev^2), 2), "%)"),
     main = "Schellhorn 2015")#    xlim=c(-4,4), ylim = c(-4,4))

abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

points(pca$x[ , this.x], pca$x[ , this.y], pch=21, col=this.colors[as.character(this.dat$habitat)], 
       bg=this.colors[as.character(this.dat$habitat)])

text(pca$x[ , this.x], pca$x[ , this.y], labels=this.dat$species, pos=4, cex=0.5, 
     col=this.colors[as.character(factor(this.dat$habitat))])

lda(this.dat[,5:20],this.dat$habitatShell)

