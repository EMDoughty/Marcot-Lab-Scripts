# # setwd("C:/Users/Evan/Dropbox/ungulate_RA/RCode")

#sources for Jon Marcot's code and specimen measurements 
source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/occFns.R")
source("~/Dropbox/code/R/common_src/sampling.R") 
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R") 

source('~/Dropbox/code/R/dentalMeasurements/src/src_dentalDataFns.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_ecologyAnalysisFns.R', chdir = TRUE)

####################################################################################################################################

occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

####################################################################################################################################

measure.mat <- getMeasureMatWithBodyMasses()

####################################################################################################################################
#### reduces matrix to just the focal order(s)
####################################################################################################################################

# focal.order <- "Artiodactyla"
# focal.order <- "Perissodactyla"
focal.order <- c("Artiodactyla", "Perissodactyla")
bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFam <- sort(unique(bigList$family[bigList$order %in% focal.order]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)

# matrix(measure.mat$taxon[!measure.mat$taxon %in% bigList$accepted_name[bigList$order %in% focal.order]], ncol=1)
# matrix(bigList$accepted_name[bigList$order %in% focal.order][! bigList$accepted_name[bigList$order %in% focal.order] %in% measure.mat$taxon], ncol=1)
measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name[bigList$order %in% focal.order], ]

ungulate.occs <- occs[occs$accepted_name %in% rownames(measure.mat),c("occurrence_no", "accepted_name", "identified_name","order", "family","max_ma", "min_ma", "early_interval", "late_interval", "lng", "lat", "formation", "collection_name", "cc","state","county", "occurrence_comments")]
ungulate.occs <- ungulate.occs[order(ungulate.occs$accepted_name),]
#write.csv(ungulate.occs, "/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/Ungulate_Occurences&Locality.csv")

ungLoc.list <- list()
state.list <- list()

unique(ungulate.occs$state)
for(xx in seq(1, length(unique(ungulate.occs$state)),1)) 
{
  ungLoc.local <- unique(ungulate.occs$collection_name[ungulate.occs$state %in% unique(ungulate.occs$state)[xx]])
  state <- unique(ungulate.occs$state)[xx]
 
  for(yy in seq(1, length(ungLoc.local),1)) 
  {
    ungLoc.species <- ungulate.occs[ungulate.occs$collection_name %in% ungLoc.local[yy],]
    state.list[[yy]] <- ungLoc.species
    names(state.list)[[yy]] <- as.character(unique(ungLoc.local[yy]))
    
    #name.folder <- paste("/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/LocalityByFile/", paste(state,"/",sep=""),sep="")
    #name.file <- paste(name.folder, paste(ungLoc.local[yy], ".csv", sep =""),sep="")
    #write.csv(ungLoc.species, name.file)
  }
  ungLoc.list[[xx]] <- state.list
}
names(ungLoc.list) <- unique(ungulate.occs$state)
#save(ungLoc.list, file="/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/Ungulate_OccurencesByLocality.Rdata")
#write(ungLoc.list, file = "/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/Ungulate_Occurences&Localities_byLocality.txt")

for(xx in unique(ungulate.occs$state)) 
{
  state.occs <- ungulate.occs[ungulate.occs$state == xx,]

name.folder <- "/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/LocalityByFile/"
#name.folder <- paste("/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/LocalityByFile/", paste(state,"/",sep=""),sep="")
name.file <- paste(name.folder, paste(xx, ".csv", sep =""),sep="")
write.csv(state.occs, name.file)
}

