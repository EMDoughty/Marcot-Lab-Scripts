extract.SpecimenNumber_old <- function(tpsArray, type = "Path", metaData=NULL)
{
  image.names <- as.character(dimnames(tpsArray)[[3]])
  #if(type = "Path")
  Split <- strsplit(image.names, "\\\\")
  name.vec <- vector()
  image.vec <- vector()
  for(xx in seq(1, length(Split),1)) 
  {
    name.vec[xx] <- Split[[xx]][max(length(Split[[xx]]))] 
    #image.vec[xx] <- Split[[xx]][max(length(Split[[xx]]))-1] 
  }
  #if there are " _ " then do substr prior to that
  specimen.vec <- substr(name.vec, start = 1, stop = regexpr(" - ", name.vec)-1)
  image.vec <- substr(name.vec, start = regexpr(" - ", name.vec)+2, stop = nchar(name.vec))
  image.vec <- gsub(" - Copy.*","\\1",image.vec)
  
  #if there are no " _ " then just have the entry be name.vec or at very end just duplicate name.vec
  for (xx in seq(1, length(specimen.vec),1)) if(specimen.vec[xx] %in% "") specimen.vec[xx] <- name.vec[xx]
  
  specimen.mat <- cbind(name.vec, specimen.vec, image.vec)
  
  specimen.mat <- as.data.frame(specimen.mat)
  specimen.mat[,"image.vec"] <- as.integer(as.character(specimen.mat[,"image.vec"]))
  specimen.mat[is.na(specimen.mat$image.vec),"image.vec"] <- 1
  
  #could move checkMissingLan here
  specimen.mat <- checkMissingLand(landfile = tpsArray, specimens = specimen.mat)
  #could  move get Binomials here
  specimen.mat <- getBinomials(specimens = specimen.mat, metafile = metaData)
  
  specimen.mat <-  getSide(specimens = specimen.mat, metafile = metaData)
  
  return(specimen.mat)
}

checkMissing4Reg <- function (spec.array)
{
  #check for and removes specimens with missing landmarks since
  #count # rows with NA's for each 
  
  if(!is.matrix(spec.array)) spec.mat.2d <- two.d.array(spec.array) else spec.mat.2d <- spec.array
  
  spec.count <- nrow(spec.mat.2d)
  
  spec.remove <- vector()
  for(xx in seq(1, nrow(spec.mat.2d),1)) 
  {
    check <- sum(is.na(spec.mat.2d[xx,]))/2
    
    if(check > 0)
    {
      spec.need <- 2*check
      
      if(spec.need > check) #how to remove specimens.....
        spec.remove[xx] <- rownames(spec.mat.2d)[xx]
    }
  }
  #print which specimens get removed
  spec.remove <- na.omit(spec.remove)
  print("Removing: \n") 
  print(spec.remove)
  
  spec.mat.2d <- spec.mat.2d[!rownames(spec.mat.2d) %in% spec.remove,]
  spec.mat.2d <- arrayspecs(A = spec.mat.2d, p = ncol(spec.mat.2d)/2, k = 2)
  
  return (spec.mat.2d)
}

getPCA <- function(specimen.array, colgroup, pointLabel, imageLabel = "", useLegend = FALSE, legend.pos.x = "topleft", legend.pos.y = NULL, this.x = 1, this.y = 2, scale.load = 1)
{
  #Procrustes
  specimen.proc <- gpagen(A = specimen.array)
  plot(specimen.proc)
  
  #PCA
  #PCA as per tutorial
  #3)plotTangentSpace.....Depricated function use gm.pcrcomp instead
  specimen.pca <- gm.prcomp(specimen.proc$coords)
  specimen.col <- rainbow(length(unique(colgroup)))
  specimen.col.vec <- as.character(colgroup)
  
  #for(xx in seq(1,length(unique(specimen.col.vec)),1)) specimen.col.vec <- gsub(as.integer(as.character(unique(specimen.col.vec)[xx])),
   #                                                                             specimen.col[xx],specimen.col.vec)
  
  for(xx in seq(1,length(specimen.col),1)) specimen.col.vec[which(specimen.col.vec %in% unique(specimen.col.vec)[xx])] <- specimen.col[xx]
  
  #dev.off()
  plot(specimen.pca)
  points(specimen.pca$x, col = specimen.col.vec, pch =16)
  text(specimen.pca$x, paste(as.character(pointLabel), as.character(imageLabel),sep=""), cex = 0.5) #, adj = -0.15)
  
  if(useLegend == TRUE) legend(x = legend.pos.x, y = legend.pos.y, legend = unique(colgroup), fill = specimen.col)
  
 # plot(specimen.pca$rotation)
  rownames(specimen.pca$rotation) <- seq(1, nrow(specimen.pca$rotation),1)
  
  #set rownames to be #.X and #.y setup
  for(xx in seq(1, nrow(specimen.pca$rotation)/2,1))
  {
    rownames(specimen.pca$rotation)[(xx*2)-1] <- paste(xx,".x",sep="")
    rownames(specimen.pca$rotation)[(xx*2)] <- paste(xx,".y",sep="")
  }
  rownames(specimen.pca$rotation)
  
  for(xx in seq(1, nrow(specimen.pca$rotation),1)) {
    lines(c(0,specimen.pca$rotation[xx,this.x]*scale.load), c(0,specimen.pca$rotation[xx,this.y]*scale.load), col = "red")
    text(x = specimen.pca$rotation[xx,this.x]*scale.load, y = specimen.pca$rotation[xx, this.y]*scale.load, 
         labels = rownames(specimen.pca$rotation)[xx], col = "red")
  }
  
  return(specimen.pca)
}

getPCA.linear <- function(specimen.array, colgroup, pointLabel, imageLabel = "", useLegend = FALSE, legend.pos.x = "topleft", legend.pos.y = NULL, this.x = 1, this.y = 2, loadings = TRUE, scale.load = 1, load.cex = 1)
{
  
  #PCA
  #PCA as per tutorial
  #3)plotTangentSpace.....Depricated function use gm.pcrcomp instead
  specimen.pca <- prcomp(specimen.array)
  specimen.col <- rainbow(length(unique(colgroup)))####
  specimen.col.vec <- as.character(colgroup)
  
  specimen.var <-   summary(specimen.pca)
  
  #for(xx in seq(1,length(unique(specimen.col.vec)),1)) specimen.col.vec <- gsub(as.integer(as.character(unique(specimen.col.vec)[xx])),
  #                                                                             specimen.col[xx],specimen.col.vec)
  
  for(xx in seq(1,length(specimen.col),1)) specimen.col.vec[which(specimen.col.vec %in% unique(specimen.col.vec)[xx])] <- specimen.col[xx]
  
  #dev.off()
  plot(specimen.pca$x[,this.x],specimen.pca$x[,this.y], 
       xlab = paste(colnames(specimen.var$importance)[this.x],specimen.var$importance[2,this.x]*100, "%",sep=" "), 
       ylab = paste(colnames(specimen.var$importance)[this.y], specimen.var$importance[2,this.y]*100, "%", sep=" "))
  points(specimen.pca$x[,this.x],specimen.pca$x[,this.y], col = specimen.col.vec, pch =16)
  text(specimen.pca$x[,c(this.x,this.y)], paste(as.character(pointLabel), as.character(imageLabel),sep=""), cex = 0.5) #, adj = -0.15)
  
  if(useLegend == TRUE) legend(x = legend.pos.x, y = legend.pos.y, legend = unique(colgroup), fill = specimen.col)
  
  # plot(specimen.pca$rotation)
  
  #set rownames to be #.X and #.y setup
  rownames(specimen.pca$rotation)
  
  if(loadings == TRUE)
  {
    for(xx in seq(1, nrow(specimen.pca$rotation),1)) {
      lines(c(0,specimen.pca$rotation[xx,this.x]*scale.load), c(0,specimen.pca$rotation[xx,this.y]*scale.load), col = "red")
      text(x = specimen.pca$rotation[xx,this.x]*scale.load, y = specimen.pca$rotation[xx, this.y]*scale.load, 
           labels = rownames(specimen.pca$rotation)[xx], col = "red", cex = load.cex)
    }
  }
  
  return(specimen.pca)
}

#need function to align figures to be on the same side
align2Side_new <- function(landfile, metafile, side = "Right", index.X = rep(0, nrow(specimens)))
{
  #specimens <- cbind(specimens, siding.y)
  meta.side <- tolower(metafile$Side)
  
  flipY <- rep(0,nrow(metafile))
  
  for(xx in seq(1,nrow(metafile),1)) 
  {
        if(!tolower(metafile[xx,"Side"]) %in% tolower(side)) flipY[xx] <-  1
  }
  
  #rotated.land <- rotate.coords(landfile, type = "flipX", index = index.X)
  #made to assume that proximal is on left and distal on the right
  ###could do two "rotateCC" to get the photo aligned 
  rotated.land <- rotate.coords(landfile, type = "flipX", index = metafile$Flip.X)
  
  rotated.land <- rotate.coords(rotated.land, type = "flipY", index = flipY)
  
  dimnames(rotated.land) <- dimnames(landfile)
  
  return(rotated.land)
}

#function to search throuhg landmark array and mark which specimens are missing landmarks and by how many.
checkMissingLand <- function(landfile, specimens)
{
  
  spec.mat <- two.d.array(landfile)
  
  miss.vec <- rep(FALSE, nrow(specimens))
  miss.count <- rep(0, nrow(specimens))
  
  for(xx in seq(1, nrow(specimens),1))
  {
    if(length(table(is.na(spec.mat[xx,]))) > 1)
    {
      miss.vec[xx] <- TRUE
      miss.count[xx] <- as.integer(summary(is.na(spec.mat[xx,]))["TRUE"])/2
      miss.count[xx] <- miss.count[xx]
    }
  }
  
  specimens <- as.data.frame(cbind(specimens, miss.count))
  
  return(specimens)
}

#need funciton to apply species list to specimens order in landmark file
getBinomials <- function(specimens, metafile)
{
  metafile <- cbind(metafile, com.catalog=as.character(paste(metafile$Institutional.Code, metafile$Catalog..,sep=" ")))
  metafile$com.catalog <- as.character(metafile$com.catalog)
  
  for(xx in seq(1,nrow(metafile),1)) 
  {
    if(metafile[xx,"com.catalog"] == "UCMP 84C57-15-1") metafile$com.catalog[xx] <- "84C57-15-1"
    if(metafile[xx,"com.catalog"] == "UCMP 1755_1") metafile$com.catalog[xx] <- "UCMP 1755"
    if(metafile[xx,"com.catalog"] == "UCMP 29878-5") metafile$com.catalog[xx] <- "UCMP 29878"
  }
  
  nrow(metafile)
  nrow(specimens)
  
  species <- rep(NA, nrow(specimens))
  
  for(zz in seq(1, nrow(specimens),1))
  {
    species.index <- which(metafile$com.catalog %in% specimens[zz,2])
    #if(length(side.index) > 1) species[zz] <- metafile$Side[side.index][metafile[side.index,]$Image.. %in% specimens$image.vec[zz]]
    
    species[zz] <- as.character(metafile$Genus[species.index])
   
    #need future versions to also grab specie sname once I add those in
  } 
  
  specimens <- cbind(specimens, species)
  
  specimens <- specimens[complete.cases(specimens),]
  
  return(specimens)
}

#side
getSide <- function(specimens, metafile)
{
  metafile <- cbind(metafile, com.catalog=as.character(paste(metafile$Institutional.Code, metafile$Catalog..,sep=" ")))
  metafile$com.catalog <- as.character(metafile$com.catalog)
  
  for(xx in seq(1,nrow(metafile),1)) 
  {
    if(metafile[xx,"com.catalog"] == "UCMP 84C57-15-1") metafile$com.catalog[xx] <- "84C57-15-1"
    if(metafile[xx,"com.catalog"] == "UCMP 1755_1") metafile$com.catalog[xx] <- "UCMP 1755"
    if(metafile[xx,"com.catalog"] == "UCMP 29878-5") metafile$com.catalog[xx] <- "UCMP 29878"
  }
  
  nrow(metafile)
  nrow(specimens)
  
  specimens$image.vec <- as.numeric(as.character(specimens$image.vec))
  metafile$Image.. <- as.numeric(as.character(metafile$Image..))
  
  side <- rep(NA, nrow(specimens))
  
 # for(zz in seq(1, nrow(specimens),1))
  #{
   # specimen.index <- which(metafile$com.catalog %in% specimens[zz,2] && metafile$Image.. %in% as.numeric(as.character(specimens[zz,3])))
    #if(length(side.index) > 1) species[zz] <- metafile$Side[side.index][metafile[side.index,]$Image.. %in% specimens$image.vec[zz]]
    
    # side[zz] <- as.character(metafile$Side[specimen.index])
    
     #metafile$side[match(metafile$com.catalog,specimens$specimen.vec)]
     
    #remove 
  #} 
  
  #specimens <-  merge(specimens, metafile, sort= FALSE, by.x = c("specimen.vec", "image.vec"), 
  #                    by.y = c("com.catalog", "Image.."), all.x = TRUE)
 
  for(xx in seq(1, nrow(specimens),1))
  {
    for(yy in seq(1, nrow(metafile),1))
    {
      if(as.character(specimens$specimen.vec[xx]) == as.character(metafile$com.catalog[yy]) && 
         as.character(specimens$image.vec[xx]) == as.character(metafile$Image..[yy]))
      {
        side[xx] <- as.character(metafile$Side[yy])
      }
    }
  }
  
  specimens <- cbind(specimens, Side = side)
  
 # specimens <- specimens[, !names(specimens) %in% c("Institutional.Code", "Catalog..", "Genus")]
  
  #specimens <- cbind(specimens, side)
  
  #specimens <- specimens[complete.cases(specimens),]
  
  return(specimens)
}

#need function(s) to clean the data of sliding landmarks if anchor points are not present
removeLand <- function(landfile, land2remove)
{
  land.seq <- seq(1, dim(landfile)[1], 1)
  land.final <- land.seq[!land2remove]
  
  land.array <- landfile[-c(land2remove),,]
  
  return(land.array)                    
}

getLoadings <- function(pca.obj)
{
  summary(pca.obj)
  pca.load <- pca.obj$rotation
  plot(pca.load)
  #pca.load
  
  return()
}

get.list.2remove <- function(meta.data, spec.direct = NA, spec.via.sp = NA)
{
  
  rm.sp <- rownames(meta.data)[meta.data$Binomial %in% spec.via.sp]
  
  output.list <- c(spec.direct,rm.sp)
  
  return(output.list[!is.na(output.list)])
}

removeSpecimen <- function(landfile, spec2remove ="Museum_Catalog#_Image", LD = 1, k=2) 
{
  landmarks <- two.d.array(landfile, sep = ".")
  
  #need to havce spec2remove be the index to remove rather than the name
  landmarks <- landmarks[!rownames(landmarks) %in% spec2remove,]
  
  landmarks <- arrayspecs(landmarks, p = LD, k=k)
  
  return(landmarks)
}

convert2TPS <- function(folder.direct = NA,filetype=".points", meta.data = NA, 
                        element = "Radius", outputfilename = "output.tps")
{
  meta.data <- meta.data[meta.data$Element %in% element,]
  
  #folder.direct <- "/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis/Radius/"
  filenames <- file.info(list.files(folder.direct, pattern = filetype, full.names=TRUE))
  
  #if(type = "Path")
  Split <- strsplit(rownames(filenames), "/")
  name.vec <- vector()
  image.vec <- vector()
  for(xx in seq(1, length(Split),1)) 
  {
    name.vec[xx] <- Split[[xx]][max(length(Split[[xx]]))] 
    #image.vec[xx] <- Split[[xx]][max(length(Split[[xx]]))-1] 
  }
  
  name.vec <- gsub(".points", "", name.vec)
  name.vec <- gsub("- ", "-", name.vec)
  name.vec <- strsplit(name.vec, " ")
  
  name.mat <- as.data.frame(matrix(nrow = length(name.vec), ncol = 3))
  
  for(yy in seq(1, nrow(name.mat),1)) 
  {
    name.mat[yy,1] <- name.vec[[yy]][1]
    name.mat[yy,2] <- name.vec[[yy]][2]
    name.mat[yy,3] <- gsub("-","", name.vec[[yy]][3])
  }
  
  colnames(name.mat) <- c("Museum", "CatalogNo", "ImageNo")
  
  #order filenames and meta.data to be the same (may be easier to get the latter in same order)
  meta.data
  nrow(name.mat)
  nrow(meta.data)
  rownames(name.mat)  <- paste(name.mat$Museum,name.mat$CatalogNo,name.mat$ImageNo,sep="_")
  rownames(meta.data) <- paste(meta.data$Museum.ID,meta.data$Catalog..,meta.data$Image,sep="_")
  meta.data$Binomial <- paste(meta.data$Accepted.Genus,meta.data$Accepted.Species, sep="_")
  
  meta.data <- meta.data[rownames(name.mat),]
  
  #convert .points files into readable format
  data.output <- vector()
  
  for(zz in seq(1,length(rownames(filenames)),1))
  {
   file.full <- readLines(rownames(filenames)[zz])
    
   lines.index <-  which(str_detect(file.full, "namedpointset"))
   file.land <- file.full[lines.index[(length(lines.index)-1)]:lines.index[(length(lines.index))]]
   
   lines.index <-  which(str_detect(file.land, "namedpointset"))
   land.data <- file.land[-lines.index]
   
   land.data <- strsplit(land.data, " " )
   
   #get # landmarks automatically here
   LM <- length(land.data)
   
   #need to import metadata file and scour it for the scale for this specimen
   spec.scale <- 1/as.numeric(meta.data$Scale..pixel.mm.[rownames(meta.data) %in% rownames(name.mat)[zz]])
   
   land.data <- lapply(land.data, function (x) x[8:9])
   
   #translate to matrix
   land.data <- matrix(unlist(land.data),nrow=LM,ncol=2, byrow = TRUE)
   
   colnames(land.data) <- c("X","Y")
    
   #remove uneccessary characters so only numbers remain (double)
   land.data <- gsub("\"","", land.data)
   land.data <- gsub("x=","", land.data)
   land.data <- gsub("y=","", land.data)
   
   land.data <- paste(land.data[,1],land.data[,2],sep=" ")
   #land.data <-  apply(land.data,c(1,2), as.double)
    
    
   spec.output <- c(paste("LM=",LM, sep=""),land.data, paste("IMAGE=",rownames(filenames)[zz],sep = ""),
     paste("ID=",paste(meta.data$Accepted.Genus[zz],meta.data$Accepted.Species[zz],sep="_"),sep=""),
     paste("SCALE=",spec.scale,sep=""))
   
   data.output <- append(data.output,spec.output)
  }
  #make the output variable into a composite and single specimen objects
  ##use single object to compile and format before appendening to the 
  
  writeLines(text = data.output,sep="\n", con=outputfilename)
  
  return(meta.data)
}

getEstimate <- function(landmarks, meta.data, method = "TPS")
{
  check <- checkMissingLand(landfile = landmarks, specimens = meta.data)
  
  if(max(check$miss.count) > 0)
  {
    landmarks <- estimate.missing(landmarks, method = method)
  } else {
    landmarks <- landmarks
  }
  
  return(landmarks)
}

get.linear.measures <- function(landmarks, metadata, meas.array = NA, meas.names = NA)
{
  meas.Mat <- matrix(meas.array, ncol=2,byrow = TRUE,
                     dimnames = list(meas.names, c("start","end")))
  element.lengths <- interlmkdist(landmarks,meas.Mat)
  
  #need to apply scaling to get final output
  meta.scale <- metadata$Scale..pixel.mm.
  
  for(xx in seq(1,nrow(element.lengths),1))
  {
    element.lengths[xx,] <- element.lengths[xx,]*meta.scale[xx]
  }
  
  name.index <- which(str_detect(colnames(element.lengths), pattern=".Length"))
  
  #Need to add a method to detect and fill NA values; might be solves using nar.rm=TRUE
  element.tot.length.Ave <- rowMeans(element.lengths[,name.index],na.rm = TRUE)
  elements.lengths <- cbind(element.lengths, element.tot.length.Ave)
  
  #make a name for the average length
  ###pull the first portion of the measure names to get element
  elem.used  <- strsplit(x = meas.names, split = "[.]")[[1]][1]
  
  new.names  <- append(meas.names, paste(elem.used, ".tot.length.Ave",sep = ""))
  
  colnames(elements.lengths) <- new.names
  
  elements.lengths <- cbind(metadata$Museum.ID ,metadata$Catalog.., metadata$Image, metadata$Binomial, as.data.frame(elements.lengths))
  # output <- aggregate(element.lengths, by = list(elements.lengths$`metadata$Catalog..`), FUN = mean)
  
  output <- elements.lengths
  
  #need to input way to triangulate/trig way to true length for femur/humerus
  
  
  return(output)
}

compile.Data <- function(elem.data = list(NA), by.merge = c("metadata$Binomial"))
{
  
  #just need to input the names of columns here and sort for thos not in those or the by input
  drop.vec <- c("metadata$Museum.ID", "metadata$Catalog..", "metadata$Image") #, "metadata$Binomial")
  
  drop.vec <- drop.vec[!drop.vec %in% by]
  
  elem.data <- lapply(X= elem.data, FUN = function(x) x[,!colnames(x) %in% drop.vec])
  
  for(xx in seq(1,length(elem.data),1))
  {
    elem.data[[xx]] <- aggregate(x = elem.data[[xx]], 
                                 by = list(Binomial = elem.data[[xx]]$`metadata$Binomial`), 
                                 FUN = mean, simplify = TRUE)
    #remove NA binomial column
    elem.data[[xx]] <- elem.data[[xx]][,-2]
  }
  #aggregate
  #remove flipped images
  
  dat.com <- merge(elem.data[[1]], elem.data[[2]], by = "Binomial", all=FALSE, no.dups = TRUE)
  if(length(elem.data) >=3)
  {
    for(xx in seq(3, length(elem.data),1))
    {
      dat.com <- merge(dat.com, elem.data[[xx]], by = "Binomial", all=FALSE, no.dups = TRUE)
    }
  }
  
  
  rownames(dat.com) <- dat.com[,1]; dat.com <- dat.com[,-1]
  
  return (dat.com)
}

get.bm <- function(this.rank = "species", focal.order = c("Artiodactyla", "Perissodactyla"))
{
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
  if (this.rank=="genus") measure.mat <- makeOneGenusMatFromSpecimenMat(measure.mat)
  
  bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
  bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
  shortFam <- sort(unique(bigList$family[bigList$order %in% focal.order]))	
  
  bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
  
  measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name[bigList$order %in% focal.order], ]
  
  
  return(measure.mat)
}

###############################################################################################################
###############################################################################################################

require(geomorph)
require(shapes)
require(stringr)

source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/occFns.R")
source("~/Dropbox/code/R/common_src/sampling.R") 
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R") 


#####
#11/11/2020 add in subfamily/tribe designation to see if there is some higher taxon trends
####can draw this from PBDB
#####

metafile <- "/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Equid Specimen Photo Metadata.csv"
meta.data <- read.csv(metafile, na.strings = c("",NA), stringsAsFactors = FALSE)
meta.data$Accepted.Species[is.na(meta.data$Accepted.Species)] <- "sp"


outputFolder <- "/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis"
#imageJ/Giji require conversion from text file to tps
#Femur
equid.femur.names <- convert2TPS(folder.direct ="/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis/Femur",
                                 filetype=".points", 
                                 meta.data = meta.data,
                                 element = "Femur",
                                 outputfilename = paste(outputFolder,"Equid_Femur.tps",sep= "/"))
#Tibia
equid.tibia.names <- convert2TPS(folder.direct ="/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis/Tibia",
                                   filetype=".points", 
                                   meta.data = meta.data,
                                   element = "Tibia",
                                   outputfilename = paste(outputFolder,"Equid_Tibia.tps",sep="/"))
#Metatarsal
equid.metatar.names <- convert2TPS(folder.direct ="/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis/Metatarsal",
                                   filetype=".points", 
                                   meta.data = meta.data,
                                   element = "Metatarsal",
                                   outputfilename = "Equid_Metatarsal.tps")
#Humerus
equid.humerus.names <- convert2TPS(folder.direct ="/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis/Humerus",
                                   filetype=".points", 
                                   meta.data = meta.data,
                                   element = "Humerus",
                                   outputfilename = "Equid_Humerus.tps")
#Radius
equid.radius.names <- convert2TPS(folder.direct ="/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis/Radius",
                                  filetype=".points", 
                                  meta.data = meta.data,
                                  element = "Radius",
                                  outputfilename = "Equid_Radius.tps")
#Metacarpal
equid.metacar.names <- convert2TPS(folder.direct ="/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis/Metacarpal",
                                   filetype=".points", 
                                   meta.data = meta.data,
                                   element = "Metacarpal",
                                   outputfilename = "Equid_Metacarpal.tps")

tpsFolder <- "/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Analysis/"

femur.tps <- readland.tps(paste(tpsFolder, "Equid_Femur.tps", sep=""), negNA = TRUE, specID = "imageID")
tibia.tps <- readland.tps(paste(tpsFolder, "Equid_Tibia.tps", sep=""), negNA = TRUE, specID = "imageID")
metatar.tps <- readland.tps(paste(tpsFolder, "Equid_Metatarsal.tps", sep=""), negNA = TRUE, specID = "imageID")
humerus.tps <- readland.tps(paste(tpsFolder, "Equid_Humerus.tps", sep=""), negNA = TRUE, specID = "imageID")
radius.tps <- readland.tps(paste(tpsFolder, "Equid_Radius.tps", sep=""), negNA = TRUE, specID = "imageID")
metacar.tps <- readland.tps(paste(tpsFolder, "Equid_Metacarpal.tps", sep=""), negNA = TRUE, specID = "imageID")

#grouping variables for factor and labeling plots
#make and import a .csv
  #need to extract image names so I can plot duplicates and tie in species
  
  #extract specimen names from tps file
  femur.specimens <- equid.femur.names
  tibia.specimens <- equid.tibia.names
  metatarsal.specimens <- equid.metatar.names
  humerus.specimens <- equid.humerus.names
  radius.specimens <- equid.radius.names
  metacarpal.specimens <- equid.metacar.names
  
  #rename pathname to specimen /image names
  dimnames(femur.tps)[[3]] <- rownames(femur.specimens)
  dimnames(tibia.tps)[[3]] <- rownames(tibia.specimens)
  dimnames(metatar.tps)[[3]] <- rownames(metatarsal.specimens)
  dimnames(humerus.tps)[[3]] <- rownames(humerus.specimens)
  dimnames(radius.tps)[[3]] <- rownames(radius.specimens)
  dimnames(metacar.tps)[[3]] <- rownames(metacarpal.specimens)
  
  #remove specimens
  #use specimens object to get index in landmark array
  species.rm <- c(NA) #c("Fortihippus_niobrarensis", "Merychippus_eohipparion", #"Kalobatippus_agatensis", 
                    #"Orohippus_sp", "Hyracotherium_sp", "Acritohippus_tertius", "Pseudhipparion_skinneri")
  
  femur.rm <- get.list.2remove(femur.specimens, spec.via.sp = species.rm)
  femur.rm <- append(femur.rm, "AMNH_99384_38")
  femur.tps <- removeSpecimen(landfile = femur.tps, spec2remove = femur.rm, LD = 6, k=2)
  femur.specimens <- femur.specimens[!rownames(femur.specimens) %in% femur.rm,]
  
  tibia.rm <- get.list.2remove(tibia.specimens, spec.via.sp = species.rm)
  tibia.tps <- removeSpecimen(landfile = tibia.tps, spec2remove = tibia.rm, LD = 5, k=2)
  tibia.specimens <- tibia.specimens[!rownames(tibia.specimens) %in% tibia.rm,]
  
  metatar.rm <- get.list.2remove(metatarsal.specimens, spec.via.sp = species.rm)
  metatar.rm <- append(metatar.rm, rownames(metatarsal.specimens[which(str_detect(rownames(metatarsal.specimens), pattern="Flip")),]))
  metatar.tps <- removeSpecimen(landfile = metatar.tps, spec2remove = metatar.rm, LD = 6, k=2)
  metatarsal.specimens <- metatarsal.specimens[!rownames(metatarsal.specimens) %in% metatar.rm,]
  
  humerus.rm <- get.list.2remove(humerus.specimens, spec.via.sp = species.rm)
  humerus.tps <- removeSpecimen(landfile = humerus.tps, spec2remove = humerus.rm, LD = 6, k=2)
  humerus.specimens <- humerus.specimens[!rownames(humerus.specimens) %in% humerus.rm,]
  
  radius.rm <- get.list.2remove(radius.specimens, spec.via.sp = species.rm)
  radius.tps <- removeSpecimen(landfile = radius.tps, spec2remove = radius.rm, LD = 6, k=2)
  radius.specimens <- radius.specimens[!rownames(radius.specimens) %in% radius.rm,]
  
  metacar.rm <- get.list.2remove(metacarpal.specimens, spec.via.sp = species.rm)
  metacar.rm <- append(metacar.rm, rownames(metacarpal.specimens[which(str_detect(rownames(metacarpal.specimens), pattern="Flip")),]))
  metacar.tps <- removeSpecimen(landfile = metacar.tps, spec2remove = metacar.rm, LD = 6, k=2)
  metacarpal.specimens <- metacarpal.specimens[!rownames(metacarpal.specimens) %in% metacar.rm,]
  
  #remove specimens
  
  #check missing landmarks
  femur.specimens <- checkMissingLand(femur.tps, femur.specimens)
  tibia.specimens <- checkMissingLand(tibia.tps, tibia.specimens)
  metatarsal.specimens <- checkMissingLand(metatar.tps, metatarsal.specimens)
  humerus.specimens <- checkMissingLand(humerus.tps, humerus.specimens)
  radius.specimens <- checkMissingLand(radius.tps, radius.specimens)
  metatarsal.specimens <- checkMissingLand(metatar.tps, metatarsal.specimens)
  metacarpal.specimens <- checkMissingLand(metacar.tps, metacarpal.specimens)
  
  #get.list of missing landmarks per species by element
  
  #rotate coords so they are on the same side
 #flip both X and Y
  
  align.femur <-  align2Side_new(landfile = femur.tps, metafile = femur.specimens, side = "Right")
  align.tibia <-  align2Side_new(landfile = tibia.tps, metafile = tibia.specimens, side = "Right")
  align.metatar <-  align2Side_new(landfile = metatar.tps, metafile = metatarsal.specimens, side = "Right")
  align.humerus <-  align2Side_new(landfile = humerus.tps, metafile = humerus.specimens, side = "Right")
  align.radius <-  align2Side_new(landfile = radius.tps, metafile = radius.specimens, side = "Right")
  align.metacar <-  align2Side_new(landfile = metacar.tps, metafile = metacarpal.specimens, side = "Right")
  
 femur.tps <- align.femur
 tibia.tps <- align.tibia
 metatar.tps <- align.metatar
 humerus.tps <- align.humerus
 radius.tps <- align.radius
 metacar.tps <- align.metacar

 #estimate or remove missing
  #estimate.missing on 3d array
    #specimens missing one end sometimes end up with extremely exaggerated engative values
  
 
  #remove landmarks 
  #femur.tps <-  removeLand(landfile = femur.tps, land2remove = c(2,8,9,12,13))
  
  # currently have too few specimens to estiamte (2 #missing+2=minimum # specimens needed to use this method reliably)
  #make a check that goes through the estim.missing() limits of 2m+2=#specimens needed and removes specimens that lack too many coord
  
  #make function to find which specimens are being estimated, mostly for double checking how they plot
  #femur.specimens <- checkMissingLand(landfile = femur.tps, specimens = femur.specimens)
 
 
  #When estimation not needed
  femur_ant.estim <- femur.tps
  tibia_ant.estim <- tibia.tps
  metatarsal_ant.estim <- metatar.tps
  humerus_ant.estim <- humerus.tps
  radius_ant.estim <- radius.tps
  metacarpal_ant.estim <- metacar.tps
 
 
 
  #estimate the remaining
  femur_ant.estim <- getEstimate(landmarks = femur.tps, meta.data = femur.specimens, method = "TPS")
  tibia_ant.estim <- getEstimate(landmarks = tibia.tps, meta.data = tibia.specimens, method = "TPS")
  metatarsal_ant.estim <- getEstimate(landmarks = metatar.tps, meta.data = metatarsal.specimens, method = "TPS")
  humerus_ant.estim <- getEstimate(landmarks = humerus.tps, meta.data = humerus.specimens, method = "TPS")
  radius_ant.estim <- getEstimate(landmarks = radius.tps, meta.data = radius.specimens, method = "TPS")
  metacarpal_ant.estim <- getEstimate(landmarks = metacar.tps, meta.data = metacarpal.specimens, method = "TPS")


  #values skewed heavily whenever missing values present in data.  worst in specimens conpletely lacking proximal head
getPCA(specimen.array = femur_ant.estim, colgroup = toupper(femur.specimens$Side), 
      pointLabel = rownames(femur.specimens), useLegend = TRUE, legend.pos.x = "topleft", 
      #pointLabel = femur.specimens$Binomial, useLegend = TRUE, legend.pos.x = "topleft", 
       legend.pos.y = NULL, scale.load = 0.10)

getPCA(specimen.array = tibia_ant.estim, colgroup = toupper(tibia.specimens$Side), 
       pointLabel = rownames(tibia.specimens),useLegend = FALSE, legend.pos.x = "topleft", 
      #pointLabel = tibia.specimens$Binomial,useLegend = FALSE, legend.pos.x = "topleft", 
       legend.pos.y = NULL,  scale.load = 0.10)

getPCA(specimen.array = metatarsal_ant.estim, colgroup = toupper(metatarsal.specimens$Side),
      pointLabel = rownames(metatarsal.specimens), useLegend = FALSE, legend.pos.x = "topleft", 
     #  pointLabel = metatarsal.specimens$Binomial, useLegend = TRUE, legend.pos.x = "topleft", 
       legend.pos.y = NULL,  scale.load = 0.05)

getPCA(specimen.array = humerus_ant.estim, colgroup = toupper(humerus.specimens$Side),
      pointLabel = rownames(humerus.specimens), useLegend = TRUE, legend.pos.x = "topleft", 
     # pointLabel = humerus.specimens$Binomial, useLegend = TRUE, legend.pos.x = "topleft", 
    legend.pos.y = NULL,  scale.load = 0.10)

getPCA(specimen.array = radius_ant.estim, colgroup = toupper(radius.specimens$Side),
       pointLabel = rownames(radius.specimens), useLegend = FALSE, legend.pos.x = "topleft", 
      # pointLabel = radius.specimens$Binomial, useLegend = TRUE, legend.pos.x = "topleft", 
       legend.pos.y = NULL,  scale.load = 0.10)

getPCA(specimen.array = metacarpal_ant.estim, colgroup = toupper(metacarpal.specimens$Side),
       pointLabel = rownames(metacarpal.specimens), useLegend = TRUE, legend.pos.x = "topleft", 
       #pointLabel = metacarpal.specimens$Binomial, useLegend = TRUE, legend.pos.x = "topleft", 
       legend.pos.y = NULL,  scale.load = 0.10)

#combine tps files to compare limbs together?
#need to discern how to get measurements from landmarks
#take length from both side and average (mets, radius,)  but other might require trig (humerus, femur, tibia)

# make this a function
####add option to procrustes then take measures...? (dont think this valid)
#####need to remember to accomidate size as factor


femur.lengths <- get.linear.measures(landmarks = femur_ant.estim, metadata = femur.specimens, 
                                     meas.array = c(1,2,3,6,1,5,1,4), 
                                     meas.names = c("fem.head.diameter","fem.dist.Width","fem.total.medial.Length","fem.total.lateral.Length"))

tibia.lengths <- get.linear.measures(landmarks = tibia_ant.estim, metadata = tibia.specimens, 
                             meas.array = c(1,5,2,4,1,2,5,3), 
                             meas.names = c("tib.prox.Width","tib.dist.Width","tib.total.medial.Length","tib.total.lateral.Length"))

metatar.lengths <- get.linear.measures(landmarks = metatarsal_ant.estim, metadata = metatarsal.specimens, 
                                       meas.array = c(1,6,2,5,3,4,1,3,4,6), 
                                       meas.names = c("mt.prox.Width","mt.distal.max.Width","mt.distal.troc.Width",
                                                     "mt.tot.medial.Length","mt.tot.lateral.Length"))

humerus.lengths <- get.linear.measures(landmarks = humerus_ant.estim, metadata = humerus.specimens, 
                                      meas.array = c(2,4,5,6,1,5,4,6), 
                                      meas.names = c("hum.prox.Width","hum.dist.Width","hum.total.medial.Length","hum.total.lateral.Length"))

radius.lengths <- get.linear.measures(landmarks = radius_ant.estim, metadata = radius.specimens, 
                                     meas.array = c(1,6,3,4,1,3,4,6), 
                                     meas.names = c("rad.prox.Width","rad.dist.Width","rad.total.medial.Length","rad.total.lateral.Length"))

metacar.lengths <- get.linear.measures(landmarks = metacarpal_ant.estim, metadata = metacarpal.specimens, 
                                       meas.array = c(1,6,2,5,3,4,1,3,4,6), 
                                       meas.names = c("mc.prox.Width","mc.distal.max.Width","mc.distal.troc.Width",
                                                      "mc.tot.medial.Length","mc.tot.lateral.Length"))
#do stepwise dicriminant
#see how a priori affects visualization

elem.data <- list(femur.lengths, tibia.lengths, metatar.lengths, humerus.lengths ,radius.lengths, metacar.lengths)

drop.vec <- c("metadata$Museum.ID", "metadata$Catalog..", "metadata$Image") #, "metadata$Binomial")

#need matrix to maintain specimen and taxon info for later reattachment
meta.elem <- matrix(ncol=4)
colnam.elem <- c(drop.vec, "metadata$Binomial"); colnames(meta.elem) <- colnam.elem

for(xx in seq(1, length(elem.data),1))
{
 # elem.data[[xx]]$Full.No...Id <- paste(paste(elem.data[[xx]]$'metadata$Museum.ID', 
  #                                    elem.data[[xx]]$'metadata$Catalog..', sep="_"),
   #                                   elem.data[[xx]]$'metadata$Image', sep = "_")
  
  elem.data[[xx]]$Full.No...Id <- paste(elem.data[[xx]]$'metadata$Museum.ID', 
                                              elem.data[[xx]]$'metadata$Catalog..', sep="_")
  
  meta.elem <- rbind(meta.elem, elem.data[[xx]][,colnames(elem.data[[xx]]) %in% colnam.elem])
  
  elem.data[[xx]] <- elem.data[[xx]][,!colnames(elem.data[[xx]]) %in% drop.vec]
}

for(xx in seq(1,length(elem.data),1))
{
  elem.data[[xx]] <- aggregate(x = elem.data[[xx]], 
                               by = list(Full.No...Id = elem.data[[xx]]$Full.No...Id), 
                               FUN = mean, simplify = TRUE)
  #remove NA binomial column
  elem.data[[xx]] <- elem.data[[xx]][,-c(2, ncol(elem.data[[xx]]))]
}

dat.com <- merge(elem.data[[1]], elem.data[[2]], by = "Full.No...Id", all= TRUE, no.dups = TRUE)
if(length(elem.data) >=3)
{
  for(xx in seq(3, length(elem.data),1))
  {
    dat.com <- merge(dat.com, elem.data[[xx]], by = "Full.No...Id", all=TRUE, no.dups = TRUE)
  }
}

#reattach taxonomic info
rownames(meta.elem) <- NULL
meta.elem <- unique(meta.elem[,c(1:2,4)])
meta.elem$Full.No...Id  <- paste(meta.elem$`metadata$Museum.ID`, meta.elem$`metadata$Catalog..`,sep="_")

dat.com$Binomial <- NA

for(xx in seq(1, nrow(dat.com),1))
{
  #use for loop to query through both objects and infill taxon name into dat.com
  for(yy in seq(1, nrow(meta.elem),1))
  {
    if(dat.com[xx,]$Full.No...Id %in% meta.elem[yy,]$Full.No...Id)
    {
      dat.com[xx,]$Binomial <- meta.elem[yy,]$`metadata$Binomial`
    }
  }
}

#purge measured onject of extra columns
#all.elem <- compile.Data(elem.data = elem.data, lit = lit.specimens, by =c("Binomial"))

all.elem <- dat.com

pull.colnames <- c("fem.total.medial.Length","fem.total.lateral.Length","tib.total.medial.Length","tib.total.lateral.Length",
                   "mt.tot.medial.Length","mt.tot.lateral.Length", "hum.total.medial.Length","hum.total.lateral.Length",
                   "hum.total.medial.Length","hum.total.lateral.Length", "rad.total.medial.Length","rad.total.lateral.Length",
                   "mc.tot.medial.Length","mc.tot.lateral.Length",
                   "mc.distal.max.Width", "mt.distal.max.Width", "")
all.elem <- all.elem[,!colnames(all.elem) %in% pull.colnames]

#read in literature measures
lit.specimens <- read.csv("/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Equid Literature Measurements 2020_12_15.csv")
rownames(lit.specimens) <- lit.specimens$Full.No...Id
lit.specimens <- lit.specimens[,-c(1:4,6:7,9:12,31)]#;lit.specimens <- lit.specimens[,-c(2:8,27)]
colnames(lit.specimens)[2] <- "Binomial"

#reorder to match lit entries
attributes(lit.specimens)$names
attributes(all.elem)$names
col.order <- c("Full.No...Id","Binomial",
               "fem.head.diameter", "fem.dist.Width","fem.tot.length.Ave",
               "tib.prox.Width","tib.dist.Width","tib.tot.length.Ave",
               "mt.prox.Width","mt.distal.troc.Width","mt.tot.length.Ave",
               "hum.prox.Width","hum.dist.Width","hum.tot.length.Ave",
               "rad.prox.Width","rad.dist.Width","rad.tot.length.Ave",
               "mc.prox.Width","mc.distal.troc.Width","mc.tot.length.Ave")

all.elem <- all.elem[,col.order]

all.elem <- rbind(all.elem, lit.specimens)

#aggregate by species
agg.elem <- aggregate(x = all.elem, by = list(Binomial = all.elem$Binomial), 
                               FUN = mean, simplify = TRUE, na.action = NULL, na.rm=TRUE)
#remove NA binomial column
agg.elem <- agg.elem[,-c(2:3)]

write.csv(agg.elem, "/Users/emdoughty/Documents/Code/R/LimbProportions/Landmarks/Complete Elements/Equid Compiled Species Measures 2021_1_19.csv")


rownames(agg.elem) <- agg.elem$Binomial
agg.elem <- agg.elem[,-1]

log.elem <- log10(agg.elem)

#remove columns that wont be useful
###remove prox femur depth=not consistent with literature

#columns
col.rem <- c("fem.head.diameter","fem.dist.Width",#"fem.tot.length.Ave",
             "tib.prox.Width","tib.dist.Width", #"tib.tot.length.Ave",
             "mt.prox.Width","mt.distal.troc.Width", #"mt.tot.length.Ave",
             "hum.prox.Width","hum.dist.Width", #"hum.tot.length.Ave",
             "rad.prox.Width","rad.dist.Width",#"rad.tot.length.Ave",
              "mc.prox.Width","mc.distal.troc.Width"#,"mc.tot.length.Ave"
           )

pca.elem <- log.elem[,!colnames(log.elem) %in% col.rem]
pca.elem <- pca.elem[complete.cases(pca.elem),]

this.x <- 2
this.y <- 3

par(mfrow = c(1,2))
getPCA.linear(specimen.array = pca.elem, colgroup = rownames(pca.elem), 
              pointLabel = rownames(pca.elem), this.x = this.x, this.y = this.y, useLegend = FALSE, legend.pos.x = "topright", 
              legend.pos.y = NULL,  loadings=FALSE,scale.load = 0)
getPCA.linear(specimen.array = pca.elem, colgroup = rownames(pca.elem), 
              pointLabel = rownames(pca.elem), this.x = this.x, this.y = this.y, useLegend = FALSE, legend.pos.x = "topright", 
              legend.pos.y = NULL,  scale.load = 0.1, load.cex=0.5)

#will need to set as lenths only to use most of literature

#all.elem$mt.fem.ratio <- all.elem$mt.tot.length.Ave/all.elem$fem.tot.length.Ave

#plot limb props through time

#avoid logging of ratios
pca.elem <- agg.elem[,!colnames(agg.elem) %in% col.rem]
pca.elem <- pca.elem[complete.cases(pca.elem),]

mt.femur <- pca.elem$mt.tot.length.Ave/pca.elem$fem.tot.length.Ave
mt.tibia <- pca.elem$mt.tot.length.Ave/pca.elem$tib.tot.length.Ave
tib.femur <- pca.elem$tib.tot.length.Ave/pca.elem$fem.tot.length.Ave

mc.hum <- pca.elem$mc.tot.length.Ave/pca.elem$hum.tot.length.Ave
mc.rad <- pca.elem$mc.tot.length.Ave/pca.elem$rad.tot.length.Ave
rad.hum <- pca.elem$rad.tot.length.Ave/pca.elem$hum.tot.length.Ave

limb.ratios <- cbind(mt.femur, mt.tibia, tib.femur, mc.hum, mc.rad, rad.hum)
rownames(limb.ratios) <- rownames(pca.elem)

#get durations from pbdb
occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- occs[occs$order %in% "Perissodactyla", ]
occs <- occs[occs$family %in% c("Equidae","Palaeotheriidae"), ]
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

occs <- occs[,c("accepted_name","genus", "max_ma","min_ma")]

occs.species.unique <- unique(occs$accepted_name)
occs.genus.unique <- unique(occs$genus)
occs.species.unique <- occs.species.unique[!occs.species.unique %in% occs.genus.unique[!occs.genus.unique %in% c("Hyracotherium","Orohippus")]]
#occs.genus.unique <- paste(occs.genus.unique, "sp.", sep="_")

#occs.species.unique <- append(occs.species.unique, occs.genus.unique)

sp.duration <- matrix(nrow = length(occs.species.unique),ncol = 3)

check.dates <- list()

for(xx in seq(1,length(occs.species.unique),1))
{
  occs.indiv <- occs[occs$accepted_name %in% occs.species.unique[xx],]
 
 # check.dates[[xx]] <- occs.indiv
  
  sp.duration[xx,] <- c(unique(occs.indiv$accepted_name),max(occs.indiv$max_ma), min(occs.indiv$min_ma))
  #print(xx)
}

sp.duration[sp.duration %in% c("Hyracotherium","Orohippus")] <- paste(c("Hyracotherium","Orohippus"),"sp",sep="_")

sp.duration.names <- sp.duration[,1]; sp.duration <- sp.duration[,-1]
sp.duration <- apply(sp.duration, 2, as.numeric)

rownames(sp.duration) <- sp.duration.names; colnames(sp.duration) <- c("max_ma","min_ma")

#attach dates to ratios/limb measures
limb.dated <- limb.ratios

#fix the dates for individual species
##Hipparion shirleyae : only one entry present in occs

##Merychippus coloradense
#check.dates[[93]]
##23.03 to 5.33
##15.97 to 13.6
##13.60 to 10.3
##13.60 to 4.9
##13.60 t0 2.588

###will set to 15.97 to 10.3 for now
sp.duration[rownames(sp.duration) %in% "Merychippus_coloradense",1] <- 15.97
sp.duration[rownames(sp.duration) %in% "Merychippus_coloradense",2] <- 10.3

##Nephipparion affine
#check.dates[[83]]
##23.3 to 5.33
##13.60 to 10.3
##13.60 to 4.9

###will set to 213.60to 10.3 for now
sp.duration[rownames(sp.duration) %in% "Neohipparion_affine",1] <- 13.60
sp.duration[rownames(sp.duration) %in% "Neohipparion_affine",2] <- 10.3

##Parahippus leonensis
#check.dates[[64]]
##20.43 to 15.970
##30.80 to 20.430
##23.03 to 5.333
##one listed as 15.97 to 13.6

###will set to 20.43 to 13.6 for now
sp.duration[rownames(sp.duration) %in% "Parahippus_leonensis",1] <- 20.43
sp.duration[rownames(sp.duration) %in% "Parahippus_leonensis",2] <- 13.6

##Protohippus supremus
#check.dates[[86]]
##23.03 to 5.33
##11.608 to 5.333
##13.60 to 10.3
##13.60 to 4.9

###will set to 13.60to 10.3 for now
sp.duration[rownames(sp.duration) %in% "Protohippus_supremus",1] <- 13.60
sp.duration[rownames(sp.duration) %in% "Protohippus_supremus",2] <- 10.3

##Psuedohipparion skinneri
#check.dates[[102]]
##23.03 to 5.33
##11.608 to 5.333
##10.30 to 4.9
##15.97 to 11.608
##13.60 to 10.3
##13.60 to 4.9

###will set to 13.60to 10.3 for now
sp.duration[rownames(sp.duration) %in% "Pseudhipparion_skinneri",1] <- 15.97
sp.duration[rownames(sp.duration) %in% "Pseudhipparion_skinneri",2] <- 10.3


sp.dur.reduce <- sp.duration[rownames(sp.duration) %in% rownames(limb.ratios),]

limb.dated <- merge(limb.dated, sp.duration, by = 0, all.x= TRUE)
rownames(limb.dated) <- limb.dated$Row.names; limb.dated <- limb.dated[,-1]      

#enter dates for modern/non-NA Equids
limb.dated[rownames(limb.dated) == 'Equus_caballus', 7:8] <- c(3.6,0.0)

habitat.Janis2012 <- read.csv("/Users/emdoughty/Downloads/Janis etal 2012 dataset.csv")
habitat.equids <- habitat.Janis2012[str_detect(habitat.Janis2012$Order.Family,c("Equidae")),]

habitat.equids$Species <- gsub(pattern = "[[:space:]]", replacement = "_", x = habitat.equids$Species)
habitat.equids$Species <- gsub(pattern = "[.]", replacement = "", x = habitat.equids$Species)
rownames(habitat.equids) <- habitat.equids$Species; habitat.equids <- habitat.equids[,-c(1:4,7:11)]

habitat.equids$genus <- gsub("([A-Za-z]+).*", "\\1", rownames(habitat.equids))
hab.gen.uniq <- unique(habitat.equids)

limb.dated <- merge(limb.dated, habitat.equids, by = 0, all.x= TRUE)
rownames(limb.dated) <- limb.dated$Row.names; limb.dated <- limb.dated[,-1]

limb.dated$genus <- gsub("([A-Za-z]+).*", "\\1", rownames(limb.dated))

for(xx in which(is.na(limb.dated$Habitat)))
{
  if(limb.dated$genus[xx] %in% hab.gen.uniq$genus)
  {
    limb.dated[xx,colnames(limb.dated) %in% c("Order.Family", "Habitat")] <- 
      hab.gen.uniq[hab.gen.uniq$genus %in% limb.dated$genus[xx],c(1:2)]
  }
}

limb.dated$Habitat <- as.character(limb.dated$Habitat)

limb.dated[is.na(limb.dated$Habitat),10] <- "Unknown"

#############################################################################################################################################
#plot ratios
ratios.plot <- limb.dated$mt.femur
ratios.label <- "mt/femer"
line.col <- "gray0"
line.trans <- 0.5
#line.col <- set to code for bodymass

plot(limb.dated$max_ma, ratios.plot, xlim=c(55,0), type="n", xaxp =c(55,5,10), xlab="Time (Ma)", ylim=c(min(ratios.plot),1.025), ylab=ratios.label, cex.axis=1.5, cex.lab=1.5)
overlayCzTimescale(do.subepochs=TRUE)
rect(xleft=24,ybottom = 0 ,xright =17,ytop = 4, col=alphaColor("green",0.5))
for (i in seq_len(nrow(limb.dated))) {
  if (is.finite(limb.dated$max_ma[i]) & is.finite(limb.dated$min_ma[i]) & limb.dated$max_ma[i] != limb.dated$min_ma[i]) lines(x=limb.dated[i,c("max_ma","min_ma")], y=c(ratios.plot[i], ratios.plot[i]), lwd=1.5, pch=21, col=alphaColor(line.col, line.trans)) #alphaColor(orderColors[i], 0.5)
  if (is.finite(limb.dated$max_ma[i]) & is.finite(limb.dated$min_ma[i]) & limb.dated$max_ma[i] != limb.dated$min_ma[i]) text(x=(limb.dated[i,"min_ma"]+limb.dated[i,"max_ma"])/2, y=ratios.plot[i]+0.002, rownames(limb.dated)[i], cex=0.5, col=alphaColor("black", 0.5)) #alphaColor(orderColors[i], 0.5)
  if (is.finite(limb.dated$max_ma[i]) & is.finite(limb.dated$min_ma[i]) & limb.dated$max_ma[i] != limb.dated$min_ma[i]) text(x=(limb.dated[i,"min_ma"]+limb.dated[i,"max_ma"])/2, y=ratios.plot[i]+0.002, rownames(limb.dated)[i], cex=0.5, col=alphaColor("black", 0.5)) #alphaColor(orderColors[i], 0.5)
}

########################################################################################################
#compare to body mass through time of equidae

###bring in body mass from ecology or via backup files

#run body mass generation in ecology/evo code
measure.mat <- get.bm()
measure.gen <- get.bm(this.rank = "genus")
rownames(measure.gen) <- paste(rownames(measure.gen),"sp",sep="_")
measure.mat <- rbind(measure.mat, measure.gen)
measure.mat[measure.mat$taxon %in% c("Hyracotherium", "Orohippus"),]$family <- "Equidae"

measure.mat[str_detect(measure.mat$taxon,c("Orohippus")),] 
measure.mat[str_detect(measure.mat$taxon,c("Hyracotherium")),] 

equid.bm <- measure.mat[measure.mat$family %in% c("Equidae","Palaeotheriidae"),]

equid.bm.dated <- merge(equid.bm,sp.duration,by=0,all.x=TRUE)
rownames(equid.bm.dated) <- equid.bm.dated$Row.names; equid.bm.dated <- equid.bm.dated[,-c(1:34)]
equid.bm.dated <- equid.bm.dated[complete.cases(equid.bm.dated),]

plot(equid.bm.dated$max_ma, equid.bm.dated$bodyMass, xlim=c(55,0), type="n", xaxp =c(55,5,10), xlab="Time (Ma)", ylim=c(0,4), ylab="log Body Mass (kg)", cex.axis=1.5, cex.lab=1.5)
overlayCzTimescale(do.subepochs=TRUE)
rect(xleft=24,ybottom = 0 ,xright =17,ytop = 4, col=alphaColor("green",0.5))
for (i in seq_len(nrow(equid.bm.dated))) {
  if (is.finite(equid.bm.dated$max_ma[i]) & is.finite(equid.bm.dated$min_ma[i]) & equid.bm.dated$max_ma[i] != equid.bm.dated$min_ma[i]) lines(x=equid.bm.dated[i,c("max_ma","min_ma")], y=c(equid.bm.dated$bodyMass[i], equid.bm.dated$bodyMass[i]), lwd=2, pch=21, col=alphaColor("gray0", 0.5)) #alphaColor(orderColors[i], 0.5)
  if (is.finite(limb.bm$max_ma[i]) & is.finite(limb.bm$min_ma[i]) & limb.bm$max_ma[i] != limb.bm$min_ma[i]) lines(x=limb.bm[i,c("max_ma","min_ma")], y=c(limb.bm$bodyMass[i], limb.bm$bodyMass[i]), lwd=2, pch=21, col=alphaColor("red", 1)) #alphaColor(orderColors[i], 0.5)
  #if (is.finite(equid.bm.dated$max_ma[i]) & is.finite(equid.bm.dated$min_ma[i]) & equid.bm.dated$max_ma[i] != equid.bm.dated$min_ma[i]) text(x=(equid.bm.dated[i,"min_ma"]+equid.bm.dated[i,"max_ma"])/2, y=equid.bm.dated$bodyMass[i]+0.002, rownames(equid.bm.dated)[i], cex=0.5, col=alphaColor("black", 0.5)) #alphaColor(orderColors[i], 0.5)
  #if (is.finite(equid.bm.dated$max_ma[i]) & is.finite(equid.bm.dated$min_ma[i]) & equid.bm.dated$max_ma[i] != equid.bm.dated$min_ma[i]) text(x=(equid.bm.dated[i,"min_ma"]+equid.bm.dated[i,"max_ma"])/2, y=equid.bm.dated$bodyMass[i]+0.002, rownames(equid.bm.dated)[i], cex=0.5, col=alphaColor("black", 0.5)) #alphaColor(orderColors[i], 0.5)
}

#plot ratios with body mass
equid.bm[rownames(equid.bm) %in% rownames(limb.dated),]
limb.bm <- merge(limb.dated,equid.bm,by=0, all.x=TRUE); limb.bm <- limb.bm[,-c(10:42)]
rownames(limb.bm) <- limb.bm$Row.names; limb.bm <- limb.bm[,-1]

#body mass for Equus caballus
limb.bm[rownames(limb.bm) %in% "Equus_caballus","bodyMass"] <- log10(mean(c(227,900))) #https://animaldiversity.org/accounts/Equus_caballus/

limb.bm <- cbind(limb.bm,limb.dated$Habitat); colnames(limb.bm)[ncol(limb.bm)] <- "Habitat"

col.hab <- limb.bm$Habitat
col.hab <- gsub("Unknown", "black", col.hab)
col.hab <- gsub("Closed", "brown", col.hab)
col.hab <- gsub("Mixed", "green", col.hab)
col.hab <- gsub("Open", "yellow", col.hab)


quartz()
par(mfrow=c(2,3))
plot(limb.bm$bodyMass, limb.bm$mt.femur, pch=16, col = col.hab)
text(limb.bm$bodyMass, limb.bm$mt.femur, rownames(limb.bm), adj=-0.1,cex=0.5)
legend("topleft",legend = unique(limb.bm$Habitat), pch=16, col = unique(col.hab))

plot(limb.bm$bodyMass, limb.bm$mt.tibia, pch=16, col = col.hab)
text(limb.bm$bodyMass, limb.bm$mt.tibia, rownames(limb.bm), adj=-0.1,cex=0.5)

plot(limb.bm$bodyMass, limb.bm$tib.femur, pch=16, col = col.hab)
text(limb.bm$bodyMass, limb.bm$tib.femur, rownames(limb.bm), adj=-0.1,cex=0.5)

plot(limb.bm$bodyMass, limb.bm$mc.hum, pch=16, col = col.hab)
text(limb.bm$bodyMass, limb.bm$mc.hum, rownames(limb.bm), adj=-0.1,cex=0.5)

plot(limb.bm$bodyMass, limb.bm$mc.rad, pch=16, col = col.hab)
text(limb.bm$bodyMass, limb.bm$mc.rad, rownames(limb.bm), adj=-0.1,cex=0.5)

plot(limb.bm$bodyMass, limb.bm$rad.hum, pch=16, col = col.hab)
text(limb.bm$bodyMass, limb.bm$rad.hum, rownames(limb.bm), adj=-0.1,cex=0.5)

#total length of fore vs hindlimb
#ratios of forelibm and hindlimb equivelent elements

#compiled lengths
comp.length <- agg.elem[,!colnames(agg.elem) %in% col.rem]
comp.length <- comp.length[complete.cases(comp.length),]

comp.length <- merge(comp.length,limb.bm,by=0, all.x = TRUE)
rownames(comp.length) <- comp.length$Row.names; comp.length <- comp.length[,-1]

comp.length$hindlimb.tot <- comp.length$fem.tot.length.Ave+comp.length$tib.tot.length.Ave+comp.length$mt.tot.length.Ave
comp.length$forelimb.tot <- comp.length$hum.tot.length.Ave+comp.length$rad.tot.length.Ave+comp.length$mc.tot.length.Ave

comp.length$foreVhind <- comp.length$forelimb.tot/comp.length$hindlimb.tot
comp.length$humVfem <- comp.length$hum.tot.length.Ave/comp.length$fem.tot.length.Ave
comp.length$radVtib <- comp.length$rad.tot.length.Ave/comp.length$tib.tot.length.Ave
comp.length$mcVmt <- comp.length$mc.tot.length.Ave/comp.length$mt.tot.length.Ave

#can do fore vs hind but need to acocunt for body mass in a way...would 3d plots work?
quartz()
par(mfrow = c(2,4))
plot(comp.length$hindlimb.tot, comp.length$forelimb.tot, pch = 16, col = col.hab)
text(comp.length$hindlimb.tot, comp.length$forelimb.tot, rownames(comp.length), adj=-0.1,cex=0.5)
legend("topleft",legend = unique(comp.length$Habitat), pch=16, col = unique(col.hab))

plot(comp.length$bodyMass, comp.length$forelimb.tot, pch = 16, col = col.hab)
text(comp.length$bodyMass, comp.length$forelimb.tot, rownames(comp.length), adj=-0.1,cex=0.5)

plot(comp.length$bodyMass, comp.length$hindlimb.tot, pch = 16, col = col.hab)
text(comp.length$bodyMass, comp.length$hindlimb.tot, rownames(comp.length), adj=-0.1,cex=0.5)

plot(comp.length$bodyMass, comp.length$foreVhind, pch = 16, col = col.hab)
text(comp.length$bodyMass, comp.length$foreVhind, rownames(comp.length), adj=-0.1,cex=0.5)

plot(comp.length$bodyMass, comp.length$humVfem, pch = 16, col = col.hab)
text(comp.length$bodyMass, comp.length$humVfem, rownames(comp.length), adj=-0.1,cex=0.5)

plot(comp.length$bodyMass, comp.length$radVtib, pch = 16, col = col.hab)
text(comp.length$bodyMass, comp.length$radVtib, rownames(comp.length), adj=-0.1,cex=0.5)

plot(comp.length$bodyMass, comp.length$mcVmt, pch = 16, col = col.hab)
text(comp.length$bodyMass, comp.length$mcVmt, rownames(comp.length), adj=-0.1,cex=0.5)

plot(comp.length$bodyMass, comp.length$mcVmt, pch = 16, col = col.hab)
#run tree code below
phylomorphospace(tree_resolv,X=reduc.comp.length[,colnames(reduc.comp.length) %in% c("bodyMass", "mcVmt")],
                 label="off", node.size=c(0.5,1), add=TRUE)
points(comp.length$bodyMass, comp.length$mcVmt, pch = 16, col = col.hab)
text(comp.length$bodyMass, comp.length$mcVmt, rownames(comp.length), adj=-0.1,cex=0.5)

require(phytools)
source("/Users/emdoughty/Dropbox/Code/R/DentalMeasurements/src/src_EvAnalysesTree.R")
source("/Users/emdoughty/Dropbox/Code/R/common_src/occFns.R", chdir = TRUE)

tree.backbone <- read.nexus("/Users/emdoughty/Dropbox/Code/R/DentalMeasurements/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
tree_resolv <- drop.tip(tree_resolv, tip=tree_resolv$tip.label[!tree_resolv$tip.label %in% rownames(comp.length)])
ranges <- sp.duration; colnames(ranges) <- c("FO","LO")
tree_resolv <- dateTreeWithRanges(this.tree = tree_resolv, strat.ranges=ranges, within.error=FALSE, min.bl=1, smooth.node.ages=TRUE, root.max=85) # internals will be flexed by extendNodeAges below

plot(tree_resolv, cex=0.25)
reduc.comp.length <- comp.length[rownames(comp.length) %in% tree_resolv$tip.label,]

length(tree_resolv$tip.label)
nrow(reduc.comp.length)
reduc.comp.length[!rownames(reduc.comp.length) %in% tree_resolv$tip.label,]


#symbols like janis
jan.sym <- habitat.Janis2012$Order.Family
jan.sym <-gsub(c("Canidae","Felidae","Ursidae"),8, jan.sym)
jan.sym <-gsub(c("Felidae"),8, jan.sym)
jan.sym <-gsub(c("Ursidae"),8, jan.sym)

jan.sym <-gsub(c("Tragulidae"),2, jan.sym)
jan.sym <-gsub(c("Cervidae"),2, jan.sym)
jan.sym <-gsub(c("Bovidae"),2, jan.sym)
jan.sym <-gsub(c("Giraffidae"),2, jan.sym)
jan.sym <-gsub(c("Antilocapridae"),2, jan.sym)
jan.sym <-gsub(c("Dromomerycidae"),2, jan.sym)

jan.sym <-gsub(c("Hippopotamidae"),19, jan.sym)
jan.sym <-gsub(c("Suidae"),19, jan.sym)

jan.sym <-gsub(c("Merycoidodontidae"),0, jan.sym)
jan.sym <-gsub(c("Camelidae"),0, jan.sym)

jan.sym <-gsub(c("Palaeotheriidae"),18, jan.sym)
jan.sym <-gsub(c("Rhinocerotidae"),18, jan.sym)
jan.sym <-gsub(c("Hyracodontidae"),18, jan.sym)

jan.sym <-gsub(c("Equidae (H)"),11, jan.sym)
jan.sym <-gsub(c("Equidae (A)"),11, jan.sym)
jan.sym <-gsub(c("Equidae (E)"),11, jan.sym)

jan.sym <-gsub(c("Chalicotheriidae"),4, jan.sym)
jan.sym <-gsub(c("Brontotheriidae"),4, jan.sym)
jan.sym <-gsub(c("Tapiroidea"),4, jan.sym)
jan.sym <-gsub(c("Tapiridae"),4, jan.sym)

jan.sym <-gsub(c("Mammutidae"),6, jan.sym)
jan.sym <-gsub(c("Gomphotheridae"),6, jan.sym)
jan.sym <-gsub(c("Elephantidae"),6, jan.sym)

#color by family
jan.col.uniq <- rainbow(length(unique(as.character(habitat.Janis2012$Order.Family))))
names(jan.col.uniq) <- unique(as.character(habitat.Janis2012$Order.Family))

jan.col.list <- merge(habitat.Janis2012,jan.col.uniq,by.x = 0, by.y= ,all.x=TRUE)

#color by habitat
par(mfrow=c(1,2))
plot(habitat.Janis2012$AFD,habitat.Janis2012$MTR..LTR, pch=as.numeric(jan.sym),
     col = habitat.Janis2012$Order.Family,
     ylim = c(-10,50))
plot.new()
legend("topleft",legend = unique(habitat.Janis2012$Order.Family), cex=0.5)


#################################################################################
#get taxon richness through time
  
bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFam <- sort(unique(bigList$family[bigList$order %in% focal.order]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)

  par(mar=c(4,4, 2.5,0.5))
  
  # taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% measure.mat$taxon[x]], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  #taxCubeG <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$genus) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
  
  # prop <- t(apply(taxCube, c(1,2), median, na.rm=TRUE))
  prop <- t(apply(taxCube, c(1,2), mean, na.rm=TRUE))
  colnames(prop)[colnames(prop)==""] <- "indeterminate"
  # dimnames(prop) <- list(rownames(intervals), shortFam)
  source("https://dl.dropbox.com/s/iy0tu983xesbig2/taxonomicEv.R")
  plotStackedRichness(this.box=prop, intervals=intervals, do.log=FALSE, overlay.labels=TRUE, numbers.only=TRUE, legend=TRUE, xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)))
  axis(side=1,at = seq(0,60,5), lwd = 2, tck = -0.03)
  axis(side=1,at = seq(1,60,1), label =FALSE)
  #med.n <- median(length(unique(unlist(sapply(repIntTaxa[[this.rep]], function(x) measure.mat$taxon[x]))))) #what is this.rep set to during this function?  variable is used in for loop in Handley
  # med.n <- median(sapply(repIntTaxa, function(x) length(unique(unlist(sapply(x, function(y) measure.mat$taxon[y]))))))
  # optList_tax <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
  #abline(v=sort(c(intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col="darkorchid4")
 # text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_tax_median[[length(optList_tax_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="darkorchid4")
  #text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0), cex=0.5, col="darkorchid4")
  box(lwd=1)
  


