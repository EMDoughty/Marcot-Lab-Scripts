sourceCall <- function(){
  
  require(phytools)
  require(mvMORPH)
  require(stringr)
  require(parallel)
  
  #sources for Jon Marcot's code and specimen measurements 
  source("https://dl.dropbox.com/s/8jy9de5owxj72p7/strat.R")
  source("https://dl.dropbox.com/s/253p4avcvb66795/occFns.R")
  source("https://dl.dropbox.com/s/9gdafsqss2b586x/phy_dateTree.R")
  source("https://dl.dropbox.com/s/9tdawj35qf502jj/amandaSrc.R")
  source("https://dl.dropbox.com/s/rlof7juwr2q4y77/blasto_Birlenbach.R")
  source("https://dl.dropbox.com/s/643op7ye4s49w8p/utils_marcot.R")
  source("https://dl.dropbox.com/s/pxnbmroe2hgcgo3/kozak_src.R")
  if(Sys.info()["sysname"] == "Darwin"){
    ################MAC
    source('~/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call cource file for analysis functions
    source('~/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R') #call cource file for data functions
    source('~/Dropbox/ungulate_RA/RCode/EvAnalysesTreeSrc.R') #call cource file for tree functions
    source('~/Dropbox/ungulate_RA/RCode/EvAnalysesPlotSrc.R') #call cource file for tree functions
    
  } else if(Sys.info()["sysname"] == "Windows"){
    #################PC   #had to swap my name with Blaires in pathnames
    source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call source for use on Evan's PC
    source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R') #call cource file for data functions
    source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysesTreeSrc.R') #call cource file for tree functions
    source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysesPlotSrc.R') #call cource file for tree functions
    
  } else print("Mac or Windows operating systems are not detected")
}

getTreeFiles <- function()
{
  # have if else statements for different computers (PC vs Mac)
  if(Sys.info()["sysname"] == "Darwin"){
    tree.backbone <- read.nexus("~/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
    clade.definitions <- read.csv("~/Dropbox/ungulate_RA/2017_3_20_Clade_species_test.csv", stringsAsFactors = FALSE)
    wildcard.positions <- read.csv("~/Dropbox/ungulate_RA/2017_4_17_MCRA_Codes.csv", stringsAsFactors = FALSE)
    regressCat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")
    setwd('~/Dropbox/ungulate_RA/RCode/Results')
  } else if(Sys.info()["sysname"] == "Windows"){
    #Need way to replace "Van Valkenburgh Lab" in automated fashion
    tree.backbone <- read.nexus("C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
    clade.definitions <- read.csv("C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/2017_3_20_Clade_species_test.csv", stringsAsFactors = FALSE)
    wildcard.positions <- read.csv("C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/2017_4_17_MCRA_Codes.csv", stringsAsFactors = FALSE)
    regressCat <- read.csv("C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")
    setwd('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/Results')
  }
  
  tree.backbone$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree.backbone$tip.label, fixed=TRUE)))
  tree.list.wildcards <- getWildcardTrees(wildcard.positions)
  
  #return list of files that will be extracted later
  return(list(tree.backbone = tree.backbone, 
              clade.definitions = clade.definitions, 
              wildcard.positions = wildcard.positions, 
              regressCat = regressCat, 
              tree.list.wildcards = tree.list.wildcards))
}

getWildcardTrees <- function(wildcard.positions)
{
  tree.list.wildcards <- list()
  for (this.wildcard in seq_len(nrow(wildcard.positions))) {
    tree.list.wildcards[[this.wildcard]] <- read.nexus(file = wildcard.positions$Filename[this.wildcard])
    tree.list.wildcards[[this.wildcard]]$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree.list.wildcards[[this.wildcard]]$tip.label, fixed=TRUE)))
    #cat(this.wildcard, ":", length(tree.list.wildcards[[this.wildcard]]$tip.label),"\n")
  }
  class(tree.list.wildcards) <- "multiPhylo"
  
  return(tree.list.wildcards)
}

getOccs <- function()
{
  occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
  occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
  occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
  #occs <- occs[!occs$order %in% c("Desmostylia", "Perissodactyla"), ]
  occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores
  
  print("Completed Occs")
  return(occs)
}

getMeasureMat <- function(occs = occs, regMat = regressCat)
{
  measureMat <- getSingleSpeciesMatrix()
  print("getSingleSpMat completed")
  rownames(measureMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = measureMat$species)
  measureMat <- appendRegressionCategories(thisMat= measureMat, regMat= regressCat)
  print("regression append completed")
  #head(measureMat)
  measureMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = measureMat$species)
  #Approxiate body masses
  measureMat <-  approxBodyMass(thisMat = measureMat)
  print("bodymass approximated")
  
  return(measureMat)
}

getDat.vec <- function(tree.list, measureMat)
{
  measureMatCrop <- matrix()
  dat.vec <- data.frame()
  if (reps == 1) { measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[1]], thisMat = measureMat)
  tree.list[[1]] <- comp_TreeSpecMat_outTree(this.tree = tree.list[[1]], thisMat = measureMat)
  } else { 
    for (ii in seq(1, length(tree.list))) { 
      measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[ii]], thisMat = measureMat) 
      tree.list[[ii]] <- comp_TreeSpecMat_outTree(this.tree = tree.list[[ii]], thisMat = measureMat)
    }
  }
  dat.vec <- as.data.frame(data.matrix(frame = measureMatCrop[,"bodyMass"]))
  dimnames(dat.vec) <- list(rownames(measureMatCrop), "body.mass")
  names(dat.vec) <- "body.mass"
  
  print("Dat vector completed")
  return(dat.vec)
}

getOptModelValues <- function(fileList = NULL, ExtractAttributes = c("sigma","break.dates"))
{
  optBM <- list()
  for(ii in seq(1, length(filenames.2Ma),1)){ #rerun and save data file for sigma and separate for break.dates
    load(as.character(filenames.2Ma[[ii]]))
    #	load(as.character(filenames.Trees[[ii]]))

    bmBest <- getOptModels(opt.list = resultsTemporalShifts, model = "BM")
    optBM[[ii]] <-list(BM=bmBest[ExtractAttributes], OU=NULL)
    
    if(ii %% 10) print(ii)
  } 
  
  return(optBM)
}