install.packages('ape')
install.packages('phytools')
install.packages('paleotree')
#install.packages('stringr')

require(ape)
require(phytools)
require(paleotree)
#require(stringr)

source('~/Dropbox/ungulate_RA/RCode/2017_1_30_CopesRule_Source_Func.R', chdir = TRUE) #call cource file for functions

#sources for Jon Marcot's code and specimen measurements 
source("https://dl.dropbox.com/s/8jy9de5owxj72p7/strat.R")
source("https://dl.dropbox.com/s/253p4avcvb66795/occFns.R")
source("https://dl.dropbox.com/s/9gdafsqss2b586x/phy_dateTree.R")
source("https://dl.dropbox.com/s/9tdawj35qf502jj/amandaSrc.R")
source("https://dl.dropbox.com/s/rlof7juwr2q4y77/blasto_Birlenbach.R")

tree_base <- read.nexus("~/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_2_1_UngulataBackboneTree")
tree_base
class(tree_base)

clade_sp <- read.csv("~/Dropbox/ungulate_RA/2017_1_26_Clade_species.csv", stringsAsFactors = FALSE)
clade_sp

#clade_sp <- read.csv("C:/Users/Evan/Dropbox/ungulate_RA/2017_1_26_Clade_species.csv",stringsAsFactors = FALSE)
#clade_sp

MCRA_Codes <- read.csv("~/Dropbox/ungulate_RA/2017_1_26_MCRA_Codes.csv", stringsAsFactors = FALSE)
MCRA_Codes

#MCRA_Codes <- read.csv("C:/Users/Evan/Dropbox/ungulate_RA/2017_1_26_MCRA_Codes.csv", stringsAsFactors = FALSE)
#MCRA_Codes <- read.csv("C:/Users/Evan/Dropbox/ungulate_RA/2017_1_26_MCRA_Codes_Laptop.csv", stringsAsFactors = FALSE) #for calling tree files on Evan's laptop
#MCRA_Codes

#place function to request number of iterations

#iter <- as.numeric(iterations()) #select number of iterations to run

reps <- 1

#tree.list <- list()
tree.list <- replicate(n = reps,expr = tree_resolution(backbone_tree = tree_base,Species_file = clade_sp,MRCA_file = MCRA_Codes))
tree.list
class(tree.list)

tree.list <- lapply(X=tree.list, FUN = drop_DupTips)


#####################################################################################################################################################
#cal3 requires tree=class phylo so line below is to test
#test.tree <- tree_resolution(backbone_tree = tree_base,Species_file = clade_sp,MRCA_file = MCRA_Codes)

#will call cal3 next
#cal3 is wrking but generating polytomies
#tree.cal<- cal3TimePaleoPhy(tree = test.tree, timeData = ranges, brRate = q_extR, extRate = q_extR, sampRate = r_sam, ntrees = 1, anc.wt = 1, node.mins = NULL, dateTreatment = "firstLast", FAD.only = FALSE, adj.obs.wt = TRUE, root.max = 200, step.size = 0.1, randres = FALSE, noisyDrop = TRUE, tolerance = 1e-04, diagnosticMode = FALSE, plot = FALSE)

#tree_cal <- bin_cal3TimePaleoPhy(tree = test.tree, timeList = ranges, brRate = q_extR, extRate = q_extR, sampRate = r_sam, ntrees = 1, anc.wt = 1, node.mins = NULL, dateTreatment = "firstLast", FAD.only = FALSE, adj.obs.wt = TRUE, root.max = 200, step.size = 0.1, randres = FALSE, noisyDrop = TRUE, tolerance = 1e-04, diagnosticMode = FALSE, plot = FALSE)

class(test.tree)


#Date tree
occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Artiodactyla,Perissodactyla&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)

tree_dated <- list()
tree_dated2 <- list()
#tree_dated <- lapply (X = tree.list, FUN = date_tree, occs) #need to swap LO and FO column headers; Check to see if the objects are being save correctly (i.e. tree dated properly)

###############
#ONE OPTION OF DO THE DATING IS HAVE CAL3 BE WITHIN an underlying function in source file
tree_dated <- lapply(X = tree.list, FUN = date_treeCal3, occs).

plot(tree_dated, font = 3, cex =0.2, label.offset = 1, adj = 0)


##############
#another option to run dating is to have cal3 call be in main program with one or more function being used to derive its inputs
rate_func <-getSamExtRate(occs)

optim(parInit(rate_func), rate_func, lower = parLower(rate_func), upper = parUpper(rate_func), method = "L-BFGS-B")
	
	#parnames(ranges_binned)
	#parbounds(ranges_binned)
	#parLower(ranges_binned)
	#parUpper(ranges_binned)
	#parInit(ranges_binned)
	
	class(rate_func)
	
	#outputs of make_durationFreqCont
red <-vector()
red <- parInit(rate_func)
class(red)
red

#q=instantaneous per-capita extinction rate; 
q_extR <- vector()
q_extR <- red[1]
q_extR

#r=instantaneous per-capita sampling rate;
r_sam<- vector()
r_sam <- red[2]
r_sam


#this method will require to use the code to extract extinction and sampling rates in main program

tree_dated2 <- lapply (X = tree.list, FUN = cal3TimePaleoPhy, ranges, q_ext, q_ext, r_sam) class(tree_dated)

plot(tree_dated, font = 3, cex =0.2, label.offset = 1, adj = 0)

#######################################################################################################################################################
