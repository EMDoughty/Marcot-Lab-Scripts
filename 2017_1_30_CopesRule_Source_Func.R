iterations <- function() 
{
  #set number of iterations for program to run
  iter <- 1
  iter<-readline(prompt="Enter number of iterations: ")
  as.integer(iter)
  return (iter)
}

tree_resolution <- function(backbone_tree,Species_file,MRCA_file)
{
#Resolve backbone Tree
tree_resolv <- multi2di(backbone_tree) #backbone_tree)
#plot(tree_resolv, font = 1, edge.width = 0.5, cex =0.2, label.offset = 1, adj = 0)
names(tree_resolv)
tree_resolv

#Identify landing edges
i <- 1 # track number of wildcard taxa completed
k <- 1 # track number of iterations completed
node <- vector()
#node <- c(0,0,0,0)#need to set so sample is not being selected as 0 so wil need to be defined in while loops near numberPos
node_choice <- 0

while(i < nrow(MRCA_file) + 1) {
	print('k')
	print(k)
	print('i')
	print(i)
	clade_tree <- read.nexus(file = MRCA_file$Filename[i])
  clade_tree <- multi2di(clade_tree)
#	plot(clade_tree, font = 0.1, cex =0.2, label.offset = 1, adj = 0)
	numberPos <- MRCA_file[i,"X._positions"]
	numberPos 	
	j <- 5 # track number of possible positions completed
	m <- 1 # save info to node list
	while(j < (5 + numberPos)) {
		print('j')
		print(j)
		position <- MRCA_file[i,j]
		print(position)
		clade_subset<- Species_file[Species_file["Family.Clade"] == position, 2:3] #c("Species.1", "Species.2")]
		clade_subset
		speciesA <- clade_subset[,"Species.1"]
		speciesA
		print(speciesA)
		speciesB <-	clade_subset[,"Species.2"]
		print(speciesB)
		print(m)
		node[m] <- fastMRCA(tree=tree_resolv, sp1=speciesA, sp2=speciesB)
		print("node")
		print(node)
		m <- m+1
		j <- j+1}

#Randomly place taxa on landing branches
	node_choice <- sample(node,size = 1, replace = FALSE, prob =NULL)
	#print("node_choice")
	#print(node_choice)
		
#Attach Subtree
	tree_resolv <- bind.tree(x=tree_resolv, y=clade_tree, position = node_choice)
		
	#plot(tree_bound, font = 3, cex =0.2, label.offset = 1, adj = 0)
	i <- i + 1
	}

#Troubleshooting: Paint Branches if issues occur
#this.map <- paintSubTree(tree=tree_resolv, node=fastMRCA(tree=tree_resolv, sp1="Hypertragulus_minor", sp2="Yumaceras_ruminalis"), state="2")
#plotSimmap(this.map,fsize=0.1, ftype = "i")

plot(tree_resolv, font = 1, cex =0.2, label.offset = 1, adj = 0)

tree_resolv <- multi2di(tree_resolv) 

#save the resolved tree to file or lists

#plot(tree_resolv, font = 0.1, cex =0.2, label.offset = 1, adj = 0)

k <- k+1

return(tree_resolv)
}

date_tree <-function(this.tree, occs)
{
	intervals <- makeIntervals(startDate=57, endDate=1, intervalLength=2)

	#occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Artiodactyla,Perissodactyla&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
	
	occs <- appendTaxonNames1.2(occs, taxonomic.level="species", keep.indet=TRUE)

	ranges <- getTaxonRangesFromOccs(occs, random = TRUE)	
	
	#rownames(ranges) <-gsub(pattern = "[[space:]]",replacement ="_", x = rownames(ranges))
	rownames(ranges) <- str_replace_all(rownames(ranges),fixed(" "), "_")

	return(this.tree)
}


getSampleRate <- function(occs)
{
	
	intervals <- makeIntervals(startDate=57, endDate=1, intervalLength=2)
	intervals

	occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Artiodactyla,Perissodactyla&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
	
	occs <- appendTaxonNames1.2(occs, taxonomic.level="species", keep.indet=TRUE)
	head(occs)
	
	ranges <- getTaxonRangesFromOccs(occs, random = TRUE)	
	
	#rownames(ranges) <-gsub(pattern = "[[space:]]",replacement ="_", x = rownames(ranges))
	rownames(ranges) <- str_replace_all(rownames(ranges),fixed(" "), "_")
	
	ranges
	class(ranges)
	
	#how do we want to have these sampled? (by epoch, stage, NALMA. etx)
	#
	
	mat_ranges <- list(ranges,intervals) # interval files needs to be comprised of persistent taxa
	mat_ranges
	
	#code from Bapst Site will likely not work as fuctions--freqRat and getSampDisc--are for discrete data
	freqRat(timeData = mat_ranges,plot=TRUE)
	freqRat()
	
	SPres <- getSampProbDisc(ranges)
	SPres
	sRate <-sProb2sRate(SPres[[0]][], int.legnth=meanInt)

	#sampleRanges(taxad = ranges, r, alpha = 1, beta = 1, rTimeRatio = 1, modern.samp.prob = 1, min.taxa = 2, ranges.only = TRUE, minInt = 0.01, merge.cryptic = TRUE, randLiveHat = TRUE, alt.method = FALSE, plot = FALSE)  ####This appears to be meant for simulations

	#ranges_binned <- binTimeData(timeData= ranges, int.length = 5, start = NA, int.times = NULL) # funtion list as Being specifically only for simulations and not real data
	
	
	#sets ranges_binned as a function that runs everytime it is called which leads to new values being output everytime
	ranges_binned <- make_durationFreqCont(timeData = ranges, groups = NULL, drop.extant = TRUE, threshold = 0.01, tol = 1e-04)
		
optim(parInit(ranges_binned), ranges_binned, lower = parLower(ranges_binned), upper = parUpper(ranges_binned), method = "L-BFGS-B")
	
	#parnames(ranges_binned)
	#parbounds(ranges_binned)
	#parLower(ranges_binned)
	#parUpper(ranges_binned)
	#parInit(ranges_binned)
	
#outputs of make_durationFreqCont
red <-vector()
red <- parInit(ranges_binned)
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

#R is per-interval taxonomic sampling probability
sRate <- sRate2sProb(r_sam, int.length = 1)
sRate

#will call cal3 next
#cal3 is wrking but generating polytomies
tree.cal<- cal3TimePaleoPhy(tree = test.tree, timeData = ranges, brRate = q_extR, extRate = q_extR, sampRate = r_sam, ntrees = 1, anc.wt = 1, node.mins = NULL, dateTreatment = "firstLast", FAD.only = FALSE, adj.obs.wt = TRUE, root.max = 200, step.size = 0.1, randres = FALSE, noisyDrop = TRUE, tolerance = 1e-04, diagnosticMode = FALSE, plot = FALSE)

plot(tree.cal, font = 3, cex =0.2, label.offset = 1, adj = 0)

#tree_cal <- bin_cal3TimePaleoPhy(tree = test.tree, timeList = ranges, brRate = q_extR, extRate = q_extR, sampRate = r_sam, ntrees = 1, anc.wt = 1, node.mins = NULL, dateTreatment = "firstLast", FAD.only = FALSE, adj.obs.wt = TRUE, root.max = 200, step.size = 0.1, randres = FALSE, noisyDrop = TRUE, tolerance = 1e-04, diagnosticMode = FALSE, plot = FALSE)

	return(rate)
}






############################################################################################################################################

#specimenMat <- getSpecimenMatFromMeasurements(filename="https://dl.dropbox.com/s/423x0zn3mpxuwc7/specimens.csv")
#specimenMat <- merge(specimenMat, getBirlenbachBlastoSpecimens(filename="https://dl.dropbox.com/s/943dq4lb1kd1h1e/blastoBirlenbach20140207.csv", info.file="https://dl.dropbox.com/s/nhrqzrclwtr5c8e/info2.csv"), all=TRUE)
#specimenMat <- merge(specimenMat, getLiteratureSpecimenMat(filename="https://dl.dropbox.com/s/qef8ts9a73j5ukb/literature.csv"), all=TRUE)

#specimenMat[sapply(specimenMat, is.nan)] <- NA
#specimenMat$species <- getCurrentTaxa(specimenMat$species)

#upLabels<-c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W") #"P2_L","P2_W",
#loLabels <- casefold(upLabels)
#thisMat <- aggregate(specimenMat[,c(upLabels, loLabels)], by=list(species=specimenMat$species), mean, na.rm=TRUE)
#thisMat <- aggregate(specimenMat[,c(upLabels, loLabels)], by=list(species=specimenMat$species), median, na.rm=TRUE)
#thisMat[,sapply(thisMat, is.numeric)] <- thisMat[,sapply(thisMat, is.numeric)] / 10  #converts mm measurements to cm for compatibility with Janis regressions
#thisMat <- transform(thisMat, p4_a=p4_l*p4_w, m1_a=m1_l*m1_w, m2_a=m2_l*m2_w, m3_a=m3_l*m3_w, M2_A=M2_L*M2_W)
#thisMat[sapply(thisMat, is.nan)] <- NA

#thisMat <- appendMissingPaleoDBSpecies(thisMat, ranges)		# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files
#thisMat[,"bodyMass"] <- getBodyMassVectorFromThisMatAllMeasures(thisMat, linked.files=TRUE)
#thisMat$bodyMass <- fillMissingBodyMasses(thisMat)				# this fills taxa missing their body mass with the average body mass of its cogeners
#thisMat[!sapply(thisMat, is.finite)] <- NA

#rownames(thisMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisMat))

############################################################################################################################################

#rownames(ranges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(ranges))

# tree <- read.nexus("~/Dropbox/code/R/NAPC2014/dat/analysis.tre")
#tree <- read.nexus("https://dl.dropbox.com/s/0f82kvl5kagef3i/analysis.tre")

#tree <- tree
#tree$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree$tip.label, fixed=TRUE)))

## randomly drop any tips that are now duplicated
	
#this.tree <- tree
#dup.taxa <- unique(this.tree$tip.label[duplicated(this.tree$tip.label)])
#for (i in seq_along(dup.taxa)) this.tree <- drop.tip(this.tree, tip=sample(which(this.tree$tip.label==dup.taxa[i]), size=sum(this.tree$tip.label==dup.taxa[i])-1))
#this.tree <- ladderize(multi2di(this.tree), right=FALSE)
#this.tree <- dateTreeWithRanges(tree=this.tree, dates=ranges, within.error=TRUE, extra=1, flex.internals=TRUE, rootAdj=5, keep.missing=FALSE) # internals will be flexed by extendNodeAges below
#this.tree <- rectifyNodeAges(this.tree, extra=1)
#this.tree <- extendNodeAges(this.tree, method="uniform", rootAdj=10, extra=1)

