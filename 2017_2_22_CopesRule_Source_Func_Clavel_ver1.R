tree_resolution_v3 <- function(tree.backbone, clade.definitions, wildcard.positions, tree.list.wildcards, occs, thisMat) {
	tree_resolv <- bindWildCards(tree.backbone, clade.definitions, wildcard.positions, tree.list.wildcards)
	#drop duplicate tips from tree
	tree_resolv <- drop_DupTips(tree_resolv)
	print("Dropping duplicate tips is complete")

	#drop tips from tree that do no coincide with entry in species data
	tree_resolv <- drop.tip(tree_resolv, tip=tree_resolv$tip.label[!tree_resolv$tip.label %in% thisMat$species])
	tree_resolv <- date_treeJon(this.tree = tree_resolv, occs = occs, root.max = 85)

	print("Tree Dating is complete")
	#paintWildcardClades_Dated(tree_dated)

	return(tree_resolv)
}

getLandingNode <- function(this.position, this.tree, clade.definitions) {
	n.members <- clade.definitions[clade.definitions["Family.Clade"] == this.position, 2]
	clade.members <- unlist(clade.definitions[clade.definitions["Family.Clade"] == this.position, 3:(2+n.members)]) #c("Species.1", "Species.2")])
	node.landing <- findMRCA(tree = this.tree, tips = as.vector(clade.members))
	node.landing
}

bindWildCards <- function(tree.backbone, clade.definitions, wildcard.positions, tree.list.wildcards) {
	tree_resolv <- multi2di(tree.backbone)
	tree_resolv$root.length <- 1
	tree_resolv$edge.length <- array(data=1, dim=nrow(tree_resolv$edge))
#	plot(tree_resolv, direction="up", cex =0.5, no.margin=TRUE)
	
	for (this.wildcard in seq_len(nrow(wildcard.positions))) {		
		#read in wildcard trees indexed in the MRCA.csv file
		# tree.wildcard <- read.nexus(file = wildcard.positions$Filename[this.wildcard]) 
	  	tree.wildcard <- multi2di(tree.list.wildcards[[this.wildcard]])
	  	tree.wildcard$edge.length <- array(data=1, dim=nrow(tree.wildcard$edge))
	  	tree.wildcard$root.edge <- 1
	  	# plot(tree.wildcard, direction="up", cex =0.5)
	  	
		this.position <- wildcard.positions[this.wildcard, 4+sample(x=wildcard.positions[this.wildcard, "X._positions"], size=1)]
		
		if (this.position == "Ruminantiamorpha") { # find the nodes between the suina and the base of the Ruminantia
			node.top <- getLandingNode("Ruminantiamorpha", this.tree=tree_resolv, clade.definitions)
			node.bottom <- getLandingNode("Ruminantiamorpha_Tayassuidae", this.tree=tree_resolv, clade.definitions)
			node.series <- getNodes(top = node.top, bottom = node.bottom,tree= tree_resolv)
			if(length(node.series) == 1) {node.landing <- node.series 
			} else node.landing <- sample(x=getNodes(top = node.top, bottom = node.bottom,tree= tree_resolv), size=1)	#### this is likely the reason for why taxa getting placed at random: sample is selecting from within the number rather tahn across the node i.e. node 1240 sample selectes 87		
		} else if (this.position == "Tylopoda") { # find the nodes between the basal node of Tylopoda and the base of the Camelidae/Oromerycidae clade
			node.top <- getLandingNode("Camelidae_Oromerycidae", this.tree=tree_resolv, clade.definitions)
			node.bottom <- getLandingNode("Tylopoda", this.tree=tree_resolv, clade.definitions)
			node.landing <- sample(x=getNodes(top = node.top, bottom = node.bottom,tree= tree_resolv), size=1)
			node.series <- getNodes(top = node.top, bottom = node.bottom,tree= tree_resolv)
			if(length(node.series) == 1) {node.landing <- node.series 
			} else node.landing <- sample(x=getNodes(top = node.top, bottom = node.bottom,tree= tree_resolv), size=1)	#### this is likely the reason for why taxa getting placed at random: sample is selecting from within the number rather tahn across the node i.e. node 1240 sample selectes 87		
		} else node.landing <- getLandingNode(this.position, this.tree=tree_resolv, clade.definitions)
			
		#Attach Subtree
		#cat("Tree_resolv", tree_resolv$root.length, "Wildcard", tree.wildcard$root.length)
		tree_resolv <- bind.tree(tree_resolv, tree.wildcard, node.landing, position=0.5)
		tree_resolv$edge.length <- array(data=1, dim=nrow(tree_resolv$edge))
	}
	# internal.nodes <- sort(unique(tree_resolv$edge[,1]))
	# sapply(internal.nodes, function(x) sum(tree_resolv$edge[,1]==x))

	tree_resolv <- reorder(tree_resolv, order = "pruningwise") # this is done to prevent graphical error when tree is plotted
	tree_resolv <- ladderize(tree_resolv, right=FALSE)
	# plot(tree_resolv, type="phylogram", use.edge.length = TRUE, direction="up", no.margin=TRUE, cex=0.5)
	tree_resolv

	# print("Binding of Wildcards is complete")
	#paintWildcardClades(tree_resolv)
}

############################################################################################################################################

#function meant as a check to make sure that the proper wildcards (i.e. oreodonts) are being placed basla to Ruminanitia but derived to suina
checkRuminantiamorpha <- function(tree, wildcardTree, clade.definitions)
{
	require(caper)
	
	loop_break <- FALSE
	
	#check to see if multi2di positioned Hypertragulidae and other ruminants more derived to ruminantiomorphs (wildcards)
	this.tree <- tree
		
	#find node of ruminants
	clade1 <- "Ruminantiamorpha"

	n <- clade.definitions[clade.definitions["Family.Clade"] == clade2, 2]
	n
	clade.members<- unlist(clade.definitions[clade.definitions["Family.Clade"] == clade1, 3:(2+n)])
	clade.members
	node_rumin <-findMRCA(tree = this.tree, tips = as.vector(clade.members), type = c("node"))
	ruminDesc <- getDescendants(this.tree, node = node_rumin)
	
	ruminTips <- clade.members(node_rumin, this.tree, tip.labels = TRUE)
	names(this.tree$edge)
	
	#find node of wildcards
	#wild clade with ruminants be best?
	wild_tips <- c(wildcardTree$tip.label, ruminTips)
	wildcard <- findMRCA(tree, tips = wild_tips)
	
	#check to see if ruminant node decendents includes node of wildcard clade
	
	wildDesc <- getDescendants(this.tree, node = wildcard)
	
	if(!wildcard %in% ruminDesc & !wildcard == node_rumin)
	{
		loop_break <- TRUE
	}
	
	#return TRUE or FALSE on whether to break the loop
	
	return(loop_break)
}

############################################################################################################################################

#function adds arbitrary species Rtip1 to base of wildcard tree: used in prior troubleshooting measures but not surrently implimented
addArbTip2Root <- function(tree.wildcard)
{
#append arbitrary tip to root of clade to prevent basal species/branches from being separated from clade when multi2di is ran on bound tree; tip will be removed after the binding 
  	tips <- tree.wildcard$tip.label
  	tips
	node_wild <- findMRCA(tree.wildcard, tips = tips)
	node_wild
	tree.wildcard <- bind.tip(tree.wildcard, tip.label="Rtip1", where = node_wild)
	tree.wildcard <- root(tree.wildcard, outgroup = "Rtip1", resolve.root = TRUE)
	tree.wildcard$edge.length <- array(data=1, dim=nrow(tree.wildcard$edge))
	plot(tree.wildcard, no.margin=TRUE)
	nodelabels()
	
	return(tree.wildcard)
}

############################################################################################################################################


iterations <- function() 
{ 
  #set number of iterations for program to run
  iter <- 1
  iter<-readline(prompt="Enter number of iterations: ")
  as.integer(iter)
  return (iter)
}

############################################################################################################################################

#Recursive function by Jon to find all nodes between two designated nodes: used to find multiple sites where tylopods and ruminantiamorphs could be attached despite ambiguity
getNodes <- function (top, bottom, tree) {
    getNodesMachine(this.node=top, bottom=bottom, tree=tree, node.vec=top)
}
getNodesMachine <- function(this.node, bottom, tree, node.vec) {
    anc.node <- tree$edge[tree$edge[,2]==this.node,1]
    if (anc.node != bottom) {
        if (anc.node == length(tree$tip.label) +1) {
            warning("\7*** in getNodesMachine: Root reached before bottom node ***")
            return(NULL)
        }
        node.vec <- c(node.vec, anc.node)
        node.vec <- getNodesMachine(this.node=anc.node, bottom=bottom, tree=tree, node.vec=node.vec)
    }
    node.vec
}

############################################################################################################################################

date_treeJon<-function(this.tree, occs, root.max)
{
	ranges <- getTaxonRangesFromOccs(occs, random =TRUE)
	rownames(ranges) <- gsub(pattern = "[[:space:]]",replacement ="_", x = rownames(ranges))
	this.tree <- dateTreeWithRanges(this.tree = this.tree, strat.ranges=ranges, within.error=FALSE, min.bl=1.0, smooth.node.ages=TRUE, root.max=root.max) # internals will be flexed by extendNodeAges below
	# tree_dated <- multi2di(this_dated)	
	# tree_dated$root.date <- root.max
	#check tree branch legnths
	
#	ttree <- tree_dated
#	par(mfrow=c(2,1))
#	hist(ttree$edge.length, breaks=seq(from=0, to=ceiling(max(ttree$edge.length)), by=0.5))
#	plot(ttree, no.margin=TRUE, cex=0.5, edge.width = 0.5)
#	max(ttree$node.date)
#	min(ttree$edge.length)
	
	return(this.tree)
}

############################################################################################################################################

drop_DupTips <- function(tree)
{
	# tree$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree$tip.label, fixed=TRUE)))
	dup.taxa <- unique(tree$tip.label[duplicated(tree$tip.label)])
	for (i in seq_along(dup.taxa)) tree <- drop.tip(tree, tip=sample(which(tree$tip.label==dup.taxa[i]), size=sum(tree$tip.label==dup.taxa[i])-1))
	tree	
		
	
	return(tree)
}

############################################################################################################################################

compile_SpecData_Old <- function(occs)
{
	occs <- appendTaxonNames1.2(occs, taxonomic.level="species", keep.indet=TRUE)
	head(occs)
	
	ranges <- getTaxonRangesFromOccs(occs, random = TRUE)	
	
	
	#compile and lable dental measurments for specimens
	specimenMat <- getSpecimenMatFromMeasurements(filename="https://dl.dropbox.com/s/423x0zn3mpxuwc7/specimens.csv")
	specimenMat <- merge(specimenMat, getBirlenbachBlastoSpecimens(filename="https://dl.dropbox.com/s/943dq4lb1kd1h1e/blastoBirlenbach20140207.csv", info.file="https://dl.dropbox.com/s/nhrqzrclwtr5c8e/info2.csv"), all=TRUE)
	specimenMat <- merge(specimenMat, getLiteratureSpecimenMat(filename="https://dl.dropbox.com/s/qef8ts9a73j5ukb/literature.csv"), all=TRUE)

	specimenMat[sapply(specimenMat, is.nan)] <- NA
	specimenMat$species <- getCurrentTaxa(specimenMat$species)

	upLabels<-c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W") #"P2_L","P2_W",
	loLabels <- casefold(upLabels)
	thisMat <- aggregate(specimenMat[,c(upLabels, loLabels)], by=list(species=specimenMat$species), mean, na.rm=TRUE)
	thisMat <- aggregate(specimenMat[,c(upLabels, loLabels)], by=list(species=specimenMat$species), median, na.rm=TRUE)
	thisMat[,sapply(thisMat, is.numeric)] <- thisMat[,sapply(thisMat, is.numeric)] / 10  #converts mm measurements to cm for compatibility with Janis regressions
	thisMat <- transform(thisMat, p4_a=p4_l*p4_w, m1_a=m1_l*m1_w, m2_a=m2_l*m2_w, m3_a=m3_l*m3_w, M2_A=M2_L*M2_W)
	thisMat[sapply(thisMat, is.nan)] <- NA

	thisMat <- appendMissingPaleoDBSpecies(thisMat, ranges)		# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files
	
	#Approximate body mass
	thisMat[,"bodyMass"] <- getBodyMassVectorFromThisMatAllMeasures(thisMat, linked.files=TRUE)
	thisMat$bodyMass <- fillMissingBodyMasses(thisMat)				# this fills taxa missing their body mass with the average body mass of its cogeners
	thisMat[!sapply(thisMat, is.finite)] <- NA

	rownames(thisMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisMat))
	
	return (thisMat)
}

############################################################################################################################################

#function is meant to be ran prior to analysis to bring all measurement datasets together into a single entity
getSingleSpeciesMatrix <- function() {
		
	#compile and lable dental measurments for specimens
	specimenMat <- getSpecimenMatFromMeasurements(filename="https://dl.dropbox.com/s/423x0zn3mpxuwc7/specimens.csv")
	specimenMat <- merge(specimenMat, getBirlenbachBlastoSpecimens(filename="https://dl.dropbox.com/s/943dq4lb1kd1h1e/blastoBirlenbach20140207.csv", info.file="https://dl.dropbox.com/s/nhrqzrclwtr5c8e/info2.csv"), all=TRUE)
	specimenMat <- merge(specimenMat, getLiteratureSpecimenMat(filename="https://dl.dropbox.com/s/qef8ts9a73j5ukb/literature.csv"), all=TRUE)

	specimenMat[sapply(specimenMat, is.nan)] <- NA
	specimenMat$species <- getCurrentTaxa(tax.vec = specimenMat$species)

	upLabels<-c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W") #"P2_L","P2_W",
	loLabels <- casefold(upLabels)

	### node that specimens are aggregated by their medians, so as to minimize the effect of outlier measurements
	thisMat <- aggregate(specimenMat[,c(upLabels, loLabels)], by=list(species=specimenMat$species), median, na.rm=TRUE)
	
	thisMat[,sapply(thisMat, is.numeric)] <- thisMat[,sapply(thisMat, is.numeric)] / 10  #converts mm measurements to cm for compatibility with Janis regressions
	thisMat <- transform(thisMat, p4_a=p4_l*p4_w, m1_a=m1_l*m1_w, m2_a=m2_l*m2_w, m3_a=m3_l*m3_w, M2_A=M2_L*M2_W)
	thisMat[sapply(thisMat, is.nan)] <- NA
	
	#thisMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = thisMat$species)
	rownames(thisMat) <- thisMat$species

	return (thisMat)
}

############################################################################################################################################

#Compile and format matrix of all measurements from multiple sources
appendRegressionCategories <- function(thisMat, regMat) {
    
    # uniqTax <- unique(occs[,c("accepted_name", "order", "family", "genus")])
    # uniqTax <- getTaxonomyForTaxa(rownames(thisMat)) # neeed to remove Squalodon and Paraceratheriidae
    # nrow(uniqTax)
    # print("got uniqTax")
   uniqTax <- lapply(c("Artiodactyla", "Perissodactyla"), FUN=getTaxonomyForOneBaseTaxon)
   uniqTax <- rbind(uniqTax[[1]], uniqTax[[2]])
   uniqTax$taxon_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$taxon_name)
    #Remove records that lack both a family and genus designation (and occurence data) from the Paleobio database
   # nrow(thisMat)
   # head(uniqTax)
   # rownames(uniqTax) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(uniqTax))
   
   # measure.matSave <- thisMat
   # #measure.mat <- measure.matSave
   # measure.matTest <- thisMat
   # nrow(measure.matTest)
    
    #delete unique taxa that lack family and genus
   thisMat$family <- uniqTax$family[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
   # thisMat$family <- sapply(thisMat$family, as.character)
   thisMat$family[thisMat$family == ""] <- NA
   
   thisMat$genus <- uniqTax$genus[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
   # thisMat$genus <- sapply(thisMat$genus, as.character)
   thisMat$genus[thisMat$genus == ""] <- NA
   
  # measure.matTest <- cbind(measure.matTest, family.names)
   # # print("cbind1")
  # measure.matTest <- cbind(measure.matTest, genus.names) 
   # # print("cbind2")
       
  # measure.matTest <- measure.matTest[!is.na(measure.matTest$family.names) | !is.na(measure.matTest$genus.names),]
   # # class(measure.matTest)
   
   # # nrow(measure.matTest) #798 species remain using occs; 806 are returned using the getTaxononomyForTaxa function
         
   # #drop additional columns for damily.names and genus.names
   # measure.matTest <- subset(measure.matTest,select= -c(family.names, genus.names))
         
    # thisMat <- measure.matTest
     
    #### 
    #append regression catagories to each species
    ####
         
    family.names <- uniqTax$family[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
    reg.vec <- regMat$reg[!is.na(regMat$family)][match(family.names, regMat$family[!is.na(regMat$family)])]  # this is the regression "labels" for the species from measure.mat in the correct order, based on family name
	   
    genus.names <- uniqTax$genus[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
    reg.vec[is.na(reg.vec)] <- regMat$reg[match(genus.names, regMat$genus)][is.na(reg.vec)]   #this is the regression "labels" for the species from measure.mat in the correct order, based on genus name; it appears that having regMat$genus[!is.na(regMat$genus)] will cause the index to improperly assign regressions
       
    thisMat$reg.vec <- reg.vec
    # print("cbindlast")
    
    #check for taxa that are not recieving a regression
 	#missingReg <- measure.mat[is.na(measure.mat $reg.vec),]
 	#missingReg <- missingReg[!grepl("sp.",missingReg$species),]
 	#missingReg <- missingReg[!grepl("indet.",missingReg$species),]
 	#missingReg <- missingReg[!grepl("cf.",missingReg$species),]
 	#write.csv(missingReg, '/Users/evandoughty/Dropbox/ungulate_RA/RCode/JonCode/2017_2_27_missingReg.csv')
 	
     # nrow(thisMat) #803 species remain after final removal 
     
     # rownames(thisMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisMat))
     
  return(thisMat)
}

############################################################################################################################################

approxBodyMass <- function(thisMat)
{
	
	#thisMat <- measureMat
	#Approximate body mass
	thisMat[,"bodyMass"] <- getBodyMassVectorFromThisMatAllMeasures(thisMat, linked.files=TRUE)
	thisMat$bodyMass <- fillMissingBodyMasses(thisMat)	# this fills taxa missing their body mass with the average body mass of its cogeners
	thisMat[!sapply(thisMat, is.finite)] <- NA

	#rownames(thisMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisMat))
	
	thisMat$species <- rownames(thisMat)
		
	return(thisMat)
}


############################################################################################################################################

#Function for checking for, isolating, and appending species entries that were not present in the initial regression catagorizing matrix.
#Entries are mearly located and appended. The reg.vec column in the regMat file will retain all NA values for the taxa missing those entries
# unless the user changes or removes them via manual or automatic means.

checkMissingReg<- function(measReg, uniqTax)
{
	missingReg <- measure.matReg[is.na(measure.matReg$reg.vec),]

	#merge/match with occs fiel to get order, genus, and family columns

	occurrences.order_name <- uniqTax$order[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.order_name)

	occurrences.family_name <- uniqTax$family[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.family_name)

	occurrences.genus_name <- uniqTax$genus[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.genus_name)

	missingReg <- missingReg[,c("occurrences.order_name","occurrences.family_name","occurrences.genus_name","species","reg.vec")]
	head(missingReg)
	colnames(missingReg)[colnames(missingReg) == "species"] <- "taxon"
	colnames(missingReg)[colnames(missingReg) == "reg.vec"] <- "reg"

	regMat <- rbind(regMat, missingReg)

	return(regMat) 
}

###############

#function for testing dating using variable root ages and other settings
dateTreeWithRanges_outMin <- function(this.tree, strat.ranges, within.error=FALSE, random.ages=TRUE, min.date=0.0, min.bl=0.0, smooth.node.ages=FALSE, root.max=max(strat.ranges), keep.missing=FALSE, missing.date=NULL, include.ranges=FALSE) {
	if (keep.missing) {
		strat.ranges <- strat.ranges[match(this.tree$tip.label, rownames(strat.ranges)),]		#reorder taxa in strat.ranges by the order in the this.tree
		rownames(strat.ranges) <- this.tree$tip.label	# the match above erases unmatched species from the rownames - this restablishes them
		if (!is.null(missing.date)) strat.ranges[!is.finite(strat.ranges)] <- missing.date
	} else if (any(!(this.tree$tip.label %in% rownames(strat.ranges)))) {
			this.tree <- drop.tip(this.tree, this.tree$tip.label[!(this.tree$tip.label %in% rownames(strat.ranges))]) # cut taxa that are not in ranges
			# strat.ranges <- strat.ranges[match(this.tree$tip.label, rownames(strat.ranges)),]		#reorder taxa in strat.ranges by the order in the this.tree
	}

	# strat.ranges <- strat.ranges[rownames(strat.ranges) %in% this.tree$tip.label,] 		#cull taxa in strat.ranges list that aren't in the this.tree
	# strat.ranges <- strat.ranges[match(this.tree$tip.label, rownames(strat.ranges)),]		#reorder taxa in strat.ranges by the order in the this.tree
	strat.ranges <- strat.ranges[this.tree$tip.label,]
	
	# set node dates of tips
	if (within.error) {
		if (random.ages) { this.tree$node.date <- getRandomRanges(strat.ranges)[,"FO"] 
		} else this.tree$node.date <- rowMeans(strat.ranges[,1:2]) #if without error, use the midpoint of FA dates
	} else this.tree$node.date <- strat.ranges[,"FO"]	# if error on the FO and LO is not provided, a column of FOs should be
		
	# make space for internal nodes
	this.tree$node.date <- c(this.tree$node.date, array(data=NA, dim=this.tree$Nnode)) #add NA dates for internal nodes
	this.tree <- setNodeDatesMachine(this.tree=this.tree, nodeID=length(this.tree$tip.label)+1, min.bl=0)
	this.tree <- makeBLFromNodeDates(this.tree)

	# apply min.bl to terminal branches ONLY
 	this.tree <- applyMinBL(this.tree, terminals.only=TRUE, min.bl=min.bl) # apply min.bl to terminal branches ONLY

 	# this determines whether to reduce min.bl to stay under max.root
 	if (any(this.tree$edge.length <= min.bl)) {
		this.min <- (root.max - max(this.tree$node.date))/getMaxShortBranchesRootToTip(this.tree, min.bl=min.bl)
		if (this.min < min.bl) {
			message <- paste0("in dateTreeWithRanges: To keep the root younger than root.max (",root.max,"), minimum length of internal branches set to ",round(this.min,3)," (shorter length than min.bl [",round(min.bl,3),"])")
			warning(message)
		}
		this.min <- min(this.min, min.bl)
		this.tree <- applyMinBL(this.tree, min.bl=this.min)
	}
	
	if (any(!is.finite(this.tree$node.date))) this.tree <- backfillMissingNodeDates(this.tree, min.bl)

	if (smooth.node.ages) this.tree <- smoothNodeAges(this.tree, method="uniform", root.max=root.max, min.bl=min.bl)

	this.tree <- makeBLFromNodeDates(this.tree) # called within smoothNodeAges, but has to be here in case smooth.node.ages is FALSE
	if (include.ranges) {
		if (within.error) { this.tree$strat.range <- cbind(this.tree$node.date[seq_along(this.tree$tip.label)], apply(strat.ranges[,c("ageLO_max", "ageLO_min")], 1, function(x) runif(1, min=x[2], max=x[1]) ))
		} else this.tree$strat.range <- cbind(this.tree$node.date[seq_along(this.tree$tip.label)], strat.ranges[,"LO"])
		this.tree$strat.range[this.tree$strat.range[,1]<this.tree$strat.range[,2],] <- this.tree$strat.range[this.tree$strat.range[,1]<this.tree$strat.range[,2],c(1,1)]
	}
	this.min
}

############################################################################################################################################

#functions compare the tree and measurement matrix with one another to verify that the same species are present: species not found in the other are culled

comp_TreeSpecMat_outTree <- function(this.tree,thisMat)
{
	thisMat <- thisMat[is.finite(thisMat$bodyMass),] 
	this.tree <- drop.tip(phy=this.tree, tip=this.tree$tip.label[!this.tree$tip.label %in% rownames(thisMat)])
	
	return(this.tree)
}

comp_TreeSpecMat_outMatrix <- function(this.tree,thisMat)
{
	thisMat <- thisMat[is.finite(thisMat$bodyMass),]
	thisMat <- thisMat[rownames(thisMat) %in% this.tree$tip.label,]
	#thisMat <- thisMat[this.tree$tip.label,]
	# cbind(this.tree$tip.label, rownames(thisMat))

	return(thisMat)	
}

############################################################################################################################################
#functions to painttree in various ways; mostly used in troubleshooting and checking that tree are constructed properly
paintWildcardClades <- function(tree)
{
	this.tree <- tree
#	this.tree <- tree_resolv
	
	this.map <-  this.tree #artio.tree

	#plotTree(this.map,node.numbers=T,pts=F,fsize=0.1, ftype = "i")
	
	#cols<-c("black","blue","red","green","cyan","magenta","yellow","gray","darkgoldenrod","hotpink","midnightblue","sienna","skyblue","seagreen","salmon","maroon","khaki","cadetblue","burlywood","chocolate","aquamarine","black","blue","red","green","cyan","magenta","yellow","gray","darkgoldenrod","hotpink","midnightblue","sienna","skyblue","seagreen","salmon","maroon","khaki","cadetblue","burlywood","chocolate","aquamarine"); names(cols)<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42)
	cols <- c("black", "cyan","magenta", "blue", "seagreen", "chocolate", "darkgray", "green", "red"); names(cols) <- c("1","5","6","7","14","20","29","30","33")
	
	cols

	#Achaenodontidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Parahyus_vagus", sp2="Achaenodon_fremdi"),state="5")
	#plotSimmap(this.map, cols,fsize=0.1, ftype = "i")

	#Agriochoeridae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Diplobunops_crassus", sp2="Protoreodon_transmontanus"),state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Diplobunops_crassus", sp2="Agriochoerus_antiquus"),state="6")
	#plotSimmap(this.map,cols,fsize=0.1, ftype = "i")

	#Merycoidontidae
	#polytomy within the Merycoidodontidae leads to multiple different configurations that require multiple different configs fpr painting
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Bathygenys_reevesi", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Aclistomycter_middletoni", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Miniochoerus_affinis", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Limnenetes_platyceps", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Aclistomycter_middletoni", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Miniochoerus_affinis", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Limnenetes_platyceps", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Limnenetes_platyceps", sp2="Miniochoerus_affinis"), state="6")
	#plotSimmap(this.map, cols, fsize=0.1, ftype = "i")

	#Anthracotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Heptacodon_yeguaensis", sp2="Aepinacodon_deflectus"),state="7")
	#plotSimmap(this.map, cols,fsize=0.1, ftype = "i")
	
	#Entelodontidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_palustris", sp2="Dinohyus_hollandi"), state="14")
	#this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Choerodon_calkinsi", sp2="Archaeotherium_palustris"), state="14")
#	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_minimus", sp2="Choerodon_caninus"), state="14")
#	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_minimus", sp2="Brachyhyops_viensis"), state="14")
#	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_minimus", sp2="Boochoerus_angustus"), state="14")
#	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_minimus", sp2="Megachoerus_latidens"), state="14")
	#plotSimmap(this.map, cols, fsize=0.1, ftype = "i")
	
	#Protoceratidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Leptotragulus_ultimus", sp2="Prosynthetoceras_orthrionanus"), state="20")
	#plotSimmap(this.map, cols,fsize=0.1, ftype = "i")
	
	#Brontotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Menops_marshi", sp2="Parvicornus_occidentalis"), state="29")
	#this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Menops_bakeri", sp2="Eotitanops_minimus"), state="29")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Megacerops_osborni", sp2="Brontops_dispar"), state="29")
	#plotSimmap(this.map, cols, fsize=0.1, ftype = "i")

	#Chalicotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Moropus_merriami", sp2="Eomoropus_anarsius"), state="30")
	#plotSimmap(this.map,cols, fsize=0.1, ftype = "i")
	
	#Heptodon
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Heptodon_posticus", sp2="Heptodon_calciculus"), state="33")
	plotSimmap(this.map, cols, fsize=0.1, ftype = "i")
	
	return(this.map)
}

paintWildcardClades_Dated <- function(tree)
{
	this.tree <- tree
	#this.tree <- tree_resolv
	#this.tree <- tree_dated
	this.map <-  this.tree #artio.tree

	#plotTree(this.map,node.numbers=T,pts=F,fsize=0.1, ftype = "i")
	
	cols<-c("black","blue","red","green","cyan","magenta","yellow","gray","darkgoldenrod","hotpink","midnightblue","sienna","skyblue","seagreen","salmon","maroon","khaki","cadetblue","burlywood","chocolate","aquamarine","black","blue","red","green","cyan","magenta","yellow","gray","darkgoldenrod","hotpink","midnightblue","sienna","skyblue","seagreen","salmon","maroon","khaki","cadetblue","burlywood","chocolate","aquamarine"); names(cols)<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42)
	cols

	#Achaenodontidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Parahyus_vagus", sp2="Achaenodon_fremdi"),state="5")
	#plotSimmap(this.map, cols,fsize=0.1, ftype = "i")

	#Agriochoeridae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Diplobunops_crassus", sp2="Protoreodon_transmontanus"),state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Diplobunops_crassus", sp2="Agriochoerus_antiquus"),state="6")
	#plotSimmap(this.map,cols,fsize=0.1, ftype = "i")

	#Merycoidontidae
	#polytomy within the Merycoidodontidae leads to multiple different configurations that require multiple different configs fpr painting
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Bathygenys_reevesi", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Aclistomycter_middletoni", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Miniochoerus_affinis", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Limnenetes_platyceps", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Aclistomycter_middletoni", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Miniochoerus_affinis", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Limnenetes_platyceps", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Limnenetes_platyceps", sp2="Miniochoerus_affinis"), state="6")
	#plotSimmap(this.map, cols, fsize=0.1, ftype = "i")

	#Anthracotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Heptacodon_yeguaensis", sp2="Aepinacodon_deflectus"),state="7")
#	plotSimmap(this.map, cols,fsize=0.1, ftype = "i")
	
	#Entelodontidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_minimus", sp2="Dinohyus_hollandi"), state="14")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_minimus", sp2="Choerodon_caninus"), state="14")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_minimus", sp2="Brachyhyops_viensis"), state="14")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Archaeotherium_minimus", sp2="Megachoerus_latidens"), state="14")
#	plotSimmap(this.map, cols, fsize=0.1, ftype = "i")
	
	#Protoceratidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Leptotragulus_ultimus", sp2="Prosynthetoceras_orthrionanus"), state="20")
	plotSimmap(this.map, cols,fsize=0.1, ftype = "i")
	
	#Brontotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Menops_bakeri", sp2="Parvicornus_occidentalis"), state="29")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Menops_bakeri", sp2="Megacerops_osborni"), state="29")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Parvicornus_occidentalis", sp2="Megacerops_osborni"), state="29")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Parvicornus_occidentalis", sp2="Brontops_dispar"), state="29")
#	plotSimmap(this.map, cols, fsize=0.1, ftype = "i")

	#Chalicotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Moropus_merriami", sp2="Eomoropus_anarsius"), state="30")
#	plotSimmap(this.map,cols, fsize=0.1, ftype = "i")
	
	#Heptodon
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Heptodon_posticus", sp2="Heptodon_calciculus"), state="33")
	plotSimmap(this.map, cols, fsize=0.1, ftype = "i")
	
	return(this.map)
}

paintOrders <- function(tree)
{
	#this.map <- tree.current
	this.map <-  tree

	print(class(this.map))

	#cols<- c(palette(rainbow(23)));names(cols)<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
	cols<-c("black","red"); names(cols)<-c("Perissodactyla","Artiodactyla")
	print(cols)

	#Artiodactyla
#	node <- findMRCA(this.map, tips = c("Achaenodon_robustus","Antiacodon_venustus"))
#	print(node)
#	artio <- clade.members(node, this.map, tip.labels = TRUE)
#	print(artio)
#	node <- findMRCA(this.map, tips = artio)
	artio <- ungulateOrderSpeciesList(artio = TRUE, perisso = TRUE, anc.clade = "perisso")
	
	print(class(this.map))

	this.map <- paintSubTree(tree=this.map, node=findMRCA(tree=this.map, tips = artio),anc.state="1", state="2")

#	print(class(this.map))
	
	#this.map <- paintSubTree(tree=this.map, node=findMRCA(this.map, tips = c("Achaenodon_robustus","Antiacodon_venustus")),anc.state="1", state="2")

	#Perissodactyla
	#this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.tree, sp1="Heptodon_posticus", sp2="Notiotitanops_mississippiensis"), state="3")

	plotSimmap(this.map, cols, fsize=0.1, ftype = "i", lwd = 1)
	# nodelabels(cex = 0.1)

	return(this.map)
}

paintTemporalRange <- function(tree, maxAge, minAge, int.length)
{
	periods <- rev(seq(minAge, maxAge, by=int.length))
  	names(periods)<-c(periods)
 	periods

	maxAge <- max(periods)

	#intervals <- makeIntervals(startDate=maxAge, endDate=minAge, intervalLength = int.length)
	#names(intervals) <- c(print(intervals))
	#intervals
	
	cols <- palette(rainbow(periods))
	
	this.map <- make.era.map(tree, limits=c(0,maxAge-periods))
	#plotSimmap(this.map, cols, fsize=0.1, ftype = "i", lwd = 1) #for some reason is not plotting phylogeny with colored intervals

	return(this.map)
}

############################################################################################################################################

#function performs BM model on given set of temporal breaks.  Failure of moddel triggers a reassignemnts of the sigma values
getLkMVMorphOneIntervalSet <- function(newBreak=NULL, oldIntBreaks=NULL, this.tree, dat, this.intervals, analysis.settings = list(do.BMsimple = FALSE, do.BMM = FALSE, do.OUsimple = FALSE, do.OUM = FALSE, do.heuristic = FALSE,adjust.date.after=10, adjust.date.increment = 0.1, this.constraint = FALSE, do.parallel = FALSE)) {
    ###############
    #this.tree <- tree.list[[1]]
    #dat <- dat.vec
    #topIntv <- min(which(intervals$ageBase > min(this.tree$node.date)))
    #baseIntv <- max(which(intervals$ageTop < max(this.tree$node.date[seq_along(this.tree$tip.label)])))
    #this.intervals <- intervals[topIntv:baseIntv,]
    #newBreak <- c(1,2,3,8)
    #oldIntBreaks <- NULL
    # newBreak <- breakList[[1]]
    ################
   
    intBreaks <- sort(c(newBreak, oldIntBreaks), decreasing=TRUE) #interval numbers
		cat("Breaks:", intBreaks, "\n")
    break.dates <- sort(as.numeric(this.intervals[intBreaks, "ageBase"]), decreasing=TRUE) #actual dates of interval/break
    this.limits <- max(this.tree$node.date) - c(max(this.tree$node.date), break.dates) # puts in time from root
    thisMap <- make.era.map(this.tree, limits=this.limits)
    #plotSimmap(thisMap,fsize=0.5)
    
    #cat("Start models","\n")
    # add for loop to loop doOneMultipleModel through all requested models
    if (analysis.settings$do.BMM) { master.BM <- doOneMultipleModel(this.tree = thisMap, dat = dat, intBreaks = intBreaks, break.dates = break.dates, this.model = "BMM", analysis.settings = analysis.settings)
  #  cat("BM complete", "\n")
    } else master.BM <- NULL
      
		if(analysis.settings$do.OUM) { master.OU <- doOneMultipleModel(this.tree = thisMap, dat = dat, intBreaks = intBreaks, break.dates = break.dates, this.model = "OUM", analysis.settings = analysis.settings)
	#	cat("BM complete", "\n")  
		} else master.OU <-  NULL 
    
    master.result <- list(bm=master.BM,ou= master.OU)
    
    return(master.result)
}

getRandomMatrix <- function(n.row, this.constraint=FALSE) {
    this <- matrix(NA, nrow=n.row, ncol=n.row)
    this[upper.tri(this)] <- runif(n=length(this[upper.tri(this)]))
    this[lower.tri(this)] <- t(this)[lower.tri(t(this))]
    diag(this) <- runif(n=n.row, min=1, max=10)
    chol_factor <- chol(this)
    if(this.constraint==FALSE) chol_factor[upper.tri(chol_factor, diag=TRUE)] else chol_factor[upper.tri(chol_factor, diag=FALSE)]
    # this[upper.tri(this, diag=TRUE)]
}

doOneMultipleModel <- function(this.tree, dat, intBreaks, break.dates, this.model= list("BMM", "OUM"), analysis.settings = list(do.BMsimple = FALSE, do.BMM = FALSE, do.OUsimple = FALSE, do.OUM = FALSE, do.heuristic = FALSE,adjust.date.after=10, adjust.date.increment = 0.1, this.constraint = FALSE, do.parallel = FALSE)) {
	 flag <- FALSE
	 this.param.BM <- list(constraint = analysis.settings$this.constraint, smean = TRUE, trend = FALSE)    #default values to begin analysis
	 this.param.OU <- list(sigma = NULL, alpha = NULL, vcv = "fixedRoot", root=TRUE, decomp = c("cholesky", "spherical", "eigen", "qr", "diagonal", "upper", "lower")) 
	 
    counter <- 0
    while(!flag) {
        counter <- counter + 1 # list of attempts
        if (counter %% analysis.settings$adjust.date.after == 0) {
            print("***** Adjusting break dates *****\n")
        		print(counter)
            break.dates <- break.dates + analysis.settings$adjust.date.increment
            this.limits <- max(this.tree$node.date) - c(max(this.tree$node.date), break.dates)
            this.tree <- make.era.map(this.tree, limits=this.limits)
        }
        result <- NULL
	    if(this.model == "BMM") {
		    try(result <- mvBM(tree=this.tree, 
                           data=dat, 
                           model= "BMM",
                           param=this.param.BM,
                           # method="rpf",
                           control=list(maxit=1e8), 
                           # optimization="subplex",
                           diagnostic=FALSE, 
                           echo=FALSE))
			result$model <- "BMM"
		} else {
			try(result <- mvOU(tree=this.tree, 
                           data=dat, 
                           model= "OUM",
                           param=this.param.OU,
                           # method="rpf",
                           control=list(maxit=1e8), 
                           # optimization="subplex",
                           diagnostic=FALSE, 
                           echo=FALSE))
			result$model <- "OUM"
		}
        if (is.null(result)) {
            print("***** NULL Result: Adjusting starting parameter values *****\n")
           	if(this.model == "BMM"){
           		print("BMM")
	            	if (analysis.settings$this.constraint==FALSE) {  this.param.BM$sigma <- replicate(n=length(intBreaks)+1,expr= getRandomMatrix(n.row=ncol(dat), this.constraint=analysis.settings$this.constraint), simplify=FALSE)
	            	} else this.param.BM$sigma <- replicate(n=length(intBreaks)+1,expr= c(getRandomMatrix(n.row=ncol(dat), this.constraint=analysis.settings$this.constraint), runif(1)), simplify=FALSE)
            } else {
	            	print("OU")
	            	this.param.OU$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint=analysis.settings$this.constraint)
         		this.param.OU$alpha <- getRandomMatrix(n.row=ncol(dat), this.constraint=analysis.settings$this.constraint)
            }
        } else if (result$convergence != 0) { 
            print("***** Failure to converge : Adjusting starting parameter values *****\n")
        		cat("result$convergence", result$convergence, "\n")
            if(this.model == "BMM") {
	            	print("BMM")
	            	if (analysis.settings$this.constraint==FALSE) { this.param.BM$sigma <- replicate(n=length(intBreaks)+1,expr= getRandomMatrix(n.row=ncol(dat), this.constraint=analysis.settings$this.constraint), simplify=FALSE)
	           		#print(this.param.BM$sigma)
	            		} else this.param.BM$sigma <- replicate(n=length(intBreaks) +1, expr= c(getRandomMatrix(n.row=ncol(dat), this.constraint=analysis.settings$this.constraint)), simplify=FALSE)	
	            	#print(this.param.BM$sigma)
            	} else {
	            	print("OUM")
	            	this.param.OU$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint=analysis.settings$this.constraint)
          		this.param.OU$alpha <- getRandomMatrix(n.row=ncol(dat), this.constraint=analysis.settings$this.constraint)
            }
        # } else if (result$hess.values != 0) { 
            # print("***** Failure to HESS: Adjusting starting parameter values *****\n")
            # if (this.constraint==FALSE) { this.param$sigma <- replicate(n=length(intBreaks)+1,expr=as.double(runif(ncol(dat)*(ncol(dat)+1)/2)), simplify=FALSE)
            # } else this.param$sigma <- replicate(n=length(intBreaks)+1,expr=as.double(runif(n=(ncol(dat)*(ncol(dat)-1)/2)+1)), simplify=FALSE)
        } else flag=TRUE
    }
   	output <- result[c("AICc")] #, "theta", "alpha", "sigma")]
    #output$intBreaks <- intBreaks
    # result$break.dates <- break.dates
    #result$map <- thisMap
	return(output)
}

testRateShiftsMVMorphIntervals <- function(tree, dat, intervals, analysis.settings = list(do.BMsimple = FALSE, do.BMM = FALSE, do.OUsimple = FALSE, do.OUM = FALSE, do.heuristic = FALSE,equal.interval = TRUE,this.constraint = FALSE,adjust.date.after=10, adjust.date.increment = 0.1, do.parallel = FALSE)) {
    ###############
    #tree <- tree.list[[1]]
    #dat <- dat.vec
    ################
    
   # rez.x <- list()
    
    # print("Drop Tip Check")
    # print(class(tree))
    # print(tree)
    if (length(dim(dat)) == 1) {  #if dat is vector (1 trait)
        this.tree <- dropTipKeepDates(tree[[1]], tree$tip.label[!tree$tip.label %in% names(dat)], min.bl=1.0)
        dat <- dat[names(dat) %in% this.tree$tip.label]
        dat <- dat[match(this.tree$tip.label, names(dat))]   #reorders dat to match the taxon order of tree - MUST stay in this function, because each tree will have different taxon order...
    } else {
        this.tree <- dropTipKeepDates(tree, tree$tip.label[!tree$tip.label %in% rownames(dat)], min.bl=1.0)
        dat <- matrix(dat[this.tree$tip.label,], ncol=1, dimnames=dimnames(dat)) # this *should* replace the following code
        
        #row.dat <- rownames(dat)[which(rownames(dat) %in% this.tree$tip.label)]
        #col.dat <- names(dat)
        #dat <- as.data.frame(dat[rownames(dat) %in% this.tree$tip.label, ])
        #rownames(dat) <- row.dat
        #row.dat <- rownames(dat)[match(this.tree$tip.label, rownames(dat))]
        #dat <- as.data.frame(dat[match(this.tree$tip.label, rownames(dat)),])  #reorders dat to match the taxon order of tree - MUST stay in this function, because each tree will have different taxon order...
        #rownames(dat) <- row.dat
        #names(dat) <- col.dat   
        # class(dat)
        #need to find way to preserve rownames
    }   
    
    print("Completed Drop Tip Check: Starts interval generation")
    
    topIntv <- min(which(intervals[,"ageBase"] > min(this.tree$node.date))) # node.date is non longer an attribute on the tree
    baseIntv <- max(which(intervals[,"ageTop"] < max(this.tree$node.date[seq_along(this.tree$tip.label)])))
    this.intervals <- intervals[topIntv:baseIntv,]
    n.intv <- nrow(this.intervals)
    optList <- list()
  	conv <- FALSE
  	#this.constraint <- FALSE
  	 
	 ######## get single rate BM
	 print("Get Single Rate BM")
	 if(analysis.settings$do.BMsimple) {    
	  	param.list <- list(constraint = analysis.settings$this.constraint)
	  	while(!conv) {
		  	x <- mvBM(tree=this.tree, data=dat, model="BM1", param=param.list, echo=FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
		    # y <- getLkMVMorphOneIntervalSet(newBreak=11, oldIntBreaks, this.tree=this.tree, dat=dat, this.intervals=this.intervals, adjust.date.after=20) #subsequent models
		    if (x$convergence != 0) {
		    		print("***** Failure to converge: Adjusting starting parameter values *****\n")
		            print(x$convergence)
		            param.list$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint= analysis.settings$this.constraint)
		            print(param.list$sigma)
		    } else conv <-TRUE
	    }
	    print("Intitial BM model complete")
	    x <- x$AICc
	    names(x) <- "AICc"
	    #x$nrates <- 1
	    #x$optBreaks <- NULL
	    #x$model <- "BM"
	}
	
	 ######## get single rate OU
	conv <- FALSE
	print("Get Single Rate OU") 
	 if(analysis.settings$do.OUsimple) { 
	 	this.param <- list(sigma = NULL, alpha = NULL, vcv = "fixedRoot", decomp = c("cholesky", "spherical", "eigen", "qr", "diagonal", "upper", "lower")) 
	 	while(!conv) {
		  	y <- mvOU(tree=this.tree, data=dat, model="OU1", param=param.list, echo=FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
		    if(y$convergence != 0) {
		    		print("***** Failure to converge: Adjusting starting parameter values *****\n")
		            print(y$convergence)
		            param.list$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint= analysis.settings$this.constraint)
		            print(param.list$sigma)
		    } else conv <-TRUE
	    }
	    print("Intitial OU model complete")
	    y <- y$AICc
	    names(y) <- "AICc"
	    #y$nrates <- 1
	    #y$optBreaks <- NULL
	    #y$model <- "OU"
	}    
    optList[[1]] <- list(bm=x, ou=y)
    rm(x, y)

######## start multiple rate analyses    

    flag <- FALSE #get thrown when AIC starts to increase; must be fewer rates than intervals
    nrates <- 2
    oldIntBreaks <- NULL 
    # rez.list <- list()
    # optList.BM <- list()
    # optList.OU <- list()
    #optSave <- optList
    # optList <- optSave
    
    print("Start Multiple Rate Loop")
    while(!flag & nrates <= n.intv) {
        cat("Beginning analysis for ", nrates, " rates...\r")
        if (analysis.settings$do.heuristic) { # builds break list; possible shift points; adds on to last nrates and optimum
            if (nrates > 2) breakList <- as.integer(seq_len(n.intv - 1)[-oldIntBreaks]) else breakList <- as.integer(seq_len(n.intv - 1)) # breaks are at the base of intervals
            breakList <- sample(breakList, size=length(breakList))   #randomly reorders breakList to avoid solidifying the "early" breaks
        } else {
            breakList <- listifyMatrixByColumn(m=combn(x=nrow(this.intervals), m=nrates-1)) # all possible combinations of shift points
        }
    		
   			#cat("Length Breaklist:", length(breakList),"\n")
   			
        if (analysis.settings$do.parallel) { 
         rez.list <- mclapply(breakList, getLkMVMorphOneIntervalSet, oldIntBreaks, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings= analysis.settings, mc.cores=4, mc.preschedule = FALSE) 
        #rez.list <- mclapply(breakList, getLkMVMorphOneIntervalSet, oldIntBreaks, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings= analysis.settings, mc.cores=detectCores(), mc.preschedule = FALSE)
		     print(rez.list)
		     } else rez.list <- lapply(breakList, getLkMVMorphOneIntervalSet, oldIntBreaks, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings= analysis.settings)

   			cat("Finished rez.list","\n")
        #need way to bypass NULL vlaues if in just one entry
        
        ### if (!is.null(BM) & lowest BM.AICc is > that last number of rates)  analysis.settings$do.BMM <- FALSE)
        ### if (!is.null(OU) & lowest OU.AICc is > that last number of rates)  analysis.settings$do.OUM <- FALSE)
        ### if(!is.null(BM) & is.null(OU)) flag <- TRUE
      #	if (min(sapply(rez.list, function(x) min(x$bm$AICc, x$ou$AICc)))  > optList[[nrates-1]]$bm$AICc) & optList[[nrates-1]]$ou$AICc) { flag <- TRUE}  #################################check function to break the loop
         
   				  if (analysis.settings$do.BMM) {
             	cat("RezList check elements BM", min(sapply(rez.list, function(x) x[["bm"]][["AICc"]])),"\n")
            	print(optList[[nrates-1]][["bm"]])
            	cat(min(sapply(rez.list, function(x) x[["bm"]][["AICc"]])) > optList[[nrates-1]][["bm"]][["AICc"]], "\n")
	         		
	         		
            	 #   if (min(sapply(rez.list, function(x) x$bm$AICc)) > optList[[nrates-1]]$bm$AICc) analysis.settings$do.BMM <- FALSE
            	if (min(sapply(rez.list, function(x) x[["bm"]][["AICc"]])) > optList[[nrates-1]][["bm"]][["AICc"]]) analysis.settings$do.BMM <- FALSE
             }
        # }
					cat("Check", analysis.settings$do.BMM,"\n")	
   				cat("Checked BM for opt rates","\n")
        # if (!is.null(rez.list$V1$ou)) {
        	#print("OU is not NULL")
   				
   			      if (analysis.settings$do.OUM) {
            	
   			      cat("RezList check elements OU", min(sapply(rez.list, function(x) x[["ou"]][["AICc"]])),"\n")
            	print(optList[[nrates-1]][["ou"]])
            
            		#   if (min(sapply(rez.list, function(x) x$ou$AICc)) > optList[[nrates-1]]$ou$AICc) analysis.settings$do.BMM <- FALSE
            		if (min(sapply(rez.list, function(x) x[["ou"]][["AICc"]])) > optList[[nrates-1]][["ou"]][["AICc"]]) analysis.settings$do.OUM <- FALSE
            	}
        # }
					cat("Checked OU for opt rates","\n")
						
        if (!analysis.settings$do.BMM & !analysis.settings$do.OUM) flag <- TRUE
        
        # cat("\n flag", flag)
        # cat("\n do.BMM", analysis.settings$do.BMM)
        # cat("\n do.OUM", analysis.settings$do.OUM)
        # print("\n")
       	cat("get opt model", "\n")
        if (!flag) {
        		cat("get opt BM","\n")  
        		if (analysis.settings$do.BMM) {
	            	optSet.BM <- which.min(sapply(rez.list, function(x) x[["bm"]][["AICc"]])) # will need to dig around for beg AIC; find model model for a given best combo of breaks
	            	#rez.list[[optSet.BM]]$bm$nrates <- nrates
	            	#rez.list[[optSet.BM]]$bm$opt.dates <- this.intervals[rez.list[[optSet.BM]]$bm$intBreaks, "ageBase"]
	            	optList.BM <- rez.list[[optSet.BM]][["bm"]]
	            	optList.BM$BreakList.index <- which(breakList[[optSet.BM]]) 
	            	if (analysis.settings$do.heuristic) oldIntBreaks[["bm"]] <- as.integer(optList[[nrates]][["intBreaks"]])
            } else optList.BM <- NULL 
          		 cat("get opt BM","\n")
            	if (analysis.settings$do.OUM) {
	            	optSet.OU <- which.min(sapply(rez.list, function(x) x[["ou"]][["AICc"]])) # will need to dig around for beg AIC; find model model for a given best combo of breaks
	            	#rez.list[[optSet.OU]]$ou$nrates <- nrates
	            	#rez.list[[optSet.OU]]$ou$opt.dates <- this.intervals[rez.list[[optSet.OU]]$ou$intBreaks, "ageBase"]
	            	optList.OU <- rez.list[[optSet.OU]][["ou"]]
	            	optList.OU$BreakList.index <- which(breakList[[optSet.OU]]) 
            	if (analysis.settings$do.heuristic) oldIntBreaks[["ou"]] <- as.integer(optList[[nrates]][["intBreaks"]])
       	  	} else optList.OU <- NULL 
	        optList[[nrates]] <- list(bm=optList.BM, ou=optList.OU)
        }
        
        # rm(rez.list)
        # gc()
        # rez.x <- append(rez.x, list(rez.list)) # will contain one extra list than optList due to the flag triggering before rez.list creation
        
         # save(rez.x, rez.list, optSet.BM, optSet.OU, optList.BM, optList.OU,optList, nrates, breakList, oldIntBreaks, topIntv, baseIntv, this.intervals, n.intv, flag, file ="/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/2017_5_24_EvoModel_BM_Param_NoNhueristic.RData")
        # cat("\n saved rate", nrates)
        optList.name <- paste("EvoAnalysis_Optlist_",paste("nrates_",nrates, sep=""),sep="")
        optList.name <- paste(paste(optList.name,timestamp(),sep=""),".RData",sep="")
        cat(optList.name,"\n")
        save(optList, breakList, rez.list, file = optList.name) #rez.list included for testing and debugging
        nrates <- nrates + 1
      }
  		# rez.nam <- paste("Nrates", c(2:(nrates-1)))
    	# names(rez.x) <- rez.nam
  
   # final.list <- list(optList, rez.x)
   # names(final.list) <- c("optList", "rez.list")
  
   return(optList)
}

getTreeOptList <- function(tree.list, dat, intervals, analysis.settings = list(do.BMsimple = FALSE, do.BMM = FALSE, do.OUsimple = FALSE, do.OUM = FALSE, do.heuristic = FALSE, adjust.date.after=10, adjust.date.increment = 0.1, do.parallel = FALSE)) {
	if (class(tree.list)=="phylo") { results <- testRateShiftsMVMorphIntervals(tree=tree.list, dat = dat, intervals = intervals, analysis.settings = analysis.settings) 
	} else results <- lapply(tree.list, FUN=testRateShiftsMVMorphIntervals, dat = dat, intervals = intervals, analysis.settings = analysis.settings)
	# for (i in seq_along(tree.list)) results[[i]] <- testRateShiftsMVMorphIntervals (tree.list[[i]], dat = dat, intervals = intervals, analysis.settings = analysis.settings)
		
	#tree.optList <- list()
	#tree.num <- 1
	
		
	#if(class(tree.list) == "phylo") { treeTotal <- 2
	#} else treeTotal <- length(tree.list) + 1
		
	#while(tree.num < treeTotal){ #loop to go through all generated trees
			# print("Tree")
			# print(tree.num)
			#cat("Tree", tree.num)
#
			#if (class(tree.list) == "phylo") { #run for single tree
				#HueristicTime <- system.time(Hueristicresults <- testRateShiftsMVMorphIntervals(tree = tree.list, dat = dat, intervals = intervals, occs = occs, do.heuristic=do.heuristic, do.parallel=do.parallel))
			#	print("Models Completed")
			#} else {
			#	HueristicTime <-system.time(Hueristicresults <- testRateShiftsMVMorphIntervals(tree = tree.list[[tree.num]], dat = dat, intervals = intervals, occs= occs, do.heuristic=do.heuristic, do.parallel=do.parallel)) # run if multiple trees
			
			#	print("Models Completed")
			#}
			#
			#print("Hueristicresults")
			#print(Hueristicresults)
			#optSetBM <- which.min(sapply(Hueristicresults[['optList']], function(x) x$AICc)) # will need to dig around for beg AIC; find model model for a given best combo of breaks
			#print("optSetBM")
			#print(optSetBM)
			##print("OUHueristicresults")
			#print(OUHueristicresults)
			#optSetOU <- which.min(sapply(OUHueristicresults[['optList']], function(x) x$AICc)) # will need to dig around for beg AIC; find model model for a given best combo of breaks
          #  print("optSetOU")
           ## print(optSetOU)
           
           # if(Hueristicresults[['optList']][[optSetBM]]$AICc < OUHueristicresults[['optList']][[optSetOU]]$AICc) {optList[[tree.num]] <- Hueristicresults[['optList']][[optSetBM]]
           # 	optList[[tree.num]]$model <- "BM" 
            #	print("BM better model")}
            
           # if(OUHueristicresults[['optList']][[optSetOU]]$AICc < Hueristicresults[['optList']][[optSetBM]]$AICc) {optList[[tree.num]] <- OUHueristicresults[['optList']][[optSetOU]]
            #	optList[[tree.num]]$model <- "OU"
            ##	print("OU better model")}
          
        
          #need to save or output list that contains all iterations of models#may need to input check for intervals to be included in filename
          # filename_save <- paste("/Users/evandoughty/Dropbox/ungulate_RA/RCode/EvoModelList_All_Results/EvoModelResultsList_Hueristics=",analysis.settings$do.heuristic,sep="")
          	# filename_save <- paste(filename_save, "_do.parallel=", sep="")
			# filename_save <- paste(filename_save,analysis.settings$do.parallel, sep="")
			# filename_save <- paste(filename_save, "_do.BMsimple=", sep = "")
			# filename_save <- paste(filename_save, analysis.settings$do.BMsimple, sep ="")
			# filename_save <- paste(filename_save, "_do.BMM=", sep="") 
			# filename_save<- paste(filename_save, analysis.settings$do.BMM, sep="")
			# filename_save <- paste(filename_save, "_do.OUsimple=", sep = "")
			# filename_save <- paste(filename_save, analysis.settings$do.OUsimple, sep ="")
			# filename_save <- paste(filename_save, "_do.OUMp=", sep="")
			# filename_save <- paste(filename_save, analysis.settings$do.OUM, sep="")
	
          # filename_save <- paste(filename_save,".RData", sep="")
          # save(results, file= filename_save)
          
         # 	tree.num <- tree.num + 1
          # }

	#return(tree.optList)
	return(results)
}

checkBreaks <- function(breakList, tree, occs, nrates, intervals)
{

######
# get species within intervals
#####
	col.dates <- getCollectionAgesFromOccs(occs=occs[, c("collection_no", "max_ma", "min_ma")], random=TRUE)
	occDates <- col.dates$collection_age[match(occs$collection_no, col.dates$collection_no)]
	intOccs <- apply(intervals, 1, function(thisIntv) occs$occurrence_no[occDates > thisIntv[1] & occDates <= thisIntv[2]])
	# intTaxa <- sapply(intOccs, function(x) unique(occs$accepted_name[occs$occurrence_no %in% x]))
	################
	int_NoTax <- matrix(nrow=3,ncol=length(rownames(intervals)))
	colnames(int_NoTax) <- rownames(intervals)
	#colnames(int_NoTax) <- intervals[,2]
	int_NoTax[1,] <- !(sapply(intOccs, length)) #if FALSE then has taxa within interval
	int_NoTax[2,] <- intervals[,1]
	int_NoTax[3,] <- intervals[,2]
	################
	
	   taxInBreaks <- int_NoTax
       rownames(taxInBreaks) <- c("TaxonMissing","ageTop","ageBase")
       
       intBreaks <- sort(unlist(breakList), decreasing=TRUE) #interval numbers; only allows atomic values when ran manually
       break.dates <- sort(c(as.numeric(intervals[intBreaks, "ageBase"])), decreasing=TRUE)
       
       break.intDates <- matrix(nrow=2,ncol=length(intBreaks))
       break.intDates[1,] <- intBreaks
       break.intDates[2,] <- break.dates
       colnames(break.intDates) <- rev(colnames(taxInBreaks))
       rownames(break.intDates) <- c("breakList","ageBase")
		
		maxDate <- tree$root.date
		#remove breaks that are older than max age of tree
		
	# Need way to compare each entry of breaklist to 
		names <- colnames(taxInBreaks)
		taxInBreaks <- rbind(taxInBreaks, as.data.frame(break.intDates))
		validBreaks <- rev(taxInBreaks[,taxInBreaks[1,] == 0])
        
        validBreakList <- lapply(breakList)
        
        if(nrates == 1) {breakList <- apply(validBreaks["breakList",], 2, as.list())}

	
return(int_NoTax)
}

test_code_Clavel <- function()
{
	# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Simulate the traits
alpha<-matrix(c(2,0.5,0.5,1),2)
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(2,3,1,1.3)
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, theta=theta, names_traits=c("head.size","mouth.size")), model="OUM", nsim=1)

## Fitting the models

# OUM - Analysis with multiple optima
mvOU(tree, data)

# with starting values provided
mvOU(tree, data, param=list(alpha=alpha, sigma=sigma)) # use the generating matrices as starting matrices
mvOU(tree, data, param=list(alpha=alpha)) # just alpha
mvOU(tree, data, param=list(sigma=sigma)) # just sigma




# Try with univariate data
alpha <- 2
sigma <- 0.5
theta <- c(0,5) # different theta for each regimes
# generate some data
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha,  theta=theta), model="OUM")

# fit the model with default starting values
mvOU(tree, data)

# try with yours starting values
mvOU(tree, data, param=list(alpha=alpha, sigma=sigma)) # use the generating matrices as starting matrices
mvOU(tree, data, param=list(alpha=alpha)) # just alpha
mvOU(tree, data, param=list(sigma=sigma)) # just sigma

sigma <- replicate(n=length(intBreaks)+1,expr= getRandomMatrix(n.row=ncol(dat), this.constraint=this.constraint), simplify=FALSE)
alpha <- replicate(n=length(intBreaks)+1,expr= getRandomMatrix(n.row=ncol(dat), this.constraint=this.constraint), simplify=FALSE)

return
}

getRefMeas <- function(measure.Mat,treeTips) #get lists of references for occurence data
{
  specimenMat <- getSpecimenMatFromMeasurements(filename="https://dl.dropbox.com/s/423x0zn3mpxuwc7/specimens.csv")
  specimenMat <- merge(specimenMat, getBirlenbachBlastoSpecimens(filename="https://dl.dropbox.com/s/943dq4lb1kd1h1e/blastoBirlenbach20140207.csv", info.file="https://dl.dropbox.com/s/nhrqzrclwtr5c8e/info2.csv"), all=TRUE)
  specimenMat <- merge(specimenMat, getLiteratureSpecimenMat(filename="https://dl.dropbox.com/s/qef8ts9a73j5ukb/literature.csv"), all=TRUE)
  
  specimenMat[sapply(specimenMat, is.nan)] <- NA
  specimenMat$species <- getCurrentTaxa(specimenMat$species)
  
  Refs <- unique(specimenMat$Reference)
  write.csv(Refs,"MeasurementReferences.csv")
  return(Refs)
}

#Function to search thisMat with crownheight/hypsodonty list
matchBM.CrownHeight <- function(thisMat)
{
	crownSp <- read.table("~/Dropbox/ungulate_RA/CrownHeight_OrigExtinct_Tests/Janis_2000_crownSp.txt")
	crownSp$V1 <- gsub(pattern = "_", replacement = " ", x = crownSp$V1)
	crownSp$V1 <- getCurrentTaxa(crownSp$V1)
	crownSp$V1 <- gsub(pattern = "[[:space:]]", replacement = "_", x = crownSp$V1)
	crownSp$V1 <- gsub(pattern="\"", replacement= "", crownSp$V1)
	bmList <- as.data.frame(thisMat$bodyMass)
	bmList[,2] <- rownames(bmList)
	names(bmList) <- c("bodyMass","taxon")
	names(crownSp) <- c("taxon","Crown Height")
	
	bmCrown<- merge(bmList,crownSp,by="taxon")
	
	missingBm <- as.data.frame(bmList[!bmList[,"taxon"] %in% bmCrown[,"taxon"],])
	missingBm[,3] <- "NA"
	names(missingBm) <- c("bodyMass","taxon","Crown Height")
	missingBm <- missingBm[,c("taxon","bodyMass","Crown Height")]
	rownames(missingBm) <- c()
	
	missingCrown <- as.data.frame(crownSp[!crownSp[,"taxon"] %in% bmCrown[,"taxon"],])
	missingCrown[,3] <- "NA"
	names(missingCrown) <- c("taxon","Crown Height", "bodyMass")
	missingCrown <- missingCrown[,c("taxon","bodyMass","Crown Height")]
	
	missingTax <- rbind(missingBm, missingCrown)
	
	missingTax <- missingTax[order(missingTax$taxon),]
	
	#write.csv(missingTax,file="~/Dropbox/ungulate_RA/CrownHeight_OrigExtinct_Tests/BmCrownHeight_MissingTaxa.csv")
	# need to impliment a check to replace names with synonyms
	return(missingTax)
}

orderEvoModel <- function(tree, occs, dat, clade = NULL, state.clade = 2,
													analysis.settings = list(do.BMsimple = FALSE, do.BMM = FALSE, 
																												do.OUsimple = FALSE, do.OUM = FALSE, 
																												do.heuristic = FALSE,this.constraint = FALSE,
																												do.parallel = FALSE))
																	#clade1 = "Artiodactyla", char.state = 2, analysis.settings)
{
	occs.cut<-occs[occs$order %in% clade,]
	occs.cut <- occs.cut[occs.cut$accepted_name %in% tree$tip.label,]
	
	this.tree <- paintSubTree(tree = tree, node = getMRCA(tree, tip = occs.cut$accepted_name), 
													 state = state.clade, anc.state = "1")
	
	cols <- c("black", "red"); names(cols) <- c("1", "2")
	
	plotSimmap(this.tree,cols, fsize = 0.05)
	optList <- list()

	conv <- FALSE
	if(analysis.settings$do.BMsimple) {    
		param.list <- list(constraint = analysis.settings$this.constraint)
		while(!conv) {
			x <- mvBM(tree=this.tree, data=dat, model="BM1", param=param.list, echo=FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
			# y <- getLkMVMorphOneIntervalSet(newBreak=11, oldIntBreaks, this.tree=this.tree, dat=dat, this.intervals=this.intervals, adjust.date.after=20) #subsequent models
			if(x$convergence != 0) {
				print("***** Failure to converge: Adjusting starting parameter values *****\n")
				print(x$convergence)
				param.list$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint= analysis.settings$this.constraint)
				print(param.list$sigma)
			} else conv <-TRUE
		}
		print("Intitial BM model complete")
		x$model <- "BM"
	}
	
	######## get single rate OU
	conv <- FALSE
	print("Get Single Rate OU") 
	if(analysis.settings$do.OUsimple) { 
		this.param <- list(sigma = NULL, alpha = NULL, vcv = "fixedRoot", decomp = c("cholesky", "spherical", "eigen", "qr", "diagonal", "upper", "lower")) 
		while(!conv) {
			y <- mvOU(tree=this.tree, data=dat, model="OU1", param=param.list, echo=FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
			if(y$convergence != 0) {
				print("***** Failure to converge: Adjusting starting parameter values *****\n")
				print(y$convergence)
				param.list$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint= analysis.settings$this.constraint)
				print(param.list$sigma)
			} else conv <-TRUE
		}
		print("Intitial OU model complete")
		y$model <- "OU"
	}    
	optList[[1]] <- list(bm=x, ou=y)
	rm(x, y)
	
	conv <- FALSE
	if(analysis.settings$do.BMM) {    
		param.list <- list(constraint = analysis.settings$this.constraint)
		while(!conv) {
			x <- mvBM(tree=this.tree, data=dat, model="BMM", param=param.list, echo=FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
			# y <- getLkMVMorphOneIntervalSet(newBreak=11, oldIntBreaks, this.tree=this.tree, dat=dat, this.intervals=this.intervals, adjust.date.after=20) #subsequent models
			if(x$convergence != 0) {
				print("***** Failure to converge: Adjusting starting parameter values *****\n")
				print(x$convergence)
				param.list$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint= analysis.settings$this.constraint)
				print(param.list$sigma)
			} else conv <-TRUE
		}
		print("Complex BM model complete")
		x$model <- "BM"
	}
		
	conv <- FALSE
	print("Get Single Rate OU") 
	if(analysis.settings$do.OUsimple) { 
	   this.param <- list(sigma = NULL, alpha = NULL, vcv = "fixedRoot", decomp = c("cholesky", "spherical", "eigen", "qr", "diagonal", "upper", "lower")) 
			while(!conv) {
				y <- mvOU(tree=this.tree, data=dat, model="OUM", param=param.list, echo=FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
				if(y$convergence != 0) {
					print("***** Failure to converge: Adjusting starting parameter values *****\n")
					print(y$convergence)
					param.list$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint= analysis.settings$this.constraint)
					print(param.list$sigma)
				} else conv <-TRUE
			}
			print("Complex OU model complete")
			y$model <- "OU"
		}    
		optList[[2]] <- list(bm=x, ou=y)
		rm(x, y)
	
	return(optList)
}

#function to go through list of optList to extract AICc, nrates, location of shifts, and parameters
getOptModels <- function(opt.list, do.plot = FALSE, folder = NULL) #function will save/check rez.list to find what rate configurations have been completed in order to continue in the event of a shutdown
{
	optTree <- list()
	
	#read in optList list or optList files
	
		#isolate optList for single tree
		for(nrates in seq(1, length(opt.list), 1)) { #print(opt.list[[nrates]])
			#acquire AICc value for BM and OU
			if(is.null(opt.list[[nrates]][["bm"]])){ opt.list[[nrates]][["bm"]][["AICc"]] <- 9999}
				#if(model$model == "OU"){
			if(is.null(opt.list[[nrates]][["ou"]])){ opt.list[[nrates]][["ou"]][["AICc"]] <- 9999}
		}
			optSet.OU <- which.min(sapply(opt.list, function(x) x[["ou"]][["AICc"]])) # will need to dig around for beg AIC; find model model for a given best combo of breaks
			optSet.BM <- which.min(sapply(opt.list, function(x) x[["bm"]][["AICc"]]))
					
			if(opt.list[[optSet.OU]][["ou"]][["AICc"]] <	opt.list[[optSet.BM]][["bm"]][["AICc"]]) optTree[[1]] <- opt.list[[optSet.OU]][["ou"]]
			
			if(opt.list[[optSet.OU]][["ou"]][["AICc"]] >	opt.list[[optSet.BM]][["bm"]][["AICc"]]) optTree[[1]] <- opt.list[[optSet.BM]][["bm"]] 
	
	return(optTree)
}
	
getOptBM <- function(opt.list) #function will save/check rez.list to find what rate configurations have been completed in order to continue in the event of a shutdown
{
	optTree <- list()
	
	#read in optList list or optList files
	
	#isolate optList for single tree
	for(nrates in seq(1, length(opt.list), 1)) { #print(opt.list[[nrates]])
		#acquire AICc value for BM and OU
		if(is.null(opt.list[[nrates]][["bm"]])){ opt.list[[nrates]][["bm"]][["AICc"]] <- 9999}
	}
	optSet.BM <- which.min(sapply(opt.list, function(x) x[["bm"]][["AICc"]]))
	
	opt.listBM <- opt.list[[optSet.BM]][["bm"]]
#	if(opt.list[[optSet.OU]]$ou$AICc <	opt.list[[optSet.BM]]$bm$AICc) optTree[[1]] <- opt.list[[optSet.OU]]$ou
	
#	if(opt.list[[optSet.OU]]$ou$AICc >	opt.list[[optSet.BM]]$bm$AICc) optTree[[1]] <- opt.list[[optSet.BM]]$bm 
	
	return(opt.listBM)
}

matrifyBMsigma <- function(opt.listBM)
{
	optTree <- list()
	
	for(tree in seq(1, length(opt.list), 1)) { #print(opt.list[[nrates]])
		#acquire AICc value for BM and OU
		opt.list[[tree]]
	}
	
	return()
}

#Function to look into extinct/origination across breaks

intervals.Radiant <- function(breakDates, MinTime, RootAge, int.Length = 1, int.Dist = 1){
	breakList <- breakDates
	if(int.Length > 0 && int.Dist > 0){
		for(ii in seq(1, int.Dist, int.Length)){
			intLow <- breakDates - ii
			intHigh <- breakDates + ii
		
			breakList <- append(breakList, intHigh)
			breakList <- append(breakList, intLow)
		}
	}
	breakList <- append(breakList, MinTime)
	breakList <- append(breakList, RootAge)
	breakList <- sort(unique(breakList))
	
	results.Int <- matrix(nrow=length(breakList)-1, ncol=2); colnames(results.Int) <- c("ageTop", "ageBase")
	
	for(jj in seq(1, length(breakList) - 1, 1)){
		results.Int[jj,] <- c(breakList[jj], breakList[jj+1])
	}
	
	rowname.Result <- vector()
	for(kk in seq(1,nrow(results.Int),1)){
		rowname.Result <- append(rowname.Result, paste(paste(paste(results.Int[kk,"ageTop"],"to",sep=""),results.Int[kk,"ageBase"],sep=""),"Ma",sep=""))
	}
	rownames(results.Int) <- rowname.Result
	results.Int <- results.Int[-(length(breakList)-1),]
	
	return(results.Int)
}

#make a function that will return the number of breaks remaining in breakList