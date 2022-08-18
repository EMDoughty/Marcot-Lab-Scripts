tree_resolution_v3 <- function(tree.backbone, clade.definitions, wildcard.positions, tree.list.wildcards, occs, thisMat, min.bl=1.0) {
	tree_resolv <- bindWildCards(tree.backbone, clade.definitions, wildcard.positions, tree.list.wildcards)
	#drop duplicate tips from tree
	tree_resolv <- drop_DupTips(tree_resolv)
	print("Dropping duplicate tips is complete")
	
	#drop tips from tree that do no coincide with entry in species data
	tree_resolv <- drop.tip(tree_resolv, tip=tree_resolv$tip.label[!tree_resolv$tip.label %in% thisMat$species])
	tree_resolv <- date_treeJon(this.tree = tree_resolv, occs = occs, root.max = 85, min.bl=min.bl)
	
	print("Tree Dating is complete")
	
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
checkRuminantiamorpha <- function(tree, wildcardTree, clade.definitions) {
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
addArbTip2Root <- function(tree.wildcard) {
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


iterations <- function() { 
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

date_treeJon<-function(this.tree, occs, root.max, min.bl=1.0) {
	ranges <- getTaxonRangesFromOccs(occs, random =TRUE)
	rownames(ranges) <- gsub(pattern = "[[:space:]]",replacement ="_", x = rownames(ranges))
	this.tree <- dateTreeWithRanges(this.tree = this.tree, strat.ranges=ranges, within.error=FALSE, min.bl=min.bl, smooth.node.ages=TRUE, root.max=root.max) # internals will be flexed by extendNodeAges below

	return(this.tree)
}

############################################################################################################################################

drop_DupTips <- function(tree) {
	# tree$tip.label <- gsub(pattern = "[[:space:]]", replacement = "_", x = getCurrentTaxa(gsub("_", " ", tree$tip.label, fixed=TRUE)))
	dup.taxa <- unique(tree$tip.label[duplicated(tree$tip.label)])
	for (i in seq_along(dup.taxa)) tree <- drop.tip(tree, tip=sample(which(tree$tip.label==dup.taxa[i]), size=sum(tree$tip.label==dup.taxa[i])-1))
	tree
}

############################################################################################################################################

compile_SpecData_Old <- function(occs){
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


#########################################################################################################

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
comp_TreeSpecMat_outTree <- function(this.tree, thisMat) {
	thisMat <- thisMat[is.finite(thisMat$bodyMass),] 
	this.tree <- drop.tip(phy=this.tree, tip=this.tree$tip.label[!this.tree$tip.label %in% rownames(thisMat)])
	return(this.tree)
}

comp_TreeSpecMat_outMatrix <- function(this.tree,thisMat) {
	thisMat <- thisMat[is.finite(thisMat$bodyMass),]
	thisMat <- thisMat[rownames(thisMat) %in% this.tree$tip.label,]
	#thisMat <- thisMat[this.tree$tip.label,]
	# cbind(this.tree$tip.label, rownames(thisMat))
	return(thisMat)	
}

############################################################################################################################################
#functions to painttree in various ways; mostly used in troubleshooting and checking that tree are constructed properly
paintWildcardClades <- function(tree) {

	this.map <-  tree #artio.tree
	
	#plotTree(this.map,node.numbers=T,pts=F,fsize=0.1, ftype = "i")
	
	#cols<-c("black","blue","red","green","cyan","magenta","yellow","gray","darkgoldenrod","hotpink","midnightblue","sienna","skyblue","seagreen","salmon","maroon","khaki","cadetblue","burlywood","chocolate","aquamarine","black","blue","red","green","cyan","magenta","yellow","gray","darkgoldenrod","hotpink","midnightblue","sienna","skyblue","seagreen","salmon","maroon","khaki","cadetblue","burlywood","chocolate","aquamarine"); names(cols)<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42)
	cols <- c("black", "cyan","magenta", "blue", "seagreen", "chocolate", "darkgray", "green", "red"); names(cols) <- c("1","5","6","7","14","20","29","30","33")
	
	#Achaenodontidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Parahyus_vagus", sp2="Achaenodon_fremdi"),state="5")
	#plotSimmap(this.map, cols,fsize=0.1, ftype = "i")
	
	#Agriochoeridae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Diplobunops_crassus", sp2="Protoreodon_transmontanus"),state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Diplobunops_crassus", sp2="Agriochoerus_antiquus"),state="6")
	#plotSimmap(this.map,cols,fsize=0.1, ftype = "i")
	
	#Merycoidontidae
	#polytomy within the Merycoidodontidae leads to multiple different configurations that require multiple different configs fpr painting
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Bathygenys_reevesi", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Aclistomycter_middletoni", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Miniochoerus_affinis", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Limnenetes_platyceps", sp2="Merychyus_minimus"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Aclistomycter_middletoni", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Miniochoerus_affinis", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Limnenetes_platyceps", sp2="Bathygenys_reevesi"), state="6")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Limnenetes_platyceps", sp2="Miniochoerus_affinis"), state="6")
	#plotSimmap(this.map, cols, fsize=0.1, ftype = "i")
	
	#Anthracotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Heptacodon_yeguaensis", sp2="Aepinacodon_deflectus"),state="7")
	#plotSimmap(this.map, cols,fsize=0.1, ftype = "i")
	
	#Entelodontidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Archaeotherium_palustris", sp2="Dinohyus_hollandi"), state="14")
	#this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Choerodon_calkinsi", sp2="Archaeotherium_palustris"), state="14")
	#	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Archaeotherium_minimus", sp2="Choerodon_caninus"), state="14")
	#	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Archaeotherium_minimus", sp2="Brachyhyops_viensis"), state="14")
	#	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Archaeotherium_minimus", sp2="Boochoerus_angustus"), state="14")
	#	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Archaeotherium_minimus", sp2="Megachoerus_latidens"), state="14")
	#plotSimmap(this.map, cols, fsize=0.1, ftype = "i")
	
	#Protoceratidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Leptotragulus_ultimus", sp2="Prosynthetoceras_orthrionanus"), state="20")
	#plotSimmap(this.map, cols,fsize=0.1, ftype = "i")
	
	#Brontotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Menops_marshi", sp2="Parvicornus_occidentalis"), state="29")
	#this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Menops_bakeri", sp2="Eotitanops_minimus"), state="29")
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Megacerops_osborni", sp2="Brontops_dispar"), state="29")
	#plotSimmap(this.map, cols, fsize=0.1, ftype = "i")
	
	#Chalicotheriidae
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Moropus_merriami", sp2="Eomoropus_anarsius"), state="30")
	#plotSimmap(this.map,cols, fsize=0.1, ftype = "i")
	
	#Heptodon
	this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Heptodon_posticus", sp2="Heptodon_calciculus"), state="33")
	plotSimmap(this.map, cols, fsize=0.1, ftype = "i")
	
	return(this.map)
}

paintOrders <- function(tree) {
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
	
	#this.map <- paintSubTree(tree=this.map, node=findMRCA(this.map, tips = c("Achaenodon_robustus","Antiacodon_venustus")),anc.state="1", state="2")
	
	#Perissodactyla
	#this.map <- paintSubTree(tree=this.map, node=fastMRCA(tree=this.map, sp1="Heptodon_posticus", sp2="Notiotitanops_mississippiensis"), state="3")
	
	plotSimmap(this.map, cols, fsize=0.1, ftype = "i", lwd = 1)
	# nodelabels(cex = 0.1)
	
	return(this.map)
}

paintTemporalRange <- function(tree, maxAge, minAge, int.length) {
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

tree_TestDate<- function(tree.backbone, clade.definitions, wildcard.positions, tree.list.wildcards, occs, thisMat, reps = 1, min.bl=1.0) {
	tree_list <- list()
	tree_resolv <- bindWildCards(tree.backbone, clade.definitions, wildcard.positions, tree.list.wildcards)
	#drop duplicate tips from tree
	tree_resolv <- drop_DupTips(tree_resolv)
	print("Dropping duplicate tips is complete")
	
	#drop tips from tree that do no coincide with entry in species data
	tree_resolv <- drop.tip(tree_resolv, tip=tree_resolv$tip.label[!tree_resolv$tip.label %in% thisMat$species])
	for(ii in seq_len(reps)){ 
		tree_list[[ii]] <- date_treeJon(this.tree = tree_resolv, occs = occs, root.max = 85, min.bl=min.bl)
	}	
	print("Tree Dating is complete")
	
	return(tree_list)
}

checkWildcards <- function(tree.list, check.wildcards = FALSE)
{
 if(check.wildcards){
 	par(mfrow=c(5,2))
 	for(ii in 1:length(tree.list)) {
 		paintWildcardClades(tree.list[[ii]])
 		}
 }
}
