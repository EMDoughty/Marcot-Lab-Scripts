source("~/Dropbox/code/R/common_src/strat.R")
# source("https://dl.dropbox.com/s/8jy9de5owxj72p7/strat.R")

makeBLFromNodeDates <- function(tree) {
	tree$edge.length <- tree$node.date[tree$edge[ , 1]] - tree$node.date[tree$edge[ , 2]]
	# tree$edge.length[tree$edge.length < extra] <- extra									extra shouldn't be added here, but in the node dates themselves...
	tree
}

lengthToRoot<-function(tip, tree) {
	lengthToRootMachine<-function(edge, tree, thisLength) {
		thisLength<-thisLength+tree$edge.length[edge]
		if (tree$edge[edge,1]!=length(tree$tip.label)+1) thisLength<-lengthToRootMachine(which(tree$edge[,2]==tree$edge[edge,1]), tree, thisLength)
		return(thisLength)
	}
	return(lengthToRootMachine(which(tree$edge[,2]==tip), tree, 0))
}
makeNodeDatesFromBL <- function(tree, rootAge=NULL, youngestTipAge=0.0) {
	tree$node.date <- array(NA, dim=tree$Nnode)
	if (is.null(rootAge)) tree$node.date[length(tree$tip.label)+1] <- youngestTipAge + max(sapply(1:length(tree$tip.label), lengthToRoot, tree)) else tree$node.date[length(tree$tip.label)+1] <- rootAge
	for (i in 1:nrow(tree$edge)) {
		tree$node.date[tree$edge[i,2]] <- tree$node.date[tree$edge[i,1]] - tree$edge.length[i]
	}
	return(tree)
}


setNodeDatesMachine<-function(tree, nodeID, extra=5e-5) {
	desc <- tree$edge[tree$edge[,1] == nodeID, 2]				#get all desc of this node
	noDates <- desc[which(!is.finite(tree$node.date[desc]))] 	#find desc without dates
	if (length(noDates) > 0) {
		for (i in noDates) {
			if (i > length(tree$tip.label)) tree <- setNodeDatesMachine(tree=tree, nodeID=i, extra) # recursively check those nodes, then get back to me
		}
	}

	if (any(is.finite(tree$node.date[desc]))) {
		tree$node.date[which(!is.finite(tree$node.date[desc]))] <- max(tree$node.date[desc], na.rm=TRUE)	# if any descendants are STILL NA, then assign them the max desc age -> STILL ALLOWS FOR NA DATES IF SISTERS ALL DESC OF A DESC ARE MISSING
		tree$node.date[nodeID] <- max(tree$node.date[desc], na.rm=TRUE) + extra
	}
	tree
}

updateNodeDatesUsingTipNodeDates<-function(tree, extra=5e-5) {
	#requires that tree$node.date for tips are correct and will not be changed
	tree$node.date <- c(tree$node.date[tree$tip.label], array(NA, dim=tree$Nnode))
	setNodeDatesMachine(tree, nodeID=length(tree$tip.label)+1, extra=extra)
}

dropTipKeepDates <- function(tree, tip, extra=0.5) {
	tree <- ape::drop.tip(phy=tree, tip=tip, rooted=FALSE)
	makeNodeDatesFromBL(tree, youngestTipAge=min(tree$node.date))
}

# scaleNodeDates<-function(cladeDuration=1, tree) {
	# tree$node.date<-tree$node.date*(cladeDuration/max(sapply(1:length(tree$tip.label), lengthToRoot, tree)))
	# return(tree)
# }

# scaleBL<-function(cladeDuration=1, tree) {
	# tree$edge.length<-tree$edge.length*(cladeDuration/max(sapply(1:length(tree$tip.label), lengthToRoot, tree)))
	# return(tree)
# }

### function ensures that anc nodes are at least [extra] older than their oldest descendant
rectifyNodeAges <- function(tree, extra=0.0) {
	for (thisNode in length(tree$node.date):(length(tree$tip.label)+1)) { # going from youngest internal node to the root node
		if (max(tree$node.date[tree$edge[tree$edge[,1] == thisNode,2]], na.rm=TRUE) >= tree$node.date[thisNode]) 
			tree$node.date[thisNode] <- max(tree$node.date[tree$edge[tree$edge[,1] == thisNode,2]], na.rm=TRUE) + extra
	}
	makeBLFromNodeDates(tree)
}

extendNodeAges <- function(tree, method=c("uniform", "exponential"), rootAdj=0.0, extra=0.0) {
	# randomly pushes internal nodes back in time, beginning with root and working toward tips 
	method <- match.arg(method)
	ntaxa <- length(tree$tip.label)
	if (rootAdj > 0) {
		thisNode <- ntaxa + 1	# this is ID of the root node
		if (method=="uniform") tree$node.date[thisNode] <- runif(n=1, min=tree$node.date[thisNode], max=(tree$node.date[thisNode] + rootAdj)) #adjust root node age
		# if (method=="exponential") tree$node.date[thisNode]<-rexp(1), min=tree$node.date[thisNode], max=(tree$node.date[thisNode]+rootAdj)) #adjust root node age
		# descBranches<-which(tree$edge[,1]==thisNode)	#find branches desc from root
		# tree$edge.length[descBranches]<-tree$node.date[thisNode]-tree$node.date[tree$edge[descBranches,2]] #lengthens branches leading FROM the root
	}
	for (thisNode in seq(from=(ntaxa + 2), to=length(tree$node.date))) { # root node (ntaxa+1) is fixed, other nodes can become older
		ancNode <- tree$edge[tree$edge[,2]==thisNode, 1]
		if (tree$node.date[ancNode] - tree$node.date[thisNode] > extra) tree$node.date[thisNode] <- runif(1, min=tree$node.date[thisNode], max=(tree$node.date[ancNode] - extra))
		# if (is.na(tree$node.date[thisNode])) print(thisNode)
		# descBranches<-which(tree$edge[,1]==thisNode)
		# tree$edge.length[descBranches]<-tree$node.date[thisNode]-tree$node.date[tree$edge[descBranches,2]] #lengthens branches leading FROM this node
		# tree$edge.length[which(tree$edge[,1]==ancNode & tree$edge[,2]==thisNode)]<-tree$node.date[ancNode]-tree$node.date[thisNode] # shortens branch leading TO this node
	}
	makeBLFromNodeDates(tree)
}

backfillMachine <- function(tree, nodeID, extra=0.0, min.date=0.0) {
	if (is.finite(tree$node.date[nodeID])) return(tree) 
	ancNode <- tree$edge[tree$edge[,2]==nodeID,1]
	if (!is.finite(tree$node.date[ancNode])) {
		if (ancNode == length(tree$tip.label) + 1) { print("backfillMachine error:	root node date is missing")
		} else tree <- backfillMachine(tree, ancNode, extra, min.date)
	} 
	
	if ((tree$node.date[ancNode] - extra) < min.date) tree$node.date[nodeID] <- min.date else tree$node.date[nodeID] <- (tree$node.date[ancNode] - extra)
	tree
}

backfillMissingNodeDates <- function(tree, extra=0.0, min.date=0.0) {
	for (i in sort(which(!is.finite(tree$node.date)), decreasing=FALSE)) tree <- backfillMachine(tree=tree, nodeID=i, extra=extra, min.date=min.date)
	tree
}

dateTreeWithRanges <- function(tree, dates, within.error=FALSE, random=TRUE, min.date=0.0, extra=0.0, flex.internals=FALSE, rootAdj=0.0, keep.missing=FALSE, missing.date=NULL, include.ranges=FALSE) {
	# tree <- thisTree
	# dates <- ranges
	dates <- dates[rownames(dates) %in% tree$tip.label,] 		#cull taxa in dates list that aren't in the tree
	dates <- dates[match(tree$tip.label, rownames(dates)),]		#reorder taxa in dates by the order in the tree
	if (keep.missing) {
		rownames(dates) <- tree$tip.label	# the match above erases unmatched species from the rownames - this restablishes them
		if (!is.null(missing.date)) dates[apply(dates, 1, function(x) all(!is.finite(x))),] <- missing.date
	} else {
		if (any(!(tree$tip.label %in% rownames(dates)))) {
			tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% rownames(dates))]) # cut taxa that are not in ranges
			dates <- dates[match(tree$tip.label, rownames(dates)),]		#reorder taxa in dates by the order in the tree
		}
	}
	
	# set node dates of tips
	if (random) tree$node.date <- getRandomRanges(dates)[,"FO"] else tree$node.date <- rowMeans(dates[,1:2]) #if without error, use the midpoint of FA dates
	# make space for internal nodes
	tree$node.date <- c(tree$node.date, array(data=NA, dim=tree$Nnode)) #add NA dates for internal nodes
	tree <- setNodeDatesMachine(tree=tree, nodeID=length(tree$tip.label)+1, extra=extra)
	tree <- backfillMissingNodeDates(tree, extra)
	if (flex.internals) tree <- extendNodeAges(tree, method="uniform", rootAdj=rootAdj, extra=extra)
	tree <- makeBLFromNodeDates(tree)
	if (include.ranges) {
		if (within.error) { tree$strat.range <- cbind(tree$node.date[seq_along(tree$tip.label)], apply(dates[,c("ageLO_max", "ageLO_min")], 1, function(x) runif(1, min=x[2], max=x[1]) ))
		} else tree$strat.range <- cbind(tree$node.date[seq_along(tree$tip.label)], dates[,"LO"])
		tree$strat.range[tree$strat.range[,1]<tree$strat.range[,2],] <- tree$strat.range[tree$strat.range[,1]<tree$strat.range[,2],c(1,1)]
	}
	tree
}

getProportionsOfBranchesWithinIntervals<-function(tree, intervals) {
	oneBranchOneInterval<-function(thisInterval, thisBranch) {
		if (thisBranch[1] > thisInterval["ageBase"] & thisBranch[2] < thisInterval["ageTop"]) {														return ((thisInterval["ageBase"] - thisInterval["ageTop"]) / (thisBranch[1]-thisBranch[2])) # bt
		} else if (thisBranch[1] <= thisInterval["ageBase"] & thisBranch[2] >= thisInterval["ageTop"]) {											return (1.0)  #FL
		} else if (thisBranch[1] <= thisInterval["ageBase"] & thisBranch[1] > thisInterval["ageTop"] & thisBranch[2] < thisInterval["ageTop"]) {	return ((thisBranch[1] - thisInterval["ageTop"]) / (thisBranch[1] - thisBranch[2]))	#Ft
		} else if (thisBranch[1] > thisInterval["ageBase"] & thisBranch[2] < thisInterval["ageBase"] & thisBranch[2] >= thisInterval["ageTop"]) {	return ((thisInterval["ageBase"]-thisBranch[2]) / (thisBranch[1] - thisBranch[2]))	#bL
		} else return(0.0)
	}
	
	# ancAge=tree$node.date[tree$edge[i,1]]
	# descAge=tree$node.date[tree$edge[i,2]]
	intvList<-listifyMatrixByRow(intervals)
	t(apply(tree$edge, 1, function(x) sapply(intvList, oneBranchOneInterval, c(tree$node.date[x[1]], tree$node.date[x[2]]))))
}

getLMAForSingleInterval<-function(thisInt, thisBranch) {
	if (thisBranch["ageDesc"]>=thisInt["ageBase"] | thisBranch["ageAnc"]<=thisInt[ "ageTop"]) { return(0)
	} else if (thisBranch["ageAnc"]>=thisInt["ageBase"] & thisBranch["ageDesc"]<=thisInt["ageTop"]) { return(thisInt["ageBase"] - thisInt["ageTop"])		#bt
	} else if (thisBranch["ageAnc"]<=thisInt["ageBase"] & thisBranch["ageDesc"]>=thisInt["ageTop"]) { return(thisBranch["ageAnc"]-thisBranch["ageDesc"] )	#FL
	} else if (thisBranch["ageAnc"]<=thisInt["ageBase"] & thisBranch["ageDesc"]<=thisInt["ageTop"]) { return(thisBranch["ageAnc"] - thisInt["ageTop"])		#Ft
	} else if (thisBranch["ageAnc"]>=thisInt["ageBase"] & thisBranch["ageDesc"]>=thisInt["ageTop"]) { return(thisInt["ageBase"] - thisBranch["ageDesc"]) }	#bL
}

getLMAForSingleBranch<-function(thisBranch, intervals) {
	apply(intervals, 1, getLMAForSingleInterval, thisBranch)
}

getBranchAges<-function(thisBranch, tree, dates, use.ranges=FALSE) { 
	if (use.ranges) {
		if (thisBranch[2]<=length(tree$tip.label)) return(c(ageAnc=tree$node.date[thisBranch[1]], ageDesc=unname(tree$strat.range[thisBranch[2],2])) )
	} #else 
	c(ageAnc=tree$node.date[thisBranch[1]], ageDesc=tree$node.date[thisBranch[2]])
}

