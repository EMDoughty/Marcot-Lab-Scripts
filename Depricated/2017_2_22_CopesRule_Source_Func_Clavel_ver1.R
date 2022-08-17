############################################################################################################################################

getRandomMatrix <- function(n.row, lower.bound=-10e6, upper.bound=10e6, this.constraint=FALSE) {
	this <- matrix(NA, nrow=n.row, ncol=n.row)
	this[upper.tri(this)] <- runif(n=length(this[upper.tri(this)]), min=lower.bound, max=upper.bound)
	this[lower.tri(this)] <- t(this)[lower.tri(t(this))]
	diag(this) <- runif(n=n.row, min=lower.bound, max=upper.bound)
	chol_factor <- chol(this)
	if(this.constraint==FALSE) chol_factor[upper.tri(chol_factor, diag=TRUE)] else chol_factor[upper.tri(chol_factor, diag=FALSE)]
	# this[upper.tri(this, diag=TRUE)]
}

############################################################################################################################################

doOneMultipleModel <- function(this.map, dat, intBreaks, break.dates, this.model= list("BMM", "OUM"), full.output=FALSE, analysis.settings) {
	flag <- FALSE
	this.param.BM <- list(constraint = analysis.settings$this.constraint, smean = TRUE, trend = FALSE)    #default values to begin analysis
	this.param.OU <- list(sigma = NULL, alpha = NULL, vcv = "fixedRoot", root=TRUE, decomp = c("cholesky", "spherical", "eigen", "qr", "diagonal", "upper", "lower")) 
	
	counter <- 0
	# while(!flag) {
		counter <- counter + 1 # counts number of attempts
		if (counter %% analysis.settings$adjust.date.after == 0) {
			print("***** Adjusting break dates *****\n")
			print(counter)
			break.dates <- break.dates + analysis.settings$adjust.date.increment
			this.limits <- max(this.map$node.date) - c(max(this.map$node.date), break.dates)
			this.map <- make.era.map(this.map, limits=this.limits)
		}
		result <- NULL
		if (this.model == "BMM") {
			########### mvMORPH
			try(result <- mvBM(tree=this.map,
												 data=dat,
												 model= "BMM",
												 param=this.param.BM,
												 # method="rpf",
												 control=list(maxit=1e8),
												 # optimization="subplex",
												 diagnostic=FALSE,
												 echo=FALSE))
			############ OUwie
			# try(result <- OUwie(phy=this.map, 
			# 										data=dat, 
			# 										model="BMS", 
			# 										simmap.tree = TRUE, 
			# 										root.age=max(this.map$node.date)))
			result$model <- "BMM"
		} else {
			########### mvMORPH
			try(result <- mvOU(tree=this.map,
												 data=dat,
												 model= "OUM",
												 param=this.param.OU,
												 # method="rpf",
												 control=list(maxit=1e8),
												 # optimization="subplex",
												 # smean=FALSE, 
												 diagnostic=FALSE,
												 echo=FALSE))
			############ OUwie
			# try(result <- OUwie(phy=this.map, 
			# 										data=dat, 
			# 										model="OUMVA", 
			# 										simmap.tree = TRUE, 
			# 										root.age=max(this.map$node.date)))
			result$model <- "OUM"
		}
		if (is.null(result) | !(all(class(result) %in% c("mvmorph", "mvmorph.bm", "mvmorph.ou")))) {	### class(result)=="try-error"
			print(result)
			print("***** NULL Result: Adjusting starting parameter values *****\n")
			if(this.model == "BMM"){
				print("BMM")
				this.param.BM$sigma <- runif(1)
				# 
				# if (analysis.settings$this.constraint==FALSE) {  this.param.BM$sigma <- replicate(n=length(intBreaks)+1,expr= getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint), simplify=FALSE)
				# } else this.param.BM$sigma <- replicate(n=length(intBreaks)+1,expr= c(getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint), runif(1)), simplify=FALSE)
			} else {
				print("OUM")
				# this.param.OU$sigma <- getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint)
				# this.param.OU$alpha <- getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint)
				dummy <- mvOU(tree=this.map, data=dat, model="OU1", echo=FALSE)
				this.param.OU$sigma <- dummy$sigma
				this.param.OU$alpha <- dummy$alpha
				rm(dummy)
			}
		} else if (result$convergence != 0) {
			print(result)
			print("***** Failure to converge : Adjusting starting parameter values *****\n")
			cat("result$convergence", result$convergence, "\n")
			if(this.model == "BMM") {
				print("BMM")
				this.param.BM$sigma <- replicate(n=length(intBreaks)+1, expr=runif(1))
				
				# if (analysis.settings$this.constraint==FALSE) {
				# 	this.param.BM$sigma <- replicate(n=length(intBreaks)+1,expr= getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint), simplify=FALSE)
				# } else this.param.BM$sigma <- replicate(n=length(intBreaks) +1, expr= c(getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint)), simplify=FALSE)
			} else {
				print("OUM")
				if (analysis.settings$this.constraint==FALSE) {
					# this.param.OU$sigma <- getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint)
					# this.param.OU$alpha <- getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint)
					
					dummy <- mvOU(tree=this.map, data=dat, model="OU1", echo=FALSE)
					this.param.OU$sigma <- dummy$sigma
					this.param.OU$alpha <- dummy$alpha
					rm(dummy)
				} else {
					# this.param.OU$sigma <- getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint)
					# this.param.OU$alpha <- getRandomMatrix(n.row=ncol(dat), lower.bound=0, upper.bound=10, this.constraint=analysis.settings$this.constraint)
					dummy <- mvOU(tree=this.map, data=dat, model="OU1", echo=FALSE)
					this.param.OU$sigma <- dummy$sigma
					this.param.OU$alpha <- dummy$alpha
					rm(dummy)
				}
			}
			# } else if (result$hess.values != 0) {
			# 	print("***** Failure to HESS: Adjusting starting parameter values *****\n")
			# 	if (this.constraint==FALSE) { this.param$sigma <- replicate(n=length(intBreaks)+1,expr=as.double(runif(ncol(dat)*(ncol(dat)+1)/2)), simplify=FALSE)
			# 	} else this.param$sigma <- replicate(n=length(intBreaks)+1,expr=as.double(runif(n=(ncol(dat)*(ncol(dat)-1)/2)+1)), simplify=FALSE)
		} else flag=TRUE
	
	if (full.output) output <- result else output <- result$AICc #, "theta", "alpha", "sigma")]
	return(output)
}

############################################################################################################################################

#function performs BM and/or OU model on given set of temporal breaks.  Failure of model triggers a reassignemnts of the sigma values
getLkMVMorphOneIntervalSet <- function(newBreak=NULL, oldIntBreaks=NULL, full.output=FALSE, this.tree, dat, this.intervals, analysis.settings) {
	intBreaks <- sort(c(newBreak, oldIntBreaks), decreasing=TRUE) #interval numbers
	cat("Breaks:", intBreaks, "\n")
	break.dates <- as.numeric(this.intervals[intBreaks, "ageBase"]) #actual dates of interval/break
	this.limits <- max(this.tree$node.date) - c(max(this.tree$node.date), break.dates) # puts in time from root
	this.map <- make.era.map(this.tree, limits=this.limits)
	
	### code below erases very short terminal segments (min.terminal in duration) 
	### within regimes which causes likelihoods to blow up. I don't know what 
	### the optimal value for min.terminal should be
	if (analysis.settings$trim.short.terminals) {
		for (this.regime in seq(from=2, to=dim(this.map$mapped.edge)[2])) {
			short.terminals <- which((this.map$edge[ ,2] <= length(this.map$tip.label)) & 
									 (this.map$mapped.edge[ ,this.regime] > 0) & 
									 (this.map$mapped.edge[ ,this.regime] < analysis.settings$min.terminal))
			if (length(short.terminals) > 0) {
				# cat("\tTrimming", length(short.terminals), "terminal branches shorter than", min.terminal, "from regime #", this.regime, ". BL:", round(this.map$edge.length[short.terminals], 2), "\n")
				this.map$mapped.edge[short.terminals, this.regime] <- 0
				this.map$maps[short.terminals] <- lapply(this.map$maps[short.terminals], function(x) x[seq_len(length(x)-1)])
			}
		}
		#plotSimmap(this.map,fsize=0.05)
	}
	
	# doOneMultipleModel through all requested models
	if (analysis.settings$do.BMM) { 
		master.BM <- doOneMultipleModel(this.map=this.map, dat=dat, intBreaks=intBreaks, break.dates=break.dates, this.model="BMM", full.output=full.output, analysis.settings = analysis.settings)
	} else master.BM <- NULL
	
	if(analysis.settings$do.OUM) { 
		master.OU <- doOneMultipleModel(this.map=this.map, dat=dat, intBreaks=intBreaks, break.dates=break.dates, this.model="OUM", full.output=full.output, analysis.settings = analysis.settings)
	} else master.OU <-  NULL 
	
	master.result <- list(BM=master.BM, OU=master.OU)
	
	return(master.result)
}

testRateShiftsMVMorphIntervals <- function(tree, dat, intervals, full.output=FALSE, do.plot.phylo.BM=FALSE, analysis.settings) {
	if (is.vector(dat)) {  #if dat is vector (1 trait)
		this.tree <- dropTipKeepDates(tree[[1]], tree$tip.label[!tree$tip.label %in% names(dat)], min.bl=1.0)
		dat <- dat[names(dat) %in% this.tree$tip.label]
		dat <- dat[match(this.tree$tip.label, names(dat))]   #reorders dat to match the taxon order of tree - MUST stay in this function, because each tree will have different taxon order...
	} else {
		this.tree <- dropTipKeepDates(tree, tree$tip.label[!tree$tip.label %in% rownames(dat)], min.bl=1.0)
		dat <- matrix(dat[this.tree$tip.label,], ncol=1, dimnames=list(rownames(dat)[rownames(dat) %in% this.tree$tip.label], colnames(dat))) # this *should* replace the following code
	}
	
	topIntv <- min(which(intervals[, "ageBase"] > min(this.tree$node.date)))
	baseIntv <- max(which(intervals[, "ageTop"] < max(this.tree$node.date[seq_along(this.tree$tip.label)])))
	this.intervals <- intervals[topIntv:baseIntv, ]
	n.intv <- nrow(this.intervals)
	
	optList <- list()

	######## get single rate BM
	x <- NULL
	if (analysis.settings$do.BMsimple) {
		# print("Get Single Rate BM")
		param.list <- list(constraint = analysis.settings$this.constraint)
		############ mvMORPH
		conv <- FALSE
		while (!conv) {
			x <- mvBM(tree = this.tree, data = dat, model = "BM1", param = param.list, echo = FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
			if (x$convergence != 0) {
				print("***** Failure to converge: Adjusting starting parameter values *****\n")
				print(x$convergence)
				# param.list$sigma <- getRandomMatrix(n.row = ncol(dat), this.constraint = analysis.settings$this.constraint)
				param.list$sigma <- runif(1)
				print(param.list$sigma)
			} else conv <- TRUE
		}
		############ OUwie
		# x <- result <- OUwie(phy=this.tree, data=dat, model="BM1", simmap.tree = FALSE, root.age=max(this.tree$node.date))
	}
	
	######## get single rate OU
	conv <- FALSE
	y <- NULL
	if (analysis.settings$do.OUsimple) {
		# print("Get Single Rate OU")
		this.param <- list(sigma = NULL, alpha = NULL, vcv = "fixedRoot", decomp = c("cholesky", "spherical", "eigen", "qr", "diagonal", "upper", "lower"))
		while (!conv) {
			y <- mvOU(tree = this.tree, data = dat, model = "OU1", param = param.list, echo = FALSE) # set 0 shift or 1 model baseline;  DOES NOT REACH CONVERGENCE OF OPTIMIZER
			if (y$convergence != 0) {
				print("***** Failure to converge: Adjusting starting parameter values *****\n")
				print(y$convergence)
				# param.list$sigma <- getRandomMatrix(n.row = ncol(dat), this.constraint = analysis.settings$this.constraint)
				param.list$sigma <- runif(1)
				print(param.list$sigma)
			} else conv <- TRUE
		}
		############ OUwie
		# y <- OUwie(phy=this.tree, data=dat, model="OU1", simmap.tree = FALSE, root.age=max(this.tree$node.date))
	}
	optList[[1]] <- list(BM=x, OU=y)
	rm(x, y)
	
	######## start multiple rate analyses    
	
	flag <- FALSE #get thrown when AIC starts to increase; must be fewer rates than intervals
	nrates <- 2
	oldIntBreaks <- NULL

	# print("Start Multiple Rate Loop")
	while (!flag & nrates <= n.intv) {
		cat("Beginning analysis for ", nrates, " rates...\r")
		### NOTE: we use "n.intv - 1" because shift is taken to bbe at the base of an interval.
		### and the oldest interval's base has no other intervals bebfore it, so can't be a shift
		### However, there can be "tree" older than the oldest interval, 
		### a rate shift at the babse of the oldest inteerval is theorestically possible
		if (analysis.settings$do.heuristic) { # builds break list; possible shift points; adds on to last nrates and optimum
			if (nrates > 2) { 
				if (analysis.settings$do.BMM) breakList.BM <- as.integer(seq_len(n.intv - 1)[-optList[[nrates-1]]$BM$breaks]) else breakList.BM <- NULL
				if (analysis.settings$do.OUM) breakList.OU <- as.integer(seq_len(n.intv - 1)[-optList[[nrates-1]]$OU$breaks]) else breakList.OU <- NULL
			} else {
				if (analysis.settings$do.BMM) breakList.BM <- as.integer(seq_len(n.intv - 1)) else breakList.BM <- NULL# breaks are at the base of intervals
				if (analysis.settings$do.OUM) breakList.OU <- as.integer(seq_len(n.intv - 1)) else breakList.OU <- NULL# breaks are at the base of intervals
			}
			# breakList <- sample(breakList, size = length(breakList)) #randomly reorders breakList to avoid solidifying the \early\" breaks"
		} else breakList <- listifyMatrixByColumn(m = combn(x = n.intv - 1, m = nrates - 1)) # all possible combinations of shift points
		
		if (analysis.settings$do.heuristic) {
			original.settings.BMM <- analysis.settings$do.BMM
			original.settings.OUM <- analysis.settings$do.OUM

			if (original.settings.BMM) {
				analysis.settings$do.OUM <- FALSE
				if (analysis.settings$do.parallel) {
					rez.list.BM <- mclapply(breakList.BM, getLkMVMorphOneIntervalSet, oldIntBreaks=oldIntBreaks$BM, full.output=FALSE, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings=analysis.settings, mc.cores = max(c(detectCores()-2, 3)), mc.preschedule = FALSE)
				} else rez.list.BM <- lapply(breakList.BM, getLkMVMorphOneIntervalSet, oldIntBreaks=oldIntBreaks$BM, full.output=FALSE, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings=analysis.settings)
				analysis.settings$do.OUM <- original.settings.OUM
			}

			if (original.settings.OUM) {
				analysis.settings$do.BMM <- FALSE
				if (analysis.settings$do.parallel) {
					rez.list.OU <- mclapply(breakList.OU, getLkMVMorphOneIntervalSet, oldIntBreaks=oldIntBreaks$OU, full.output=FALSE, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings=analysis.settings, mc.cores = max(c(detectCores()-2, 3)), mc.preschedule = FALSE)
				} else rez.list.OU <- lapply(breakList.OU, getLkMVMorphOneIntervalSet, oldIntBreaks=oldIntBreaks$OU, full.output=FALSE, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings=analysis.settings)
				analysis.settings$do.BMM <- original.settings.BMM
			}
			
			rez.list <- list()
			if (original.settings.BMM & original.settings.OUM) {
				for (i in seq_along(rez.list.BM)) rez.list[[i]] <- list(BM=rez.list.BM[[i]]$BM, OU=rez.list.OU[[i]]$OU)
			} else if (original.settings.BMM) {
				rez.list <- rez.list.BM
			} else rez.list <- rez.list.OU
			
		} else {	#do.heuristic is FALSE = do exhaustive search
			if (analysis.settings$do.parallel) {
				rez.list <- mclapply(breakList, getLkMVMorphOneIntervalSet, oldIntBreaks, full.output=FALSE, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings=analysis.settings, mc.cores = max(c(detectCores()-2, 3)), mc.preschedule = FALSE)
			} else rez.list <- lapply(breakList, getLkMVMorphOneIntervalSet, oldIntBreaks, full.output=FALSE, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings=analysis.settings)
		}
		
		# cbind(BM=sapply(rez.list, function(x) x$BM), OU=sapply(rez.list, function(x) x$OU))
		
		# if the BM AICc of this number of rates is not lower than the previous numnber of rates, then stop checking BM in higher numbers of rates
		if (analysis.settings$do.BMM) {
			#print(rez.list)
			#cat(min(sapply(rez.list, function(x) x$BM)),"\n")
			if (min(sapply(rez.list, function(x) x$BM)) > optList[[nrates - 1]][["BM"]]$AICc) analysis.settings$do.BMM <- FALSE
		}
		# if the OU AICc of this number of rates is not lower than the previous numnber of rates, then stop checking OU in higher numbers of rates
		if (analysis.settings$do.OUM) {
			if (min(sapply(rez.list, function(x) x$OU)) > optList[[nrates - 1]][["OU"]]$AICc) analysis.settings$do.OUM <- FALSE
		}
		# if we stopped checking both models, abort the while loop
		if (!analysis.settings$do.BMM & !analysis.settings$do.OUM) flag <- TRUE
		
		if (!flag) {
			if (analysis.settings$do.BMM) {
				optSet.BM <- which.min(sapply(rez.list, function(x) x$BM))
				this.breaks <- sort(c(oldIntBreaks$BM, breakList.BM[[optSet.BM]]), decreasing = TRUE)
				optList.BM <- getLkMVMorphOneIntervalSet(newBreak=this.breaks, oldIntBreaks=NULL, full.output=TRUE, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings=analysis.settings)$BM
				optList.BM$breaks <- this.breaks
				optList.BM$break.dates <- this.intervals[optList.BM$breaks, "ageBase"]
				optList.BM$model <- "BM"
				if (analysis.settings$do.heuristic) oldIntBreaks$BM <- as.integer(optList.BM$breaks)
			} else optList.BM <- NULL
			
			if (analysis.settings$do.OUM) {
				optSet.OU <- which.min(sapply(rez.list, function(x) x$OU))
				this.breaks <- sort(c(oldIntBreaks$OU, breakList.OU[[optSet.OU]]), decreasing = TRUE)
				optList.OU <- getLkMVMorphOneIntervalSet(newBreak=breakList.OU[[optSet.OU]], oldIntBreaks=oldIntBreaks$OU, full.output=TRUE, this.tree=this.tree, dat=dat, this.intervals=this.intervals, analysis.settings=analysis.settings)$OU
				optList.OU$breaks <- this.breaks
				optList.OU$break.dates <- this.intervals[optList.OU$breaks, "ageBase"]
				optList.OU$model <- "OU"
				if (analysis.settings$do.heuristic) oldIntBreaks$OU <- as.integer(optList.OU$breaks)
			} else optList.OU <- NULL
			optList[[nrates]] <- list(BM=optList.BM, OU=optList.OU)
			if (do.plot.phylo.BM) {
				plot.phylo(ladderize(this.tree, right=FALSE), direction="upwards", no.margin=TRUE, cex=0.3)
				plotRateShiftsBM(this.rez=optList[[nrates]], this.tree = this.tree)
			}
		}
		
		nrates <- nrates + 1
	}
	# print(getOptModelStatisticsFromOptList (optList))
	return(optList)
}

getOptModelStatisticsFromOptList <- function(optList) {
	max.rates.BM <- max(which(!sapply(lapply(optList, function(x) x$BM), is.null))) 
	max.rates.OU <- max(which(!sapply(lapply(optList, function(x) x$OU), is.null)))
	BM <- data.frame(breaks=c(rev(optList[[max.rates.BM]]$BM$breaks), "-"), 
						 break.dates=c(paste(rev(optList[[max.rates.BM]]$BM$break.dates), "Ma", ""), "-"),
						 sigma=round(rev(as.vector(optList[[max.rates.BM]]$BM$sigma)), 3))
	
	if (max.rates.OU>1) {
		OU <- data.frame(breaks=c(rev(optList[[max.rates.OU]]$OU$breaks), "-", "Anc"), 
						 break.dates=c(paste(rev(optList[[max.rates.OU]]$OU$break.dates), "Ma", ""), "-", "-"),
						 sigma=round(rev(as.vector(optList[[max.rates.OU]]$OU$sigma)), 3), 
						 alpha=round(rev(as.vector(optList[[max.rates.OU]]$OU$alpha)), 3), 
						 theta=rev(round(as.vector(optList[[max.rates.OU]]$OU$theta), 3)))
	} else OU <- data.frame(breaks="_", break.dates="_", sigma=optList[[max.rates.OU]]$OU$sigma, alpha=optList[[max.rates.OU]]$OU$alpha, theta=optList[[max.rates.OU]]$OU$theta)
	list(BM=BM, OU=OU)
}

getTreeOptList <- function(tree.list, dat, intervals, analysis.settings = list(do.BMsimple = FALSE, do.BMM = FALSE, do.OUsimple = FALSE, do.OUM = FALSE, do.heuristic = FALSE, adjust.date.after=10, adjust.date.increment = 0.1, do.parallel = FALSE)) {
	if (class(tree.list)=="phylo") { results <- testRateShiftsMVMorphIntervals(tree=tree.list, dat = dat, intervals = intervals, analysis.settings = analysis.settings) 
	} else results <- lapply(tree.list, FUN=testRateShiftsMVMorphIntervals, dat = dat, intervals = intervals, analysis.settings = analysis.settings)
	# for (i in seq_along(tree.list)) results[[i]] <- testRateShiftsMVMorphIntervals (tree.list[[i]], dat = dat, intervals = intervals, analysis.settings = analysis.settings)

	return(results)
}

checkBreaks <- function(breakList, tree, occs, nrates, intervals) {
	####### get species within intervals
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
				# param.list$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint= analysis.settings$this.constraint)
				param.list$sigma <- runif(1)
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
				# param.list$sigma <- getRandomMatrix(n.row=ncol(dat), this.constraint= analysis.settings$this.constraint)
				param.list$sigma <- runif(1)
				print(param.list$sigma)
			} else conv <-TRUE
		}
		print("Intitial OU model complete")
		y$model <- "OU"
	}    
	optList[[1]] <- list(BM=x, OU=y)
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
	print("Get Complex Rate OU") 
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
	optList[[2]] <- list(BM=x, OU=y)
	rm(x, y)
	
	return(optList)
}

#function to go through list of optList to extract AICc, nrates, location of shifts, and parameters
getBestOptModels <- function(opt.list, do.plot = FALSE, folder = NULL) #function will save/check rez.list to find what rate configurations have been completed in order to continue in the event of a shutdown
{
	optTree <- list()
	
	#read in optList list or optList files
	
	#isolate optList for single tree
	for(nrates in seq(1, length(opt.list), 1)) { #print(opt.list[[nrates]])
		#acquire AICc value for BM and OU
		if(is.null(opt.list[[nrates]][["BM"]])){ opt.list[[nrates]][["BM"]][["AICc"]] <- 9999}
		#if(model$model == "OU"){
		if(is.null(opt.list[[nrates]][["OU"]])){ opt.list[[nrates]][["OU"]][["AICc"]] <- 9999}
	}
	optSet.OU <- which.min(sapply(opt.list, function(x) x[["OU"]][["AICc"]])) # will need to dig around for beg AIC; find model model for a given best combo of breaks
	optSet.BM <- which.min(sapply(opt.list, function(x) x[["BM"]][["AICc"]]))
	
	if(opt.list[[optSet.OU]][["OU"]][["AICc"]] <	opt.list[[optSet.BM]][["BM"]][["AICc"]]) optTree[[1]] <- opt.list[[optSet.OU]][["OU"]]
	
	if(opt.list[[optSet.OU]][["OU"]][["AICc"]] >	opt.list[[optSet.BM]][["BM"]][["AICc"]]) optTree[[1]] <- opt.list[[optSet.BM]][["BM"]] 
	
	return(optTree)
}

getOptModels <- function(opt.list, model = c("BM", "OU", "BMOU")) #function will save/check rez.list to find what rate configurations have been completed in order to continue in the event of a shutdown
{
	optTree <- list()
	
	#read in optList list or optList files
	
	if(model == "BM" | model == "BMOU"){
		#isolate optList for single tree
		for(nrates in seq(1, length(opt.list), 1)) { #print(opt.list[[nrates]])
			#acquire AICc value for BM and OU
			if(is.null(opt.list[[nrates]][["BM"]])){ opt.list[[nrates]]$BM$AICc <- 9999}
		}
		optSet.BM <- which.min(sapply(opt.list, function(x) x[["BM"]][["AICc"]]))
		
		opt.listResultBM <- opt.list[[optSet.BM]][["BM"]]
		#	if(opt.list[[optSet.OU]]$OU$AICc <	opt.list[[optSet.BM]]$BM$AICc) optTree[[1]] <- opt.list[[optSet.OU]]$OU
		
		#	if(opt.list[[optSet.OU]]$OU$AICc >	opt.list[[optSet.BM]]$BM$AICc) optTree[[1]] <- opt.list[[optSet.BM]]$BM 
	}
	
	if(model == "OU" | model == "BMOU"){
		#isolate optList for single tree
		for(nrates in seq(1, length(opt.list), 1)) { #print(opt.list[[nrates]])
			#acquire AICc value for BM and OU
			if(is.null(opt.list[[nrates]][["OU"]])){ opt.list[[nrates]][["OU"]][["AICc"]] <- 9999}
		}
		optSet.OU <- which.min(sapply(opt.list, function(x) x[["OU"]][["AICc"]]))
		
		opt.listResultOU <- opt.list[[optSet.OU]][["OU"]]
		#	if(opt.list[[optSet.OU]]$OU$AICc <	opt.list[[optSet.BM]]$BM$AICc) optTree[[1]] <- opt.list[[optSet.OU]]$OU
		
		#	if(opt.list[[optSet.OU]]$OU$AICc >	opt.list[[optSet.BM]]$BM$AICc) optTree[[1]] <- opt.list[[optSet.BM]]$BM 
	}
	
	if(model == "BMOU")
	{
		opt.listResult = list(BM=opt.listResultBM, OU=opt.listResultOU)
	} else if(model == "BM"){
		opt.listResult <- opt.listResultBM
	} else if(model == "OU"){
		opt.listResult <- opt.listResultOU	
	}
	
	return(opt.listResult)
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
	
	return(as.data.frame(results.Int))
}

#######################
#look into setting a seed
######################

#runSubclade <- FALSE
#Set the code below to subsample tree for specific clades
#clade_select <- list(runSubclade=TRUE, clade_level = "order", cladename = "Perissodactyla")
#if(runSubclade == TRUE) {
#	clade_sp <- unique(occs[which(occs[clade_select$clade_level] == clade_select$cladename),c("accepted_name", "accepted_rank")])
#	clade_sp <- clade_sp[which(clade_sp$accepted_rank == "species"),]
#	clade_sp <- clade_sp[clade_sp$accepted_name %in% tree.list[[1]]$tip.label,]

#	tree_clade <- extract.clade(phy=tree.list[[1]],node = getMRCA(phy = tree.list[[1]], clade_sp$accepted_name))
#	tree.list[[1]] <- tree_clade
#make sure that dat.vec is truncated before running evo model
#}
