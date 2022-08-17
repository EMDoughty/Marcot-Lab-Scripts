regimeHist_countBox <- function(countBox = NULL, breaks = NULL, optList, thisMat, netFreq = TRUE, regimeFreq=FALSE,
                       netPlotType = "absolute", plot.together = FALSE, grayscale = TRUE)
{
  if(class(countBox) == "matrix"){
      if(is.null(breaks)) breaks <- optList_bm_allReps[[1]][[length(optList_bm_allReps)]]$optBreaks
    
      regimeBM.list <- list()
      regimeSp.List <- list()
    
    colnames(countBox) <- str_remove(colnames(countBox), " Ma")
    
    regimeBM <- countBox[,which(as.double(colnames(countBox)) > breaks[1])]
    regimeBM.list[[1]] <- apply(regimeBM, c(1), sum)
    rm(regimeBM)
    
    for(mm in seq(2,length(breaks),1))
    {
      #get regimes for remaining sections
      
      #seq(min(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),
      #max(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),1))
      regimeBM <- countBox[,which(as.double(colnames(countBox)) < breaks[mm-1] 
                                    & as.double(colnames(countBox)) > breaks[mm])]
      regimeBM.list[[mm]] <- apply(regimeBM, c(1), sum)
      
      rm(regimeBM)
    }  
    
    regimeBM <- countBox[,which(as.double(colnames(countBox)) < breaks[length(breaks)])]
    regimeBM.list[[length(breaks)+1]] <- apply(regimeBM, c(1),sum)
    
    rm(regimeBM)
    
    regimeNames <- vector()
    regimeNames[1] <- paste(">",max(breaks,sep=""))
    
    for(ii in seq(2,length(regimeBM.list)-1,1))
    {
      regimeNames[ii] <- paste(paste(breaks[ii-1]," to ",sep=""),breaks[ii],sep="")
    }
    regimeNames[length(regimeBM.list)] <- paste(breaks[length(regimeBM.list)-1],">",sep="")
    names(regimeBM.list) <- regimeNames 
    bmBreaks <- c(0, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, 4.0)
    
    btwnBreaks <- vector()
    for(ii in seq(1, length(bmBreaks)-1,1)) btwnBreaks[ii] <- (bmBreaks[ii]+bmBreaks[ii+1])/2 
    
    if(grayscale == TRUE) breakCol <- gray.colors(length(bmBreaks), start=0, end=1)
    else breakCol <- rainbow(length(bmBreaks))
    
    #make vector that conbtains the btwBreak values at qunatity listed by countBox
    for(hh in seq(1, length(regimeBM.list),1))
    {
      
      regimeBM.list[[hh]] <- c(rep(btwnBreaks[1],regimeBM.list[[hh]][1]),rep(btwnBreaks[2],regimeBM.list[[hh]][2]),
                               rep(btwnBreaks[3],regimeBM.list[[hh]][3]),rep(btwnBreaks[4],regimeBM.list[[hh]][4]),
                               rep(btwnBreaks[5],regimeBM.list[[hh]][5]),rep(btwnBreaks[6],regimeBM.list[[hh]][6]))
    }
    
    hist.list <- list()
    for(jj in seq(1, length(regimeBM.list),1))
    {
      
      hist.list[[jj]] <- hist(regimeBM.list[[jj]], breaks = bmBreaks, plot = FALSE)
    }
    
    net.plot <- hist.list
    net.plot[[length(regimeNetChange)]] <- NULL
    names(net.plot) <- breaks
    
    for(tt in seq(2,length(hist.list),1))
    {
      net.plot[[tt-1]]$counts <- hist.list[[tt]]$counts - hist.list[[tt-1]]$counts
      #tt <- tt+1
    }
    
    max(sapply(net.plot, function(x) max(x$counts)))
    min(sapply(net.plot, function(x) min(x$counts)))
    #get limits for plot dimensions and scale
    hist.plot_ylimMax <- max(sapply(hist.list, function(x) max(x$density)))
    net.plot_ylimMax <- max(sapply(net.plot, function(x) max(x$counts)))
    net.plot_ylimMin <- min(sapply(net.plot, function(x) min(x$counts)))
  
    quartz(width = 11, height = 4)
    par(mfrow = c(1, length(hist.list)),mar=c(2,3,0,0))
    for(ii in seq(1, length(hist.list),1))
    {
      if(ii == 1) axesCheck <- TRUE 
      else axesCheck <- FALSE
      plot(hist.list[[ii]], col = breakCol, axes = axesCheck, ylim = c(0,1), main = NULL)
    }
    quartz(width = 11, height = 4)
    par(mfrow = c(1, length(net.plot)+1),mar=c(2,3,0,0)) #have +1 to plot lenght so plots are same size as those for hist.plot
    for(ii in seq(1, length(net.plot),1))
    {
      if(ii == 1) axesCheck <- TRUE 
      else axesCheck <- FALSE
      plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = axesCheck, ylim = c(net.plot_ylimMin,net.plot_ylimMax), 
           main = NULL)
    }
  }
  
  return()
}
#############################################################################################################################
regimeHist_HistMedian<- function(repIntSp = NULL, breaks = NULL, optList, thisMat, netFreq = TRUE, regimeFreq=FALSE,
                                netPlotType = "absolute", plot.together = FALSE)
{
  repIntSp_regimes <- repIntSp
  if(is.null(breaks)) print("Error: Input a vector of break dates")
  
  regimeBM.list <- list()
  regimeSp.List <- list()
  all.hist <- list()
  
  for(ii in seq(1, length(repIntSp_regimes),1))
  {
    names(repIntSp_regimes[[ii]]) <- str_remove(names(repIntSp_regimes[[ii]]), " Ma")
  }
  
  for (ii in seq(1, length(repIntSp_regimes),1)) {
    regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) > breaks[1])]))
    regimeSp <- regimeSp[regimeSp %in% thisMat$species]
    regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
    
    regimeSp.List[[1]] <- regimeSp
    regimeBM.list[[1]] <- regimeBM
    
    rm(regimeSp, regimeBM)
  
    for(mm in seq(2,length(breaks),1))
    {
    #get regimes for remaining sections
    
    #seq(min(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),
    #max(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),1))
      regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] 
                                                            & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])]))
      regimeSp <- regimeSp[regimeSp %in% thisMat$species]
      regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
    
      regimeSp.List[[mm]] <- regimeSp
      regimeBM.list[[mm]] <- regimeBM
    
      rm(regimeSp, regimeBM)
    }
  
    regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[length(breaks)])]))
    regimeSp <- regimeSp[regimeSp %in% thisMat$species]
    regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
  
    regimeSp.List[[length(breaks)+1]] <- regimeSp
    regimeBM.list[[length(breaks)+1]] <- regimeBM
  
    rm(regimeSp, regimeBM)
  
    regimeNames <- vector()
    regimeNames[1] <- paste(">",max(breaks,sep=""))
    for(nn in seq(2,length(regimeBM.list)-1,1))
    {
      regimeNames[nn] <- paste(paste(breaks[nn-1]," to ",sep=""),breaks[nn],sep="")
    }
    regimeNames[length(regimeBM.list)] <- paste(breaks[length(regimeBM.list)-1],">",sep="")
    names(regimeBM.list) <- regimeNames 
    bmBreaks <- c(0, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, 4.0)
    breakCol <- rainbow(length(bmBreaks))
  
    hist.list <- list()
    #quartz()
    #par(mfrow=c(2,length(regimeBM.list)/2+1))
    for(jj in seq(1, length(regimeBM.list),1))
    {
      hist.list[[jj]] <- hist(regimeBM.list[[jj]], main=names(regimeBM.list[[jj]]), breaks = bmBreaks, 
                              las=1, col = breakCol, freq = regimeFreq, xlab="log Body Mass (kg)", plot = FALSE)
    }
    all.hist[[ii]] <- hist.list 
  }
  
  #get median across all hist counts and test to get median density
  median.density <- list()
  median.counts <- list()
  
  for(rr in seq(1, length(all.hist[[1]]),1)) #to go through each regime
    {
    regimeDensityMat <- matrix(nrow = length(all.hist),ncol= length(bmBreaks)-1)
    regimeCountMat <- matrix(nrow = length(all.hist),ncol= length(bmBreaks)-1)
    
    for(gg in seq(1, length(all.hist),1))
    {
      regimeDensityMat[gg,] <- all.hist[[gg]][[rr]]$density 
      regimeCountMat[gg,] <- all.hist[[gg]][[rr]]$counts
    }
    
    median.density[[rr]] <- apply(regimeDensityMat,2,median)
    median.counts[[rr]] <- apply(regimeCountMat,2,median)
  }
  hist.plot <- hist.list
  
  #overwite with median density and counts
  quartz(width = 11, height = 4)
  par(mfrow = c(1, length(hist.plot)),mar=c(0,0.5,0,0))
  for(pp in seq(1, length(hist.plot),1))
  {
    hist.plot[[pp]]$density <- median.density[[pp]]
    hist.plot[[pp]]$counts <- median.counts[[pp]]
    #plot(hist.plot[[pp]], col = breakCol, axes = FALSE)
  }
  
  #get net change between regimes
  net.plot <- hist.plot; net.plot[[length(hist.plot)]] <- NULL
  for(ii in seq(1, length(net.plot),1))
  {
    net.plot[[ii]]$counts <- hist.plot[[ii+1]]$counts-hist.plot[[ii]]$counts
    #plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = FALSE)
  }
  max(sapply(net.plot, function(x) max(x$counts)))
  min(sapply(net.plot, function(x) min(x$counts)))
  #get limits for plot dimensions and scale
  hist.plot_ylimMax <- max(sapply(hist.plot, function(x) max(x$density)))
  net.plot_ylimMax <- max(sapply(net.plot, function(x) max(x$counts)))
  net.plot_ylimMin <- min(sapply(net.plot, function(x) min(x$counts)))
  
  #plot both the hist and net changes
  quartz(width = 11, height = 4)
  par(mfrow = c(1, length(hist.plot)),mar=c(0,0.5,0,0))
  for(ii in seq(1, length(hist.plot),1))
  {
    plot(hist.plot[[ii]], col = breakCol, axes = TRUE, ylim = c(0,hist.plot_ylimMax), main = NULL)
  }
  quartz(width = 11, height = 4)
  par(mfrow = c(1, length(net.plot)),mar=c(0,0.5,0,0))
  for(ii in seq(1, length(net.plot),1))
  {
    plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = TRUE, ylim = c(net.plot_ylimMin,net.plot_ylimMax), 
         main = NULL)
  }
 return() 
}


