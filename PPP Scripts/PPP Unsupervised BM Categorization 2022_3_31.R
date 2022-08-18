mom.data <- read.csv("/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/MOM_v4.1.csv")
#mom.data <- read.csv("/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/MOM_v4.1_herb.csv")


mom.data <- mom.data[,c("Continent", "Species.Distribution", "Status", 
                        "Order", "FAMILY", "Genus", "Species", 
                        "Combined.Mass..g.","LogMass..g.", "Combined.Mass..kg.", "LogMass..kg."),]

mom.data <- mom.data[mom.data$Order %in% c("Artiodactyla", "Perissodactyla"),]

### set groups using Janis 2000 breaks
bmBreaks_herb <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
sizecateg <- c("<5 kg", 
               "5-25 kg",
               "25-150kg", 
               "150-500 kg", 
               ">500 kg")
mom.data$Janis.cat <- NA

for(ii in seq(1, length(sizecateg),1)) 
{
  mom.data[mom.data$LogMass..kg. > bmBreaks_herb[ii] & mom.data$LogMass..kg. < bmBreaks_herb[ii+1],]$Janis.cat <- sizecateg[ii]
}

colors.Janis <- rainbow(length(sizecateg)); names(colors.Janis) <- sizecateg

kmeans_results <- list()
kmeans_breaks <- list()

wss <- vector(length = 9)
sil <- vector(length = 9)

for(ii in seq(2,10,1))
{
  kmeans_results[[ii]] <- kmeans(mom.data$LogMass..kg., centers = ii, iter.max = 10000, nstart = 50)
  wss[ii-1] <- kmeans_results[[ii]]$tot.withinss
  sil[ii-1] <- mean(silhouette(kmeans_results[[ii]]$cluster, dist(mom.data$LogMass..kg.))[,3])
  
  kmean_Onebreaks <- matrix(ncol = 3, nrow = ii)
  colnames(kmean_Onebreaks) <- c("k", "min.kg", "max.kg")
  
  for(jj in seq(1, max(unique(kmeans_results[[ii]]$cluster)),1)) 
  {
    min.categ <- min(mom.data[kmeans_results[[ii]]$cluster == jj,"LogMass..kg."])
    max.categ <- max(mom.data[kmeans_results[[ii]]$cluster == jj,"LogMass..kg."])
    
    kmean_Onebreaks[jj,] <- c(jj, min.categ, max.categ)
  }
  kmeans_breaks[[ii]] <- kmean_Onebreaks 
}

#good resource and code to follow
##https://uc-r.github.io/kmeans_clustering

#Determine Optimal K value
par(mfrow=c(1,2))
##Elbow Method
plot(2:10, wss, type = "b", pch = 19,
     xlab = "Number of clusters K", ylab = "Total Within Sum of Square")
abline(v=5)
##Silhoettes
require(cluster)
plot(2:10, sil, type = "b", pch = 19,
     xlab = "Number of clusters K", ylab = "Average Silhoettes")
abline(v=5)

#Plot PPP pred Categories vs Optimal K Means
par(mfrow=c(1,3))
plot(mom.data$LogMass..kg., pch = 19, col = colors.Janis[mom.data$Janis.cat], main = "Janis 2000")
abline(h = bmBreaks_herb)

plot(mom.data$LogMass..kg., pch = 19, col = rainbow(5)[kmeans_results[[5]]$cluster], main = "K_5")
#abline(h= kmeans_results[[5]]$centers)
abline(h = c((kmeans_breaks[[5]][4,"max.kg"] + kmeans_breaks[[5]][3,"min.kg"])/2,
             (kmeans_breaks[[5]][3,"max.kg"] + kmeans_breaks[[5]][2,"min.kg"])/2,
             (kmeans_breaks[[5]][2,"max.kg"] + kmeans_breaks[[5]][1,"min.kg"])/2,
             (kmeans_breaks[[5]][1,"max.kg"] + kmeans_breaks[[5]][5,"min.kg"])/2)
       )

plot(mom.data$LogMass..kg., pch = 19, col = rainbow(6)[kmeans_results[[6]]$cluster], main = "K_6")
#abline(h= kmeans_results[[6]]$centers)
#abline(h = bmBreaks_herb)
abline(h = c((kmeans_breaks[[6]][1,"max.kg"] + kmeans_breaks[[6]][4,"min.kg"])/2,
             (kmeans_breaks[[6]][4,"max.kg"] + kmeans_breaks[[6]][3,"min.kg"])/2,
             (kmeans_breaks[[6]][3,"max.kg"] + kmeans_breaks[[6]][6,"min.kg"])/2,
             (kmeans_breaks[[6]][6,"max.kg"] + kmeans_breaks[[6]][5,"min.kg"])/2,
             (kmeans_breaks[[6]][5,"max.kg"] + kmeans_breaks[[6]][2,"min.kg"])/2)
       )

#################################################################################################################

####### Predators

#removed taxa that were under scrutiny and those with a weight of -999
mom.data <- read.csv("/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/MOM_v4.1_pred.csv")
mom.data <- mom.data[!mom.data$FAMILY %in% c("Otariidae","Odobenidae","Phocidae"),]

mom.data <- mom.data[,colnames(mom.data) %in% c("Status", 
                        "Order", "FAMILY", "Genus", "Species", 
                        "Combined.Mass..g.","LogMass..g.", "Combined.Mass..kg.", "LogMass..kg.")]

mom.data <- mom.data[mom.data$Order %in% c("Carnivora"),]

### set groups using PPP prey breaks
bmBreaks_pred <- c(-Inf, 0, 0.845098, 1.322219, 2, Inf)  #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
sizecateg <- c("<1kg",
               "1-7kg",
               "7-21kg", 
               "21-100kg",
               ">100kg")
mom.data$Pred.cat <- NA

for(ii in seq(1, length(sizecateg),1)) 
{
  mom.data[mom.data$LogMass..kg. > bmBreaks_pred[ii] & mom.data$LogMass..kg. < bmBreaks_pred[ii+1],]$Pred.cat <- sizecateg[ii]
}

colors.pred <- rainbow(length(sizecateg)); names(colors.pred) <- sizecateg

kmeans_results <- list()
kmeans_breaks <- list()

wss <- vector(length = 9)
sil <- vector(length = 9)

for(ii in seq(2,10,1))
{
  kmeans_results[[ii]] <- kmeans(mom.data$LogMass..kg., centers = ii, iter.max = 10000, nstart = 50)
  wss[ii-1] <- kmeans_results[[ii]]$tot.withinss
  sil[ii-1] <- mean(silhouette(kmeans_results[[ii]]$cluster, dist(mom.data$LogMass..kg.))[,3])
  
  kmean_Onebreaks <- matrix(ncol = 3, nrow = ii)
  colnames(kmean_Onebreaks) <- c("k", "min.kg", "max.kg")
  
  for(jj in seq(1, max(unique(kmeans_results[[ii]]$cluster)),1)) 
  {
    min.categ <- min(mom.data[kmeans_results[[ii]]$cluster == jj,"LogMass..kg."])
    max.categ <- max(mom.data[kmeans_results[[ii]]$cluster == jj,"LogMass..kg."])
    
    kmean_Onebreaks[jj,] <- c(jj, min.categ, max.categ)
  }
  kmeans_breaks[[ii]] <- kmean_Onebreaks 
}

#good resource and code to follow
##https://uc-r.github.io/kmeans_clustering

#Determine Optimal K value
par(mfrow=c(1,2))
##Elbow Method
plot(2:10, wss, type = "b", pch = 19,
     xlab = "Number of clusters K", ylab = "Total Within Sum of Square")
abline(v=6)
##Silhoettes
plot(2:10, sil, type = "b", pch = 19,
     xlab = "Number of clusters K", ylab = "Average Silhoettes")
abline(v=6)

#Plot PPP pred Categories vs Optimal K Means
par(mfrow=c(1,3))
plot(mom.data$LogMass..kg., pch = 19, col = colors.pred[mom.data$Pred.cat], main = "PPP Categories")
abline(h = bmBreaks_pred)

plot(mom.data$LogMass..kg., pch = 19, col = rainbow(2)[kmeans_results[[2]]$cluster], main = "K_2")
#abline(h= kmeans_results[[5]]$centers)
abline(h = c((kmeans_breaks[[2]][2,"max.kg"] + kmeans_breaks[[2]][1,"min.kg"])/2))

plot(mom.data$LogMass..kg., pch = 19, col = rainbow(6)[kmeans_results[[6]]$cluster], main = "K_6")
#abline(h= kmeans_results[[6]]$centers)
#abline(h = bmBreaks_herb)
abline(h = c((kmeans_breaks[[6]][6,"max.kg"] + kmeans_breaks[[6]][1,"min.kg"])/2,
             (kmeans_breaks[[6]][1,"max.kg"] + kmeans_breaks[[6]][4,"min.kg"])/2,
             (kmeans_breaks[[6]][4,"max.kg"] + kmeans_breaks[[6]][3,"min.kg"])/2,
             (kmeans_breaks[[6]][3,"max.kg"] + kmeans_breaks[[6]][2,"min.kg"])/2,
             (kmeans_breaks[[6]][2,"max.kg"] + kmeans_breaks[[6]][5,"min.kg"])/2))



