install.packages('ape')
install.packages('phytools')

require(ape)
require(phytools)
tree_base <- read.nexus("~/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_1_25_UngulataBackBone")

#clade_sp <- read.csv("~/Dropbox/ungulate_RA/2017_1_26_Clade_species.csv", stringsAsFactors = FALSE)
#clade_sp

clade_sp <- read.csv("C:/Users/Evan/Dropbox/ungulate_RA/2017_1_26_Clade_species.csv",stringsAsFactors = FALSE)
clade_sp

#MCRA_Codes <- read.csv("~/Dropbox/ungulate_RA/2017_1_26_MCRA_Codes.csv", stringsAsFactors = FALSE)
#MCRA_Codes

MCRA_Codes <- read.csv("C:/Users/Evan/Dropbox/ungulate_RA/2017_1_26_MCRA_Codes.csv", stringsAsFactors = FALSE)
MCRA_Codes

iterations <- function()
{
#set number of iterations for program to run
iter <- 1
iter<-readline(prompt="Enter number of iterations: ")
as.integer(iter)
return (iter)
}

iter <- iterations()

#set this to be a function
tree_resolution <- function(backbone_tree,Species_file,MRCA_file,iterations)
{
#Resolve backbone Tree
tree_resolv <- multi2di(backbone_tree)
plot(tree_resolv, font = 1, edge.width = 0.5, cex =0.2, label.offset = 1, adj = 0)
names(tree_resolv)
tree_resolv

k <- 1 #track number of iterations completed

while(k < iterations+1) {
print("Iteration: ")
print(k)

#Identify landing edges
i <- 1
node <- vector()
#node <- c(0,0,0,0)#need to set so sample is not being selected as 0 so wil need to be defined in while loops near numberPos
nodechoice <- 0

while(i < nrow(MRCA_file) + 1) {
	print('i')
	print(i)
	clade_tree <- read.nexus(file = MRCA_file$Filename[i])
  clade_tree <- multi2di(clade_tree)
	plot(clade_tree, font = 0.1, cex =0.2, label.offset = 1, adj = 0)
	numberPos <- MRCA_file[i,"X._positions"]
	numberPos 	
	j <- 5
	m <- 1
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
		typeof(speciesA)
		print(m)
		node[m] <- fastMRCA(tree=tree_resolv, sp1=speciesA, sp2=speciesB) #fastMRCA breaking at x<match(...) and y<-match(..) prior to passing to getAncestor() internal function keeps breaking
		node
		m <- m+1
		j <- j+1}
#Randomly place taxa on landing branches
	node_choice <- sample(node,size = 1, replace = FALSE, prob =NULL)
	node_choice
	typeof(node_choice)
#Attach Subtree
	tree_resolv <- bind.tree(x=tree_resolv, y=clade_tree, position = node_choice)
		
	plot(tree_resolv, font = 3, cex =0.2, label.offset = 1, adj = 0)
	i <- i + 1
	}

#Troubleshooting: Paint Branches if issues occur
#this.map <- paintSubTree(tree=tree_resolv, node=fastMRCA(tree=tree_resolv, sp1="Hypertragulus_minor", sp2="Yumaceras_ruminalis"), state="2")
#plotSimmap(this.map,fsize=0.1, ftype = "i")

plot(tree_resolv, font = 3, cex =0.2, label.offset = 1, adj = 0)

return(tree_resolv)
}

save_tree <- tree_resolution(backbone_tree = tree_base,Species_file = clade_sp,MRCA_file = MCRA_Codes,iterations = iter)
save_tree

#one tree for the function

#save tree to variable and file
name <- paste("tree_resolve",k)
#filename <- paste("~/Dropbox/ungulate_RA/NAUngulata_Trees/ResolvedTrees/",name, collapse=" ") #pathing on mac
filename <- paste("C:/Users/Evan/Dropbox/ungulate_RA/NAUngulata_Trees/ResolvedTrees/",name, collapse=" ") #for pathing on personal laptop
filename
filename <- paste(filename,".nex")
filename
filename <-gsub(" ","", filename)
filename
write.nexus(save_tree, file=filename, translate = TRUE)

k<-k+1
}



