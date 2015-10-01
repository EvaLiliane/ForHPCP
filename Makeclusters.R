rm(list = ls())

# Arguments 
args <- commandArgs(TRUE) # Should be 4 arguments.

if(length(args) == 0){
  print("No arguments supplied.")
  ##supply default values
  infile <- "/home/evaliliane/Documents/PhD/HPCP/Output/ChildrenIndFitData.csv"
  outfile <- "ChildrenClusters.csv"  
  outfig <- "Children"
} else {
  infile <- eval( parse(text=args[1]))
  outfile <- eval( parse(text=args[2]))
  outfig <- eval( parse(text=args[3]))
}

# Read data in
indout <- read.csv(infile, header = T)
#plot(indout)

setwd("/home/evaliliane/Documents/PhD/HPCP/Output")

namefig1 <- paste("NewScatterPlot_Estimates_",outfig,".png",sep="")
png(namefig1, width = 800, height = 600)
plot(indout[3:8], main="Scatterplot of Estimated Parameters")
dev.off()

namefig2 <- paste("NewBoxPlot_Estimates_",outfig,".png",sep="")
png(namefig2, width = 800, height = 600)
boxplot(indout[3:8], main="Boxplot of Estimated Parameters")
dev.off()

######################### CLUSTERS ########################################################################
# ----------------------------------------
# HC clustering. 
# ----------------------------------------
# It is defined to work on a large dataset. 
d <- dist(as.matrix(indout[,c("patient","RSS")]))
system.time(hc <- hclust(d, "centroid") )

# # Look at the plot to choose number of patient
plot(hc,labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = NULL, xlab = NULL, ylab = "Height")


cutree(hc, h = ) # cut tree into k clusters
# # Draw dendogram with red borders around the k clusters
# rect.hclust(hc, k=3, border="red") 
# # divide rows into k clusters


library(cluster)
pam(indout[,c("patient","RSS")])
silhouette(hc)
