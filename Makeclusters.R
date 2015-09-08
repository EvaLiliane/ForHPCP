rm(list = ls())

# Read data in
#indout <- read.csv("~/Documents/PhD/Codes/IndividualData10-New.csv")
plot(indout)


setwd("/home/evaliliane/Documents/PhD/Codes/IndMod")

png("NewScatterPlot_Estimates.png", width = 800, height = 600)
plot(indout[3:6], main="Scatterplot of Estimated Parameters")
dev.off()

png("NewBoxPlot_Estimates.png", width = 800, height = 600)
boxplot(indout[3:6], main="Boxplot of Estimated Parameters")
dev.off()



pat <- dat$patient
newtest <- mydata[mydata$patient %in% pat,]
plot(newtest)

######################### CLUSTERS ########################################################################
# ----------------------------------------
# HC clustering. 
# ----------------------------------------
# It is defined to work on a large dataset. 
d <- dist(as.matrix(indout[1:2]))
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
