#HackBio Transcriptomics Project

# Load required libraries
library(readr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(devtools)
library(ggpubr)
library(tidyr)
library(dplyr)
library(data.table)
library(tidyverse)
library (useful)

# Import .csv file in R
asthma <- read_delim("Asthma.txt", delim = "\t", 
                     escape_double = FALSE, trim_ws = TRUE)

# Set scale
# before the PCA analysis
res.pca <- prcomp(asthma,  scale = TRUE)
  
# Default plot
fviz_pca_ind(res.pca)

# Change title and axis labels
fviz_pca_ind(res.pca) +
  labs(title ="PCA", x = "PC1", y = "PC2")

# Change axis limits by specifying the min and max
fviz_pca_ind(res.pca) +
  xlim(-200, 50) + ylim(-40, 40)

# Use text only
fviz_pca_ind(res.pca, geom="text") +
  xlim(-200, 50) + ylim(-40, 40)
  

# Use points only
fviz_pca_ind(res.pca, geom="point") + 
  xlim(-200, 50) + ylim(-40, 40)

# # Change the size of points
fviz_pca_ind(res.pca, geom="point", pointsize = 2) +
  xlim(-200, 50) + ylim(-40, 40)

# Change point color and theme
fviz_pca_ind(res.pca, col.ind = "blue")+
  xlim(-200, 50) + ylim(-40, 40)
  theme_minimal()

# Control automatically the color of individuals
# using the cos2 or the contributions
# cos2 = the quality of the individuals on the factor map
fviz_pca_ind(res.pca, col.ind="cos2") +
  xlim(-200, 50) + ylim(-40, 40)

# Gradient color
fviz_pca_ind(res.pca, col.ind="cos2") +
  xlim(-200, 50) + ylim(-40, 40)
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)

# Change the theme and use only points
fviz_pca_ind(res.pca, col.ind="cos2", geom = "point") +
  xlim(-200, 50) + ylim(-40, 40) +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)+ theme_minimal()

# Color by the contributions
fviz_pca_ind(res.pca, col.ind="contrib") +
  xlim(-200, 50) + ylim(-40, 40) +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=4)

# Control the transparency of the color by the
# contributions
fviz_pca_ind(res.pca, alpha.ind="contrib") +
  xlim(-200, 50) + ylim(-40, 40) +
  theme_minimal()

# Select and visualize individuals with cos2 > 0.96
fviz_pca_ind(res.pca, select.ind = list(cos2 = 0.96))



######################## fviz_pca_var(): Graph of variables #################
# Default plot
fviz_pca_var(res.pca) +
  xlim(-2, 1)

# Use points and text
fviz_pca_var(res.pca, geom = c("point", "text")) +
  xlim(-2, 1) 

# Change color and theme
fviz_pca_var(res.pca, col.var="steelblue") +
  xlim(-2, 1)
  theme_minimal()  
  

############################# fviz_pca_biplot(): Biplot of individuals of variables #######################
fviz_pca_biplot(res.pca) +
    xlim(-200, 50) + ylim(-40, 40)
  
# Keep only the labels for variables
fviz_pca_biplot(res.pca, label ="var") +
  xlim(-200, 50) + ylim(-40, 40)

# Keep only labels for individuals
fviz_pca_biplot(res.pca, label ="ind") +
  xlim(-200, 50) + ylim(-40, 40)

# Control automatically the color of individuals using the cos2
fviz_pca_biplot(res.pca, label ="var", col.ind="cos2") +
  xlim(-200, 50) + ylim(-40, 40) +
  theme_minimal()




######################################################################################
######################################################################################
############################# METHYLATION ANALYSIS ###################################
######################################################################################
######################################################################################
# Import .csv file in R
Methylation <- read_delim("Asthma_Methylation.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

# For the purpose of visualization, I will like to use any two columns that correlate
cor(Methylation)
plot (Methylation$patient1, Methylation$patient2)
plot (Methylation$patient1, Methylation$patient3)

# So we can have same results
set.seed(102)

# Choosing the right number of clusters
MethylationBEST <- FitKMeans(Methylation, max.clusters=10) 
PlotHartigan(MethylationBEST)

# Perform KMEAN Clustering
MethylationK3 <- kmeans (x = Methylation, centers = 8)

# Well, I can also add this cluster information to my dataset
Methylation$clusters <- c(MethylationK3$cluster)

# To visualize this cluster information and compare patients
ggplot(Methylation, aes(x = patient1, y = patient2, color = factor(clusters))) + geom_point() + theme_bw()
ggplot(Methylation, aes(x = patient2, y = patient5, color = factor(clusters))) + geom_point() + theme_bw()
