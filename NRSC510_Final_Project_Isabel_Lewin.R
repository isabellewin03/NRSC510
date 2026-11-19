# NRSC510 final project in R 

# NRSC510 class project 
## Isabel Lewin 

# installing new packages that I will need for PCA + data visualization
#install.packages("factoextra") # for some PCA visualizations
#install.packages("MASS")
#install.packages("plotly") # for 3D scatter plot 
#install.packages("ggbiplot") # for creating PCA biplots using ggplot
#install.packages("clValid") # for calculating the Dunn's index
#install.packages("fastcluster") # for faster heirarchical clustering
#install.packages("uwot") # for umap 


# loading packages that I will use for this project
library(vroom)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(factoextra)
library(MASS)
library(cluster) 
library(fastcluster)
library(uwot)

# ------------------------------------------------------------------------------
# DATA CLEANING AND SET-UP

# to copy file path for working directory: 
## option + command + c
setwd('/Users/isabellewin/Documents/UBC MSc GPN/NRSC510/NRSC510 final project')

# reading in dataset (alex's static microglia morphology data) using the vroom function
static_microglia <- vroom(delim = ',',
                          #header = TRUE,
                          'final-dataframe.csv')

# checking dataset format to ensure that it is correct 
View(static_microglia)

# checking the variable type of each column using glimpse()
glimpse(static_microglia)
# the majority of microglia characteristics are double-precision floating-point number <dbl>

# using the summary function to examine characteristics of each variable more closely 
summary(static_microglia)

# removing the 'cluster' column because this was created by a previous student in my lab after conducting K-means clustering 
## the 'cluster' column was column #35 so I used the indexing brackets to remove it and replaced the old dataset with a new one that doesn't contain #35
### I also removed any row from the BrainRegion column that were not HYPO (hypothalamus)
#### this is because the student who generated these data also looked at specific hypothalamic nuclei, which I don't want to examine
static_microglia_filtered <- static_microglia[,-35]%>%
  filter(BrainRegion == 'HYPO')%>%
  rename(Hull.pixels = Density.of.foreground.pixels.in.hull.area,
         Ratio.hull.axes = Span.ratio.of.hull..major.minor.axis., 
         Max.hull.span = Maximum.span.across.hull,
         Width.bound.rect = Width.of.bounding.rectangle, 
         Ht.bound.rect = Height.of.bounding.rectangle,
         Max.hull.radius = Maximum.radius.from.hull.s.center.of.mass,
         Max.min.hull.radius = Max.min.radii.from.hull.s.center.of.mass, 
         CV.hull.radius = Relative.variation..CV..in.radii.from.hull.s.center.of.mass, 
         Diameter.bound.circle = Diameter.of.bounding.circle,
         Max.circle.radius = Maximum.radius.from.circle.s.center.of.mass,
         Max.min.circle.radius = Max.min.radii.from.circle.s.center.of.mass,
         CV.circle.radius = Relative.variation..CV..in.radii.from.circle.s.center.of.mass,
         Mean.circle.radius = Mean.radius.from.circle.s.center.of.mass, 
         Branches = X..of.branches, 
         Junctions = X..of.junctions, 
         End.pt.vox = X..of.end.point.voxels,
         Junction.vox = X..of.junction.voxels,
         Slab.vox = X..of.slab.voxels, 
         Mean.branch.length = Average.branch.length,
         Triple.pts = X..of.triple.points, 
         Quad.pts = X..of.quadruple.points, 
         Max.branch.length = Maximum.branch.length)

# checking to make sure indexing and filtering step worked
glimpse(static_microglia_filtered) 
# these steps worked! number of columns reduced by 1, rows reduced by ~2000

# ------------------------------------------------------------------------------
# INITIAL PCA STEP

# i am going to use PCA for the first part of my analysis for this dataset, and before i do that, i need to standardize my data
## this will prevent any variables with large values from skewing the PCA results

# my first step is to prevent the categorical variables from being standardized: 
## to do this i am using the sapply() function to isolate the numeric columns and make them into a vector
static_microglia_num_cols <- sapply(static_microglia_filtered, is.numeric)

# now I will use the lapply() function to standardize the numeric columns in my dataset 
## i will perform z-standardization, using the following operation: z = (value - mean)/standard deviation, where z is the new, standardized value 
### the scale() function automatically performs z-scoring so I will use it for this case
static_microglia_filtered[static_microglia_num_cols] <- lapply(static_microglia_filtered[static_microglia_num_cols], 
                                                               function(x) as.numeric(scale(x)))

# checking to make sure my standardized dataset looks correct
View(static_microglia_filtered)
# looks good! all numeric variables have been standardized by z-scoring

# now I will run a PCA using the function prcomp()
## here, I excluded columns 1-7 from analysis because they are categorical variables and I don't want them in my analysis
static_microglia_pca <- prcomp(static_microglia_filtered[,-c(1:7)], 
                               scale = FALSE)

# checking out my data to make sure it looks right after doing PCA 
summary(static_microglia_pca)
# looks good!

# this is going to print out the names of more data saved in static_microglia_pca
names(static_microglia_pca)
# looking at each one individually: 
static_microglia_pca$sdev # this is the standard deviation of the components 
static_microglia_pca$rotation # these are the eigenvectors 
static_microglia_pca$center # this is the standard deviation of the variables 
static_microglia_pca$scale # this is the mean of the variables 
static_microglia_pca$x # these are the PCA scores for all of the components 

# generating a scree plot of to display the Eigenvalues for each principal component in my PCA 
## this will help me decide how many components to retain in my analysis 
### the function fviz_eig() is a part of the factoextra package 
fviz_eig(static_microglia_pca, 
         addlabels = TRUE, 
         ylim = c(0, 60),
         barfill = 'skyblue',
         barcolor = 'black',
         main = "Microglia Morphology PCA Scree Plot", 
         ggtheme = theme_classic())
# from this plot, it seems as though PCs 1-3 explain most of the variance (~82%)
## PC1 and PC2 alone explain 75.3% 

# i am also going to make an elbow plot with the Kaiser criterion cutoff (red line)
## this cutoff states that you should only retain PCs with a variance percentage (Eigenvalue) > 1
pca_scree <- data.frame(
  PC = 1:length(static_microglia_pca$sdev),
  Eigenvalue = static_microglia_pca$sdev^2)

ggplot(pca_scree, aes(x = PC, y = Eigenvalue)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 3) +
  # Add Kaiser Criterion line at y = 1
  geom_hline(yintercept = 1, linterp = "dashed", color = "red", linewidth = 1) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Eigenvalue") +
  theme_minimal()
# this shows me that I need to retain PCs 1-3 

# since these are the PCs I'm mainly interested in, I can make a component matrix table with them alone for easier viewing
## component matrix for PCs 1-3
component_matrix_PC123 <- static_microglia_pca$rotation[,-c(4:27)] # these are the eigenvectors
View(component_matrix_PC123)

# now I am going to make biplots to visualize the spread of my data with PCA
# creating a biplot that overlays PC1 and PC2 
ggbiplot(static_microglia_pca,
         groups = static_microglia_filtered$Treatment,
         ellipse = TRUE,
         var.axes = TRUE,
         varname.size = 2,
         varname.adjust = 2,
         alpha = 0.4,
         repel = TRUE, 
         vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
  scale_fill_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC1 and PC2",
       color = "Treatment", fill = "Treatment") +
  coord_fixed()
# I can't get the labels not to overlap, so I will have to save the plot as an EPS file and edit the locations of the labels

# # creating a biplot that overlays PC2 and PC3, again visualizing only the top 8 contributing variables 
ggbiplot(static_microglia_pca,
         choices = c(2, 3), # specifying PCs 2 and 3
         groups = static_microglia_filtered$Treatment,
         ellipse = TRUE,
         var.axes = TRUE,
         varname.size = 2,
         varname.adjust = 2,
         alpha = 0.4,
         repel = TRUE, 
         vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
  scale_fill_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC1 and PC2",
       color = "Treatment", fill = "Treatment") +
  coord_fixed()

# i can use these plots (which both use functions in the same family as above) to visualize the data points and variables separately
## data points
fviz_pca_ind(static_microglia_pca, geom = "point",
             habillage = static_microglia_filtered$Sex, addEllipses = TRUE)
## variables 
fviz_pca_var(static_microglia_pca, col.var = "contrib",
             select.var = list(contrib = 8), 
             gradient.cols = c("blue", "orange", "red"), repel = TRUE)


# now I can also create a 3D scatterplot to look at the spread of PCs 1-3
## i will do this using the plotly package 
library(plotly)

# identifying the scores for the first 3 PCs, plus the grouping variable (treatment)
scores <- as.data.frame(static_microglia_pca$x[, 1:3])
scores$Treatment <- static_microglia_filtered$Treatment

# generating the variance explained, for the axis labels
ve <- round(summary(static_microglia_pca)$importance[2, ] * 100, 1)

# creating the 3D scatter plot 
fig <- plot_ly(scores,
               x = ~PC1, y = ~PC2, z = ~PC3,
               color = ~Treatment,
               colors = c(CON = "#E69F00", COLD = "#0072B2"),
               type = "scatter3d",
               mode = "markers",
               marker = list(size = 1, opacity = 0.8))

fig <- fig %>% layout(
  title = "Microglia Morphology PCA (PC1-PC3)",
  scene = list(
    xaxis = list(title = paste0("PC1 (", ve[1], "%)")),
    yaxis = list(title = paste0("PC2 (", ve[2], "%)")),
    zaxis = list(title = paste0("PC3 (", ve[3], "%)"))
  )
)

fig

# if i want to save this as an html file, i can use this code
## library(htmlwidgets)
## saveWidget(fig, "pca_3d.html", selfcontained = TRUE)

# ------------------------------------------------------------------------------
# PROCESSING MY DATA MORE TO REMOVE REDUNDANT VARIABLES 

# after completing PCA on my data, it is clear that there are several redundant/correlating variables within my dataset
# now I am going to remove any redundant variables and examine the resulting PCA to see if this can help clarify the analysis

## first, I am going to check variable correlation empirically using cor(), which uses Pearson correlation
library(corrplot)
# isolating numbered columns from my data
static_microglia_num_cols <- sapply(static_microglia_filtered, is.numeric)

# using the cor() function to see which variables correlate with each other 
static_microglia_cor <- cor(static_microglia_filtered[static_microglia_num_cols], use = "pairwise.complete.obs")
# plotting this as a correlation matrix
## a Pearson correlation coefficient of |r| = 1 means exact correlation 
### r = +1 means a positive correlation, but r = -1 means a negative correlation
corrplot(static_microglia_cor, method = "color", type = "upper", tl.cex = 0.6) 

# to make this simpler, I am going to list variable pairs with |r| > 0.9, which I am defining as very high correlation
high <- which(abs(static_microglia_cor) > 0.90 & upper.tri(static_microglia_cor), arr.ind = TRUE)
data.frame(var1 = rownames(static_microglia_cor)[high[, 1]],
           var2 = colnames(static_microglia_cor)[high[, 2]],
           r = static_microglia_cor[high])
# from this list and the above plot, I can see that many variables correlate empirically
## but it is still hard to determine what variables need to be removed from this list 

# for a better visual, I am going to use a cluster dendrogram to visualize the different blocks of variables that are redundant 
## this method converts the Pearson correlation coefficients from above into a distance clustering format (by calculating cluster score = 1-|r|), and then creates a dendrogram 
### again, I want to look at variables with high correlation (|r| > 0.9), which in this format translates to cluster score > 0.1 
dist_mat <- as.dist(1 - abs(static_microglia_cor))
hc <- hclust(dist_mat, method = "complete")
plot(hc, cex = 0.7)
abline(h = 0.1, col = "red", lty = 2)   # this red line marks |r| = 0.9

# I am also going to print out the names of the variables in each of these clusters below |r| = 0.9
groups <- cutree(hc, h = 0.1)
split(names(groups), groups) 

# now I can make decisions about what variables to cut based on empirical evidence and my knowledge about how the variables were collected
# ie. I will remove variables if they are both highly correlated with others AND are measuring similar elements

# from group '3', I will keep both "Ratio.of.hull.axes" and "Circularity", because they are negatively correlated and correlation is on the low side |r| = 0.90       
## from group '4', I will keep "Perimeter", "Max.span.across.hull", "Mean.radius", "Mean.radius.from.circle.center", and remove the rest
### from group '5', I will keep both "Area" and "Number.of.slab.voxels", because they are not redundant measurements
#### from group '12', I will keep "Number.of.branches" and "Number.of.triple.points" remove "Number.of.junctions" and "Number.of.junction.voxels", which are very similar measurements to "Number.of.branches" 

# this brings my total number of (numeric) variables from 27 to __
keep <- c("BrainRegion", "MouseID", "Treatment", "Sex", "SlideNumber", "ID", "UniqueID", 
          "Foreground.pixels", "Hull.pixels", "Ratio.hull.axes", "Circularity", "Perimeter",
          "Mean.radius", "Mean.circle.radius", "Area", "Slab.vox",
          "Width.bound.rect", "Ht.bound.rect", "Max.min.hull.radius", 
          "CV.hull.radius", "Max.min.circle.radius",
          "CV.circle.radius", "Branches", "Triple.pts", "End.pt.vox", 
          "Mean.branch.length", "Quad.pts", "Max.branch.length") 
# creating a new dataframe that only keeps the above 18 numeric variables 
static_microglia_non_red  <- static_microglia_filtered[keep]

# checking to make sure the data were processed correctly 
View(static_microglia_non_red)
glimpse(static_microglia_non_red)
# looks good!

# now that I have removed these highly correlated variables from my dataset, I am going to perform PCA again using the steps from above 
static_microglia_non_red_pca <- prcomp(static_microglia_non_red[,-c(1:7)], 
                                   scale = FALSE)

# generating a scree plot to help me decide how many components to retain in my analysis 
fviz_eig(static_microglia_non_red_pca, 
         addlabels = TRUE, 
         ylim = c(0, 60),
         barfill = 'skyblue',
         barcolor = 'black',
         main = "Microglia Morphology PCA Scree Plot", 
         ggtheme = theme_classic())

# i am also going to make an elbow plot 
pca_scree <- data.frame(
  PC = 1:length(static_microglia_non_red_pca$sdev),
  Eigenvalue = static_microglia_non_red_pca$sdev^2)
ggplot(pca_scree, aes(x = PC, y = Eigenvalue)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Eigenvalue") +
  theme_minimal()
# based on the scree and elbow plots, PCs 1-4 explain most of the variance (83.6% of total variance)


# now I am going to make biplots to visualize the spread of my data with PCA
# biplot that overlays PC1 and PC2 
biplot_PCA_12 <- ggbiplot(static_microglia_non_red_pca,
                           groups = static_microglia_non_red$Treatment,
                           ellipse = TRUE,
                           var.axes = TRUE,
                           varname.size = 2,
                           varname.adjust = 2,
                           alpha = 0.4,
                           repel = TRUE, 
                           vector.args = list(linewidth = 0.3)) +
                    scale_color_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
                    scale_fill_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
                    guides(shape = "none") +
                    theme_minimal() +
                    labs(title = "Biplot of PC1 and PC2",
                         color = "Treatment", fill = "Treatment") +
                    coord_fixed()
biplot_PCA_12

# saving biplot as an svg file so I can edit the positioning of the variable labels which are currently overlapping in R 
## when I need to do edit the labels on the other plots, I will use these same lines of code
#install.packages("svglite")
library(svglite)
ggsave("biplot_PC1_PC2.svg", biplot_PCA_12, width = 4, height = 6)

# biplot that overlays PC2 and PC3
ggbiplot(static_microglia_non_red_pca,
         choices = c(2, 3), # specifying PCs 2 and 3
         groups = static_microglia_non_red$Treatment,
         ellipse = TRUE,
         var.axes = TRUE,
         varname.size = 2,
         varname.adjust = 2,
         alpha = 0.4,
         repel = TRUE, 
         vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
  scale_fill_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC2 and PC3",
       color = "Treatment", fill = "Treatment") +
  coord_fixed()

# out of curiosity, creating a biplot of PC3 vs PC4 
ggbiplot(static_microglia_non_red_pca,
         choices = c(3, 4), # specifying PCs 3 and 4
         groups = static_microglia_non_red$Treatment,
         ellipse = TRUE,
         var.axes = TRUE,
         varname.size = 2,
         varname.adjust = 2,
         alpha = 0.4,
         repel = TRUE, 
         vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
  scale_fill_manual(values = c(CON = "darkorange", COLD = "cornflowerblue")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC3 and PC4",
       color = "Treatment", fill = "Treatment") +
  coord_fixed()

# now I can also create a 3D scatterplot to look at the spread of PCs 1-3
## i will do this using the plotly package 
library(plotly)

# identifying the scores for the first 3 PCs, plus the grouping variable (treatment)
scores_non_red <- as.data.frame(static_microglia_non_red_pca$x[, 1:3])
scores_non_red$Treatment <- static_microglia_filtered$Treatment

# generating the variance explained, for the axis labels
ve <- round(summary(static_microglia_non_red_pca)$importance[2, ] * 100, 1)

# creating the 3D scatter plot 
fig_non_red <- plot_ly(scores_non_red,
               x = ~PC1, y = ~PC2, z = ~PC3,
               color = ~Treatment,
               colors = c(CON = "darkorange", COLD = "cornflowerblue"),
               type = "scatter3d",
               mode = "markers",
               marker = list(size = 1, opacity = 0.8))

fig_non_red <- fig_non_red %>% layout(
  title = "Microglia Morphology PCA (PC1-PC3)",
  scene = list(
    xaxis = list(title = paste0("PC1 (", ve[1], "%)")),
    yaxis = list(title = paste0("PC2 (", ve[2], "%)")),
    zaxis = list(title = paste0("PC3 (", ve[3], "%)"))
  )
)
fig_non_red

# ------------------------------------------------------------------------------
# CONDUCTING K-MEANS CLUSTERING 

# I am going to keep PCs 1-3 for this analysis, as they represent the majority of the variance 

# setting the seed to keep my results reproducible each time I run the code 
set.seed(123) # I need to do this because I will be taking random subsamples 

# determining the optimal number of clusters

# first I am going to use the elbow method  
### here I am taking a random subsample of 10k cells from my dataset to run this on because my computer does not have enough memory to run on the whole dataset
static_microglia_subsample <- sample(nrow(static_microglia_pc1_3), 
                                     min(10000, nrow(static_microglia_pc1_3)), 
                                     replace = FALSE)
fviz_nbclust(static_microglia_pc1_3[static_microglia_subsample,], FUNcluster = kmeans, 
             method = "wss", k.max = 10, nstart = 25, iter.max = 200) +
             labs(subtitle = "Elbow Method")
# this method suggests that either k = 2 or k = 4 is optimal (those are where slight plateaus occur)

# now I am going to use the silhouette method 
fviz_nbclust(static_microglia_pc1_3[static_microglia_subsample, ], kmeans, 
             method = "silhouette", k.max = 10, nstart = 25, iter.max = 200) +
             labs(subtitle = "Silhouette Method")
# this method suggests that k = 2 is optimal 

# these two methods likely suggest that fewer clusters are better because the points in my dataset were in one tightly packed cloud after PCA 
## this is because the morphologies in my dataset are very continuos, rather than discrete/dinstinct groups
### for this reason, I will test several k-values ranging from 2-6

# setting k = 3 to start 
## running kmeans() on my entire dataset of PC1-PC3
static_microglia_kmeans_3 <- kmeans(static_microglia_pc1_3, centers = 3, nstart = 25)
## printing out the values obtained from this computation 
### size = number of cells per cluster
str(static_microglia_kmeans_3)

# visualizing clusters using fviz_cluster()
fviz_cluster(static_microglia_kmeans_3, data = static_microglia_pc1_3)

# calculating the averages of each variable by cluster to try to interpret what each cluster means 
## creating a copy of my df and adding k=3 cluster column 
### making a copy of the df
static_microglia_non_red_k3 <- static_microglia_non_red
static_microglia_non_red_k3$Cluster <- factor(static_microglia_kmeans_3$cluster)
View(static_microglia_non_red_k3) # it worked!

## parsing out the numeric columns from the new df using sapply() again
static_microglia_non_red_k3_num_cols <- static_microglia_non_red_k3[sapply(static_microglia_non_red_k3, is.numeric)]

# averaging the values for each variable from each cluster (note: I am taking the means of the SCALED values)
## this step divides each clusters and takes the mean value of each variable for that specific cluster 
static_microglia_cluster_means <- aggregate(static_microglia_non_red_k3_num_cols,
                           by = list(Cluster = static_microglia_non_red_k3$Cluster),
                           FUN = mean)
static_microglia_cluster_means[-1] <- round(static_microglia_cluster_means[-1], 2)
# printing a table with the mean values for each variable from each cluster
static_microglia_cluster_means
# since the values were z-standardized, each value is in standard deviation units relative to all cells
## if the value is positive, it is above average, negative is below average, 0 = exactly average 

# visualizing this as a heatmap for convenience 
#install.packages('pheatmap')
library(pheatmap)
pheatmap(t(static_microglia_cluster_means[, -1]), labels_col = static_microglia_cluster_means$Cluster,
         cluster_cols = FALSE, scale = "none")

# based on these results, each cluster represents the following morphological group: 
## cluster 1: ameboid (highest circularity, low branch length, area, number of branches, etc.)
## cluster 2: rod-like (low circularity, high branch length, ratio of hull axes)
## cluster 3: ramified (high number of branches, triple points, area)

# plotting clusters on biplot of PC1 and PC2
ggbiplot(static_microglia_non_red_pca,
                          groups = static_microglia_non_red_k3$Cluster,
                          #ellipse = TRUE,
                          var.axes = TRUE,
                          varname.size = 2,
                          varname.adjust = 2,
                          alpha = 0.4,
                          repel = TRUE, 
                          vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c('1' = "darkorange", '2' = "cornflowerblue", '3' = "mediumslateblue")) +
  scale_fill_manual(values = c('1' = "darkorange", '2' = "cornflowerblue", '3' = "mediumslateblue")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC1 and PC2",
       color = "Cluster", fill = "Cluster") +
  coord_fixed()

# now I will perform statistical tests to evaluate the differences between clusters
## I will examine how similar each cluster is based on different variables 
### I will choose variables that have large loadings and correlate with a certain cluster 

# calculating the pca loadings to have an empirical measure when choosing variabels 
pca_loadings <- static_microglia_non_red_pca$rotation %*% diag(static_microglia_non_red_pca$sdev)
pca_loadings


## testing variable: FOREGROUND PIXELS - large loading magnitude, aligns with cluster 1
# I am using kruskall wallis, which is the non-parametric version of anova
kruskal.test(Foreground.pixels ~ Cluster, data = static_microglia_non_red_k3)
# results: chi-squared = 16811, df = 2, p-value < 2.2e-16
# indicates that there is a significant difference between clusters, so I will follow up with a post-hoc test 
# post-hoc test: wilcox/mann-whitney u-test 
pairwise.wilcox.test(static_microglia_non_red_k3$Foreground.pixels, 
                     static_microglia_non_red_k3$Cluster, 
                     p.adjust.method = "BH") # Bonferroni correction for multiple testing
# results: p<2e-16 for each comparison of two clusters 
# indicates that each cluster has significant differences in foreground pixels 

## testing variable: MAX/MIN RATIO OF HULL AXES (measure of cell shape - circularity/length) - large loading magnitude, aligns with cluster 2
# I am using kruskall wallis, which is the non-parametric version of anova
kruskal.test(Max.min.hull.radius  ~ Cluster, data = static_microglia_non_red_k3)
# results: chi-squared = 14532, df = 2, p-value < 2.2e-16
# indicates that there is a significant difference between clusters, so I will follow up with a post-hoc test 
# post-hoc test: wilcox/mann-whitney u-test 
pairwise.wilcox.test(static_microglia_non_red_k3$Max.min.hull.radius, 
                     static_microglia_non_red_k3$Cluster, 
                     p.adjust.method = "BH") # Bonferroni correction for multiple testing
# results: p<2e-16 for each comparison of two clusters 
# indicates that each cluster has significant differences in cell max/min hull radius (measure of cell shape - circularity/length)

## testing variable: CIRCULARITY - large loading magnitude, aligns with cluster 3
# I am using kruskall wallis, which is the non-parametric version of anova
kruskal.test(Circularity ~ Cluster, data = static_microglia_non_red_k3)
# results: chi-squared = 16076, df = 2, p-value < 2.2e-16
# indicates that there is a significant difference between clusters, so I will follow up with a post-hoc test 
# post-hoc test: wilcox/mann-whitney u-test 
pairwise.wilcox.test(static_microglia_non_red_k3$Circularity, 
                     static_microglia_non_red_k3$Cluster, 
                     p.adjust.method = "BH") # Bonferroni correction for multiple testing
# results: p<2e-16 for each comparison of two clusters 
# indicates that each cluster has significant differences in cell circularity

# after running these stats, I can confirm that my data 

# ------------------------------------------------------------------------------
# CONDUCTING K-MEANS CLUSTERING (k = 4, 5, 6, and 2)

# setting k = 4
static_microglia_kmeans_4 <- kmeans(static_microglia_pc1_3, centers = 4, nstart = 25)
fviz_cluster(static_microglia_kmeans_4, data = static_microglia_pc1_3)
# creating new df with added cluster column 
static_microglia_non_red_k4 <- static_microglia_non_red
static_microglia_non_red_k4$Cluster <- factor(static_microglia_kmeans_4$cluster)
# plotting clusters on biplot of PC1 and PC2
ggbiplot(static_microglia_non_red_pca,
         groups = static_microglia_non_red_k4$Cluster,
         #ellipse = TRUE,
         var.axes = TRUE,
         varname.size = 2,
         varname.adjust = 2,
         alpha = 0.4,
         repel = TRUE, 
         vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c('1' = "darkorange", '2' = "cornflowerblue", '3' = "mediumslateblue", '4' = "darkgreen")) +
  scale_fill_manual(values = c('1' = "darkorange", '2' = "cornflowerblue", '3' = "mediumslateblue", '4' = "darkgreen")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC1 and PC2",
       color = "Cluster", fill = "Cluster") +
  coord_fixed()


# setting k = 5 
static_microglia_kmeans_5 <- kmeans(static_microglia_pc1_3, centers = 5, nstart = 25)
fviz_cluster(static_microglia_kmeans_5, data = static_microglia_pc1_3)
# creating new df with added cluster column 
static_microglia_non_red_k5 <- static_microglia_non_red
static_microglia_non_red_k5$Cluster <- factor(static_microglia_kmeans_5$cluster)
# plotting clusters on biplot of PC1 and PC2
ggbiplot(static_microglia_non_red_pca,
         groups = static_microglia_non_red_k5$Cluster,
         #ellipse = TRUE,
         var.axes = TRUE,
         varname.size = 2,
         varname.adjust = 2,
         alpha = 0.4,
         repel = TRUE, 
         vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c('1' = "darkorange", '2' = "cornflowerblue", '3' = "mediumslateblue", '4' = "darkgreen", '5' = "magenta")) +
  scale_fill_manual(values = c('1' = "darkorange", '2' = "cornflowerblue", '3' = "mediumslateblue", '4' = "darkgreen", '5' = "magenta")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC1 and PC2",
       color = "Cluster", fill = "Cluster") +
  coord_fixed()

# setting k = 6
static_microglia_kmeans_6 <- kmeans(static_microglia_pc1_3, centers = 6, nstart = 25)
fviz_cluster(static_microglia_kmeans_6, data = static_microglia_pc1_3)
# creating new df with added cluster column 
static_microglia_non_red_k6 <- static_microglia_non_red
static_microglia_non_red_k6$Cluster <- factor(static_microglia_kmeans_6$cluster)
# plotting clusters on biplot of PC1 and PC2
ggbiplot(static_microglia_non_red_pca,
         groups = static_microglia_non_red_k6$Cluster,
         #ellipse = TRUE,
         var.axes = TRUE,
         varname.size = 2,
         varname.adjust = 2,
         alpha = 0.4,
         repel = TRUE, 
         vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c('1' = "darkorange", '2' = "cornflowerblue", '3' = "mediumslateblue", '4' = "darkgreen", '5' = "magenta", '6' = "yellow")) +
  scale_fill_manual(values = c('1' = "darkorange", '2' = "cornflowerblue", '3' = "mediumslateblue", '4' = "darkgreen", '5' = "magenta", '6' = "yellow")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC1 and PC2",
       color = "Cluster", fill = "Cluster") +
  coord_fixed()

# setting k = 2
static_microglia_kmeans_2 <- kmeans(static_microglia_pc1_3, centers = 2, nstart = 25)
fviz_cluster(static_microglia_kmeans_2, data = static_microglia_pc1_3)
# creating new df with added cluster column 
static_microglia_non_red_k2 <- static_microglia_non_red
static_microglia_non_red_k2$Cluster <- factor(static_microglia_kmeans_2$cluster)
# plotting clusters on biplot of PC1 and PC2
ggbiplot(static_microglia_non_red_pca,
         groups = static_microglia_non_red_k2$Cluster,
         #ellipse = TRUE,
         var.axes = TRUE,
         varname.size = 2,
         varname.adjust = 2,
         alpha = 0.4,
         repel = TRUE, 
         vector.args = list(linewidth = 0.3)) +
  scale_color_manual(values = c('1' = "darkorange", '2' = "cornflowerblue")) +
  scale_fill_manual(values = c('1' = "darkorange", '2' = "cornflowerblue")) +
  guides(shape = "none") +
  theme_minimal() +
  labs(title = "Biplot of PC1 and PC2",
       color = "Cluster", fill = "Cluster") +
  coord_fixed()

# ------------------------------------------------------------------------------
# RUNNING SILHOUETTE TEST AND DUNN INDEX FOR K-MEANS CLUSTERS

# these tests will help me determine which k value is optimal 

# starting with the average silhouette width method 
## first I am creating a list of all the k-means products I have created 
all_kmeans <- list(k2 = static_microglia_kmeans_2, k3 = static_microglia_kmeans_3,
                k4 = static_microglia_kmeans_4, k5 = static_microglia_kmeans_5,
                k6 = static_microglia_kmeans_6)
# now I am going to calculate the average silhouette widths (ASW) for each k-value 
static_microglia_pc1_3_dist <- dist(static_microglia_pc1_3)   # distances between ALL cells
sil_avg <- sapply(all_kmeans, function(m) {
  silhouette_width <- silhouette(m$cluster, static_microglia_pc1_3_dist) # this calculates the AWS
  mean(silhouette_width[, "sil_width"])
})
round(sil_avg, 3)

##### NOTE: I've included this section of commented out code because I had a lot of issues dealing with a large dataset and having enough memory to do so
# this code here uses a subsample of the objects within all_kmeans (see below) to work around this, however it is not ideal because it does not use the entire dataset
# however, I wanted to include it as a backup option in case anyone references this script and needs another method that will use less memory to run
#sil_avg <- sapply(all_kmeans, function(m) {
  # sample 3000 cells (without replacement) from each cluster for each object in all_kmeans
  # I need to use a subsample because my dataset is too large to run the analysis on all of my cells  
#  cluster_subsample <- unlist(lapply(split(seq_along(m$cluster), m$cluster),
#                                     function(i) i[sample.int(length(i), min(5000, length(i)))]))
#  d <- dist(static_microglia_pc1_3[cluster_subsample, ])
#  s <- silhouette(m$cluster[cluster_subsample], d)
#  mean(s[, "sil_width"])
#})
#round(sil_avg, 3)

# now I will plot these values as I have done before
sil_df <- data.frame(k = 2:6, sil = as.numeric(sil_avg))
best_k <- sil_df$k[which.max(sil_df$sil)] # adds line at the optimal k-value 
# using ggplot to recreate the style of silhouette plots in fviz_nbclust for the above analysis 
ggplot(sil_df, aes(x = k, y = sil, group = 1)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  geom_vline(xintercept = best_k, linetype = "dashed", color = "steelblue") +
  scale_x_continuous(breaks = sil_df$k) +
  labs(title = "Optimal number of clusters (method = \"silhouette\")",
       subtitle = "Silhouette Method",
       x = "Number of clusters k",
       y = "Average silhouette width") +
  theme_classic()
# k = 2 has the highest value, which means 2 clusters are the most optimal. ASW decreases as k-value increases 

# the second method I will use is to calculate the Dunn index
## the Dunn index assesses how compact and well-separated clusters are in cluster analysis
### at first I tried to use the clValid Dunn function, but this took too much memory, so I had an LLM write a faster function for me that uses the same math
fast_dunn <- function(x, cl_list, chunk = 500) {
  x <- as.matrix(x)
  sq <- rowSums(x^2)
  within_max  <- rep(0, length(cl_list))
  between_min <- rep(Inf, length(cl_list))
  for (s in seq(1, nrow(x), by = chunk)) {
    rows <- s:min(s + chunk - 1, nrow(x))
    # distances from this chunk of cells to ALL cells
    dmat <- sqrt(pmax(outer(sq[rows], sq, "+") - 2 * x[rows, , drop = FALSE] %*% t(x), 0))
    for (j in seq_along(cl_list)) {
      cl <- cl_list[[j]]
      same <- outer(cl[rows], cl, "==")
      within_max[j]  <- max(within_max[j],  max(dmat[same]))
      between_min[j] <- min(between_min[j], min(dmat[!same]))
    }
  }
  setNames(between_min / within_max, names(cl_list))
}
# printing out Dunn index values 
km_cl_list <- lapply(all_kmeans, function(m) m$cluster)
round(fast_dunn(static_microglia_pc1_3, km_cl_list), 4)
# these values are VERY small, which is to because the minimum distance between clusters is very small 
## I can still interpret them, but the differences between should be taken with a grain of salt, since they are also very slight/small
### the Dunn index suggests that k = 2 is optimal and k = 4 is the next best option


# ------------------------------------------------------------------------------
# CONDUCTING WARDD2 CLUSTERING 

# since the k-means algorithm is best-suited for round, spherical clusters and my data is very compact and groups aren't spread out, I wanted to attempt another method 
## my data are very spread out because they are from embryonic samples, where microglial morphologies are less distinct and more continuous between subtypes

# wardD2 is a hierarchical clustering method suggested by Margo that works well with data that is tightly packed together 
## website that I referenced (URL): https://rstudio-pubs-static.s3.amazonaws.com/599072_93cf94954aa64fc7a4b99ca524e5371c.html#Hierarchical_clustering

# I am going to be using wardD2 clustering from the fastcluster package, which will allow me to run the whole dataset on my computer 
static_microglia_pc1_3_wardd2 <- hclust.vector(as.matrix(static_microglia_pc1_3), method = "ward") #in this package, the method ward = wardD2 

# creating a dendrogram to represent the clustering data
plot(static_microglia_pc1_3_wardd2, cex = 0.4, labels = FALSE)

# now I need to determine the optimal number of clusters from wardD2 clustering 
## I will first use the 'vertical distance method'
### in this method, you identify the longest vertical line that isn't intersected by horizontal lines and the number of vertical lines parallel to it is the number of clusters 

# printing out the 
dendrogram_heights <- sort(static_microglia_pc1_3_wardd2$height, decreasing = TRUE)[1:10]

dendrogram_heights
# length of the vertical line that corresponds to each number of clusters k
gaps <- data.frame(k = 2:10,
                   gap = dendrogram_heights[-length(dendrogram_heights)] - dendrogram_heights[-1])
gaps

View(gaps)
# k with the longest vertical line
gaps$k[which.max(gaps$gap)]

barplot(gaps$gap, names.arg = gaps$k,
        xlab = "Number of clusters k", ylab = "Length of vertical line (height gap)")

#### in this case, it is ___ clusters


# now I am going to conduct average silhouette width analysis 
## creating a list of wardD2 clustering outputs for k-values from 2-6 
k_values2_6 <- 2:6
ward_list <- lapply(k_values2_6, function(k) cutree(static_microglia_pc1_3_wardd2, k = k))
names(ward_list) <- paste0("k", k_values2_6)
# using the silhouette function to run analysis 
sil_avg_wardD2 <- sapply(ward_list, function(cl) mean(silhouette(cl, static_microglia_pc1_3_dist)[, "sil_width"]))
round(sil_avg, 3)
# plotting average silhouette width for wardD2 clustering 
sil_df_wardD2 <- data.frame(k = 2:6, sil = as.numeric(sil_avg_wardD2))
best_k_wardD2 <- sil_df_wardD2$k[which.max(sil_df_wardD2$sil)] # adds line at the optimal k-value 
# using ggplot to recreate the style of silhouette plots in fviz_nbclust for the above analysis 
ggplot(sil_df_wardD2, aes(x = k, y = sil, group = 1)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  geom_vline(xintercept = best_k_wardD2, linetype = "dashed", color = "steelblue") +
  scale_x_continuous(breaks = sil_df_wardD2$k) +
  labs(title = "Optimal number of clusters (method = \"silhouette\")",
       subtitle = "Silhouette Method",
       x = "Number of clusters k",
       y = "Average silhouette width") +
  theme_classic()
# based on this analysis, k = 2 is the best number of clusters, followed by k = 4 

# I will also calculate the Dunn index for each number of wardD2 clusters 
## here I am using the function created when calculating Dunn index for K-means clustering, but subbing in the wardD2 data
ward_dunn <- fast_dunn(static_microglia_pc1_3, ward_list)
round(ward_dunn, 4)
# this analysis suggests that ____

# after conducting the above analyses, I can conclude that k = 2 and k = 4 are the best k-values 
## k = 2 is optimal, however, it likely has less biological relevance than k = 4 
### to determine which is more relevant in terms of morphological features from my dataset, I am going to create heatmaps for clustering with both k-values and see which features contribute to them

# starting with wardD2 clustering with k = 2 
static_microglia_non_red_ward_k2 <- static_microglia_non_red
static_microglia_non_red_ward_k2$Cluster <- factor(cutree(static_microglia_pc1_3_wardd2, k = 2))
# parsing out the numeric columns (Cluster is a factor, so it is left out)
static_microglia_non_red_ward_k2_num_cols <- static_microglia_non_red_ward_k2[sapply(static_microglia_non_red_ward_k2, is.numeric)]
# averaging the (scaled) values of each variable within each cluster
static_microglia_ward_cluster_means_k2 <- aggregate(static_microglia_non_red_ward_k2_num_cols,
                                                 by = list(Cluster = static_microglia_non_red_ward_k2$Cluster),
                                                 FUN = mean)
static_microglia_ward_cluster_means_k2[-1] <- round(static_microglia_ward_cluster_means_k2[-1], 2)
static_microglia_ward_cluster_means_k2
# creating my heatmap
library(pheatmap)
pheatmap(t(static_microglia_ward_cluster_means_k2[, -1]),
         labels_col = static_microglia_ward_cluster_means_k2$Cluster,
         cluster_cols = FALSE, scale = "none")
# for k = 2: 
## cluster 1 cells have a smaller area and perimeter, fewer branches, a smaller radius, low circularity, low max branch length
## cluster 2 cells have a larger area and perimeter, more branches, a larger radius, low circularity, relatively low max branch length
### even though this k = 2 is the better cluster number via empirical methods, it makes less sense when looking at the implicated features 
#### cluster 1 cells have low area and less branches that cluster 2, this likely puts cells that are more ameboid and rod-like in the same cluster (1)

# moving on to wardD2 clustering with k = 4
## making a copy of the df and adding the Ward cluster column 
static_microglia_non_red_ward_k4 <- static_microglia_non_red
static_microglia_non_red_ward_k4$Cluster <- factor(cutree(static_microglia_pc1_3_wardd2, k = 4))
# parsing out the numeric columns 
static_microglia_non_red_ward_k4_num_cols <- static_microglia_non_red_ward_k4[sapply(static_microglia_non_red_ward_k4, is.numeric)]
# averaging the (scaled) values of each variable within each cluster
static_microglia_ward_cluster_means_k4 <- aggregate(static_microglia_non_red_ward_k4_num_cols,
                                                 by = list(Cluster = static_microglia_non_red_ward_k4$Cluster),
                                                 FUN = mean)
static_microglia_ward_cluster_means_k4[-1] <- round(static_microglia_ward_cluster_means_k4[-1], 2)
static_microglia_ward_cluster_means_k4
# creating the heatmap
library(pheatmap)
pheatmap(t(static_microglia_ward_cluster_means_k4[, -1]),
         labels_col = static_microglia_ward_cluster_means_k4$Cluster,
         cluster_cols = FALSE, scale = "none")
# for k = 4: 
## cluster 1: average area, perimeter, number of branches and circularity, low hull radius
## cluster 2: small area, perimeter, few branches, highest circularity (of the 4), many hull pixels, low max branch length
## cluster 3: high max branch length, ratio of hull axes and hull radius, low circularity and hull pixels 
## cluster 4: large area, perimeter, many branches, low circularity, low mean branch length 






# ------------------------------------------------------------------------------
# UMAP DIMENSIONALITY REDUCTION

# this is another method of dimensionality reduction that I want to use 






# ------------------------------------------------------------------------------


# from margo to extend running memory:
#options(future.globals.maxSize = 200000 * 1024^2)

