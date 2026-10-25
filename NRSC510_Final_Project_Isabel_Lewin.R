# NRSC510 final project in R 

# NRSC510 class project 
## Isabel Lewin 

# installing new packages that I will need for PCA + data visualization
#install.packages("factoextra") # for PCA 
#install.packages("MASS")
#install.packages("plotly") # for 3D scatter plot 
#install.packages("ggbiplot") # for creating PCA biplots using ggplot

# libraries 
library(vroom)
library(dplyr)
library(ggplot2)
library(factoextra)
library(MASS)
library(ggbiplot)

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
static_microglia_filtered[static_microglia_num_cols] <- lapply(static_microglia_filtered[static_microglia_num_cols], function(x) as.numeric(scale(x)))

# checking to make sure my standardized dataset looks correct
View(static_microglia_filtered)
# looks good! all numeric variables have been standardized by z-scoring

# now I will run a PCA using the function prcomp()
## here, I excluded columns 1-7 from analysis because they are categorical variables and I don't want them in my analysis
static_microglia_pca <- prcomp(static_microglia_filtered[,-c(1:7)], 
                               scale = TRUE)

# checking out my data to make sure it looks right after doing PCA 
summary(static_microglia_pca)
# looks good!

# this is going to print out the names of more data saved in static_microglia_pca
names(static_microglia_pca)
# looking at each one individually: 
static_microglia_pca$sdev # this is the standard deviation of the components 
static_microglia_pca$rotation # these are the Eigenvectors 
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
component_matrix_PC123 <- static_microglia_pca$rotation[,-c(4:27)] # these are the Eigenvectors 
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
                                   scale = TRUE)

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
ggbiplot(static_microglia_non_red_pca,
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
  labs(title = "Biplot of PC1 and PC2",
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
  labs(title = "Biplot of PC1 and PC2",
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






# plots 
# PC1 vs 2
# colour scatterplots, make it heatmap 
# ggplot 


# tips from margo:

# if data looks odd:
# remove outliers?
# make sure to use variables that are unique -- no summaries of the same variable (ie. mean and std)
# check to see if log normalization or another norm method might help remove noise 
# try visualizing pc1-4 to check how many pcs you have 
# if more than 3, could use umap 



