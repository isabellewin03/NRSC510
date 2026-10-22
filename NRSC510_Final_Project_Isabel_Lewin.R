# NRSC510 final project in R 

# NRSC510 class project 
## Isabel Lewin 

# installing new packages that I will need for PCA + data visualization
##install.packages("factoextra") # for PCA 
##install.packages("MASS")
#install.packages("plotly") # for 3D scatter plot 

# libraries 
library(vroom)
library(dplyr)
library(ggplot2)
library(factoextra)
library(MASS)

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
  rename(Foreground.hull.pixel.density = Density.of.foreground.pixels.in.hull.area,
         Ratio.of.hull.axes = Span.ratio.of.hull..major.minor.axis., 
         Max.span.across.hull = Maximum.span.across.hull,
         Max.radius.from.hull.center = Maximum.radius.from.hull.s.center.of.mass,
         Max.min.radius.from.hull.center = Max.min.radii.from.hull.s.center.of.mass, 
         CV.in.radii.from.hull.center = Relative.variation..CV..in.radii.from.hull.s.center.of.mass, 
         Max.radius.from.circle.center = Maximum.radius.from.circle.s.center.of.mass,
         Max.min.radius.from.circle.center = Max.min.radii.from.circle.s.center.of.mass,
         CV.in.radii.from.circle.center = Relative.variation..CV..in.radii.from.circle.s.center.of.mass,
         Mean.radius.from.circle.center = Mean.radius.from.circle.s.center.of.mass, 
         Number.of.branches = X..of.branches, 
         Number.of.junctions = X..of.junctions, 
         Number.of.end.point.voxels = X..of.end.point.voxels,
         Number.of.junction.voxels = X..of.junction.voxels,
         Number.of.slab.voxels = X..of.slab.voxels, 
         Number.of.triple.points = X..of.triple.points, 
         Number.of.quadruple.points = X..of.quadruple.points)

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
# add new column for cumulative explained vairance -- sum

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
## this biplot only shows the top 8 contributing variables for simplification
pca_12 <- fviz_pca_biplot(static_microglia_pca,
                geom = c("point"),
                label = "var",
                habillage = static_microglia_filtered$Treatment,
                col.var = "black",
                labelsize = 4,
                repel = TRUE,
                ggrepel.max.overlaps = Inf,
                alpha.var = "contrib", 
                alpha.legend = "none",
                select.var = list(contrib = 8),   # top 8 contributors only
                title = "Microglia Morphology PCA Biplot",
                xlab = "PC1 (55.3%)",
                ylab = "PC2 (20.0%)",
                ggtheme = theme_minimal())
pca_12 <- pca_12 + guides(alpha = "none")
pca_12

# # creating a biplot that overlays PC2 and PC3, again visualizing only the top 8 contributing variables 
pca_23 <- fviz_pca_biplot(static_microglia_pca,
                axes = c(2, 3),   
                geom = c("point"),
                label = "var",
                habillage = static_microglia_filtered$Treatment,
                col.var = "black",
                labelsize = 4,
                repel = TRUE,
                ggrepel.max.overlaps = Inf,
                alpha.var = "contrib", 
                select.var = list(contrib = 8),   # top 8 contributors only
                title = "Microglia Morphology PCA Biplot",
                xlab = "PC1 (55.3%)",
                ylab = "PC2 (20.0%)",
                ggtheme = theme_minimal())
pca_23 <- pca_23 + guides(alpha = "none")
pca_23

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

# now I am going to remove any redundant variables from my analysis (ie. variables that are summary statistics of one another)
## in addition to examining the names of the variables, I am going to check empirically using cor() to see which variables correlate with each other 
library(corrplot)

# isolating numbered columns from my data
static_microglia_num_cols <- sapply(static_microglia, is.numeric)

# using the cor() function to see which variables correlate with each other 
cor_mat <- cor(static_microglia[static_microglia_num_cols], use = "pairwise.complete.obs")
corrplot(cor_mat, method = "color", type = "upper", tl.cex = 0.6) # plotting this 

# listing variable pairs with |r| > 0.9
high <- which(abs(cor_mat) > 0.95 & upper.tri(cor_mat), arr.ind = TRUE)
data.frame(var1 = rownames(cor_mat)[high[, 1]],
           var2 = colnames(cor_mat)[high[, 2]],
           r = cor_mat[high])

print(colnames(static_microglia))

# using a cluster dendrogram to visualize the different blocks of variables that are redundant 
## blocks of variables below the red line (0.1) 
dist_mat <- as.dist(1 - abs(cor_mat))
hc <- hclust(dist_mat, method = "complete")
plot(hc, cex = 0.7)
abline(h = 0.05, col = "red", lty = 2)   # cut at |r| = 0.9

groups <- cutree(hc, h = 0.05)
split(names(groups), groups) # this prints my variable names listed by block 

keep <- c("Foreground.pixels", "Density.of.foreground.pixels.in.hull.area", "Span.ratio.of.hull..major.minor.axis.", 
          "Perimeter", "Area", "Circularity", "Width.of.bounding.rectangle", "Height.of.bounding.rectangle", 
          "Max.min.radii.from.hull.s.center.of.mass", "Relative.variation..CV..in.radii.from.hull.s.center.of.mass", 
          "Max.min.radii.from.circle.s.center.of.mass", "Relative.variation..CV..in.radii.from.circle.s.center.of.mass", 
          "X..of.branches", "X..of.end.point.voxels", "X..of.slab.voxels", "Average.branch.length", "X..of.quadruple.points", 
          "Maximum.branch.length")         
cor_new <- cor(static_microglia[keep], use = "pairwise.complete.obs")
max(abs(cor_new[upper.tri(cor_new)]))     # should now be < 0.9


static_microglia_redundancy <- static_microglia[keep]
static_microglia_pca_red <- prcomp(static_microglia_redundancy[,-c(1:7)], 
                                   scale = TRUE)
fviz_eig(static_microglia_pca_red, 
         addlabels = TRUE, 
         ylim = c(0, 60),
         barfill = 'skyblue',
         barcolor = 'black',
         main = "Microglia Morphology PCA Scree Plot", 
         ggtheme = theme_classic())
fviz_pca_biplot(static_microglia_pca_red, 
                label = 'var', 
                habillage = static_microglia$Sex,
                col.var = "black", 
                repel = TRUE,
                title = "Microglia Morphology PCA Biplot", 
                xlab = "PC1 (55.3%)", 
                ylab = "PC2 (20.0%)", 
                ggtheme = theme_minimal())


#static_microglia_edited <- static_microglia%>%filter(BrainRegion == 'HYPO')
# checking dataframe and it looks good!
View(static_microglia) # now my dataframe only contains 
glimpse(static_microglia) # now my dataframe only has 34 columns 




#----------------------------
# can't get this to work - is it worth it?


# try plotly for 3d scatter -- 



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



