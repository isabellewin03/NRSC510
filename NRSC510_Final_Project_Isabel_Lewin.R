# NRSC510 final project in R 

# NRSC510 class project 
## Isabel Lewin 

# installing new packages that I will need for PCA + data visualization
##install.packages("factoextra")
##install.packages("MASS")

# libraries 
library(vroom)
library(dplyr)
library(ggplot2)
library(factoextra)
library(MASS)

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
# the majority of micrglia characteristics are double-precision floating-point number <dbl>

# using the summary function to examine characteristics of each variable more closely 
summary(static_microglia)

# removing the 'cluster' column because this was created by a previous student in my lab after conducting K-means clustering 
## the 'cluster' column was column #35 so I used the indexing brackets to remove it and replaced the old dataset with a new one that doesn't contain #35
### I also removed any row from the BrainRegion column that were not HYPO (hypothalamus)
#### this is because the student who generated these data also looked at specific hypothalamic nuclei, which I don't want to examine
static_microglia <- static_microglia[,-35]%>%filter(BrainRegion == 'HYPO')

#static_microglia_edited <- static_microglia%>%filter(BrainRegion == 'HYPO')
# checking dataframe and it looks good!
View(static_microglia) # now my dataframe only contains 
glimpse(static_microglia) # now my dataframe only has 34 columns 


# i am going to use PCA for the first part of my analysis for this dataset, and before i do that, i need to standardize my data
## this will prevent any variables with large values from skewing the PCA results

# my first step is to prevent the categorical variables from being standardized: 
## to do this i am using the sapply() function to isolate the numeric columns 
static_microglia_num_cols <- sapply(static_microglia, is.numeric)

# now I will use the lapply() function to standardize the numeric columns in my dataset 
## i will perform z-standardization, using the following operation: z = (value - mean)/standard deviation, where z is the new, standardized value 
### the scale() function automatically performs z-scoring so I will use it for this case
static_microglia[static_microglia_num_cols] <- lapply(static_microglia[static_microglia_num_cols], function(x) as.numeric(scale(x)))


# checking to make sure my standardized dataset looks correct
View(static_microglia)
# looks good! all numeric variables have been z-scored

# now I will run a PCA using the function prcomp()
## here, I excluded columns 1-7 from analysis because they are categorical variables and I don't want them in my analysis
static_microglia_pca <- prcomp(static_microglia[,-c(1:7)], 
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

# creating a biplot 
fviz_pca_biplot(static_microglia_pca, 
                label = 'var', 
                habillage = static_microglia$Sex,
                title = "Microglia Morphology PCA Biplot")






# tips from margo:

# if data looks odd:
# remove outliers?
# make sure to use variables that are unique -- no summaries of the same variable (ie. mean and std)
# check to see if log normalization or another norm method might help remove noise 
# try visualizing pc1-4 to check how many pcs you have 
# if more than 3, could use umap 



