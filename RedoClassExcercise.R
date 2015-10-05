#clear all previous variable/objects
rm(list=ls())

#Set the working directory to my class folder which contains the downloaded files
setwd('C:/Users/Emil/Documents/Bio720')

#Data file info:
#M172 and M200 are sm
#M180 and M257 are lg

#Read in the data we had before into 4 tables names M172.data for file M172*.xprs, M180.data for file M180*.xprs, etc.
M172.data <- read.table("M172_sm_male_hdhorn_ATCACG_L003_results.xprs", header=TRUE)
M180.data <- read.table("M180_lg_male_hdhorn_CAGATC_L001_results.xprs", header=TRUE)
M200.data <- read.table("M200_sm_male_hdhorn_ACAGTG_L008_results.xprs", header=TRUE)
M257.data <- read.table("M257_lg_male_hdhorn_ATGTCA_L002_results.xprs", header=TRUE)

#View dimmensions of all tables
dim(M172.data)
dim(M180.data)
dim(M200.data)
dim(M257.data)

#Look at the first 5 targer_id entries of each table to see original order (can then compare after sort to make sure table is sorted)
M172.data[1:5,2]
M180.data[1:5,2]
M200.data[1:5,2]
M257.data[1:5,2]

#Sort the tables based on column 2 ("target_id") and put into new tables called s*.data (for sorted[Data_ID].data))
sM172.data <- M172.data[order(M172.data$target_id),]
sM180.data <- M180.data[order(M180.data$target_id),]
sM200.data <- M200.data[order(M200.data$target_id),]
sM257.data <- M257.data[order(M257.data$target_id),]

#Look at the first five entries in column "target_id" to see that the tables are now properly sorted
sM172.data[1:5,2]
sM180.data[1:5,2]
sM200.data[1:5,2]
sM257.data[1:5,2]

#Combine the TPM measures (column 15) of each data set into a new table called comb_data
comb_data <- cbind(sM172.data[,15], sM180.data[,15], sM200.data[,15], sM257.data[,15])

#Give the columns and rows the right names
colnames(comb_data) <- c("M172_sm", "M180_lg", "M200_sm", "M257_lg")
#rownames(comb_data) <- sM172.data[,2] Commented this out cause it looks like garbage

#Look at the top few lines of comb_data to make sure all the information is there
head(comb_data)

#Also check the dimension to make sure they make sense
dim(comb_data)

#Look at the correlations between measures of TPM and store it in a table called Ccomb_data
Ccomb_data <- cor(comb_data)

#Production of scatterplots of all data files against all data files

#The below example will plot column 1 of comb_data (M172_sm) versus column 2 of comb_data (M180_lg).
#This could then also be repeated for column 1 (C1) vs column 3 (C3), C1 vs C4, C2 vs C3, C2 vs C4, C3 vs C4
#Doing it this way will produce seperate scatterplots where the first column is always along the X-axis
#plot(comb_data[,1], comb_data[,2])

#The scatterplots can also be visualized using pairs(comb_data). 
#This automatically plots every column against every other column in a single figure. Data sets appear on the X-axis in their row and on the y-axis in their column.
pairs(comb_data)
pairs(Ccomb_data)

#Make a vector containing all the row means for comb_data
row_means <- rowMeans(comb_data)

#Add the row_means to the comb_data matrix, calling the new matrix gene_means. Add gene names as well.
gene_means <- cbind(comb_data, row_means)
rownames(gene_means) <- sM172.data$target_id

#Check the dimensions of gene_means to make sure everything makes sense
dim(gene_means)

#Create a subset of genes which have a row means greater than 10
G10_means <- subset(gene_means,row_means>10)

#-----------------------Making an MA Plot-----------------------
#MA plots (microarray plots?) show differences in gene expression between a control group and a test group
#Check the first few minutes of this for a reminder: https://www.youtube.com/watch?v=46-t2jOYsyY
#The plot is M vs A where:
#A is log2 average of (xy); ie. .5*(log2(y)+log2(x))
#M is y-x; ie. M=log2(y)-log2(x)
#When doing an MA plot on biological replications should first average them

#Removing row names from G10_means for clarity's sake
NN_means <- G10_means
rownames(NN_means) <- NULL

#Removing row_means as its useless at this point
NN_means <- NN_means[,-5]

#Looping through NN_means to average sm and lg observations
sm_means <- NULL
lg_means <- NULL
for (n in 1:nrow(NN_means))
{
  x <- 0
  y <- 0
  x <- (NN_means[n,1] + NN_means[n,3])/2
  y <- (NN_means[n,2] + NN_means[n,4])/2
  sm_means <- append(sm_means,x)
  lg_means <- append(lg_means,y)
}

#Converts from whatever stores named numbers to just numbers. Makes the output a little nicer and easier to read.
sm_means <- as.vector(sm_means)
lg_means <- as.vector(lg_means)

#Make a matrix containing sm_means and lg_means
SMvLG <- cbind(sm_means, lg_means)

#Calculating the means of the logbase2 sm and lg means (ie. A)
A <- rowMeans(log2(SMvLG))

#Calculating the difference in the log base 2 values for lg vs sm
M <- NULL
for (i in 1:nrow(SMvLG))
{
 x <- NULL
 x <- (log2(lg_means[i])-log2(sm_means[i]))
 M <- append(M,x)
}

#Making my MA plot
plot(A,M,main="LG versus SM Data")