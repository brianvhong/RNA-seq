
## Library
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

################## Introduction and Data Import 
# Read the data into R
seqdata <- read.delim("./data/GSE60450_LactationGenewiseCounts.txt", stringsAsFactors = FALSE)

# Read the sample information into R
sampleinfo <- read.delim("data/SampleInfo.txt", stringsAsFactors = TRUE)

head(seqdata)

################# Format data 
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]

colnames(countdata)

# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
head(countdata)

# Check to see if column names == sample name
table(colnames(countdata)==sampleinfo$SampleName)


################## Convert counts to DGEList object
y <- DGEList(countdata) ## edgeR package to store object
# have a look at y
y

# See what slots are stored in y
names(y)

# Library size information is stored in the samples slot
y$samples

# We can also store the groups for the samples in the DGEList object.
group <- paste(sampleinfo$CellType, sampleinfo$Status,sep=".")
# Take a look
group

# Convert to factor
group <- factor(group)
# Take another look.
group

# Add the group information into the DGEList
y$samples$group <- group
y$samples

################## Adding Annotation (We will demonstrate how to do this using the org.Mm.eg.db package)

# First we need to decide what information we want. In order to see what we can extract we can run the columns function on the annotation database.
columns(org.Mm.eg.db)

# We definitely want gene symbols and perhaps the full gene name. Let’s build up our annotation information in a separate data frame using the select function.
ann <- select(org.Mm.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))

# Have a look at the annotation
head(ann)

# Let’s double check that the ENTREZID column matches exactly to our y$counts rownames.
table(ann$ENTREZID==rownames(y$counts))

#We can slot in the annotation information into the genes slot of y. (Please note that if the select function returns a 1:many mapping then you can’t just append the annotation to the y object. An alternative way to get annotation will be discussed during the analysis of the second dataset.)
y$genes <- ann

################## Filtering lowly expressed genes
#  In this dataset, we choose to retain genes if they are expressed at a counts-per-million (CPM) above 0.5 in at least two samples.
# We’ll use the cpm function from the edgeR library (M. D. Robinson, McCarthy, and Smyth 2010) to generate the CPM values and then filter. Note that by converting to CPMs we are normalising for the different sequencing depths for each sample

# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)

# As a general rule, a good threshold can be chosen by identifying the CPM that corresponds to a count of 10, which in this case is about 0.5. 
# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

######################## Challenge #########################
# 1) Plot the counts-per-million versus counts for the second sample.
plot(myCPM[,2],countdata[,2],ylim=c(0,50),xlim=c(0,3))

# 2) Add a vertical line at 0.5 and a horizontal line at 10.
plot(myCPM[,2],countdata[,2], ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, h = 10)

# 3) Add the lines again, colouring them blue
plot(myCPM[,2],countdata[,2], ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, h = 10, col = "blue")

# Now that we’ve checked our filtering method we will filter the DGEList object.
y <- y[keep, keep.lib.sizes=FALSE]


################## Quality Control
# Library size and distribution plots
# First, we can check how many reads we have for each sample in the y.
y$samples$lib.size

# We can also plot the library sizes as a barplot to see whether there are any major discrepancies between the samples more easily.
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# we can also adjust the labelling if we want
barplot(y$samples$lib.size/1e06, names=colnames(y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

# Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts
# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Do any samples appear to be different compared to the others?
# ANSwER: Some samples have smaller medians

# Multidimensional scaling plots
plotMDS(y)

# We can colour the samples according to the grouping information. We can also change the labels, or instead of labels we can have points.
# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType)

## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(sampleinfo$Status)

col.status <- c("blue","red","black")[sampleinfo$Status]
col.status

plotMDS(y, col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

## Discussion : Look at the MDS plot coloured by cell type. Is there something strange going on with the samples? Identify the two samples that don’t appear to be in the right place.
# ANSWER: MCL1.LA and MCl1.LB, CL1.DG and CL1.DH

# There is a sample info corrected file in your data directory
# Old sample info
sampleinfo

# I'm going to write over the sampleinfo object with the corrected sample info
sampleinfo <- read.delim("./data/SampleInfo_Corrected.txt", stringsAsFactors = TRUE)
sampleinfo

# We need to correct the info for the groups
group <- factor(paste(sampleinfo$CellType,sampleinfo$Status,sep="."))
y$samples$group <- group

# Redo the MDSplot with corrected information
par(mfrow=c(1,2))
col.cell <- c("purple","orange")[sampleinfo$CellType]
col.status <- c("blue","red","black")[sampleinfo$Status]
plotMDS(y,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")


# Discussion: What is the greatest source of variation in the data (i.e. what does dimension 1 represent)? What is the second greatest source of variation in the data?
# Answer: The greatest source of variation are between groups. Dimenstion 1 represents root-mean-square of the largest 500 log2-fold changes between that pair of samples. The 2nd greatest variation is status

## CHALLENGE
# 1) Redo plots choosing your own color
# 2) Change the plotting character to a symbol instead of the column names HINT: use pch argument. Try pch=16 and see what happens.

par(mfrow=c(1,2))
col.cell <- c("red","lightblue")[sampleinfo$CellType]
col.status <- c("pink","orange","purple")[sampleinfo$Status]
plotMDS(y,col=col.cell, pch=16)
legend("topleft",fill=c("red","lightblue"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(y,col=col.status, pch=16)
legend("topleft",fill=c("pink","orange","purple"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

# Another alternative is to generate an interactive MDS plot using the Glimma package. This allows the user to interactively explore the different dimensions.
labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Status)
glMDSPlot(y, labels=labels, groups=group, folder="mds")


## Hierarchical clustering with heatmaps

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")

# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()


#### Challenge
# 1) Change the colour scheme to “PiYG” and redo the heatmap. Try ?RColorBrewer and see what other colour schemes are available.
# 2) Change the sample names to group using the labCol argument
# 3) Redo the heatmap using the top 500 LEAST variable genes.

# Get the gene names for the top 500 least variable genes
select_var_least <- names(sort(var_genes, decreasing=FALSE))[1:500]
head(select_var_least)
# Subset logcounts matrix
least_variable_lcpm <- logcounts[select_var_least,]
dim(least_variable_lcpm)
## Get some nicer colours
mypalette <- brewer.pal(11,"PiYG")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]
# Plot the heatmap
heatmap.2(least_variable_lcpm, col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row", labCol = group)
###################################################

### Add annotation for both cell type and status using aheatmap function from NMF package
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
aheatmap(highly_variable_lcpm,col=rev(morecols(50)),main="Top 500 most variable genes across samples",annCol=sampleinfo[, 3:4],labCol=group, scale="row")


################## Normalisation for composition bias
# TMM normalization is performed to eliminate composition biases between libraries 
# The calcNormFactors function calculates the normalization factors between libraries. TMM normalisation (and most scaling normalisation methods) scale relative to one sample

# Apply normalisation to DGEList object
y <- calcNormFactors(y)
y$samples

# The last two samples have much smaller normalisation factors, and MCL1.LA and MCL1.LB have the largest. If we plot mean difference plots using the plotMD function for these samples, we should be able to see the composition bias problem. We will use the logcounts, which have been normalised for library size, but not for composition bias.
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")

# The mean-difference plots show average expression (mean: x-axis) against log-fold-changes (difference: y-axis). Because our DGEList object contains the normalisation factors, if we redo these plots using y, we should see the composition bias problem has been solved.
par(mfrow=c(1,2))
plotMD(y,column = 7)
abline(h=0,col="grey")
plotMD(y,column = 11)
abline(h=0,col="grey")

### Challenge
#Plot the biased and unbiased MD plots side by side for the same sample to see the before and after TMM normalisation effect.
par(mfrow=c(2,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(y,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")
plotMD(y,column = 11)
abline(h=0,col="grey")

# We need to save a few data objects to use for Day 2 so we don’t have to rerun everything
save(group,y,logcounts,sampleinfo, seqdata, file="day1objects.Rdata")


