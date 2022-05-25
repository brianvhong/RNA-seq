## Install Packages
if (!requireNamespace("BiocManager"))
        install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))

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
