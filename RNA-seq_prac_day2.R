################## Differential expression with limma-voom (Day 2)
## Library
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
# Load the objects into the workspace that we created yesterday
load("day1objects.Rdata")
objects()


## Create the design matrix
# Look at group variable again
group

# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design

## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design

## Voom transform the data: Voom will automatically adjust the library sizes using the norm.factors already calculate
#Check mean-variance plot
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
v

# What is contained in this object?
names(v)

## Challenge
# What is in the targets slot of v and what does it correspond to in y?
# What are the dimensions of the weights slot in v?
# Answer: Target in V indicate the group, library size and norm factor. It correspond to sample in y. 15804x12 for dimensions

# We can repeat the box plots for the normalised data to compare to before normalisation. The expression values in v$E are already log2 values so we don’t need to log-transform.
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

## Testing for differential expression
#  First we fit a linear model for each gene using the lmFit function in limma. lmFit needs the voom object and the design matrix that we have already specified, which is stored within the voom object.
# Fit the linear model
fit <- lmFit(v)
names(fit)

# we are interested in knowing which genes are differentially expressed between the pregnant and lactating group in the basal cells. This is done by defining the null hypothesis as basal.pregnant - basal.lactate = 0 for each gene. Note that the group names must exactly match the column names of the design matrix.
cont.matrix <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate, L.PregVsLac = luminal.pregnant - luminal.lactate, levels=design)
cont.matrix
# Now we can apply the contrasts matrix to the fit object to get the statistics and estimated parameters of our comparison that we are interested in. Here we call the contrasts.fit function in limma.
fit.cont <- contrasts.fit(fit, cont.matrix)
# The final step is to call the eBayes function, which performs empirical Bayes shrinkage on the variances, and estimates moderated t-statistics and the associated p-values.
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

# We can use the limma decideTests function to generate a quick summary of DE genes for the contrasts.
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

vennDiagram(summa.fit, include=c("up", "down"),
            counts.col=c("purple", "black"))

# Plots after testing for DE
# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"B.PregVsLac"], values = c(-1, 1), hl.col=c("blue","red"))

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL, main="B.PregVsLac")




### Do it for L.PregvsLac
par(mfrow=c(1,2))

plotMD(fit.cont, coef=1,status=summa.fit[,"L.PregVsLac"], values = c(-1, 1), hl.col=c("blue","red"))

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=200,names=fit.cont$genes$SYMBOL, main="L.PregVsLac")

########################

# Before following up on the DE genes with further lab work, it is recommended to have a look at the expression levels of the individual samples for the genes of interest. 
# We can quickly look at grouped expression using stripchart. We can use the normalised log expression values in the voom object (v$E).
par(mfrow=c(1,3))
# Let's look at the first gene in the topTable, Wif1, which has a rowname 24117
stripchart(v$E["24117",]~group)
# This plot is ugly, let's make it better
stripchart(v$E["24117",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,col=1:6,method="jitter")
# Let's use nicer colours
nice.col <- brewer.pal(6,name="Dark2")
stripchart(v$E["24117",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="Wif1")


## CHALLENGE
# Take the top gene from the L.PregVsLactate comparison and make a stripchart of grouped expression as above. (Don’t forget to change the title of the plot.)
par(mfrow=c(1,3))

stripchart(v$E["12992",]~group)
# This plot is ugly, let's make it better
stripchart(v$E["12992",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,col=1:6,method="jitter")
# Let's use nicer colours
nice.col <- brewer.pal(6,name="Dark2")
stripchart(v$E["12992",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="Csn1s2b")




# An interactive version of the volcano plot above that includes the raw per sample values in a separate panel is possible via the glXYPlot function in the Glimma package.
group2 <- group
levels(group2) <- c("basal.lactate","basal.preg","basal.virgin","lum.lactate", "lum.preg", "lum.virgin")
glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
         xlab="logFC", ylab="B", main="B.PregVsLac",
         counts=v$E, groups=group2, status=summa.fit[,1],
         anno=fit.cont$genes, side.main="ENTREZID", folder="volcano")