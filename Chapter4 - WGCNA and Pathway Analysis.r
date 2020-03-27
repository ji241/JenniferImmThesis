##Script to perform the WGCNA to look for odules of co-methylated probes as iPSCs differentiate into neurons
####Date:18/10/2019
##Author: Jennifer Imm

library(WGCNA)
library(cluster)

#load in data
data<-read.csv("/mnt/data1/Jenny/EPIC QC/iPSCsThroughDiff/wo2samples/Dasen normalised data_wo2sample.csv", row.names=1)
dim(data)
#[1] 847103     14


#REMOVE CROSS-HYBRIDISING PROBES
cross<-read.table("/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
crosslist<-cross[,1]
length(crosslist)
#[1]44210  
data.1<-data[ ! rownames(data) %in% crosslist,]
dim(data)  	# 847103     14
dim(data.1) 	# 803777     14


##REMOVE PROBES with SNPs
snpProbes<-read.table("/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
dim(snpProbes)
# 340327     18
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]
dim(snpProbes)
# 10888    18
SNPlist<-snpProbes[,1]
data.2 <-data.1 [! rownames(data.1) %in% SNPlist,]

##REMOVE bad probes
badProbes<-read.csv("/mnt/data1/EPIC_reference/EPICArrayProbesToFilter.csv", stringsAsFactors = FALSE, header = TRUE)
dim(badProbes)
# 977    48
badlist<-badProbes[,1]
data.3 <- data.2[ ! rownames(data.2) %in% badlist,]
dim(data.1) #
#  803777     14
dim(data.2) #
#  794501     14
dim(data.3) #
# 793851     14


#Load in phenotype file
pheno<-read.csv ("/mnt/data1/Jenny/EPIC QC/Phenotypes.csv", stringsAsFactors = FALSE, header = TRUE) ##pheno file

#Remove two failed samples from file
pheno2<- pheno[c(1,4:16),]

#Set column names of data file to be the sample IDs
colnames(data.3)<-pheno2$Sample.ID

#calculate the variance and median of variance
var.jenny<- (apply(data.3, 1,function(x){var(x, na.rm=TRUE)}))#calculate the variance
med <- median(var.jenny)
med
# 0.001137884


##remove non variable probes (probes with variance <=0.001137884 (= 396926 probes)

#create list of non-variable probes
non_var_probes<-which(var.jenny <= med)
length(non_var_probes)
#396926

#Create filtered data file that has non-variable probes removed
JENNY_METHYL_filtered<-data.3[-non_var_probes,]
dim(JENNY_METHYL_filtered)
#[1] 396925     14


#Transpose data
datExpr<-data.frame(t(JENNY_METHYL_filtered))



##Call the network topology analysis function, depending on whether you want a signed or unsigned network this can be changed
# Will help you determine power to use later on
powers = 1:20
options(stringsAsFactors = FALSE);
enableWGCNAThreads(32)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "unsigned", blockSize = 10000)


##plot results
pdf("Unsigned_WGCNA_NetworkTopolgy_wo2samples.pdf",width=10, height=5)
par(mfrow=c(1,2))
## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, Unsigned R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels= powers, cex = 0.6, col = "red")
##this line corresponds to using an R^2 cut off of h
abline(h = 0.9, col = "red")
##Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex=0.6, col = "red")
dev.off()

#Pick soft power based on which value gives mean connectivity closest to 0.9
enableWGCNAThreads(64)
softPower = 20;

#Run the WGCNA
##use this function for large dataset. NOTE maxBlockSize depends on computer memory. Also depending on what your minimum module size required is.
WGCNA = blockwiseModules(datExpr, maxBlockSize = 5000, power = softPower, TOMType = "unsigned", minModuleSize = 500, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs= FALSE, verbose = 3, nThreads = 32)

#save output
save(WGCNA, file = "WGCNA.Rdata")

#Look at how many probes are in each module
table(WGCNA$colors) ## label 0 is reserved for genes outside of all modules
##output 
 0      1      2      3      4      5      6      7      8      9     10
139815  81651  22732  18399  17345  16456  13576  13365  10685   9757   9513
    11     12     13     14     15     16     17     18     19     20     21
  8439   7119   5130   4055   3891   2396   2363   2162   2032   1414   1254
    22     23     24     25
  1129    992    744    511

#Assign the modules a colour
bwModuleLabels<- (WGCNA$colors)
bwModuleColors <- labels2colors((WGCNA$colors)) ## convert labels to colors for plotting
table(bwModuleColors)

#output
bwModuleColors
     black          blue         brown          cyan     darkgreen
        13365         22732         18399          4055          1129
     darkgrey       darkred darkturquoise         green   greenyellow
          744          1254           992         16456          8439
         grey        grey60     lightcyan    lightgreen   lightyellow
       139815          2363          2396          2162          2032
      magenta  midnightblue        orange          pink        purple
         9757          3891           511         10685          9513
          red     royalblue        salmon           tan     turquoise
        13576          1414          5130          7119         81651
       yellow
        17345



#create data frame of data before WGCNA was run and modules
x<-data.frame(bwModuleColors)
y<-data.frame(colnames(datExpr))
probeCol<-cbind(x,y)


#Determine module eigengenes
options(stringsAsFactors = FALSE)
datME = moduleEigengenes(datExpr, bwModuleColors)$eigengenes
x<-signif(cor(datME, use="p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )



nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

## recalculate MEs with color labels
MEs0<-datME
datME = orderMEs(MEs0)

#Create objects containing day of differentiation and timepoint information
timepoint<-pheno2$DayofDifferentiation
day<-pheno2$DayofDifferentiation

#Create matrices for module trait correlations and p-values, to be populated later using loop
moduleTraitCor <- matrix(data = NA, ncol = 1, nrow = ncol(datME))
rownames(moduleTraitCor) <- colnames(datME)
colnames(moduleTraitCor)<-"Days"

moduleTraitPvalue <- matrix(data = NA, ncol = 1, nrow = ncol(datME))
rownames(moduleTraitPvalue) <- colnames(datME)
colnames(moduleTraitPvalue)<-"Days"


#load and intstall packages
library(nlme)
install.packages("QuantPsyc")
library("QuantPsyc")

dat<-data.frame(datME, timepoint, day)

#use loop to work out module trait correlations and p-values
for(i in 1:26){
module<-as.numeric(datME[,i])
	try(res1<-lm(module ~ timepoint, dat=dat), silent = TRUE)
	if(class(res1) != "try-error")
		moduleTraitCor[i,1]<- as.numeric(lm.beta(res1)[1])
		moduleTraitPvalue[i,1]<-coef(summary(res1))["timepoint",4]
	}

#create plot showing module eigengene for each module over time
pdf("boxplot.modules_wo2samples.pdf", width=20, height=10)
par(mfrow = c(2,5))
	for(i in 1:26){
	boxplot(datME[,i]~ timepoint, main=as.character(colnames(datME[i])))
	}
dev.off()	


#create heatmap of module trait correlations and p-values for each module
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf(file = "Module-trait-relationships_wo2samples.pdf", height=10)
par(mar = c(8, 10, 3, 3))
## display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()




#calculate the module membership for blue module 

MEblue<-probeCol[which(probeCol$bwModuleColors=="blue"),2]
length(MEblue)
#22732
allcpgs<-colnames(datExpr)

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

#calculate module membership
MM.blue<-data.frame(datKME$MM.blue, row.names=row.names(datKME))
MM.blue.filtered<-abs(data.frame(MM.blue[MEblue,], row.names=row.names(datKME[MEblue,])))


#Calculate the Gene Significance levels 
GS1=as.numeric(cor(timepoint,datExpr, use="p"))
length(GS1)
#396925
GeneSignificance= abs(GS1)
GS<-data.frame(rownames=colnames(datExpr),GeneSignificance, row.names=1)

summary(GS)
 summary(GS)
 GeneSignificance
 Min.   :0.0000029
 1st Qu.:0.2083207
 Median :0.4013196
 Mean   :0.4078796
 3rd Qu.:0.5917121
 Max.   :0.9779917


#calculate correlation between module membership and gene significance
GS.blue.filtered<-data.frame(GS[MEblue,])
cor.test(MM.blue.filtered[,1],GS.blue.filtered[,1])

###output
#Pearson's product-moment correlation

        Pearson's product-moment correlation

data:  MM.blue.filtered[, 1] and GS.blue.filtered[, 1]
t = 224.16, df = 22730, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.8256900 0.8337887
sample estimates:
     cor
0.829783


#Repeat for black module

#calculate the module membership
MEblack<-probeCol[which(probeCol$bwModuleColors=="black"),2]
length(MEblack)
#13365
allcpgs<-colnames(datExpr)

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

MM.black<-data.frame(datKME$MM.black, row.names=row.names(datKME))
MM.black.filtered<-abs(data.frame(MM.black[MEblack,], row.names=row.names(datKME[MEblack,])))


#Calculate the Gene Significance levels 
GS1=as.numeric(cor(timepoint,datExpr, use="p"))
length(GS1)
#396925
GeneSignificance= abs(GS1)
GS<-data.frame(rownames=colnames(datExpr),GeneSignificance, row.names=1)

summary(GS)
GeneSignificance
 Min.   :0.0000029
 1st Qu.:0.2083207
 Median :0.4013196
 Mean   :0.4078796
 3rd Qu.:0.5917121
 Max.   :0.9779917

#calculate correlation
GS.black.filtered<-data.frame(GS[MEblack,])
cor.test(MM.black.filtered[,1],GS.black.filtered[,1])
#output Pearson's product-moment correlation

       
        Pearson's product-moment correlation

data:  MM.black.filtered[, 1] and GS.black.filtered[, 1]
t = 117.51, df = 13363, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.7044360 0.7211146
sample estimates:
      cor
0.7128761


#Repeat for greenyellow module

#calculate the module membership
MEgreenyellow<-probeCol[which(probeCol$bwModuleColors=="greenyellow"),2]
length(MEgreenyellow)
#8439
allcpgs<-colnames(datExpr)

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

MM.greenyellow<-data.frame(datKME$MM.greenyellow, row.names=row.names(datKME))
MM.greenyellow.filtered<-abs(data.frame(MM.greenyellow[MEgreenyellow,], row.names=row.names(datKME[MEgreenyellow,])))


#Calculate the Gene Significance levels 
GS1=as.numeric(cor(timepoint,datExpr, use="p"))
length(GS1)
#396925
GeneSignificance= abs(GS1)
GS<-data.frame(rownames=colnames(datExpr),GeneSignificance, row.names=1)


summary(GS)
GeneSignificance
 Min.   :0.0000029
 1st Qu.:0.2083207
 Median :0.4013196
 Mean   :0.4078796
 3rd Qu.:0.5917121
 Max.   :0.9779917

#calculate correlation
GS.greenyellow.filtered<-data.frame(GS[MEgreenyellow,])
cor.test(MM.greenyellow.filtered[,1],GS.greenyellow.filtered[,1])
#output Pearson's product-moment correlation

        Pearson's product-moment correlation

data:  MM.greenyellow.filtered[, 1] and GS.greenyellow.filtered[, 1]
t = 149.89, df = 8437, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.8467021 0.8583560
sample estimates:
      cor
0.8526351


#Repeat for red module

#calculate the module membership
MEred<-probeCol[which(probeCol$bwModuleColors=="red"),2]
length(MEred)
#13576
allcpgs<-colnames(datExpr)

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

MM.red<-data.frame(datKME$MM.red, row.names=row.names(datKME))
MM.red.filtered<-abs(data.frame(MM.red[MEred,], row.names=row.names(datKME[MEred,])))


#Calculate the Gene Significance levels 
GS1=as.numeric(cor(timepoint,datExpr, use="p"))
length(GS1)
#396925
GeneSignificance= abs(GS1)
GS<-data.frame(rownames=colnames(datExpr),GeneSignificance, row.names=1)

summary(GS)
GeneSignificance
 Min.   :0.0000029
 1st Qu.:0.2083207
 Median :0.4013196
 Mean   :0.4078796
 3rd Qu.:0.5917121
 Max.   :0.9779917

#Calculate correlation
GS.red.filtered<-data.frame(GS[MEred,])
cor.test(MM.red.filtered[,1],GS.red.filtered[,1])
#output Pearson's product-moment correlation

        Pearson's product-moment correlation

data:  MM.red.filtered[, 1] and GS.red.filtered[, 1]
t = 186.29, df = 13574, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.8430404 0.8525018
sample estimates:
      cor
0.8478385



#Save data
save(MEblack,MEblue,MEgreenyellow, MEred, allcpgs, file="GOWGCNAMODULES_wo2samples.Rdata")



#############################################
PATHWAY ANALYSIS ON WGCNA MODULES 
############################################

#Only run pathway analysis on top 15% probes from each module based on module memnbership

#filter based on module membership
bluefiltered <- MM.blue.filtered[abs(MM.blue.filtered[,1])>quantile(abs(MM.blue.filtered[,1]),0.85),,drop=F]
dim(bluefiltered)
#[1] 3410    1

#####Pathway Analysis for significatnt modules identified in WGCNA

options(stringsAsFactors = FALSE) 

source("/mnt/data1/Jenny/WGCNA/crystalmeth.R")

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")


###Blue Module###

#set sig cpgs to be probe IDs of blue module
sigcpg<-rownames(bluefiltered)
#set all cpgs to be probes IDs from data before WGCNA was run ie normalised anf filtered betas
allcpg <- allcpgs

sigcpg <- as.character(sigcpg)

projectName <- "iPSC through Diff WGCNA Blue filtered"

scriptVersion <- "-v1"

##Run Go pathway Analysis
go <- crystalmeth(sig.cpg = sigcpg, 
                  all.cpg = allcpg,
                  collection = "GO",array.type = "EPIC" )

go<-go[order(go$FDR),]


#save output
write.csv(go, "GO Gene Analysis Output WGCNA Blue Filtered.csv")
save(go, file = "GO Gene Analysis Output WGCNA Blue Filtered.Rdata")


##################################################

#filter and run pathway analysis on black module 

blackfiltered <- MM.black.filtered[abs(MM.black.filtered[,1])>quantile(abs(MM.black.filtered[,1]),0.85),,drop=F]
dim(blackfiltered)
#[1] 2005    1


###Black Module
sigcpg<-rownames(blackfiltered)
allcpg <- allcpgs

sigcpg <- as.character(sigcpg)

projectName <- "iPSC through Diff WGCNA Black filtered"

scriptVersion <- "-v1"

##Go pathway Analysis
go <- crystalmeth(sig.cpg = sigcpg, 
                  all.cpg = allcpg,
                  collection = "GO",array.type = "EPIC" )

go<-go[order(go$FDR),]


#save output
write.csv(go, "GO Gene Analysis Output WGCNA black Filtered.csv")
save(go, file = "GO Gene Analysis Output WGCNA black Filtered.Rdata")

##############################################

#filter and run pathway analysis on greenyellow module 

greenyellowfiltered <- MM.greenyellow.filtered[abs(MM.greenyellow.filtered[,1])>quantile(abs(MM.greenyellow.filtered[,1]),0.85),,drop=F]
dim(greenyellowfiltered)
#[1] 1266    1

##
###Greenyellow Module
sigcpg<-rownames(greenyellowfiltered)
allcpg <- allcpgs

sigcpg <- as.character(sigcpg)

projectName <- "iPSC through Diff WGCNA Greenyellow filtered"

scriptVersion <- "-v1"

##Go pathway analysis Analysis
go <- crystalmeth(sig.cpg = sigcpg, 
                  all.cpg = allcpg,
                  collection = "GO",array.type = "EPIC" )

go<-go[order(go$FDR),]


#save output
write.csv(go, "GO Gene Analysis Output WGCNA Greenyellow Filtered.csv")
save(go, file = "GO Gene Analysis Output WGCNA Greenyellow Filtered.Rdata")

#############################################################

#filter and run pathway analysis on red module 

redfiltered <- MM.red.filtered[abs(MM.red.filtered[,1])>quantile(abs(MM.red.filtered[,1]),0.85),,drop=F]
dim(redfiltered)
#[1] 2037    1

##
###red Module
sigcpg<-rownames(redfiltered)
allcpg <- allcpgs

sigcpg <- as.character(sigcpg)

projectName <- "iPSC through Diff WGCNA red filtered"

scriptVersion <- "-v1"

##Go Analysis
go <- crystalmeth(sig.cpg = sigcpg, 
                  all.cpg = allcpg,
                  collection = "GO",array.type = "EPIC" )

go<-go[order(go$FDR),]


#save output
write.csv(go, "GO Gene Analysis Output WGCNA Red Filtered.csv")
save(go, file = "GO Gene Analysis Output WGCNA Red Filtered.Rdata")



