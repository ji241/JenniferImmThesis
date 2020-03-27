##Script to assess quality of Illumina EPIC data of iPSCs as they differentiate into neurons
####Date:10.12.18
##Author: Jennifer Imm


#set working directory
setwd("")

#load required packages
library(methylumi)
library(wateRmelon)
require(gdata)

### load phenotype data, DNA methylation data and probe annotation information
pheno<-read.csv ("Phenotypes.csv", stringsAsFactors = FALSE, header = TRUE) ##pheno file
pheno$Basename<-as.character(pheno$Basename) ##changing the Basenames to character strings
epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7) ##probe annotation file


#Set the idat file path so the HPC knows where to find them when they are called upon
idatPath<-c("idats/")

#Create an mset file which matches each idat file to the relevant sample, this is done by BeadChip barcode.
#This also converts the idats to methylumi objects so the fluorescence intensities can be interpreted in the later steps
msetEPIC <- readEPIC(idatPath=idatPath, barcodes = pheno$Basename, parallel = FALSE)

#save file
save(msetEPIC, file = "/mnt/data1/Jenny/EPIC QC/iPSCsThroughDiff/wo2samples/MsetEPIC_wo2samples.Rdata")

#To reload next time
#load("/mnt/data1/Jenny/EPIC QC/iPSCsThroughDiff/wo2samples/MsetEPIC_wo2samples.Rdata", verbose = FALSE)


##matches betas msetEPIC file to epic manifest column Name
epicManifest<-epicManifest[match(rownames(betas(msetEPIC)), epicManifest$Name),]

######################################
####QC Pipeline#######################
######################################

#Step 1. Assessing the fluorescence intensities

#1a. Checking median intensities

### extract sample intensities and plot
m_intensities<-methylated(msetEPIC)
u_intensities<-unmethylated(msetEPIC)

M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

#M.median output
200821610056_R01C01 200821610056_R02C01 200821610057_R01C01 200821610057_R02C01
               2320                1136                 992                1781
200821610056_R03C01 200821610056_R04C01 200821610057_R03C01 200821610057_R04C01
               4369                3769                3645                3105
200821610056_R05C01 200821610056_R06C01 200821610057_R05C01 200821610057_R06C01
               3394                3298                3637                3633
200821610056_R07C01 200821610056_R08C01 200821610057_R07C01 200821610057_R08C01
               3315                3476                2748                3407

#U.median output 
200821610056_R01C01 200821610056_R02C01 200821610057_R01C01 200821610057_R02C01
               1429                1032                 811                1146
200821610056_R03C01 200821610056_R04C01 200821610057_R03C01 200821610057_R04C01
               2199                2001                1667                1632
200821610056_R05C01 200821610056_R06C01 200821610057_R05C01 200821610057_R06C01
               1767                1572                1496                1403
200821610056_R07C01 200821610056_R08C01 200821610057_R07C01 200821610057_R08C01
               1819                1789                1232                1480

#save as tables
write.table(M.median,file="M.median.txt",quote=F,col.names=F)
write.table(U.median,file="U.median.txt",quote=F,col.names=F)

#create output file that will have information on each QC metric
QCmetrics<-cbind(pheno, M.median, U.median) ## create a table to store output of QC pipeline

create Chip ID column in pheno file
ChipID <- c(200821610056, 200821610056, 200821610057, 200821610057, 200821610056, 200821610056, 200821610057, 200821610057, 200821610056, 200821610056, 200821610057, 200821610057, 200821610056, 200821610056, 200821610057, 200821610057)
pheno <- cbind(pheno, ChipID)


intens.Thres<-1000 ## change this to adjust the threshold at which you filter 

#create scatterplot and histogram of median intensities without cut offs
pdf("Scatterplot_SampleIntensity.pdf")
par(mfrow = c(1,2))
hist(M.median, xlab = "Median Methylated Intensity")
hist(U.median, xlab = "Median Unmethylated Intensity")
par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median Methylated Intensity", ylab = "Median Unmethylated Intensity", main="Median Sample Intensities", cex=2)
text(M.median, U.median, labels=pheno$Sample.ID, cex=0.9, pos=1)
dev.off()

#create same plots with cut offs
pdf("Scatterplot_SampleIntensity_Abline.pdf")
par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity")
hist(U.median, xlab = "Median U intensity")
par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity", xlim = c(500,3000))
abline(v = intens.Thres, col = "red")
abline(h = intens.Thres, col = "red")
dev.off()


#1b. Assessing if BeadChip or position on BeadChip affect median intensities.

## calculate a summary statistic for each chip
chip.M.median<-aggregate(M.median, by = list(unlist(strsplit(colnames(m_intensities), "_"))[seq(from = 1, to = 2*ncol(m_intensities), by = 2)]), FUN = median)
chip.U.median<-aggregate(U.median, by = list(unlist(strsplit(colnames(m_intensities), "_"))[seq(from = 1, to = 2*ncol(u_intensities), by = 2)]), FUN = median)


##boxplot of median intensities coloured by BeadChip ID
pdf("iPSCsThroughDiff/wo2samples/Boxplot_SampleIntensity_ByPlate.pdf")
par(mfrow = c(1,2))
par(mar = c(8, 4, 1, 1))
nCol<-length(unique(pheno$ChipID))
boxplot(M.median ~ pheno$ChipID, ylab = "Median Methylated Intensity", xlab = "", las = 2, col = rainbow(nCol)) 
mtext(text="Chip ID", side=1, line=6.8) 
boxplot(U.median ~ pheno$ChipID, ylab = "Median Unmethylated Intensity", xlab = "", las = 2, col = rainbow(nCol))
mtext(text="Chip ID", side=1, line=6.8) 
mtext("Median Sample Intensities by Chip ID", outer=TRUE, line=-1)
dev.off()

##scatterplot of median intensities against one another coloured by BeadChip ID
pdf("iPSCsThroughDiff/wo2samples/Scatterplot_SampleIntensity_ColByPlate2.pdf")
nCol<-length(unique(pheno$ChipID))## assumes there is a column called Plate in your phenotype file
plot(M.median, U.median, pch = 16, cex=2, xlab = "Median Methylated Intensity", ylab = "Median Unmethylated Intensity", ylim= c(0,2700), col = rainbow(nCol)[factor(pheno$ChipID)])
legend("topleft", levels(factor(pheno$ChipID)), col = rainbow(nCol), pch = 16, title="Chip ID")
mtext("Median Sample Intensities by Chip ID", outer=TRUE, line=-4)
dev.off()

## boxplot of median intensities coloured by their position on the BeadChip
pdf("iPSCsThroughDiff/wo2samples/Boxplot_SampleIntensity_ByPos.pdf")
par(mfrow = c(1,2))
par(mar = c(8, 4, 1, 1))
nCol<-length(unique(pheno$Meth.Chip.Position))## assumes there is a column called Plate in your phenotype file
boxplot(M.median ~ pheno$Meth.Chip.Position, ylab = "Median Methylated Intensity", las = 2, col = rainbow(nCol), xlab="")
mtext(text="Position on Chip", side=1, line=5) 
boxplot(U.median ~ pheno$Meth.Chip.Position, ylab = "Median Unmethylated Intensity", las = 2, col = rainbow(nCol), xlab="")
mtext(text="Position on Chip", side=1, line=5) 
mtext("Median Sample Intensities by Position", outer=TRUE, line=-1)
dev.off()

#Scatterplot of median intensities colouring by position on BeadChip
pdf("iPSCsThroughDiff/wo2samples/Scatterplot_SampleIntensity_ColByPos2.pdf")
nCol<-length(unique(pheno$Meth.Chip.Position))
plot(M.median, U.median, pch = 16, cex=2, xlab = "Median Methylated Intensity", ylab = "Median Unmethylated Intensity", ylim= c(0,2500), col = rainbow(nCol)[factor(pheno$Meth.Chip.Position)])
legend("topleft", levels(factor(pheno$Meth.Chip.Position)), col = rainbow(nCol), pch = 16, title="Position on Chip")
mtext("Median Sample Intensities by Position", outer=TRUE, line=-4)
dev.off()

##2. Calculate bisulphite conversion statistic

bs<-bscon(msetEPIC)
#bs output
200821610056_R01C01 200821610056_R02C01 200821610057_R01C01 200821610057_R02C01
           93.50013            89.76363            90.78108            92.42073
200821610056_R03C01 200821610056_R04C01 200821610057_R03C01 200821610057_R04C01
           94.24837            92.37949            93.85473            91.10159
200821610056_R05C01 200821610056_R06C01 200821610057_R05C01 200821610057_R06C01
           90.80429            92.41864            93.11873            90.62292
200821610056_R07C01 200821610056_R08C01 200821610057_R07C01 200821610057_R08C01
           91.40210            92.60794            92.95649            93.60908
#Save as table
write.table(bs,file="/mnt/data1/Jenny/EPIC QC/iPSCsThroughDiff/wo2samples/bscon.txt",quote=F,col.names=F)
		   
##Create histogram showing the BScon statistic for each sample
pdf("iPSCsThroughDiff/wo2samples/HistogramBisulphiteConversionStatistics.pdf")
hist(bs, xlab = "Median % Bisulfite Conversion", main = "Bisulfite Conversion Efficiency", col = c("#87ceeb", "#eba487"))
dev.off()


##3. Beta Density plot

pheno<-pheno[match(colnames(betas), pheno$Basename),]	   
		   
#transpose beta values
t(betas)->betast
rownames(pheno)<-pheno$Basename
#merge pheno file with transposed betas
pheno_betas <- merge(betast, pheno, by = "row.names")

##create density plot of beta values from rows 2 to 866905
##create density line for each sample individually using denisty() and then plot them all on one plot using lines function
plot(density(pheno_betas[,2:866905]), main = "Density plot of beta values", col = "blue")

density(betas[,1], na.rm=T)->TEST1
density(betas[,4], na.rm=T)->TEST4
density(betas[,5], na.rm=T)->TEST5
density(betas[,6], na.rm=T)->TEST6
density(betas[,7], na.rm=T)->TEST7
density(betas[,8], na.rm=T)->TEST8
density(betas[,9], na.rm=T)->TEST9
density(betas[,10], na.rm=T)->TEST10
density(betas[,11], na.rm=T)->TEST11
density(betas[,12], na.rm=T)->TEST12
density(betas[,13], na.rm=T)->TEST13
density(betas[,14], na.rm=T)->TEST14
density(betas[,15], na.rm=T)->TEST15
density(betas[,16], na.rm=T)->TEST16

#create plot
pdf("iPSCsThroughDiff/wo2samples/betas_density_plot.pdf")
plot(TEST1, col="red", main = "Density Plot of Beta Values", ylim = c(0, 4.5))
lines(TEST4, col="red")
lines(TEST5, col="purple")
lines(TEST6, col="purple")
lines(TEST7, col="purple")
lines(TEST8, col="purple")
lines(TEST9, col="blue")
lines(TEST10, col="blue")
lines(TEST11, col="blue")
lines(TEST12, col="blue")
lines(TEST13, col="green")
lines(TEST14, col="green")
lines(TEST15, col="green")
lines(TEST16, col="green")
legend("topleft", c("iPSC", "NPC", "Neuron 37", "Neuron 58"), col = c("red", "purple", "blue", "green"), pch = 18, title="Cell Stage")
dev.off()



### 4. P filter
msetEPIC.pf <- pfilter(msetEPIC)
pFilterPass<-colnames(betas(msetEPIC)) %in% colnames(betas(msetEPIC.pf))
#output
2 samples having 1 % of sites with a detection p-value greater than 0.05 were removed
Samples removed: 200821610056_R02C01 200821610057_R01C01
12962 sites were removed as beadcount <3 in 5 % of samples
7170 sites having 1 % of samples with a detection p-value greater than 0.05 were removed


###5. Mitotic age

#Remove two failed samples from p-filter stage before running mito age calculator
#Subset pheno file to remove samples
pheno2<- pheno[c(1,4:16),]

#Load in function needed to run mito age calculator
source("/mnt/data1/reference_files/MiAge/function_library.r") ### library of all functions used in the calculation of mitotic ages                                                                                 


clocksitesID=as.matrix(read.csv("/mnt/data1/reference_files/MiAge/Additional_File1.csv",header=T))[,1]   ###clock CpG sites

#ensure betas are in correct format
betas2 <- betas(msetEPIC)


dim(betas2)
#866895     16

#subset betas to remove two failed samples
betas2 <-  betas2[,c(1,4:16)]

#select clock specific cpg sites from betas
beta=  betas2[ match(clocksitesID,rownames(betas2)),]

#check correct dimensions
dim(beta)
     268       14
	 
#load site specific parameters
load("/mnt/data1/reference_files/MiAge/site_specific_parameters.Rdata") #### site-specific parameters

#Run the mitotic age clock
b=methyl.age[[1]];c=methyl.age[[2]];d=methyl.age[[3]]

n=mitotic.age(beta,b,c,d) ### estimated mitotic age
names(n)=colnames(beta)
write.table(n,file="/mnt/data1/Jenny/EPIC QC/iPSCsThroughDiff/wo2samples/mitotic_age.txt",quote=F,col.names=F)

#bind output to QCmetrics object
QCmetrics=cbind(QCmetrics, n)
 
#Create a plot showing the mitotic age of each time point
QCmetrics$Day.of.Differentiation=as.factor(QCmetrics$Day.of.Differentiation)
pdf("/mnt/data1/Jenny/EPIC QC/iPSCsThroughDiff/wo2samples/mitotic_age_prediction.pdf")
plot(QCmetrics$CellStage,QCmetrics$n,ylab="Estimated no. of Cell Divisions",xlab="Cell Stage",main="Mitotic Age Prediction", col = c("#ff00ff", 
"#ffcc00", "#9650a1", "#45b449"))
legend("topleft", c("iPSC", "NPC", "Neuron-D37", "Neuron-D58"), col = c("#ff00ff", "#ffcc00", "#9650a1", "#45b449"), pch = 18, title="Cell Stage")
dev.off()


##Save QC file
write.csv(QCmetrics, file="/mnt/data1/Jenny/EPIC QC/iPSCsThroughDiff/wo2samples/QC_without_2_samples.csv")

##save pfiltered betas
save(msetEPIC.pf, file = "iPSCsThroughDiff/without 2 sample/QC'd betas_wo2_Samples")



#################################################################################
#NOMRALISATION USING DASEN
#################################################################################

#Pull out type I and II probes
EPIC.data <- msetEPIC.pf@featureData
onetwod <- EPIC.data$DESIGN
proId <- which(onetwod == "I")
proIId <- which(onetwod == "II")
EPIC.data_I  <- msetEPIC.pf[proId,]
EPIC.data_II <- msetEPIC.pf[proIId,]

#calculate betas of pfiltered data
betas.pf   <- betas(msetEPIC.pf)
betas.pfI  <- betas(EPIC.data_I)
betas.pfII <- betas(EPIC.data_II)

#Run dasen normalisation on p-filtered data
EPIC.dasen<-dasen(msetEPIC.pf)


#pull out the methylated and unmethylated intensities of normalised data
m_intensities_dasen<-methylated(EPIC.dasen)
u_intensities_dasen<-unmethylated(EPIC.dasen)
betas_d<-betas(EPIC.dasen)

#calculate median methylated and unmethylated intensities after normalisation
M.median_dasen<-apply(m_intensities_dasen, 2, median)
U.median_dasen<-apply(u_intensities_dasen, 2, median)

##M.median_dasen output
200821610056_R01C01 200821610057_R02C01 200821610056_R03C01 200821610056_R04C01
           3340.000            3340.000            3340.000            3340.214
200821610057_R03C01 200821610057_R04C01 200821610056_R05C01 200821610056_R06C01
           3339.943            3340.214            3340.357            3340.357
200821610057_R05C01 200821610057_R06C01 200821610056_R07C01 200821610056_R08C01
           3340.429            3339.786            3340.429            3340.357
200821610057_R07C01 200821610057_R08C01
           3340.214            3340.429

		   
##U.median_dasen output
200821610056_R01C01 200821610057_R02C01 200821610056_R03C01 200821610056_R04C01
           1636.643            1636.643            1636.643            1636.643
200821610057_R03C01 200821610057_R04C01 200821610056_R05C01 200821610056_R06C01
           1636.643            1637.286            1636.500            1637.214
200821610057_R05C01 200821610057_R06C01 200821610056_R07C01 200821610056_R08C01
           1637.214            1637.500            1637.168            1636.500
200821610057_R07C01 200821610057_R08C01
           1636.500            1636.643

#subset pheno file so it doesnt contain two failed samples
pheno_3 <- pheno[-c(2,3),]


#Create a scatterplot and histogram of normalised median intensities, label each point with the sample ID
pdf("iPSCsThroughDiff/wo2samples/Scatterplot_MEDIANIntensities_EPIC.Dasen.pdf")		   
par(mfrow = c(1,2))
hist(M.median_dasen, xlab = "Median Methylated intensity")
hist(U.median_dasen, xlab = "Median Unmethyalted intensity")
par(mfrow = c(1,1))
plot(M.median_dasen, U.median_dasen, pch = 16, cex=2, xlab = "Median Methylated Intensity", ylab = "Median Unmethylated Intensity", main="Median Sample Intensities after Normalisation")
text(M.median_dasen, U.median_dasen, labels=pheno_3$Sample.ID, cex=0.9, pos=1)
dev.off()


##Create a beta density plot of normalised data
#transpose normalised betas
t(betas_d)->betast
rownames(pheno)<-pheno$Basename
#merge normalised betas with pheno file
pheno_betas_dasen <- merge(betast, pheno, by = "row.names")

dim(pheno_betas_dasen)
## 14 847113

###density plot of Dasen normalised data, coloured by cell stage
plot(density(pheno_betas_dasen[,2:847113]), main = "Density plot of beta values after Dasen", col = "blue")

density(betas_d[,1], na.rm=T)->TEST1
density(betas_d[,2], na.rm=T)->TEST2
density(betas_d[,3], na.rm=T)->TEST3
density(betas_d[,4], na.rm=T)->TEST4
density(betas_d[,5], na.rm=T)->TEST5
density(betas_d[,6], na.rm=T)->TEST6
density(betas_d[,7], na.rm=T)->TEST7
density(betas_d[,8], na.rm=T)->TEST8
density(betas_d[,9], na.rm=T)->TEST9
density(betas_d[,10], na.rm=T)->TEST10
density(betas_d[,11], na.rm=T)->TEST11
density(betas_d[,12], na.rm=T)->TEST12
density(betas_d[,13], na.rm=T)->TEST13
density(betas_d[,14], na.rm=T)->TEST14

#create the plot
pdf("iPSCsThroughDiff/wo2samples/betas_density_plot_dasen.pdf")
plot(TEST1, col="red", main = "Density Plot of Beta Values after Normalisation", ylim = c(0, 4.5))
lines(TEST2, col="red")
lines(TEST3, col="purple")
lines(TEST4, col="purple")
lines(TEST5, col="purple")
lines(TEST6, col="purple")
lines(TEST7, col="blue")
lines(TEST8, col="blue")
lines(TEST9, col="blue")
lines(TEST10, col="blue")
lines(TEST11, col="green")
lines(TEST12, col="green")
lines(TEST13, col="green")
lines(TEST14, col="green")
legend("topleft", c("iPSC", "NPC", "Neuron 37", "Neuron 58"), col = c("red", "purple", "blue", "green"), pch = 18, title="Cell Stage")
dev.off()

#Pull out type I and II probes
data.dasen <- EPIC.dasen@featureData
onetwod <- data.dasen$DESIGN
proId <- which(onetwod == "I")
proIId <- which(onetwod == "II")
EPIC.dasen_I  <- EPIC.dasen[proId,]
EPIC.dasen_II <- EPIC.dasen[proIId,]

#calculate betas of dasen data of each probe type
betas.dasen   <- betas(EPIC.dasen)
betas.dasenI  <- betas(EPIC.dasen_I)
betas.dasenII <- betas(EPIC.dasen_II)  

##create CSV files for normalised data
write.csv(betas.dasen, "iPSCsThroughDiff/without 2 sample/Dasen normalised data_wo2sample.csv")
write.csv(betas.dasenI, "iPSCsThroughDiff/without 2 sample/Dasen normalised data type I probes_wo2sample.csv")
write.csv(betas.dasenII, "iPSCsThroughDiff/without 2 sample/Dasen normalised data type II probes_wo2sample.csv")



###Create beta density plot for only the type I probes 
#transpose
t(betas.dasenI)->betast
rownames(pheno)<-pheno$Basename
#merge
pheno_betas_dasen <- merge(betast, pheno, by = "row.names")

dim(pheno_betas_dasen)
## 14 138055

###density plot of Dasen normalised data, coloured by cell stage
plot(density(betas.dasenI[,2:138055]), main = "Type I Probes after Normalisation", col = "blue")

density(betas.dasenI[,1], na.rm=T)->TEST1
density(betas.dasenI[,2], na.rm=T)->TEST2
density(betas.dasenI[,3], na.rm=T)->TEST3
density(betas.dasenI[,4], na.rm=T)->TEST4
density(betas.dasenI[,5], na.rm=T)->TEST5
density(betas.dasenI[,6], na.rm=T)->TEST6
density(betas.dasenI[,7], na.rm=T)->TEST7
density(betas.dasenI[,8], na.rm=T)->TEST8
density(betas.dasenI[,9], na.rm=T)->TEST9
density(betas.dasenI[,10], na.rm=T)->TEST10
density(betas.dasenI[,11], na.rm=T)->TEST11
density(betas.dasenI[,12], na.rm=T)->TEST12
density(betas.dasenI[,13], na.rm=T)->TEST13
density(betas.dasenI[,14], na.rm=T)->TEST14

#create plot
pdf("iPSCsThroughDiff/wo2samples/betas_density_plot_dasenI.pdf")
plot(TEST1, col="red", main = "Type I Probes after Normalisation", ylim = c(0, 6))
lines(TEST2, col="red")
lines(TEST3, col="purple")
lines(TEST4, col="purple")
lines(TEST5, col="purple")
lines(TEST6, col="purple")
lines(TEST7, col="blue")
lines(TEST8, col="blue")
lines(TEST9, col="blue")
lines(TEST10, col="blue")
lines(TEST11, col="green")
lines(TEST12, col="green")
lines(TEST13, col="green")
lines(TEST14, col="green")
legend("topleft", c("iPSC", "NPC", "Neuron 37", "Neuron 58"), col = c("red", "purple", "blue", "green"), pch = 18, title="Cell Stage")
dev.off()


##beta density plot showing dasen normalised data but only type II probes
#transpose
t(betas.dasenII)->betast
rownames(pheno)<-pheno$Basename
#merge
pheno_betas_dasen <- merge(betast, pheno, by = "row.names")

dim(pheno_betas_dasen)
## 14 709068

###density plot of Dasen normalised data, coloured by cell stage
plot(density(betas.dasenI[,2:709068]), main = "Type II Probes after Normalisation", col = "blue")

density(betas.dasenII[,1], na.rm=T)->TEST1
density(betas.dasenII[,2], na.rm=T)->TEST2
density(betas.dasenII[,3], na.rm=T)->TEST3
density(betas.dasenII[,4], na.rm=T)->TEST4
density(betas.dasenII[,5], na.rm=T)->TEST5
density(betas.dasenII[,6], na.rm=T)->TEST6
density(betas.dasenII[,7], na.rm=T)->TEST7
density(betas.dasenII[,8], na.rm=T)->TEST8
density(betas.dasenII[,9], na.rm=T)->TEST9
density(betas.dasenII[,10], na.rm=T)->TEST10
density(betas.dasenII[,11], na.rm=T)->TEST11
density(betas.dasenII[,12], na.rm=T)->TEST12
density(betas.dasenII[,13], na.rm=T)->TEST13
density(betas.dasenII[,14], na.rm=T)->TEST14

#create plot
pdf("iPSCsThroughDiff/wo2samples/betas_density_plot_dasenII.pdf")
plot(TEST1, col="red", main = "Type II Probes after Normalisation", ylim = c(0, 6))
lines(TEST2, col="red")
lines(TEST3, col="purple")
lines(TEST4, col="purple")
lines(TEST5, col="purple")
lines(TEST6, col="purple")
lines(TEST7, col="blue")
lines(TEST8, col="blue")
lines(TEST9, col="blue")
lines(TEST10, col="blue")
lines(TEST11, col="green")
lines(TEST12, col="green")
lines(TEST13, col="green")
lines(TEST14, col="green")
legend("topleft", c("iPSC", "NPC", "Neuron 37", "Neuron 58"), col = c("red", "purple", "blue", "green"), pch = 18, title="Cell Stage")
dev.off()


###############################
#####Hemi Methylation##########
###############################

#Assess whether there is any difference in hemi-methylated probes between timepoints 

##subsetted out the beta values for non cpg probes and looked at dimension
betas.chr <- betas_d[grep("ch", rownames(betas_d)),]
dim(betas.chr)
[1] 2854   14

##created a subset of only cpg probes doing same thing as above using -grep, so everything but anything that starts with ch
betas.chr_WO <- betas_d[-grep("ch", rownames(betas_d)),]
dim(betas.chr_WO)
[1] 844249     14
dim(betas_d)
[1] 847103     14


##subset the non cpg probes by cell type or timepoint
betas.chr.1<-betas.chr[,c(1:2)]
betas.chr.2<-betas.chr[,c(3:6)]
betas.chr.3<-betas.chr[,c(7:10)]
betas.chr.4<-betas.chr[,c(11:14)]

#subset cpg probes by cell type
betas.chr_WO.1<-betas.chr_WO[,c(1:2)]
betas.chr_WO.2<-betas.chr_WO[,c(3:6)]
betas.chr_WO.3<-betas.chr_WO[,c(7:10)]
betas.chr_WO.4<-betas.chr_WO[,c(11:14)]

#pulled out specifically the hemi methylated non cpg probes, from 20% to 80% methylated
hemimeth_non.1<- which(rowMeans(betas.chr.1)>0.2 & rowMeans(betas.chr.1)<0.8)
length(hemimeth_non.1)
#2450
hemimeth_non.2<- which(rowMeans(betas.chr.2)>0.2 & rowMeans(betas.chr.2)<0.8)
length(hemimeth_non.2)
#28
hemimeth_non.3<- which(rowMeans(betas.chr.3)>0.2 & rowMeans(betas.chr.3)<0.8)
length(hemimeth_non.3)
#80
hemimeth_non.4<- which(rowMeans(betas.chr.4)>0.2 & rowMeans(betas.chr.4)<0.8)
length(hemimeth_non.4)
#85

#pulled out specifically the hemi methylated cpg probes, from 20% to 80% methylated
hemimeth_cpg.1<- which(rowMeans(betas.chr_WO.1)>0.2 & rowMeans(betas.chr_WO.1)<0.8)
length(hemimeth_cpg.1)
#448990
hemimeth_cpg.2<- which(rowMeans(betas.chr_WO.2)>0.2 & rowMeans(betas.chr_WO.2)<0.8)
length(hemimeth_cpg.2)
#400328
hemimeth_cpg.3<- which(rowMeans(betas.chr_WO.3)>0.2 & rowMeans(betas.chr_WO.3)<0.8)
length(hemimeth_cpg.3)
#405154
hemimeth_cpg.4<- which(rowMeans(betas.chr_WO.4)>0.2 & rowMeans(betas.chr_WO.4)<0.8)
length(hemimeth_cpg.4)
#393191


#Average each of the hemi methylation values for non-cpgs and cpgs for each cell type
hemimeth_non.1_beta <- betas.chr.1[hemimeth_non.1, ]
hemimeth_non.2_beta <- betas.chr.2[hemimeth_non.2, ]
hemimeth_non.3_beta <- betas.chr.3[hemimeth_non.3, ]
hemimeth_non.4_beta <- betas.chr.4[hemimeth_non.4, ]

hemimeth_cpg.1_beta <- betas.chr_WO.1[hemimeth_cpg.1, ]
hemimeth_cpg.2_beta <- betas.chr_WO.2[hemimeth_cpg.2, ]
hemimeth_cpg.3_beta <- betas.chr_WO.3[hemimeth_cpg.3, ]
hemimeth_cpg.4_beta <- betas.chr_WO.4[hemimeth_cpg.4, ]

hemi_non.1_mean <- mean(hemimeth_non.1_beta)
hemi_non.2_mean <- mean(hemimeth_non.2_beta)
hemi_non.3_mean <- mean(hemimeth_non.3_beta)
hemi_non.4_mean <- mean(hemimeth_non.4_beta)

hemi_cpg.1_mean <- mean(hemimeth_cpg.1_beta)
hemi_cpg.2_mean <- mean(hemimeth_cpg.2_beta)
hemi_cpg.3_mean <- mean(hemimeth_cpg.3_beta)
hemi_cpg.4_mean <- mean(hemimeth_cpg.4_beta)

###Plots####

#Plot showing Hemi−Methylation of Non−CpG Probes
pdf("iPSCsThroughDiff/without 2 sample/mean_noncpg_hemimet_wo2sample.pdf")
boxplot(rowMeans(hemimeth_non.1_beta),rowMeans(hemimeth_non.2_beta), rowMeans(hemimeth_non.3_beta), rowMeans(hemimeth_non.4_beta), main="Hemi-Methylation of Non-CpG Probes", xlab="Cell Type", ylab="Methylation Beta Value", names= c("iPSC", "NPC", "Neuron 37", "Neuron 58"), col=c("red", "purple", "blue", "green"))
dev.off()

#Plot showing Hemi−Methylation of Non−CpG Probes
pdf("iPSCsThroughDiff/without 2 sample/mean_cpg_hemimet_wo2sample.pdf")
boxplot(rowMeans(hemimeth_cpg.1_beta), rowMeans(hemimeth_cpg.2_beta), rowMeans(hemimeth_cpg.3_beta), rowMeans(hemimeth_cpg.4_beta), main="Hemi-Methylation of CpG Probes", xlab="Cell Type", ylab="Methylation Beta Value", names= c("iPSC", "NPC", "Neuron 37", "Neuron 58"), col=c("red", "purple", "blue", "green"))
dev.off()

##plot both of the above graphs together side by side on same graph
pdf("iPSCsThroughDiff/without 2 sample/mean_hemimet_wo2sample.pdf")
boxplot(rowMeans(hemimeth_non.1_beta),rowMeans(hemimeth_non.2_beta), rowMeans(hemimeth_non.3_beta), rowMeans(hemimeth_non.4_beta), rowMeans(hemimeth_cpg.1_beta), rowMeans(hemimeth_cpg.2_beta), rowMeans(hemimeth_cpg.3_beta), rowMeans(hemimeth_cpg.4_beta), main="Hemi-Methylation of CpG Probes", xlab="Cell Type", ylab="Methylation Beta Value", names= c("iPSC", "NPC", "Neuron 37", "Neuron 58", "iPSC", "NPC", "Neuron 37", "Neuron 58"), col=c("red", "purple", "blue", "green", "red", "purple", "blue", "green"))
dev.off()



########################################
#####Epigenetic Age Calculator##########
########################################

#Date: 14/11/19


#set working directory
setwd("")

#Load in required packages
library(minfi)
library(minfiData)
library(sva)
library(wateRmelon)
require(gdata)
library(methylumi)

#load in Pheno file and epic manifest
pheno<-read.csv ("Phenotypes.csv", stringsAsFactors = FALSE, header = TRUE) ##pheno file
pheno$Basename<-as.character(pheno$Basename) ##changing the Basenames to character strings
epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B2.csv", skip = 7) ##probe annotation file


#subset pheno file so it doesnt contain two failed samples
pheno<- pheno[c(1,4:16),]

#load in msetEPIC file
load("/mnt/data1/Jenny/EPIC QC/iPSCsThroughDiff/wo2samples/MsetEPIC_wo2samples.Rdata", verbose = FALSE)
epicManifest<-epicManifest[match(rownames(betas(msetEPIC)), epicManifest$Name),]

#apply pfilter to msetEPIC file
msetEPIC.pf <- pfilter(msetEPIC)

#normalise
EPIC.dasen<-dasen(msetEPIC.pf)

#Load in 27k annotation file 
##data for age calculator - online age calculator created by Steve Horvath to predict biol age based on methylation
probeAnnotation27k=read.csv("/mnt/data1/reference_files/AgeCalculator/datMiniAnnotation27k.csv")  ##read annotation file


##Remove cross-hyb/SNP/Bad probes

#load in txt file containing info on cross hybridising probes
cross<-read.table("/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
crosslist<-cross[,1]
length(crosslist)
#[1]44210  
#remove these cross-hyb probes from normalised object
data.1<-EPIC.dasen[ ! rownames(EPIC.dasen) %in% crosslist,]
dim(EPIC.dasen) 
#Features  Samples 
#847103       14 
dim(data.1)
#Features  Samples 
#803777       14

#load in txt file containing information on SNP probes
snpProbes<-read.table("/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
dim(snpProbes)
# 340327     18
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]
dim(snpProbes)
# 10888    18
SNPlist<-snpProbes[,1]
#remove SNP probes from data
data.2 <-data.1 [! rownames(data.1) %in% SNPlist,]

#load in txt file containing information on bad probes
badProbes<-read.csv("/mnt/data1/EPIC_reference/EPICArrayProbesToFilter.csv", stringsAsFactors = FALSE, header = TRUE)
dim(badProbes)
# 977    48
badlist<-badProbes[,1]
#remove bad probes from data
data.3 <- data.2[ ! rownames(data.2) %in% badlist,]
dim(data.1) #
#  803777     14
dim(data.2) #
#  794501     14
dim(data.3)
# 793851     14

#extract beta values from normalised dataset
betas.upload_dasenRM<-betas(data.3)

#pull out the 27k probes from the dataset
betas.upload_dasenRM<-betas.upload_dasenRM[match(probeAnnotation27k$Name, rownames(betas.upload_dasenRM)),] ##match your betas to annotation file
rownames(betas.upload_dasenRM)<-probeAnnotation27k$Name ##make row names of betas.upload same as rownames in Name of annotation file

#create a csv of this object, this is what you upload online 
write.csv(betas.upload_dasenRM, "iPSCsThroughDiff/wo2samples/Data_upload_DasenRM.csv")  ### download this file to your computer and compress

#create a simplified pheno file to upload with betas
#to contain information on gender of samples, a tissue column filled with NAs, the age of samples
pheno_tmp<-cbind(pheno$Basename, NA)
female.indictor<-rep(0, nrow(pheno))
#female.indictor[which(pheno$SEX == "male")]<-1 ##if genders provided this can be altered
tissue<-rep("NA", nrow(pheno_tmp))
pheno_tmp<-cbind(pheno_tmp, female.indictor, tissue)
colnames(pheno_tmp)<-c("SampleID", "Age", "Female", "Tissue")

#create csv of this
write.csv(pheno_tmp, "iPSCsThroughDiff/wo2samples/Pheno_AgeCalc_Noob.csv") ##update with where you want file to go and name

### upload the compressed betas file and this new phenotype file to http://dnamage.genetics.ucla.edu/new and save resulting output 


##Create plot to show epigenetic age associated with each timepoint or cell stage

###Plots

#load in output file from age calculator 
AgeCalc <- read.csv("Noob_Dasen_Age.csv")

#Make sure cell stage column is a factor so plots in the correct order not alphabetically 
AgeCalc$Cell.Stage <- factor(AgeCalc$Cell.Stage, levels = c("iPSC", "NPC", "Neuron-D37", "Neuron-D58"))

#create plot showing epigenetic age with each timepoint/cell stage
#add stripchart to show individual points on plot
pdf("Age Calc Noob and Dasen.pdf", width=12, height=6)
plot(AgeCalc$Cell.Stage,AgeCalc$Dasen_Age, ylim=c(-0.8, -0.4), ylab="Estimated Epigenetic Age (years)",xlab="Cell Stage",main="Biological Age Prediction using Dasen", col = c("#ff00ff", "#ffcc00", "#9650a1", "#45b449"))
stripchart(AgeCalc$Dasen_Age~AgeCalc$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.5)
dev.off()






