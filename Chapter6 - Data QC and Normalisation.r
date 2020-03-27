##Script to perform the quality control and normalisation for SV40 microglial cells that have been treated with LPS.
####Date:04/10/2019
##Author: Jennifer Imm

SV40 EPIC Analysis

#Rsession Information
R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-redhat-linux-gnu (64-bit)

#Set the working directory
setwd("/mnt/data1/Jenny/EPIC_QC/SV40 EPIC/")


#load required packages
library(methylumi)
library(wateRmelon)
require(gdata)

### load phenotype data, DNA methylation data and probe annotation information
pheno<-read.csv ("/mnt/data1/SV40_EPIC/SV40pheno.csv", stringsAsFactors = FALSE, header = TRUE) ##pheno file
epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B2.csv", skip = 7) ##probe annotation file


#Look at pheno file to check all is in order
head(pheno)
 Sample.ID Group Gender Sex Meth.Chip.ID Meth.Chip.Position
1   Control_1  ctrl   Male   0 202073200065             R01C01
2       LPS_1   lps   Male   0 202073200065             R02C01
3 LPS_Recov_1  lpsR   Male   0 202073200065             R03C01
4 LPS_Recov_2  lpsR   Male   0 202073200080             R01C01
5   Control_2  ctrl   Male   0 202073200080             R02C01
6       LPS_2   lps   Male   0 202073200080             R03C01
             Basename
1 202073200065_R01C01
2 202073200065_R02C01
3 202073200065_R03C01
4 202073200080_R01C01
5 202073200080_R02C01
6 202073200080_R03C01




#Set the idat file path so the HPC knows where to find them when they are called upon
idatPath<-c("/mnt/data1/Jenny/EPIC QC/SV40 EPIC/idats")

#Create an mset file which matches each idat file to the relevant sample, this is done by BeadChip barcode.
#This also converts the idats to methylumi objects so the fluorescence intensities can be interpreted in the later steps
msetEPIC <- readEPIC(idatPath=idatPath, barcodes = pheno$Basename, parallel = FALSE)



##Match the mset file just created to the EPIC manifest file so that we have all relevant info such as chromosomal location etc matched to our sample information
epicManifest<-epicManifest[match(rownames(betas(msetEPIC)), epicManifest$Name),]

#save output
save(msetEPIC, file = "MsetEPIC_SV40.Rdata")

#Code to reload next time
#load("/mnt/data1/Jenny/EPIC QC/SV40 EPIC/MsetEPIC_SV40.Rdata", verbose = FALSE)


######################################
####QC Pipeline#######################
######################################

#Step 1. Assessing the fluorescence intensities

#1a. Checking median intensities

### extract the median methylated and unmethylated intensities from the mset file 
m_intensities<-methylated(msetEPIC)
u_intensities<-unmethylated(msetEPIC)

#calculate the median methylated and unnmethylated intensities
M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

#output showing median methylated intensities 
202073200065_R01C01 202073200065_R02C01 202073200065_R03C01 202073200080_R01C01
               2565                1821                1271                2176
202073200080_R02C01 202073200080_R03C01 202073200082_R01C01 202073200082_R02C01
               1714                1601                1766                1564
202073200082_R03C01 202073200083_R01C01 202073200083_R02C01 202073200083_R03C01
               1694                1705                1861                2059

#output showing median unmethylated intensities 
202073200065_R01C01 202073200065_R02C01 202073200065_R03C01 202073200080_R01C01
               1860                1359                 997                1643
202073200080_R02C01 202073200080_R03C01 202073200082_R01C01 202073200082_R02C01
               1296                1193                1344                1221
202073200082_R03C01 202073200083_R01C01 202073200083_R02C01 202073200083_R03C01
               1281                1669                1778                1852



#Create a file to save all QC data to together with the pheno file
QCmetrics<-cbind(pheno, M.median, U.median)


#Create plots to visually assess median intensities

#Set an intensity cut off below which we should think about excluding samples.
intens.Thres<-1000 ## change this to adjust the threshold at which you filter 

#create plot without cut offs
pdf("Scatterplot_SampleIntensity.pdf")
par(mfrow = c(1,2))
hist(M.median, xlab = "Median Methylated intensity")
hist(U.median, xlab = "Median Unmetylated intensity")
par(mfrow = c(1,1))
plot(M.median, U.median, xlim=c(1000,2700), pch = 16, xlab = "Median Methylated Intensity", ylab = "Median Unmethylated Intensity", main="Median Sample Intensities")
text(M.median, U.median, labels=pheno$Sample.ID, cex=0.8, pos=1)
dev.off()

#create plot with cut offs
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

#Create plots to visually assess whether BeadChip or position on BeadChip affect median intensity
## plot each BeadChip as a boxplot
pdf("Boxplot_SampleIntensity_ByPlate.pdf")
par(mfrow = c(1,2))
par(mar = c(8, 4, 1, 1))
nCol<-length(unique(pheno$Meth.Chip.ID))## assumes there is a column called Plate in your phenotype file
boxplot(M.median ~ pheno$Meth.Chip.ID, ylab = "Median Methylated Intensity", xlab = "", las = 2, col = rainbow(nCol)) 
mtext(text="Chip ID", side=1, line=6.8) 
boxplot(U.median ~ pheno$Meth.Chip.ID, ylab = "Median Unmethylated Intensity", xlab = "", las = 2, col = rainbow(nCol))
mtext(text="Chip ID", side=1, line=6.8) 
dev.off()

## alternatively colour points in original scatterplot by BeadChip ID
pdf("Scatterplot_SampleIntensity_ColByPlate.pdf")
nCol<-length(unique(pheno$Meth.Chip.ID))## assumes there is a column called Plate in your phenotype file
plot(M.median, U.median, cex=2 pch = 16, xlab = "Median Methylated Intensity", ylab = "Median Unmethylated Intensity", ylim= c(0,3000), col = rainbow(nCol)[factor(pheno$Meth.Chip.ID)])
legend("topleft", levels(factor(pheno$Meth.Chip.ID)), col = rainbow(nCol), pch = 16, title="Chip ID", cex=1.2)
dev.off()

## plot each position on BeadChip as a boxplot
pdf("Boxplot_SampleIntensity_ByPos.pdf")
par(mfrow = c(1,2))
par(mar = c(8, 4, 1, 1))
nCol<-length(unique(pheno$Meth.Chip.Position))## assumes there is a column called Plate in your phenotype file
boxplot(M.median ~ pheno$Meth.Chip.Position, ylab = "Median Methylated Intensity", las = 2, col = rainbow(nCol), xlab="")
mtext(text="Position on Chip", side=1, line=5) 
boxplot(U.median ~ pheno$Meth.Chip.Position, ylab = "Median Unmethylated intensity", las = 2, col = rainbow(nCol), xlab="")
mtext(text="Position on Chip", side=1, line=5) 
dev.off()

#colour by position on chip as scatterplot
pdf("Scatterplot_SampleIntensity_ColByPos.pdf")
nCol<-length(unique(pheno$Meth.Chip.Position))
plot(M.median, U.median, cex=2, pch = 16, xlab = "Median Methylated Intensity", ylab = "Median Unmethylated Intensity", ylim= c(0,3000), col = rainbow(nCol)[factor(pheno$Meth.Chip.Position)])
legend("topleft", levels(factor(pheno$Meth.Chip.Position)), col = rainbow(nCol), pch = 16, title="Position on Chip")
dev.off()


####2. Check bisulfite conversion worked using fully methylated control probes. 

#There are a 8 fully methylated control probes included which can be used to generate a score to estimate the success of the bisulfite conversion step. As they are fully methylated these should have DNA methylation values of ~1. The bisulfite conversion score is essentially the median of these probes, and value < 80 is taken as a failure. 
#Use bs function to calculate bisulfite conversion efficiency of each sample
bs<-bscon(msetEPIC)

#create a histogram to show distribution of conversion efficiencies
pdf("HistogramBisulphiteConversionStatistics.pdf")
hist(bs, xlab = "Median % BS conversion", main = "Bisulphite Conversion Effiency", col=c("lightskyblue2","mediumpurple1"))
dev.off()

#Also create a table containing output
write.table(bs,file="SV40_bscon.txt",quote=F,col.names=F)

#Add this to QC output file 
QCmetrics<-cbind(QCmetrics, bs)


#Also worthwhile to have column detailing which treatment group each sample belongs to 
Treatment <- c("Control", "LPS", "Recovery", "Recovery", "Control", "LPS", "LPS", "Recovery", "Control", "Recovery", "LPS", "Control")
QCmetrics<-cbind(QCmetrics, Treatment)



#####3. beta density plot

#transpose beta values
t(betas)->betast

#make sure rows of pheno file are the Basename of each sample
rownames(pheno)<-pheno$Basename

#merge pheno file with transposed beta values
pheno_betas <- merge(betast, pheno, by = "row.names")

#Save output
write.csv(pheno_betas, "Pheno_betas.csv")

dim(pheno_betas)
# 12 866158


##create density plot of beta values 
##create density line for each sample individually using denisty() and then plot them all on one plot using lines function
plot(density(pheno_betas[,2:866158]), main = "Density plot of beta values", col = "blue")


density(betas[,1], na.rm=T)->TEST1
density(betas[,2], na.rm=T)->TEST2
density(betas[,3], na.rm=T)->TEST3
density(betas[,4], na.rm=T)->TEST4
density(betas[,5], na.rm=T)->TEST5
density(betas[,6], na.rm=T)->TEST6
density(betas[,7], na.rm=T)->TEST7
density(betas[,8], na.rm=T)->TEST8
density(betas[,9], na.rm=T)->TEST9
density(betas[,10], na.rm=T)->TEST10
density(betas[,11], na.rm=T)->TEST11
density(betas[,12], na.rm=T)->TEST12

pdf("betas_density_plot.pdf")
plot(TEST1, col="red", main = "Density Plot of Beta Values", ylim = c(0, 4.5))
lines(TEST2, col="blue")
lines(TEST3, col="green")
lines(TEST4, col="green")
lines(TEST5, col="red")
lines(TEST6, col="blue")
lines(TEST7, col="blue")
lines(TEST8, col="green")
lines(TEST9, col="red")
lines(TEST10, col="green")
lines(TEST11, col="blue")
lines(TEST12, col="red")

legend("topleft", c("Control", "LPS", "LPS + Recovery"), col = c("red", "blue", "green"), pch = 18, title="Treatment")
dev.off()

##create density plots for type I and II probes separately
##Pull out probe type information from mset file
data<- msetEPIC@featureData
onetwod <- data$DESIGN

proId <- which(onetwod == "I")
proIId <- which(onetwod == "II")


##subset data based on whether type I or II probe
EPIC.data_I  <- msetEPIC[proId,]
EPIC.data_II <- msetEPIC[proIId,]


##generate beta values 
betasI <- betas(EPIC.data_I)
betasII <- betas(EPIC.data_II)


#####Density plot of type I probes
t(betasI)->betast
rownames(pheno)<-pheno$Basename
pheno_betas <- merge(betast, pheno, by = "row.names")

density(betasI[,1], na.rm=T)->TEST1
density(betasI[,2], na.rm=T)->TEST2
density(betasI[,3], na.rm=T)->TEST3
density(betasI[,4], na.rm=T)->TEST4
density(betasI[,5], na.rm=T)->TEST5
density(betasI[,6], na.rm=T)->TEST6
density(betasI[,7], na.rm=T)->TEST7
density(betasI[,8], na.rm=T)->TEST8
density(betasI[,9], na.rm=T)->TEST9
density(betasI[,10], na.rm=T)->TEST10
density(betasI[,11], na.rm=T)->TEST11
density(betasI[,12], na.rm=T)->TEST12

pdf("betas_density_I.pdf")
plot(TEST1, col="red", main = "Type I Probes Before Normalisation", ylim = c(0, 5.5), xlab="Beta Value")
lines(TEST2, col="blue")
lines(TEST3, col="green")
lines(TEST4, col="green")
lines(TEST5, col="red")
lines(TEST6, col="blue")
lines(TEST7, col="blue")
lines(TEST8, col="green")
lines(TEST9, col="red")
lines(TEST10, col="green")
lines(TEST11, col="blue")
lines(TEST12, col="red")

legend("topleft", c("Control", "LPS", "LPS + Recovery"), col = c("red", "blue", "green"), pch = 18, title="Treatment")
dev.off()



#####Density plot of type II probes
t(betasII)->betast
rownames(pheno)<-pheno$Basename
pheno_betas <- merge(betast, pheno, by = "row.names")

plot(density(pheno_betas[,2:723956]), main = "Type II Probes Before Normalisation", col = "blue")

density(betasII[,1], na.rm=T)->TEST1
density(betasII[,2], na.rm=T)->TEST2
density(betasII[,3], na.rm=T)->TEST3
density(betasII[,4], na.rm=T)->TEST4
density(betasII[,5], na.rm=T)->TEST5
density(betasII[,6], na.rm=T)->TEST6
density(betasII[,7], na.rm=T)->TEST7
density(betasII[,8], na.rm=T)->TEST8
density(betasII[,9], na.rm=T)->TEST9
density(betasII[,10], na.rm=T)->TEST10
density(betasII[,11], na.rm=T)->TEST11
density(betasII[,12], na.rm=T)->TEST12

pdf("betas_density_II.pdf")
plot(TEST1, col="red", main = "Type II Probes Before Normalisation", ylim = c(0, 4.5), xlab="Beta Value")
lines(TEST2, col="blue")
lines(TEST3, col="green")
lines(TEST4, col="green")
lines(TEST5, col="red")
lines(TEST6, col="blue")
lines(TEST7, col="blue")
lines(TEST8, col="green")
lines(TEST9, col="red")
lines(TEST10, col="green")
lines(TEST11, col="blue")
lines(TEST12, col="red")

legend("topleft", c("Control", "LPS", "LPS + Recovery"), col = c("red", "blue", "green"), pch = 18, title="Treatment")
dev.off()


###4. Calculate mitotic age

#load in source information containing list of probes needed by algorithm
source("/mnt/data1/reference_files/MiAge/function_library.r") ### library of all functions used in the calculation of mitotic ages                                                                                 


clocksitesID=as.matrix(read.csv("/mnt/data1/reference_files/MiAge/Additional_File1.csv",header=T))[,1]   ### epigenetic clock CpG sites

#ensure betas are in correct format
betas2 <- betas(msetEPIC)

##select and pull out CpG sites needed to run the clock
beta=  betas2[ match(clocksitesID,rownames(betas2)),]##select clock CpG sites 

#check you have pulled out the right number of probes
dim(beta)
     268       12

#Load in site specific parameters file
load("/mnt/data1/reference_files/MiAge/site_specific_parameters.Rdata") #### site-specific parameters

#Run the mitotic age clock
b=methyl.age[[1]];c=methyl.age[[2]];d=methyl.age[[3]]

n=mitotic.age(beta,b,c,d) ### estimated mitotic age
names(n)=colnames(beta)
write.table(n,file="SV40_mitotic_age.txt",quote=F,col.names=F)

#put this output in QC output file
QCmetrics=cbind(QCmetrics,n)

QCmetrics$Treatment=as.factor(QCmetrics$Treatment)

#create plot showing mitotic age split by treatment
pdf("SV40_mitotic_age_prediction.pdf")
plot(QCmetrics$Treatment,QCmetrics$n,ylab="Estimated no. of Cell Divisions",xlab="Treatment",main="Mitotic Age Prediction", col = c("red", "blue", "green"))
dev.off()


###5. P filter

#default perc output 
pfilter(msetEPIC, perc=1)

#Output
0 samples having 1 % of sites with a detection p-value greater than 0.05 were removed
Samples removed:
13207 sites were removed as beadcount <3 in 5 % of samples
9180 sites having 1 % of samples with a detection p-value greater than 0.05 were removed



msetEPIC.pf <- pfilter(msetEPIC)
pFilterPass<-colnames(betas(msetEPIC)) %in% colnames(betas(msetEPIC.pf))

QCmetrics<-cbind(QCmetrics, pFilterPass)


#Create csv file containing QC information 
write.csv(QCmetrics, file="QC_06062019_JI.csv")

##save pfiltered betas
save(msetEPIC.pf, file = "QC'd betas_06062019_JI.rdata")


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

#Calculate median methylated and unmethylated signal intensities for each sample
M.median_dasen<-apply(m_intensities_dasen, 2, median)
U.median_dasen<-apply(u_intensities_dasen, 2, median)

#Median methylated intensity output
M.median_dasen
202073200065_R01C01 202073200065_R02C01 202073200065_R03C01 202073200080_R01C01
           1863.891            1863.583            1863.833            1863.667
202073200080_R02C01 202073200080_R03C01 202073200082_R01C01 202073200082_R02C01
           1863.583            1863.583            1863.583            1864.432
202073200082_R03C01 202073200083_R01C01 202073200083_R02C01 202073200083_R03C01
           1863.583            1863.583            1864.250            1863.583


#Median unmethylated intensity output
U.median_dasen
202073200065_R01C01 202073200065_R02C01 202073200065_R03C01 202073200080_R01C01
           1482.333            1482.333            1482.833            1482.000
202073200080_R02C01 202073200080_R03C01 202073200082_R01C01 202073200082_R02C01
           1482.083            1482.667            1481.958            1482.750
202073200082_R03C01 202073200083_R01C01 202073200083_R02C01 202073200083_R03C01
           1482.417            1482.500            1482.467            1482.500
		   


#Create scatterplots and histograms to look at median intensities
pdf("Scatterplot_MEDIANIntensities_EPIC.Dasen.pdf")		   
par(mfrow = c(1,2))
hist(M.median_dasen, xlab = "Median M intensity")
hist(U.median_dasen, xlab = "Median U intensity")
par(mfrow = c(1,1))
plot(M.median_dasen, U.median_dasen, xlim=c(1863.4,1864.6) , ylim=c(1481.8,1483), pch = 16, xlab = "Median Methylated Intensity", ylab = "Median Unmethylated Intensity", main="Dasen Normalised Sample Intensities")
text(M.median_dasen, U.median_dasen, labels=pheno$Sample.ID, cex=0.8, pos=1)
dev.off()


##Create beta density plot to look at distribution of beta values
#transpose data
t(betas_d)->betast
rownames(pheno)<-pheno$Basename

#merge with phenotype file
pheno_betas_dasen <- merge(betast, pheno, by = "row.names")  

dim(pheno_betas_dasen)
12 844194

##Create density plot using density for each sample separately and then plot each sample using lines function
plot(density(pheno_betas_dasen[,2:844194]), main = "Density plot of beta values", col = "blue")


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

pdf("betas_density_plot_dasen.pdf")
plot(TEST1, col="red", main = "Density Plot of Dasen Normalised Beta Values", ylim = c(0, 4.5), xlab="Beta Value")
lines(TEST2, col="blue")
lines(TEST3, col="green")
lines(TEST4, col="green")
lines(TEST5, col="red")
lines(TEST6, col="blue")
lines(TEST7, col="blue")
lines(TEST8, col="green")
lines(TEST9, col="red")
lines(TEST10, col="green")
lines(TEST11, col="blue")
lines(TEST12, col="red")

legend("topleft", c("Control", "LPS", "LPS + Recovery"), col = c("red", "blue", "green"), pch = 18, title="Treatment")
dev.off()


###Pull out type I and type II probes
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
write.csv(betas.dasen, "Dasen normalised data.csv")
write.csv(betas.dasenI, "Dasen normalised data type I probes.csv")
write.csv(betas.dasenII, "Dasen normalised data type II probes.csv")

save(betas.dasen, file = "Dasen normalised data.rdata")
save(betas.dasenI, file = "Dasen normalised data type I probes.rdata")
save(betas.dasenII, file = "Dasen normalised data type II probes.rdata")


#####Density plot of type I probes
t(betas.dasenI)->betast
rownames(pheno)<-pheno$Basename
pheno_betas_dasen <- merge(betast, pheno, by = "row.names")

dim(pheno_betas_dasen)
#12 138449

plot(density(pheno_betas_dasen[,2:138449]), main = "Type I Probes After Normalisation", col = "blue")

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

pdf("betas_density_plot_dasenI.pdf")
plot(TEST1, col="red", main = "Type I Probes after Normalisation", ylim = c(0, 5.5), xlab="Beta Value")
lines(TEST2, col="blue")
lines(TEST3, col="green")
lines(TEST4, col="green")
lines(TEST5, col="red")
lines(TEST6, col="blue")
lines(TEST7, col="blue")
lines(TEST8, col="green")
lines(TEST9, col="red")
lines(TEST10, col="green")
lines(TEST11, col="blue")
lines(TEST12, col="red")

legend("topleft", c("Control", "LPS", "LPS + Recovery"), col = c("red", "blue", "green"), pch = 18, title="Treatment")
dev.off()




#####Density plot of type II probes
t(betas.dasenII)->betast
rownames(pheno)<-pheno$Basename
pheno_betas_dasen <- merge(betast, pheno, by = "row.names")

dim(pheno_betas_dasen)
#12 705754

plot(density(pheno_betas_dasen[,2:705754]), main = "Type II Probes After Normalisation", col = "blue")

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

pdf("betas_density_plot_dasenII.pdf")
plot(TEST1, col="red", main = "Type II Probes After Normalisation", ylim = c(0, 4.5), xlab="Beta Value")
lines(TEST2, col="blue")
lines(TEST3, col="green")
lines(TEST4, col="green")
lines(TEST5, col="red")
lines(TEST6, col="blue")
lines(TEST7, col="blue")
lines(TEST8, col="green")
lines(TEST9, col="red")
lines(TEST10, col="green")
lines(TEST11, col="blue")
lines(TEST12, col="red")

legend("topleft", c("Control", "LPS", "LPS + Recovery"), col = c("red", "blue", "green"), pch = 18, title="Treatment")
dev.off()



###############################
#####Hemi Methylation##########
###############################

#Assess whether there is any difference in hemi-methylated probes between treatments 

##subsetted out the beta values for non cpg probes and looked at dimension
betas.chr <- betas_d[grep("ch", rownames(betas_d)),]
dim(betas.chr)
[1] 2809   12

##created a subset of only cpg probes doing same thing as above using -grep, so everything but anything that starts with ch
betas.chr_WO <- betas_d[-grep("ch", rownames(betas_d)),]
dim(betas.chr_WO)
[1] 841376     12
dim(betas_d)
[1] 844185     12


colnames(betas_d)
 [1] "202073200065_R01C01" "202073200065_R02C01" "202073200065_R03C01"
 [4] "202073200080_R01C01" "202073200080_R02C01" "202073200080_R03C01"
 [7] "202073200082_R01C01" "202073200082_R02C01" "202073200082_R03C01"
[10] "202073200083_R01C01" "202073200083_R02C01" "202073200083_R03C01"


##subset the non cpg probes by treatment
betas.chr.1<-betas.chr[,c(1, 5, 9, 12)]
betas.chr.2<-betas.chr[,c(2, 6, 7, 11)]
betas.chr.3<-betas.chr[,c(3, 4, 8, 10)]

#subset cpg probes by treatment
betas.chr_WO.1<-betas.chr_WO[,c(1, 5, 9, 12)]
betas.chr_WO.2<-betas.chr_WO[,c(2, 6, 7, 11)]
betas.chr_WO.3<-betas.chr_WO[,c(3, 4, 8, 10)]


#pulled out specifically the hemi methylated non cpg probes, from 20% to 80% methylated
hemimeth_non.1<- which(rowMeans(betas.chr.1)>0.2 & rowMeans(betas.chr.1)<0.8)
length(hemimeth_non.1)
#20
hemimeth_non.2<- which(rowMeans(betas.chr.2)>0.2 & rowMeans(betas.chr.2)<0.8)
length(hemimeth_non.2)
#31
hemimeth_non.3<- which(rowMeans(betas.chr.3)>0.2 & rowMeans(betas.chr.3)<0.8)
length(hemimeth_non.3)
#29


#pulled out specifically the hemi methylated cpg probes, from 20% to 80% methylated
hemimeth_cpg.1<- which(rowMeans(betas.chr_WO.1)>0.2 & rowMeans(betas.chr_WO.1)<0.8)
length(hemimeth_cpg.1)
#427508
hemimeth_cpg.2<- which(rowMeans(betas.chr_WO.2)>0.2 & rowMeans(betas.chr_WO.2)<0.8)
length(hemimeth_cpg.2)
#439563
hemimeth_cpg.3<- which(rowMeans(betas.chr_WO.3)>0.2 & rowMeans(betas.chr_WO.3)<0.8)
length(hemimeth_cpg.3)
#440051



#Average each of the hemi methylation values for each treatment
hemimeth_non.1_beta <- betas.chr.1[hemimeth_non.1, ]
hemimeth_non.2_beta <- betas.chr.2[hemimeth_non.2, ]
hemimeth_non.3_beta <- betas.chr.3[hemimeth_non.3, ]


hemimeth_cpg.1_beta <- betas.chr_WO.1[hemimeth_cpg.1, ]
hemimeth_cpg.2_beta <- betas.chr_WO.2[hemimeth_cpg.2, ]
hemimeth_cpg.3_beta <- betas.chr_WO.3[hemimeth_cpg.3, ]


hemi_non.1_mean <- mean(hemimeth_non.1_beta)
hemi_non.2_mean <- mean(hemimeth_non.2_beta)
hemi_non.3_mean <- mean(hemimeth_non.3_beta)


hemi_cpg.1_mean <- mean(hemimeth_cpg.1_beta)
hemi_cpg.2_mean <- mean(hemimeth_cpg.2_beta)
hemi_cpg.3_mean <- mean(hemimeth_cpg.3_beta)



#plot as box plot
#Plot showing Hemi−Methylation of Non−CpG Probes
pdf("mean_noncpg_hemimet.pdf")
boxplot(rowMeans(hemimeth_non.1_beta),rowMeans(hemimeth_non.2_beta), rowMeans(hemimeth_non.3_beta), main="Hemi-Methylation of Non-CpG Probes", xlab="Treatment", ylab="Methylation Beta Value", names= c("Control", "LPS", "LPS + Recovery"), col=c("red", "blue", "green"))
dev.off()

#Plot showing Hemi−Methylation of CpG Probes
pdf("mean_cpg_hemimet.pdf")
boxplot(rowMeans(hemimeth_cpg.1_beta), rowMeans(hemimeth_cpg.2_beta), rowMeans(hemimeth_cpg.3_beta), main="Hemi-Methylation of CpG Probes", xlab="Treatment", ylab="Methylation Beta Value", names= c("Control", "LPS", "Recovery"), col=c("red", "blue", "green"))
dev.off()

##plot both of the above graphs together side by side on same graph
pdf("mean_hemimet.pdf")
boxplot(rowMeans(hemimeth_non.1_beta),rowMeans(hemimeth_non.2_beta), rowMeans(hemimeth_non.3_beta), rowMeans(hemimeth_cpg.1_beta), rowMeans(hemimeth_cpg.2_beta), rowMeans(hemimeth_cpg.3_beta), main="Hemi-Methylation of Probes", xlab="Treatment", ylab="Methylation Beta Value", names= c("Control", "LPS", "Recovery", "Control", "LPS", "Recovery"), col=c("red", "blue", "green", "red", "blue", "green"))
dev.off()




##Before we can analyse data need to remove cross hybridising, bad and SNP probes as outlined in McCartney et al.



############################################################
#REMOVE CROSS-HYBRIDISING PROBES
############################################################


betas.dasen<-read.csv("/mnt/data1/Jenny/EPIC QC/SV40 EPIC/Dasen normalised data.csv", row.names=1)
dim(data)
#844185     12

cross<-read.table("/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
crosslist<-cross[,1]
length(crosslist)
#[1]44210  
data.1<-betas.dasen[ ! rownames(betas.dasen) %in% crosslist,]
dim(betas.dasen)  	# 844185     12
dim(data.1) 	# 800932     12

############################################################
##REMOVE PROBES with SNPs
############################################################

snpProbes<-read.table("/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
dim(snpProbes)
# 340327     18
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]
dim(snpProbes)
# 10888    18
SNPlist<-snpProbes[,1]
data.2 <-data.1 [! rownames(data.1) %in% SNPlist,]
dim(data.2) # 791748     12


############################################################
##REMOVE bad probes
############################################################

badProbes<-read.csv("/mnt/data1/EPIC_reference/EPICArrayProbesToFilter.csv", stringsAsFactors = FALSE, header = TRUE)
dim(badProbes)
badlist<-badProbes[,1]
data.3 <- data.2[ ! rownames(data.2) %in% badlist,]

dim(data.3) 
#  791549     12    



##Save output to take on for analysis
write.csv(data.3, "DasenNormData_badprobesrm.csv")
save.image(file="DasenNormData_badprobesrm.rdata")








