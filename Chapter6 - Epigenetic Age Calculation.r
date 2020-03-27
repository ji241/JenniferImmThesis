##Script to extract probes necessary to run epigenetic age calculator and the visualisation of this data
####Date:18/12/2019
##Author: Jennifer Imm


##Normalisse data using Dasen before running epigen age calculator and using the latest version of Horvath clock

##load in packages
library(methylumi)
library(wateRmelon)
require(gdata)

##load in normalised data and EPIC manifest 

#EPIC manifest file
epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B2.csv", skip = 7) ##probe annotation file

#load in mset file
load("/mnt/data1/Jenny/EPIC QC/SV40 EPIC/MsetEPIC_SV40.Rdata", verbose = FALSE)

#27K array annotation file
probeAnnotation27k=read.csv("/mnt/data1/reference_files/AgeCalculator/datMiniAnnotation27k.csv")  ##read annotation file

#Match the epic manifest to the mset file
epicManifest<-epicManifest[match(rownames(betas(msetEPIC)), epicManifest$Name),]

#run pfilter to make sure failing probes are removed
msetEPIC.pf <- pfilter(msetEPIC)

#normalise data
EPIC.dasen<-dasen(msetEPIC.pf)

##remove cross-hyb, SNP and bad probes

cross<-read.table("/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
crosslist<-cross[,1]
length(crosslist)
#[1]44210  
data.1<-EPIC.dasen[ ! rownames(EPIC.dasen) %in% crosslist,]
dim(EPIC.dasen) 
#Features  Samples 
#844185       12 
dim(data.1)
#Features  Samples 
#800932       12 

snpProbes<-read.table("/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
dim(snpProbes)
# 340327     18
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]
dim(snpProbes)
# 10888    18
SNPlist<-snpProbes[,1]
data.2 <-data.1 [! rownames(data.1) %in% SNPlist,]

badProbes<-read.csv("/mnt/data1/EPIC_reference/EPICArrayProbesToFilter.csv", stringsAsFactors = FALSE, header = TRUE)
dim(badProbes)
# 977    48
badlist<-badProbes[,1]
data.3 <- data.2[ ! rownames(data.2) %in% badlist,]
dim(data.1) #
#Features  Samples 
#800932       12
dim(data.2) #
#  791748       12
dim(data.3)
# 7791549       12 


#rename data file
betas.upload_dasen<-betas(data.3)

#exatract the 27K array CpGs Gs from the dataset, these are the ones needed to be uploaded online
betas.upload_data<-betas.upload_dasen[match(probeAnnotation27k$Name, rownames(betas.upload_dasenRM)),]

write.csv(betas.upload_data, "Data_upload_DasenRM2.csv")  ### download this file to your computer and compress

##Upload this file to http://dnamage.genetics.ucla.edu/new and save resulting output 


##Create plot to show epigenetic age associated with each treatment

#Load in the epigenetic age output 
EpiAge<- read.csv("Age Calculator output.csv", header = T)

#Subset output to contain only columns of interest
EpiAge<-EpiAge[,c(1:2)]

#ensure pheno and epi age file have same rownames
rownames(pheno)<-pheno$Basename
rownames(EpiAge)<- EpiAge$SampleID

#Create new object that has matched the phenotype file to the epigenetic age output file
DNAmAge<-cbind(pheno, EpiAge$DNAmAge[match(rownames(pheno),rownames(EpiAge))])

#Assign column names
colnames(DNAmAge)<-c("Sample.ID", "Gender", "Sex", "Meth.Chip.ID", "Meth.Chip.Position", "Basename", "treatment", "treatment_cod", "MetAge")

#Change treatment column to be a factor containing three levels
DNAmAge$treatment <- factor(DNAmAge$treatment, levels = c("Control", "LPS", "Recovery"))


##Create plots showing epigenetic age with different treatments
pdf("Epigenetic Age after Dasen.pdf", width=12, height=5)
par(mfrow = c(1,2))
hist(DNAmAge$MetAge, col=c("#ffa500","#005aff"), main="Estimated Epigenetic Age", xlab="Estimated Epigenetic Age (years)")
plot(DNAmAge$treatment,DNAmAge$MetAge, col=c("red","blue","green"), ylab="Estimated Epigenetic Age (years)", xlab="Treatment", main="Epigenetic Age Prediction")
dev.off()

#create output file
write.csv(DNAmAge, "Pheno+DNAmAge.csv") 

#Perform ANOVA and post-hoc Tukey to see if there is any difference in the epigenetic ages between treatments. 
res.aov<- aov(DNAmAge$MetAge~DNAmAge$treatment)
summary(res.aov)
#                   Df Sum Sq Mean Sq F value Pr(>F)
#DNAmAge$treatment  2   65.0   32.50   1.218   0.34
#Residuals          9  240.2   26.68  


TukeyHSD(res.aov)
#                      diff        lwr       upr     p adj
#LPS-Control      -5.2318699 -15.429987  4.966247 0.3659213
#Recovery-Control -4.5766659 -14.774783  5.621452 0.4541852
#Recovery-LPS      0.6552041  -9.542913 10.853321 0.9824520



