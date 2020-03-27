##Script to extract probes necessary to run epigenetic age calculator and the visualisation of this data
####Date:04/12/2019
##Author: Jennifer Imm and Ehsan Pishva

#Load in and format data
betas<-read.csv("/mnt/data1/Jenny/EPIC QC/SV40 EPIC/DasenNormData_badprobesrm.csv", sep=",", header=T)
rownames(betas)<-betas[,1]
betas[,1]<-NULL

#load in EPIC manifest
epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7)
rownames(epicManifest)<-epicManifest$IlmnID

#match manifest to betas 
x<-intersect(rownames(betas), rownames(epicManifest))

#subset manifest and betas using matched betas/manifest object
annot<-epicManifest[x,]
betas<-betas[x,]
identical(rownames(betas), rownames(annot))
#[1] TRUE

#load in phenotype file
pheno<-read.csv("/mnt/data1/Jenny/EPIC QC/SV40 EPIC/SV40Pheno.csv")

#change colnames of betas to match phenotype file
colnames(betas)<-substr(colnames(betas),2,20)
betas<-betas[,as.character(pheno$Basename)]

identical(colnames(betas),as.character(pheno$Basename))
#[1] TRUE

#calculate median absolute derivation (ie mad) and pull out top 20% most variable probes
mad.jenny<- (apply(betas, 1,mad))
a<-quantile(mad.jenny,0.8)
var_probes<-which(mad.jenny > a) ##taking the top 20% most variable probes

#Label the three treatments from 1-3
pheno$condition<-NULL
pheno[which((substr(pheno$Sample.ID,1,7))=="Control"),"condition"]=1 ## giving the various conditions numbers 1, 2 and 3
pheno[which((substr(pheno$Sample.ID,1,3))=="LPS"),"condition"]=2
pheno[which((substr(pheno$Sample.ID,1,9))=="LPS_Recov"),"condition"]=3

#save file
save(betas,pheno,var_probes,file="Jenny_SV40.Rdata") ##file with most variable probes, all betas and pheno in

####################################################
###perform epigenome wide associaion study (EWAS)###
####################################################

#Assign data object to contain most variable probes in control and LPS conditions
data<-betas[var_probes,which(pheno$condition==1|pheno$condition==2)]

#define condition as Control and LPS
condition<-as.factor(pheno[which(pheno$condition==1|pheno$condition==2),]$condition)

#This is the function
EWAS <- function(x, condition){
  try(fit<-lm(x ~ condition), silent = T)
  if(inherits(fit,'try-error')) return(rep(NA,4))
  return(coef(summary(fit))[2,])
} 

#define how much of HPC you want to use (ie how many clusters)
cl<- makeCluster(32)

#run the EWAS for control vs LPS and save to results object
results<-t(parApply(cl,data , 1, EWAS, condition)) 

#name column of output (estimate, standard error, t-value and p-value)
colnames(results) <- c("Est2vs1","SE2vs1","t2vs1","P2vs1") 

#change to data frame
res12<-data.frame(results)

#pull out probes that have p-value of less than 1e-3
length(which(res12$P2vs1< 1.e-3)) ##arbitrary p value
#pull out probe IDs of these significant probes
x.1<-rownames(res12[which(res12$P2vs1< 1.e-3),])

##Define betas and condition for LPS vs LPS+recovery EWAS
data<-betas[x.1,which(pheno$condition==2|pheno$condition==3)]
condition<-as.factor(pheno[which(pheno$condition==2|pheno$condition==3),]$condition)

#Run LPS vs LPS+recovery EWAS
results<-t(parApply(cl,data , 1, EWAS, condition))
head(results)
#rename output columns
colnames(results) <- c("Est3vs2","SE3vs2","t3vs2","P3vs2")
res23<-data.frame(results)

Define betas and condition for Control vs LPS+recovery EWAS
data<-betas[x.1,which(pheno$condition==1|pheno$condition==3)]
condition<-as.factor(pheno[which(pheno$condition==1|pheno$condition==3),]$condition)

#Run EWAS for control vs LPS+recovery
results<-t(parApply(cl,data , 1, EWAS, condition))
head(results)
#rename columns in output
colnames(results) <- c("Est3vs1","SE3vs1","t3vs1","P3vs1")
res13<-data.frame(results)

#combine the results of the three EWAS'
res<-cbind(res12[x.1,], res23[x.1,], res13[x.1,]) ##results of 3 EWAS

#save file and stop clusters
save(res, file="res.SV40.rda")
stopCluster(cl)

#Pull out rows in manifest that are in the results file
x<-intersect(rownames(res), rownames(epicManifest))
annot<-epicManifest[x,]
res<-res[x,]
identical(rownames(res), rownames(annot))

#bind together results and annotation and save
res.annot<-cbind(res, annot)
save(res.annot, file="res.annot.SV40.rda")


####################################
###Define groups and create plots###
####################################


#Date: 12.12.2019


##Load in data
load("Jenny_SV40.Rdata")
load("res.annot.SV40.rda")
load("res.SV40.rda")
load("SV40.gam.rdata")


#Add treatment column to pheno file
pheno$treatment <- c("Control","Control","Control","Control", "LPS", "LPS", "LPS", "LPS", "Recovery", "Recovery", "Recovery", "Recovery")

#Define four treatment groups

#eg hypo_short pull out that are signif between control and LPS (ie P2vs1<1.e-3), that are becoming less methylated (estimate < 0),
#that have increasing methylation between LPS and recovery (estimate >0) and also are not significant between control and recovery (ie P3vs1 > 0.05)
hypo_short<-res[which(res$P2vs1<1.e-3 & res$Est2vs1 < 0 &  res$Est3vs2 > 0 & res$P3vs1 > 0.05),]
#order by p-value
hypo_short_ord <- hypo_short[order(hypo_short$P2vs1),]
#pull out top ten DMPs
hypo_short_ord[1:10,]
#subset dat so only contains hypo_short probes
hypo_short_plot <- dat[rownames(hypo_short_ord),]

#Do similar for other three groupings
hypo_long<-res[which(res$P2vs1<1.e-3 & res$Est2vs1 < 0 &  res$Est3vs2 < 0 & res$P3vs1 < 0.05),]
hypo_long_ord <- hypo_long[order(hypo_long$P2vs1),]
hypo_long_ord[1:8,]
hypo_long_plot <- dat[rownames(hypo_long_ord),]

hyper_short<-res[which(res$P2vs1<1.e-3 & res$Est2vs1 > 0 &  res$Est3vs2 < 0 & res$P3vs1 > 0.05),]
hyper_short_ord <- hyper_short[order(hyper_short$P2vs1),]
hyper_short_ord[1:10,]
hyper_short_plot <- dat[rownames(hyper_short_ord),]

hyper_long<-res[which(res$P2vs1<1.e-3 & res$Est2vs1 > 0 &  res$Est3vs2 > 0 & res$P3vs1 < 0.05),]
hyper_long_plot <- dat[rownames(hyper_long),]


###Calculate delta values for 1vs2 (control vs LPS) and 1vs3 (control and LPS+recovery)

res$D2vs1 <- res$Est2vs1*100
res$D3vs1 <- res$Est3vs1*100

#Pull out rows that are specific to each methylation change category
D_hyposhort <- res[rownames(hypo_short_plot),]
D_hypolong <- res[rownames(hypo_long_plot),]
D_hypershort<- res[rownames(hyper_short_plot),]
D_hyperlong<- res[rownames(hyper_long_plot),]

#Create csv files of files
write.csv(D_hyposhort, file="res_HypoShort.csv")
write.csv(D_hypolong, file="res_HypoLong.csv")
write.csv(D_hypershort, file="res_HyperShort.csv")
write.csv(D_hyperlong, file="res_HyperLong.csv")

#########
##PLOTS##
#########

#plot all short term hypomethylating probes
pdf("hypo_short.pdf", height=10, width=10)
par(mfrow=c(4,4))
for(i in 1:97){
  boxplot(as.numeric(hypo_short_plot[i,])~pheno$treatment, ylab="Beta Value", xlab="Treatment", main=rownames(hypo_short_plot[i,]), col=c("#ebebff", "#8989ff", "#4e4eff"))
  stripchart(as.numeric(hypo_short_plot[i,])~pheno$treatment, add=TRUE, vertical=TRUE, method="jitter", pch=16)
}
dev.off()

#plot top 4 probes
pdf("hypo_short_top4.pdf", height=10, width=10)
par(mfrow=c(2,2))
for(i in 1:4){
  boxplot(as.numeric(hypo_short_plot[i,])~pheno$treatment, ylab="Beta Value", xlab="Treatment", main=rownames(hypo_short_plot[i,]), col=c("#ebebff", "#8989ff", "#4e4eff"))
  stripchart(as.numeric(hypo_short_plot[i,])~pheno$treatment, add=TRUE, vertical=TRUE, method="jitter", pch=16, cex=1.5)
}
dev.off()

#plot all long term hypomethylating probes
pdf("hypo_long.pdf", height=10, width=10)
par(mfrow=c(3,3))
for(i in 1:8){
  boxplot(as.numeric(hypo_long_plot[i,])~pheno$treatment, ylab="Beta Value", xlab="Treatment", main=rownames(hypo_long_plot[i,]), col=c("#1e8bb7", "#135873", "#0e3e52"))
  stripchart(as.numeric(hypo_long_plot[i,])~pheno$treatment, add=TRUE, vertical=TRUE, method="jitter", pch=16)
}
dev.off()

#plot top 4 long term hypomethylating probes
pdf("hypo_long_top4.pdf", height=10, width=10)
par(mfrow=c(2,2))
for(i in 1:2){
  boxplot(as.numeric(hypo_long_plot[i,])~pheno$treatment, ylab="Beta Value", xlab="Treatment", main=rownames(hypo_long_plot[i,]), col=c("#1e8bb7", "#135873", "#0e3e52"))
  stripchart(as.numeric(hypo_long_plot[i,])~pheno$treatment, add=TRUE, vertical=TRUE, method="jitter", pch=16, cex=1.5)
}
dev.off()


#plot all short term hypermethylating probes
pdf("hyper_short.pdf", height=10, width=10)
par(mfrow=c(4,4))
for(i in 1:67){
  boxplot(as.numeric(hyper_short_plot[i,])~pheno$treatment, ylab="Beta Value", xlab="Treatment", main=rownames(hyper_short_plot[i,]), col=c("#ffebeb", "#ff8989", "#ff4e4e"))
  stripchart(as.numeric(hyper_short_plot[i,])~pheno$treatment, add=TRUE, vertical=TRUE, method="jitter", pch=16)
}
dev.off()

#plot top 4
pdf("hyper_short_top4.pdf", height=10, width=10)
par(mfrow=c(2,2))
for(i in 1:4){
  boxplot(as.numeric(hyper_short_plot[i,])~pheno$treatment, ylab="Beta Value", xlab="Treatment", main=rownames(hyper_short_plot[i,]), col=c("#ffebeb", "#ff8989", "#ff4e4e"))
  stripchart(as.numeric(hyper_short_plot[i,])~pheno$treatment, add=TRUE, vertical=TRUE, method="jitter", pch=16,cex=1.5)
}
dev.off()


pdf("hyper_long.pdf")
  boxplot(as.numeric(hyper_long_plot[1,])~pheno$treatment, ylab="Beta Value", xlab="Treatment", main=rownames(hyper_long_plot[1,]), col=c("#ff0000", "#b10000", "#620000"))
  stripchart(as.numeric(hyper_long_plot[1,])~pheno$treatment, add=TRUE, vertical=TRUE, method="jitter", pch=16)
dev.off()


#####Match with Epic manifest annotation and save files
epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7)
rownames(epicManifest)<-epicManifest$IlmnID

hypo_short_cgs<-merge(hypo_short_plot, epicManifest, by="row.names")
hypo_long_cgs<-merge(hypo_long_plot, epicManifest, by="row.names")

hyper_short_cgs<-merge(hyper_short_plot, epicManifest, by="row.names")
hyper_long_cgs<-merge(hyper_long_plot, epicManifest, by="row.names")

write.csv(hypo_short_cgs, file="HypoShort_Probes.csv")
write.csv(hypo_long_cgs, file="HypoLong_Probes.csv")
write.csv(hyper_short_cgs, file="HyperShort_Probes.csv")
write.csv(hyper_long_cgs, file="HyperLong_Probes.csv")


write.csv(hypo_short_plot, file="Ordered_hypoShort.csv")
write.csv(hypo_long_plot, file="Ordered_hypoLong.csv")
write.csv(hyper_short_plot, file="Ordered_hyperShort.csv")


####Working out whether spread/outlier in recovery group is same sample

#give treatents in pheno file codes
pheno$treatment_cod <- c(1,1,1,1,2,2,2,2,3,3,3,3)

##test that the plot and labels work (labelling each point with sample ID)
plot(pheno$treatment_cod, as.numeric(hyper_long_plot[1,]), pch=16)
text(pheno$treatment_cod, as.numeric(hyper_long_plot[1,]), labels=pheno$Sample.ID, cex=0.5, pos=2)

##Create labelled plots
pdf("hypo_short_labelled.pdf", height=10, width=10)
par(mfrow=c(3,3))
for(i in 1:97){
  plot(pheno$treatment_cod, as.numeric(hypo_short_plot[i,]), pch=16, xlab="Treatment", ylab="Beta Value", main=rownames(hypo_short_plot[i,]))
  text(pheno$treatment_cod, as.numeric(hypo_short_plot[i,]), labels=pheno$Sample.ID, cex=0.5, pos=2)
}
dev.off()

pdf("hypo_long_labelled.pdf", height=10, width=10)
par(mfrow=c(3,3))
for(i in 1:8){
  plot(pheno$treatment_cod, as.numeric(hypo_long_plot[i,]), pch=16, xlab="Treatment", ylab="Beta Value", main=rownames(hypo_long_plot[i,]))
  text(pheno$treatment_cod, as.numeric(hypo_long_plot[i,]), labels=pheno$Sample.ID, cex=0.5, pos=2)
}
dev.off()

pdf("hyper_short_labelled.pdf", height=10, width=10)
par(mfrow=c(3,3))
for(i in 1:67){
  plot(pheno$treatment_cod, as.numeric(hyper_short_plot[i,]), pch=16, xlab="Treatment", ylab="Beta Value", main=rownames(hyper_short_plot[i,]))
  text(pheno$treatment_cod, as.numeric(hyper_short_plot[i,]), labels=pheno$Sample.ID, cex=0.5, pos=2)
}
dev.off()


pdf("hyper_long_labelled.pdf", height=10, width=10)
plot(pheno$treatment_cod, as.numeric(hyper_long_plot[1,]), pch=16, xlab="Treatment", ylab="Beta Value", main=rownames(hyper_long_plot[1,]))
text(pheno$treatment_cod, as.numeric(hyper_long_plot[1,]), labels=pheno$Sample.ID, cex=0.5, pos=2)
dev.off()



######################
###Pathway Analysis###
######################


###Merge all probes together to run pathway analysis

#load missMethyl package into environment
library(missMethyl)

#merge cpgs
hypo_cpgs<- rbind(hypo_short_plot, hypo_long_plot)
hyper_cpgs <- rbind(hyper_short_plot, hyper_long_plot)
sig_cpgs <- rbind(hypo_cpgs, hyper_cpgs)

#make sure sig_cpgs is just a list of cg IDs
sig_cpgs<-row.names(sig_cpgs)

##load in illumina annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


#read in normalised QC'd betas
betas<-read.csv("/mnt/data1/Jenny/EPIC QC/SV40 EPIC/DasenNormData_badprobesrm.csv", sep=",", header=T)

#set all CpGs to be the cgIDs
all_cpgs<- betas$X



#Doing missmethyl pathway analysis but pulling out gene names as well in pathway
#using function (crystalmeth_fun) provided by Ehsan Pishva

source("/mnt/data1/Ehsan/References/crystalmeth_fun.R")


##Run GO pathway analysis, order and pull out significant pathways
GO<- crystalmeth(sig.cpg = sig_cpgs, 
                 all.cpg = all_cpgs)


ord_paths <- GO[order(GO$P.DE),]
top_paths<-ord_paths[which(ord_paths$P.DE<0.05),]

dim(ord_paths)
dim(top_paths)

#save files
write.csv(GO, file="SV40_Pathway_Crystalmeth.csv")
write.csv(top_paths, file="SV40_SigPathways_Crystalmeth.csv")


##Plotting pathway analysis data as Treemaps using Revigo

#GO pathway IDs and p-values were loaded into Revigo tool online

#Revigo Parameters
#Database = Homo Sapien
#Semantic similarity mesaure = Resnik
#http://revigo.irb.hr/


#Plot for biological pathways (this code is generated by the Revigo tool so you can visualise/adapt it)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0010737","protein kinase A signaling",0.179,2.4170,0.862,0.000,"protein kinase A signaling"),
                     c("GO:0006399","tRNA metabolic process",1.016,2.0238,0.902,0.011,"tRNA metabolism"),
                     c("GO:0008033","tRNA processing",0.664,1.5400,0.898,0.373,"tRNA metabolism"),
                     c("GO:0050873","brown fat cell differentiation",0.231,2.0073,0.873,0.027,"brown fat cell differentiation"),
                     c("GO:0010470","regulation of gastrulation",0.219,1.8183,0.861,0.119,"brown fat cell differentiation"),
                     c("GO:0050672","negative regulation of lymphocyte proliferation",0.369,1.7787,0.779,0.055,"negative regulation of lymphocyte proliferation"),
                     c("GO:0032945","negative regulation of mononuclear cell proliferation",0.369,1.7787,0.781,0.603,"negative regulation of lymphocyte proliferation"),
                     c("GO:0070664","negative regulation of leukocyte proliferation",0.387,1.7446,0.792,0.483,"negative regulation of lymphocyte proliferation"),
                     c("GO:1903055","positive regulation of extracellular matrix organization",0.098,2.3275,0.812,0.055,"positive regulation of extracellular matrix organization"),
                     c("GO:0046324","regulation of glucose import",0.369,1.7814,0.585,0.245,"positive regulation of extracellular matrix organization"),
                     c("GO:0046323","glucose import",0.444,1.7216,0.606,0.512,"positive regulation of extracellular matrix organization"),
                     c("GO:1904659","glucose transmembrane transport",0.121,1.3914,0.595,0.512,"positive regulation of extracellular matrix organization"),
                     c("GO:0015850","organic hydroxy compound transport",1.246,1.4173,0.765,0.202,"positive regulation of extracellular matrix organization"),
                     c("GO:0015749","monosaccharide transport",0.917,1.3325,0.618,0.491,"positive regulation of extracellular matrix organization"),
                     c("GO:1903053","regulation of extracellular matrix organization",0.167,1.7841,0.825,0.434,"positive regulation of extracellular matrix organization"),
                     c("GO:0032370","positive regulation of lipid transport",0.317,1.8646,0.705,0.127,"positive regulation of extracellular matrix organization"),
                     c("GO:0032368","regulation of lipid transport",0.577,1.3425,0.717,0.433,"positive regulation of extracellular matrix organization"),
                     c("GO:0010827","regulation of glucose transport",0.600,1.5913,0.589,0.512,"positive regulation of extracellular matrix organization"),
                     c("GO:0034219","carbohydrate transmembrane transport",0.162,1.3259,0.604,0.491,"positive regulation of extracellular matrix organization"),
                     c("GO:0008645","hexose transport",0.906,1.3394,0.612,0.509,"positive regulation of extracellular matrix organization"),
                     c("GO:0016999","antibiotic metabolic process",0.012,1.3040,0.886,0.056,"antibiotic metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

#Install and load package to plot treemap
install.packages("treemap")
library(treemap) 

#Plot for biological pathways (BP)
pdf( file="revigo_treemap_BP.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "REVIGO GO Biological Pathways Treemap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none",
  fontsize.labels = 20
)
dev.off()


#Do the same for cellular components (CC)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0005811","lipid particle",0.384,1.6720,0.747,0.000,"lipid particle"),
                     c("GO:0015630","microtubule cytoskeleton",5.906,1.3515,0.796,0.239,"lipid particle"),
                     c("GO:0072686","mitotic spindle",0.373,1.4110,0.713,0.168,"lipid particle"),
                     c("GO:0031594","neuromuscular junction",0.303,1.3903,0.849,0.000,"neuromuscular junction"),
                     c("GO:1902554","serine/threonine protein kinase complex",0.481,1.4549,0.860,0.000,"serine/threonine protein kinase complex"),
                     c("GO:0016528","sarcoplasm",0.373,1.3303,0.821,0.030,"sarcoplasm"),
                     c("GO:0031301","integral component of organelle membrane",0.946,1.6168,0.787,0.036,"integral component of organelle membrane"),
                     c("GO:0031300","intrinsic component of organelle membrane",1.001,1.4922,0.788,0.203,"integral component of organelle membrane"),
                     c("GO:0005783","endoplasmic reticulum",9.140,1.5895,0.870,0.055,"endoplasmic reticulum"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

#CC Plot
pdf( file="revigo_treemap_CC.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "REVIGO GO Cellular Component Treemap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none",
  fontsize.labels = 20
)
dev.off()

#Finally do the same for the molecula functions (MF)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0004867","serine-type endopeptidase inhibitor activity",0.560,1.7406,0.738,0.000,"serine-type endopeptidase inhibitor activity"),
                     c("GO:0005179","hormone activity",0.705,1.6130,0.749,0.000,"hormone activity"),
                     c("GO:0005201","extracellular matrix structural constituent",0.451,2.4386,0.728,0.000,"extracellular matrix structural constituent"),
                     c("GO:0016866","intramolecular transferase activity",0.156,2.8605,0.633,0.000,"intramolecular transferase activity"),
                     c("GO:0016853","isomerase activity",0.901,1.4013,0.711,0.115,"intramolecular transferase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

#MF Plot
pdf( file="revigo_treemap_MF.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "REVIGO GO Molecular Function Treemap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none",
  fontsize.labels = 20
)

dev.off()





