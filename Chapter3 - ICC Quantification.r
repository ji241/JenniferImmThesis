##Script to analyse the protein changes occuring throughout iPSC to neuronal differentiation
####Date: 13/11/19
##Author: Jennifer Imm


##R session information
R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-redhat-linux-gnu (64-bit)

##Set working directory
setwd("")


##Load in data for OCT4 and B-tubulin
Oct_Btub<-read.csv("Btub_Oct 4_KW.csv")


##Load in relevant packages
install.packages("PMCMRplus")
library("PMCMRplus")
require(PMCMR)

#Kruskal Wallis test to look for differences in B-tubulin across at different cell stages
#Post hoc Dunn test 
#including bonferroni correction
posthoc.kruskal.dunn.test(Oct_Btub$Beta_Tubulin, Oct_Btub$Cell_Stage, p.adjust.method = "bonferroni")

##output
Pairwise comparisons using Dunn's-test for multiple	
                         comparisons of independent samples 

data:  Oct_Btub$Beta_Tubulin and Oct_Btub$Cell_Stage 

       iPSC    Neuron 
Neuron 7.7e-16 -      
NPC    0.012   3.3e-07

P value adjustment method: bonferroni 

########################################################################################################
#Repeat above test for OCT4 levels

posthoc.kruskal.dunn.test(Oct_Btub$OCT_4, Oct_Btub$Cell_Stage, p.adjust.method = "bonferroni")


##output
Pairwise comparisons using Dunn's-test for multiple	
comparisons of independent samples 

data:  Oct_Btub$OCT_4 and Oct_Btub$Cell_Stage 

        iPSC    Neuron
Neuron  4.7e-08 -     
  NPC   3.2e-16 0.025 

########################################################################################################
##Repeat above test for SOX2 levels
posthoc.kruskal.dunn.test(sox_nest$SOX_2, sox_nest$Cell_Stage, p.adjust.method = "bonferroni")


##output
Pairwise comparisons using Dunn's-test for multiple	
                         comparisons of independent samples 

data:  sox_nest$SOX_2 and sox_nest$Cell_Stage 

       iPSC    NPC   
NPC    0.0005  -     
Neuron 2.7e-11 0.0067

P value adjustment method: bonferroni 


########################################################################################################
##Repeat above test for Nestin levels
posthoc.kruskal.dunn.test(sox_nest$Nestin, sox_nest$Cell_Stage, p.adjust.method = "bonferroni")


##output
Pairwise comparisons using Dunn's-test for multiple	
comparisons of independent samples 

data:  sox_nest$Nestin and sox_nest$Cell_Stage 

            iPSC    NPC    
  NPC    1.3e-06 -      
  Neuron 0.24224 0.00066

P value adjustment method: bonferroni

########################################################################################################
##Repeat above test for NANOG levels
posthoc.kruskal.dunn.test(nan_sse$NANOG, nan_sse$Cell_Stage, p.adjust.method = "bonferroni")

##output
Pairwise comparisons using Dunn's-test for multiple	
                         comparisons of independent samples 

data:  nan_sse$NANOG and nan_sse$Cell_Stage 

    iPSC   
NPC 2.9e-13

########################################################################################################
##Repeat above test for SSEA4 levels
posthoc.kruskal.dunn.test(nan_sse$SSEA4, nan_sse$Cell_Stage, p.adjust.method = "bonferroni")


##output
Pairwise comparisons using Dunn's-test for multiple	
comparisons of independent samples 

data:  nan_sse$SSEA4 and nan_sse$Cell_Stage 

    iPSC   
NPC 2.9e-13


#########################################################################################################
##Create plots showing levels of different proteins over time

pdf("Btub and Oct4 ICC.pdf", height=5, width=10)
par(mfrow=c(1,2))
plot(Oct_Btub$Cell_Stage, Oct_Btub$Beta_Tubulin, ylim=c(0,500), col=c("#ff8989", "#e200e2", "#8989ff"), main="Beta-tubulin levels throughout differentiation", xlab="Cell Stage", ylab="Fluorescence Intensity (a.u)")
stripchart(Oct_Btub$Beta_Tubulin~Oct_Btub$Cell_Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=0.7)
plot(Oct_Btub$Cell_Stage, Oct_Btub$OCT_4, ylim=c(0,500), col=c("#ff8989", "#e200e2", "#8989ff"), main="Oct4 levels throughout differentiation", xlab="Cell Stage", ylab="Fluorescence Intensity (a.u)")
stripchart(Oct_Btub$OCT_4~Oct_Btub$Cell_Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=0.7)
dev.off()


pdf("Nestin and SOX2 ICC.pdf", height=5, width=10)
par(mfrow=c(1,2))
plot(sox_nest$Cell_Stage, sox_nest$Nestin, ylim=c(0,200), col=c("#ff8989", "#e200e2", "#8989ff"), main="Nestin levels throughout differentiation", xlab="Cell Stage", ylab="Fluorescence Intensity (a.u)")
stripchart(sox_nest$Nestin~sox_nest$Cell_Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=0.7)
plot(sox_nest$Cell_Stage, sox_nest$SOX_2, ylim=c(0,400), col=c("#ff8989", "#e200e2", "#8989ff"), main="Sox2 levels throughout differentiation", xlab="Cell Stage", ylab="Fluorescence Intensity (a.u)")
stripchart(sox_nest$SOX_2~sox_nest$Cell_Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=0.7)
dev.off()


pdf("NANOG and SSEA4 ICC.pdf", height=5, width=10)
par(mfrow=c(1,2))
plot(nan_sse$Cell_Stage, nan_sse$NANOG, ylim=c(0,80), col=c("#ff8989", "#e200e2"), main="Nanog levels throughout differentiation", xlab="Cell Stage", ylab="Fluorescence Intensity (a.u)")
stripchart(nan_sse$NANOG~nan_sse$Cell_Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=0.7)
plot(nan_sse$Cell_Stage, nan_sse$SSEA4, ylim=c(0,200), col=c("#ff8989", "#e200e2"), main="SSEA4 levels throughout differentiation", xlab="Cell Stage", ylab="Fluorescence Intensity (a.u)")
stripchart(nan_sse$SSEA4~nan_sse$Cell_Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=0.7)
dev.off()











