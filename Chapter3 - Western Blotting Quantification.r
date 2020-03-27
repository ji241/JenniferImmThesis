##ANOVA analysis of Western Blotting data to look at protein level changes as iPSCs differentiate into neurons
##Date: 26/11/19
##Author: Jenny Imm

#session information
R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-redhat-linux-gnu (64-bit)

##Set working directory
setwd("")

##Load in data for OCT4 and Btub
WB_dat<-read.csv("WesternBlot_Anova.csv")

################################################################

#ANOVA to look for differences in B-tubulin across at different cell stages
#Post hoc Tukey Test


summary(Btub<-aov(WB_dat$Btub~WB_dat$Cell.Stage))
##output
                  Df Sum Sq Mean Sq F value Pr(>F)  
WB_dat$Cell.Stage  2 152397   76199  8.6227 0.0172 *
      Residuals          6  53022    8837                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


################################################################
##ANOVA and post-hoc Tukey for Sox2 levels


summary(SOX2<-aov(WB_dat$SOX2~WB_dat$Cell.Stage))
##output
                  Df Sum Sq Mean Sq F value Pr(>F)
WB_dat$Cell.Stage  2 0.3461  0.1730   0.385  0.696
Residuals          6 2.6984  0.4497

################################################################
##ANOVA and post-hoc Tukey for SSEA4 levels


summary(SSEA4<- aov(WB_dat$SSEA4~WB_dat$Cell.Stage))
##output
                  Df  Sum Sq Mean Sq F value Pr(>F)
WB_dat$Cell.Stage  2 0.03235 0.01618   0.476  0.643
Residuals          6 0.20392 0.03399  

################################################################
##ANOVA and post-hoc Tukey for NESTIN levels


summary(NESTIN <- aov(WB_dat$NESTIN~WB_dat$Cell.Stage))
##output
                  Df   Sum Sq  Mean Sq F value Pr(>F)
WB_dat$Cell.Stage  2 0.003834 0.001917   0.811  0.488
Residuals          6 0.014186 0.002364


################################################################


#OCT4 data was not normally distributed so a Kruskal-Wallis test with post hoc Dunns test was performed

install.packages("dunn.test")
library("dunn.test")

attach(WB_dat)
dunn.test(OCT4, Cell.Stage, kw=TRUE, method="bonferroni")
detach(WB_dat)



##output
data: OCT4 and Cell.Stage
Kruskal-Wallis chi-squared = 7.2, df = 2, p-value = 0.03


Comparison of OCT4 by Cell.Stage                        
(Bonferroni)                                  
Col Mean-|
Row Mean |       iPSC     Neuron
---------+----------------------
  Neuron |   2.683281
         |    0.0109*
         |
     NPC |   1.341640  -1.341640
         |     0.2696     0.2696





################################################################
#Create plots showing levels of different proteins over time


WB_dat$Cell.Stage <- factor(WB_dat$Cell.Stage, levels = c("iPSC", "NPC", "Neuron"))

pdf("WesternBlot_Neuron.pdf")
par(mfrow=c(3,2))
plot(WB_dat$Cell.Stage, WB_dat$Btub, col=c("#ff8989", "#e200e2", "#8989ff"), ylim=c(0,600), main="Beta-tubulin levels throughout differentiation", xlab="Cell Stage", ylab="Relative Expression")
stripchart(WB_dat$Btub~WB_dat$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.2)
plot(WB_dat$Cell.Stage, WB_dat$OCT4, col=c("#ff8989", "#e200e2", "#8989ff"), ylim=c(0,7), main="Oct4 levels throughout differentiation", xlab="Cell Stage", ylab="Relative Expression")
stripchart(WB_dat$OCT4~WB_dat$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.2)
plot(WB_dat$Cell.Stage, WB_dat$SOX2, col=c("#ff8989", "#e200e2", "#8989ff"), ylim=c(0,2.5), main="Sox2 levels throughout differentiation", xlab="Cell Stage", ylab="Relative Expression")
stripchart(WB_dat$SOX2~WB_dat$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.2)
plot(WB_dat$Cell.Stage, WB_dat$SSEA4, col=c("#ff8989", "#e200e2", "#8989ff"), ylim=c(-0.2,0.8), main="SSEA4 levels throughout differentiation", xlab="Cell Stage", ylab="Relative Expression")
stripchart(WB_dat$SSEA4~WB_dat$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.2)
plot(WB_dat$Cell.Stage, WB_dat$NESTIN, col=c("#ff8989", "#e200e2", "#8989ff"), ylim=c(-0.1,0.2), main="Nestin levels throughout differentiation", xlab="Cell Stage", ylab="Relative Expression")
stripchart(WB_dat$NESTIN~WB_dat$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.2)
dev.off()


################################################################
##Embryoid Bodies


#As there was only two samples in each group for the embryoid bodies no statistics were performed, therefore just plots showing the protein changes were created


EB_dat<-read.csv("/mnt/data1/Jenny/Western/WesternBlot_Anova_EB.csv")

EB_dat$Cell.Stage <- factor(EB_dat$Cell.Stage, levels = c("iPSC", "EB"))

pdf("WesternBlot_EB.pdf", width=10, height=5)
par(mfrow=c(1,3))
plot(EB_dat$Cell.Stage, EB_dat$Btub, col=c("#ff8989", "#89ffff"), ylim=c(0,40), main="Beta-tubulin levels throughout differentiation", xlab="Cell Stage", ylab="Relative Expression")
stripchart(EB_dat$Btub~EB_dat$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.2)
plot(EB_dat$Cell.Stage, EB_dat$SMA, col=c("#ff8989", "#89ffff"), ylim=c(0,15), main="SMA levels throughout differentiation", xlab="Cell Stage", ylab="Relative Expression")
stripchart(EB_dat$SMA~EB_dat$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.2)
plot(EB_dat$Cell.Stage, EB_dat$AFP, col=c("#ff8989", "#89ffff"), ylim=c(0,0.3), main="AFP levels throughout differentiation", xlab="Cell Stage", ylab="Relative Expression")
stripchart(EB_dat$AFP~EB_dat$Cell.Stage, method="jitter", vertical=TRUE, add=TRUE, pch=16, cex=1.2)
dev.off()
