# Code to accompany ""
# By:

getwd()

#### Load packages & color schemes####
library(ggplot2)
library(vegan)
library(reshape2)
library(tidyverse)

mycols<-c("#d73027","#fc8d59","#4575b4","#91bfdb")

#### Import strain frequency data ####
# This is a 119 rows each of which represent a set of pooled nodules sampled from Medicago
# There are 69 Columns the first is a concatonated sample name
# Community_Host_TrtRep (e.g C68_A17_L1)
# Ensifer community is one of six communities: C68, C8, C3, MK,KH, HM
# Host is one of two genotypes: A17 or R108
# Treatments: X (high density, No Nitrogen), L (low density, No Nitrogen), N (high density + Nitrogen), LN (low density + Nitrogen)

freqsC68 = read.csv('data/2018-10-19_core_strain_freqs/C68freq.tsv',
                    sep='\t', header=TRUE, as.is=TRUE)

# Store initial frequencies for standardizing and remove from data set
initial<-freqsC68[grep(pattern = "initial",x = freqsC68$pool),]
initial_means = apply(initial[,-1], 2, mean)
freqs<-freqsC68[-grep(pattern = "initial",x = freqsC68$pool),]

# Define relevent columns for downstream analysis from the sample names#
trts<-data.frame(do.call('rbind',strsplit(freqs$pool,"[_]")))
colnames(trts)<-c("Community","Host","trt_rep")
trts$trt<-gsub('([A-Z]+)[0-9]+','\\1',trts$trt_rep)
trts$rep<-as.numeric(gsub('[A-Z]+([0-9]+)','\\1',trts$trt_rep))
trts$Density <- "L"
trts[trts$trt=="X"|trts$trt=="N",]$Density <- "H"
trts$Nitrogen <- "X"
trts[trts$trt=="LN"|trts$trt=="N",]$Nitrogen <- "N"
trts$trt<-factor(trts$trt,levels = c("X","N","L","LN"),labels=c("HiD","HiD+N","LoD","LoD+N"))
trts$Host_trt<-paste(trts$Host,trts$trt,sep="_")
trts$Host_trt<-factor(trts$Host_trt)
trts$Host_trt_rep<-paste(trts$Host_trt,trts$rep,sep="_")

#### RDA, Diversity, and Predicted Benefit (Figure 2) ####
# Calculate relative fitness #
#mine<-as.matrix(freqs[,-1])
#mine<-log2(mine/initial_means)
mine<-log2(mapply("/", freqs[,-1], initial_means))
mine[mine< (-8)]<- (-8)
fitness<-data.frame(trts,mine)  #### Relative fitness values with metadata
freqs<-data.frame(trts,as.matrix(freqs[,-1]))  ### Raw frequency values w/ metadata

##### RDA Analysis of Strain Communities (Table 1 and Figure 2a & S4) #####

# Produce Table 1
# Note that because these tests are based on permutations p-values will vary slightly each time the model is run

## RDA_all_H*D*N model
rdaHDN<-rda(fitness[,c(10:77)]~Host*Density*Nitrogen,fitness, scale=TRUE)
modelHDN<-anova(rdaHDN, step=1000, perm.max=1000, by= "terms")

## RDA_all_H*D model
rdaHD<-rda(fitness[,c(10:77)]~Host*Density,fitness, scale=TRUE)
modelHD<-anova(rdaHD, step=1000, perm.max=1000, by= "terms")

## RDA_A17_D*N model
fitness_A17<- fitness[fitness$Host=="A17",]
rdaADN<-rda(fitness_A17[,c(10:77)]~Density*Nitrogen,fitness_A17, scale=TRUE)
modelADN<-anova(rdaADN, step=1000, perm.max=1000, by= "terms")

## RDA_R108_D*N model
fitness_R108<- fitness[fitness$Host=="R108",]
rdaRDN<-rda(fitness_R108[,c(10:77)]~Density*Nitrogen,fitness_R108, scale=TRUE)
modelRDN<-anova(rdaRDN, step=1000, perm.max=1000, by= "terms")

#Organize model results
HDN<-data.frame(Dataset="All",Term=row.names(modelHDN),
                Df=modelHDN$Df,
                Prop.Var=round(modelHDN$Variance/sum(modelHDN$Variance),3),
                Fstat=round(modelHDN$F,2),
                Pvalue=round(modelHDN$`Pr(>F)`,3),
                Radj=round(as.numeric(RsquareAdj(rdaHDN)[2]),3))
HD<-data.frame(Dataset="All",Term=row.names(modelHD),
                Df=modelHD$Df,
                Prop.Var=round(modelHD$Variance/sum(modelHD$Variance),3),
                Fstat=round(modelHD$F,2),
                Pvalue=round(modelHD$`Pr(>F)`,3),
                Radj=round(as.numeric(RsquareAdj(rdaHD)[2]),3))

ADN<-data.frame(Dataset="A17",Term=row.names(modelADN),
           Df=modelADN$Df,
           Prop.Var=round(modelADN$Variance/sum(modelADN$Variance),3),
           Fstat=round(modelADN$F,2),
           Pvalue=round(modelADN$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(rdaADN)[2]),3))

RDN<-data.frame(Dataset="R108",Term=row.names(modelRDN),
                Df=modelRDN$Df,
                Prop.Var=round(modelRDN$Variance/sum(modelRDN$Variance),3),
                Fstat=round(modelRDN$F,2),
                Pvalue=round(modelRDN$`Pr(>F)`,3),
                Radj=round(as.numeric(RsquareAdj(rdaRDN)[2]),3))

write.table(x = rbind(HDN,HD,ADN,RDN),file = "tables/Table1_RDAModelResultsnew.tsv",sep="\t",row.names=FALSE)

# Figure 2a: Visualize the RDA results

pdf(file="figures/Figure2a_RDAPlot.pdf",width = 4,height=4,useDingbats = FALSE)
scale<-1
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
rdaplot <- ordiplot(rdaHDN, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-2,2),ylim=c(3,-6),scaling=scale, xlab=paste("RDA 1 (",round(summary(rdaHDN)$cont$importance[1,1],2),"% var.)",sep=""), 
                    ylab=paste("RDA 2 (",round(summary(rdaHDN)$cont$importance[1,2],2),"% var.)",sep=""))
points(rdaHDN,"wa", cex=0.8,pch=16,col=rep(mycols,2)[as.numeric(fitness$Host_trt)])
ordiellipse(rdaHDN, fitness$Host_trt, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
ordispider(rdaHDN, fitness$Host_trt, lwd=1,label =TRUE,col=paste(mycols),cex=.5)
dev.off()

pdf(file="figures/FigureS4b_RDAPlotR108.pdf",width = 4,height=4,useDingbats = FALSE)
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
rdaplot <- ordiplot(rdaRDN, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-3,2),ylim=c(-4,4),scaling=scale, xlab=paste("RDA 1 (",round(summary(rdaRDN)$cont$importance[1,1],2),"% var.)",sep=""), 
                    ylab=paste("RDA 2 (",round(summary(rdaRDN)$cont$importance[1,2],2),"% var.)",sep=""))
points(rdaRDN,"wa", cex=0.8,pch=16,col=rep(mycols,2)[as.numeric(fitness_R108$Host_trt)])
ordiellipse(rdaRDN, fitness_R108$Host_trt, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
ordispider(rdaRDN, fitness_R108$Host_trt, lwd=1,label =TRUE,col=paste(mycols),cex=.5)
dev.off()

pdf(file="figures/FigureS4b_RDAPlotA17.pdf",width = 4,height=4,useDingbats = FALSE)
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
rdaplot <- ordiplot(rdaADN, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-2,2),ylim=c(4,-4),scaling=scale, xlab=paste("RDA 1 (",round(summary(rdaADN)$cont$importance[1,1],2),"% var.)",sep=""), 
                    ylab=paste("RDA 2 (",round(summary(rdaADN)$cont$importance[1,2],2),"% var.)",sep=""))
points(rdaADN,"wa", cex=0.8,pch=16,col=rep(mycols,2)[as.numeric(fitness_A17$Host_trt)])
ordiellipse(rdaADN, fitness_A17$Host_trt, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
ordispider(rdaADN, fitness_A17$Host_trt, lwd=1,label =TRUE,col=paste(mycols),cex=.5)
dev.off()

#### Diversity analysis (Figure 2b & Table S1) ######

# Calculate diversity metrics
mydiver<-data.frame(freqs[,1:9],renyi(as.matrix(freqs[,-1:-9])))

# Run an ANOVA
write.table(x = anova(lm(X1~Host*Density*Nitrogen,data=mydiver)),file="tables/TableS1_Diversity_BothHosts.tsv",sep="\t",row.names = TRUE,col.names = NA)
write.table(x = anova(lm(X1~Density*Nitrogen,data=mydiver[mydiver$Host=="A17",])),file="tables/TableS1_Diversity_A17Host.tsv",sep="\t",row.names = TRUE,col.names = NA)
write.table(x = anova(lm(X1~Density*Nitrogen,data=mydiver[mydiver$Host=="R108",])),file="tables/TableS1_Diversity_R108Host.tsv",sep="\t",row.names = TRUE,col.names = NA)

# Make Diversity Figure
pdf(file = "figures/Figure2c_DiversityBoxplot.pdf",height=4,width=4)
boxplot(mydiver$X1~mydiver$Host_trt,col=rep(mycols,2),ylab="Shannon's Diversity",xlab="",lty=1,axes=FALSE,outline=FALSE)# Shannons diversity
axis(side = 1,at = c(1:8),labels = paste(levels(mydiver$Host_trt)),las=2,cex.axis=.75)
axis(side = 2 )
points(mydiver$X1~jitter(as.numeric(mydiver$Host_trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
abline(v = 4.5,col="black",lty=2)
dev.off()


#### Predicted Plant Benefit analysis (Figure 2c and Table S3) ####

# Load in data (this is the same data used in original Burghardt 2018 PNAS paper)
benefit<-read.csv(file = "data/SingleStrain_phenotype_summary.tsv",sep = '\t')
benefit$strain<-paste("X",benefit$strain,sep="")
mystrains<-colnames(fitness)[-1:-9]

#Subset down to strains used in this experiment
A17benefit<-as_tibble(benefit) %>% filter(plant_genotype=="A17",strain %in% !!mystrains) %>% select(c(strain,weight))
R108benefit<-as_tibble(benefit) %>% filter(plant_genotype=="R108",strain %in% !!mystrains) %>% select(c(strain,weight))

freq<-data.frame(t(freqs[,c(-1:-9)]))
colnames(freq)<- paste(trts$Host_trt_rep)
freq<-cbind(freq,initial=initial_means)
freq$strain<-rownames(freq)
R108freq<-merge(x=R108benefit,y=freq,by="strain")
A17freq<-merge(x=A17benefit,y=freq,by="strain")

Predsize<-data.frame(R108=(colSums(R108freq[,-1:-2]*R108freq$weight,na.rm = TRUE)/colSums(R108freq[,-1:-2],na.rm = TRUE)-min(R108freq$weight,na.rm=TRUE))/(max(R108freq$weight,na.rm=TRUE)-min(R108freq$weight,na.rm=TRUE)),
                     A17=(colSums(A17freq[,-1:-2]*A17freq$weight,na.rm = TRUE)/colSums(A17freq[,-1:-2],na.rm = TRUE)-min(A17freq$weight,na.rm=TRUE))/(max(A17freq$weight,na.rm=TRUE)-min(A17freq$weight,na.rm=TRUE)))

initialA17<-Predsize["initial","A17"]
initialR108<-Predsize["initial","R108"]

Predsize<-Predsize[-49,]

Predsize$Host_trt<-factor(freqs$Host_trt)
Predsize$Host<-factor(freqs$Host)
Predsize$Density<-factor(freqs$Density)
Predsize$Nitrogen<-factor(freqs$Nitrogen)

# Focus in on A17 benefits for A17 hosts and R108 benefits for R108
Predsize$Benefit <- 0
Predsize[Predsize$Host=="A17",]$Benefit <- Predsize[Predsize$Host=="A17",]$A17
Predsize[Predsize$Host=="R108",]$Benefit <- Predsize[Predsize$Host=="R108",]$R108

# Create a Boxplot for Figure 2c
pdf("figures/Figure2c_PredictedPlantBenefit.pdf",width = 4,height=4,useDingbats = FALSE)
boxplot(Predsize$Benefit~Predsize$Host_trt,col=rep(mycols,2),xlab="Treatment",ylab="Predicted Benefit",ylim=c(.25,.75),outline=FALSE,axes=FALSE,lty=1)
axis(side = 1,at = c(1:8),labels = paste(levels(Predsize$Host_trt)),las=2,cex.axis=.5)
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
segments(x0=1,x1=4,y0=initialA17,col="grey40",lty=1)
segments(x0=5,x1=8,y0=initialR108,col="grey40",lty=1)
points(Predsize$Benefit~jitter(as.numeric(Predsize$Host_trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
dev.off()

# Table S2: ANOVA's models for Host Benefit
modelHDN<-anova(lm(Benefit~Host*Density*Nitrogen,data=Predsize))
modelADN<-anova(lm(Benefit~Density*Nitrogen,data=Predsize[Predsize$Host=="A17",]))
modelRDN<-anova(lm(Benefit~Density*Nitrogen,data=Predsize[Predsize$Host=="R108",]))

HDN<-data.frame(Host="All",Term=row.names(modelHDN),
                Df=modelHDN$Df,
                SumSq=round(modelHDN$`Sum Sq`,3),
                Fstat=round(modelHDN$`F value`,2),
                Pvalue=round(modelHDN$`Pr(>F)`,3))

ADN<-data.frame(Host="All",Term=row.names(modelADN),
                Df=modelADN$Df,
                SumSq=round(modelADN$`Sum Sq`,3),
                Fstat=round(modelADN$`F value`,2),
                Pvalue=round(modelADN$`Pr(>F)`,3))

RDN<-data.frame(Host="All",Term=row.names(modelRDN),
                Df=modelRDN$Df,
                SumSq=round(modelRDN$`Sum Sq`,3),
                Fstat=round(modelRDN$`F value`,2),
                Pvalue=round(modelRDN$`Pr(>F)`,3))

write.table(x = rbind(HDN,ADN,RDN),file = "tables/TableS2_BenefitModels.tsv",sep="\t",row.names=FALSE)

####  Load in Data for Plant Pheotypic Trait Analysis ####
# Graphs for Nod Num, weight, Veg, Root
# ANOVA with complexity as continuous variable H*C

#Load in Phenotypic Data and organize
harvest<-read.csv(file="data/R108_SR3_HarvestData_Final.txt",sep="\t",header = TRUE)
harvest$Plants_harvested<-harvest$Plant_num
all<-harvest[,c(1:10,17,11:16)]
A17<-read.csv(file="data/A17_SR3_HarvestData_Final.txt",sep="\t",header = TRUE)
all<-rbind(all,A17)
all$Treatment<- factor(x = all$Treatment,levels = c("X","N","L","LN"),labels=c("HiD","HiD+N","LoD","LoD+N"))
all$Host<-all$Plant_genotype
all$Density<- "H"
all[all$Treatment=="LoD"|all$Treatment=="LoD+N",]$Density<- "L"
all$Nitrogen<- "X"
all[all$Treatment=="HiD+N"|all$Treatment=="LoD+N",]$Nitrogen<- "N"
all$Complexity <- 68
all[all$Community=="C8",]$Complexity<- 8
all[all$Community=="C3",]$Complexity<- 3
all[all$Community=="HM"|all$Community=="HK"|all$Community=="MK",]$Complexity<- 3

#Calculate variables of interest
all$weight_per_Nod<-all$Nod_weight_grams/all$Nod_num
all$Veg_per_plant<- all$Veg_dryweight_g/all$Plant_num
all$Root_per_plant<- all$Root_dryweight_g/all$Plant_num
all$Nodnum_per_plant<- all$Nod_num/all$Plants_harvested
all$Nodweight_per_plant<- (all$weight_per_Nod*all$Nodnum_per_plant)
all$RootShoot<- (all$Root_per_plant/all$Veg_per_plant)

# Grab relevent columns
all<-all[,c("Host","Community","Treatment","Density","Nitrogen","Complexity","Harvest_date","Nod_num","Nodnum_per_plant","weight_per_Nod","Nodweight_per_plant","Veg_per_plant","Root_per_plant","RootShoot")]

# Create two datasets: Nitrogen*Density & Community Complexity
Pheno_HND<-all[all$Community=="C68",]
Pheno_Complexity<-rbind(all[all$Community=="C68"&all$Treatment=="HiD",],all[!all$Community=="C68",])

#### Nitrogen*Density Traits (Figures 3, S2a, S7, & Table S3) ####
# Run Factorial ANOVA models for all six focal Nodule and Host traits
modNN<-anova(lm(Nodnum_per_plant~Host*Density*Nitrogen, data=Pheno_HND)) # Host*Density
modNW<-anova(lm(weight_per_Nod*1000~Host*Density*Nitrogen, data=Pheno_HND)) # Host*Density*Nitrogen
modNWP<-anova(lm(Nodweight_per_plant*1000~Host*Density*Nitrogen, data=Pheno_HND)) # NONE
modV<-anova(lm(Veg_per_plant*1000~Host*Density*Nitrogen, data=Pheno_HND)) #### Only Host
modR<-anova(lm(Root_per_plant*1000~Host*Density*Nitrogen, data=Pheno_HND)) ### Only Host
modRS<-anova(lm(RootShoot~Host*Density*Nitrogen, data=Pheno_HND)) ### Host*Nitrogen and marginal Host*Density

# Run submodels divided by host for those traits that had a significant interaction
modNN.A<-anova(lm(Nodnum_per_plant~Density*Nitrogen, data=Pheno_HND[Pheno_HND$Host=="A17",]))
modNN.R<-anova(lm(Nodnum_per_plant~Density*Nitrogen, data=Pheno_HND[Pheno_HND$Host=="R108",]))
modNW.A<-anova(lm(weight_per_Nod*1000~Density*Nitrogen, data=Pheno_HND[Pheno_HND$Host=="A17",]))
modNW.R<-anova(lm(weight_per_Nod*1000~Density*Nitrogen, data=Pheno_HND[Pheno_HND$Host=="R108",]))
modRS.A<-anova(lm(RootShoot~Density*Nitrogen, data=Pheno_HND[Pheno_HND$Host=="A17",]))
modRS.R<-anova(lm(RootShoot~Density*Nitrogen, data=Pheno_HND[Pheno_HND$Host=="R108",]))

# Organize the data for saving
full<-data.frame(Terms=row.names(modNN),DF=modNN$Df,
           NN.SS=round(modNN$`Sum Sq`,1),NN.P=round(modNN$`Pr(>F)`,3),
           NW.SS=round(modNW$`Sum Sq`,1),NW.P=round(modNW$`Pr(>F)`,3),
           NWP.SS=round(modNWP$`Sum Sq`,1),NWP.P=round(modNWP$`Pr(>F)`,3),
           V.SS=round(modV$`Sum Sq`,1),V.P=round(modV$`Pr(>F)`,3),
           R.SS=round(modR$`Sum Sq`,1),R.P=round(modR$`Pr(>F)`,3),
           RS.SS=round(modRS$`Sum Sq`,1),RS.P=round(modRS$`Pr(>F)`,3))

sub<-data.frame(Terms=row.names(modNN.A),DF=modNN.A$Df,
                 NN.A.SS=round(modNN.A$`Sum Sq`,1),NN.A.P=round(modNN.A$`Pr(>F)`,3),
                 NN.R.SS=round(modNN.R$`Sum Sq`,1),NN.R.P=round(modNN.R$`Pr(>F)`,3),
                 NW.A.SS=round(modNW.A$`Sum Sq`,1),NW.A.P=round(modNW.A$`Pr(>F)`,3),
                 NW.R.SS=round(modNW.R$`Sum Sq`,1),NW.R.P=round(modNW.R$`Pr(>F)`,3),
                 RS.A.SS=round(modRS.A$`Sum Sq`,1),RS.A.P=round(modRS.A$`Pr(>F)`,3),
                 RS.R.SS=round(modRS.R$`Sum Sq`,1),RS.R.P=round(modRS.R$`Pr(>F)`,3))
               
write.table(x = cbind(full),file = "tables/TableS3_PlantPhenotypes_FullModels.tsv",sep="\t",row.names=FALSE)
write.table(x = cbind(sub),file = "tables/TableS3_PlantPhenotypes_HostSubsetModels.tsv",sep="\t",row.names=FALSE)

# Graph the Nodule traits
Pheno_HND$Trt<-paste(Pheno_HND$Host,Pheno_HND$Treatment,sep="_")
Pheno_HND$Trt<-factor(Pheno_HND$Trt,levels=c("A17_HiD","A17_HiD+N","A17_LoD","A17_LoD+N","R108_HiD","R108_HiD+N","R108_LoD","R108_LoD+N"))

pdf(file = "figures/Figure3_NoduleBoxplot_N*D.pdf",height=5,width=5,useDingbats = FALSE)
boxplot(Nodnum_per_plant~Trt,col=mycols[c(1:4,1:4)],xlab="Treatment",ylab="Nodule number per plant",outline=FALSE,axes=FALSE,lty=1, data=Pheno_HND)
axis(side = 1,at = c(1:8),labels = paste(levels(Pheno_HND$Trt)))
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
points(Pheno_HND$Nodnum_per_plant~jitter(as.numeric(Pheno_HND$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)

boxplot(weight_per_Nod*1000~Trt,col=mycols[c(1:4,1:4)],xlab="Treatment",ylab="Mean weight per nodule (mg)",outline=FALSE,axes=FALSE,lty=1, data=Pheno_HND)
axis(side = 1,at = c(1:8),labels = paste(levels(Pheno_HND$Trt)))
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
points(Pheno_HND$weight_per_Nod*1000~jitter(as.numeric(Pheno_HND$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)

boxplot(Nodweight_per_plant*1000~Trt,col=mycols[c(1:4,1:4)],xlab="Treatment",ylab="Nodule weight per plant (mg)",outline=FALSE,axes=FALSE,lty=1, data=Pheno_HND)
axis(side = 1,at = c(1:8),labels = paste(levels(Pheno_HND$Trt)))
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
points(Pheno_HND$Nodweight_per_plant*1000~jitter(as.numeric(Pheno_HND$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)

dev.off()

pdf(file = "figures/FigureS2a_PlantBoxplot_N*D.pdf",height=5,width=6.5, useDingbats = FALSE)
boxplot(Nod_num~Trt,col=mycols[c(1:4,1:4)],xlab="Treatment",ylab="Nodule number sampled",outline=FALSE,axes=FALSE,lty=1, data=Pheno_HND)
axis(side = 1,at = c(1:8),labels = c(paste(levels(Pheno_HND$Treatment)),paste(levels(Pheno_HND$Treatment))))
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
points(Pheno_HND$Nod_num~jitter(as.numeric(Pheno_HND$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
text(x=c(2.5,6.5),y=c(550,550),labels=c("A17","R108"))
dev.off()

pdf(file = "figures/FigureS7_PlantBoxplot_N*D.pdf",height=5,width=6.5, useDingbats = FALSE)
boxplot(Veg_per_plant~Trt,col=mycols[c(1:4,1:4)],xlab="Treatment",ylab="Veg biomass per plant (g)",outline=FALSE,axes=FALSE,lty=1, data=Pheno_HND)
axis(side = 1,at = c(1:8),labels = c(paste(levels(Pheno_HND$Treatment)),paste(levels(Pheno_HND$Treatment))))
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
points(Pheno_HND$Veg_per_plant~jitter(as.numeric(Pheno_HND$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
text(x=c(2.5,6.5),y=c(.7,.7),labels=c("A17","R108"))

boxplot(Root_per_plant~Trt,col=mycols[c(1:4,1:4)],xlab="Treatment",ylab="Root biomass per plant (g)",outline=FALSE,axes=FALSE,lty=1, data=Pheno_HND)
axis(side = 1,at = c(1:8),labels = c(paste(levels(Pheno_HND$Treatment)),paste(levels(Pheno_HND$Treatment))))
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
points(Pheno_HND$Root_per_plant~jitter(as.numeric(Pheno_HND$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
text(x=c(2.5,6.5),y=c(.7,.7),labels=c("A17","R108"))

boxplot(RootShoot~Trt,col=mycols[c(1:4,1:4)],xlab="Treatment",ylab="Root:Shoot ratio",outline=FALSE,axes=FALSE,lty=1, data=Pheno_HND)
axis(side = 1,at = c(1:8),labels = c(paste(levels(Pheno_HND$Treatment)),paste(levels(Pheno_HND$Treatment))))
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
abline(h = 1,col="grey",lty=1)
points(Pheno_HND$RootShoot~jitter(as.numeric(Pheno_HND$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
text(x=c(2.5,6.5),y=c(2,2),labels=c("A17","R108"))

boxplot(RootShoot~Trt,col=mycols[c(1:4,1:4)],xlab="Treatment",ylab="Root:Shoot ratio",outline=FALSE,axes=FALSE,lty=1, data=Pheno_HND)
axis(side = 1,at = c(1:8),labels = c(paste(levels(Pheno_HND$Treatment)),paste(levels(Pheno_HND$Treatment))))
axis(side = 2 )
abline(v = 4.5,col="black",lty=2)
abline(h = 1,col="grey",lty=1)
points(Pheno_HND$RootShoot~jitter(as.numeric(Pheno_HND$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
text(x=c(2.5,6.5),y=c(2,2),labels=c("A17","R108"))

dev.off()

#### Community Complexity Traits (Figures S2b, S8 & Table S4) #####
Pheno_Complexity$Host<- factor(Pheno_Complexity$Host,levels=c("R108","A17"))

# Run anova Models with Complexity as a quantitative predictor
modNN<-anova(lm(Nodnum_per_plant~Host*Complexity, data=Pheno_Complexity)) # Host
modNW<-anova(lm(weight_per_Nod*1000~Host*Complexity, data=Pheno_Complexity)) # Host*Complexity
modV<-anova(lm(Veg_per_plant~Host*Complexity, data=Pheno_Complexity)) #### Only Host
modR<-anova(lm(Root_per_plant~Host*Complexity, data=Pheno_Complexity)) ### Host, Complexity marginal
modRS<-anova(lm(RootShoot~Host*Complexity, data=Pheno_Complexity)) ### none
modNWP<-anova(lm(Nodweight_per_plant~Host*Complexity, data=Pheno_Complexity))### none

# Host submodels for Nodweight and root per plant showed inklings of interactions
modNW.A<-anova(lm(weight_per_Nod*1000~Complexity, data=Pheno_Complexity[Pheno_Complexity$Host=="A17",]))
modNW.R<-anova(lm(weight_per_Nod*1000~Complexity, data=Pheno_Complexity[Pheno_Complexity$Host=="R108",])) 
modR.A<-anova(lm(Root_per_plant~Complexity, data=Pheno_Complexity[Pheno_Complexity$Host=="A17",]))
modR.R<-anova(lm(Root_per_plant~Complexity, data=Pheno_Complexity[Pheno_Complexity$Host=="R108",])) 

# Compile the results for writing to a table
full<-data.frame(Terms=row.names(modNN),DF=modNN$Df,
                 NN.SS=round(modNN$`Sum Sq`,3),NN.P=round(modNN$`Pr(>F)`,3),
                 NW.SS=round(modNW$`Sum Sq`,3),NW.P=round(modNW$`Pr(>F)`,3),
                 NWP.SS=round(modNWP$`Sum Sq`,3),NWP.P=round(modNWP$`Pr(>F)`,3),
                 V.SS=round(modV$`Sum Sq`,3),V.P=round(modV$`Pr(>F)`,3),
                 R.SS=round(modR$`Sum Sq`,3),R.P=round(modR$`Pr(>F)`,3),
                 RS.SS=round(modRS$`Sum Sq`,3),RS.P=round(modRS$`Pr(>F)`,3))

write.table(x = cbind(full),file = "tables/TableS4_CommunityComplexityModels.tsv",sep="\t",row.names=FALSE)

#Figure Graphs
mycols.x<-rev(c("#ffffcc","#c2e699","#78c679","#238443"))
Pheno_Complexity$Community<-factor(Pheno_Complexity$Community,levels=c("C68","C8","C3","HK","HM","MK"),labels= c("C68","C8","C3","H:K","M:H","M:K"))
Pheno_Complexity$Trt<-factor(paste(Pheno_Complexity$Host,Pheno_Complexity$Community,sep="_"),levels=c("A17_C68",  "A17_C8","A17_C3", "A17_H:K" ,  "A17_M:H" ,  "A17_M:K" ,  "R108_C68",  "R108_C8",
                                                                                                      "R108_C3",  "R108_H:K" , "R108_M:H",  "R108_M:K"))
pdf(file = "figures/FigureS8_ComplexityPhenotypes.pdf",height=5,width=7)
  boxplot(Nodnum_per_plant~Trt,col=mycols.x[c(1:4,4,4,1:4,4,4)],xlab="Community",ylab="Nodule per plant",outline=FALSE,axes=FALSE,lty=1, data=Pheno_Complexity)
    axis(side = 1,at = c(1:12),labels = c(paste(levels(Pheno_Complexity$Community)),paste(levels(Pheno_Complexity$Community))))
    axis(side = 2 )
    abline(v = 6.5,col="black",lty=2)
    points(Pheno_Complexity$Nodnum_per_plant~jitter(as.numeric(Pheno_Complexity$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
    text(x=c(3.5,9.5),y=c(135,135),labels=c("A17","R108"))

  boxplot(weight_per_Nod*1000~Trt,col=mycols.x[c(1:4,4,4,1:4,4,4)],xlab="Community",ylab="Mean nodule wieght (mg)",outline=FALSE,axes=FALSE,lty=1, data=Pheno_Complexity)
    axis(side = 1,at = c(1:12),labels = c(paste(levels(Pheno_Complexity$Community)),paste(levels(Pheno_Complexity$Community))))
    axis(side = 2 )
    abline(v = 6.5,col="black",lty=2)
    points(Pheno_Complexity$weight_per_Nod*1000~jitter(as.numeric(Pheno_Complexity$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
    text(x=c(3.5,9.5),y=c(5.8,5.8),labels=c("A17","R108"))

  boxplot(Nodweight_per_plant*1000~Trt,col=mycols.x[c(1:4,4,4,1:4,4,4)],xlab="Community",ylab="Nodule wieght per plant (mg)",outline=FALSE,axes=FALSE,lty=1, data=Pheno_Complexity)
    axis(side = 1,at = c(1:12),labels = c(paste(levels(Pheno_Complexity$Community)),paste(levels(Pheno_Complexity$Community))))
    axis(side = 2 )
    abline(v = 6.5,col="black",lty=2)
    points(Pheno_Complexity$Nodweight_per_plant*1000~jitter(as.numeric(Pheno_Complexity$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
    text(x=c(3.5,9.5),y=c(105,105),labels=c("A17","R108"))

  boxplot(Veg_per_plant*1000~Trt,col=mycols.x[c(1:4,4,4,1:4,4,4)],xlab="Community",ylab="Vegetative biomass per plant (mg)",outline=FALSE,axes=FALSE,lty=1, data=Pheno_Complexity)
    axis(side = 1,at = c(1:12),labels = c(paste(levels(Pheno_Complexity$Community)),paste(levels(Pheno_Complexity$Community))))
    axis(side = 2 )
    abline(v = 6.5,col="black",lty=2)
    points(Pheno_Complexity$Veg_per_plant*1000~jitter(as.numeric(Pheno_Complexity$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
    text(x=c(3.5,9.5),y=c(550,550),labels=c("A17","R108"))

  boxplot(Root_per_plant*1000~Trt,col=mycols.x[c(1:4,4,4,1:4,4,4)],xlab="Community",ylab="Root biomass per plant (mg)",outline=FALSE,axes=FALSE,lty=1, data=Pheno_Complexity)
    axis(side = 1,at = c(1:12),labels = c(paste(levels(Pheno_Complexity$Community)),paste(levels(Pheno_Complexity$Community))))
    axis(side = 2 )
    abline(v = 6.5,col="black",lty=2)
    points(Pheno_Complexity$Root_per_plant*1000~jitter(as.numeric(Pheno_Complexity$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
    text(x=c(3.5,9.5),y=c(650,650),labels=c("A17","R108"))

  boxplot(RootShoot~Trt,col=mycols.x[c(1:4,4,4,1:4,4,4)],xlab="Community",ylab="Root:Shoot ratio",outline=FALSE,axes=FALSE,lty=1, data=Pheno_Complexity)
    axis(side = 1,at = c(1:12),labels = c(paste(levels(Pheno_Complexity$Community)),paste(levels(Pheno_Complexity$Community))))
    axis(side = 2 )
    abline(v = 6.5,col="black",lty=2)
    points(Pheno_Complexity$RootShoot~jitter(as.numeric(Pheno_Complexity$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
    text(x=c(3.5,9.5),y=c(1.9,1.9),labels=c("A17","R108"))
    abline(h=1,col="grey")
dev.off()

pdf(file = "figures/FigureS2b_ComplexityNodulesSampled.pdf",height=5,width=7)
boxplot(Nod_num~Trt,col=mycols.x[c(1:4,4,4,1:4,4,4)],xlab="Community",ylab="Nodule pool",outline=FALSE,axes=FALSE,lty=1, data=Pheno_Complexity)
axis(side = 1,at = c(1:12),labels = c(paste(levels(Pheno_Complexity$Community)),paste(levels(Pheno_Complexity$Community))))
axis(side = 2 )
abline(v = 6.5,col="black",lty=2)
points(Pheno_Complexity$Nod_num~jitter(as.numeric(Pheno_Complexity$Trt)),col=c("black"), pch=19 , ylim=c(.25,.75),cex=.75)
text(x=c(1.5,7.5),y=c(650,650),labels=c("A17","R108"))
abline(h=1,col="grey")
dev.off()

#### Community Complexity Strain Frequencies (Figure 4, Figure S3, S6)  ####
library(colorspace)
library(scales)
library(Hmisc)

#Load in the Data
freqsHM = read.csv('data/2018-10-19_core_strain_freqs/HMfreq.tsv',sep='\t', header=TRUE, as.is=TRUE)
freqsKH = read.csv('data/2018-10-19_core_strain_freqs/KHfreq.tsv',sep='\t', header=TRUE, as.is=TRUE)
freqsMK = read.csv('data/2018-10-19_core_strain_freqs/MKfreq.tsv',sep='\t', header=TRUE, as.is=TRUE)

freqsC3 = read.csv('data/2018-10-19_core_strain_freqs/C3freq.tsv',sep='\t', header=TRUE, as.is=TRUE)
freqsC3<-freqsC3[grep(pattern = "C3",x = freqsC3$pool),]

freqsC8 = read.csv('data/2018-10-19_core_strain_freqs/C8freq.tsv',sep='\t', header=TRUE, as.is=TRUE)
freqsC8sub3<-freqsC8[,colnames(freqsC8) %in% colnames(freqsC3)]
freqsC68sub<-rbind(freqsC68[grep(pattern = "A17_X",x = freqsC68$pool),], freqsC68[grep(pattern = "R108_X",x = freqsC68$pool),],freqsC68[grep(pattern = "initial",x = freqsC68$pool),])

freqsC68sub8<-freqsC68sub[,colnames(freqsC68sub) %in% colnames(freqsC8)]
freqsC68sub3<-freqsC68sub[,colnames(freqsC68sub) %in% colnames(freqsC3)]

freqsC68subHM<-freqsC68sub[,colnames(freqsC68sub) %in% colnames(freqsHM)]
freqsC8subHM<-freqsC8[,colnames(freqsC8) %in% colnames(freqsHM)]
freqsC3subHM<-freqsC3[,colnames(freqsC3) %in% colnames(freqsHM)]

freqsC68subKH<-freqsC68sub[,colnames(freqsC68sub) %in% colnames(freqsKH)]
freqsC8subKH<-freqsC8[,colnames(freqsC8) %in% colnames(freqsKH)]
freqsC3subKH<-freqsC3[,colnames(freqsC3) %in% colnames(freqsKH)]

freqsC68subMK<-freqsC68sub[,colnames(freqsC68sub) %in% colnames(freqsMK)]
freqsC8subMK<-freqsC8[,colnames(freqsC8) %in% colnames(freqsMK)]
freqsC3subMK<-freqsC3[,colnames(freqsC3) %in% colnames(freqsMK)]

FrequencyAll<-function(freqs,title,mycols) {
  mycols<-data.frame(Strain=c("HM006.1","M162","KH46C","KH35C","M270","T073","USDA1021", "USDA1157"),
                     colors=c("darkblue","darkgoldenrod","olivedrab","darkred","black","gold","darkorange","lightblue"))
  initial_pools = data.frame(Host="initial",freqs[grep(pattern = "initial",x = freqs$pool),])
  A17_pools=data.frame(Host="A17",freqs[grep(pattern = "A17",x = freqs$pool),])
  R108_pools=data.frame(Host="R108",freqs[grep(pattern = "R108",x = freqs$pool),])
  Prop=rbind(initial_pools,A17_pools,R108_pools)
  
  mdfr <- melt(Prop, id.vars = c("Host","pool"))
  colnames(mdfr)<-c("Environment","pool","Strain","value")
  mycol<-mycols[mycols$Strain %in% mdfr$Strain,]
  mdfr$Strain<-factor(mdfr$Strain, levels = paste(mycol$Strain))
  #ggplot(mdfr, aes(pool, value, fill = Strain)) +
    ggplot(mdfr, aes(Environment, value, fill = Strain)) +
    geom_bar(position = "fill", stat = "identity")+
    scale_y_continuous(labels = percent)+
    scale_fill_manual(values = paste(mycol$colors))+
    ggtitle(paste(title))+
    theme_minimal()
}

pdf(file="figures/FigureS6_RelativeAbundanceAcrossComplexity.pdf",height = 4,width=3,fonts = "Helvetica",useDingbats = FALSE)
FrequencyAll(freqsHM, "HM Pairwise")
FrequencyAll(freqsKH,"KH Pairwise")
FrequencyAll(freqsMK, "MK Pairwise")
FrequencyAll(freqsC3, "C3 Community")
FrequencyAll(freqsC8sub3, "C3 Strains in C8")
FrequencyAll(freqsC68sub3, "C3 Strains in C68")
FrequencyAll(freqsC8, "C8 Community")
FrequencyAll(freqsC68sub8, "C8 Strains in C68")
dev.off()

##### Community Composition for N*D Treatments (Figure S3)  #####

#Create the DataFrame

Prop<-rbind(freqsC68[,c(-1:-7,-9)],cbind(Host_trt="initial",initial[,-1]))
mdfr <- melt(Prop, id.vars = c("Host_trt"))
colnames(mdfr)<-c("Environment","Strain","value")
mdfr$value <-as.numeric(mdfr$value)
mdfr$Environment<- factor(mdfr$Environment,levels=c("A17_HiD","A17_HiD+N","A17_LoD","A17_LoD+N","initial","R108_HiD","R108_HiD+N","R108_LoD","R108_LoD+N"))

# Get the colors and orders right
set.seed(seed = 2)
mycols<-data.frame(Strain=c("HM006.1","M162","KH46C","KH35C","M270","T073","USDA1021", "USDA1157"),
                   colors=c("darkblue","darkgoldenrod","olivedrab","darkred","black","gold","darkorange","lightblue"))
mdfr$Strain<-factor(mdfr$Strain, levels = c(as.character(mycols$Strain),levels(mdfr$Strain)[!levels(mdfr$Strain) %in% mycols$Strain]))
mine<-data.frame(Strain=levels(mdfr$Strain)[!levels(mdfr$Strain) %in% mycols$Strain],colors=sample(gray.colors(n = 60,start = .05,end=1),size = 60))
mycols<-rbind(mycols,mine)

A17_pools$Trt <- NULL
myorder<-colnames(A17_pools[-c(1:2)])[order(colMeans(A17_pools[,-c(1:2)]),decreasing=TRUE)]
mycols<-mycols[match(myorder,mycols$Strain),]
mdfr$Strain<-factor(mdfr$Strain, levels = c(myorder))

#Create Figure S3
pdf(file="figures/FigureS3_NoduleComposition_N*D.pdf",height = 6,width=7,fonts = "Helvetica",useDingbats = FALSE)
  ggplot(mdfr, aes(Environment, value, fill = Strain)) +
    geom_bar(position = "fill", stat = "identity")+
    scale_y_continuous(labels = percent)+
    scale_fill_manual(values = paste(mycols$colors))+
    theme_minimal()
dev.off()

#### Community Complexity Rank Order (Figure 4) ####

C8<-as.tibble(rbind(freqsC8,freqsC68sub8)) %>% separate(pool,into=c("Community","Host",NA),sep = "_") %>% 
  filter(Host!="initial") %>% unite(Combo,c(Host,Community)) %>% group_by(Combo) %>% summarize_if(is.numeric, mean) %>% 
  pivot_longer(USDA1157:KH46C,'Strain',values_to = 'freq') %>% pivot_wider(names_from=Combo,values_from=freq) %>% mutate_if(is.double,function(x) dense_rank(-x))

C3<-as.tibble(rbind(freqsC3,freqsC68sub3)) %>% separate(pool,into=c("Community","Host",NA),sep = "_") %>% 
  filter(Host!="initial") %>% unite(Combo,c(Host,Community)) %>% group_by(Combo) %>% summarize_if(is.numeric, mean) %>% 
  pivot_longer(HM006.1:KH46C,'Strain',values_to = 'freq') %>% pivot_wider(names_from=Combo,values_from=freq) %>% mutate_if(is.double,function(x) dense_rank(-x))

mycols<-data.frame(Strain=c("HM006.1","M162","KH46C","KH35C","M270","T073","USDA1021", "USDA1157"),
                   colors=c("darkblue","darkgoldenrod","olivedrab","darkred","black","gold","darkorange","lightblue"))

C8<-merge(C8,mycols,by="Strain")
C3<-merge(C3,mycols,by="Strain",all.x = TRUE)

  plot(jitter(C8$A17_C68)~C8$A17_C8,pch=19,col=paste(C8$colors))
  points(jitter(C8$R108_C68)~C8$R108_C8,pch=17,col=paste(C8$colors))


pdf(file="figures/Figure4b_rankC8inC68.pdf", height=5,width=5,useDingbats = FALSE)
  plot(jitter(C8$A17_C68,factor=.5)~C8$A17_C8,pch=19,col=paste(C8$colors),xlab="Rank in C8",ylab="Rank in C68",axes=FALSE)
    points(jitter(C8$R108_C68,factor=.5)~C8$R108_C8,pch=17,col=paste(C8$colors))  
    axis(side = 1,at = c(1:8))
    axis(side = 2 )
    legend("topleft",legend=c("A17","R108"),pch = c(19,17))
dev.off()

pdf(file="figures/Figure4a_rankC3inC68.pdf", height=4,width=4,useDingbats = FALSE)
  plot(jitter(C3$A17_C68,factor=.5)~C3$A17_C3,pch=19,col=paste(C3$colors),xlab="Rank in C3",ylab="Rank in C68",axes=FALSE)
    points(jitter(C3$R108_C68,factor=.5)~C3$R108_C3,pch=17,col=paste(C3$colors))  
    axis(side = 1,at = c(1:3))
    axis(side = 2, at = c(1:3))
    legend("topleft",legend=c("A17","R108"),pch = c(19,17))
dev.off()


##### Initial Community CFU's (Figure S1) #####
library(plyr)
Initial<-read.csv(file="data/ColonyCounts_InitialCommunities_2Feb2018.txt",sep="\t",header = TRUE)
Initial<-Initial[!Initial$Strain %in% c("SM1021"),]

meso.i<-ddply(Initial, .(Strain, Letter,Dilution), summarize,
              mean = mean(NumColonies))
meso.i$Strain<-factor(meso.i$Strain,levels=c("C68","C8","C3","HK","HM","MK","C68L"))

pdf(file="figures/FigureS1_InitialCommunityCFU.pdf",height=3,width=4, useDingbats = FALSE)
ggplot(meso.i,aes(x=Strain,y=mean*10^Dilution,color=Strain))+
  geom_point()+
  scale_y_log10()+
  xlab("Strain Community")+
  ylab("Colony Forming Units (CFU) per ml innoculum")+
  theme_minimal()
dev.off()

##### PCA's (Figure S4) ########
library("FactoMineR")
library("factoextra")

pdf(file="figures/FigureS5_PCAFitness.pdf",height=3,width=4, useDingbats = FALSE)
myPCA<-PCA(fitness[,-1:-9], scale.unit = FALSE, graph = FALSE)
x<-fviz_pca_ind(myPCA, geom.ind = "point", col.ind = fitness$Host_trt, axes = c(1,2),
                addEllipses = TRUE, ellipse.type = "confidence",
                legend.title = "Groups",mean.point=FALSE)
x + scale_fill_manual(values = rep(mycols,2))+ scale_color_manual(values = rep(mycols,2))

myPCA<-PCA(fitness_A17[,-1:-9], scale.unit = FALSE, graph = FALSE)
x<-fviz_pca_ind(myPCA, geom.ind = "point", col.ind = fitness_A17$Host_trt, axes = c(1,2),
                addEllipses = TRUE, ellipse.type = "confidence",
                legend.title = "Groups",mean.point=FALSE)
x+ scale_fill_manual(values = mycols)+ scale_color_manual(values = mycols)

myPCA<-PCA(fitness_R108[,-1:-9], scale.unit = FALSE, graph = FALSE)
x<-fviz_pca_ind(myPCA, geom.ind = "point", col.ind = fitness_R108$Host_trt, axes = c(1,2),
                addEllipses = TRUE, ellipse.type = "confidence",
                legend.title = "Groups",mean.point=FALSE)
x+ scale_shape_manual(values = c(5:8)) + scale_fill_manual(values = mycols)+ scale_color_manual(values = mycols)
dev.off()

#### Calculating Sequencing Coverages cited in methods of Manuscript ####

coverage<-read.table(file="data/winter_2018_pools_USDA1106_read_counts_2019-09-24.tsv",header = TRUE,sep = "\t")

coverage_A17_Eco<-rbind(coverage[grep(pattern = "A17_X",x = coverage$pool),], coverage[grep(pattern = "A17_N",x = coverage$pool),], coverage[grep(pattern = "A17_L",x = coverage$pool),])
coverage_A17_Eco<-coverage_A17_Eco[-grep(pattern = "_D",x = coverage_A17_Eco$pool),]

coverage_R108_Eco<-rbind(coverage[grep(pattern = "R108_X",x = coverage$pool),], coverage[grep(pattern = "R108_N",x = coverage$pool),], coverage[grep(pattern = "R108_L",x = coverage$pool),])
coverage_R108_Eco<-coverage_R108_Eco[-grep(pattern = "_D",x = coverage_R108_Eco$pool),]

summary(c((coverage_A17_Eco$aln_n_reads*125)/6716230,(coverage_A17_Eco$aln_n_reads*125)/6716230))

summary(c(coverage_A17_Eco$raw_n_pairs,coverage_A17_Eco$raw_n_pairs))

######## Heatmap + Phylogeny ########

### Phylogeny Ordered Heatmap @ 24 wks

# Prepare the data-Calculate median strain fitness for each treatment @ 24wks and create it into a single dataframe with a column for each treatment and a row for each strain. 

fitness_long<- as_tibble(fitness) %>% pivot_longer(cols = USDA1157:X1719,names_to= "strain",values_to = "fitness")

fitness_med<- fitness_long %>% group_by(Host_trt,strain) %>% summarise_if(is.numeric,median) %>% pivot_wider(names_from = Host_trt,values_from = fitness)
fitness_med$strain<- as.character(fitness_med$strain)
fitness_med$strain <- gsub("X","",fitness_med$strain) # get rid of the X's in front of the strains
fitness_med$strain <- gsub("USDA","",fitness_med$strain) # get rid of the USDA's in front of the strains



#library (gplots)
library(ape)
## Read in the tree file created by Brendan
tree = read.tree('data/tree.nw')

### there are some discrepancies in strain names that need to be fixed...
tree$tip.label[tree$tip.label=="KH46c"]<-"KH46C"
tree$tip.label[tree$tip.label=="HM006-1"]<-"HM006.1"
tree$tip.label[tree$tip.label=="KH35c"]<-"KH35C"
tree$tip.label[tree$tip.label=="USDA1021"]<-"1021"
tree$tip.label[tree$tip.label=="USDA1157"]<-"1157"

# load in teh function for making a heatmap with the tree #
heatmap.phylo <- function(x, Rowp, Colp, ...) {
  l = length(seq(-8, 3.9, 0.1))
  pal = colorRampPalette(c('#2166ac','#92c5de', '#f7f7f7', '#b2182b'))(l)
  row_order = Rowp$tip.label[Rowp$edge[Rowp$edge[, 2] <= Ntip(Rowp), 2]] 
  col_order = Colp$tip.label[Colp$edge[Colp$edge[, 2] <= Ntip(Colp), 2]] 
  x <- x[row_order, col_order]
  xl <- c(0.5, ncol(x)+0.5)
  yl <- c(0.5, nrow(x)+0.5)
  layout(matrix(c(0,1,0, 2,3,4, 0,5,0), nrow=3, byrow=TRUE),
         width=c(3.5,    4.5, 1),
         height=c(0.2, 3, 0.18))
  par(mar=rep(0,4))
  plot(Colp, direction="downwards", show.tip.label=FALSE,
       xaxs="i", x.lim=xl)
  par(mar=rep(0,4))
  plot(Rowp, direction="rightwards", show.tip.label=FALSE, 
       yaxs="i", y.lim=yl)
  lpp = .PlotPhyloEnv$last_plot.phylo 
  segments(lpp$xx[1:Ntip(Rowp)], lpp$yy[1:Ntip(Rowp)], par('usr')[2],
           lpp$yy[1:Ntip(Rowp)], lty=3, col='grey50')
  par(mar=rep(0,4), xpd=TRUE)
  image((1:ncol(x))-0.5, (1:nrow(x))-0.5, t(x), col=pal,
        xaxs="i", yaxs="i", axes=FALSE, xlab="",ylab="",breaks=seq(-8,4,0.1))
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", yaxs="i", xlim=c(0,2), ylim=yl)
  text(rep(0,nrow(x)),1:nrow(x), row_order, pos=4, family='Helvetica',
       cex=1, xpd=NA)
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", xaxs="i", ylim=c(0,2), xlim=xl)
  text(1:ncol(x),rep(2,ncol(x)), col_order, srt=90, adj=c(1,.5), family='Helvetica',
       cex=1.5)
}

# Create the matrix and get the column dendrogram for the heatmap from it.
m = structure(as.matrix(fitness_med[, -1:-2]),
              dimnames=list(unlist(fitness_med[, 'strain']),
                            unlist(names(fitness_med)[-1:-2])))

col_dendro = as.dendrogram(hclust(dist(t(m))))

# And make the plot....
pdf(file="figures/FigX_heatmap_Phylo.pdf",width = 5,height=8, useDingbats=FALSE)
heatmap.phylo(x = m, Rowp = tree, Colp = as.phylo(as.hclust(col_dendro)))
dev.off()

l = length(seq(-8, 3.9, 0.1))
pal = colorRampPalette(c('#2166ac','#92c5de', '#f7f7f7', '#b2182b'))(l)

pdf(file="figures/FigX_heatmap_legend.pdf",height = 3, width = 6)
legend_image <- as.raster(matrix(pal, nrow=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
#text(y=-0.2, x = seq(0,2,l=4), labels = seq(-8,4,l=4))
rasterImage(legend_image, 0, 0, 2,2)
dev.off()
