## This is a script to perform the analysis of data presented in Quick et al., 2019

## Load useful libaries
library(ggplot2)
library(lme4)
library(plyr)
library(dplyr)
library(data.table)
library(multcomp) ## to use glht()
library(car)
library(Rfast)
library(plotrix)
library(MASS)
library(lmtest)

################################################


## The first thing is to set your R session to the base directory you just downloaded from github
## insert path below...

setwd("/Users/max.feldman/Documents/github/r.quick_castle_trv_cleanse_2019")

setwd()


##### CREATE DIRECTORY PATHS ##### 

## Make the directory of the folder you downloaded the current working directory
home.dir<-getwd()
setwd(home.dir)

## Read in data
data<-read.csv("r.quick_trv_castle_cleanse_final.csv")


## All nematode counts are based upon the number of nematodes present in 250 cc of soil
## Samples are taken from 10 liter pots

## Calculation of Rf (Reproduction factor) for nematodes
## Rf = # of nematodes now / # of nematodes to start with

## Experiment 1 began with 60 nematodes / 10,000 cm of soil 
## 60 / 40 = 1.5 nematodes / 250 cm of soil

## Experiment 2 began with 1060 / 10,000 cm of soil
## 1060 / 40 = 26.5 nematodes / 250 cm of soil

## Lets add this column to data
data$Rf<-rep('NA', nrow(data))
data[data$experiment == 1, 'Rf']<-((data[data$experiment == 1, 'nematode'])/(60/40))
data[data$experiment == 2, 'Rf']<-((data[data$experiment == 2, 'nematode'])/(1060/40))

## lets look at nematode counts on non-bait plants
data.non.bait<-data[data$bait == 'no',]
## Convert Rf to numeric from character
data.non.bait$Rf<-as.numeric(as.character(data.non.bait$Rf))

## Lets get a mean for each variable category
data.non.bait.ag <- aggregate(data.non.bait[,c('plant', 'experiment', 'nematode','month','Rf')], by=list(data.non.bait$plant, data.non.bait$experiment, data.non.bait$month), mean)
data.non.bait.ag<-data.non.bait.ag[,-c(2:4)]
colnames(data.non.bait.ag)[1]<-c("plant")

################################################

## Lets model nematode count (first using poission)
nematode.poi.mdl<-glm(nematode ~ plant + month + experiment, data=data.non.bait, family=poisson)
summary(nematode.poi.mdl)

## Model does not fit
1 - pchisq(summary(nematode.poi.mdl)$deviance, 
           summary(nematode.poi.mdl)$df.residual
)

## P-value is < 0.05

Anova(nematode.poi.mdl, type=c("III"), test.statistic="LR")
nematode_anova.poi.table<-as.data.frame(Anova(nematode.poi.mdl, type=c("III"), test.statistic="LR"))

## Now with negative binomial
nematode.nb.mdl<-glm.nb(nematode ~ plant + month + experiment, data=data.non.bait)
summary(nematode.nb.mdl)

## Model fits less poorly cannot reject the null
1 - pchisq(summary(nematode.nb.mdl)$deviance, 
           summary(nematode.nb.mdl)$df.residual
)

## p-value = 0.08821647 which is > 0.05

Anova(nematode.nb.mdl, type=c("III"), test.statistic="LR")
nematode_anova.nb.table<-as.data.frame(Anova(nematode.nb.mdl, type=c("III"), test.statistic="LR"))

write.csv(nematode_anova.nb.table, file="Table_1a.csv", quote=F)

## lets try without tobacco
no.tobacco.nematode.mdl<-glm(nematode ~ plant + month + experiment, data=data.non.bait[data.non.bait$plant != 'tobacco',], family=poisson)
summary(no.tobacco.nematode.mdl)

Anova(no.tobacco.nematode.mdl, type=c("III"), test.statistic="LR")

1 - pchisq(summary(no.tobacco.nematode.mdl)$deviance, 
           summary(no.tobacco.nematode.mdl)$df.residual
)

## Lets model Rf value (derived data is not count so we will not use poisson distribution)
## data is not normal
hist(data.non.bait$Rf)

## log transformation makes the data ~normal
hist(log(data.non.bait$Rf))

## Get an offset value that enables you to use log transformation
## The value is 1/2 the minimum value
offset<-min(data.non.bait[data.non.bait$Rf > 0,"Rf"])/2

Rf.mdl<-lm(log(Rf + offset) ~ plant + month + experiment, data=data.non.bait)
summary(Rf.mdl)


Anova(Rf.mdl, type=c("III"), test.statistic="LR")
rf_anova.table<-as.data.frame(Anova(Rf.mdl, type=c("III"), test.statistic="LR"))

write.csv(rf_anova.table, file="Table_1b.csv", quote=F)

## Lets make some plots of the data

## Nematode count
p<-ggplot(data.non.bait, aes(x=month, y=nematode, color=factor(experiment))) + geom_point() + facet_wrap(~plant) 
q<-p + geom_point(data=data.non.bait.ag, aes(x=month, y=nematode,color=factor(experiment), size=4))
s<-q + geom_line(data=data.non.bait.ag, aes(x=month, y=nematode,color=factor(experiment)))
s


t<-ggplot(data.non.bait.ag, aes(x=month, y=nematode, color=factor(experiment))) + geom_point() + facet_wrap(~plant) 
t


p<-ggplot(data.non.bait, aes(x=month, y=Rf, color=factor(experiment))) + geom_point() + facet_wrap(~plant) 
q<-p + geom_point(data=data.non.bait.ag, aes(x=month, y=Rf,color=factor(experiment), size=4))
s<-q + geom_line(data=data.non.bait.ag, aes(x=month, y=Rf,color=factor(experiment)))
s


## Lets make a summary table that describes the distribtuion (mean, median, and standard error of the mean)
data.non.bait.mean <- aggregate(data.non.bait[,c('plant', 'experiment', 'nematode','month','Rf')], by=list(data.non.bait$plant, data.non.bait$experiment, data.non.bait$month), mean)
data.non.bait.mean<-data.non.bait.mean[,-c(2:4)]
data.non.bait.mean<-data.non.bait.mean[,c(1,2,4,3,5)]
colnames(data.non.bait.mean)<-c("Plant", "Experiment","Month","Nematode_count","Rf_mean" )

data.non.bait.median <- aggregate(data.non.bait[,c('plant', 'experiment', 'nematode','month','Rf')], by=list(data.non.bait$plant, data.non.bait$experiment, data.non.bait$month), med)
data.non.bait.median<-data.non.bait.median[,-c(2:4)]
data.non.bait.median<-data.non.bait.median[,c(1,2,4,3,5)]
colnames(data.non.bait.median)<-c("Plant", "Experiment","Month","Nematode_median","Rf_median" )

data.non.bait.stderr <- aggregate(data.non.bait[,c('plant', 'experiment', 'nematode','month','Rf')], by=list(data.non.bait$plant, data.non.bait$experiment, data.non.bait$month), std.error)
data.non.bait.stderr<-data.non.bait.stderr[,-c(4,5,7)]
colnames(data.non.bait.stderr)<-c("Plant", "Experiment","Month","Nematode_stderr","Rf_stderr" )

data.non.bait.summary<-merge(data.non.bait.mean, data.non.bait.median, by=c("Plant", "Experiment", "Month"))
data.non.bait.summary<-merge(data.non.bait.summary, data.non.bait.stderr, by=c("Plant", "Experiment", "Month"))

data.non.bait.summary$Plant<-as.character(data.non.bait.summary$Plant)
data.non.bait.summary[data.non.bait.summary$Plant == 'alfalfa', "Plant"]<-c("Alfalfa")
data.non.bait.summary[data.non.bait.summary$Plant == 'burbank', "Plant"]<-c("Burbank")
data.non.bait.summary[data.non.bait.summary$Plant == 'tobacco', "Plant"]<-c("Tobacco")
data.non.bait.summary[data.non.bait.summary$Plant == 'castle', "Plant"]<-c("Castle")

data.non.bait.summary<-data.non.bait.summary[with(data.non.bait.summary, order(Plant, Experiment, Month)),]
data.non.bait.summary<-data.non.bait.summary[,c(1:4,6,8,5,7,9)]

write.csv(data.non.bait.summary, file="Table_S1.csv", quote=F, row.names=F)

data.non.bait$plant<-as.character(data.non.bait$plant)
data.non.bait[data.non.bait$plant == "alfalfa", "plant"]<-c("Alfalfa")
data.non.bait[data.non.bait$plant == "castle", "plant"]<-c("Castle")
data.non.bait[data.non.bait$plant == "burbank", "plant"]<-c("Burbank")
data.non.bait[data.non.bait$plant == "tobacco", "plant"]<-c("Tobacco")
colnames(data.non.bait)[c(2:4,7)]<-c("Plant", "Experiment", "Nematode_count", "Month")


## Figure 1

p<-ggplot(data.non.bait, aes(x=Month, y=Nematode_count, color=factor(Experiment))) + geom_point(size=0.2) + facet_wrap(~Plant, scales="free_y") + scale_color_manual(values=c("Orange", "Red"),name = "Experiment")
q<-p+geom_point(data=data.non.bait.summary, aes(x=Month, y=Nematode_count, color=factor(Experiment)), shape=4, size=2)
s<-q+geom_errorbar(data=data.non.bait.summary, aes(ymin=Nematode_count-Nematode_stderr, ymax=Nematode_count+Nematode_stderr), width=.4,position=position_dodge(0.05))
t<-s+theme_bw() + ylab("Nematode Count (per 250 cc)")
t 

pdf("Figure_1.pdf")
print(t)
dev.off()


## Figure S1

## Need to change name from Rf_mean to Rf so we can make the plot
colnames(data.non.bait.summary)[7]<-c("Rf")

p<-ggplot(data.non.bait, aes(x=Month, y=Rf, color=factor(Experiment))) + geom_point(size=0.2) + facet_wrap(~Plant, scales="free_y") + scale_color_manual(values=c("Orange", "Red"),name = "Experiment")
q<-p+geom_point(data=data.non.bait.summary, aes(x=Month, y=Rf, color=factor(Experiment)), shape=4, size=2)
s<-q+geom_errorbar(data=data.non.bait.summary, aes(ymin=Rf - Rf_stderr, ymax=Rf + Rf_stderr), width=.4, position=position_dodge(0.05))
t<-s+theme_bw() + ylab("Rf")
t

pdf("Figure_S1.pdf")
print(t)
dev.off()




################################################

# lets look at nematode counts on non-bait plants
data.bait<-data[data$bait == 'yes' | data$plant == "tobacco",]
## Convert Rf to numeric from character
data.bait$Rf<-as.numeric(as.character(data.bait$Rf))

## Lets get a mean for each variable category
data.bait.ag <- aggregate(data.bait[,c('plant', 'experiment', 'nematode','month','Rf')], by=list(data.bait$plant, data.bait$experiment, data.bait$month), mean)
data.bait.ag<-data.bait.ag[,-c(2:4)]
colnames(data.bait.ag)[1]<-c("plant")

## get RT-PCR count data
root_count<-as.data.frame(table(data[,c("plant","month","bait", "experiment", "root")]))
shoot_count<-as.data.frame(table(data[,c("plant","month","bait", "experiment", "foliage")]))
colnames(root_count)[which(colnames(root_count) == 'Freq')]<-c("root_count")
colnames(shoot_count)[which(colnames(shoot_count) == 'Freq')]<-c("shoot_count")

## lets keep data.frame with all entries
root_count_all<-root_count
shoot_count_all<-shoot_count

## lets get a total number of plants (positive and negative) so we can run some correlation statistics

## root first
root_count_pos<-root_count[root_count$root == "positive",]
root_count_pos<-root_count_pos[,-c(5)]
colnames(root_count_pos)[5]<-c("root_count_pos")

root_count_neg<-root_count[root_count$root == "negative",]
root_count_neg<-root_count_neg[,-c(5)]
colnames(root_count_neg)[5]<-c("root_count_neg")

root_count_total<-merge(root_count_pos, root_count_neg, by=c("plant", "month", "experiment", "bait"))
root_count_total$root_total<-root_count_total$root_count_pos + root_count_total$root_count_neg
root_count_total<-root_count_total[root_count_total$root_total != 0,]
root_count_total$root_percent_pos<-root_count_total$root_count_pos/root_count_total$root_total

## now shoot
shoot_count_pos<-shoot_count[shoot_count$foliage == "positive",]
shoot_count_pos<-shoot_count_pos[,-c(5)]
colnames(shoot_count_pos)[5]<-c("shoot_count_pos")

shoot_count_neg<-shoot_count[shoot_count$foliage == "negative",]
shoot_count_neg<-shoot_count_neg[,-c(5)]
colnames(shoot_count_neg)[5]<-c("shoot_count_neg")
shoot_count_total<-merge(shoot_count_pos, shoot_count_neg, by=c("plant", "month", "experiment", "bait"))
shoot_count_total$shoot_total<-shoot_count_total$shoot_count_pos + shoot_count_total$shoot_count_neg
shoot_count_total<-shoot_count_total[shoot_count_total$shoot_total != 0,]
shoot_count_total$shoot_percent_pos<-shoot_count_total$shoot_count_pos/shoot_count_total$shoot_total

rtpcr_counts<-merge(root_count_total, shoot_count_total, by=c("plant", "month", "experiment", "bait"))


root_count_all$tissue<-rep("root", nrow(root_count_all))
colnames(root_count_all)[5]<-c("rt.pcr")
colnames(root_count_all)[6]<-c("count")

shoot_count_all$tissue<-rep("shoot", nrow(shoot_count_all))
colnames(shoot_count_all)[5]<-c("rt.pcr")
colnames(shoot_count_all)[6]<-c("count")

count_data_all<-rbind(root_count_all, shoot_count_all)
count_data_all$plant<-factor(count_data_all$plant, levels=c("tobacco","alfalfa","burbank","castle"))
count_data_all$month<-factor(count_data_all$month, levels=c("1","2","3","4"))
count_data_all$bait<-factor(count_data_all$bait, levels=c("no","yes"))
count_data_all$experiment<-factor(count_data_all$experiment, levels=c("1","2"))

################################################
## fit model

## Recode as numbers for logistic anova
tester<-data.bait
tester$foliage<-as.character(tester$foliage)
tester$root<-as.character(tester$root)

tester[tester$foliage == "positive", "foliage"]<-c(1)
tester[tester$foliage == "negative", "foliage"]<-c(0)
tester[tester$root == "positive", "root"]<-c(1)
tester[tester$root == "negative", "root"]<-c(0)
tester$root<-as.numeric(as.character(tester$root))
tester$foliage<-as.numeric(as.character(tester$foliage))

trv.root.mdl<-glm(root ~ plant + month + experiment, data=tester, family="binomial")
summary(trv.root.mdl)
Anova(trv.root.mdl, type=c("III"), test.statistic="LR")

trv.root_anova.table<-as.data.frame(Anova(trv.root.mdl, type=c("III"), test.statistic="LR"))

write.csv(trv.root_anova.table, file="Table_2a.csv", quote=F)

trv.shoot.mdl<-glm(foliage ~ plant + month + experiment, data=tester, family="binomial")
summary(trv.shoot.mdl)
Anova(trv.shoot.mdl, type=c("III"), test.statistic="LR")

trv.shoot_anova.table<-as.data.frame(Anova(trv.shoot.mdl, type=c("III"), test.statistic="LR"))

write.csv(trv.shoot_anova.table, file="Table_2b.csv", quote=F)


################################################
## Lets make some plots of the data

count_data_all[count_data_all$rt.pcr == "positive", 'count']<-count_data_all[count_data_all$rt.pcr == "positive", 'count'] * -1

## Plots
count_data_all_bait<-count_data_all[count_data_all$bait == 'yes' | count_data_all$plant == 'tobacco',]

## Lets capitalize the names so the plot looks pretty
count_data_all_bait$plant<-as.character(count_data_all_bait$plant)
count_data_all_bait[count_data_all_bait$plant == "tobacco", "plant"]<-c("Tobacco")
count_data_all_bait[count_data_all_bait$plant == "alfalfa", "plant"]<-c("Alfalfa")
count_data_all_bait[count_data_all_bait$plant == "burbank", "plant"]<-c("Burbank")
count_data_all_bait[count_data_all_bait$plant == "castle", "plant"]<-c("Castle")

count_data_all_bait$rt.pcr<-as.character(count_data_all_bait$rt.pcr)
count_data_all_bait[count_data_all_bait$rt.pcr == "negative", "rt.pcr"]<-c("Negative")
count_data_all_bait[count_data_all_bait$rt.pcr == "positive", "rt.pcr"]<-c("Positive")

## Experiment 1 shoot tissue
p<-ggplot(count_data_all_bait[count_data_all_bait$experiment == 1 & count_data_all_bait$tissue == 'shoot' & as.numeric(as.character(count_data_all_bait$month)) > 1,], aes(x=month, y=count, color=rt.pcr, fill=rt.pcr)) + geom_bar(stat="identity") + facet_grid(.~plant)
q<-p + scale_color_manual(values=c("darkgreen", "green"), guide="none") + scale_fill_manual(values=c("darkgreen", "green"), name="RT-PCR") + theme_bw() + ylab("RT-PCR Count") + xlab("Months") + ggtitle("Experiment 1")

pdf("Figure_2a.pdf")
print(q)
dev.off()

## Experiment 1 root tissue
p<-ggplot(count_data_all_bait[count_data_all_bait$experiment == 1 & count_data_all_bait$tissue == 'root' & as.numeric(as.character(count_data_all_bait$month)) > 1,], aes(x=month, y=count, color=rt.pcr, fill=rt.pcr)) + geom_bar(stat="identity") + facet_grid(.~plant)
q<-p + scale_color_manual(values=c("brown", "orange"), guide="none") + scale_fill_manual(values=c("brown", "orange"), name="RT-PCR") + theme_bw() + ylab("RT-PCR Count") + xlab("Months") + ggtitle("Experiment 1")

pdf("Figure_2b.pdf")
print(q)
dev.off()

## Experiment 2 shoot tissue
p<-ggplot(count_data_all_bait[count_data_all_bait$experiment == 2 & count_data_all_bait$tissue == 'shoot' & as.numeric(as.character(count_data_all_bait$month)) > 1,], aes(x=month, y=count, color=rt.pcr, fill=rt.pcr)) + geom_bar(stat="identity") + facet_grid(.~plant)
q<-p + scale_color_manual(values=c("darkgreen", "green"), guide="none") + scale_fill_manual(values=c("darkgreen", "green"), name="RT-PCR") + theme_bw() + ylab("RT-PCR Count") + xlab("Months") + ggtitle("Experiment 2")

pdf("Figure_2c.pdf")
print(q)
dev.off()

## Experiment 2 root tissue
p<-ggplot(count_data_all_bait[count_data_all_bait$experiment == 2 & count_data_all_bait$tissue == 'root' & as.numeric(as.character(count_data_all_bait$month)) > 1,], aes(x=month, y=count, color=rt.pcr, fill=rt.pcr)) + geom_bar(stat="identity") + facet_grid(.~plant)
q<-p + scale_color_manual(values=c("brown", "orange"), guide="none") + scale_fill_manual(values=c("brown", "orange"), name="RT-PCR") + theme_bw() + ylab("RT-PCR Count") + xlab("Months") + ggtitle("Experiment 2")

pdf("Figure_2d.pdf")
print(q)
dev.off()

## Lets examine the relationship between SRN number and TRV bait plant infection
data$Rf<-as.numeric(as.character(data$Rf))
data.nematode.ag <- aggregate(data[,c('plant', 'experiment', 'nematode','month','Rf')], by=list(data$plant, data$experiment, data$month, data$bait), mean)

data.nematode.ag<-data.nematode.ag[,-c(2,3,5)]
colnames(data.nematode.ag)[1:2]<-c('plant', 'bait')

## Here we speccify tobacco as a bait plant 
data.nematode.ag<-data.nematode.ag[data.nematode.ag$bait == "no" | data.nematode.ag$plant == "tobacco",]


## Because nematode counts are only performed months 1-3 and RT-PCR is done on bait plants in months 1-4 
## Lets add the initial nematode concentrations as month 1 and shift all nematode counts to one month later so they 
## are matching

#nematode_count_initial<-data.frame("plant"=c("alfalfa", "burbank", "castle", "tobacco","alfalfa", "burbank", "castle", "tobacco"),"experiment"=c(1,1,1,1,2,2,2,2), "nematode"=c(4.5, 4.5, 4.5, 4.5, 79.5, 79.5, 79.5, 79.5), "month"=rep(1, 8), 'Rf'=rep(1,8))
data.nematode.ag<-data.nematode.ag[,-c(2)]         
data.nematode.ag$month<-data.nematode.ag$month + 1
#data.nematode.ag<-rbind(nematode_count_initial,data.nematode.ag)
data.nematode.ag<-data.nematode.ag[data.nematode.ag$month < 5,]

rtpcr_counts_bait<-rtpcr_counts[rtpcr_counts$bait == "yes" | rtpcr_counts$plant == "tobacco",]
rtpcr_counts_bait<-rtpcr_counts_bait[,-c(4)]

rtpcr_counts_bait<-rtpcr_counts_bait[rtpcr_counts_bait$month != 1,]

## Merge both
count_and_nematode<-merge(rtpcr_counts_bait, data.nematode.ag,  by=c("plant", "month", "experiment"))

## check correlation
cor.test(count_and_nematode$root_percent_pos, count_and_nematode$nematode)
cor.test(count_and_nematode$shoot_percent_pos, count_and_nematode$nematode)
cor.test(count_and_nematode$root_percent_pos, count_and_nematode$Rf)
cor.test(count_and_nematode$shoot_percent_pos, count_and_nematode$Rf)

## No significant correlation between nematode count, Rf value and % of bait plants testing positive by RT-PCR

ggplot(data=count_and_nematode, aes(y=root_percent_pos, x=log(nematode))) + geom_point()
