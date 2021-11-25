library(readxl)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(lme4)
library(glmmTMB)
library(AER)
library(AICcmodavg)
library(glmmML)
require(MuMIn)

### The analyses here require several different datasets. PART ONE deals with transect data on the basic ecological characeristics of A. melanantha (Q1). PART TWO addressess patterns of seeds on the fur of geladas (Q2-Q4). PART THREE uses data from focal follows on short-term fluxes of seeds (Q5).


## PART ONE. This portion of the script deals with the transect performed in September and October 2018 to determine basic aspects of A. melanantha ecology. These correspond to Q1 in the manuscript ##

trans <- read_excel("~/Desktop/transect.xlsx",col_types = c("date", "text", "numeric","numeric", "numeric", "numeric", "numeric", "text"))
trans<-as.data.frame(trans)
head(trans)

length(trans$Num_Agro) # 180 rows
length(which(trans$Num_Agro==0)) # 135 rows with zero agrocharis

length(which(trans$Num_Agro==0))/length(trans$Num_Agro) # 75% (135/180) of plots had zero Agrocharis, and those that did have between 2 and 8 plants in them 

#### analysis by microhabitat ####
habitat<-trans %>% group_by(Hab) %>% summarize(numb=sum(Num_Agro),prop=numb/135) # 87% of plants were in predominately short grass habitat, although these habitats were often complex mixes of grasses and herbs. 9% were found in tall grass habitat. 
habitat
sum(habitat$numb) # 135 is the total number of plants counted across all the transects

trans$flow<-trans$Num_Flow_Agro/trans$Num_Agro
hist(trans$flow)
flowering<-trans %>% group_by(Hab) %>% summarize(numb=mean(flow,na.rm=T))
flowering # fewer than a quarter of the plants are actually flowering at this time

## This is the dataframe with burr counts on the umbels of the plants

bur <- read_excel("~/Desktop/burr.xlsx",col_types = c("date", "text", "numeric", "numeric", "numeric", "numeric",   "numeric", "numeric", "numeric","text"))
bur<-as.data.frame(bur)
head(bur)

# get average heights etc. each row here is a plant in a plot. the columns are averaged values across all the umbels on that plant, so it's average plant height. number of seeds is summed up
summ1<-bur %>% group_by(Plot, Plant) %>% summarize(ht=mean(Tot_Ht), ht_can=mean(Ht_Canopy),sds=sum(Num_Sds))
head(summ1,20)

### Absolute height
summary(summ1$ht) # min 5cm, max 37cm, mean 16.1 cm.
length(summ1$ht) # 33 plots, each of which is  plant
sd(summ1$ht) # 7.1 sd

## Relative canopy height 
summary(summ1$ht_can) # min is -7cm, max is 13.8 cm, mean is 2.0cm
length(summ1$ht_can) # n=33, each of which is a plant
sd(summ1$ht_can) #sd is 4.8

### Number of seeds per plant
summary(summ1$sds) # min 1, mean 25, max 95
length(summ1$sds) # 33
sd(summ1$sds) # sd = 24
median(summ1$sds) # median is 16

### Average seeds per umbel
mean(bur$Num_Sds) # mean of 5.4
length(bur$Num_Sds) # 153 umbels
summary(bur$Num_Sds) # min 1, mean 5.4, max 16 umbels
median(bur$Num_Sds) # median is 5
sd(bur$Num_Sds) # SD is 2.9

### PART TWO. This portion of the analysis deals with the seeds on the fur of geladas (Q2-Q4) ##

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # grey

'%ni%' <- Negate('%in%')

sheet<-read.csv("/Users/vivekvasivenkataraman/Desktop/agrosheet.csv")
sheet<-as.data.frame(sheet)
sheet$SNAME<-as.factor(sheet$SNAME)
sheet<-sheet %>% mutate(newdate=dmy(DATE)) %>% mutate(month=month(newdate),yr=year(newdate),dy=day(newdate))# add column with more standard notation for date and add columns that pull out month, year and date# new variable to make new age categories based on size. Only unintuitive assignment is LJs being included as adult females.

sheet<-sheet %>% filter(SNAME!="LL1")
sheet<-sheet %>% filter(SNAME!="LL2")
sheet<-sheet %>% filter(SNAME!="LTD")

sheet$newage<-0 # new column for re-labeling age-sex classes
sheet$newage[which(sheet$AGE=="SJ" | sheet$AGE=="MJ")]<-"Small"
sheet$newage[which(sheet$AGE=="BRI" | sheet$AGE=="BLI")]<-"Very small"
sheet$newage[which(sheet$AGE=="A" | sheet$AGE=="SAM"| sheet$AGE=="SA" & sheet$SEX=="M")]<-"Large"
sheet$newage[which(sheet$AGE=="A" & sheet$SEX=="F")]<-"Medium"
sheet$newage[which(sheet$AGE=="LJ" & sheet$SEX=="M")]<-"Medium"
sheet$newage[which(sheet$AGE=="LJ" & sheet$SEX=="F")]<-"Medium"
sheet$newage<-as.factor(sheet$newage)
sheet$newage <- factor(sheet$newage, order=TRUE, levels = c("Very small", "Small","Medium","Large")) # order factor for body size

### for body size as continuous variable
sheet$size<-0
sheet$size[sheet$newage=="Very small"]<-3.5
sheet$size[sheet$newage=="Small"]<-7
sheet$size[sheet$newage=="Medium"]<-10
sheet$size[sheet$newage=="Large"]<-19
sheet$size<-as.numeric(sheet$size)
sheet<-sheet %>% drop_na(newage) # for forthcoming GLMM, drop unknown age-sex class

### Import dataset on OMU unit names and merge with sheet

units<-read.csv("~/Desktop/units.csv")

sheet <- left_join(sheet,units,by="SNAME")

sheet<-sheet %>% drop_na(OMU) # drop unknown OMUs
sheet[is.na(sheet)] <- 0 #change NA's to zeros. these are real zeros.

## A few basic attributes of the dataset
length(sheet$TOTAL) # 5344 rows
reps<-sheet %>% group_by(SNAME,newdate) %>% summarize(number=n()) # show repeated sampling of individuals
table(reps$number) # 5100 samples were once per day, 116 instances of twice per day, 4 instances of 3x per day. So, in 98% of cases observations were made on animals once per day
max(sheet$TOTAL) # 511 burrs on a single individual in a day
length(unique(sheet$SNAME)) #  225 unique inds
sum(sheet$TOTAL) # 12649 seeds counted overall

#look at how many of each age sex class were sampled
agesex<-sheet %>% group_by(month,newage) %>% summarize(number=n(),avg=mean(TOTAL),maxim=max(TOTAL)) # show repeated sampling of individuals
head(agesex)
agesex<-as.data.frame(agesex)
b<-agesex %>% group_by(newage) %>% summarize(ss=mean(number))
b<-as.data.frame(b)
b # avg sample sizes across months

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# Generalized linear mixed models (GLMMs) to examine determinants of seed counts
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

## First, need to examine issue of zero-inflation in the dataset
sum(sheet$TOTAL > 0) # 1539 nonzero
sum(sheet$TOTAL == 0) # 3805 are zero, so, 60% of dataset is zeros
# plot(table(sheet$TOTAL)) # visual representation of zero-inflation

sheet<- sheet %>% filter(month!=3) 
sheet <- sheet %>% filter(month!=4)

#turn some variables into factors for analysis
sheet$yr<-as.factor(sheet$yr)
sheet$OMU<-as.factor(sheet$OMU)
sheet$month<-as.factor(sheet$month)

# List of candidate models
cand.models <- list()

cand.models[[1]] <- glmer(TOTAL ~ size + yr + month + (1|OMU/SNAME), data = sheet,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[2]] <- glmer(TOTAL ~ size + month + (1|OMU/SNAME), data = sheet,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[3]] <- glmer(TOTAL ~ size + yr + (1|OMU/SNAME), data = sheet,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[4]] <- glmer(TOTAL ~ yr + month + (1|OMU/SNAME), data = sheet,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[5]] <- glmer(TOTAL ~ yr + (1|OMU/SNAME), data = sheet,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[6]] <- glmer(TOTAL ~ month  + (1|OMU/SNAME), data = sheet,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[7]] <- glmer(TOTAL ~ size + (1|OMU/SNAME), data = sheet,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[8]] <- glmer(TOTAL ~ (1|OMU/SNAME), data = sheet,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

Modnames <- paste("mod", 1:length(cand.models), sep = " ")

##generate AICc table
aictab(cand.set = cand.models, modnames = Modnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)
      
# R2 values for candidate models, report marginal R2 below      
r.squaredGLMM(cand.models[[1]]) # m = 0.71
r.squaredGLMM(cand.models[[2]])# m = 0.66
r.squaredGLMM(cand.models[[3]])# m = 0.22
r.squaredGLMM(cand.models[[4]])# m = 0.63
r.squaredGLMM(cand.models[[5]])# m = 0.02
r.squaredGLMM(cand.models[[6]])# m = 0.60
r.squaredGLMM(cand.models[[7]])# m = 0.11
r.squaredGLMM(cand.models[[8]])# m = 0.00

summary(cand.models[[1]]) # summary of best model

############### FIGURE 3 . Raw seed counts in log scale as a function of month ##################
# add jeff's code here

# plot on single annual scale coded by year
ggplot(sheet,aes(x=as.factor(month),y=log(TOTAL+1),color=as.factor(yr))) + geom_jitter(shape=19,alpha=0.3,size=2,width=0.2) + theme_classic(base_size=14) + ylab("log(number of seeds+1)") + theme(legend.title=element_blank()) + theme(legend.position = c(0.15, 0.8)) + xlab("Month") + theme(legend.text = element_text(colour="black", size=12))
ggsave("fig3.png",device="png",width=7)

##### FIGURE 4. Average number on individuals across months as a function of age and sex ################
avg_tot1<-sheet %>% group_by(newage,month,SNAME) %>% summarise(avg=mean(TOTAL)) 
avg_tot2<-avg_tot1 %>% group_by(newage,month) %>% summarise(average=mean(avg),stdev=sd(avg),stderr=stdev/sqrt(length(avg)),max=max(avg)) 
avg_tot2<-as.data.frame(avg_tot2)
 
avg_tot2<-avg_tot2 %>% filter(month %in% c(6,7,8,9,10,11,12)) # only plot high seed months

ggplot(avg_tot2,aes(x=as.factor(month),y=average,fill=newage)) + geom_bar(stat="identity",color="black",size=0.4,position=position_dodge()) + geom_errorbar(aes(ymin = average, ymax = average + stderr),position=position_dodge())  + theme_classic(base_size=18) + ylab("Mean number of seeds per individual") + theme(legend.position = c(0.15, 0.75)) + labs(fill='Gelada body size') + xlab("") + theme(legend.title=element_text(size=11,face="bold")) + theme(legend.text=element_text(size=10))+
theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+  scale_fill_manual(values=cbPalette) 

ggsave("fig4_new.png",device="png",width=8)

####################################################################################################
### For body region figures, need to make new data frame t that transforms 'sheet' from wide to long ###
####################################################################################################
t<-sheet %>% gather(Location,Count,Face:DistalTail) # 
t<-as.data.frame(t)
head(t)

################################################################################
# New dataframe t1 builds on dataframe t and preserves location and takes the average count for individuals who were resampled during a day.
t1<-t %>% group_by(SNAME,newdate,newage,Location,month) %>% summarize(tot1=mean(Count,na.rm=T))
t1<-as.data.frame(t1)
head(t1, 50)

################ FIGURE 5. Seeds across different body regions ############################

t4<-t1 %>% group_by(newage,Location,month) %>% summarize(tot3=mean(tot1,na.rm=T)) # now summarize by date. average across individuals and dates to arrive at average for body location and size class
t4<-as.data.frame(t4)
head(t4, 50)

t4<-t4 %>% filter(month %in% c(8,9,10)) # Keep it to just the months when seeds were around.

t4$newage <- factor(t4$newage, levels = c("Very small", "Small", "Medium","Large"))

t4$newloc<-0

t4$newloc[which(t4$Location=="LowerLeg")]<-"Hindlimb"
t4$newloc[which(t4$Location=="ThighKnee")]<-"Hindlimb"
t4$newloc[which(t4$Location=="Hindquarter")]<-"Hindlimb"
t4$newloc[which(t4$Location=="HindLeg")]<-"Hindlimb"
t4$newloc[which(t4$Location=="ProximalTail")]<-"Hindlimb"
t4$newloc[which(t4$Location=="MidTail")]<-"Hindlimb"
t4$newloc[which(t4$Location=="DistalTail")]<-"Hindlimb"
t4$newloc[which(t4$Location=="Hip")]<-"Hindlimb"
t4$newloc[which(t4$Location=="AnkleFoot")]<-"Hindlimb"
t4$newloc[which(t4$Location=="WristHand")]<-"Forelimb"
t4$newloc[which(t4$Location=="UpperArmElbow")]<-"Forelimb"
t4$newloc[which(t4$Location=="Forearm")]<-"Forelimb"
t4$newloc[which(t4$Location=="LowerBack")]<-"Mid-upper body"
t4$newloc[which(t4$Location=="Shoulder")]<-"Mid-upper body"
t4$newloc[which(t4$Location=="UpperBack")]<-"Mid-upper body"
t4$newloc[which(t4$Location=="Ventrum")]<-"Mid-upper body"
t4$newloc[which(t4$Location=="Flank")]<-"Mid-upper body"
t4$newloc[which(t4$Location=="MidBack")]<-"Mid-upper body"
t4$newloc[which(t4$Location=="Head")]<-"Face&Head"
t4$newloc[which(t4$Location=="Face")]<-"Face&Head"
t4$newloc[which(t4$Location=="Cape")]<-"Cape"


t4$newloc <- factor(t4$newloc, levels = c("Hindlimb", "Mid-upper body", "Forelimb","Face&Head","Cape"))

t5 <- t4 %>% group_by(newage,newloc) %>% summarize(sums=mean(tot3,na.rm=T))

ggplot(t5,aes(x=newage,y=log(sums+1),fill=newloc)) + geom_col(color="black",position="dodge") + theme_classic(base_size=15) + ylab("Log(mean number of seeds + 1)")+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = c(0.2, 0.7))+xlab("")+theme(legend.title=element_blank())+  scale_fill_manual(values=cbPalette)+ labs(fill='Body location') + theme(legend.title=element_text(size=11,face="bold")) + theme(legend.text=element_text(size=10))

ggsave("fig5_oct18.png",device="png",width=8)

####### Figure S4 with more detailed body part locations #####
t3<-t1 %>% group_by(newage,Location,month) %>% summarize(tot3=mean(tot1,na.rm=T)) # now summarize by date. average across individuals and dates to arrive at average for body location and size class
t3<-as.data.frame(t3)
t3<-t3 %>% filter(month %in% c(6,7,8,9,10,11,12)) # Keep it to just the months when seeds were around.
t3$newage <- factor(t3$newage, levels = c("Very small", "Small", "Medium","Large"))
t3<-t3 %>% group_by(newage,Location) %>% summarize(tot3new=mean(tot3,na.rm=T),stdvnew=sd(tot3),stdernew=stdvnew/length(tot3))

###### PART THREE. FIGURE 6. Seed flux ###########

Rate<-read_excel("~/Desktop/seed paper/Triana/20210510/Agrocharis_Rate_20210510.xlsx")

#Convert count variable to numeric from factor
Rate$COUNT<-as.numeric(Rate$COUNT)

#Consider only focals where there was a change
Change_Data<-subset(Rate,COUNT_CHANGE == "Y")

#Colours for lines
colours <- c("#999999", "#E69F00", "#56B4E9", "#CC79A7","#009E73")

#Figure 6
Figure_6<-ggplot(data=Change_Data, aes(x=EVENT_CODE, y=COUNT, colour=BODY_REGION, group=FOCAL))+
  geom_point(size=2)+geom_line(size=1.25)+
  facet_wrap(~SNAME, scales="free_y")+
  labs(x="Focal Interval Number", y="Number of seeds", colour="Body Region")+
  scale_colour_manual(values=colours, labels=c("Back","Forearm","Hindlimb","Cape", "Ventrum"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.x=element_text(face="bold", margin=margin(0.1,0,0.1,0, "cm")),
        axis.title=element_text(size=16), axis.text=element_text(face="bold"),
        legend.title=element_text(size=24), legend.text=element_text(size=16), 
        legend.background=element_rect(color="black",linetype="solid"),
        legend.position="bottom")
Figure_6
