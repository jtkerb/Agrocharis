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
library(gridExtra)
library(quantreg)
library(olsrr)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # grey

'%ni%' <- Negate('%in%')

### Read in agrocharis dataset

sheet<-read.csv("data/agrosheet.csv", , fileEncoding="UTF-8-BOM")
sheet<-as.data.frame(sheet)
sheet$SNAME<-as.factor(sheet$SNAME)
sheet<-sheet %>% mutate(newdate=dmy(DATE)) %>% mutate(month=month(newdate),yr=year(newdate),dy=day(newdate))# add column with more standard notation for date and add columns that pull out month, year and date# new variable to make new age categories based on size. Only unintuitive assignment is LJs being included as adult females.

sheet$newage<-0 # new column for re-labeling age-sex classes
sheet$newage[which(sheet$AGE=="SJ" | sheet$AGE=="MJ")]<-"Small"
sheet$newage[which(sheet$AGE=="BRI" | sheet$AGE=="BLI")]<-"Very small"
sheet$newage[which(sheet$AGE=="A" | sheet$AGE=="SAM"| sheet$AGE=="SA" & sheet$SEX=="M")]<-"Large"
sheet$newage[which(sheet$AGE=="A" & sheet$SEX=="F")]<-"Medium"
sheet$newage[which(sheet$AGE=="LJ" & sheet$SEX=="M")]<-"Medium"
sheet$newage[which(sheet$AGE=="LJ" & sheet$SEX=="F")]<-"Medium"
sheet$newage<-as.factor(sheet$newage)
sheet$newage <- factor(sheet$newage, order=TRUE, levels = c("Very small", "Small","Medium","Large")) # order factor for body size

sheet[is.na(sheet)] <- 0 #change NA's to zeros
sheet<-sheet %>% drop_na(newage)

#Read in and prepare the NDVI dataset
#Note this is the NDVI file as directly downloaded from the 
#Earth engine script called 'GuassaNDVI_Agrocharis_Paper'
df1<-read.csv("data/NDVI_Guassacore_00-21_QA1.csv", header=T)

df2<-df1 %>%
  na_if("") %>%   # re-code empty strings "" based on NAs
  na.omit # remove NAs

df3 <- df2 %>% rename(full_date = system.time_start, rawNDVI = undefined)

#preapre columns
df3$full_date<-as.character(df3$full_date)
df3$formatted_full_date<-as.Date(df3$full_date, format="%b %d, %Y")
df3$MonthDay<-format(df3$formatted_full_date, format = "%m-%d")
df3$Month<-format(df3$formatted_full_date, format = "%m")
df3$DOY<-yday(df3$formatted_full_date)
df3$YR<-year(df3$formatted_full_date)
df3$NDVI <- as.numeric(gsub(",","",df3$rawNDVI))/10000

#subset to study period years
df3_SP<-df3 %>% filter(YR >= 2017) %>% filter (YR <= 2019)


## A few basic attributes of the agrocharis dataset
length(sheet$TOTAL) # 5344 rows
reps<-sheet %>% group_by(SNAME,newdate) %>% summarize(number=n()) # show repeated sampling of individuals
table(reps$number) # rather minor. 5123 samples were once per day, 116 instances of twice per day, 4 instances of 3x per day. So, in 98% of cases observations were made on animals once per day
max(sheet$TOTAL) # 511 burrs on a single individual in a day
length(unique(sheet$SNAME)) #188 unique individuals throughout study, new2021 says 232 inds
sum(sheet$TOTAL) # 12552, 12563 burrs counted overall
#look at how many of each age sex class were sampled
agesex<-sheet %>% group_by(month,newage) %>% summarize(number=n(),avg=mean(TOTAL),maxim=max(TOTAL)) # show repeated sampling of individuals
head(agesex)
agesex<-as.data.frame(agesex)
agesex # not very balanced sampling, need to explain this or deal with analytically 

b<-agesex %>% group_by(newage) %>% summarize(ss=mean(number))
b<-as.data.frame(b)
b # avg sample sizes across months

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# Generalized linear mixed models (GLMMs) to examine determinants of seed counts
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

## First, need to examine issue of zero-inflation in the dataset
sum(sheet$TOTAL > 0) # 1539 nonzero
sum(sheet$TOTAL == 0) # 3805 are zero, so, 60% of dataset is zeros
plot(table(sheet$TOTAL)) # visual representation of zero-inflation

# new datasheet for models only
sheet_mods <- sheet %>% filter(month!=3) # remove months 3 and 4 for analysis because no seeds available in environment
sheet_mods <- sheet %>% filter(month!=4)
sheet_mods$month<-as.factor(sheet_mods$month) #turn month into factor for analysis
sheet_mods$yr<-as.factor(sheet_mods$yr)

write.csv(sheet_mods, 'data/sheet_mods.csv')

# List of candidate models
cand.models <- list()
cand.models[[1]] <- glmer(TOTAL ~ newage + yr + month + (1|SNAME:month), data = sheet_mods,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[2]] <- glmer(TOTAL ~ newage + month + (1|SNAME:month), data = sheet_mods,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[3]] <- glmer(TOTAL ~ newage + yr + (1|SNAME:month), data = sheet_mods,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[4]] <- glmer(TOTAL ~ yr + month + (1|SNAME:month), data = sheet_mods,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[5]] <- glmer(TOTAL ~ yr + (1|SNAME:month), data = sheet_mods,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[6]] <- glmer(TOTAL ~ month  + (1|SNAME:month), data = sheet_mods,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[7]] <- glmer(TOTAL ~ newage + (1|SNAME:month), data = sheet_mods,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
cand.models[[8]] <- glmer(TOTAL ~ (1|SNAME:month), data = sheet_mods,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

Modnames <- paste("mod", 1:length(cand.models), sep = " ")

##generate AICc table
aictab(cand.set = cand.models, modnames = Modnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)
      
# R2 values for candidate models, report marginal R2 below      
r.squaredGLMM(cand.models[[1]]) # m = 0.64
r.squaredGLMM(cand.models[[2]])# m = 0.62
r.squaredGLMM(cand.models[[3]])# m = 0.07
r.squaredGLMM(cand.models[[4]])# m = 0.58
r.squaredGLMM(cand.models[[5]])# m = 0.007
r.squaredGLMM(cand.models[[6]])# m = 0.57
r.squaredGLMM(cand.models[[7]])# m = 0.04
r.squaredGLMM(cand.models[[8]])# m = 0.00

summary(cand.models[[1]]) # summary of best model

#####################################################################################################
#####################################################################################################
# FIGURE 3 . Raw seed counts in log scale as a function of month ##########################################################################################################################################################################################################################

# remove year, plot on single annual scale
# NDVI through year, all years
df3$Month<-as.numeric(df3$Month)




all<-left_join(df3, sheet, by = c("Month" = "month"))

top<- all %>%                                    
  arrange(desc(TOTAL)) %>% 
  group_by(Month) %>%
  slice(1:500)

head(all)
(Fig_compare<-ggplot(all,aes(x=NDVI,y=log(TOTAL+1),color=as.factor(yr))) +
    geom_jitter(shape=19,alpha=0.02,size=2,width=0.2) + 
    theme_classic(base_size=14) + 
    ylab("log(number of seeds+1)") + 
    theme(legend.title=element_blank()) + 
    #theme(legend.position = c(0.15, 0.8)) + 
    xlab("NDVI") + 
    theme(legend.text = element_text(colour="black", size=12)) +
    geom_quantile(quantiles = 0.9, size=4))

(Fig_TOP500_compare<-ggplot(top,aes(x=NDVI,y=log(TOTAL+1))) +
    geom_jitter(shape=19,alpha=0.3,size=2,width=0.02) + 
    theme_classic(base_size=14) + 
    ylab("log(number of seeds+1)") + 
    theme(legend.title=element_blank()) + 
    #theme(legend.position = c(0.15, 0.8)) + 
    xlab("NDVI") + 
    theme(legend.text = element_text(colour="black", size=12)) +
    geom_smooth(method='lm'))
    #geom_quantile(quantiles = 0.5, size=4))

(Fig_TOP500_compare2<-ggplot(top,aes(y=NDVI,x=log(TOTAL+1))) +
    geom_jitter(shape=19,alpha=0.1,size=2,width=0.05) + 
    theme_classic(base_size=14) + 
    xlab("log(number of seeds+1)") + 
    theme(legend.title=element_blank()) + 
    #theme(legend.position = c(0.15, 0.8)) + 
    ylab("NDVI") + 
    theme(legend.text = element_text(colour="black", size=12)) +
    #geom_smooth(method='lm', aes(color="Linear Model"), se=FALSE, lwd=2) +
    #geom_smooth(method="lm", aes(color="Exp Model"), formula= (y ~ exp(-x)), se=FALSE, linetype = 1)+
    geom_smooth(method="loess", se=FALSE, lwd=2, col=2))
    #geom_smooth(method="loess",aes(color="Loess"), se=FALSE, lwd=2))
#geom_quantile(quantiles = 0.5, size=4))

(Fig_compare2<-ggplot(all,aes(x=NDVI,y=log(TOTAL+1))) +
    geom_jitter(shape=19,alpha=0.02,size=2,width=0.02) + 
    theme_classic(base_size=14) + 
    ylab("log(number of seeds+1)") + 
    theme(legend.title=element_blank()) + 
    #theme(legend.position = c(0.15, 0.8)) + 
    xlab("NDVI") + 
    theme(legend.text = element_text(colour="black", size=12)) +
    geom_quantile(quantiles = 0.9, size=4))

(Fig1a<-ggplot(df3, aes(x=as.numeric(Month), y=NDVI)) +
  theme_classic(base_size=14) + 
  geom_smooth(method = "loess", size = 1.5, col='darkgreen') +
  xlab("Month") +
  geom_point(data=df3, aes(x=Month, y=NDVI), color='black', alpha=0.1)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_point(data=df3_SP, aes(as.numeric(Month), y=NDVI,color=as.factor(YR)), alpha=0.3, size =2)+
  theme(legend.position='none'))


(Fig1b<-ggplot(sheet,aes(x=as.factor(month),y=log(TOTAL+1),color=as.factor(yr))) +
  geom_jitter(shape=19,alpha=0.3,size=2,width=0.2) + 
  theme_classic(base_size=14) + 
  ylab("log(number of seeds+1)") + 
  theme(legend.title=element_blank()) + 
  theme(legend.position = c(0.15, 0.8)) + 
  xlab("Month") + 
  theme(legend.text = element_text(colour="black", size=12)))

f1A <- ggplotGrob(Fig1a)
f1B <- ggplotGrob(Fig1b)
f1C <- ggplotGrob(Fig_TOP500_compare2)

maxWidth = grid::unit.pmax(f1A$widths[2:5], f1B$widths[2:5])
f1A$widths[2:5] <- as.list(maxWidth)
f1B$widths[2:5] <- as.list(maxWidth)
Fig1ab<-arrangeGrob(f1A, f1B,f1C, ncol=1)
ggsave(Fig1ab, filename ="figs/fig1ab_new2.png",device="png",width=6, height = 14)

####################################################################################################
####################################################################################################
### FIGURE 4. Average number on individuals across months as a function of age and sex ###############################################################################################
####################################################################################################
avg_tot1<-sheet %>% group_by(newage,month,SNAME) %>% summarise(avg=mean(TOTAL)) 
avg_tot2<-avg_tot1 %>% group_by(newage,month) %>% summarise(average=mean(avg),stdev=sd(avg),stderr=stdev/sqrt(length(avg)),max=max(avg)) 
avg_tot2<-as.data.frame(avg_tot2)
#avg_tot2<-avg_tot2  filter(month==)

avg_tot2<-avg_tot2 %>% filter(month %in% c(6,7,8,9,10,11,12)) # only plot high seed months

ggplot(avg_tot2,aes(x=as.factor(month),y=average,fill=newage)) + geom_bar(stat="identity",color="black",size=0.4,position=position_dodge()) + geom_errorbar(aes(ymin = average, ymax = average + stderr),position=position_dodge())  + theme_classic(base_size=18) + ylab("Average number of seeds") + theme(legend.position = c(0.15, 0.75)) + labs(fill='Gelada body size') + xlab("") + theme(legend.title=element_text(size=11,face="bold")) + theme(legend.text=element_text(size=10))+
theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+  scale_fill_manual(values=cbPalette) 

ggsave("figs/fig4new.png",device="png",width=8)

####################################################################################################
### make new data frame t that makes 'sheet' from wide to long ###
####################################################################################################
t<-sheet %>% gather(Location,Count,Face:DistalTail) # 
t<-as.data.frame(t)
head(t)

####################################################################################################
####################################################################################################
# New dataframe t1 builds on dataframe t and preserves location and takes the average count for individuals who were resampled during a day.
t1<-t %>% group_by(SNAME,newdate,newage,Location,month) %>% summarize(tot1=mean(Count,na.rm=T))
t1<-as.data.frame(t1)

############################################################################################################
############################################################################################################
################# Figure 5. Examine the location of burrs in different places ############################
############################################################################################################
############################################################################################################

t3<-t1 %>% group_by(newage,Location,month) %>% summarize(tot3=mean(tot1,na.rm=T)) # now summarize by date. average across individuals and dates to arrive at average for body location and size class
t3<-as.data.frame(t3)

t3<-t3 %>% filter(month %in% c(6,7,8,9,10,11,12)) # Keep it to just the months when seeds were around.

t3$newage <- factor(t3$newage, levels = c("Very small", "Small", "Medium","Large"))

t3<-t3 %>% group_by(newage,Location) %>% summarize(tot3new=mean(tot3,na.rm=T),stdvnew=sd(tot3),stdernew=stdvnew/length(tot3))

### Single plot like above but with all size categories in same plot ###
ggplot(t3,aes(x=reorder(Location,-tot3new),y=tot3new,fill=newage)) + geom_col(color="black") + theme_classic(base_size=15) + ylab("Average number of seeds")+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = c(0.85, 0.7))+xlab("")+theme(legend.title=element_blank())+  scale_fill_manual(values=cbPalette)+ labs(fill='Gelada body size') + theme(legend.title=element_text(size=11,face="bold")) + theme(legend.text=element_text(size=10))

ggsave("figs/fig5new.png",device="png",width=8)



