#NOTE - this is just a temporary working script. All finalized data analyses and plotting are in 'disp.R'
library(ggplot2)
library(dplyr)
library(lubridate)
library(gridExtra)

#Date	RAW	DOY	NDVI	Year
#Note this raw NDVI file is directly downloaded, unformatted from the 
#Earth engine script called 'GuassaNDVI_Agrocharis_Paper'
Guassa_raw_NDVI<-read.csv("data/NDVI_Guassacore_00-21_QA1.csv", header=T)

df2<-Guassa_raw_NDVI %>%
  # recode empty strings "" by NAs
  na_if("") %>%
  # remove NAs
  na.omit

df3 <- df2 %>% rename(full_date = system.time_start, rawNDVI = undefined)
df3$full_date<-as.character(df3$full_date)
df3$formatted_full_date<-as.Date(df3$full_date, format="%b %d, %Y")

df3$MonthDay<-format(df3$formatted_full_date, format = "%m-%d")
df3$DOY<-yday(df3$formatted_full_date)
df3$YR<-year(df3$formatted_full_date)
df3$NDVI <- as.numeric(gsub(",","",df3$rawNDVI))/10000
df3_SP<-df3 %>% filter(YR >= 2017) %>% filter (YR <= 2019)

# NDVI through year, all years
ggplot(df3, aes(x=MonthDay, y=NDVI)) +
  geom_smooth(method = "loess", size = 1.5) +
  geom_point(data=df3, aes(x=MonthDay, y=NDVI), color='black', alpha=0.3)+
  #geom_point(data=GuassaNDVIc18, aes(x=MonthDay, y=NDVI), color='red', cex =5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# NDVI through year, all years
NDVI_22yrs<-ggplot(df3, aes(x=DOY, y=NDVI)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_rect(data=NULL,aes(xmin=305,xmax=365,ymin=-Inf,ymax=Inf),
            fill="lightyellow")+
  geom_rect(data=NULL,aes(xmin=0,xmax=151,ymin=-Inf,ymax=Inf),
            fill="lightyellow")+
  geom_rect(data=NULL,aes(xmin=151,xmax=305,ymin=-Inf,ymax=Inf),
            fill="lightgreen")+
  geom_smooth(method = "loess", size = 1.5, col='darkgreen') +
  geom_point(data=df3, aes(x=DOY, y=NDVI), color='black', alpha=0.3)+
  #geom_point(data=GuassaNDVIc18, aes(x=DOY, y=NDVI), color='red', cex =5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# NDVI through year, study years
NDVI_studyperiod<-ggplot(df3_2017_2019, aes(x=DOY, y=NDVI)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_rect(data=NULL,aes(xmin=305,xmax=365,ymin=-Inf,ymax=Inf),
            fill="lightyellow")+
  geom_rect(data=NULL,aes(xmin=0,xmax=151,ymin=-Inf,ymax=Inf),
            fill="lightyellow")+
  geom_rect(data=NULL,aes(xmin=151,xmax=305,ymin=-Inf,ymax=Inf),
            fill="lightgreen")+
  geom_smooth(method = "loess", size = 1.5) +
  geom_point(data=df3_2017_2019, aes(x=DOY, y=NDVI), color='black', alpha=0.3)+
  #geom_point(data=GuassaNDVIc18, aes(x=DOY, y=NDVI), color='red', cex =5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Ag_data<-read.csv('data/sheet_mods.csv', header=TRUE)

Ag_data$DOY<-yday(Ag_data$newdate)
Ag_data$LogTOTAL.1<-log(Ag_data$TOTAL+1)

LogSeeds_studyperiod<-ggplot(Ag_data, aes(x=DOY, y=LogTOTAL.1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_rect(data=NULL,aes(xmin=305,xmax=365,ymin=-Inf,ymax=Inf),
            fill="lightyellow")+
  geom_rect(data=NULL,aes(xmin=0,xmax=151,ymin=-Inf,ymax=Inf),
            fill="lightyellow")+
  geom_rect(data=NULL,aes(xmin=151,xmax=305,ymin=-Inf,ymax=Inf),
            fill="lightgreen")+
  #geom_smooth(method = "loess", size = 1.5, col = 'black') +
  geom_jitter(data=Ag_data, aes(x=DOY, y=LogTOTAL.1), color='black', alpha=0.3, width = 0.5)+
  #geom_point(data=GuassaNDVIc18, aes(x=DOY, y=NDVI), color='red', cex =5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




grid.arrange(NDVI_studyperiod, LogSeeds_studyperiod, ncol=1)
grid.arrange(NDVI_22yrs, LogSeeds_studyperiod, ncol=1)

gB <- ggplotGrob(LogSeeds_studyperiod)
gA <- ggplotGrob(NDVI_22yrs)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1)

df3_2018<-df3 %>% filter(YR == 2018)

NDVI_18.lo <- loess(NDVI ~ DOY, df3_2018)
predict(NDVI_18.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)
dailyNDVI_18<-predict(NDVI_18.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)
plot(dailyNDVI_18$fit)
points(df3_2018$DOY, df3_2018$NDVI)

NDVI_17_19.lo <- loess(NDVI ~ DOY, df3_2017_2019)
predict(NDVI_17_19.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)
dailyNDVI_17_19<-predict(NDVI_17_19.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)
plot(dailyNDVI_17_19$fit)
points(df3_2017_2019$DOY, df3_2017_2019$NDVI)
points(dailyNDVI$fit)

NDVI.lo <- loess(NDVI ~ DOY, df3)
predict(NDVI.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)
dailyNDVI<-predict(NDVI.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)
plot(dailyNDVI$fit)
points(df3$DOY, df3$NDVI)

Ag_data<-read.csv('data/sheet_mods.csv', header=TRUE)

Ag_data$DOY<-yday(Ag_data$newdate)

plot(Ag_data$DOY, log(Ag_data$TOTAL+1))

hist(log(Ag_data$TOTAL+1))

Ag.lo <- loess(log(TOTAL+1) ~ DOY, Ag_data)
predict(Ag.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)
dailyAg<-predict(Ag.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)
plot(dailyAg$fit)

abline(h=0)
points(log(Ag_data$TOTAL+1))


#Create dataframe with daily NDVI estimates and daily Agrocharis data

Ag_NDVI_df<-NULL
Ag_NDVI_df$DOY<-seq(1,365, 1)
Ag_NDVI_df$LogTotal.90ile<-NULL
for (i in 1:365){
  Ag_NDVI_df$LogTotal.90ile[Ag_NDVI_df$DOY==i] <-quantile(log(Ag_data$TOTAL[Ag_data$DOY==i]+1), 0.9)
}
plot(Ag_NDVI_df$LogTotal.90ile)

Ag_log.lo <- loess(LogTotal.90ile ~ DOY, Ag_NDVI_df)
predict(Ag_log.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)

Ag_NDVI_df$Ag_log.90ile.lo<-as.vector(predict(Ag_log.lo, data.frame(DOY = seq(1, 365, 1)), se = TRUE)$fit)

Ag_NDVI_df$NDVI<-as.vector(dailyNDVI$fit)

Ag_NDVI_df<-data.frame(Ag_NDVI_df)

plot(Ag_log.90ile.lo~ DOY, data=Ag_NDVI_df)
