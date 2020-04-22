#loadlibraries
install.packages(c("tidyverse", "magrittr", "stringr"))
library(tidyverse)
library(magrittr)
library(readr)
library(stringr)
library(data.table)

#start here, change path and change path on csv
allfiles <- dir(path = "C:/Users/koush/Documents/RawDataTERT/M/M8", pattern = "*.csv")
setwd("~/RawDataTERT/M/M8")
alldata <- NULL
for (file in allfiles) {temp <- read_csv(file)
temp$id <- sub(".csv", "", file)
alldata <- rbind(alldata, temp)
}

alldata <- alldata %>%
  dplyr::mutate(id=stringr::str_replace(id, pattern = "_Amplitude", replacement = ""))

alldata= mutate(alldata, rowid=substr(id, 4,4)) 
alldata= mutate(alldata, wellid=substr(id, 4,6)) 
names(alldata)[names(alldata)=="Ch1 Amplitude"] <- "one"
names(alldata)[names(alldata)=="Ch2 Amplitude"] <- "two"
alldata = subset(alldata, select = -c(Cluster,id) )

#gating strategy adopted from reading algorithm page of Lettalli ddpCR package for the definition of the empty droplets.
#The gate that is chosen is per plate.

ch1gate4 <- sum(mean(alldata$one), 4*sd(alldata$one))
ch1gate4.5 <- sum(mean(alldata$one), 4.5*sd(alldata$one))
ch1gate5 <- sum(mean(alldata$one), 5*sd(alldata$one))
ch1gate5.5 <- sum(mean(alldata$one), 5.5*sd(alldata$one))
ch1gate6 <- sum(mean(alldata$one), 6*sd(alldata$one))
ch1gate6.5 <- sum(mean(alldata$one), 6.5*sd(alldata$one))
ch1gate7 <- sum(mean(alldata$one), 7*sd(alldata$one))
ch1gate7.5 <- sum(mean(alldata$one), 7.5*sd(alldata$one))
ch1gate8 <- sum(mean(alldata$one), 8*sd(alldata$one))

ch2gate2 <- sum(mean(alldata$two), 2*sd(alldata$two))
ch2gate2.5 <- sum(mean(alldata$two), 2.5*sd(alldata$two))
ch2gate3 <- sum(mean(alldata$two), 3*sd(alldata$two))
ch2gate3.5 <- sum(mean(alldata$two), 3.5*sd(alldata$two))
ch2gate4 <- sum(mean(alldata$two), 4*sd(alldata$two))
ch2gate4.5 <- sum(mean(alldata$two), 4.5*sd(alldata$two))
ch2gate5 <- sum(mean(alldata$two), 5*sd(alldata$two))
ch2gate5.5 <- sum(mean(alldata$two), 5.5*sd(alldata$two))
ch2gate6 <- sum(mean(alldata$two), 6*sd(alldata$two))

mergedA <- subset(alldata, rowid == "A")
mergedB <- subset(alldata, rowid == "B")
mergedC <- subset(alldata, rowid == "C")
mergedD <- subset(alldata, rowid == "D")
mergedE <- subset(alldata, rowid == "E")
mergedF <- subset(alldata, rowid == "F")
mergedG <- subset(alldata, rowid == "G")
mergedH <- subset(alldata, rowid == "H")
mergedI <- subset(alldata, rowid == "I")
mergedJ <- subset(alldata, rowid == "J")

MAFA <- 100*sum(mergedA$one > ch1gate6)/sum(mergedA$two > ch2gate3)
MAFB <- 100*sum(mergedB$one > ch1gate6)/sum(mergedB$two > ch2gate3)
MAFC <- 100*sum(mergedC$one > ch1gate6)/sum(mergedC$two > ch2gate3)
MAFD <- 100*sum(mergedD$one > ch1gate6)/sum(mergedD$two > ch2gate3)
MAFE <- 100*sum(mergedE$one > ch1gate6)/sum(mergedE$two > ch2gate3)
MAFF <- 100*sum(mergedF$one > ch1gate6)/sum(mergedF$two > ch2gate3)
MAFG <- 100*sum(mergedG$one > ch1gate6)/sum(mergedG$two > ch2gate3)
MAFH <- 100*sum(mergedH$one > ch1gate6)/sum(mergedH$two > ch2gate3)
MAFI <- 100*sum(mergedI$one > ch1gate6)/sum(mergedI$two > ch2gate3)
MAFJ <- 100*sum(mergedJ$one > ch1gate6)/sum(mergedJ$two > ch2gate3)

MAFA2 <- 100*sum(mergedA$one > ch1gate6.5)/sum(mergedA$two > ch2gate3)
MAFB2 <- 100*sum(mergedB$one > ch1gate6.5)/sum(mergedB$two > ch2gate3)
MAFC2 <- 100*sum(mergedC$one > ch1gate6.5)/sum(mergedC$two > ch2gate3)
MAFD2 <- 100*sum(mergedD$one > ch1gate6.5)/sum(mergedD$two > ch2gate3)
MAFE2 <- 100*sum(mergedE$one > ch1gate6.5)/sum(mergedE$two > ch2gate3)
MAFF2 <- 100*sum(mergedF$one > ch1gate6.5)/sum(mergedF$two > ch2gate3)
MAFG2 <- 100*sum(mergedG$one > ch1gate6.5)/sum(mergedG$two > ch2gate3)
MAFH2 <- 100*sum(mergedH$one > ch1gate6.5)/sum(mergedH$two > ch2gate3)
MAFI2 <- 100*sum(mergedI$one > ch1gate6.5)/sum(mergedI$two > ch2gate3)
MAFJ2 <- 100*sum(mergedJ$one > ch1gate6.5)/sum(mergedJ$two > ch2gate3)

MAFA3 <- 100*sum(mergedA$one > ch1gate6)/sum(mergedA$two > ch2gate3.5)
MAFB3 <- 100*sum(mergedB$one > ch1gate6)/sum(mergedB$two > ch2gate3.5)
MAFC3 <- 100*sum(mergedC$one > ch1gate6)/sum(mergedC$two > ch2gate3.5)
MAFD3 <- 100*sum(mergedD$one > ch1gate6)/sum(mergedD$two > ch2gate3.5)
MAFE3 <- 100*sum(mergedE$one > ch1gate6)/sum(mergedE$two > ch2gate3.5)
MAFF3 <- 100*sum(mergedF$one > ch1gate6)/sum(mergedF$two > ch2gate3.5)
MAFG3 <- 100*sum(mergedG$one > ch1gate6)/sum(mergedG$two > ch2gate3.5)
MAFH3 <- 100*sum(mergedH$one > ch1gate6)/sum(mergedH$two > ch2gate3.5)
MAFI3 <- 100*sum(mergedI$one > ch1gate6)/sum(mergedI$two > ch2gate3.5)
MAFJ3 <- 100*sum(mergedJ$one > ch1gate6)/sum(mergedJ$two > ch2gate3.5)

MAFA4 <- 100*sum(mergedA$one > ch1gate6.5)/sum(mergedA$two > ch2gate3.5)
MAFB4 <- 100*sum(mergedB$one > ch1gate6.5)/sum(mergedB$two > ch2gate3.5)
MAFC4 <- 100*sum(mergedC$one > ch1gate6.5)/sum(mergedC$two > ch2gate3.5)
MAFD4 <- 100*sum(mergedD$one > ch1gate6.5)/sum(mergedD$two > ch2gate3.5)
MAFE4 <- 100*sum(mergedE$one > ch1gate6.5)/sum(mergedE$two > ch2gate3.5)
MAFF4 <- 100*sum(mergedF$one > ch1gate6.5)/sum(mergedF$two > ch2gate3.5)
MAFG4 <- 100*sum(mergedG$one > ch1gate6.5)/sum(mergedG$two > ch2gate3.5)
MAFH4 <- 100*sum(mergedH$one > ch1gate6.5)/sum(mergedH$two > ch2gate3.5)
MAFI4 <- 100*sum(mergedI$one > ch1gate6.5)/sum(mergedI$two > ch2gate3.5)
MAFJ4 <- 100*sum(mergedJ$one > ch1gate6.5)/sum(mergedJ$two > ch2gate3.5)

MAFA5 <- 100*sum(mergedA$one > ch1gate6)/sum(mergedA$two > ch2gate4)
MAFB5 <- 100*sum(mergedB$one > ch1gate6)/sum(mergedB$two > ch2gate4)
MAFC5 <- 100*sum(mergedC$one > ch1gate6)/sum(mergedC$two > ch2gate4)
MAFD5 <- 100*sum(mergedD$one > ch1gate6)/sum(mergedD$two > ch2gate4)
MAFE5 <- 100*sum(mergedE$one > ch1gate6)/sum(mergedE$two > ch2gate4)
MAFF5 <- 100*sum(mergedF$one > ch1gate6)/sum(mergedF$two > ch2gate4)
MAFG5 <- 100*sum(mergedG$one > ch1gate6)/sum(mergedG$two > ch2gate4)
MAFH5 <- 100*sum(mergedH$one > ch1gate6)/sum(mergedH$two > ch2gate4)
MAFI5 <- 100*sum(mergedI$one > ch1gate6)/sum(mergedI$two > ch2gate4)
MAFJ5 <- 100*sum(mergedJ$one > ch1gate6)/sum(mergedJ$two > ch2gate4)

MAFA6 <- 100*sum(mergedA$one > ch1gate6.5)/sum(mergedA$two > ch2gate4)
MAFB6 <- 100*sum(mergedB$one > ch1gate6.5)/sum(mergedB$two > ch2gate4)
MAFC6 <- 100*sum(mergedC$one > ch1gate6.5)/sum(mergedC$two > ch2gate4)
MAFD6 <- 100*sum(mergedD$one > ch1gate6.5)/sum(mergedD$two > ch2gate4)
MAFE6 <- 100*sum(mergedE$one > ch1gate6.5)/sum(mergedE$two > ch2gate4)
MAFF6 <- 100*sum(mergedF$one > ch1gate6.5)/sum(mergedF$two > ch2gate4)
MAFG6 <- 100*sum(mergedG$one > ch1gate6.5)/sum(mergedG$two > ch2gate4)
MAFH6 <- 100*sum(mergedH$one > ch1gate6.5)/sum(mergedH$two > ch2gate4)
MAFI6 <- 100*sum(mergedI$one > ch1gate6.5)/sum(mergedI$two > ch2gate4)
MAFJ6 <- 100*sum(mergedJ$one > ch1gate6.5)/sum(mergedJ$two > ch2gate4)

#MAFreadout <- c(MAFA, MAFB, MAFC, MAFD, MAFE, MAFF, MAFG)
MAFreadout <- c(MAFA, MAFB, MAFC, MAFD, MAFE, MAFF, MAFG, MAFH, MAFI, MAFJ)
MAFreadout2 <- c(MAFA2, MAFB2, MAFC2, MAFD2, MAFE2, MAFF2, MAFG2, MAFH2, MAFI2, MAFJ2)
MAFreadout3 <- c(MAFA3, MAFB3, MAFC3, MAFD3, MAFE3, MAFF3, MAFG3, MAFH3, MAFI3, MAFJ3)
MAFreadout4 <- c(MAFA4, MAFB4, MAFC4, MAFD4, MAFE4, MAFF4, MAFG4, MAFH4, MAFI4, MAFJ4)
MAFreadout5 <- c(MAFA5, MAFB5, MAFC5, MAFD5, MAFE5, MAFF5, MAFG5, MAFH5, MAFI5, MAFJ5)
MAFreadout6 <- c(MAFA6, MAFB6, MAFC6, MAFD6, MAFE6, MAFF6, MAFG6, MAFH6, MAFI6, MAFJ6)


AlleleFreq <- round(MAFreadout, digits = 2)
AlleleFreq2 <- round(MAFreadout2, digits = 2)
AlleleFreq3 <- round(MAFreadout3, digits = 2)
AlleleFreq4 <- round(MAFreadout4, digits = 2)
AlleleFreq5 <- round(MAFreadout5, digits = 2)
AlleleFreq6 <- round(MAFreadout6, digits = 2)

CallID <- ifelse(AlleleFreq >= 0.50, "Mut", "WT")
CallID2 <- ifelse(AlleleFreq2 >= 0.50, "Mut", "WT")
CallID3 <- ifelse(AlleleFreq3 >= 0.50, "Mut", "WT")
CallID4 <- ifelse(AlleleFreq4 >= 0.50, "Mut", "WT")
CallID5 <- ifelse(AlleleFreq5 >= 0.50, "Mut", "WT")
CallID6 <- ifelse(AlleleFreq6 >= 0.50, "Mut", "WT")

PlateRow <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#PlateRow <- c("A", "B", "C", "D", "E", "F", "G")

chartedMAF <- data.frame(PlateRow, AlleleFreq, CallID, AlleleFreq2, CallID2, AlleleFreq3, CallID3, AlleleFreq4, CallID4, AlleleFreq5, CallID5, AlleleFreq6, CallID6)

chartedMAF

write_csv(chartedMAF, path = "M8_HeatMapReport.csv")