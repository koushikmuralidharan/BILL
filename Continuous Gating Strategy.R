#loadlibraries ===============================================================================================================
install.packages(c("tidyverse", "magrittr", "stringr"))
library(tidyverse)
library(magrittr)
library(readr)
library(stringr)
library(data.table)

#start here, change path and change path on csv======================================================================================================================
allfiles <- dir(path = "D:/Bill.R/RawDataTERT/M/M7", pattern = "*.csv")
setwd("~/RawDataTERT/M/M7")
alldata <- NULL
for (file in allfiles) {temp <- read_csv(file)
temp$id <- sub(".csv", "", file)
alldata <- rbind(alldata, temp)
}

#imports data into alldata===========================================================================================================================================
alldata <- alldata %>%
  dplyr::mutate(id=stringr::str_replace(id, pattern = "_Amplitude", replacement = ""))

alldata= mutate(alldata, rowid=substr(id, 4,4)) 
alldata= mutate(alldata, wellid=substr(id, 4,6)) 
names(alldata)[names(alldata)=="Ch1 Amplitude"] <- "one"
names(alldata)[names(alldata)=="Ch2 Amplitude"] <- "two"
alldata = subset(alldata, select = -c(Cluster,id) )

#defines standard deviations with intervals of 0.05 =================================================================================================================
deviations1 <- seq(4, 8, by = 0.05)
deviations2 <- seq(2, 6, by = 0.05)

#function that creates gates, along with objects that contain gate vectors ==========================================================================================
Ch1gatecreator <- function(deviations) {
      Ch1gate <- sum(mean(alldata$one) + deviations*sd(alldata$one))
    return(Ch1gate)
}
Ch1gatelist <- vector(mode = "double", length = 81)
for (i in seq(deviations1)) {
  Ch1gatelist[[i]] <- Ch1gatecreator(deviations1[[i]])
}
Ch2gatecreator <- function(deviations) {
  Ch2gate <- sum(mean(alldata$two) + deviations*sd(alldata$two))
  return(Ch2gate)
}
Ch2gatelist <- vector(mode = "double", length = 81)
for (i in seq(deviations2)) {
  Ch2gatelist[[i]] <- Ch2gatecreator(deviations2[[i]])
}

#table with every possible combination of gate ----------------------------------------------------------------------------------------------------------------------
GateCombo <- data.frame(expand.grid(Ch1gatelist, Ch2gatelist))
names(GateCombo)[names(GateCombo)=="Var1"] <- "Ch1gate"
names(GateCombo)[names(GateCombo)=="Var2"] <- "Ch2gate"

#subsetting alldata into row-based variables ------------------------------------------------------------------------------------------------------------------------
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

MutantA <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantA[[i]] <- sum(mergedA$one >= Ch1gatelist[[i]])
}
WildtypeA <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeA[[i]] <- sum(mergedA$two >= Ch2gatelist[[i]])
}
MutantB <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantB[[i]] <- sum(mergedB$one >= Ch1gatelist[[i]])
}
WildtypeB <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeB[[i]] <- sum(mergedB$two >= Ch2gatelist[[i]])
}
MutantC <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantC[[i]] <- sum(mergedC$one >= Ch1gatelist[[i]])
}
WildtypeC <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeC[[i]] <- sum(mergedC$two >= Ch2gatelist[[i]])
}
MutantD <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantD[[i]] <- sum(mergedD$one >= Ch1gatelist[[i]])
}
WildtypeD <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeD[[i]] <- sum(mergedD$two >= Ch2gatelist[[i]])
}
MutantE <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantE[[i]] <- sum(mergedE$one >= Ch1gatelist[[i]])
}
WildtypeE <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeE[[i]] <- sum(mergedE$two >= Ch2gatelist[[i]])
}
MutantF <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantF[[i]] <- sum(mergedF$one >= Ch1gatelist[[i]])
}
WildtypeF <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeF[[i]] <- sum(mergedF$two >= Ch2gatelist[[i]])
}
MutantG <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantG[[i]] <- sum(mergedG$one >= Ch1gatelist[[i]])
}
WildtypeG <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeG[[i]] <- sum(mergedG$two >= Ch2gatelist[[i]])
}
MutantH <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantH[[i]] <- sum(mergedH$one >= Ch1gatelist[[i]])
}
WildtypeH <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeH[[i]] <- sum(mergedH$two >= Ch2gatelist[[i]])
}
MutantI <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantI[[i]] <- sum(mergedI$one >= Ch1gatelist[[i]])
}
WildtypeI <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeI[[i]] <- sum(mergedI$two >= Ch2gatelist[[i]])
}
MutantJ <- vector("double", length = 81)
for (i in seq(Ch1gatelist)) {
  MutantJ[[i]] <- sum(mergedJ$one >= Ch1gatelist[[i]])
}
WildtypeJ <- vector("double", length = 81)
for (i in seq(Ch2gatelist)) {
  WildtypeJ[[i]] <- sum(mergedJ$two >= Ch2gatelist[[i]])
}

GateUno <- pull(GateCombo, "Ch1gate")
GateDos <- pull(GateCombo, "Ch2gate")
MAFACombo <- data.frame(expand.grid(MutantA, WildtypeA))
MAFBCombo <- data.frame(expand.grid(MutantB, WildtypeB))
MAFCCombo <- data.frame(expand.grid(MutantC, WildtypeC))
MAFDCombo <- data.frame(expand.grid(MutantD, WildtypeD))
MAFECombo <- data.frame(expand.grid(MutantE, WildtypeE))
MAFFCombo <- data.frame(expand.grid(MutantF, WildtypeF))
MAFGCombo <- data.frame(expand.grid(MutantG, WildtypeG))
MAFHCombo <- data.frame(expand.grid(MutantH, WildtypeH))
MAFICombo <- data.frame(expand.grid(MutantI, WildtypeI))
MAFJCombo <- data.frame(expand.grid(MutantJ, WildtypeJ))

GateUno <- pull(GateCombo, "Ch1gate")
GateDos <- pull(GateCombo, "Ch2gate")
MutPosA <- pull(MAFACombo, "Var1")
WTPosA <- pull(MAFACombo, "Var2")
MutPosB <- pull(MAFBCombo, "Var1")
WTPosB <- pull(MAFBCombo, "Var2")
MutPosC <- pull(MAFCCombo, "Var1")
WTPosC <- pull(MAFCCombo, "Var2")
MutPosD <- pull(MAFDCombo, "Var1")
WTPosD <- pull(MAFDCombo, "Var2")
MutPosE <- pull(MAFECombo, "Var1")
WTPosE <- pull(MAFECombo, "Var2")
MutPosF <- pull(MAFFCombo, "Var1")
WTPosF <- pull(MAFFCombo, "Var2")
MutPosF <- pull(MAFFCombo, "Var1")
WTPosF <- pull(MAFFCombo, "Var2")
MutPosG <- pull(MAFGCombo, "Var1")
WTPosG <- pull(MAFGCombo, "Var2")
MutPosH <- pull(MAFHCombo, "Var1")
WTPosH <- pull(MAFHCombo, "Var2")
MutPosI <- pull(MAFICombo, "Var1")
WTPosI <- pull(MAFICombo, "Var2")
MutPosJ <- pull(MAFJCombo, "Var1")
WTPosJ <- pull(MAFJCombo, "Var2")


MAFA <- 100*MutPosA/WTPosA
MAFB <- 100*MutPosB/WTPosB
MAFC <- 100*MutPosC/WTPosC
MAFD <- 100*MutPosD/WTPosD
MAFE <- 100*MutPosE/WTPosE
MAFF <- 100*MutPosF/WTPosF
MAFG <- 100*MutPosG/WTPosG
MAFH <- 100*MutPosH/WTPosH
MAFI <- 100*MutPosI/WTPosI
MAFJ <- 100*MutPosJ/WTPosJ

#Readout <- data.frame(GateUno, GateDos, MAFA, MAFB, MAFC, MAFD, MAFE, MAFF, MAFG, MAFH, MAFI, MAFJ)
Readout <- data.frame(GateUno, GateDos, MAFA, MAFB, MAFC, MAFD, MAFE, MAFF, MAFG, MAFH)

setwd("D:/Bill.R/")
write.csv(Readout, file = "M7_cont.csv")