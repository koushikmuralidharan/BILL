#loadlibraries ===============================================================================================================================================
install.packages("OptimalCutpoints")
library(optimalcutpoints)
install.packages("pROC")
library(pROC)
install.packages("ROCR")
library(ROCR)
install.packages("cutpointr")
library(cutpointr)
install.packages(c("tidyverse", "magrittr", "stringr"))
library(tidyverse)
library(magrittr)
library(readr)
library(stringr)
library(data.table)
install.packages("ggsci")
library(ggsci)

#sets up table for analysis ------------------------------------------------------------------------------------------------------
setwd("D:/GatingTests")
allvalues <- read_csv("Combined.csv", )
mutantvalues <- select(allvalues, contains("Mut"))
wildtypevalues <- select(allvalues, contains("WT"))

flippedvalues <- data.frame(t(allvalues))
workablevalues <- rownames_to_column(flippedvalues, var = "SNAPSHOT")
Cutpoint <- workablevalues[-c(1, 2, 3),]
Status <- ifelse(grepl("M", Cutpoint$SNAPSHOT), 1, 0)
CutpointFormat <- data.frame(add_column(Cutpoint, Status))

#sets up Mutant and WT tables as separate dataframes for analysis --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MutantSamples <- data.frame(subset(CutpointFormat, CutpointFormat$Status == 1))
WTSamples <- data.frame(subset(CutpointFormat, CutpointFormat$Status == 0))
TP <- vector(mode = "double", length = 19683)
FN <- vector(mode = "double", length = 19683)
FP <- vector(mode = "double", length = 19683)
TN <- vector(mode = "double", length = 19683)
Gate <- vector(mode = "double", length = 19683)
Cutoffvalue <- vector(mode = "double", length = 19683)

#for loop that dictates that for every gate in selectedgates generated previously, calculates contingency table values---------------------------------------------------------------------------------------------------------------------------------------
for (number in seq_along(1:19683)) {
  criterion <- data.frame(SelectedGatesSummary[number,])
  Gatenumber <- criterion$Gatenumber
  SNAPMut <- as.vector(select(MutantSamples, num_range("X", Gatenumber)))
  SNAPWT <- as.vector(select(WTSamples, num_range("X", Gatenumber)))
  TP[[number]] <- sum(SNAPMut >= criterion$Cutoff)
  FN[[number]] <- sum(SNAPMut < criterion$Cutoff)
  FP[[number]] <- sum(SNAPWT >= criterion$Cutoff)
  TN[[number]] <- sum(SNAPWT < criterion$Cutoff)
  Gate[[number]] <- criterion$Gatenumber
  Cutoffvalue[[number]] <- criterion$Cutoff
}

RawTestingTable <- data.frame(TP, FN, FP, TN, Gate)
Sensitivity <- (100*TP/(TP+FN))
Specificity <- (100*TN/(TN+FP))
AllGatesTest3 <- data.frame(Gate, SelectedGatesSummary$Gatenumber, Cutoffvalue, TP, FN, FP, TN, Sensitivity, Specificity, SelectedGatesSummary$Ch1gate, SelectedGatesSummary$Ch2gate)

View(AllGatesTest1)
View(AllGatesTest2)

Emptydataset <- NULL

#Attempt to Create a Function to Filter Gates in a query format ===========================================================================================
ReporterOld <- function(value, thres){
  A <- data.frame(AllGatesTest1[value,], "1")
  B <- data.frame(AllGatesTest2[value,], "2")
  C <- data.frame(AllGatesTest3[value,], "3")
  OutputReportPre <- bind_rows(A, B)
  OutputReport <- bind_rows(OutputReportPre, C)
  Emptydataset <- bind_rows(OutputReport)
}

Reporter <- function(value, thres){
  A <- data.frame(filter(AllGatesTest1, Gate == value & Cutoffvalue == thres), "1")
  B <- data.frame(filter(AllGatesTest2, Gate == value & Cutoffvalue == thres), "2")
  C <- data.frame(filter(AllGatesTest3, Gate == value & Cutoffvalue == thres), "3")
  OutputReportPre <- bind_rows(A, B)
  OutputReport <- bind_rows(OutputReportPre, C)
}

#Filtering gates and results by LB criteria "Atleast 60/90" ==============================================================================================================================

shortlist <- filter(AllGatesTest3, AllGatesTest3$Sensitivity >= 60, AllGatesTest3$Specificity >= 90)

View(shortlist)
shortlisttraining <- subset(AllGatesTest1, Gate %in% shortlist$Gate)
shortlisttraining <- subset(shortlisttraining, Cutoffvalue %in% shortlist$Cutoffvalue)
shortlistvalidation <- subset(AllGatesTest2, Gate %in% shortlist$Gate)
shortlistvalidation <- subset(shortlistvalidation, Cutoffvalue %in% shortlist$Cutoffvalue)

#CodeforHeatmapGeneration ------------------------------------------------------------------------------------------------------------
df_total <- data.frame(seq(1:158)) 
for (number in seq(1:146)) {
  qualifiers <- data.frame(shortlist[number, ])
  GateChoice <- qualifiers$SelectedGatesSummary.Gatenumber
  CutoffGate <- qualifiers$Cutoffvalue
  MAFList <- CutpointFormat[[1+GateChoice]]
  NormalMAF <- MAFList - CutoffGate
  CallStatus <- data.frame((if_else(MAFList >= CutoffGate, 1, 0)))
  df_total <- cbind(df_total, NormalMAF)
}

names(df_total) <- c(seq(1:147))
HeatmapOutput <- subset(df_total, select = -1)
SampleName <- c(CutpointFormat$SNAPSHOT)
HeatmapReport <- cbind(SampleName, HeatmapOutput)
View(HeatmapReport)

write.csv(HeatmapReport, file = "NormalMAFHeatmapReportFinal.csv")

#CreatingFinalReport of Results --------------------------------------------------------------------------------------------------------
colnames(shortlist)[8] = "Sensitivity_Combined"
colnames(shortlist)[9] = "Specificity_Combined"
colnames(shortlisttraining)[8] = "Sensitivity_Training"
colnames(shortlisttraining)[9] = "Specificity_Training"
colnames(shortlistvalidation)[8] = "Sensitivity_Validation"
colnames(shortlistvalidation)[9] = "Specificity_Validation"

common_column_names <- c("Gate", "Cutoffvalue")
coleworld <- merge(shortlist, shortlisttraining, by=common_column_names, all.x = TRUE)
kdotwhereyouat <- merge(coleworld, shortlistvalidation, by=common_column_names, all.x = TRUE)
yeezus <- data.frame(kdotwhereyouat$Gate, kdotwhereyouat$SelectedGatesSummary.Ch1gate, kdotwhereyouat$SelectedGatesSummary.Ch2gate.y, kdotwhereyouat$Cutoffvalue, kdotwhereyouat$Sensitivity_Combined, kdotwhereyouat$Specificity_Combined, kdotwhereyouat$Sensitivity_Training, kdotwhereyouat$Specificity_Training, kdotwhereyouat$Sensitivity_Validation, kdotwhereyouat$Specificity_Validation)
final_report <- setnames(yeezus, c("Gate","Ch1gate","Ch2gate","Cutoffvalue", "Comb_Sensitivity", "Comb_Specificity", "Train_Sensitivity", "Train_Specificity", "Val_Sensitivity", "Cal_Specificity"))
SensitivitySum <- final_report$Comb_Sensitivity + final_report$Train_Sensitivity + final_report$Val_Sensitivity
FinalReport <- data.frame(final_report, SensitivitySum)

#putting out the final report sensitivities.=====================================================================================================================================================================================================================================================================================================================================================================
write.csv(FinalReport, file = "FinalReportSensitivities.csv")
View(FinalReport)

m <- ggplot(FinalReport, aes(FinalReport$Ch1gate, FinalReport$Ch2gate))
m + geom_bar()
 