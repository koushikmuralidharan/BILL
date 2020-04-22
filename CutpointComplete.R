#load libraries ====================================================================================================================
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

#start here, change path and change path on csv -----------------------------------------------------------------------------------
setwd("D:/GatingTests")
allvalues <- read_csv("Training.csv")

mutantvalues <- select(allvalues, contains("Mut"))
wildtypevalues <- select(allvalues, contains("WT"))


#loads data in workable fashion ----------------------------------------------------------------------------------------------------
flippedvalues <- data.frame(t(allvalues))
workablevalues <- rownames_to_column(flippedvalues, var = "SNAPSHOT")
Cutpoint <- workablevalues[-c(1, 2, 3),]
Status <- ifelse(grepl("M", Cutpoint$SNAPSHOT), 1, 0)
CutpointFormat <- data.frame(add_column(Cutpoint, Status))


#putstogether cutoff calculations ======================================================================================
MaxSP_Sensitivity <- vector("double", length = 6561)
MaxSP_Specificity <- vector("double", length = 6561)
MaxSP_Threshold <- vector("double", length = 6561)
SP90_Sensitivity <- vector("double", length = 6561)
SP90_Specificity <- vector("double", length = 6561)
SP90_Threshold <- vector("double", length = 6561)
ROCO1_Sensitivity <- vector("double", length = 6561)
ROCO1_Specificity <- vector("double", length = 6561)
ROCO1_Threshold <- vector("double", length = 6561)
MaxedSum_Sensitivity <- vector("double", length = 6561)
MaxedSum_Specificity <- vector("double", length = 6561)
MaxedSum_Threshold <- vector("double", length = 6561)

for (number in seq_along(1:6561)) {
  data <- data.frame(CutpointFormat$Status, CutpointFormat[[1+number]]) 
  names(data)[1] <- "Status"
  names(data)[2] <- "data"
  MinValueSp <- 0.9
  cutpoint <- optimal.cutpoints(X = "data", status = "Status", tag.healthy = 0, methods = c("MaxSp", "MinValueSp", "ROC01", "MaxSpSe"), data = data, pop.prev = NULL, categorical.cov = NULL, control = control.cutpoints(), ci.fit = TRUE)
  MaxSP_Sensitivity[[number]] <- cutpoint[["MaxSp"]][["Global"]][["optimal.cutoff"]][["Se"]][[1]]
  MaxSP_Specificity[[number]] <- cutpoint[["MaxSp"]][["Global"]][["optimal.cutoff"]][["Sp"]][[1]]
  SP90_Sensitivity[[number]] <- cutpoint[["MinValueSp"]][[1]][["optimal.cutoff"]][["Se"]][[1]]
  SP90_Specificity[[number]] <- cutpoint[["MinValueSp"]][[1]][["optimal.cutoff"]][["Sp"]][[1]]
  ROCO1_Sensitivity[[number]] <- cutpoint[["ROC01"]][[1]][["optimal.cutoff"]][["Se"]][[1]]
  ROCO1_Specificity[[number]] <- cutpoint[["ROC01"]][[1]][["optimal.cutoff"]][["Sp"]][[1]]
  MaxedSum_Sensitivity[[number]] <- cutpoint[["MaxSpSe"]][[1]][["optimal.cutoff"]][["Se"]][[1]]
  MaxedSum_Specificity[[number]] <- cutpoint[["MaxSpSe"]][[1]][["optimal.cutoff"]][["Sp"]][[1]]
}


ReportCutpoint <- data.frame(MaxSP_Specificity, MaxSP_Sensitivity, SP90_Specificity, SP90_Sensitivity, ROCO1_Specificity, ROCO1_Sensitivity, MaxedSum_Specificity, MaxedSum_Sensitivity)
names(ReportCutpoint) [names(ReportCutpoint)=="Value"] <- "MAXSP Specificity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.1"] <- "MAXSP Sensitivity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.2"] <- "SP90 Specificity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.3"] <- "SP90 Sensitivity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.4"] <- "ROCO1 Specificity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.5"] <- "ROCO1 Sensitivity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.6"] <- "MaxedSum Specificity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.7"] <- "MaxedSum Sensitivity"

ThresholdReport <- ReportCutpoint*100
deviations1 <- seq(4, 8, by = 0.05)
deviations2 <- seq(2, 6, by = 0.05)
gatenumbers <- seq(1:6561)
GateCombo <- data.frame(expand.grid(deviations1, deviations2))
GateListCombo <- data.frame(GateCombo, gatenumbers)

BillBoardThreshold <- data.frame(MaxSP_Threshold, SP90_Threshold, ROCO1_Threshold, MaxedSum_Threshold)

BillBoard <- data.frame(GateListCombo, ThresholdReport)

View(BillBoard)
View(BillBoardThreshold)

setwd("D:/Bill.R/")
write.csv(BillBoard, file = "Sensitivities_Spec_Graphable.csv")

