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
setwd("D:/Bill.R/Bootstrapping")
allvalues <- read_csv("R.csv", )
mutantvalues <- select(allvalues, contains("Mut"))
wildtypevalues <- select(allvalues, contains("WT"))


#loads data in workable fashion ----------------------------------------------------------------------------------------------------
flippedvalues <- data.frame(t(allvalues))
workablevalues <- rownames_to_column(flippedvalues, var = "SNAPSHOT")
Cutpoint <- workablevalues[-c(1, 2, 3),]
Status <- ifelse(grepl("M", Cutpoint$SNAPSHOT), 1, 0)
CutpointFormat <- data.frame(add_column(Cutpoint, Status))


#want to convert the following into a function ======================================================================================
data <- data.frame(CutpointFormat$Status, CutpointFormat$X1)
names(data)[1] <- "Status"
names(data)[names(data)=="CutpointFormat.X1"] <- "data"
MinValueSp <- 0.9

cutpoint <- optimal.cutpoints(X = "data", status = "Status", tag.healthy = 0, methods = c("MaxSp", "MinValueSp", "ROC01", "MaxSpSe"), data = data, pop.prev = NULL, categorical.cov = NULL, control = control.cutpoints(), ci.fit = TRUE)
summary(cutpoint)

#want these to be vector outputs for X1-X6562 =======================================================================================
MaxSP_Sensitivity <- cutpoint[["MaxSp"]][["Global"]][["optimal.cutoff"]][["Se"]][[1]]
MaxSP_Specificity <- cutpoint[["MaxSp"]][["Global"]][["optimal.cutoff"]][["Sp"]][[1]]
SP90_Sensitivity <- cutpoint[["MinValueSp"]][[1]][["optimal.cutoff"]][["Se"]][[1]]
SP90_Specificity <- cutpoint[["MinValueSp"]][[1]][["optimal.cutoff"]][["Sp"]][[1]]
ROCO1_Sensitivity <- cutpoint[["ROC01"]][[1]][["optimal.cutoff"]][["Se"]][[1]]
ROCO1_Specificity <- cutpoint[["ROC01"]][[1]][["optimal.cutoff"]][["Sp"]][[1]]
MaxedSum_Sensitivity <- cutpoint[["MaxSpSe"]][[1]][["optimal.cutoff"]][["Se"]][[1]]
MaxedSum_Specificity <- cutpoint[["MaxSpSe"]][[1]][["optimal.cutoff"]][["Sp"]][[1]]  
  
  
ReportCutpoint <- data.frame(MaxSP_Specificity, MaxSP_Sensitivity, SP90_Specificity, SP90_Sensitivity, ROCO1_Specificity, ROCO1_Sensitivity, MaxedSum_Sensitivity, MaxedSum_Specificity)
names(ReportCutpoint) [names(ReportCutpoint)=="Value"] <- "MAXSP Specificity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.1"] <- "MAXSP Sensitivity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.2"] <- "SP90 Specificity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.3"] <- "SP90 Sensitivity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.4"] <- "ROCO1 Specificity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.5"] <- "ROCO1 Sensitivity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.6"] <- "MaxedSum Specificity"
names(ReportCutpoint) [names(ReportCutpoint)=="Value.7"] <- "Maxed Sensitivity"


