#loadlibraries ===============================================================================================================================================
install.packages("OptimalCutpoints")
library(OptimalCutpoints)
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

#plotting sensitivities and specificities as a function of changing gate ========================================================================
Spec100Plot <- ggplot(data = BillBoard, aes(x = Var1, y = MaxSP_Sensitivity, color = "Sens")) + scale_fill_continuous() + geom_line() + facet_wrap(BillBoard$Var2) + geom_line(data = BillBoard, aes(x = Var1, y = MaxSP_Specificity, color = "Spec"))
Specto90Plot <- ggplot(data = BillBoard, aes(x = Var1, y = SP90_Sensitivity, color = "Sens")) + scale_fill_continuous() + geom_line() + facet_wrap(BillBoard$Var2) + geom_line(data = BillBoard, aes(x = Var1, y = SP90_Specificity, color = "Spec"))
ROCO01Plot <- ggplot(data = BillBoard, aes(x = Var1, y = ROCO1_Sensitivity, color = "Sens")) + scale_fill_continuous() + geom_line() + facet_wrap(BillBoard$Var2) + geom_line(data = BillBoard, aes(x = Var1, y = ROCO1_Specificity, color = "Spec"))
MaxedSumSensitivity <- ggplot(data = BillBoard, aes(x = Var1, y = MaxedSum_Sensitivity, color = "Sens")) + scale_fill_continuous() + geom_line() + facet_wrap(BillBoard$Var2) + geom_line(data = BillBoard, aes(x = Var1, y = MaxedSum_Specificity, color = "Spec"))

#loop that calculates thresholds for 3/4 =============================================================================================================================================================================================================================================
Spec100Threshold <- vector("double", length = 6561)
Spec90Threshold <- vector("double", length = 6561)
ROCO1Threshold <- vector("double", length = 6561)
for (number in seq_along(1:6561)) {
  data <- data.frame(CutpointFormat$Status, CutpointFormat[[1+number]])
  names(data)[1] <- "Status"
  names(data)[2] <- "data"
  MinValueSp <- 0.9
  cutpoint <- optimal.cutpoints(X = "data", status = "Status", tag.healthy = 0, methods = c("MaxSp", "MinValueSp", "ROC01", "MaxSpSe"), data = data, pop.prev = NULL, categorical.cov = NULL, control = control.cutpoints(), ci.fit = TRUE)
  Spec100Threshold[[number]] <- summary(cutpoint)[["MaxSp"]][["Global"]][["optimal.cutoff"]][["cutoff"]]
  Spec90Threshold[[number]] <- summary(cutpoint)[["MinValueSp"]][[1]][["optimal.cutoff"]][["cutoff"]]
  ROCO1Threshold[[number]] <- summary(cutpoint)[["ROC01"]][[1]][["optimal.cutoff"]][["cutoff"]]
  }

#tables with results from BillBoard and from threshold calculating loop ==============================================================================================================================================================================================================
Spec100Frame <- data.frame(Spec100Threshold, BillBoard$gatenumbers, BillBoard$Var1, BillBoard$Var2, BillBoard$MaxSP_Sensitivity, BillBoard$MaxSP_Specificity)
SelectedGatesSpec100 <- subset(Spec100Frame, BillBoard.MaxSP_Sensitivity >= 0)
Spec90Frame <- data.frame(Spec90Threshold, BillBoard$gatenumbers, BillBoard$Var1, BillBoard$Var2, BillBoard$SP90_Sensitivity, BillBoard$SP90_Specificity)
SelectedGatesSpec90 <- subset(Spec90Frame, BillBoard.SP90_Specificity >= 0)
ROCO1Frame <- data.frame(ROCO1Threshold, BillBoard$gatenumbers, BillBoard$Var1, BillBoard$Var2, BillBoard$ROCO1_Sensitivity, BillBoard$ROCO1_Specificity)
SelectedGatesROCO1 <- subset(ROCO1Frame, BillBoard.ROCO1_Specificity >= 0)
MaxSumFrame <- data.frame(BillBoard$gatenumbers, BillBoard$Var1, BillBoard$Var2, BillBoard$MaxedSum_Sensitivity, BillBoard$MaxedSum_Specificity)
SelectedGatesMaxSum <- subset(MaxSumFrame, BillBoard.MaxedSum_Specificity >= 0)

names(SelectedGatesSpec100) <- c("Cutoff", "Gatenumber", "Ch1gate", "Ch2gate", "Sensitivity", "Specificity")
names(SelectedGatesSpec90) <- c("Cutoff", "Gatenumber", "Ch1gate", "Ch2gate", "Sensitivity", "Specificity")
names(SelectedGatesROCO1) <- c("Cutoff", "Gatenumber", "Ch1gate", "Ch2gate", "Sensitivity", "Specificity")

Table1 <- bind_rows(SelectedGatesSpec100, SelectedGatesSpec90)
SelectedGatesSummary <- bind_rows(Table1, SelectedGatesROCO1)
