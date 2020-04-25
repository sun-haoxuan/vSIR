rm(list = ls())
library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales)
library(pinyin)
library(smooth)
library(Mcomp)

source('Code/basic_functions.R')
source('Code/code_for_vSIR_model.R')
source('Code/code_for_prediction.R')

## Results of vSIR model
lst = vSIR_Alldata('Data/nCoV_Pinyin0220/', smooth=TRUE, w=5, delta=1)

## Prediction results 
l1 = lst$dat0; l2 =lst$dat; l3 = lst$vbs
l4 = merge(l2, l3, by = c("province", "date"))
pred.result = vSIR_Predict(l1, l4)
