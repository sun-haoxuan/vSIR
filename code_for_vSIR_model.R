## 2020-02-19
## Haoxuan Sun
## 
## Functions for vSIR Model
library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales)

## Moving Average
mav <- function(x){
  w=c(0.3,0.4,0.3)
  y=stats::filter(x, w, sides=2)
  y[1]=(3/10)*x[2]+(7/10)*x[1]
  y[length(y)]=(3/10)*x[length(x)-1]+(7/10)*x[length(x)]
  y=as.double(y)
  return(y)
}

## get Data
vSIR_Data <- function(opendir, smooth=TRUE){
  A = read.csv(opendir, head = TRUE, fileEncoding = 'UTF-8')
  if(A$city[1] == "total"){
    A1 = A[A$infected != 0, ] %>% select(province, date, infected, dead, recovered)
  }else{
    A1 = A[A$infected != 0, ] %>% select(city, date, infected, dead, recovered)
    A1$date = c(sub(".csv", "", A1$date))
    colnames(A1)[1] <- "province"
  }
  A1 = A1 %>% mutate(date = as.character(date) %>% mdy()) %>% mutate(RD = dead + recovered)  # Recoverd+Dead
  A1$infected = A1$infected-A1$RD
  
  if(smooth){
    A1$infected=mav(A1$infected)
    A1$dead=mav(A1$dead)
    A1$recovered=mav(A1$recovered)
    A1$RD=mav(A1$RD)
  }
  
  return(A1)
}

## Time-varying coefficient SIR model
vSIR_Vary <- function(A1, w=5, delta=1){

  A1 = A1 %>% mutate(infected.fit = 0)
  
  ## beta, gamma  
  beta.temp = c()
  SE.temp = c()
  gamma.temp = c()
  for (j in w : dim(A1)[1]){
    index = c((j - w + 1) : j)
    Data.temp = A1[index, ]
    
    ## gamma
    infected.diff = as.numeric(diff(Data.temp$infected,lag=delta))
    infected.lag = Data.temp$infected[1:(nrow(Data.temp)-delta)]*delta
    RD.diff = as.numeric(diff(Data.temp$RD,lag=delta))
    fit.infected.diff = lm(infected.diff ~ 0 + infected.lag)
    fit.RD.diff = lm(RD.diff ~ 0 + infected.lag)
    gamma.temp[j - w + 1]=coef(fit.RD.diff)[1]
    
    if(Data.temp$infected[length(Data.temp$infected)]==0){
      beta.temp[j - w + 1] = NA
      SE.temp[j - w + 1] = NA
      A1$infected.fit[j] = 0
    }else{
      fit = lm(log(Data.temp$infected) ~ c(1 : w))
      beta.temp[j - w + 1] = coef(fit)[2] + coef(fit.RD.diff)[1]
      SE.temp[j - w + 1] = (summary(fit)$coefficients[2, 2]^2+
                              summary(fit.RD.diff)$coefficients[2]^2)^(1/2)
      if (j == w){
        A1$infected.fit[1 : w] = exp(predict(fit))
      } else {A1$infected.fit[j] = exp(predict(fit)[w])}
    }
    
  }
  cov.temp = cov(beta.temp,gamma.temp)
  coef.one = data.frame(date = A1$date[w : dim(A1)[1]], 
                        beta = beta.temp,
                        beta.CI.Low = beta.temp-1.96*(SE.temp^2+cov.temp^2)^(1/2),
                        beta.CI.Up = beta.temp+1.96*(SE.temp^2+cov.temp^2)^(1/2),
                        beta.SE = (SE.temp^2+cov.temp^2)^(1/2),
                        gamma = gamma.temp) %>% mutate(date = as.character(date))
  coef.one$beta[which(coef.one$beta<0)]=0
  coef.one$beta.CI.Low[which(coef.one$beta.CI.Low<0)]=0
  coef.one$beta.CI.Up[which(coef.one$beta.CI.Up<0)]=0
  
  return(list('data'=A1, 'beta'=coef.one))
}

## Merge All Results
vSIR_Alldata <- function(opendata, smooth=TRUE, w=5, delta=1){
  files <- dir(opendata)
  l1 <- data.frame(); l2 <- data.frame(); l3 <- data.frame();
  for (ct in files){
    opendir <- paste0(opendata,ct)
    
    dat0 = vSIR_Data(opendir, smooth=FALSE)
    l1 = rbind(l1,dat0)
    
    dat = vSIR_Data(opendir)
    l2 = rbind(l2,dat)
    
    vbs = vSIR_Vary(dat,w,delta)$beta
    vbs$province=dat$province[1]
    l3 = rbind(l3,vbs)
  }
  
  l1 <- l1[which(!l1$date %in% c(as.Date("2020-01-11"),as.Date("2020-01-16"),
                                 as.Date("2020-01-17"),as.Date("2020-01-18"))),]
  l2 <- l2[which(!l2$date %in% c(as.Date("2020-01-11"),as.Date("2020-01-16"),
                                 as.Date("2020-01-17"),as.Date("2020-01-18"))),]
  colnames(l3)=c('date','betas','betas.L','betas.U','betas.SE','gamma','province')
  
  return(list('dat0'=l1,'dat'=l2,'vbs'=l3))

}

