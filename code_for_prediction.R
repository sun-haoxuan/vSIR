## Created on March 24 2020 by Han Yan

## Code for making prediction of the peak time, ending time and final number et.al.

library(tidyverse)
library(lubridate)
library(pinyin)
library(smooth)
library(Mcomp)
library(ggplot2)

# ---------- Part 1. Data preparation ---------- #

path0 = "Data/nCov_pinyin0220/"

mav <- function(x){
  w=c(0.3,0.4,0.3)
  y=stats::filter(x, w, sides=2)
  y[1]=(3/10)*x[2]+(7/10)*x[1]
  y[length(y)]=(3/10)*x[length(x)-1]+(7/10)*x[length(x)]
  y=as.double(y)
  return(y)
}
# province name inculde 30 provinces + 15 Hubei cities
All.data = read.csv("Data/Ets.csv")
province.name = as.character(All.data$pv)
province.name= unique(province.name)
cityname = c("Wuhan","Huangshi","Shiyan","Yichang","Xiangyang","Ezhou","Jingmen","Xiaogan","Jingzhou",
             "Huanggang","Xianning","Suizhou","Enshizhou","Xiantao","Tianmen")

Infect.data = data.frame()
Infect.data.2 = data.frame()
infected.all = 0
dead.all = 0
recovered.all = 0
infected.all1 = 0
dead.all1 = 0
recovered.all1 = 0

for (i in 1 : length(province.name)){
  path = paste(path0, province.name[i], ".csv", sep = "")
  A = read.csv(path, head = TRUE)
  # A$infected = A$infected - A$dead - A$recovered
  date1 = A$date
  infected.all = infected.all + A$infected
  dead.all = dead.all + A$dead
  recovered.all = recovered.all + A$recovered
  if (i > 1){
    infected.all1 = infected.all1
    dead.all1 = dead.all1 + A$dead
    recovered.all1 = recovered.all1 + A$recovered
  }
  A1 = A[A$infected != 0, ]  %>% select(province, date, infected, dead, recovered)
  A1$infected = A1$infected - A1$dead - A1$recovered
  A1$province = province.name[i]
  A1 = A1 %>% mutate(date = as.character(date) %>% mdy()) %>% mutate(RD = dead + recovered) %>% mutate(It = infected - RD)
  #A1 = A1 %>% filter(date < as.Date("2020-02-11"))
  #A1 = A1 %>% mutate(province = py(as.character(province), dic = pydic((method = "quanpin"), dic = "pinyin2")))
  Infect.data = rbind(Infect.data, A1)
  
  A2 = A1
  A1$smoothed.infected.2 = mav(A1$infected)
  A1$smoothed.RD.2 = mav(A1$RD)
  A1$smooth.It.2 = mav(A1$It)
  
  ## copy smoothed data back to the matrix
  A2$infected = A1$smoothed.infected.2
  A2$RD = A1$smoothed.RD.2
  A2$It = A1$smooth.It.2
  Infect.data.2 = rbind(Infect.data.2, A2)
}

Infect.data$province = gsub("_","",Infect.data$province)
Infect.data.2$province = gsub("_","",Infect.data.2$province)


# ---------- Part 2. Parameter estimation ---------- #

w = 5 # window size
Beta.varycoef = data.frame(date = as.character(Infect.data.2$date[1 : 10])) %>% mutate(date = as.character(date))
province.name = unique(Infect.data.2$province)

for (i in 1 : length(province.name)){
  A4 = Infect.data.2[Infect.data.2$province == province.name[i], ]
  
  beta.temp = c()
  gamma.temp = c()
  
  for (j in w : dim(A4)[1]){
    index = c((j - w + 1) : j)
    Data.temp = A4[index, ]
    
    RD.diff2 = diff(Data.temp$RD, lag = 1)
    infected.lag2 = Data.temp$infected[-dim(Data.temp)[1]]
    fit.gamma = lm(RD.diff2 ~ 0 + infected.lag2)
    gamma.temp[j - w + 1] = coef(fit.gamma)
    #gamma.temp[j-w+1] = (Data.temp$RD[w]-Data.temp$RD[1])/sum(Data.temp$infected[2:w])
    fit = lm(log(Data.temp$infected + 0.00001) ~ c(1 : w))
    beta.temp[j - w + 1] = coef(fit)[2] + gamma.temp[j- w + 1]
  }
  beta.one = data.frame(date = A4$date[w : dim(A4)[1]], beta = beta.temp, gamma = gamma.temp) %>% mutate(date = as.character(date))
  names(beta.one) = c("date", province.name[i], paste(province.name[i], paste(province.name[i], "gamma", sep = "-")))
  Beta.varycoef = full_join(Beta.varycoef, beta.one, by = 'date')
}


Beta.varycoef.reshape=data.frame()
for (j in 1 : length(province.name)){
  index.j = 2 * (j - 1) + 2
  Beta.varycoef.reshape = rbind(Beta.varycoef.reshape, data.frame( province = names(Beta.varycoef)[index.j], date = Beta.varycoef[, 1], beta = Beta.varycoef[, index.j], gamma = Beta.varycoef[, index.j + 1] ) )
}
Beta.varycoef.reshape=Beta.varycoef.reshape[!is.na(Beta.varycoef.reshape$beta),]
Beta.varycoef.reshape$province = as.character(Beta.varycoef.reshape$province)

Beta.varycoef.reshape$date=as.Date(Beta.varycoef.reshape$date)


# Merge record and beta hat

Infect.data.3 =  full_join(Infect.data.2, Beta.varycoef.reshape, by = c("province","date"))
Infect.data.3 = Infect.data.3 %>% select(province,date,infected, RD, beta,gamma) %>%filter(!is.na(beta))

# ---------- Part 3. Functions for the reciprocal fitting model ---------- #

SSE = function(x, y, power, para){
  a = para[1]
  b = para[2]
  c = power
  res = sum((y - b / (x^c - a))^2)
  return(res)
}

least.square = function(x, y, power, para0){
  res = nlm(SSE, para0, x = x, y = y, power = power, stepmax = 10^5)
  return(c(res$estimate, power, res$code, res$minimum))
}

min.power = function(x, y, power, para0){
  min.value = c()
  for (i in 1 : length(power)){
    res.temp = least.square(x, y, power[i], para0)
    min.value[i] = res.temp[5]
  }
  power.min = power[order(min.value)[1]]
  res.final = least.square(x, y, power.min, para0)
  return(res.final)
}


# ------ Part 4. Simulation function built with the ODE ------ #

simulator <- function(gamma, a, b, eta, t1, I1, R1, max_t, N){
  
  ODEfun <- function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      dS <- - (b /((t + t1 )^eta -a)) *I*S/N
      dI <- (b /((t + t1 )^eta -a) )*I*S/N-gamma*I
      dR <- gamma*I
      list(c(dS,dI,dR))
    })
  } 
  
  parameters <- c(a = a, b = b, eta = eta, gamma = gamma, t1 = t1)
  state <- c(S = N-I1,I = I1, R = R1)
  times <- 0:(max_t-1)
  
  # Simulation trajectories of SIE model
  sim <- deSolve::ode(y = state, times = times, func = ODEfun, parms = parameters)
  list(t = sim[,1], S = sim[,2], I = sim[,3], R = sim[,4])
  
}

# ------ Other functions ------ #

quantile.L = function(i){
  L = boot.predict[,i]
  l = quantile(L, probs = 0.05)
}
quantile.U = function(i){
  U = boot.predict[,i]
  u = quantile(U, probs = 0.95)
}

thresh.0 = function(x){
  return(max(x,0))
}

# Population data
province.pop <- read.csv("/Data/province_population.csv")
province.pop$province <- sapply(province.pop$province, as.character)

remove.province = c("Xianggang", "Aomen","Taiwan","All provinces")
province.name = province.name[! province.name %in% remove.province]

province.name.2 = c("Hubei","All provinces","All without Hubei")


# ---------- Part 5. Prediction based on above results ---------------- #

t.length = 1000

B = 100

peak_time_L = c()
peak_time_U = c()
end_time_L = c()
end_time_U = c()

stop1_time_L = c()
stop1_time_U = c()
stop14_time_L = c()
stop14_time_U = c()


final_infected_L = c()
final_infected_U = c()
current_infected = c()
R0_province = c()
predict.df.3 = data.frame()
peak_time = c()
end_time= c()
stop_time1 = c()
stop_time14= c()
final_num = c()
#province.name = province.name


for(j in 1:length(province.name)){ 
  
  pop = province.pop[which(province.pop$province == province.name[j]),2]*10000
  
  temp.3 = Infect.data.3 %>% filter(province == province.name[j])
  temp.nosmooth = Infect.data %>% filter(province == province.name[j]) # nosmooth(Infect.data) is the original data
  temp.nosmooth = temp.nosmooth[-(1:4),] # 减去一个window(没有beta的估计)
  temp.3.dim = dim(temp.3)
  
  I.ns = temp.nosmooth$infected
  R.ns = temp.nosmooth$RD
  S.ns = pop - I.ns - R.ns # We will use the last day of them as the initial value of the ODE
  
  I.origin = temp.3$infected # We use the smoothed record to do bootstrap
  R.origin = temp.3$RD
  S.origin = pop - I.origin - R.origin
  
  beta.origin = temp.3$beta # the estimated varying beta
  beta.origin = sapply(beta.origin,thresh.0 ) # Some beta are negative and should be turned to be zero 
  gamma.origin =temp.3$gamma # the estimated varying gamma
  t.range = temp.3.dim[1]
  
  boot.I <- matrix(0, nrow = t.range, ncol = B) # boot.I(R.S) are matrices which record the bootstrap samples
  boot.R <- matrix(0, nrow = t.range, ncol = B) 
  boot.S <- matrix(0, nrow = t.range, ncol = B) 
  boot.S[1, ] <- rep(S.origin[1], B) 
  boot.I[1, ] <- rep(I.origin[1], B) 
  boot.R[1, ] <- rep(R.origin[1], B)
  
  
  for (b in 1 : B){
    for (i in 1 : (t.range - 1)){
      be = beta.origin[i]
      ga = gamma.origin[i]
      boot.S[(i + 1), b] = boot.S[i, b] - rpois(1, lambda = be* boot.S[i, b] * boot.I[i, b] / pop+0.00001) 
      # here + 0.00001 is to aviod lamba = 0 
      boot.R[(i + 1), b] = boot.R[i, b] + rpois(1, lambda = boot.I[i, b] * ga+0.00001)
      boot.I[(i + 1), b] = max(pop - boot.S[(i + 1), b] - boot.R[(i + 1), b],0)
    }
  }
  
  boot.end.day = c()
  boot.stop1.day = c()
  boot.stop7.day = c()
  boot.peak.day = c()
  boot.total.num = c()
  boot.peak.num = c()
  
  w = 5
  varycoef.boot = matrix(0, dim(temp.3)[1] - (w - 1), B)
  varygamma.boot = matrix(0, dim(temp.3)[1] - (w - 1), B)
  
  for(b in 1:B){
    A4 = data.frame(infected = boot.I[, b], RD = boot.R[, b])
    A4$date = temp.3$date
    
    beta.temp = c()
    gamma.temp = c()
    
    for (jj in w : t.range){
      index = c((jj - w + 1) : jj)
      Data.temp = A4[index, ]
      
      RD.diff2 = diff(Data.temp$RD, lag = 1)
      infected.lag2 = Data.temp$infected[-dim(Data.temp)[1]]
      fit.gamma = lm(RD.diff2 ~ 0 + infected.lag2, na.action = na.omit)
      gamma.temp[jj - w + 1] = coef(fit.gamma)
      #gamma.temp[jj -w +1] = (Data.temp$RD[w] - Data.temp$RD[1])/sum(Data.temp$infected[2:w])
      
      fit = lm(log(Data.temp$infected+0.00001) ~ c(1 : w))
      beta.temp[jj - w + 1] = coef(fit)[2] + gamma.temp[jj - w+1]
    }
    varycoef.boot[, b] = beta.temp
    varygamma.boot[, b] = gamma.temp
  }
  
  bias = rowMeans(varycoef.boot) - temp.3$beta[-c(1 : (w - 1))]
  beta.correction = varycoef.boot - bias
  
  bias.gamma = rowMeans(varygamma.boot) - temp.3$gamma[-c(1 : (w - 1))]
  gamma.correction = varygamma.boot - bias.gamma
  
  t.range.new = dim(beta.correction)[1]
  
  boot.predict = matrix(0.0, nrow = B, ncol = t.length+t.range-1)
  
  
  for(b in 1:B){
    beta.temp.correction = beta.correction[,b]
    beta.temp.correction=beta.temp.correction[length(beta.temp.correction-7):length(beta.temp.correction)]
    # Use the last 8 days to fit the reciprocal model 
    beta.temp.correction = sapply(beta.temp.correction, thresh.0)
    
    time.index = 1 : length(beta.temp.correction)
    fit.reciprocal = min.power(time.index, beta.temp.correction, seq(0.5, 5, 0.1), c(0, 1))
    a_new = fit.reciprocal[1]
    b_new = fit.reciprocal[2]
    b_new = thresh.0(b_new)
    eta_new = fit.reciprocal[3]  
    t.current.new = length(beta.temp.correction)
    
    gamma.current = gamma.correction[dim(gamma.correction)[1],b]
    gamma.current = max(gamma.current, 1/14,na.rm = TRUE)
    # gamma.current should not below 1/14. 
    
    sim = simulator(gamma = gamma.current, a = a_new, b = b_new, eta= eta_new,
                    t1 = t.current.new, I1 = I.ns[t.range],
                    R1 = R.ns[t.range], max_t = t.length, N = pop)
    
    N.process = temp.nosmooth$infected + temp.nosmooth$RD
    N.process[t.range:(t.range + t.length - 1)] = sim$I + sim$R
    N.diff1 = diff(N.process,lag = 1)
    #N.diff7 = diff(N.process,lag = 7)
    N.diff14 = diff(N.process,lag = 14)
    I.process = temp.3$infected 
    I.process[t.range:(t.range+t.length-1)] = sim$I
    
    #boot.end.day[b] = min(which(I.process < (max(I.process))/1000),na.rm = TRUE)
    boot.end.day[b] = min(which(I.process < 1),na.rm = TRUE)
    boot.stop1.day[b] = min(which(N.diff1 <1 ),na.rm = TRUE)
    boot.stop14.day[b] = min(which(N.diff14 <1),na.rm = TRUE)+13
    
    boot.peak.day[b] <- which.max(I.process)
    boot.total.num[b] <- sim$R[t.length]
    boot.peak.num[b] = max(I.process)
    boot.predict[b,] = I.process  
    
  }  
  
  quantile.L = function(i){
    L = boot.predict[,i]
    l = quantile(L, probs = 0.025)
  }
  quantile.U = function(i){
    U = boot.predict[,i]
    u = quantile(U, probs = 0.975)
  }
  
  
  #predict.temp = data.frame(date = seq(temp.3$date[1],temp.3$date[t.range]+t.length-1, "days"),
  #                        I.mean = colMeans(boot.predict), 
  #                        I.L = sapply(1:(t.length+t.range-1), quantile.L), 
  #                        I.U = sapply(1:(t.length+t.range-1), quantile.U),
  #                        province = rep(province.name[j],t.length+t.range-1))
  #predict.temp = predict.temp[1:200,]    
  
  #predict.df.3 = rbind(predict.df.3, predict.temp)
  
  
  # get the confidence intervals and standard errors
  peak.CI <- ceiling(quantile(boot.peak.day, probs = c(0.025, 0.975),na.rm = TRUE))
  peak.sd <- sd(boot.peak.day)
  end.CI <- ceiling(quantile(boot.end.day, probs = c(0.025, 0.975),na.rm = TRUE))
  stop1.CI = ceiling(quantile(boot.stop1.day, probs = c(0.025, 0.975),na.rm = TRUE))
  stop14.CI = ceiling(quantile(boot.stop14.day, probs = c(0.025, 0.975),na.rm = TRUE))
  peak.sd <- sd(boot.peak.day)
  end.sd <- sd(boot.end.day)
  total.CI <- ceiling(quantile(boot.total.num, probs = c(0.025, 0.975),na.rm = TRUE))
  total.sd <- sd(boot.total.num)
  
  #R0_province[j] = temp.3$R0_2_wks[temp.3.dim[1]]
  peak_time_L[j] = as.character(peak.CI[1]+temp.3$date[1]-1)
  peak_time_U[j] = as.character(peak.CI[2]+temp.3$date[1]-1)
  end_time_L[j] = as.character(end.CI[1]+temp.3$date[1]-1)
  end_time_U[j] = as.character(end.CI[2]+temp.3$date[1]-1)
  stop1_time_L[j] = as.character(stop1.CI[1] + temp.3$date[1]-1)
  stop1_time_U[j] = as.character(stop1.CI[2] + temp.3$date[1]-1)
  stop14_time_L[j] = as.character(stop14.CI[1] + temp.3$date[1]-1)
  stop14_time_U[j] = as.character(stop14.CI[2] + temp.3$date[1]-1)
  
  final_infected_L[j] = total.CI[1]
  final_infected_U[j] = total.CI[2]
  current_infected[j] = ceiling(temp.nosmooth$infected[temp.3.dim[1]]+temp.nosmooth$RD[temp.3.dim[1]])
  
  
  peak_time[j] = as.character(sprintf("%s/%s - %s/%s",month(peak_time_L[j]),day(peak_time_L[j]),
                                      month(peak_time_U[j]),day(peak_time_U[j])))
  
  end_time[j] = as.character(sprintf("%s/%s - %s/%s",month(end_time_L[j]),day(end_time_L[j]),
                                     month(end_time_U[j]),day(end_time_U[j])))
  stop_time1[j] = as.character(sprintf("%s/%s - %s/%s",month(stop1_time_L[j]),day(stop1_time_L[j]),
                                       month(stop1_time_U[j]),day(stop1_time_U[j])))
  stop_time14[j] = as.character(sprintf("%s/%s - %s/%s",month(stop14_time_L[j]),day(stop14_time_L[j]),
                                        month(stop14_time_U[j]),day(stop14_time_U[j])))
  
  final_num[j] = sprintf("%s - %s", final_infected_L[j],final_infected_U[j])
}


result= data.frame(province = province.name, peak.time = peak_time, end.time = end_time,  
                     stop.time.14 = stop_time14, final.num = final_num, current_infected = current_infected)

result.hubei = result %>% filter(province %in% cityname)
result.nohubei = result %>% filter(!province %in% cityname)

