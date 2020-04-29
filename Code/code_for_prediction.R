## 2020-04-24
## Han Yan
## 
## Prediction function for vSIR Model


#' vSIR_Predict: function of prediction
#' @param origin.df: origin data of all provinces
#' @param est.df: estimation results of all provinces
#' @param w: regression window
#' @param D.R recover time
#' @param t.length: time range of the projection
#' @param B: times of the bootstrap

vSIR_Predict = function(origin.df, est.df, w = 5, D.R = 17.5,  t.length = 1000, B = 100){
  pop.df = read.csv("Data/province_population.csv", head = TRUE, fileEncoding = 'UTF-8')
  pv.names = as.character(unique(est.df$province))
  pred.sum = matrix(0, nrow = length(pv.names), ncol = 5)
  for(p in 1 : length(pv.names)){
    pv = pv.names[p]
    pop = pop.df$population[which(pop.df$province == pv)] * 10000
    temp.3 = est.df %>% filter(province == pv)
    temp.nosmooth = origin.df %>% filter(province == pv) 
    temp.nosmooth = temp.nosmooth[-(1 : (w - 1)),] 
    temp.3.dim = dim(temp.3)
    
    I.ns = temp.nosmooth$infected
    R.ns = temp.nosmooth$RD
    S.ns = pop - I.ns - R.ns # We will use the last day of them as the initial value of the ODE
    
    I.origin = temp.3$infected # We use the smoothed record to do bootstrap
    R.origin = temp.3$RD
    S.origin = pop - I.origin - R.origin
    
    beta.origin = temp.3$betas # the estimated varying beta
    beta.origin = sapply(beta.origin,thresh.0 ) # Some beta are negative and should be turned to be zero 
    gamma.origin = temp.3$gamma # the estimated varying gamma
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
    
    boot.peak.day = c()
    boot.end.day = c()
    boot.total.num = c()
    
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
    
    bias = rowMeans(varycoef.boot) - beta.origin[-c(1 : (w - 1))]
    beta.correction = varycoef.boot - bias
    
    bias.gamma = rowMeans(varygamma.boot) - gamma.origin[-c(1 : (w - 1))]
    gamma.correction = varygamma.boot - bias.gamma
    
    t.range.new = dim(beta.correction)[1]
    
    boot.predict = matrix(0.0, nrow = B, ncol = t.length+t.range-1)
    
    
    for(b in 1:B){
      beta.temp.correction = beta.correction[,b]
      beta.temp.correction=beta.temp.correction[(length(beta.temp.correction)-7):length(beta.temp.correction)]
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
      gamma.current = max(gamma.current, 1/D.R,na.rm = TRUE)
      # gamma.current should not below 1/14. 
      sim = simulator(gamma = gamma.current, a = a_new, b = b_new, eta= eta_new,
                      t1 = t.current.new, I1 = I.ns[length(I.ns)],
                      R1 = R.ns[length(I.ns)], max_t = t.length, N = pop)
      I.process = temp.3$infected 
      I.process[t.range:(t.range+t.length-1)] = sim$I
      
      #boot.end.day[b] = min(which(I.process < (max(I.process))/1000),na.rm = TRUE)
      boot.end.day[b] = min(which(I.process < 1),na.rm = TRUE)
      
      boot.peak.day[b] <- which.max(I.process)
      boot.total.num[b] <- sim$R[t.length]
      # boot.peak.num[b] = max(I.process)
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
    
    
    
    # get the confidence intervals and standard errors
    peak.CI <- ceiling(quantile(boot.peak.day, probs = c(0.025, 0.975),na.rm = TRUE))
    end.CI <- ceiling(quantile(boot.end.day, probs = c(0.025, 0.975),na.rm = TRUE))
    total.CI <- ceiling(quantile(boot.total.num, probs = c(0.025, 0.975),na.rm = TRUE))
    
    peak_time_L = as.character(peak.CI[1]+temp.3$date[1]-1)
    peak_time_U = as.character(peak.CI[2]+temp.3$date[1]-1)
    end_time_L = as.character(end.CI[1]+temp.3$date[1]-1)
    end_time_U = as.character(end.CI[2]+temp.3$date[1]-1)
    final_infected_L = total.CI[1]
    final_infected_U = total.CI[2]
    current_infected = ceiling(temp.nosmooth$infected[temp.3.dim[1]]+temp.nosmooth$RD[temp.3.dim[1]])
    
    
    peak_time = as.character(sprintf("%s/%s - %s/%s",month(peak_time_L),day(peak_time_L),
                                     month(peak_time_U),day(peak_time_U)))
    end_time = as.character(sprintf("%s/%s - %s/%s",month(end_time_L),day(end_time_L),
                                    month(end_time_U),day(end_time_U)))
    final_num = sprintf("%s - %s", final_infected_L,final_infected_U)
    
    pred.sum[p, 1] = pv; pred.sum[p, 2] = peak_time; pred.sum[p, 3] = end_time; 
    pred.sum[p, 4] = final_num; pred.sum[p, 5] = current_infected
    #print(pv) # for debugging
  }
  colnames(pred.sum) = c("province", "peak.time", "end.time", "final.size", "current infected")
  pred.sum = as.data.frame(pred.sum)
  return(pred.sum)
}


