## 2020-04-24
## Haoxuan Sun; Han Yan
## 
## Basic Functions


#' mav: function for smoothing
#' @param x: a time series to smooth


mav <- function(x){
  w=c(0.3,0.4,0.3)
  y=stats::filter(x, w, sides=2)
  y[1]=(3/10)*x[2]+(7/10)*x[1]
  y[length(y)]=(3/10)*x[length(x)-1]+(7/10)*x[length(x)]
  y=as.double(y)
  return(y)
}

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


## Simulation function built with the ODE
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


thresh.0 = function(x){
  return(max(x,0))
}

