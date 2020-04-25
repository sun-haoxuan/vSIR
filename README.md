# COVID-19-vSIR-model
This repository includes the R code for vSIR model in our paper "Tracking Reproductivity of COVID-19 Epidemic in China with
Varying Coefficient SIR Model".

* basic_functions.R

The code provide functions for smoothing, reciprocal fitting and ODE simulation. Here we use a three point smoothing with a weight $3:4:3$ and $7:3$ and $3:7$ for the first and last day. The reciprocal function $\hat{\beta_t} = b/t^{\eta} - a$ fitted by non-linear minimization.

* code_for_vSIR_model.R

The code is for reproducing the estimation of $(\beta(t)$, $\gamma(t)$ and $R^D(t))$. The estimates are based on the infected and removed cases each day. See *Sun et al.(2020)* for method and theory details. 

* code_for_prediction.R

The code is for reproducing the prediction result. The prediction for the peak time, ending time and final infected number are based on reciprocal fitting model and simulation function built with the ODE. The prediction intervals are constructed with the Bootstrap method. See *Sun et al.(2020)* for method and theory details.
