#' Expected Number of Events Since First Subject Randomized
#' 
#' This function calculates the expected number of events under alternative hypothesis
#' at calendar time t, which is calculated from first subject randomized. The function
#' returns the expected number of events for each arm at time t, based on the provided
#' enrollment distribution function and random lost-to-followup distribution if applicable.
#' If the total sample size is not provided, then only the corresponding probability of event
#' for each arm is provided.
#' 
#' @param DCO Analysis time calculated from first subject randomization date.
#' @param r  Randomization ratio of experimental arm : control arm as r:1. When r = 1, it is equal allocation. Default r = 1.
#' @param h0 Hazard function of control arm. h0(t) = log(2)/m0 means T~exponential distribution with median m0.
#' @param S0 Survival function of control arm. In general, S0(t) = exp(- integral of h0(u) for u from 0 to t).
#'           but providing S0(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f0(t) = h0(t) * S0(t).
#' @param h1 Hazard function of experimental arm. h1(t) = log(2)/m1 means T~exponential distribution with median m0.
#' @param S1 Survival function of experimental arm. In general, S1(t) = exp(- integral of h1(u) for u from 0 to t).
#'           but providing S1(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f1(t) = h1(t) * S1(t).
#' @param Lambda Distribution function of enrollment. For uniform enrollment, 
#' Lambda(t) = (t/A) where A is the enrollment period, i.e., Lambda(t) = t/A for 0<=t<=A, and 
#' Lambda(t) = 1 when t > A. For more general non-uniform enrollment with weight psi, 
#' Lambda(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). Default Lambda is uniform distribution function.
#' @param G0 Distribution function of lost-to-follow-up censoring process for control arm.
#' @param G1 Distribution function of lost-to-follow-up censoring process for experimental arm.
#' @param n Total sample size for two arms.  
#'
#' @return An object with a dataframe below.
#'  \describe{
#'       \itemize{
#'       \item ne0: number of events for control group
#'       \item ne1: number of events for experimental group
#'       \item ne: total number of events for two groups
#'       }
#'  }
#'  
#' @examples 
#' #Example (1) Trial scenario: 1:1 randomization, n = 450, enrollment follows non-uniform 
#' enrollment distribution with weight 1.5 and enrollment period is 18 months. 
#' Control arm ~ exponential distribution with median 12 months, and 
#' Experimental arm ~ exponential distribution (Proportional Hazards) with median 12 / 0.7 months.
#' Assuming no lost-to-followup. Find the expected number of events at calendar time 24 months, i.e.
#' 6 months after last patient randomized.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' Hazard function and survival function for experimental arm
#' lambda1 = lambda0 * HR
#' h1 = function(t){lambda1}; S1= function(t){exp(-lambda1 * t)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' max.DCO = 30
#' nE = matrix(NA, nrow = max.DCO, ncol=3)
#' for (DCO in 1:max.DCO) {
#'   o = fe(DCO = DCO, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450)
#'   nE[DCO, 1] = o$e0; nE[DCO, 2] = o$e1; nE[DCO, 3] = o$e
#' }
#' plot(1:max.DCO, nE[,3], type="n", xlab="Months", ylab = "Num of events")   
#' lines(1:max.DCO, nE[, 1], lty = 1, col=1)
#' lines(1:max.DCO, nE[, 2], lty = 2, col=2)
#' lines(1:max.DCO, nE[, 3], lty = 3, col=3)
#' legend(0, max(nE), c("Control", "Experimental", "Total"), col=1:3, lty=1:3, bty="n", cex=0.8)
#' 
#' #Example (2) Same trial set up as example (1) but assuming delayed effect for 
#' experimental arm. The delayed period is assumed 6 months, and after delay the
#' hazard ratio is assumed 0.65.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' Hazard function and survival function for experimental arm
#' h1 = function(t){lambda0*as.numeric(t < delay)+HR*lambda0*as.numeric(t >= delay)}
#' c = exp(-delay*lambda0*(1-HR)); 
#' S1 = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' max.DCO = 30
#' nE = matrix(NA, nrow = max.DCO, ncol=3)
#' for (DCO in 1:max.DCO) {
#'   o = fe(DCO = DCO, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450)
#'   nE[DCO, 1] = o$e0; nE[DCO, 2] = o$e1; nE[DCO, 3] = o$e
#' }
#' plot(1:max.DCO, nE[,3], type="n", xlab="Months", ylab = "Num of events")   
#' lines(1:max.DCO, nE[, 1], lty = 1, col=1)
#' lines(1:max.DCO, nE[, 2], lty = 2, col=2)
#' lines(1:max.DCO, nE[, 3], lty = 3, col=3)
#' legend(0, max(nE), c("Control", "Experimental", "Total"), col=1:3, lty=1:3, bty="n", cex=0.8)
#' 
#' @export

fe = function(n = 450, DCO = 24, r = 1, 
              h0 = function(t){log(2)/12}, 
              S0=function(t){exp(-log(2)/12 * t)}, 
              h1 = function(t){log(2)/12*0.70}, 
              S1= function(t){exp(-log(2)/12 * 0.7 * t)},
              Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
              G0 = function(t){0}, G1 = function(t){0}){

  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 
  
  #Density function
  f0 = function(t) {return(h0(t) * S0(t))}
  f1 = function(t) {return(h1(t) * S1(t))}
  
  #Integrand of control arm and experimental arm for calculation of prob. of event
  I0 = function(t){Lambda(DCO-t) * (1 - G0(t)) * f0(t)}
  I1 = function(t){Lambda(DCO-t) * (1 - G1(t)) * f1(t)}
  
  #prob. of event for control and experimental arm
  pe0 = integrate(I0, lower=0, upper=DCO)$value
  pe1 = integrate(I1, lower=0, upper=DCO)$value

  pe = r0 * pe0 + r1 * pe1
  
  #expected number of events
  if(!is.null(n)){
    e0 = n * Lambda(DCO) * r0 * pe0
    e1 = n * Lambda(DCO) * r1 * pe1
    e = e0 + e1
  }
  
  e = data.frame(cbind(e0, e1, e, pe0, pe1, pe))
  
  return(e)
}
