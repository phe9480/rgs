#' Max-combo Test (Stratified) At A Single Time Point
#' 
#' Max-combo is a class of multiple combination tests. This function considers
#' the maximum of multiple weighted log-rank tests, i.e. z_maxcombo = max(z1, z2, ..., zk),
#' where zi is the weighted log-rank test statistic with weight function wi. 
#' This function calculates the p-value and associated statistics when performing
#' the max-combo test. The weight function wi can be the flexible Stabilized 
#' Fleming-Harrington class sFH(rho, gamma, tau, s.tau), or any user-defined 
#' weight function. 
#' 
#' Refer to Karrinson (2016) for the method of max-combo when 
#' it is defined based on FH(0, 0), FH(1, 0), FH(0, 1) and FH(1, 1). However, 
#' this function extends the concept to include any type of weighted log-rank tests.
#' The kth weighted log-rank test has parameters of (rho_k, gamma_k, tau_k, s.tau_k)
#' defined in the stabilized Fleming-Harrington class (He et al 2021) or user-defined 
#' weight function f.ws_k, for k = 1, ..., K. 
#' For the stabilized Fleming-Harrington class, specify either tau or s.tau, 
#' which are thresholds in survival time and survival rate, respectively. In addition,
#' the weight function is estimated within each strata. One-sided test is more appropriate
#' in maxcombo test.
#' 
#' @param  time  Survival time
#' @param  event Event indicator; 1 = event, 0 = censor
#' @param  group Treatment group; 1 = experimental group, 0 = control
#' @param  rho A vector of rho parameters for Fleming-Harrington (rho, gamma) weighted log-rank tests.
#' @param  gamma A vector of gamma parameters for Fleming-Harrington (rho, gamma) weighted log-rank tests.
#'         For log-rank test, set rho = gamma = 0.
#' @param  tau  A vector of cut points for stabilized FH test, sFH(rho, gamma, tau); with weight
#'       function defined as w(t) = s_tilda^rho*(1-s_tilda)^gamma, where
#'       s_tilda = max(s(t), s.tau) or max(s(t), s(tau)) if s.tau = NULL
#'       tau = Inf reduces to regular Fleming-Harrington test(rho, gamma)
#' @param  s.tau  A vector of survival rate cuts S(tau) at t = tau; default 0.5, ie. cut at median.
#'       s.tau = 0 reduces to regular Fleming-Harrington test(rho, gamma)
#'       
#' @param  f.ws  Self-defined weight functions of survival rate. 
#'         For example, f.ws=list(
#'          lr = function(s){return(1)}, 
#'          fh005=function(s){sqrt(1-s)},
#'          sfh01=function(s){1-apply(cbind(s, 0.5),MARGIN=1,FUN=max)}
#'          ) 
#'         defines 3 weighted log-rank tests: log-rank, FH(0, 0.5), and sFH(0, 1, s.tau=0.5)
#'         When f.ws is specified, the weight function takes it as priority.
#'         
#'         Note: The actual math formula for weighting function is based on 
#'         the left of each event time t, i.e., w(t) = f.ws(s(t-)). 
#'         For FH(rho, gamma) test, when gamma = 0, then the first event time has 
#'         weight 1;  when gamma > 0, the first event weight is 0. This can 
#'         ensure consistency with FH(0,0) = logrank; and FH(1,0) = generalized Wilcoxon.
#'         
#' @param  strata1 Stratification variable 1
#' @param  strata2 Stratification variable 2
#' @param  strata3 Stratification variable 3
#' @param  side  Type of test. one.sided or two.sided. default = two.sided
#' 
#' @return 
#' \describe{
#' \item{data}{dataframe including the following variables:}
#'          \itemize{
#'          \item time: Survival time 
#'          \item event: Event indicator; 0 = censor; 1 = event
#'          \item group: Group indicator; 0 = control; 1 = experimental treatment
#'          \item strata1; Strata 1 
#'          \item strata2; Strata 2 
#'          \item strata3; Strata 3 
#'          }
#' \item{corr}{Correlation matrix of weighted log-rank tests}        
#' \item{test.results}{dataframe including the following variables:}
#'          \itemize{
#'          \item rho:     Stabilized Fleming-Harrington test parameter
#'          \item gamma:   Stabilized Fleming-Harrington test parameter
#'          \item tau:     Stabilized Fleming-Harrington test parameter if provided
#'          \item s.tau:   Stabilized Fleming-Harrington test parameter if provided
#'          \item test.side: 1 or 2
#'          \item z:       Normalized weighted log-rank z statistics
#'          \item z.max:   maximum of weighted log-rank test statistics(1-sided test)
#'          \item abs.z.max maximum of absolute values of weighted log-rank test statistics (2-sided tests)
#'          \item p:       P value
#'          }
#' \item{wt}{Weight function(s) used in the calculation}
#' }
#' @examples
#' #Example (1) Stratified max-combo test including 3 weighted log-rank test statistics:
#' #log-rank, FH(0, 0.5), and sFH(0, 1, tau=0.5).
#' time=rexp(100); event=sample(c(0,1), 100, replace = TRUE); 
#' group=c(rep(0, 50), rep(1, 50)); 
#' strata1=sample(c(1,2), 100, replace = TRUE)
#' strata2=sample(c(1,2), 100, replace = TRUE) 
#' strata3=sample(c(3,4), 100, replace = TRUE)
#' rho = c(0,0,0); gamma=c(0,0.5,1); tau = NULL; s.tau=c(0, 0, 0.5);
#' f.ws=list(
#'          lr = function(s){return(1)}, 
#'          fh005=function(s){sqrt(1-s)},
#'          sfh01=function(s){1-apply(cbind(s, 0.5),MARGIN=1,FUN=max)}
#'          )
#' side = c("one.sided")
#' 
#' wlr.maxcombo(time=time, event=event, group=group, strata1=strata1, 
#' strata2=strata2, strata3=strata3, rho=rho, gamma=gamma, 
#' tau = tau, s.tau=s.tau, f.ws=f.ws, side = side)
#'
#' #Equivalent to 
#' wlr.maxcombo(time=time, event=event, group=group, strata1=strata1, 
#' strata2=strata2, strata3=strata3, rho=rho, gamma=gamma, 
#' tau = tau, s.tau=s.tau, f.ws=NULL, side = side)
#' 
#' #Example (2) If there is only 1 weighted log-rank test, max-combo reduces to weighted log-rank test
#' #FH01 test
#' wlr.maxcombo(time=time, event=event, group=group, strata1=strata1, 
#' strata2=strata2, strata3=strata3, rho=0, gamma=1, 
#' tau = NULL, s.tau=0, f.ws=NULL, side = side)
#' 
#' #equivalent to
#' wlr(time=time, event=event, group=group, strata1=strata1, 
#' strata2=strata2, strata3=strata3, rho=0, gamma=1, 
#' tau = NULL, s.tau=0, f.ws=NULL, side = side)
#' 
#' #Example (3) maxcombo of (logrank, FH01, FH11)
#' wlr.maxcombo(time=time, event=event, group=group, strata1=strata1, 
#' strata2=strata2, strata3=strata3, rho=c(0,0,1), gamma=c(0,1,1), 
#' tau = NULL, s.tau=c(0,0,0), f.ws=NULL, side = side)
#' 
#' #Example (4) maxcombo of (logrank, FH01)
#' wlr.maxcombo(time=time, event=event, group=group, strata1=strata1, 
#' strata2=strata2, strata3=strata3, rho=c(0,0), gamma=c(0,1), 
#' tau = NULL, s.tau=c(0,0), f.ws=NULL, side = side)
#' 
#' #Example (5) maxcombo of (logrank, modestly logrank)
#' wlr.maxcombo(time=time, event=event, group=group, strata1=strata1, 
#' strata2=strata2, strata3=strata3, rho=NULL, gamma=NULL, 
#' tau = NULL, s.tau=NULL, f.ws=list(
#' lr = function(s){1},
#' mlr = function(s){1/apply(cbind(s, 0.5),MARGIN=1,FUN=max)}), side = side)
#' 
#' @export
wlr.maxcombo = function(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1),
    group=c(0,1,0,1,0,1,0,1), strata1=NULL, strata2=NULL, strata3=NULL, 
    rho = c(0,0,0), gamma=c(0,0.5,1), tau = NULL, s.tau=c(0, 0, 0.5),
    f.ws=list(
      lr = function(s){return(1)}, 
      fh005=function(s){sqrt(1-s)},
      sfh01=function(s){1-apply(cbind(s, 0.5),MARGIN=1,FUN=max)}
    ), 
    side = c("two.sided", "one.sided")) {

  #CALCULATE the z-statistic for each of WLR tests. 
  #Number of weighted log-rank tests K.
  if (!is.null(f.ws)){K = length(f.ws)} else {K=length(rho)}
  z = rep(NA, K)
  for (k in 1:K) {
    wlr.k = wlr(time=time, event=event, group=group, strata1=strata1, 
                strata2=strata2, strata3=strata3, rho=rho[k], gamma=gamma[k], 
                tau = tau[k], s.tau=s.tau[k], f.ws=f.ws[[k]], side = side)
    z[k] = as.numeric(wlr.k$test.results$z)
  }
  
  if(K > 1){
    #correlation matrix of (z1, ..., zK)
    corr = matrix(1, nrow=K, ncol=K)
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        corr[i, j] = wlr.cov(time=time, event=event, group=group, strata1=strata1,
                  strata2=strata2, strata3=strata3, 
              rho1=rho[i], gamma1=gamma[i], tau1 = tau[i], s.tau1=s.tau[i],
              rho2=rho[j], gamma2=gamma[j], tau2 = tau[j], s.tau2=s.tau[j],
              f.ws1=f.ws[[i]], f.ws2=f.ws[[j]])$corr
        corr[j, i] = corr[i, j]
      }
    }

    if (side[1] == "one.sided"){
      test.side = 1
      z.max = max(z)
      p = 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(z.max, K), corr = corr, abseps = 1e-8, maxpts=100000)
      test.results =  data.frame(cbind(z, z.max, p, test.side))
    
    } else {
      test.side = 2
      abs.z.max = max(abs(z))    
      p = 2*(1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(abs.z.max, K), corr = corr, abseps = 1e-8, maxpts=100000))
      test.results =  data.frame(cbind(z, abs.z.max, p, test.side))
    }
    o=list()
    o$corr = corr
    o$test.results = test.results
    if(!is.null(f.ws)){wt = f.ws} else{
      wt = data.frame(cbind(rho, gamma, tau, s.tau))
    }
    o$wt = wt
  } else{
      #output consistent format as maxcombo
      tmp = wlr.k$test.results
      z = tmp$z; z.max = z; p = tmp$p; test.side = tmp$test.side 
      test.results = data.frame(cbind(z, z.max, p, test.side))
      o = wlr.k; o$test.results = test.results
    }
  return(o)     
}
