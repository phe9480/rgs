#' Streamlined Trial Data Simulations and Analysis Using Weighted Log-rank Test For Mixture Cure Rate Distribution
#'
#' Simulate Randomized two-arm trial data with the following characteristics:  
#' (1) randomization time (entry time) is generated according to the specified non-uniform accrual pattern, 
#' i.e. the cumulative recruitment at calendar time t is (t/A)^w with weight w and enrollment complete in A months.
#' w = 1 means uniform enrollment, which is usually not realistic due to graduate sites activation process.
#' (2) Survival time follows piece-wise exponential distribution for each arm.
#' (3) N total patients with r:1 randomization ratio
#' (4) Random drop off can be incorporated into the censoring process.
#' (5) Data cutoff dates are determined by specified vector of target events for all analyses.
#' (6) A dataset is generated for each analysis according to the specified number of target events. 
#'     Multiple analyses can be specified according to the vector of targetEvents, eg, targetEvents = c(100, 200, 300) 
#'     defines 3 analyses at 100, 200, and 300 events separately.
#' (7) Weighted log-rank test is then performed for each simulated group sequential dataset.
#'
#' @param  nSim Number of trials
#' @param  N Total number patients in two arms.
#' @param  A Total accrual period in months
#' @param  w Weight parameter in cumulative enrollment pattern. 
#' The cumulative enrollment at month t is (t / A)^w, eg, at month 6, 
#'   the enrollment is N*(6/24)^2 = N/16 for 24 months planned accrual period.
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' @param  lambda0 Hazard rates for control arm of intervals defined by cuts; for exponential(lambda0) distribution,
#'         lambda0 = log(2) / median;
#' @param  lambda1 Hazard rates for experimental arm for intervals; for exponential(lambda1) distribution,
#'         lambda1 = log(2) / median; For delayed effect under H1, lambda1 is a vector (below).
#' @param  cuts Timepoints to form intervals for piecewise exponential distribution. For example,
#   \itemize{
#   \item Proportional hazards with hr = 0.65. Then lambda0 = log(2)/m0, lambda1 = log(2)/m0*hr, cuts = NULL. 
#   \item Delayed effect at month 6, and control arm has constant hazard (median m0) and 
#'       experimental arm has hr = 0.6 after delay, then cuts = 6, and 
#'       lamda0 = log(2) / m0 or lambda0 = rep(log(2) / m0, 2), 
#'       lamda1 = c(log(2)/m0, log(2)/m0*hr). 
#   \item Delayed effect at month 6, and control arm has crossover to subsequent IO 
#'       treatment after 24 mo, so its hazard decreases 20%. Then, 
#'       lambda0 = c(log(2)/m0, log(2)/m0, log(2)/m0*0.8), 
#'       lambda1 = c(log(2)/m0, log(2)/m0*hr, log(2)/m0*hr), and
#'       cuts = c(6, 24), which forms 3 intervals (0, 6), (6, 24), (24, infinity)
#       }
#' @param dropOff0 Drop Off rate per month, eg, 1%, for control arm
#' @param dropOff1 Drop Off rate per month, eg, 1%, for experimental arm
#' @param targetEvents A vector of target events is used to determine DCOs. For example, 
#'              397 target events are used to determine IA DCO; and 496 events are used 
#'              to determine the FA cutoff.
#' @param logrank Indicator whether log-rank test is requested besides the weighted logrank tests. "Y" or "N". Default "Y".
#'                 If "Y", the traditional log-rank test will be used based on survdiff() function.
#'                 If "N", the weighted log-rank test with weighting function specified in fws will be used.
#' @fws.options   Weighting strategies in the following format as examples. If fws.options is provided,
#' then the nphDesign object's weighting strategy will be ignored. fws can contain multiple weighting strategies.
#'        For example, fws = list(fws1, fws2) means 2 weighting strategies are evaluated, where
#'        fws1 = list(IA = list(lr), FA=list(lr, fh01)); fws2 = list(IA = list(lr), FA=list(fh01)). Each
#'        weighting strategy is specified as following examples for illustration.
#'        \itemize{
#'        \item (1) Single-time analysis using log-rank test: fws1 = list(FA = list(lr));
#'        \item (2) Two interim analyses and 1 final analysis with logrank at IA1, 
#'        max(logrank, fleming-harrington(0,1) at IA2, and 
#'        max(logrank, fleming-harrington(0,1), fleming-harrington(1,1)) at final): 
#'        fws2 = list(IA1 = list(lr), IA2=list(lr, fh01), FA=list(lr,fh01,fh11)).
#'        \item (3) One IA and one FA: stabilized Fleming-Harrington (0,1) at IA, 
#'        and max(logrank, stabilized Fleming-Harrington (0, 1)) at FA. 
#'        fw3 = list(IA = list(sfh01), FA=list(lr, sfh01)).
#'        \item General format of weighting strategy specification is: (a) the weight
#'        functions for each analysis must be provided in list() object even there is only
#'        1 weight function for that an analysis. (b) The commonly used functions 
#'        are directly available including lr: log-rank test; fh01: Fleming-Harrington (0,1) test;
#'        fh11: Fleming-Harrington(1,1) test; fh55: Fleming-Harrington(0.5, 0.5) test.
#'        stabilized versions at median survival time: sfh01, sfh11, sfh55. Modestly 
#'        log-rank test: mlr.
#'        (c) User-defined weight function can also be handled, but the weight function
#'        must be defined as a function of survival rate, i.e., fws = function(s){...}, 
#'        where s is the survival rate S(t-). 
#'        \item Options of weighting strategies for exploration: fws.options=list(fws1, fws2, fws3).
#'        }
#' @param overall.alpha One-sided overall type I error, default 0.025.       
#' @param type1err Incremental type I error. sum(type1err) = overall.alpha. When targetEvents is NULL, type1err is required.
#' @param H0 "Y" or "N" to indicate whether the simulation is for type I error
#' @param Cox "Y" or "N" to indicate whether Cox regression model is required to produce HR estimate
#' @param Median "Y" or "N" to indicate whether median is requested to estimate.
#' @param scanfreq0   A vector of scan frequency for control arm, eg, every 6 weeks for 2 scans, 
#' then every 9 weeks thereafter. scanfreq = c(rep(6, 2), rep(9, 1000))*7. Applicable for PFS data.
#' @param scanfreq1   A vector of scan frequency for experimental arm, eg, every 6 weeks for 2 scans, 
#' then every 9 weeks thereafter. scanfreq = c(rep(6, 2), rep(9, 1000))*7. Applicable for PFS data.
#' 
#' @return An object with a dataframe for each analysis including the following variables:
#' \describe{
#' \item{power}{Power for each analysis}
#' \item{overall.power}{Overall power of the group sequential design}
#' \item{wlr.simulations}{Simulation results for each simulated study data. 
#' An array with dimensions: nSim simulations X M testing strategies 
#' (fws.options) X K analyses X 5 variables:
#'      \itemize{ 
#'      \item z value
#'      \item p value
#'      \item analysis
#'      \item rejection boundary
#'      \item testing result (1 = positive; 0 = negative)
#'      }}
#' \item{lr.power}{Power for each analysis using log-rank test. Available if logrank ="Y"}
#' \item{lr.overall.power}{Overall power of the group sequential design using logrank test}
#' \item{lr.simulations}{Simulation results for each simulated study data using logrank test.
#'      nSim simulations X K analyses X 5 variables as above}
#' \item{CutOff}{A matrix with dimension (nSim, K). Date cutoff date for each analysis in each simulated trial.}
#' }
#' @examples
#' lr = function(s){1}
#' fh01 = function(s){1-s}
#' fh11 = function(s){s*(1-s)}
#' sfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1-s1)}
#' 
#' #Example (1): Simulate 10 samples from proportional hazards scenario. 
#' fws1 = list(IA1 = list(lr), FA = list(lr))
#' fws2 = list(IA1 = list(lr), FA = list(fh11))
#' fws3 = list(IA1 = list(lr), FA = list(lr, fh11))
#' 
#' 
#' #7 weighting strategies for exploration 
#' fws = list(fws1, fws2, fws3)
#' 
#' set.seed(2022)
#' 
#' #Simulations for exploring 3 weighting strategies
#' #medians 12.8 vs 17.9
#' m0 = qmcr(0.5, p = 0.12, alpha = log(2)/10, beta=1, gamma=1, lambda=0, tau=0, psi=1)
#' 
#' m1 = qmcr(0.5, p = 0.12, alpha = log(2)/10, beta=1, gamma=1, lambda=0, tau=6, psi=0.6)
#' 
#' sim = simulation.rgs.mcr(nSim=10, N = 600, A = 18, w=1.5, r=1, p=c(0.12, 0.12),
#' alpha = c(log(2)/10,log(2)/10), beta=c(1,1), gamma=c(1,1), 
#' lambda=c(0,0), tau=c(0,6), psi=c(1,0.6), drop=c(0,0),
#' targetEvents = c(300, 420), DCO = NULL,
#' sf = "LDOF", overall.alpha = 0.025, type1err = NULL,
#' logrank="N", fws.options=fws, H0 = "N")
#' 
#' sim = simulation.rgs.mcr(nSim=1, N = 500, A = 18, w=1.5, r=1, p=c(0.12, 0.12),
#' alpha = c(log(2)/10,log(2)/10), beta=c(1,1), gamma=c(1,1), 
#' lambda=c(0,0), tau=c(0,6), psi=c(1,0.6), drop=c(0,0),
#' targetEvents = c(300, 420), DCO = NULL,
#' sf = "LDOF", overall.alpha = 0.025, type1err = NULL,
#' Cox = "Y", scanfreq0=NULL, scanfreq1=NULL,
#' logrank="N", fws.options=fws, H0 = "N")
#' 
#' @export 
simulation.rgs.mcr = function(nSim=10, N = 600, A = 18, w=1.5, r=1, p=c(0.1,0.1), 
                              alpha = c(log(2)/12,log(2)/12), beta=c(1,1), gamma=c(1,1), 
                              lambda=c(0,0), tau=c(0, 6), psi=c(1, 0.6), drop=c(0,0),
                              targetEvents = NULL, DCO = NULL,
                              sf = "LDOF", overall.alpha = 0.025, type1err = NULL,
                              logrank="N", fws.options=NULL, H0 = "N", Cox = "N", Median="N",
                              scanfreq0=NULL, scanfreq1=NULL, seed=2022) {
  set.seed(seed)
  
  ##############################
  #Check if the simulation is for type I error
  ##############################
  if (H0 == "Y"){
    alpha = rep(alpha[1],2)
    beta = rep(beta[1],2)
    gamma = rep(gamma[1],2)
    lambda = rep(lambda[1], 2)
    tau = rep(tau[1],2)
    psi = rep(psi[1],2)
    drop = rep(drop[1], 2)
  } 
  
  ##############################
  #M Options of test strategies
  ##############################
  M = length(fws.options)
  
  ##############################
  #K analyses
  ##############################
  K = ifelse(is.null(DCO), length(targetEvents), length(DCO))
  if (is.null(type1err)*is.null(targetEvents)){print("Error: When targetEvents is NULL, type1err is required.")}
  if (is.null(type1err)){timing = targetEvents / targetEvents[K]}
  
  ##############################
  #Alpha spending
  ##############################
  side = 1 #always 2-sided

  #if alpha is not provided, use sf to derive alpha. 
  #if alpha is provided, then sf is ignored.
  if(is.null(type1err) && !is.null(overall.alpha)){
    ld.obf = function(s){
      if (side == 1){a = 2*(1 - pnorm(qnorm(1-overall.alpha/2)/sqrt(s)))}
      if (side == 2){a = 2*2*(1 - pnorm(qnorm(1-overall.alpha/4)/sqrt(s)))}
      return(a)
    }
    ld.pk = function(s){overall.alpha * log(1 + (exp(1)-1)*s)}
    
    if (K > 1){
      if (sf == "LDOF"){gs.alpha = ld.obf(s = timing)}
      if (sf == "LDPK") {gs.alpha = ld.pk(s = timing)}
      type1err[1] = gs.alpha[1]
      for(i in 2:K){type1err[i] = gs.alpha[i] - gs.alpha[i-1]}      
    } else {
      type1err = overall.alpha
    }
  }
  
  ##############################
  #Simulation
  ##############################
  if (M > 0) {wlr.sim = array(NA, dim=c(nSim, M, K, 5))}
  if (logrank == "Y") {lr.sim = array(NA, dim=c(nSim, K, 5))}
  CutOff = hr = matrix(NA, nrow=nSim, ncol = K)
  
  #(1). Generate data

  for (i in 1:nSim) {
    if (is.null(scanfreq0) || is.null(scanfreq1)) {
      dati = simulation.mcr(nSim=1, N = N, A = A, w=w, r=r, p=p, 
                           alpha = alpha, beta=beta, gamma=gamma, 
                           lambda=lambda, tau=tau, psi=psi, drop=drop,
                           targetEvents = targetEvents, DCO = DCO) 
    } else{
      dati = simulation.mcr.pfs(nSim=1, N = N, A = A, w=w, r=r, p=p, 
                               alpha = alpha, beta=beta, gamma=gamma, 
                               lambda=lambda, tau=tau, psi=psi, drop=drop,
                               targetEvents = targetEvents, DCO = DCO,
                               scanfreq0=scanfreq0, scanfreq1=scanfreq1)
    }

    #(2). Add group variable
    for (j in 1:K) {
      dati[[j]]$group = as.numeric(dati[[j]]$treatment == "experimental")
      
      dati[[j]]$calendarCutoff = as.numeric(dati[[j]]$calendarCutoff)
      dati[[j]]$survTimeCut = as.numeric(dati[[j]]$survTimeCut)
      dati[[j]]$cnsrCut = as.numeric(dati[[j]]$cnsrCut)
      CutOff[i,j] =dati[[j]]$calendarCutoff[1]
    }  
    
    #(3). Testing strategy m
    if (M > 0) {
      for (m in 1:M){    
        #Perform weighted log-rank test for each analysis in strategy m
        wlri = wlr.inference(data=dati, alpha = type1err, f.ws=fws.options[[m]])$test.results
        wlri = wlri[!duplicated(wlri$analysis),]
        wlri$result = as.numeric(wlri$inference=="Positive")
        wlr.sim[i, m, , ] = as.matrix(wlri[,c(2,3,5,6,8)])
      }
    }
    #(4). Standard log-rank test for all analyses if requested
    if (logrank=="Y"){
      if(K > 1){
        #GSD boundary for each analysis
        if(side == 1) {
          z.bd <- gsDesign::gsDesign(k=K,  alpha=overall.alpha,timing=timing[1:(K-1)], 
                                     sfu=gsDesign::sfLDOF)$upper$bound 
        } else{
          z.bd <- gsDesign::gsDesign(k=K,  alpha=overall.alpha/2,timing=timing[1:(K-1)], 
                                     sfu=gsDesign::sfLDOF)$upper$bound 
        }
      } else {
        if(side == 1) {z.bd = qnorm(1-overall.alpha)} else{z.bd = qnorm(1-overall.alpha/2)}
      }
      for (j in 1:K){
        lr.test = survival::survdiff(survival::Surv(as.numeric(survTimeCut), 1-as.numeric(cnsrCut)) ~ group, data = dati[[j]])
        
        #convert to z value in correct direction: z>0 means better experimental arm.
        better = as.numeric(lr.test$obs[1] > lr.test$obs[2])
        sign = 2*better - 1
        z = sqrt(lr.test$chisq) * sign
        
        #count power
        lr.sim[i, j, 1] = z
        if (side == 1){ p = 1 - pnorm(z)} else {p = 2*(1 - pnorm(z))}
        
        lr.sim[i, j, 2] = p
        lr.sim[i, j, 3] = j
        lr.sim[i, j, 4] = z.bd[j]
        lr.sim[i, j, 5] = as.numeric(z > z.bd[j])
      }
    }
    
    #(5)Cox proportional hazards model
    if (Cox == "Y"){
      for (j in 1:K){
        cox.fit = survival::coxph(survival::Surv(survTimeCut,1-cnsrCut) ~ group, data=dati[[j]])
        hr[i,j] = summary(cox.fit)$coefficients[2]
      }
    }
    
    #(6) Median
    if (Median == "Y"){
      for (j in 1:K){
        km = survival::survfit(survival::Surv(survTimeCut,1-cnsrCut) ~ group, data=dati[[j]])
        med0[i,j] = summary(km)$table[,7][1]
        med1[i,j] = summary(km)$table[,7][2]
      }
    }
    
  }
  
  o=list()
  ##############################
  #Power
  ##############################
  if (M > 0) {
    pow = matrix(NA, nrow=M, ncol=K)
    overall.pow = rep(0, M)
  
    for(m in 1:M){for(j in 1:K){pow[m, j] = sum(wlr.sim[,m,j,5])/nSim}}
    for(m in 1:M){for(i in 1:nSim){
      overall.pow[m] =  overall.pow[m] + as.numeric(sum(wlr.sim[i,m,,5])>0)
    }}
    overall.pow = overall.pow / nSim
    
    o$power = pow; o$overall.power = overall.pow
    o$wlr.simulations = wlr.sim
  }
  ##############################
  #Output
  ##############################
  if(logrank=="Y"){
    lr.pow = rep(NA, K)
    for (j in 1:K) {lr.pow[j] = sum(lr.sim[,j,5])/nSim}
    
    lr.overall.pow = 0
    for (i in 1:nSim){lr.overall.pow = lr.overall.pow + as.numeric(sum(lr.sim[i,,5])>0)}
    lr.overall.pow = lr.overall.pow / nSim
    
    o$lr.overall.power = lr.overall.pow
    o$lr.power = lr.pow
    o$lr.simulations = lr.sim
  }
  o$CutOff = CutOff
  if (Cox == "Y"){ o$hr = hr}
  if (Median == "Y") {o$med0 = med0; o$med1 = med1}
  
  return(o)
}

