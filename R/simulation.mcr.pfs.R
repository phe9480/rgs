#' Simulation of PFS Data Following Mixture Cure Rate Distribution after incorporating  scan intervals
#' 
#' Consider the survival function of mixture cure rate (MCR) distribution: 
#' `S(t) = p + (1-p)S0(t)`, where `S0(t)` is a survival distribution 
#' for susceptible subject, i.e., S0(0)=1 and S0(t) -> 0 as t -> Infinity.
#' S0(t) can be any proper survival function. For generality of S0(t), consider 
#' the generalized modified Weibull (GMW) distribution with parameters 
#' (alpha, beta, gamma, lambda) Martinez et al(2013). 
#' 
#' alpha: scale parameter
#' beta and gamma: shape parameters
#' lambda: acceleration parameter
#' 
#' S0(t) = 1-(1-exp(-alpha*t^gamma*exp(lambda*t)))^beta
#' 
#' Special cases:
#' (1) Weibull dist: lambda = 0, beta = 1. Beware of the parameterization difference.
#' (2) Exponential dist: lambda = 0, beta = 1, gamma = 1. The shape parameter (hazard rate) is alpha.
#' (3) Rayleigh dist: lambda = 0, beta = 1, gamma = 2.
#' (4) Exponentiated Weibull dist (EW): lambda = 0
#' (5) Exponentiated exponential dist (EE): lambda = 0 and gamma = 1
#' (6) Generalized Rayleigh dist (GR): lambda = 0, gamma = 2
#' (7) Modified Weibull dist (MW): beta = 1
#' 
#' Let T be the survival time according to survival function S1(t) below.
#' Denote T's distribution as MCR(p, alpha, beta, gamma, lambda, tau, psi), 
#' where tau is the delayed effect and psi is the hazard ratio after delayed effect,
#' ie. proportional hazards to S(t) after delay tau.
#' `S(t) = p + (1-p)*S0(t)`
#' `S1(t) = S(t)I(t<tau) + S(tau)^(1-psi)*S(t)^psi`
#' In brief, `S0(t)` is the proper GMW distribution `(alpha, beta, gamma, lambda)`;
#' `S(t)` is MCR with additional cure rate parameter `p`;
#' In reference to `S(t)`, `S1(t)` is a delayed effect distribution and proportional hazards after delay.
#' 
#' @param  nSim Number of trials
#' @param  N Total number patients in two arms.
#' @param  A Total accrual period in months. 
#' @param  w Weight parameter in cumulative enrollment pattern. The 
#'           cumulative enrollment at month t is (t / A)^w, eg, 
#'           at month 6, the enrollment is N*(6/24)^2 = N/16 for 
#'           24 months planned accrual period.
#' @param  r Randomization ratio `r:1`, where r refers to the experimental arm, 
#' eg, `r=2` in 2:1 ratio
#' @param  p Cure rate parameter. When p = 0, it reduces to `GMW(alpha, beta, gamma, lambda)` 
#' distribution. Each of the following parameters is a vector of 2 components for 
#' two treatment groups (control, experimental arm): 
#' `p`, `alpha`, `beta`, `gamma`, `lambda`, `tau`, `psi`, `drop`.
#' @param  alpha Generalized modified Weibull (GMW) distribution parameters. `alpha > 0`
#' @param  beta Generalized modified Weibull (GMW) distribution parameters. `beta > 0`
#' @param  gamma Generalized modified Weibull (GMW) distribution parameters. `gamma >= 0` but `gamma` and `lambda` cannot be both 0.
#' @param  lambda Generalized modified Weibull (GMW) distribution parameters. `lambda >= 0` but `gamma` and `lambda`` cannot be both 0.
#' @param  tau  Threshold for delayed effect period. `tau = 0` reduces no delayed effect.
#' @param  psi  Hazard ratio after delayed effect. `psi = 1` reduces to the survival function without proportional hazards after delayed period (`tau`).
#' @param  drop Drop-off rate per unit time. For example, if 3% drop off in 1 year of followup, then `drop = 0.03/12`.
#' @param  targetEvents A vector of target events is used to determine DCOs. For example, 
#'              397 target events are used to determine IA DCO; and 496 events are used 
#'              to determine the FA cutoff.
#' @param DCO   A vector of data cut-off time in months, calculated from first subject in. 
#'              Default NULL. The date cut-off times will be determined by targetEvents.
#'              If provided, then the targetEvents will be ignored.
#' @param scanfreq0   A vector of scan frequency for control arm, eg, every 6 weeks for 2 scans, 
#' then every 9 weeks thereafter. scanfreq = c(rep(6, 2), rep(9, 1000))*7. Applicable for PFS data.
#' @param scanfreq1   A vector of scan frequency for experimental arm, eg, every 6 weeks for 2 scans, 
#' then every 9 weeks thereafter. scanfreq = c(rep(6, 2), rep(9, 1000))*7. Applicable for PFS data.
#' 
#' @return A dataframe for each analysis including the following variables:
#' \describe{
#' \item{sim}{sequence number of simulated dataset;}
#' \item{treatment}{treatment group with values of "control" and "experimental"}
#' \item{enterTime}{Time of randomization in calendar time}
#' \item{calendarTime}{the time when event/censoring occurred in calendar time}
#' \item{survTime}{Survival time for analysis, = calendarTime - enterTime}
#' \item{cnsr}{censor status (0=event; 1=censor) before administrative censoring due to data cut}
#' \item{calendarCutOff}{Data CutOff Time (DCO);}
#' \item{survTimeCut}{Survival time after cut}
#' \item{cnsrCut}{Censor status after cut}
#' }
#' 
#' @examples
#' #Example. Delayed effect 6 months, proportional hazards HR = 0.6
#' #Total 600 pts, 1:1 randomization, control median OS 12 mo; 
#' #HR = 0.65, enrollment 24 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 400 and 500 events respectively.
#' 
#' set.seed(2022)
#' data = simulation.mcr.pfs(nSim=1, N = 600, A = 18, w=1.5, r=1, p=c(0.1,0.1), 
#' alpha = c(log(2)/12,log(2)/12), beta=c(1,1), gamma=c(1,1), 
#' lambda=c(0,0), tau=c(0,6), psi=c(1,0.6), drop=c(0,0),
#' targetEvents = c(350, 500), DCO = NULL, scanfreq0=c(rep(6, 2), rep(9, 100))*7/(365.25/12), 
#' scanfreq1=c(rep(6, 2), rep(9, 100))*7/(365.25/12)) 
#' 
#' data.IA = data[[1]][data[[1]]$sim==1,]
#' data.FA = data[[2]][data[[2]]$sim==1,]
#' m0 = qmcr(0.5, p=0.1, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=0, psi=1)
#' m1 = qmcr(0.5, p=0.1, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=6, psi=0.6)
#' 
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-as.numeric(cnsrCut))~treatment, data=data.IA)
#' plot(km.IA, ylab="Survival", xlim=range(km.FA$time))
#' abline(h = c(0.1, 0.5), col="gray80", lty=2)
#' abline(v=c(6, m1, m0), col="gray80", lty=2)
#'  
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-as.numeric(cnsrCut))~treatment, data=data.FA)
#' plot(km.FA, ylab="Survival", xlim=range(km.FA$time)) 
#' abline(h = c(0.1, 0.5), col="gray80", lty=2)
#' abline(v=c(6, m1, m0), col="gray80", lty=2)
#' 
#' #No consideration of scan intervals
#' set.seed(2022)
#' data0 = simulation.mcr.pfs(nSim=1, N = 600, A = 18, w=1.5, r=1, p=c(0.1,0.1), 
#' alpha = c(log(2)/12,log(2)/12), beta=c(1,1), gamma=c(1,1), 
#' lambda=c(0,0), tau=c(0,6), psi=c(1,0.6), drop=c(0,0),
#' targetEvents = c(350, 500), DCO = NULL, scanfreq0=NULL, 
#' scanfreq1=NULL) 
#' 
#' data0.IA = data0[[1]][data0[[1]]$sim==1,]
#' data0.FA = data0[[2]][data0[[2]]$sim==1,]
#' m0 = qmcr(0.5, p=0.1, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=0, psi=1)
#' m1 = qmcr(0.5, p=0.1, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=6, psi=0.6)
#' 
#' km.IA0<-survival::survfit(survival::Surv(survTimeCut,1-as.numeric(cnsrCut))~treatment, data=data0.IA)
#' km.FA0<-survival::survfit(survival::Surv(survTimeCut,1-as.numeric(cnsrCut))~treatment, data=data0.FA)
#' plot(km.IA0, ylab="Survival", xlim=range(km.FA0$time))
#' abline(h = c(0.1, 0.5), col="gray80", lty=2)
#' abline(v=c(6, m1, m0), col="gray80", lty=2)
#'  
#' plot(km.FA0, ylab="Survival", xlim=range(km.FA0$time)) 
#' abline(h = c(0.1, 0.5), col="gray80", lty=2)
#' abline(v=c(6, m1, m0), col="gray80", lty=2)
#' 
#' @export
#' 
simulation.mcr.pfs = function(nSim=100, N = 600, A = 18, w=1.5, r=1, p=c(0.1,0.1), 
                  alpha = c(log(2)/12,log(2)/12), beta=c(1,1), gamma=c(1,1), 
                  lambda=c(0,0), tau=c(0,0), psi=c(1,1), drop=c(0,0),
                  targetEvents = c(300, 420), DCO = NULL, 
                  scanfreq0=c(rep(6, 2), rep(9, 100))*7/(365.25/12), 
                  scanfreq1=c(rep(6, 2), rep(9, 100))*7/(365.25/12)){

  ###################
  #Number of analyses
  ###################
  if(is.null(DCO)){ L = length(targetEvents)} else{L = length(DCO)}
  
  ###########################
  #Recruitment data utility
  ###########################
  f.nEachMonth = function (N=600, A=18, w=1.5, r=1){
    #f.nEachMonth(N=600, A=18, w=1.5, r=1)
    N1 = N * (r/(r+1)); N0 = N - N1
    #When r > 1, the control arm has smaller number of pts.
    #Just need to determine enrollment for control arm per month,
    #then to obtain enrollment for experimental arm by n1i = n0i * r.
    n1 = n0 = rep(NA, A) #enrollment by month
    randdt1 = rep(NA, N1) #randomization date
    randdt0 = rep(NA, N0)
    #Determine number of pts per month for control arm
    #(i-1)th month cumulative enrolled pts
    cLastN0 = 0
    for (i in 1:A) {
      #ith month: cumulative #pts
      cN0i = max(round((i/A)^w * N0), 1)
      n0[i] = max(cN0i - cLastN0, 1)
      if (i == A) {n0[i] = N0 - sum(n0[1:(A-1)]) }
      cLastN0 = cN0i  
    }
    n1 = n0 * r
    
    #Patch for extreme rare scenarios that 0 enrollment in the last month
    if(n0[A] == 0 && n0[A-1] > 1){n0[A-1] = n0[A-1]-1; n0[A]=1}
    if(n1[A] == 0 && n1[A-1] > 1){n1[A-1] = n1[A-1]-1; n1[A]=1}
    
    o = list()
    o$n0 = n0; o$n1 = n1
    return(o)
  }
  ###########################
  #Data cut utility
  ###########################
  f.dataCut = function(data, targetEvents = 397, DCO = NULL) {
    data0 = data
    data0.order <- data0[order(as.numeric(data0$calendarTime)), ] #order by calendar time
    data.event <- data0.order[data0.order$cnsr == 0,] #Events Only
    
    data.event$event.seq <- seq.int(nrow(data.event)) #event sequence number
    if(is.null(DCO)){
      #Data cutoff in calendar time added to the original dataframe as a variable
      data0$calendarCutoff = as.numeric(data.event$calendarTime[data.event$event.seq == targetEvents])
    } else {
      data0$calendarCutoff = DCO
    }
    data0$survTimeCut = ifelse(as.numeric(data0$calendarTime) <= as.numeric(data0$calendarCutoff), as.numeric(data0$survTime), as.numeric(data0$calendarCutoff) - as.numeric(data0$enterTime))
    data0$cnsrCut = ifelse(as.numeric(data0$calendarTime) <= as.numeric(data0$calendarCutoff), data0$cnsr, 1)
    
    return(data0)
  }
  
  out = list(NULL)
  if (L > 1){for (k in 2:L){out = c(out, list(NULL))}}
  
  for (i in 1:nSim){
    nEachMonth = f.nEachMonth(N=N, A=A, w=w, r=r)
    n0 = sum(nEachMonth$n0); n1 = sum(nEachMonth$n1)
    
    ########################
    #MCR data for each arm
    ########################
    T0 = rmcr(n=n0, p=p[1], alpha = alpha[1], beta=beta[1], 
              gamma=gamma[1], lambda=lambda[1], tau=tau[1], psi=psi[1])
    T1 = rmcr(n=n0, p=p[2], alpha = alpha[2], beta=beta[2], 
              gamma=gamma[2], lambda=lambda[2], tau=tau[2], psi=psi[2])
    
    #Permutation of the original ordered samples
    T0 = sample(T0); T1 = sample(T1)
    
    ############################
    #Drop Off data for each arm
    ############################
    if (drop[1] > 0) {W0 = rexp(n0, rate=drop[1])} else {W0 = rep(Inf, n0)}
    if (drop[2] > 0) {W1 = rexp(n1, rate=drop[2])} else {W1 = rep(Inf, n1)}
    
    ############################
    #Censor data from Drop-off
    ############################
    Y0 = apply(cbind(T0, W0), 1, min)
    Y1 = apply(cbind(T1, W1), 1, min)
    
    event0 = as.numeric(T0 < Inf)
    event0[W0 < T0] = 0
    event1 = as.numeric(T1 < Inf)
    event1[W1 < T1] = 0
    
    ############################
    #EnterTime, CalendarTime
    ############################
    enterTime0 = rep(NA, n0)
    enterTime0[1:nEachMonth$n0[1]] = runif(nEachMonth$n0[1], min=0, max=1)
    if (A > 1) {for (m in 2:A){
      LL = sum(nEachMonth$n0[1:(m-1)])+1
      UU = sum(nEachMonth$n0[1:m])
      enterTime0[LL:UU] = runif(nEachMonth$n0[m], min=m-1, max=m)
    }}
    
    enterTime1 = rep(NA, n1)
    enterTime1[1:nEachMonth$n1[1]] = runif(nEachMonth$n1[1], min=0, max=1)
    if (A > 1) {for (m in 2:A){
      LL = sum(nEachMonth$n1[1:(m-1)])+1
      UU = sum(nEachMonth$n1[1:m])
      enterTime1[LL:UU] = runif(nEachMonth$n1[m], min=m-1, max=m)
    }}
    
    ############################
    #Assemble data before cut
    ############################
    sim = rep(i, N)
    treatment = c(rep("control", n0), rep("experimental", n1))
    enterTime = c(enterTime0, enterTime1)
    survTime = as.numeric(c(Y0, Y1))
    
    #trick infinity
    survTime[survTime > 1e6] = 1e6
    calendarTime = as.numeric(enterTime) + as.numeric(survTime)
    cnsr = c(1-event0, 1-event1)
    
    #PFS scan interval
    survTime = pfs_scan(time = survTime, cnsr = cnsr, 
                   treatment = treatment, 
                   scanfreq0=scanfreq0, scanfreq1=scanfreq1)
    
    dati = data.frame(cbind(sim, treatment, enterTime, calendarTime, survTime, cnsr))
    
    ############################
    #Cut data
    ############################
    dati.cut = NULL
    for (ii in 1:L){
      dati.cut[[ii]] = f.dataCut(data=dati, targetEvents=targetEvents[ii], DCO = DCO[ii])
      out[[ii]] = rbind(out[[ii]], dati.cut[[ii]])
    }
  }
  return(out)
}

