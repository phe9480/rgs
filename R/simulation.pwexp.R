#' Simulate data following piece-wise exponential distribution with non-uniform accrual pattern and lost to follow-up
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
#'
#' @param  nSim Number of trials
#' @param  N Total number patients in two arms.
#' @param  A Total accrual period in months
#' @param  w Weight parameter in cumulative enrollment pattern. 
#' The cumulative enrollment at month t is (t / A)^w, eg, at month 6, 
#'   the enrollment is N*(6/24)^2 = N/16 for 24 months planned accrual period.
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' @param  lam0 Hazard rates for control arm of intervals defined by cuts; for exponential(lam0) distribution,
#'         lam0 = log(2) / median;
#' @param  lam1 Hazard rates for experimental arm for intervals; for exponential(lam1) distribution,
#'         lam1 = log(2) / median; For delayed effect under H1, lam1 is a vector (below).
#' @param  cuts Timepoints to form intervals for piecewise exponential distribution. For example,
#   \itemize{
#   \item Proportional hazards with hr = 0.65. Then lam0 = log(2)/m0, lam1 = log(2)/m0*hr, cuts = NULL. 
#   \item Delayed effect at month 6, and control arm has constant hazard (median m0) and 
#'       experimental arm has hr = 0.6 after delay, then cuts = 6, and 
#'       lam0 = log(2) / m0 or lam0 = rep(log(2) / m0, 2), 
#'       lam1 = c(log(2)/m0, log(2)/m0*hr). 
#   \item Delayed effect at month 6, and control arm has crossover to subsequent IO 
#'       treatment after 24 mo, so its hazard decreases 20%. Then, 
#'       lam0 = c(log(2)/m0, log(2)/m0, log(2)/m0*0.8), 
#'       lam1 = c(log(2)/m0, log(2)/m0*hr, log(2)/m0*hr), and
#'       cuts = c(6, 24), which forms 3 intervals (0, 6), (6, 24), (24, infinity)
#       }
#' @param drop0 Drop Off rate per month, eg, 3% every year for control arm, then drop0=0.03/12
#' @param drop1 Drop Off rate per month, eg, 1%, for experimental arm
#' @param targetEvents A vector of target events is used to determine DCOs. For example, 
#'              397 target events are used to determine IA DCO; and 496 events are used 
#'              to determine the FA cutoff.
#' @param DCO   A vector of data cut-off time in months, calculated from first subject in. 
#'              Default NULL. The date cut-off times will be determined by targetEvents.
#'              If provided, then the targetEvents will be ignored.
#' @return An object with a dataframe for each analysis including the following variables:
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
#' @examples
#' #Example (1): Simulate 10 samples from proportional hazards scenario. 
#' #Total 600 pts, 1:1 randomization, control median OS 12 mo; 
#' #HR = 0.65, enrollment 24 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 400 and 500 events respectively.
#' 
#' sim.ph = simulation.pwexp(nSim=10, N = 600, A = 24, w=1.5, r=1, lam0=log(2)/12, lam1= log(2)/12*0.65, cuts=NULL, drop0= 0, drop1= 0, targetEvents = c(400, 500))
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.ph[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.ph[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' #Example (2): Simulate 10 samples with delayed effect at month 6;
#' #Total 600 pts, 1:1 randomization, control median OS 12 mo; 
#' #HR = 0.65, enrollment 24 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 400 and 500 events respectively.
#' 
#' sim.delay6 = simulation.pwexp(nSim=10, N = 600, A = 24, w=1.5, r=1, lam0=log(2)/12, lam1= c(log(2)/12,log(2)/12*0.65), cuts=6, drop0= 0, drop1= 0, targetEvents = c(400, 500))
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.delay6[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data= sim.delay6[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' #Example (3): Simulate 10 samples with delayed effect at month 6 
#' #Control arm has crossover to subsequent IO after 24 mo, so its hazard decreases 20%.
#' #control arm has constant hazard (median 11.7 mo) and experimental arm has 
#' #hr = 1 and 0.65 at intervals (0, 6) and (6, 24) respectively.
#' #HR = 0.65, enrollment 24 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 400 and 500 events respectively.
#' 
#' crossEffect = 0.8 #Hazard ratio in control arm (after crossover vs before crossover)
#' lam0 = log(2)/12*c(1, 1, crossEffect); lam1 = log(2)/12*c(1, hr, hr)
#' sim.delay6crs=simulation.pwexp(nSim=10,N=600,A=24,w=1.5,r=1,lam0=lam0, lam1=lam1,cuts=c(6, 24),drop0=0,drop1=0, targetEvents=c(400, 500))
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.delay6crs[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data= sim.delay6crs[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' @export 
simulation.pwexp = function(nSim=100, N = 600, A = 21, w=1.5, r=1, 
                            lam0=log(2)/12, lam1=log(2)/12*0.65, cuts=NULL, 
                            drop0=0, drop1=0, targetEvents = c(400, 500), 
                            DCO = NULL) {
  
  f.nEachMonth = function (N=600, A=24, w=2, r=2) {
    
    N1 = N * (r/(r+1))
    N0 = N - N1
    
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
    o$n0 = n0
    o$n1 = n1
    return(o)
  }
  f.dataCut = function(data, targetEvents = 397, DCO = NULL) {
    data0 = data
    data0.order <- data0[order(data0$calendarTime), ] #order by calendar time
    data.event <- data0.order[data0.order$cnsr == 0,] #Events Only
    
    data.event$event.seq <- seq.int(nrow(data.event)) #event sequence number
    if(is.null(DCO)){
      #Data cutoff in calendar time added to the original dataframe as a variable
      data0$calendarCutoff = data.event$calendarTime[data.event$event.seq == targetEvents] 
    } else {
      data0$calendarCutoff = DCO
    }
    data0$survTimeCut = ifelse(data0$calendarTime <= data0$calendarCutoff, data0$survTime, data0$calendarCutoff-data0$enterTime)
    data0$cnsrCut = ifelse(data0$calendarTime <= data0$calendarCutoff, data0$cnsr, 1)
    
    return(data0)
  }
  
  nEachMonth = f.nEachMonth(N=N, A=A, w=w, r=r)
  
  gamma = nEachMonth$n0 + nEachMonth$n1
  eta0 = -log(1-drop0)
  eta1 = -log(1-drop1)
  
  o = nphsim::nphsim(nsim=nSim,lambdaC=lam0,lambdaE=lam1, ssC=N/(r+1), ssE=N*r/(r+1),
             intervals=cuts, gamma=gamma, R=rep(1, A),eta=eta0, etaE=eta1,fixEnrollTime = FALSE)
  dat = o$simd
  data.out = NULL
  #number of analyses
  if(is.null(DCO)){
    L = length(targetEvents)
  } else{
    L = length(DCO)
  }
  
  for (j in 1:nSim) {
    
    dataj = dat[dat$sim == j,]
    
    #rename variables
    dataj <- dplyr::rename(dataj, enterTime = enterT, calendarTime = ct, survTime= survival)

    #cut data according to the specified target events    
    dataj.cut <- lapply(1:L, function(k) {
      f.dataCut(data=dataj, targetEvents[k], DCO = DCO[k])
    })
    
    data.out=lapply(1:L, function(k) {
      rbind(data.out[[k]], dataj.cut[[k]])
    })
    
  }
  return(data.out)
}