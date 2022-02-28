#' Cut simulated PFS data according to the pre-specified scan intervals
#'
#' This function cuts the progression data according to the scan intervals. 
#' The observed PD date is determined by the immediate scan following the unobservable 
#' true PD date. The scan interval should be considered only for events after 
#' incorporating drop-off process. Data cutoff process, if needed, should follow
#' the scan interval consideration. The event status should not be changed after
#' considering the scan interval.
#'
#' @param PFS  Time to event variable
#' @param cnsr Censor status variable (0: event; 1: censored)
#' @param treatment Treatment arm code (0: control arm; 1: experimental arm)
#' @param scanfreq0   A vector of scan frequency for control arm, eg, every 6 weeks for 2 scans, 
#' then every 9 weeks thereafter. scanfreq = c(rep(6, 2), rep(9, 1000))*7
#' @param scanfreq1   A vector of scan frequency for experimental arm, eg, every 6 weeks for 2 scans, 
#' then every 9 weeks thereafter. scanfreq = c(rep(6, 2), rep(9, 1000))*7
#' 
#' @examples
#' #Example: Simulate 1 sample from proportional hazards scenario. 
#' #Total 600 pts, 1:1 randomization, control median PFS 8 mo; 
#' #HR = 0.6
#' 
#' time0 = rexp(300, rate=log(2)/12)
#' time1 = rexp(300, rate=log(2)/12*0.6)
#' time = c(time0, time1)
#' cnsr = rep(0, 600)
#' treatment = c(rep("control", 300), rep("experiment", 300))
#' s = sample(1:600)
#' data = data.frame(cbind(time, cnsr))
#' data$treatment = treatment
#' data = data[s,]
#' 
#' scanfreq0 = scanfreq1 = c(rep(6, 2), rep(9, 100))*7/(365.25/12)
#' km<-survival::survfit(survival::Surv(time, 1-cnsr)~treatment, data=data)
#' plot(km,xlab="Month Since Randomization",ylab="Survival",
#' main="True unobservable PFS", lty=1:2,xlim=c(0,36))
#' 
#' data$PFS = pfs_scan(time = data$time, cnsr = data$cnsr, 
#'            treatment = data$treatment, scanfreq0 = scanfreq0, 
#'            scanfreq1=scanfreq1)
#' 
#' km2<-survival::survfit(survival::Surv(PFS, 1-cnsr)~treatment, data=data)
#' plot(km2,xlab="Month Since Randomization",ylab="Survival",
#' main="Observed PFS", lty=1:2,xlim=c(0,36))
#' 
#' @export 
#' 
#' 
pfs_scan = function(time = data$survTimeCut, cnsr = data$cnsrCut, 
                treatment = data$treatment, 
                scanfreq0=c(rep(6, 2), rep(9, 100))*7/(365.25/12), 
                scanfreq1=c(rep(6, 2), rep(9, 100))*7/(365.25/12)) {
  if (is.null(scanfreq0) || is.null(scanfreq1)) {
    PFS = time
  } else {
    #Convert the scan frequency into scan intervals
    scanint0 = scanfreq0
    for (i in 1:length(scanfreq0)){
      scanint0[i] = sum(scanfreq0[1:i])
    }
    scanint1 = scanfreq1
    for (i in 1:length(scanfreq1)){
      scanint1[i] = sum(scanfreq1[1:i])
    }
  
    #Find the scan interval (L_s, R_s) such that T in (L_s, R_s)
    PFS = rep(NA, length(time))
    for (j in 1:length(time)){
      timej = time[j]
      if (cnsr[j] == 0){
        #Consider scan intervals ONLY for events
        if (treatment[j] == "control"){
          PFS[j] = scanint0[timej < scanint0][1]
        } else {
          PFS[j] = scanint1[timej < scanint1][1]
        }
      } else {
        PFS[j] = timej
      }
    }
  }
  return(PFS)
}