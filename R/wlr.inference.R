#'  Calculation of Actual Rejection Boundary and Statistical Inference for Group 
#'  Sequential Design Using Weighted Log-rank Tests
#' 
#'  This function calculates the actual rejection boundary for the next analysis 
#'  using weighted log-rank test based on the datasets from previous analyses. 
#'  In standard log-rank test, the rejection boundaries can be determined based 
#'  on the number of events and the alpha spending function, because the asymptotic
#'  distribution can be approximated by the number of events. However, when using weighted 
#'  log-rank test, the asymptotic distribution is usually associated with the 
#'  pooled survival curve estimated from the actual data, eg, Fleming-Harrington class.
#'  As a result, the actual rejection boundary depends on the asymptotic correlation 
#'  estimated from the actual data. Refer to Tsiatis (1982) for the consistent 
#'  estimator of the asymptotic correlation matrix.
#'  
#' @param data  The dataframe used for the analyses including all data in previous analyses.
#'              For example, for the first analysis, data = list(IA1 = data1). 
#'              For the 2nd analysis, data = list(IA1 = data1, IA2=data2), where
#'              data1 is the data used for the first analysis. If IA2 is the final
#'              analysis, data = list(IA1 = data1, FA=data2). For the 3rd analysis,
#'              data = list(IA1 = data1, IA2=data2, IA3=data3) or 
#'              data = list(IA1 = data1, IA2=data2, FA=data3).
#'              Every dataset must be sorted by subject id and is required have 
#'              a format of 1 subject 1 record. In the current version, every 
#'              dataset must have the same subjects. Furthermore, each dataset 
#'              must include the following variables: (1) survTimeCut - Survival time 
#'              for each analysis; (2) cnsrCut - Censoring status for each analysis 
#'              (0 = event; 1 = censored); (3) group - Treatment group (0 = control, 
#'              1 = experimental treatment); (4) strata1, strata2, strata3: 
#'              Stratification variables. Up to 3 stratification factors are allowed
#'              in the current version.                  
#' @param alpha Incremental type I errors (1-sided) allocated to the analyses in group 
#' sequential tests. Must be 1-sided type I error. The sum of alpha is the overall
#' alpha allocated to the hypothesis test. For example, IA has type I error 0.01,
#' and the overall type I error is 0.025, then alpha = c(0.01, 0.015). If alpha is
#' specified, sf is ignored. For another example, if the current analysis is IA2 
#' and the cumulative type I error spending up to IA2 is 0.02, and IA1 has type I error
#' spending of 0.005, then alpha = c(0.005, 0.015). 
#' alpha_i = P(accept H0 for first i-1 analyses and reject H0 at ith analysis | H0). 
#' If the cumulative alpha spending by ith analysis is ai, then alpha_i = ai - a_{i-1}.
#' @param  f.ws  Weight functions in the weighted logrank tests for all previous analyses
#'               and the current analysis. For example, IA1 uses standard log-rank test, 
#'               IA2 uses a maxcombo test of (logrank, FH01), and FA uses a maxcombo 
#'               test of (logrank, FH01, FH11). Then f.ws = list(IA1=list(lr), 
#'               IA2=list(lr, fh01), FA=list(lr, fh01, fh11)), where define 
#'               lr = function(s){1}; fh01=function(s){1-s}; fh11 = function(s){s*(1-s)}.
#'               For another example, if only fh01 is used for all three analyses, 
#'               then f.ws = list(IA1=list(fh01), IA2=list(fh01), FA1=list(fh01)).
#'
#' @return 
#'  \itemize{
#'  \item  test.results:  Including actual rejection boundary and statistical inference.
#'  \item  corr: Correlation matrix among the weighted log-rank tests at different analyses.
#'  \item  wt: Weight functions used for each analysis
#'  \item  alpha: Incremental type I error used for each analysis
#'  }
#'  
#' @examples 
#' #Simulate a dataset: 600 subjects, control arm median 12 months following exponential
#' distribution, 1:1 randomization, recruitment period 24 months with weight parameter 1.5.
#' Experimental arm has delayed effect of 6 months. The HR after delay is 0.65.
#' Three analyses are performed at 300, 400, and 500 events.
#' data0 = simulation.pwexp(nSim=1, N = 600, A = 24, w=1.5, r=1, lam0=log(2)/12, 
#'       lam1=c(log(2)/12, log(2)/12*0.65), cuts=6, drop0= 0, drop1= 0, 
#'       targetEvents = c(300, 400, 500))
#' data1 = data0[[1]][sim==1,]; data2 = data0[[2]][sim==1,]; data3 = data0[[3]][sim==1,]
#' #Add strata variables       
#' data1$strata1 = data1$strata2 =data1$strata3 =sample(c(1,2), 600, replace = TRUE);
#' data2$strata1 = data2$strata2 =data2$strata3 =sample(c(1,2), 600, replace = TRUE);
#' data3$strata1 = data3$strata2 =data3$strata3 =sample(c(1,2), 600, replace = TRUE);
#' data1$group = as.numeric(data1$treatment == "experimental")
#' data2$group = as.numeric(data2$treatment == "experimental")
#' data3$group = as.numeric(data3$treatment == "experimental")
#' 
#' km.IA1<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=data1)
#' plot(km.IA1,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.IA2<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=data2)
#' plot(km.IA2,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=data3)
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' #Define weight functions for weighted log-rank tests
#' lr = function(s){1}; fh01 = function(s){(1-s)}; fh11 = function(s){s*(1-s)}
#' sfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1-s1)} 
#' 
#' #Example (1). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and FH01);
#' #               FA uses max(log-rank, FH01, FH11).
#' 
#' wlr.inference(data=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.005, 0.01, 0.01), 
#'     strata1 = data1$strata1, strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), IA2=list(lr, fh01), FA=list(lr, fh01, fh11)))
#'   
#' #Example (2a). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and sFH01);
#' #             FA uses max(log-rank, FH01). Unstratified analysis
#'               
#' wlr.inference(data=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.005, 0.01, 0.01), strata1 = NULL, strata2 = NULL, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), IA2=list(lr, sfh01), FA=list(lr, fh01)))
#'          
#' #Example (2b). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and sFH01);
#' #              FA uses max(log-rank, FH01). 
#' #              Stratified analysis of 2 stratification factors
#'               
#' wlr.inference(data=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.005, 0.01, 0.01), strata1 = data1$strata1, 
#'     strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), IA2=list(lr, sfh01), FA=list(lr, fh01)))
#'          
#' #Example (2c). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and sFH01);
#' #              FA uses max(log-rank, FH01). 
#' #              Stratified analysis of 3 stratification factors
#'               
#' wlr.inference(data=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.005, 0.01, 0.01),strata1 = data1$strata1, 
#'     strata2 = data1$strata2, strata3 = data1$strata3,
#'     f.ws=list(IA1=list(lr), IA2=list(lr, sfh01), FA=list(lr, fh01)))
#'          
#' #Example (3). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and sFH01);
#' #               FA uses max(log-rank, FH01). 
#'               
#' wlr.inference(data=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.005, 0.01, 0.01), strata1 = data1$strata1, 
#'     strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), IA2=list(lr, sfh01), FA=list(lr, fh01)))
#'          
#' #Example (4). 2 IAs and FA. All analyses use max(logrank, sfh01).
#'               
#' wlr.inference(data=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.005, 0.01, 0.01), strata1 = data1$strata1, 
#'     strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr, sfh01), IA2=list(lr, sfh01), FA=list(lr, sfh01)))
#'          
#'  #For each analysis, the wlr.maxcombo() function can perform
#'  #stratified maxcombo test and produce the equivalent p value.
#'  #However, wlr.maxcombo() function cannot provide statistical inference
#'  #without considering the actual rejection boundary.
#'          
#' wlr.maxcombo(time=data1$survTimeCut, event=1-data1$cnsrCut, group=data1$group, 
#'  rho=NULL, gamma=NULL, tau = NULL, s.tau=NULL,
#'  strata1 = data1$strata1, strata2=data1$strata2, strata3=data1$strata3,
#'  f.ws=list(lr, sfh01), side = 1)
#' 
#' #Example (5). 2 IAs and FA. All analyses use logrank.
#'               
#' wlr.inference(data=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.005, 0.01, 0.01), strata1 = data1$strata1, 
#'     strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), IA2=list(lr), FA=list(lr)))
#' @export
#' 
#' 
#'  
wlr.inference = function(data=list(IA1=data1, IA2=data2, FA=data3), 
        alpha = c(0.005, 0.01, 0.01), 
        strata1 = NULL, strata2 = NULL, strata3 = NULL,
        f.ws=list(IA1=list(lr), IA2=list(lr, fh01), FA=list(lr, fh01, fh11))){

  #K analyses
  K=length(f.ws)
  
  #Number of test components in each analysis
  J = lengths(f.ws)
  
  #Find the complete correlation matrix among all weighted log-rank 
  #tests J[1]+J[2]+...+J[K] dimensions
  corr = matrix(1, nrow=sum(J), ncol=sum(J))
  #calculate the correlation between Zij and Z_i'j'
  for (i in 1:K){
    for (j in 1:J[i]){
      for (ip in i:K){
        for (jp in 1:J[ip]){
          row = as.numeric(i>=2)*sum(J[1:(i-1)])+j #row location of the corr matrix
          #incremental location pointer for column compared to row of the corr matrix
          incr = (as.numeric(ip>=2)*sum(J[1:(ip-1)])+jp)-(as.numeric(i>=2)*sum(J[1:(i-1)])+j)
          col = row + incr
          #incr controls the computation only limited to upper right corner
          if(incr > 0){
            #information matrix for Zij and Zi'j'
            datai = data[[i]]; dataip = data[[ip]]
            corr[row, col] = wlr.cov2t(time1=as.numeric(datai$survTimeCut), event1=1-as.numeric(datai$cnsrCut), 
                            time2=as.numeric(dataip$survTimeCut), event2=1-as.numeric(dataip$cnsrCut), 
                            group=datai$group, 
                            strata1=strata1, strata2=strata2, strata3=strata3,
                            f.ws1=f.ws[[i]][[j]], f.ws2=f.ws[[ip]][[jp]])$corr
            corr[col, row] = corr[row, col]
          }
        }
      }
    }
  }
  #Rejection boundary: recursively solve for the rejection boundary for each analysis
  b = rep(NA, K)
  
  #First Analysis
  if (J[1] >= 2){
      #1-sided test
      f.b = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
            mean=rep(0, J[1]), corr = corr[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }
      b[1] = uniroot(f=f.b, interval=c(1, 20), tol = 1e-8)$root
  } else{
    b[1] = qnorm(1-alpha[1])
  }
  #Recursively solve other boundaries from 2nd analysis to Kth analysis
  if(K > 1){
    for(i in 2:K){
        f.b = function(x){
          LL1 = rep(-Inf, sum(J[1:(i-1)]))
          UU1 = rep(b[1], J[1])
          if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b[j], J[j]))}}
          idx1 = 1:sum(J[1:(i-1)])
          if (length(LL1) == 1) {
            P1 = pnorm(UU1)
          } else {
            P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, mean=rep(0, length(idx1)), 
                                  corr = corr[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
          }
          
          LL2 = rep(-Inf, sum(J[1:i]))
          UU2 = c(UU1, rep(x, J[i]))
          
          idx = 1:sum(J[1:i]) #index of the corr matrix
          P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, mean=rep(0, length(idx)), 
                                corr = corr[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(P1 - P2 - alpha[i])
        }
        b[i] = uniroot(f=f.b, interval=c(1, 20), tol = 1e-8)$root
    }
  }
  #Calculate p value for each analysis
  test.results = NULL
  for (i in 1:K){
    testi = wlr.maxcombo(time=data[[i]]$survTimeCut, event=1-data[[i]]$cnsrCut,
            group=data[[i]]$group, 
            strata1=strata1, strata2=strata2, strata3=strata3, 
            rho = NULL, gamma=NULL, tau = NULL, s.tau=NULL,
            f.ws=f.ws[[i]], side = "one.sided")$test.results
    testi$analysis = i; testi$z.bound = b[i]; 
    
    testi$inference = ifelse(testi$z.max > b[i], "Positive", "Negative")
    test.results = rbind(test.results, testi)
  }
  
  o = list()
  o$corr = corr
  o$test.results = test.results 
  o$wt = f.ws
  o$alpha = alpha
  return(o)
}
