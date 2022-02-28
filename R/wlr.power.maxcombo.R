#'  Power Calculation for Group Sequential Design Using A Max-combo Test
#' 
#' @param n    Total sample size for two arms.
#' @param DCO  Analysis time, calculated from first subject in.
#' @param r  Randomization ratio of experimental arm : control arm as r:1. 
#'           When r = 1, it is equal allocation. Default r = 1.
#' @param alpha Allocated one-sided alpha levels. sum(alpha) is the total type I error.
#'           If alpha spending function a(t) is used for information time c(t1, ..., tK),
#'           then alpha1 = a(t1), alpha2 = a(t2)-a(t1), ..., alphaK = a(tK)-a(t_{K-1}),
#'           and the total alpha for all analyses is a(tK).
#' @param b  Rejection boundary in normalized Z. If b is NULL, then the boundaries 
#'           will be calculated based on alpha at each analysis time T. Default b = NULL.
#'           If b is provided, then alpha is ignored. Default NULL.
#' @param h0 Hazard function of control arm. h0(t) = log(2)/m0 means T~exponential distribution with median m0.
#' @param S0 Survival function of control arm. In general, S0(t) = exp(- integral of h0(u) for u from 0 to t).
#'           but providing S0(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f0(t) = h0(t) * S0(t).
#' @param h1 Hazard function of experimental arm. h1(t) = log(2)/m1 means T~exponential distribution with median m0.
#' @param S1 Survival function of experimental arm. In general, S1(t) = exp(- integral of h1(u) for u from 0 to t).
#'           but providing S1(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f1(t) = h1(t) * S1(t).
#' @param cuts A vector of cut points to define piecewise distributions. 
#' If cuts is not specified or incorrectly specified, it might occasionally have numerical integration issue.
#' @param  rho Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#' @param  gamma Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#'         For log-rank test, set rho = gamma = 0.
#' @param  tau  Cut point for stabilized FH test, sFH(rho, gamma, tau); with weight
#'       function defined as w(t) = s_tilda^rho*(1-s_tilda)^gamma, where
#'       s_tilda = max(s(t), s.tau) or max(s(t), s(tau)) if s.tau = NULL
#'       tau = Inf reduces to regular Fleming-Harrington test(rho, gamma)
#' @param  s.tau  Survival rate cut S(tau) at t = tau1; default 0.5, ie. cut at median.
#'         s.tau = 0 reduces to regular Fleming-Harrington test(rho, gamma)
#' @param f.ws  Self-defined weight function of survival rate, eg, f.ws = function(s){1/max(s, 0.25)}
#'         When f.ws is specified, sFH parameter will be ignored.
#' @param Lambda Cumulative distribution function of enrollment. 
#' @param G0 Cumulative distribution function of drop-off for control arm, eg, G.ltfu=function(t){1-exp(-0.03/12*t)}
#'               is the distribution function for 3 percent drop-off in 12 months of followup.
#' @param G1 Cumulative distribution function of drop-off for experimental arm, eg, G.ltfu=function(t){1-exp(-0.03/12*t)}
#'               is the distribution function for 3 percent drop-off in 12 months of followup.
#' @param mu.method Method for mean of weighted logrank test Z statistic. 
#' "Schoenfeld" or "H1"
#' @param cov.method Method for covariance matrix calculation in power calculation: 
#'              "H0", "H1", "H1.LA" for null hypothesis H0, H1, H1 under local alternative
#' @return An object with dataframes below.
#'  \itemize{
#'  \item  n:      Total number of subjects for two arms
#'  \item  DCO:    Expected analysis time
#'  \item  targetEvents: Expected number of events
#'  \item  power:  Power of the max-combo test at each analysis
#'  \item  overall.power Overall power of the study
#'  \item  incr.power Incremental power for each analysis. 
#'                  The sum of all incremental powers is the overall power.
#'  \item  medians: Median of each treatment group
#'  \item  b:    Expected rejection boundary in z value
#'  \item  Expected_HR: Expected HR
#'  }
#'  Omega0: Covariance matrix under H0; 
#'  Omega1: Covariance matrix for power calculation request in cov.method.
#'  
#' @examples 
#' #Distributions for both arms
#' m0 = 12; #median OS for control arm
#' lambda0 = log(2) / m0
#' h0 = function(t){lambda0}; 
#' S0 = function(t){exp(-lambda0 * t)}
#' HRd = 0.60 #hazard ratio after delay
#' 
#' h.D3=function(t){lambda0*as.numeric(t<3)+HRd*lambda0*as.numeric(t>=3)}
#' c3 = exp(-3*lambda0*(1-HRd)); 
#' S.D3 = function(t){S0(t)*as.numeric(t<3)+c3*exp(-HRd*lambda0*t)*as.numeric(t>=3)}
#' 
#' #Define weight functions for weighted log-rank tests
#' lr = function(s){1}
#' fh01 = function(s){(1-s)}
#' fh11 = function(s){s*(1-s)}
#' 
#' #Enrollment
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)};
#' G0 = function(t){0}; G1 = function(t){0}; 
#' 
#' #Schoenfeld method with power based on covariance matrix under H0
#' wlr.power.maxcombo(DCO = c(24, 36),  
#'   alpha=c(0.01, 0.04)/2, 
#'   r = 1, n = 500, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, 
#'   f.ws = list(IA1 = list(lr), FA=list(fh01)), 
#'   Lambda=Lambda, G0=G0, G1=G1,
#'   mu.method = "Schoenfeld", cov.method = "H0")
#'   
#' #Schoenfeld method with power based on covariance matrix under H1
#' wlr.power.maxcombo(DCO = c(24, 36),  
#'   alpha=c(0.01, 0.04)/2, 
#'   r = 1, n = 500, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, 
#'   f.ws = list(IA1 = list(lr), FA=list(fh01)), 
#'   Lambda=Lambda, G0=G0, G1=G1,
#'   mu.method = "Schoenfeld", cov.method = "H1")
#'   
#' #Schoenfeld method with power based on covariance matrix under H1 in Local Alternative (simplified)
#' wlr.power.maxcombo(DCO = c(24, 36),  
#'   alpha=c(0.01, 0.04)/2, 
#'   r = 1, n = 500, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, 
#'   f.ws = list(IA1 = list(lr), FA=list(fh01)), 
#'   Lambda=Lambda, G0=G0, G1=G1,
#'   mu.method = "Schoenfeld", cov.method = "H1.LA")   
#'      
#' #Mean(Z) under H1 with power based on covariance matrix under H0
#' wlr.power.maxcombo(DCO = c(24, 36),  
#'   alpha=c(0.01, 0.04)/2, 
#'   r = 1, n = 500, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, 
#'   f.ws = list(IA1 = list(lr), FA=list(fh01)), 
#'   Lambda=Lambda, G0=G0, G1=G1,
#'   mu.method = "H1", cov.method = "H0")
#'   
#' #Mean(Z) under H1 with power based on covariance matrix under H1
#' wlr.power.maxcombo(DCO = c(24, 36),  
#'   alpha=c(0.01, 0.04)/2, 
#'   r = 1, n = 500, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, 
#'   f.ws = list(IA1 = list(lr), FA=list(fh01)), 
#'   Lambda=Lambda, G0=G0, G1=G1,
#'   mu.method = "H1", cov.method = "H1")
#'   
#' #Mean(Z) under H1 with power based on covariance matrix under H1 in Local Alternative (simplified)
#' wlr.power.maxcombo(DCO = c(24, 36),  
#'   alpha=c(0.01, 0.04)/2, 
#'   r = 1, n = 500, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, 
#'   f.ws = list(IA1 = list(lr), FA=list(fh01)), 
#'   Lambda=Lambda, G0=G0, G1=G1,
#'   mu.method = "H1", cov.method = "H1.LA")   
#'
#' #max-combo(logrank, FH11) at FA   
#' #Mean(Z) under H1 with power based on covariance matrix under H1 in Local Alternative (simplified)
#' wlr.power.maxcombo(DCO = c(24, 36),  
#'   alpha=c(0.01, 0.04)/2, 
#'   r = 1, n = 500, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, 
#'   f.ws = list(IA1 = list(lr), FA=list(lr, fh11)), 
#'   Lambda=Lambda, G0=G0, G1=G1,
#'   mu.method = "H1", cov.method = "H1.LA")   
#'   
#' #max-combo(logrank, FH11) at FA   
#' #Mean(Z) under H1 with power based on covariance matrix under H1
#' wlr.power.maxcombo(DCO = c(24, 36),  
#'   alpha=c(0.01, 0.04)/2, 
#'   r = 1, n = 500, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, 
#'   f.ws = list(IA1 = list(lr), FA=list(lr, fh11)), 
#'   Lambda=Lambda, G0=G0, G1=G1,
#'   mu.method = "H1", cov.method = "H1")   
#'         
#' @export
#' 
#' 
wlr.power.maxcombo = function(n = 600, r = 1, DCO = c(24, 36), 
            alpha=c(0.02, 0.03)/2,   
            h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
            h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
            f.ws = list(IA1=list(function(s){1}), FA=list(function(s){1}, function(s){s*(1-s)})), cuts=NULL,
            Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
            G0 = function(t){0}, G1 = function(t){0}, mu.method = "Schoenfeld", cov.method = "H0"){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  
  r1 = r / (r + 1); r0 = 1 - r1 

  #K analyses
  K=length(f.ws)

  #Number of test components in each analysis
  J = lengths(f.ws)
  
  ###################################################
  #(1)Calculate the target events according to DCO
  ###################################################
  targetEvents0 = targetEvents1 = rep(NA, K)
  for (i in 1:K){
    e.tmp = fe(n=n,DCO=DCO[i],r=r,h0=h0,S0=S0,h1=h1,S1=S1,Lambda=Lambda,G0=G0,G1=G1)
    targetEvents0[i] = e.tmp$e0
    targetEvents1[i] = e.tmp$e1
  }
    
  #Mean(Z)
  mu = matrix(NA, nrow=K, ncol=max(J))
  for(i in 1:K){
    for(j in 1:J[[i]]){
      mu[i, j] = wlr.mu(n=n, DCO = DCO[i], r=r, h0 = h0, S0=S0, h1 = h1, S1=S1,
                        rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, 
                        f.ws = f.ws[[i]][[j]], Lambda = Lambda, 
                        G0 = G0, G1 = G1, mu.method = mu.method)
    }
  }
  
  ###################################################
  #(2)Covariance matrix of all Z's: J[1]+J[2]+...+J[K] dimensions
  ###################################################
  
  #Omega0: Covariance matrix under H0
  #Omega1: Covariance matrix under H1
  Omega0 = Omega1 = matrix(NA, nrow=sum(J), ncol=sum(J))
  #Calculate the correlation between Zij and Z_i'j'
  for (i in 1:K){
    for (j in 1:J[[i]]){
      row = as.numeric(i>=2)*sum(J[1:(i-1)])+j #row location of the cov matrix
      Omega0[row, row] = 1
      Omega1[row, row] = wlr.Zcov(DCO = c(DCO[i], DCO[i]), r=r, 
                                  h0 = h0, S0=S0, h1 = h1, S1=S1,
                                  cuts = cuts, rho=NULL, gamma=NULL, 
                                  tau=NULL, s.tau=NULL,
                                  f.ws=c(f.ws[[i]][[j]],f.ws[[i]][[j]]), 
                                  Lambda = Lambda, G0 = G0, G1 = G1, 
                                  Hypo = cov.method)
      for (ip in i:K){
        for (jp in 1:J[ip]){
          #incremental location pointer for column compared to row of the cov matrix
          incr = (as.numeric(ip>=2)*sum(J[1:(ip-1)])+jp)-(as.numeric(i>=2)*sum(J[1:(i-1)])+j)
          col = row + incr
          #incr controls the computation only limited to upper right corner
          if(incr > 0){
            Omega0[row, col] = wlr.Zcov(DCO = c(DCO[i], DCO[ip]), r=r, 
                                        h0 = h0, S0=S0, h1 = h1, S1=S1,
                                        cuts = cuts, rho=NULL, gamma=NULL, 
                                        tau=NULL, s.tau=NULL,
                                        f.ws=c(f.ws[[i]][[j]],f.ws[[ip]][[jp]]), 
                                        Lambda = Lambda, G0 = G0, G1 = G1, 
                                        Hypo = "H0")
            Omega0[col, row] = Omega0[row, col]
            Omega1[row, col] = wlr.Zcov(DCO = c(DCO[i], DCO[ip]), r=r, 
                                        h0 = h0, S0=S0, h1 = h1, S1=S1,
                                        cuts = cuts, rho=NULL, gamma=NULL, 
                                        tau=NULL, s.tau=NULL,
                                        f.ws=c(f.ws[[i]][[j]],f.ws[[ip]][[jp]]), 
                                        Lambda = Lambda, G0 = G0, G1 = G1, 
                                        Hypo = cov.method)
            Omega1[col, row] = Omega1[row, col]            
          }
        }
      }
    }
  }
  
  ###################################################
  #(3)Rejection boundary: recursively solve for 
  #   the rejection boundary for each analysis
  ###################################################
  
  b = rep(NA, K)
  
  #First Analysis
  #maxcombo has at least 2 components
  if (J[[1]] >= 2){
      f.b = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
            mean=0, corr = Omega0[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }
      b[1] = uniroot(f=f.b, interval=c(1, 20), tol = 1e-8)$root
  } else{  b[1] = qnorm(1-alpha[1])  }
  
  #Recursively solve other boundaries from 2nd analysis to Kth analysis
  if(K > 1){
    for(i in 2:K){
        f.b = function(x){
          LL1 = rep(-Inf, sum(J[1:(i-1)]))
          UU1 = rep(b[1], J[1])
          if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b[j], J[j]))}}

          idx1 = 1:sum(J[1:(i-1)])
          if (length(UU1) == 1) {
            P1 = pnorm(UU1)
          } else {
            P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, 
                    mean=rep(0, length(idx1)),corr = Omega0[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
          }
          
          LL2 = rep(-Inf, sum(J[1:i]))
          UU2 = c(UU1, rep(x, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, 
                  mean=rep(0, length(idx)), corr = Omega0[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(P1 - P2 - alpha[i])
        }
        b[i] = uniroot(f=f.b, interval=c(1, 20), tol = 1e-8)$root
    }
  }
  
  ###################################################
  #(4)Power
  ###################################################

  power = rep(NA, K); #Power for each analysis individually

  for (i in 1:K){
    if(J[i] >= 2){
      mui = c(mu[i,]); mui = mui[!is.na(mui)]
      ix = (as.numeric(i>1)*sum(J[1:(i-1)])+1) : sum(J[1:i])
      power[i] = 1 - mvtnorm::pmvnorm(lower=rep(-Inf,J[i]),upper=rep(b[i], J[i]), 
                                      mean=mui, sigma = Omega1[ix, ix], abseps = 1e-8, maxpts=100000)[1]
    } else {power[i] = 1 - pnorm(b[i], mean=mu[i, 1])}  
  }
  
  #Overall power and incremental power
  overall.power = power[1]; incr.power = rep(0, K)
  
  if (K > 1){
      for (i in 2:K){
        LL1 = rep(-Inf, sum(J[1:(i-1)]))
        UU1 = rep(b[1], J[1])
        if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b[j], J[j]))}}
        mu.i1 = mu[1,] 
        if (i > 2) {for (j in 2:(i-1)){mu.i1 = c(mu.i1, mu[j,])}}
        mu.i1 = mu.i1[!is.na(mu.i1)] #mean vector up to (i-1)th analysis
        
        idx1 = 1:sum(J[1:(i-1)])
        if (length(LL1) == 1) {
          P1 = pnorm(UU1, mean=mu.i1)
        } else {
          P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, mean=mu.i1, 
                                sigma = Omega1[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
        }

        LL2 = rep(-Inf, sum(J[1:i]))
        UU2 = rep(b[1], J[1])
        for (j in 2:i){UU2 = c(UU2, rep(b[j], J[j]))}
        
        mu.i = mu[1,] 
        for (j in 2:i){mu.i = c(mu.i, mu[j,])}
        mu.i = mu.i[!is.na(mu.i)] #mean vector up to ith analysis
        
        idx = 1:sum(J[1:i]) #indices of the corr matrix
        P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, mean=mu.i, 
                              sigma = Omega1[idx, idx], abseps = 1e-8, maxpts=100000)[1]
        
        incr.power[i] = P1 - P2
        overall.power = overall.power + incr.power[i]
      }
  }
  ###################################################
  #(5)Misc
  ###################################################
  #Calculate the medians
  f.m0 = function(t){S0(t) - 0.5}
  f.m1 = function(t){S1(t) - 0.5}
  
  m0 = uniroot(f.m0, interval= c(1, 100), tol = 1e-8)$root
  m1 = uniroot(f.m1, interval= c(1, 100), tol = 1e-8)$root
  medians = c(m0, m1)
  targetEvents = targetEvents0 + targetEvents1
  
  #######Critical Values (Outside the scope of the manuscript)######
  CV.HR.H0 = CV.HR.H1 = rep(NULL, K)
  if (sum(J) == K){
    lr0 = 0
    for (i in 1:K){for (j in 1:J[i]){lr0 = f.ws[[i]][[j]](rnorm(1)) + lr0}}
    if (lr0 == K){
      #logrank test in all analyses
      for (i in 1:K){
        CV.HR.H0[i] = exp(-b[i]/sqrt(Lambda(DCO[i])*r0*r1*targetEvents[i]))
        r0.H1 = targetEvents0[i] / targetEvents[i]
        r1.H1 = 1 - r0.H1
        CV.HR.H1[i] = exp(-b[i]/sqrt(Lambda(DCO[i])*r0*r1*targetEvents[i]))
      }
    }
  }
  
  #Expected HR using un-weighted Cox regression
  Expected_HR = rep(NULL, K)
  for (i in 1:K){
    Expected_HR[i] = 0.001
  }
  
  o = list()
  o$Omega0 = Omega0
  o$Omega1 = Omega1
  o$mu = mu
  Analysis = 1:K
  o$design = data.frame(cbind(n, Analysis, DCO, targetEvents, power, incr.power, overall.power, b, medians, CV.HR.H0, CV.HR.H1, Expected_HR))
  return(o)
}

