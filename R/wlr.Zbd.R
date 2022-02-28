#'  Rejection Boundary (z) for Weighted Logrank Test in Group Sequential Design
#' 
#'  Consider the general testing framework: At analysis i, the test is max(WLR_1, ..., WLR_ni).
#'  This function calculates the rejection boundary in z value for each analysis 
#'  with the alpha allocation is given. The calculation is based on the asymptotic
#'  multivariate normal distribution. So it can be a challenge if the number of dimensions is large, 
#'  in which case one can use simulation approach with a large number of draws like 100 million 
#'  in order to achieve high precision.
#'  
#'  
#' @param DCO  Analysis time. 
#' @param alpha A vector of allocated alpha levels for interim and final analyses. 
#' sum(alpha) = 0.025 (1-sided).
#' @param r  Randomization ratio of experimental arm : control arm as r:1. When r = 1, it is equal allocation. Default r = 1.
#' @param h0 Hazard function of control arm. 
#' @param S0 Survival function of control arm.
#' @param h1 Hazard function of experimental arm. 
#' @param S1 Survival function of experimental arm. 
#' @param f.ws  Weight functions used for each analysis.
#'              For example (1), if there are 3 analyses planned including IA1, IA2, and FA.
#'              IA1 uses log-rank test, IA2 uses a maxcombo test of (logrank, FH01),
#'              and FA uses a maxcombo test of (logrank, FH01, FH11). Then specify as
#'              f.ws = list(IA1=list(lr), IA2=list(lr, fh01), FA=list(lr, fh01, fh11)),
#'              where define lr = function(s){1}; fh01=function(s){1-s}; fh11 = function(s){s*(1-s)};
#'              For example (2), if only fh01 is used for all three analyses IA1, IA2 and FA.
#'              f.ws = list(IA1=list(fh01), IA2=list(fh01), FA1=list(fh01)).
#'              For example (3), if logrank is used in a study design with only one analysis, then
#'              f.ws = list(FA=list(lr)).
#' @param Lambda Distribution function of enrollment. For example, 
#'        Lambda(t) = (t/A)^psi*I(0<=t<=A) + I(t>A) for enrollment in (0, A) with weight psi.
#' @param G0 Cumulative distribution function of drop-off for control arm, eg, G.ltfu=function(t){1-exp(-0.03/12*t)}
#'               is the distribution function for 3 percent drop-off in 12 months of followup.
#' @param G1 Cumulative distribution function of drop-off for experimental arm.
#'
#' @return 
#'  \itemize{
#'  \item  cov: Asymptotic covariance matrix, also correlation matrix under H0.
#'  \item  bd:  A vector of rejection boundaries for all analyses.
#'  }
#'  
#' @examples 
#' m0 = 12; #median RFS for control arm
#' lam0 = log(2) / m0
#' h0 = function(t){lam0}; 
#' S0 = function(t){exp(-lam0 * t)}
#' HRd = 0.60 #hazard ratio after delay
#' 
#' h.D6=function(t){lam0*as.numeric(t<6)+HRd*lam0*as.numeric(t>=6)}
#' c6 = exp(-6*lam0*(1-HRd)); 
#' S.D6 = function(t){exp(-lam0*t)*as.numeric(t<6)+c6*exp(-HRd*lam0*t)*as.numeric(t>=6)}
#' f.logHR.D6 = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*HRd)}
#' 
#' #Define weight functions for weighted log-rank tests
#' lr = function(s){1}
#' fh01 = function(s){(1-s)}
#' fh11 = function(s){s*(1-s)}
#' #modestly log-rank
#' mfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1/s1)}
#' 
#' #Define the testing strategy
#' test.IA1 = list(lr)
#' test.IA2 = list(lr, fh11)
#' test.IA3 = list(lr, fh01, fh11)
#' 
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' G0 = G1 = function(t){0} 
#'   
#' wlr.Zbd(DCO = c(24, 36, 48), r = 1, alpha=c(0.01, 0.02, 0.02)/2,
#' h0 = h0, S0 = S0, h1 = h.D6, S1 = S.D6, 
#' f.ws = list(IA1=test.IA1, IA2=test.IA2, FA=test.IA3),
#' Lambda = Lambda, G0 = G0, G1 = G1)
#'   
#' @export
#' 
wlr.Zbd = function(DCO = c(24, 36, 48), r = 1, alpha=c(0.01, 0.02, 0.02)/2,
                   h0 = function(t){log(2)/12}, 
                   S0= function(t){exp(-log(2)/12 * t)},
                   h1 = function(t){log(2)/12*0.70}, 
                   S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
                   f.ws = list(IA1=list(lr), IA2=list(lr, fh11), FA=list(lr, fh01, fh11)),
                   Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                   G0 = function(t){0}, G1 = function(t){0}){
  
  r1 = r / (r + 1)  #Proportion of subjects in experimental arm
  r0 = 1 - r1       #Proportion of subjects in control arm
  K=length(f.ws)    #K analyses
  J = lengths(f.ws) #Number of test components in each analysis

  ###############################
  #Correlation matrix under H0
  ###############################
  
  #Complete correlation matrix structure: J[1]+J[2]+...+J[K] dimensions
  Omega = matrix(1, nrow=sum(J), ncol=sum(J))
  #calculate the correlation between Zij and Z_i'j'
  for (i in 1:K){
    for (j in 1:J[[i]]){
      for (ip in i:K){
        for (jp in 1:J[ip]){
          row = as.numeric(i>=2)*sum(J[1:(i-1)])+j #row location of the corr matrix
          
          #incremental location pointer for column compared to row of the corr matrix
          incr = (as.numeric(ip>=2)*sum(J[1:(ip-1)])+jp)-(as.numeric(i>=2)*sum(J[1:(i-1)])+j)
          col = row + incr
          #incr controls the computation only limited to upper right corner
          if(incr > 0){
            #information matrix for Zij and Zi'j'
            zcov = wlr.Zcov(DCO = c(DCO[i], DCO[ip]), r = r, 
                            h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
                            rho=NULL, gamma=NULL, tau=NULL, s.tau=NULL, 
                            f.ws = c(f.ws[[i]][[j]],f.ws[[ip]][[jp]]),
                            Lambda = Lambda, G0 = G0, G1 = G1, Hypo = "H0")
            
            Omega[col, row] = Omega[row, col] = zcov
          }
        }
      }
    }
  }
  
  ######################################
  #Rejection boundary:
  #Recursively solve for each analysis
  ######################################
  bd = rep(NA, K) #Rejection boundary in z value
  
  #-----------------
  #First Analysis
  #-----------------
  
  #If at least 2 components in 1st analysis
  if (J[[1]] >= 2){
    f.bd1 = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]), upper=rep(x, J[1]), 
                           mean=0, corr = Omega[1:J[1],1:J[1]], 
                           abseps = 1e-8, maxpts=100000)[1] - alpha[1]
    }
    bd[1] = uniroot(f=f.bd1, interval=c(1, 20), tol = 1e-8)$root
  } else { bd[1] = qnorm(1-alpha[1])}
  
  #-----------------
  #Further Analysis
  #-----------------
  
  #Recursively solve other boundaries from 2nd analysis to Kth analysis
  if(K > 1){
    for(i in 2:K){
      f.bd = function(x){
        LL1 = rep(-Inf, sum(J[1:(i-1)]))
        UU1 = rep(bd[1], J[1])
        if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(bd[j], J[j]))}}
        
        idx1 = 1:sum(J[1:(i-1)])
        if (length(UU1) == 1) {
          P1 = pnorm(UU1)
        } else {
          P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, 
                                mean=rep(0, length(idx1)),
                                corr = Omega[idx1, idx1], 
                                abseps = 1e-8, maxpts=100000)[1]
        }
          
        LL2 = rep(-Inf, sum(J[1:i]))
        UU2 = c(UU1, rep(x, J[i]))
          
        idx = 1:sum(J[1:i]) #indices of the corr matrix
        P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, 
                              mean=rep(0, length(idx)), 
                              corr = Omega[idx, idx], 
                              abseps = 1e-8, maxpts=100000)[1]
        return(P1 - P2 - alpha[i])
      }
      bd[i] = uniroot(f=f.bd, interval=c(1, 20), tol = 1e-8)$root
    }
  }
  
  o = list()
  o$cov = Omega; o$bd = bd
  
  return(o)
}

