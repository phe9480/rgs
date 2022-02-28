#' Asymptotic Distribution of Weighted Log-rank Statistic
#' 
#'  This function calculates the asymptotic distribution (mean and variance) at a DCO.
#'  
#' @param T  Analysis time, calculated from first subject in.
#' @param r  Proportion of subjects in experimental arm. For 1:1 randomization, r = 1/2.
#' @param n  Total sample size for two arms. 
#' @param h0 Hazard function of control arm, eg, h0(t) = log(2)/m0 for the exponential distribution with median m0.
#' @param S0 Survival function of control arm. 
#' @param h1 Hazard function of experimental arm. 
#' @param S1 Survival function of experimental arm.
#' @param cuts Piece wise interval cut points of the curresponding survival functions. 
#' For delayed effect model with 6 months of delay, cuts = 6. Incorrect 
#' specification of cuts may lead to numerical integration issues.
#' @param f.ws  Self-defined weight function of survival rate, eg, f.ws = function(s){1/max(s, 0.25)}
#'         When f.ws is specified, sFH parameter will be ignored.
#' @param Lambda Cumulative distribution function of enrollment. 
#' @param G0 Cumulative distribution function of drop-off for control arm, eg, G.ltfu=function(t){1-exp(-0.03/12*t)}
#'               is the distribution function for 3 percent drop-off in 12 months of followup.
#' @param G1 Cumulative distribution function of drop-off for experimental arm, eg, G.ltfu=function(t){1-exp(-0.03/12*t)}
#'               is the distribution function for 3 percent drop-off in 12 months of followup.
#' @param method Non-centrality parameter options: "SIGWA", "Schoenfeld", "YL 2020", "LM 2019". Default "SIGWA". 
#' "YL 2020" is based on Yung and Liu (2020) and "LM 2019" is based on Luo et al (2019).
#'
#' @return A dataframe with variables including
#'  \itemize{
#'  \item  mu:     Mean of the weighted logrank test statistic Z = U/sqrt(V). E(Z) = mu = sqrt(n)*delta/sigma
#'  \item  mu_SF:  Mean of the weighted logrank test statistic Z using Schoenfeld method. E(Z) = mu_SF = sqrt(n)*eta/sigma
#'  \item  delta:  Mean of U/n. E(U/n)=delta = r0 * delta0 + r1 * delta1
#'  \item  delta0: Mean of U0/n0 for control arm. E(U0/n0) = delta0
#'  \item  delta1: Mean of U1/n1 for experimental arm E(U1/n1) = delta1
#'  \item  eta:    Mean of U/n in weighted logrank test in Schoenfeld method. E(U/n)=eta
#'  \item  sigma2: Estimate of V/n --> sigma2 for n large. Also, under H0, var(n^(-1/2)U) = sigma2 
#'  Under H1, var(n^(-1/2)U) = sigma2 + O(n^(-1/2)). So sigma2 is also an estimate of var(n^(-1/2)U)
#'  \item  sigma2_b: Var(n^(-1/2)U) using Yung and Liu (2020)'s method Z ~ N(mu, sigma2_b/sigma2)
#'  \item  sigma2_s: Var(n^(-1/2)U) using Luo et al (2019)'s method Z ~ N(mu, sigma2_s/sigma2)
#'  \item  sigma2_h: Var(n^(-1/2)U) using He et al (2022)'s method Z ~ N(mu, sigma2_h/sigma2)
#'  \item  R:  Location Factor. Under H1, mu = sqrt(n)*R. For a fixed sample size n, R represents the magnitude of followup.
#'  \item  var_YL: var_YL = sigma2_b/sigma2
#'  \item  var_LM: var_LM = sigma2_s/sigma2
#'  \item  var_HK: var_HK = sigma2_h/sigma2
#'  \item  pow_YL: Power for a test with rejection bound 1.96: P(Z > 1.96)=PHI((mu-1.96)*sigma/sigma_x). 
#'  2.5% alpha with only 1 analysis. 
#'  \item  pow_LM: See pow_YL 
#'  \item  pow_HK: See pow_YL
#'  \item  pow_SF: Power of Schoenfeld method: PHI(mu_SF-1.96)
#'  \item  pow_mu: Power of mu method with simpler approximation: PHI(mu-1.96)
#'  }
#'  Another dataframe including variables: method, mean, var, pow
#' @examples 
#' #Example: Assume HR 0.65 after 6 months of delayed effect. Control arm follows
#' exponential distribution with median 12 months.
#' 
#' HR = 0.4; delay = 6; lam0 = log(2) / 12; 
#' h0 = function(t){lam0}; S0 = function(t){exp(-lam0 * t)}
#' h1.D6 = function(t){lam0*as.numeric(t < delay)+HR*lam0*as.numeric(t >= delay)}
#' c = exp(-delay*lam0*(1-HR)); 
#' S1.D6 = function(t){exp(-lam0*t)*as.numeric(t<delay) + c*exp(-HR*lam0*t)*as.numeric(t>=delay)}
#' f.logHR.D6 = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*HR)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' G0.ltfu = G1.ltfu = function(t){0}
#' 
#' wlr.mu(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = Lambda)
#' wlr.asym(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, cuts=c(6))
#' wlr.asym(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, cuts=c(6))
#' wlr.asym(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, cuts=c(6))
#' wlr.mu.schoenfeld(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = Lambda)
#'      
#' wlr.asym(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, cuts=c(6))
#' wlr.asym(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, cuts=c(6))
#' wlr.asym(DCO = 30, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, cuts=NULL)
#'      
#' 
#' @export
#' 
wlr.asym = function(DCO = 24, r = 1, n = 450, 
       h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
       h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
       cuts=c(6), 
       rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
       Lambda = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
       G0 = function(t){0}, G1 = function(t){0}){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 
  
  #Density functions
  f0 = function(t) {return(h0(t) * S0(t))}
  f1 = function(t) {return(h1(t) * S1(t))}  
  
  #Censoring mechanism
  H0 = function(t) {Lambda(DCO-t) * (1-G0(t))}
  H1 = function(t) {Lambda(DCO-t) * (1-G1(t))}

  #Pooled survival function for weight function
  S.bar = function(t){r0 * S0(t) + r1 * S1(t)} 
  
  #At-risk function and density function adjusted by censoring mechanism
  S.tilda = function(t){r0 * S0(t)*H0(t) + r1 * S1(t)*H1(t)}
  f.tilda = function(t){r0 * f0(t)*H0(t) + r1 * f1(t)*H1(t)}
  
  #Proportion of at-risk at each arm
  pi1 = function(t){r1 * S1(t) * H1(t) / S.tilda(t)}
  pi0 = function(t){1 - pi1(t)}
  
  #Weight function
  f.w = function(t, f.S = S0, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma){
    s = f.S(t)
    #First priority: f.ws
    if(!is.null(f.ws)){
      w = f.ws(s)
    }else {
      #Second priority: s.tau
      if (!is.null(s.tau)){
        s.til = apply(cbind(s, s.tau), MARGIN=1,FUN=max);
      } else {
        s.til = apply(cbind(s, f.S(tau)), MARGIN=1,FUN=max);        
      }
      w = s.til^rho*(1-s.til)^gamma
    }
    return(w)
  }
  
  ##Intervals for piecewise integration
  at = unique(c(0, cuts[cuts <= DCO], DCO))
  
  ## delta ##
  I.delta = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w *pi0(t)* pi1(t)*S.tilda(t)*(h0(t)-h1(t)))
  }
  delta = 0
  for (j in 1:(length(at)-1)){
    deltaj = integrate(I.delta, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
    delta = delta + deltaj
  }
  
  ## sigma2 ##
  I.sig2 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w^2*pi0(t)*pi1(t)*f.tilda(t))
  }
  sigma2 = 0
  for (j in 1:(length(at)-1)){
    sigma2j = integrate(I.sig2, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
    sigma2 = sigma2 + sigma2j
  }
  
  ## E[Z]: mu for all methods except for Schoenfeld ##
  R = delta / sqrt(sigma2)
  mu = sqrt(n * Lambda(DCO)) * R
  
  ## E[Z]: Schoenfeld method ##
  I.eta = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * (log(h0(t)/h1(t)))*pi0(t)* pi1(t)* f.tilda(t))
  }
  eta = 0
  for (j in 1:(length(at)-1)){
    etaj = integrate(I.eta, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
    eta = eta + etaj
  }
  mu_SF = sqrt(n * Lambda(DCO)) * eta / sqrt(sigma2)
  
  ## var(n^(-1/2)U): Y&L method ##
  I.delta0 = function(t){
      w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
      - w * pi1(t) * (f0(t)*H0(t) - 1/r0*pi0(t)*f.tilda(t))
  }
  delta0 = 0
  for (j in 1:(length(at)-1)){
    delta0j = integrate(I.delta0, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
    delta0 = delta0 + delta0j
  } 
  #Reverse the sign for consistency with the manuscript
  delta0 = - delta0
   
  I.delta1 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    w * pi0(t) * (f1(t)*H1(t) - 1/r1*pi1(t)*f.tilda(t))
  }
  delta1 = 0
  for (j in 1:(length(at)-1)){
    delta1j = integrate(I.delta1, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
    delta1 = delta1 + delta1j
  }
  delta1 = - delta1
    
  ## sigma2j: j = 0 or 1 ##
  f.sigma2j = function(j){
      
    #Inside integral
    inner <- function(x) {
      tmp = function(x1){
        I.inner = function(x2){
          w = f.w(x2, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
          w * (j - pi1(x2)) * f.tilda(x2) / S.tilda(x2)
        }
        integrate(I.inner, lower=0, upper=x1)$value
      }  
      sapply(x, tmp)
    }
      
    I.sigma2jA = function(t){
      w=f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma) 
      fj.tilda = f0(t)*H0(t)*as.numeric(j==0) + f1(t)*H1(t)*as.numeric(j==1)
      w^2*(j - pi1(t))^2 * fj.tilda 
    }
    sigma2jA = 0
    for (k in 1:(length(at)-1)) {
      sigma2jA = sigma2jA + integrate(I.sigma2jA,lower=at[k], upper=at[k+1])$value
    }
      
    sigma2jB <- -1 * (delta0*as.numeric(j==0) + delta1*as.numeric(j==1))^2
    
    I.sigma2jC = function(t){
      w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma) 
      Sj.tilda = S0(t) * H0(t) * as.numeric(j==0) + S1(t) * H1(t) * as.numeric(j==1)
      fj.tilda = f0(t) * H0(t) * as.numeric(j==0) + f1(t) * H1(t) * as.numeric(j==1)
      w*(j - pi1(t)) * (f.tilda(t) / S.tilda(t) * Sj.tilda - fj.tilda) * inner(t)
    }
      
    sigma2jC = 0
    for (k in 1:(length(at)-1)){
      sigma2jC = sigma2jC + 2 * integrate(I.sigma2jC, lower=at[k], upper=at[k+1])$value
    }
    sigma2jA + sigma2jB + sigma2jC
  }
  sigma2_0 = f.sigma2j(0)
  sigma2_1 = f.sigma2j(1)
  
  sigma2_b = r0 * sigma2_0 + r1 * sigma2_1
  var_YL = sigma2_b / sigma2
  
  ## var(n^(-1/2)U): L&M method ##
  sigma2_s = r0*r1*(delta0 - delta1)^2 + sigma2_b
  var_LM = sigma2_s / sigma2
  
  ## var(n^(-1/2)U): H&K method ##
  ## sigma2_h ##
  I.sig2h = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w^2*pi0(t)*pi1(t)*(r1*h0(t)*S1(t)*H1(t)+r0*h1(t)*S0(t)*H0(t)))
  }
  sigma2_h = 0
  for (j in 1:(length(at)-1)){
    sigma2hj = integrate(I.sig2h, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
    sigma2_h = sigma2_h + sigma2hj
  }
  var_HK = sigma2_h / sigma2
  
  ## POWER CALCULATION ASSUMING 2.5% ALPHA ##
  za = qnorm(1-0.025)
  pow_SF = pnorm(mu_SF - za)
  pow_YL = pnorm((mu-za)/sqrt(var_YL))
  pow_LM = pnorm((mu-za)/sqrt(var_LM))
  pow_HK = pnorm((mu-za)/sqrt(var_HK))
  pow_mu = pnorm(mu-za)

  o1 = data.frame(cbind(mu, mu_SF, delta, delta0, delta1, 
         eta, sigma2, sigma2_b, sigma2_s, sigma2_h, R, 
         var_YL, var_LM, var_HK, 
         pow_YL, pow_LM, pow_HK, pow_SF, pow_mu))
  o2 = NULL
  o2$method = c("Schoenfeld", "YL", "LM", "HK", "mu")
  o2$mean = c(mu_SF, mu, mu, mu, mu)
  o2$var  = c(1, var_YL, var_LM, var_HK, 1)
  o2$pow  = c(pow_SF, pow_YL, pow_LM, pow_HK, pow_mu)
  o = list(); o$par = o1; o$dist = data.frame(o2)
  
  return(o)
}
