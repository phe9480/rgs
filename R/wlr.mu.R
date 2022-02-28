#'  The Expectation of Weighted Logrank Test Z = U(t) / sqrt(V(t)) at analysis time t
#' 
#'  Delta = E(U(t)/n); V(t)/n --> sigma2 
#'  Z = n^(-1/2)U(t)/sqrt(V(t)/n), so E(Z) = mu = n^(1/2)Delta/sigma
#'  
#' @param DCO  Analysis time, calculated from first subject in.
#' @param r  Randomization ratio of experimental arm : control arm as r:1. 
#'           When r = 1, it is equal allocation. Default r = 1.
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
#' @param mu.method Method of mean(Z): "Schoenfeld" or "H1".
#'  
#' 
#' @examples 
#' #1:1 randomization, control~exp(median=12);
#' #experiment~delayed effect 6 mo and HR 0.6 after delay
#' #Fleming-Harrington (1,1) test; Enrollment: 18 mo with weight 1.5
#' #drop-off rate: 3% in 1-year of follow-up in control arm; no drop-off in exp arm
#' HR = 0.6; delay = 6; lam0 = log(2) / 12; 
#' h0 = function(t){lam0}; S0 = function(t){exp(-lam0 * t)}
#' h1.D6 = function(t){lam0*as.numeric(t < delay)+HR*lam0*as.numeric(t >= delay)}
#' c = exp(-delay*lam0*(1-HR)); 
#' S1.D6 = function(t){exp(-lam0*t)*as.numeric(t<delay) + c*exp(-HR*lam0*t)*as.numeric(t>=delay)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' drop0 = 0.03/12; drop1 = 0
#' 
#' Delta = wlr.Delta(DCO = 24, r = 1, h0 = h0, S0 = S0, h1 = h1.D6, S1 = S1.D6, 
#' cuts=c(6), rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0})
#' #-0.02315818
#' 
#' sig2=wlr.sigma2(DCO = 24, r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0})
#' 
#' #sqrt(600)*Delta/sqrt(sig2)
#' 
#' wlr.mu(DCO = 24, r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0},
#' mu.method="H1")
#' 
#' wlr.mu(DCO = 24, r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0},
#' mu.method="Schoenfeld")
#' 
#' 
#' @export
 
wlr.mu = function(DCO = 24, n=600, r = 1, 
                  h0 = function(t){log(2)/12}, 
                  S0=function(t){exp(-log(2)/12 * t)},
                  h1 = function(t){log(2)/12*0.70}, 
                  S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
                  cuts=c(6), 
                  rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
                  Lambda = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                  G0 = function(t){0}, G1 = function(t){0}, mu.method="H1"){

  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 
  
  #mu = sqrt(n) * Delta / sqrt(sigma2)
  
  Delta = wlr.Delta(DCO = DCO, r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
                    cuts=cuts, rho = rho, gamma = gamma, tau = tau, 
                    s.tau = s.tau, f.ws = f.ws, 
                    Lambda = Lambda, G0 = G0, G1 = G1)
  sigma2 = wlr.sigma2(DCO = DCO, r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
                      cuts=cuts, rho = rho, gamma = gamma, tau = tau, 
                      s.tau = s.tau, f.ws = f.ws, 
                      Lambda = Lambda, G0 = G0, G1 = G1, Hypo="H1")
  
  if(mu.method == "H1") {
    mu = -sqrt(n*Lambda(DCO))*Delta / sqrt(sigma2)
  } else if (mu.method == "Schoenfeld") {
    #Density functions
    f0 = function(t) {return(h0(t) * S0(t))}
    f1 = function(t) {return(h1(t) * S1(t))}  
    f0.tilda = function(t){f0(t)*(1-G0(t))}
    f1.tilda = function(t){f1(t)*(1-G1(t))}
    
    f.tilda = function(t){r0 * f0.tilda(t) + r1 * f1.tilda(t)}
    S0.tilda = function(t){S0(t)*(1-G0(t))}
    S1.tilda = function(t){S1(t)*(1-G1(t))}
    S.tilda = function(t){r0 * S0.tilda(t) + r1 * S1.tilda(t)}
    S.bar = function(t){r0 * S0(t) + r1 * S1(t)} 
    
    #Proportion of at-risk at each arm
    pi1 = function(t){r1 * S1.tilda(t) / S.tilda(t)}
    pi0 = function(t){1 - pi1(t)}
    
    #Weight function
    f.w = function(t, f.S, f.ws, tau, s.tau, rho, gamma){
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
    
    ##Intervals for piecewise integration to avoid occasional numerical accuracy issues
    at = unique(c(0, cuts[cuts <= DCO], DCO))
    
    ## E[Z]: Schoenfeld method ##eta == "Delta" in H1 scenario
    I.eta = function(t){
      w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
      return(w * (log(h0(t)/h1(t)))*pi0(t)*pi1(t)*Lambda(DCO-t)*f.tilda(t))
    }
    eta = 0
    for (j in 1:(length(at)-1)){
      etaj = integrate(I.eta, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
      eta = eta + etaj
    }
    mu = sqrt(n * Lambda(DCO)) * eta / sqrt(sigma2)
  }
  
  return(mu)  
}

