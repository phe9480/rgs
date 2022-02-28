#'  Asymptotic Covariance of n1^(-1/2)U1(t) and n2^(-1/2)U2(t) in Two Weighted Logrank Tests
#' 
#'  This function calculates the asymptotic covariance of n1**(-1/2)*U1(t) and 
#'  n2**(-1/2)*U2(t) in two weighted logrank tests Z1 and Z2,
#'  where Zj = nj**(-1/2)*Uj(t) / sqrt(Vj(t)/n).
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
#' @param G1 Cumulative distribution function of drop-off for experimental arm.
#' @param Hypo Hypothesis: "H0", "H1", "H1.LA". Under H0, the pooled hazard h.bar(t) = f.tilda(t)/S.tilda(t).
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
#' wlr.Ucov(DCO = c(12, 24), r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = c(0,0), gamma = c(0,0), tau = c(NULL,NULL), s.tau = c(0,0), f.ws = list(f.ws1=NULL, f.ws2=NULL),
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0},
#' Hypo = "H1")
#' 
#' wlr.Ucov(DCO = c(24, 24), r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = c(0,0), gamma = c(0,1), tau = c(NULL,NULL), s.tau = c(0,0), f.ws = list(f.ws1=NULL, f.ws2=NULL),
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0},
#' Hypo = "H0")
#' 
#' @export
#' 
#' 
wlr.Ucov = function(DCO = c(24, 32), r = 1, 
                    h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
                    h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
                    cuts=c(6), 
                    rho = c(0,0), gamma = c(0,0), tau = c(NULL,NULL), s.tau = c(0,0), 
                    f.ws = list(f.ws1=NULL, f.ws2=NULL),
                    Lambda = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                    G0 = function(t){0}, G1 = function(t){0},
                    Hypo = "H0"){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 
  
  #Density functions
  f0 = function(t) {return(h0(t) * S0(t))}
  f1 = function(t) {return(h1(t) * S1(t))}  
  #Pooled survival function for weight function
  f.bar = function(t){r0 * f0(t) + r1 * f1(t)} 
  S.bar = function(t){r0 * S0(t) + r1 * S1(t)} 
  
  #At-risk function and density function adjusted by censoring mechanism
  S0.tilda = function(t){S0(t)*(1-G0(t))}
  S1.tilda = function(t){S1(t)*(1-G1(t))}
  f0.tilda = function(t){f0(t)*(1-G0(t))}
  f1.tilda = function(t){f1(t)*(1-G1(t))}
  
  S.tilda = function(t){r0 * S0.tilda(t) + r1 * S1.tilda(t)}
  f.tilda = function(t){r0 * f0.tilda(t) + r1 * f1.tilda(t)}
  
  #Proportion of at-risk at each arm
  pi1 = function(t){r1 * S1.tilda(t) / S.tilda(t)}
  pi0 = function(t){1 - pi1(t)}
  #pi0 and pi1 under H0 (ie f0 = f1 = f.bar)
  pi1.H0 = function(t){r1 * (1-G1(t)) / (r0 * (1-G0(t)) + r1 * (1-G1(t)))}
  pi0.H0 = function(t){1 - pi1.H0(t)}
  f.tilda.H0 = function(t){(r0 * (1-G0(t)) + r1 * (1-G1(t)))*f.bar(t)}
  
  #f.delta function
  f.delta = function(t){
    pi0(t)*pi1(t)*S.tilda(t)*(h1(t)-h0(t))
  }
  
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
  at = unique(c(0, cuts[cuts <= DCO[1]], DCO[1]))
  
  if (Hypo =="H0"){
    #Cv1(t1, t2; w1, w2)
    f.Cv1a = function(t){
      w1 = f.w(t, f.S = S.bar, f.ws=f.ws[[1]], tau=tau[1], s.tau=s.tau[1], rho=rho[1], gamma=gamma[1])
      w2 = f.w(t, f.S = S.bar, f.ws=f.ws[[2]], tau=tau[2], s.tau=s.tau[2], rho=rho[2], gamma=gamma[2])
      tmp = pi0.H0(t)*pi1.H0(t)*f.tilda.H0(t)
      return(w1*w2*Lambda(DCO[1]-t)*tmp)
    }
    
    Cv1a = 0
    for (j in 1:(length(at)-1)){
      Cv1aj = integrate(f.Cv1a, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
      Cv1a = Cv1a + Cv1aj
    }
    Cv1 = Cv1a
    Delta1 = Delta2 = 0
  } else if (Hypo =="H1"){
    #Cv1(t1, t2; w1, w2)
    f.Cv1a = function(t){
      w1 = f.w(t, f.S = S.bar, f.ws=f.ws[[1]], tau=tau[1], s.tau=s.tau[1], rho=rho[1], gamma=gamma[1])
      w2 = f.w(t, f.S = S.bar, f.ws=f.ws[[2]], tau=tau[2], s.tau=s.tau[2], rho=rho[2], gamma=gamma[2])
      tmp = pi0(t)*pi1(t)*S.tilda(t)*(pi0(t)*h1(t)+pi1(t)*h0(t))
        
      return(w1*w2*Lambda(DCO[1]-t)*tmp)
    }
    
    Cv1a = 0
    for (j in 1:(length(at)-1)){
      Cv1aj = integrate(f.Cv1a, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
      Cv1a = Cv1a + Cv1aj
    }
    
    #Delta1 and Delta2
    Delta1 = wlr.Delta(DCO = DCO[1], r = r, h0 = h0, S0 = S0,h1 = h1, S1 = S1, 
                       cuts=cuts, rho = rho[1], gamma = gamma[1], tau = tau[1], s.tau = s.tau[1], f.ws=f.ws[[1]],
                       Lambda = Lambda, G0 = G0, G1 = G1)
    
    Delta2 = wlr.Delta(DCO = DCO[2], r = r, h0 = h0, S0 = S0,h1 = h1, S1 = S1, 
                       cuts=cuts, rho = rho[2], gamma = gamma[2], tau = tau[2], s.tau = s.tau[2], f.ws=f.ws[[2]],
                       Lambda = Lambda, G0 = G0, G1 = G1)
    Cv1 = Cv1a - Delta1*Delta2
  } else if (Hypo =="H1.LA"){ 
    #simplified covariance under Local alternative
    #Cv1(t1, t2; w1, w2)
    f.Cv1a = function(t){
      w1 = f.w(t, f.S = S.bar, f.ws=f.ws[[1]], tau=tau[1], s.tau=s.tau[1], rho=rho[1], gamma=gamma[1])
      w2 = f.w(t, f.S = S.bar, f.ws=f.ws[[2]], tau=tau[2], s.tau=s.tau[2], rho=rho[2], gamma=gamma[2])
      tmp = pi0(t)*pi1(t)*f.tilda(t)
      
      return(w1*w2*Lambda(DCO[1]-t)*tmp)
    }
    
    Cv1a = 0
    for (j in 1:(length(at)-1)){
      Cv1aj = integrate(f.Cv1a, lower=at[j], upper=at[j+1], abs.tol=1e-8)$value
      Cv1a = Cv1a + Cv1aj
    }
    Cv1 = Cv1a
    Delta1 = Delta2 = 0
  }
  
  ####Cv2(t1, t2; w1, w2)
  #s1 = s' in the range of (0, DCO1) and s2 in (0, DCO2)
  #s[1] = s' and s[2] = s for the notations in the manuscript
  f.Cv2 = function(s){
    I = as.numeric(s[1] >= s[2])
    #I = ifelse(s[1] >= s[2], 1, 0)
    w1 = f.w(s[1], f.S = S.bar, f.ws=f.ws[[1]], tau=tau[1], s.tau=s.tau[1], rho=rho[1], gamma=gamma[1])
    w2 = f.w(s[2], f.S = S.bar, f.ws=f.ws[[2]], tau=tau[2], s.tau=s.tau[2], rho=rho[2], gamma=gamma[2])
    Lam = Lambda(min(DCO[1]-s[1], DCO[2]-s[2]))
    fa = f.tilda(s[2])/S.tilda(s[2])*(pi1(s[2])-pi0(s[1]))*f.delta(s[1])
    fb = f.tilda(s[1])/S.tilda(s[1])*(pi1(s[1])-pi0(s[2]))*f.delta(s[2])    

    return(w1*w2*Lam*(I*fa+(1-I)*fb)*100)
  }
  
  #Vectorize f.Cv2 for double integration
  f.Cv2.v <- function(s) {
    matrix(apply(s, 2, function(z) f.Cv2(z)), ncol = ncol(s))
  }

  Cv2 = 0 #For both H0 and H1.LA, Cv2 = 0
  if(Hypo == "H1"){
    Cv2 = cubature::adaptIntegrate(f.Cv2.v, rep(0,2), DCO, vectorInterface = TRUE)$integral/100
  }
  cov = sqrt(Lambda(DCO[1])/Lambda(DCO[2]))*(Cv1 + Cv2)
  
  o = list()
  o$cov=cov; o$Delta = c(Delta1, Delta2); 
  o$Cv1a = Cv1a; o$Cv1 = Cv1; o$Cv2 = Cv2
  return(o)
}


