#'  Asymptotic Covariance of Z1 and Z2 in Two Weighted Logrank Tests
#' 
#'  This function calculates the asymptotic covariance of two weighted logrank 
#'  tests Z1 and Z2, where Zj = nj**(-1/2)*Uj(t) / sqrt(Vj(t)/n).
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
#' wlr.Zcov(DCO = c(12, 24), r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = c(0,0), gamma = c(0,0), tau = c(NULL,NULL), s.tau = c(0,0), f.ws = list(f.ws1=NULL, f.ws2=NULL),
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0},
#' Hypo = "H0")
#' 
#' 
#' wlr.Zcov(DCO = c(24, 36), r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = c(0,0), gamma = c(0,1), tau = c(NULL,NULL), s.tau = c(0,0), f.ws = list(f.ws1=NULL, f.ws2=NULL),
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0},
#' Hypo = "H0")
#' 
#' wlr.Zcov(DCO = c(24, 36), r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = c(0,0), gamma = c(0,0), tau = c(NULL,NULL), s.tau = c(0,0), f.ws = list(f.ws1=NULL, f.ws2=NULL),
#' Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0},
#' Hypo = "H0")
#' 
#' wlr.Zcov(DCO = c(36, 36), r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=c(6), rho = c(0,0), gamma = c(1,1), tau = c(NULL,NULL), s.tau = c(0,0), f.ws = list(f.ws1=NULL, f.ws2=NULL),
#' Lambda = Lambda, G0 = function(t){0}, G1 = function(t){0},
#' Hypo = "H0")
#' 
#' eIA = fe(DCO = 24, r = 1, h0 = h0, S0 = S0, h1 = h1.D6, S1 = S1.D6, 
#' Lambda = Lambda, n = 600, G0 = G0, G1 = G1)$e
#' 
#' eFA = fe(DCO = 36, r = 1, h0 = h0, S0 = S0, h1 = h1.D6, S1 = S1.D6,
#' Lambda = Lambda, n = 600, G0 = G0, G1 = G1)$e
#' 
#' sigma2IA = wlr.sigma2(DCO = 24, r = 1, h0 = h0, S0 = S0, 
#' h1 = h1.D6, S1 = S1.D6, cuts=NULL, 
#' rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws=NULL,
#' Lambda = Lambda, G0 = G0, G1 = G1, Hypo="H0")
#' 
#' sigma2FA = wlr.sigma2(DCO = 36, r = 1, h0 = h0, S0 = S0, 
#' h1 = h1.D6, S1 = S1.D6, cuts=NULL, 
#' rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws=NULL,
#' Lambda = Lambda, G0 = G0, G1 = G1, Hypo="H0")
#' 
#' sigma2IA - eIA/(4*600)
#' sigma2FA - eFA/(4*600)
#' 
#' wlr.Zcov(DCO = c(24, 36), r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
#' cuts=NULL, rho = c(0,0), gamma = c(0,0), tau = c(NULL,NULL), s.tau = c(0,0), 
#' f.ws = list(f.ws1=NULL, f.ws2=NULL),
#' Lambda = Lambda, G0 = function(t){0}, G1 = function(t){0},
#' Hypo = "H0")
#' 
#' 
#' @export
#' 
#' 
wlr.Zcov = function(DCO = c(24, 32), r = 1, 
                    h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
                    h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
                    cuts=c(6), 
                    rho = c(0,0), gamma = c(0,0), tau = c(NULL,NULL), s.tau = c(0,0), 
                    f.ws = list(f.ws1=NULL, f.ws2=NULL),
                    Lambda = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                    G0 = function(t){0}, G1 = function(t){0},
                    Hypo = "H0"){
  
  #U1 and U2 covariance
  Ucov = wlr.Ucov(DCO = DCO, r = r, h0 = h0, S0 = S0, h1 = h1,S1 = S1, 
         cuts=cuts, rho = rho, gamma = gamma, tau = tau, s.tau = s.tau, f.ws = f.ws,
         Lambda = Lambda, G0 = G0, G1 = G1, Hypo = Hypo)$cov
  sigma1_2 = wlr.sigma2(DCO = DCO[1], r = r, h0 = h0, S0 = S0,h1 = h1, S1 = S1, 
                        cuts=cuts, rho = rho[1], gamma = gamma[1], tau = tau[1], s.tau = s.tau[1], f.ws=f.ws[[1]],
                        Lambda = Lambda, G0 = G0, G1 = G1, Hypo = Hypo)
  sigma2_2 = wlr.sigma2(DCO = DCO[2], r = r, h0 = h0, S0 = S0,h1 = h1, S1 = S1, 
                        cuts=cuts, rho = rho[2], gamma = gamma[2], tau = tau[2], s.tau = s.tau[2], f.ws=f.ws[[2]],
                        Lambda = Lambda, G0 = G0, G1 = G1, Hypo = Hypo)

  return(Ucov/sqrt(sigma1_2*sigma2_2))
}


