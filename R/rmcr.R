#' Random Number Generator for Mixture Cure Rate Distribution
#' 
#' Consider the survival function of mixture cure rate (MCR) distribution: 
#' S(t) = p + (1-p)S0(t), where S0(t) is a survival distribution 
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
#' S(t) = p + (1-p)*S0(t)
#' S1(t) = S(t)I(t<tau) + S(tau)^(1-psi)*S(t)^psi
#' In brief, S0(t) is the proper GMW distribution (alpha, beta, gamma, lambda);
#' S(t) is MCR with additional cure rate parameter p;
#' In reference to S(t), S1(t) is a delayed effect distribution and proportional hazards after delay.
#' 
#' MCR(p, alpha, beta, gamma, lambda, tau, psi=1) reduces to S(t)
#' MCR(p, alpha, beta, gamma, lambda, tau=0, psi) reduces to S(t)^psi, i.e. proportional hazards.
#' MCR(p=0, alpha, beta, gamma, lambda, tau=0, psi=1) reduces to S0(t)
#' MCR(p=0, alpha, beta=1, gamma=1, lambda=0, tau=0, psi=1) reduces to exponential dist.
#' 
#' Draw u from U(0, 1); 
#'   if u >= 1 - S(tau)^(1-psi)*p^psi, let t = Infinity; 
#'   otherwise, let v = ((1-u1)^(1/psi) / (S(tau)^(1-psi)) - p) / (1-p),
#'   and t = inv.S0(v), where inv.S0 is the inverse function of S0(t).
#' Repeat N times.
#' 
#' Then the N samples {t_1, ..., t_N} follow the MCR dist above.
#' 
#' @param  n Number of samples following the MCR distribution
#' @param  p Cure rate parameter. When p = 0, it reduces to GMW(alpha, beta, gamma, lambda) distribution.
#' @param  alpha Generalized modified Weibull (GMW) distribution parameters. alpha > 0
#' @param  beta Generalized modified Weibull (GMW) distribution parameters. beta > 0
#' @param  gamma Generalized modified Weibull (GMW) distribution parameters. gamma >= 0 but gamma and lambda cannot be both 0.
#' @param  lambda Generalized modified Weibull (GMW) distribution parameters. lambda >= 0 but gamma and lambda cannot be both 0.
#' @param  tau  Threshold for delayed effect period. tau = 0 reduces no delayed effect.
#' @param  psi  Hazard ratio after delayed effect. psi = 1 reduces to the survival function without proportional hazards after delayed period (tau).
#' 
#' 
#' @examples
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#' set.seed(2022)
#' 
#' #(1) Exponential distribution mixture with cure rate 0.15
#' p = 0.15; alpha = log(2)/12
#' x = rmcr(n=10000, p=0.15, alpha = log(2)/12, beta=1, gamma=1, lambda=0)
#' m = qmcr(pct=0.5, p=0.15, alpha = log(2)/12, beta=1, gamma=1, lambda=0)
#' median(x)
#' km<-survival::survfit(survival::Surv(x,rep(1,length(x)))~rep(1,length(x)))
#' plot(km, ylab="Survival")
#' abline(h = c(p, 0.5), col="gray80", lty=2)
#' 
#' #(2) Weibull distribution mixture with cure rate 0.15
#' p = 0.15; alpha = log(2)/12
#' x = rmcr(n=10000, p=0.15, alpha = log(2)/12, beta=1, gamma=0.9, lambda=0)
#' m = qmcr(pct=0.5, p=0.15, alpha = log(2)/12, beta=1, gamma=0.9, lambda=0)
#' median(x)
#' 
#' #(3) Rayleigh distribution mixture with cure rate 0.15
#' p = 0.15; alpha = log(2)/100
#' x = rmcr(n=10000, p=0.15, alpha = alpha, beta=1, gamma=2, lambda=0)
#' m = qmcr(pct=0.5, p=0.15, alpha = alpha, beta=1, gamma=2, lambda=0)
#' median(x)
#' 
#' #(4) GMW distribution mixture with cure rate 0.15
#' p = 0.15; alpha = log(2)/12
#' x = rmcr(n=10000, p=0.15, alpha = alpha, beta=1, gamma=1, lambda=0.2)
#' m = qmcr(pct=0.5, p=0.15, alpha = alpha, beta=1, gamma=1, lambda=0.2)
#' median(x)
#' 
#' #(5) No cure rate
#' x = rmcr(n=100000, p=0, alpha = log(2)/12, beta=1, gamma=1, lambda=0)
#' median(x)
#' 
#' #(6) Delayed effect 6 months, proportional hazards HR = 0.6
#' set.seed(2022)
#' #control arm:
#' x0 = rmcr(n=1000, p=0.15, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=0, psi=1)
#' median(x0)
#' m0 = qmcr(pct=0.5, p=0.15, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=0, psi=1)
#' m0
#' 
#' #Experimental arm:
#' x1 = rmcr(n=1000, p=0.15, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=6, psi=0.5)
#' median(x1)
#' m1 = qmcr(pct=0.5, p=0.15, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=6, psi=0.5)
#' m1
#' x = c(x0, x1); group = c(rep(0,length(x0)),rep(1,length(x1)))
#' event = rep(1, length(x)) #all events
#' 
#' km<-survival::survfit(survival::Surv(x,event)~group)
#' plot(km, ylab="Survival")
#' abline(h = c(p, 0.5), col="gray80", lty=2)
#' abline(v=c(6, m1, m0), col="gray80", lty=2)
#' 
#' @export
#' 
rmcr = function(n=1, p=0.3, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=0, psi=1) {
  
  S0 = function(t){1-(1-exp(-alpha*t^gamma*exp(lambda*t)))^beta}
  S  = function(t){p + (1-p)*S0(t)}
  
  #Inverse function of S0(t)
  inv.S0 = function(s){
    phi = -log(1-(1-s)^(1/beta))/alpha
    if (lambda == 0){
      t = phi^{1/gamma}
    } else if(gamma == 0){
      t = phi^{1/lambda}
    } else{
      g = function(x){gamma*log(x) + lambda*x - log(phi)}
      t = uniroot(f=g, interval = c(0.0001, 1e10), tol=1e-3)$root
    }
  }
  s.tau = S(tau)
  u.UL = 1 - s.tau^(1-psi)*p^psi
  
  u = runif(n)
  u1 = u[u < u.UL]
  u2 = u[u >= u.UL]
  
  #For 1-u1 < S(tau), t = inv.S(1-u) = inv.S0((1-u-p)/(1-p))

  u1a = u1[1-u1 >= s.tau]
  u1b = u1[1-u1 < s.tau]
  
  va = (1 - u1a - p)/(1-p)
  t1a = sapply(va, inv.S0)
  
  sb = exp(log(1-u1b)/psi-(1-psi)/psi*log(s.tau))
  vb = (sb - p) / (1-p)
  #vb = ((1-u1b)^(1/psi) / (s.tau^((1-psi)/psi)) - p) / (1-p)
  t1b = sapply(vb, inv.S0)
  
  t2 = rep(Inf, length(u2))
  
  ans=NULL
  if (length(t1a) > 0){ans = c(ans, t1a)}
  if (length(t1b) > 0){ans = c(ans, t1b)}
  if (length(t2) > 0) {ans = c(ans, t2)}

  return(ans)
}

