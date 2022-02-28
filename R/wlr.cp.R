#'  Conditional Power When IA Negative And Benefit-Risk Assessment Compared to Standard Logrank Test
#' 
#' @param b   A vector of rejection boundary in Z for IA and FA
#' @param blr Rejection boundary at FA using logrank test at IA and FA
#' @param mu1 Mean(Z1): IA Z
#' @param cv1 Cov(Z1)
#' @param mu  Mean(Z1, Zf): IA Z and FA Z
#' @param cv  Cov(Z1, Zf)
#' @param mu1lr Mean(Z1, Zlr_f): IA Z1 and logrank Z at FA
#' @param cv1lr Cov(Z1, Zlr_f)
#' @param mulr  Mean(Z1, Zlr_f, Zf)
#' @param cvlr  Cov(Z1, Zlr_f, Zf)
#' @param lrInFA FA includes logrank test component. "Y", "N"
#' @return 
#' (1) Conditional power when IA is negative;
#' (2) Benefit: P(logrank neg at FA; but weighted logrank test negative, IA negative)
#' (3) Risk: P(logrank pos if used at FA; but weighted logrank #' (3) Risk: #' (3) Risk: )#'   

#' @examples 
#' 
#' DCO = c(24, 36); 
#' h0 = function(t){log(2)/12}; S0= function(t){exp(-log(2)/12 * t)};
#' h1 = function(t){log(2)/12*0.70}; S1= function(t){exp(-log(2)/12 * 0.7 * t)}; 
#' primary.test = list(IA1=list(function(s){1}), FA=list(function(s){1}, function(s){s*(1-s)})); 
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)};
#' G0 = function(t){0}; G1 = function(t){0}; 
#'
#' pow = wlr.power.maxcombo(DCO = DCO, alpha=c(0.008, 0.01699),
#' r = 1, n = 450, h0 = h0, S0=S0, 
#' h1 = h1, S1 = S1, f.ws = primary.test, 
#' Lambda=Lambda, G0=G0, G1=G1,
#' mu.method = "H1", cov.method = "H1.LA")
#' 
#' b = pow$design$b
#' mu1 = pow$mu[1,!is.na(pow$mu[1,])]
#' cv1 = pow$Omega1[1,1]
#' mu = NULL
#' for(i in 1:nrow(pow$mu)){mu=c(mu,pow$mu[i,!is.na(pow$mu[i,])])}
#' cv = pow$Omega1
#' 
#' lr.test = list(IA1=list(function(s){1}), FA=list(function(s){1})); 
#' 
#' pow.lr = wlr.power.maxcombo(DCO = DCO, alpha=c(0.008, 0.01699),
#' r = 1, n = 450, h0 = h0, S0=S0, 
#' h1 = h1, S1 = S1, f.ws = lr.test, 
#' Lambda=Lambda, G0=G0, G1=G1,
#' mu.method = "H1", cov.method = "H1.LA")
#' 
#' b.lrf = lr.pow$design$b[2]
#' mu1lr = c(pow$mu[1,1],lr.pow$mu[2])   #Mean(Z1, Zlr_f)
#' cv1lr = lr.pow$Omega1                 #cov(Z1, Zlr_f)
#' mulr  = mu    #Mean(Z1, Zlr_f, Zf)
#' cvlr  = cv    #cov(Z1, Zlr_f, Zf). Note Zf already includes Zlr_f in this case
#' wlr.cp (b = b, b.lrf=b.lrf, mu1=mu1, cv1=cv1, mu=mu, cv=cv, 
#'         mu1lr=mu1lr, cv1lr=cv1lr, mulr=mulr, cvlr=cvlr)
#'                  
#' @export
#' 
#' 
wlr.cp = function(b = c(2.3, 2.0), b.lrf, mu1, cv1, mu, cv, 
                  mu1lr, cv1lr, mulr, cvlr, lrInFA="Y"){
  
  #(1) CP(FA positive | IA negative) = P(Zf>bf|Z1<b1)
  #   = 1 - P(Z1<b1, Zf<bf)/P(Z1<b1)
  L1 = length(mu1)
  Lf = length(mu) - length(mu1)
  
  ##P(Z1<b1) 
  if(L1 == 1){
    z1 = pnorm(b[1], mean=mu1, sd=sqrt(cv1))
  } else {
    z1 = mvtnorm::pmvnorm(lower=rep(-Inf,L1), upper=rep(b[1], L1),
                            mean=mu1, sigma = cv1, 
                            abseps = 1e-8, maxpts=100000)[1]
  }
  ##Prob. P(Z1<b1, Zf<bf)
  z1.zf = mvtnorm::pmvnorm(lower=rep(-Inf,L1+Lf),upper = c(rep(b[1],L1),rep(b[2],Lf)),
                           mean=mu, sigma = cv, 
                           abseps = 1e-8, maxpts=100000)[1]
  ##Prob. FA Pos
  cp = 1 - z1.zf / z1
  
  #(2) P(FA positive, logrank neg, IA negative) 
  #   = P(Zf > bf, Zlr_f < b.lrf, Z1 < b1)

  #Need covariance matrix of (Z1, Zlr_f, Zf)
  #eg, IA: lr; FA: max(lr, G01, G11), then the covariance matrix 
  #produced by wlr.power.maxcombo is adequate. However, if 
  #IA: lr; FA: G11; then needs to construct the covariance matrix by creating
  #a test strategy: test.tmp = list(IA=list(lr), FA=list(lr, G11)), 
  #then call wlr.power.maxcomb to create the covariance matrix.

  #P(Z1<b1, Zlr_f < b.lrf)
  z1.zlrf = mvtnorm::pmvnorm(lower=rep(-Inf,L1+1), 
                               upper=c(rep(b[1], L1), b.lrf),
                               mean=mu1lr, sigma = cv1lr,
                               abseps = 1e-8, maxpts=100000)[1]    
  
  #P(Z1<b1, Zlr_f<b.lrf, Zf<bf)
  if (lrInFA == "Y"){
    z1.zlrf.zf = mvtnorm::pmvnorm(lower=rep(-Inf,L1+Lf), 
                                  upper=c(rep(b[1], L1), min(b.lrf,b[2]), rep(b[2], Lf-1)),
                                  mean=mulr, sigma = cvlr,
                                  abseps = 1e-8, maxpts=100000)[1]  
  } else {
    z1.zlrf.zf = mvtnorm::pmvnorm(lower=rep(-Inf,L1+Lf+1), 
                                  upper=c(rep(b[1], L1), b.lrf, rep(b[2], Lf)),
                                  mean=mulr, sigma = cvlr,
                                  abseps = 1e-8, maxpts=100000)[1]    
  }
  p.benefit = (z1.zlrf - z1.zlrf.zf)

  #(3) P(FA neg, logrank pos if logrank used for FA, IA negative) 
  #  = P(Z1<b1, Zlr_f<blr_f, Zf<bf)
  
  p.risk = (z1.zf - z1.zlrf.zf)
  o = list()
  o$cp = cp; o$p.benefit = p.benefit; o$p.risk = p.risk
  
  return(o)
}

  
  
  