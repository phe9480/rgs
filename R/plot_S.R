#' Display of Survival Curves per Study Design
#' 
#' This function plots the survival curves per study design
#' 
#' @param t  A sequence of survival time for the plot (x-axis)
#' @param S  An object containly a list of survival functions to draw
#' @param ... Other graphic parameters passed to the plot
#'
#' @return Display of the graph
#'  
#' @examples 
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' Hazard function and survival function for experimental arm
#' lambda1 = lambda0 * HR
#' h1 = function(t){lambda1}; S1= function(t){exp(-lambda1 * t)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' plotS(Tmax = 50, S = list(S0, S1), param=list(xlab="Time (mo)", ylab="PFS", main="PFS: HR = 0.65; Design Assumptions"))
#' 
#' @export
#' 
#' 
plot_S = function(S = list(S0, S1), Tmax = 50,
                  leg=list(x=30, y=1, txt=c("Control", "Exp. Arm")),
                  param = list(xlab="Survival Time (mo)", 
                             ylab="PFS",
                             main="PFS Survival Curve Per Study Design")){
  
  #number of survival curves
  g = length(S)
  t = seq(0, Tmax, by = 0.1)

  s0 = apply(cbind(t), MARGIN=1,FUN=S[[1]])
  
  plot(t, s0, type = "n", xlab = param$xlab, 
       ylab=param$ylab, main = param$main)
  for (i in 1:g){
    si = apply(cbind(t), MARGIN=1,FUN=S[[i]])
    lines(t, si, lwd = 3, col = i, lty = i)
  }
  abline(h = seq(0, 1, 0.1), col="gray80", lty=3)
  abline(v=seq(0, Tmax, by=2), col="gray80", lty=3)
  
  legend(leg$x, leg$y, leg$txt, col=1:g, lwd=2, lty=1:g, bty="n", cex=0.8)
}
