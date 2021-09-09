rdsensitivity <- function(y,x,c=0,fuzzy){
  
  Trt <- factor(x > c, exclude = NA)  
  
  require(rdrobust,knitr)
  
  if(missing(fuzzy)){
  
    if(FALSE){  
    b  <-  rdbwselect(y = y, x = x)
    b_opt  <-  b$bws[1]
    
    tri_opt  <-  pmax(0, 1 - abs(x) / b$bws[1])
    big_tri <- 2*tri_opt
    small_tri <- .5*tri_opt
    rect_opt  <-  ifelse(- b$bws[1] < x-c & x-c < b$bws[1], 1, 0)
    big_rect <- 2*rect_opt
    small_rect <- .5*rect_opt
    
    fit_tri_opt_linear  <-  lm(y ~ x+ Trt + x:Trt, 
                             subset = tri_opt > 0, weight = tri_opt)
    fit_tri_opt_quad  <-  lm(y ~ x + I(x^2)+ Trt + x:Trt + I(x^2):Trt, 
                           subset = tri_opt > 0, weight = tri_opt)
    fit_tri_opt_cubic  <-  lm(y ~ x + I(x^2) + I(x^3) + 
                              Trt + x:Trt + I(x^2):Trt + I(x^3):Trt, 
                            subset = tri_opt > 0, weight = tri_opt)
    
    fit_big_tri_linear  <-  lm(y ~ x+ Trt + x:Trt, 
                             subset = big_tri > 0, weight = big_tri)
    fit_big_tri_quad  <-  lm(y ~ x + I(x^2)+ Trt + x:Trt + I(x^2):Trt, 
                           subset = big_tri > 0, weight = big_tri)
    fit_big_tri_cubic  <-  lm(y ~ x + I(x^2) + I(x^3) + 
                              Trt + x:Trt + I(x^2):Trt + I(x^3):Trt, 
                            subset = big_tri > 0, weight = big_tri)
    
    fit_small_tri_linear  <-  lm(y ~ x+ Trt + x:Trt, 
                               subset = small_tri > 0, weight = small_tri)
    fit_small_tri_quad  <-  lm(y ~ x + I(x^2)+ Trt + x:Trt + I(x^2):Trt, 
                             subset = small_tri > 0, weight = small_tri)
    fit_small_tri_cubic  <-  lm(y ~ x + I(x^2) + I(x^3) + 
                                Trt + x:Trt + I(x^2):Trt + I(x^3):Trt, 
                              subset = small_tri > 0, weight = small_tri)
    
    fit_rect_opt_linear  <-  lm(y ~ x+ Trt + x:Trt, 
                              subset = rect_opt > 0, weight = rect_opt)
    fit_rect_opt_quad  <-  lm(y ~ x + I(x^2)+ Trt + x:Trt + I(x^2):Trt, 
                            subset = rect_opt > 0, weight = rect_opt)
    fit_rect_opt_cubic  <-  lm(y ~ x + I(x^2) + I(x^3) + 
                               Trt + x:Trt + I(x^2):Trt + I(x^3):Trt, 
                             subset = rect_opt > 0, weight = rect_opt)
    
    fit_big_rect_linear  <-  lm(y ~ x+ Trt + x:Trt, 
                              subset = big_rect > 0, weight = big_rect)
    fit_big_rect_quad  <-  lm(y ~ x + I(x^2)+ Trt + x:Trt + I(x^2):Trt, 
                            subset = big_rect > 0, weight = big_rect)
    fit_big_rect_cubic  <-  lm(y ~ x + I(x^2) + I(x^3) + 
                               Trt + x:Trt + I(x^2):Trt + I(x^3):Trt, 
                             subset = big_rect > 0, weight = big_rect)
    
    fit_small_rect_linear  <-  lm(y ~ x+ Trt + x:Trt, 
                                subset = small_rect > 0, weight = small_rect)
    fit_small_rect_quad  <-  lm(y ~ x + I(x^2)+ Trt + x:Trt + I(x^2):Trt, 
                              subset = small_rect > 0, weight = small_rect)
    fit_small_rect_cubic  <-  lm(y ~ x + I(x^2) + I(x^3) + 
                                 Trt + x:Trt + I(x^2):Trt + I(x^3):Trt, 
                               subset = small_rect > 0, weight = small_rect)
    
    results <- rbind(
      coeftest(fit_tri_opt_linear, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_tri_opt_quad, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_tri_opt_cubic, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_big_tri_linear, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_big_tri_quad, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_big_tri_cubic, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_small_tri_linear, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_small_tri_quad, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_small_tri_cubic, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_rect_opt_linear, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_rect_opt_quad, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_rect_opt_cubic, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_big_rect_linear, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_big_rect_quad, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_big_rect_cubic, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_small_rect_linear, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_small_rect_quad, vcov = vcovHC, type = "HC2")["TrtTRUE",],
      coeftest(fit_small_rect_cubic, vcov = vcovHC, type = "HC2")["TrtTRUE",])
    
    info <- cbind(c(rep("Triangular",9),rep("Rectangular",9)),rep(c(rep(b_opt,3),rep(2*b_opt,3),rep(.5*b_opt,3)),2),rep(c("Linear","Quadratic","Cubic"),6))
    
    colnames(info) <- c("Kernel","Bandwidth","Order")  
    info[,2] <- round(as.numeric(info[,2]),3)
    results <- round(results,3)
    knitr::kable(cbind(info,results[,1:4]),digits = 2,align='ccccccc')}
    
    fit_tri_opt_linear <- rdrobust(y, x,c ,p = 1,kernel = "triangular")
    fit_tri_opt_quad <- rdrobust(y=y, x=x,c=c ,p = 2,kernel = "triangular")
    fit_tri_opt_cubic <- rdrobust(y=y, x=x,c=c ,p = 3,kernel = "triangular")
    
    fit_big_tri_linear <- rdrobust(y=y, x=x,c=c ,p = 1,kernel = "triangular",h= 2*fit_tri_opt_linear$bws[1,])
    fit_big_tri_quad <- rdrobust(y=y, x=x,c=c ,p = 2,kernel = "triangular", h= 2* fit_tri_opt_quad$bws[1,])
    fit_big_tri_cubic <- rdrobust(y=y, x=x,c=c ,p = 3,kernel = "triangular", h= 2* fit_tri_opt_cubic$bws[1,])
    
    fit_small_tri_linear <- rdrobust(y=y, x=x,c=c ,p = 1,kernel = "triangular", h= .5* fit_tri_opt_linear$bws[1,])
    fit_small_tri_quad <- rdrobust(y=y, x=x,c=c ,p = 2,kernel = "triangular", h= .5* fit_tri_opt_quad$bws[1,])
    fit_small_tri_cubic <- rdrobust(y=y, x=x,c=c ,p = 3,kernel = "triangular", h= .5* fit_tri_opt_cubic$bws[1,])
    
    fit_rect_opt_linear <- rdrobust(y=y, x=x,c=c ,p = 1,kernel = "uniform")
    fit_rect_opt_quad <- rdrobust(y=y, x=x,c=c ,p = 2,kernel = "uniform")
    fit_rect_opt_cubic <- rdrobust(y=y, x=x,c=c ,p = 3,kernel = "uniform")
    
    fit_big_rect_linear <- rdrobust(y=y, x=x,c=c ,p = 1,kernel = "uniform", h= 2* fit_rect_opt_linear$bws[1,])
    fit_big_rect_quad <- rdrobust(y=y, x=x,c=c ,p = 2,kernel = "uniform", h= 2* fit_rect_opt_quad$bws[1,])
    fit_big_rect_cubic <- rdrobust(y=y, x=x,c=c ,p = 3,kernel = "uniform", h= 2* fit_rect_opt_cubic$bws[1,])
    
    fit_small_rect_linear <- rdrobust(y=y, x=x,c=c ,p = 1,kernel = "uniform", h= .5* fit_rect_opt_linear$bws[1,])
    fit_small_rect_quad <- rdrobust(y=y, x=x,c=c ,p = 2,kernel = "uniform",, h= .5* fit_rect_opt_quad$bws[1,])
    fit_small_rect_cubic <- rdrobust(y=y, x=x,c=c ,p = 3,kernel = "uniform", h= .5* fit_rect_opt_cubic$bws[1,])
    
    results <- rbind(fit_tri_opt_linear$Estimate,
                   fit_tri_opt_quad$Estimate,
                   fit_tri_opt_cubic$Estimate,
                   fit_big_tri_linear$Estimate,
                   fit_big_tri_quad$Estimate,
                   fit_big_tri_cubic$Estimate,
                   fit_small_tri_linear$Estimate,
                   fit_small_tri_quad$Estimate,
                   fit_small_tri_cubic$Estimate,
                   fit_rect_opt_linear$Estimate,
                   fit_rect_opt_quad$Estimate,
                   fit_rect_opt_cubic$Estimate,
                   fit_big_rect_linear$Estimate,
                   fit_big_rect_quad$Estimate,
                   fit_big_rect_cubic$Estimate,
                   fit_small_rect_linear$Estimate,
                   fit_small_rect_quad$Estimate,
                   fit_small_rect_cubic$Estimate)
    
    resultsp1 <- rbind(fit_tri_opt_linear$pv[1],
                     fit_tri_opt_quad$pv[1],
                     fit_tri_opt_cubic$pv[1],
                     fit_big_tri_linear$pv[1],
                     fit_big_tri_quad$pv[1],
                     fit_big_tri_cubic$pv[1],
                     fit_small_tri_linear$pv[1],
                     fit_small_tri_quad$pv[1],
                     fit_small_tri_cubic$pv[1],
                     fit_rect_opt_linear$pv[1],
                     fit_rect_opt_quad$pv[1],
                     fit_rect_opt_cubic$pv[1],
                     fit_big_rect_linear$pv[1],
                     fit_big_rect_quad$pv[1],
                     fit_big_rect_cubic$pv[1],
                     fit_small_rect_linear$pv[1],
                     fit_small_rect_quad$pv[1],
                     fit_small_rect_cubic$pv[1])
    
    resultsp2 <- rbind(fit_tri_opt_linear$pv[2],
                     fit_tri_opt_quad$pv[2],
                     fit_tri_opt_cubic$pv[2],
                     fit_big_tri_linear$pv[2],
                     fit_big_tri_quad$pv[2],
                     fit_big_tri_cubic$pv[2],
                     fit_small_tri_linear$pv[2],
                     fit_small_tri_quad$pv[2],
                     fit_small_tri_cubic$pv[2],
                     fit_rect_opt_linear$pv[2],
                     fit_rect_opt_quad$pv[2],
                     fit_rect_opt_cubic$pv[2],
                     fit_big_rect_linear$pv[2],
                     fit_big_rect_quad$pv[2],
                     fit_big_rect_cubic$pv[2],
                     fit_small_rect_linear$pv[2],
                     fit_small_rect_quad$pv[2],
                     fit_small_rect_cubic$pv[2])
    
    resultsp3 <- rbind(fit_tri_opt_linear$pv[3],
                     fit_tri_opt_quad$pv[3],
                     fit_tri_opt_cubic$pv[3],
                     fit_big_tri_linear$pv[3],
                     fit_big_tri_quad$pv[3],
                     fit_big_tri_cubic$pv[3],
                     fit_small_tri_linear$pv[3],
                     fit_small_tri_quad$pv[3],
                     fit_small_tri_cubic$pv[3],
                     fit_rect_opt_linear$pv[3],
                     fit_rect_opt_quad$pv[3],
                     fit_rect_opt_cubic$pv[3],
                     fit_big_rect_linear$pv[3],
                     fit_big_rect_quad$pv[3],
                     fit_big_rect_cubic$pv[3],
                     fit_small_rect_linear$pv[3],
                     fit_small_rect_quad$pv[3],
                     fit_small_rect_cubic$pv[3])
    
    results <- cbind(results,resultsp1,resultsp2,resultsp3)
    
    info <- cbind(c(rep("Triangular",9),rep("Rectangular",9)),rep(c(rep(b_opt,3),rep(2*b_opt,3),rep(.5*b_opt,3)),2),rep(c("Linear","Quadratic","Cubic"),6))
    
    colnames(info) <- c("Kernel","Bandwidth","Order")  
    info[,2] <- round(as.numeric(info[,2]),3)
    results <- round(results,3)
    colnames(results) <- c("Estimate", "Estimate (bc)", "SE","SE (robust)","P>|z|","P>|z| (bc)","P>|z| (robust)")
    knitr::kable(cbind(info,results),digits = 2, align='cccccccccc')
    
  }
  
  else {
    
    
    fit_tri_opt_linear <- rdrobust(y, x,c,fuzzy ,p = 1,kernel = "triangular")
    fit_tri_opt_quad <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 2,kernel = "triangular")
    fit_tri_opt_cubic <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 3,kernel = "triangular")
    
    fit_big_tri_linear <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 1,kernel = "triangular",h= 2*fit_tri_opt_linear$bws[1,])
    fit_big_tri_quad <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 2,kernel = "triangular", h= 2* fit_tri_opt_quad$bws[1,])
    fit_big_tri_cubic <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 3,kernel = "triangular", h= 2* fit_tri_opt_cubic$bws[1,])
    
    fit_small_tri_linear <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 1,kernel = "triangular", h= .5* fit_tri_opt_linear$bws[1,])
    fit_small_tri_quad <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 2,kernel = "triangular", h= .5* fit_tri_opt_quad$bws[1,])
    fit_small_tri_cubic <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 3,kernel = "triangular", h= .5* fit_tri_opt_cubic$bws[1,])
    
    fit_rect_opt_linear <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 1,kernel = "uniform")
    fit_rect_opt_quad <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 2,kernel = "uniform")
    fit_rect_opt_cubic <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 3,kernel = "uniform")
    
    fit_big_rect_linear <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 1,kernel = "uniform", h= 2* fit_rect_opt_linear$bws[1,])
    fit_big_rect_quad <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 2,kernel = "uniform", h= 2* fit_rect_opt_quad$bws[1,])
    fit_big_rect_cubic <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 3,kernel = "uniform", h= 2* fit_rect_opt_cubic$bws[1,])
    
    fit_small_rect_linear <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 1,kernel = "uniform", h= .5* fit_rect_opt_linear$bws[1,])
    fit_small_rect_quad <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 2,kernel = "uniform",, h= .5* fit_rect_opt_quad$bws[1,])
    fit_small_rect_cubic <- rdrobust(y=y, x=x,c=c,fuzzy ,p = 3,kernel = "uniform", h= .5* fit_rect_opt_cubic$bws[1,])
    
    results <- rbind(fit_tri_opt_linear$Estimate,
                   fit_tri_opt_quad$Estimate,
                   fit_tri_opt_cubic$Estimate,
                   fit_big_tri_linear$Estimate,
                   fit_big_tri_quad$Estimate,
                   fit_big_tri_cubic$Estimate,
                   fit_small_tri_linear$Estimate,
                   fit_small_tri_quad$Estimate,
                   fit_small_tri_cubic$Estimate,
                   fit_rect_opt_linear$Estimate,
                   fit_rect_opt_quad$Estimate,
                   fit_rect_opt_cubic$Estimate,
                   fit_big_rect_linear$Estimate,
                   fit_big_rect_quad$Estimate,
                   fit_big_rect_cubic$Estimate,
                   fit_small_rect_linear$Estimate,
                   fit_small_rect_quad$Estimate,
                   fit_small_rect_cubic$Estimate)
    
    resultsp1 <- rbind(fit_tri_opt_linear$pv[1],
                     fit_tri_opt_quad$pv[1],
                     fit_tri_opt_cubic$pv[1],
                     fit_big_tri_linear$pv[1],
                     fit_big_tri_quad$pv[1],
                     fit_big_tri_cubic$pv[1],
                     fit_small_tri_linear$pv[1],
                     fit_small_tri_quad$pv[1],
                     fit_small_tri_cubic$pv[1],
                     fit_rect_opt_linear$pv[1],
                     fit_rect_opt_quad$pv[1],
                     fit_rect_opt_cubic$pv[1],
                     fit_big_rect_linear$pv[1],
                     fit_big_rect_quad$pv[1],
                     fit_big_rect_cubic$pv[1],
                     fit_small_rect_linear$pv[1],
                     fit_small_rect_quad$pv[1],
                     fit_small_rect_cubic$pv[1])
    
    resultsp2 <- rbind(fit_tri_opt_linear$pv[2],
                     fit_tri_opt_quad$pv[2],
                     fit_tri_opt_cubic$pv[2],
                     fit_big_tri_linear$pv[2],
                     fit_big_tri_quad$pv[2],
                     fit_big_tri_cubic$pv[2],
                     fit_small_tri_linear$pv[2],
                     fit_small_tri_quad$pv[2],
                     fit_small_tri_cubic$pv[2],
                     fit_rect_opt_linear$pv[2],
                     fit_rect_opt_quad$pv[2],
                     fit_rect_opt_cubic$pv[2],
                     fit_big_rect_linear$pv[2],
                     fit_big_rect_quad$pv[2],
                     fit_big_rect_cubic$pv[2],
                     fit_small_rect_linear$pv[2],
                     fit_small_rect_quad$pv[2],
                     fit_small_rect_cubic$pv[2])
    
    resultsp3 <- rbind(fit_tri_opt_linear$pv[3],
                     fit_tri_opt_quad$pv[3],
                     fit_tri_opt_cubic$pv[3],
                     fit_big_tri_linear$pv[3],
                     fit_big_tri_quad$pv[3],
                     fit_big_tri_cubic$pv[3],
                     fit_small_tri_linear$pv[3],
                     fit_small_tri_quad$pv[3],
                     fit_small_tri_cubic$pv[3],
                     fit_rect_opt_linear$pv[3],
                     fit_rect_opt_quad$pv[3],
                     fit_rect_opt_cubic$pv[3],
                     fit_big_rect_linear$pv[3],
                     fit_big_rect_quad$pv[3],
                     fit_big_rect_cubic$pv[3],
                     fit_small_rect_linear$pv[3],
                     fit_small_rect_quad$pv[3],
                     fit_small_rect_cubic$pv[3])
    
    results <- cbind(results,resultsp1,resultsp2,resultsp3)
    
    info <- cbind(c(rep("Triangular",9),rep("Rectangular",9)),rep(c(rep(b_opt,3),rep(2*b_opt,3),rep(.5*b_opt,3)),2),rep(c("Linear","Quadratic","Cubic"),6))
    
    colnames(info) <- c("Kernel","Bandwidth","Order")  
    info[,2] <- round(as.numeric(info[,2]),3)
    results <- round(results,3)
    colnames(results) <- c("Estimate", "Estimate (bc)", "SE","SE (robust)","P>|z|","P>|z| (bc)","P>|z| (robust)")
    knitr::kable(cbind(info,results),digits = 2, align='cccccccccc')
  }
}