library(mvtnorm)
Datageneration = function(n, beta, corr, Dist = 'mzNormal', link = 'mod1')
{
   d <- length(beta) - 2
   sigmax <- matrix(corr, d-1, d-1) + diag(2-corr, d-1)
   if( Dist == 'mzNormal' ){X  <- rmvnorm(n, rep(0, d-1), sigmax)}
   if( Dist == 'T5' ){X  <- rmvt(n, sigma = sigmax, df = 5)}
   if( Dist == 'T3' ){X  <- rmvt(n, sigma = sigmax, df = 3)}
   X =cbind(1, X,  runif(n, 0, 0.5), runif(n, 0, 0.5))
   
   #X = cbind(1, X)
   p <- d + 1
   if(link == 'mod1')
   {
      prob = 1 - 1 / (1 + exp(c(X %*% beta)))  
   }
   if(link == 'mod2')
   {
     prob = 1 - 1 / (1 + exp(c(X %*% beta) + 0.5* sin(0.5*X[,2]) - 0.5* sin(0.5*X[,3]) + 0.2* sin(0.2*X[,8]))) 
   }
   if(link == 'mod3')
   {
     prob = 1 - 1 / (1 + exp(c(X %*% beta) + 0.5 * X[,6]^2 - 0.5 * X[,8]^2 + exp( X[,2:3]%*%rep(0.5,2) ))) 
   }
   Y = rbinom(n,1,prob)
   list(Y = Y, X = X)
}
##############active entropy####
SES <- function(pool.x, pool.y, sub.x, sub.y, sub.p, batch, c.beta, i.model, setup = FALSE, method = 'Self',link = 'mod1') 
{
   if(setup == TRUE)
   {
      pool.x = pool.x; pool.y = pool.y; sub.x = sub.x; sub.y=sub.y
      sub.p = sub.p; c.beta = c.beta; pi=NULL; i.model = i.model
   }else{
      full.x = rbind(sub.x, pool.x)
      if(link == 'mod2')
      {
         imp.v = predict(i.model,newdata=data.frame(x2=pool.x[,2],x3=pool.x[,3],x4=pool.x[,4],x5=pool.x[,5],x6=pool.x[,6],
                                               x7=pool.x[,7],x8=pool.x[,8],x9=sin(0.5*pool.x[,2]),
                                               x10=sin(0.5*pool.x[,3]),x11=sin(0.2*pool.x[,8])),type="response")
      }else if(link == 'mod3'){
         imp.v = predict(i.model,newdata=data.frame(x2=pool.x[,2],x3=pool.x[,3],x4=pool.x[,4],x5=pool.x[,5],x6=pool.x[,6],
                                               x7=pool.x[,7],x8=pool.x[,8],x9=pool.x[,6]^2,
                                               x10=pool.x[,8]^2,x11=exp(pool.x[,2]*0.3),x12=exp( pool.x[,3]*0.3 )),type="response")
      }else if(link == 'Single'){
         imp.v = predict(i.model,newdata=data.frame(x2=pool.x[,2],x3=pool.x[,3],x4=pool.x[,4],x5=pool.x[,5],x6=pool.x[,6],
                                                    x7=pool.x[,7],x8=pool.x[,8]),type="response")
      }else{ 
        imp.v = predict(i.model,newdata=data.frame(x2=pool.x[,2],x3=pool.x[,3],x4=pool.x[,4],x5=pool.x[,5],x6=pool.x[,6],
                                     x7=pool.x[,7],x8=pool.x[,8]),type="response")
      }  
      c.pool.n = length(pool.y); t.n = c.pool.n + length(sub.y)
      score = full.x %*% c.beta
      p_s = c(1 - 1 / (1 + exp(score)))
      w.s <- p_s * (1 - p_s)
      W.s <- solve(t(full.x) %*% (full.x * w.s))
      #if(method == 'Mean')
      #{
      #   wegt.s = sqrt( ((1-imp.v)*imp.v)^2 * rowSums((pool.x%*%W.s)^2))
      #}   
      wegt.s = sqrt((imp.v - p_s[-(1:length(sub.y))])^2 * rowSums((pool.x%*%W.s)^2))
      if( sum(is.na(wegt.s)) != 0 )
      {
        wegt.s[is.na(wegt.s)] = min(wegt.s[!is.na(wegt.s)] ) 
      }
      ses.prob = wegt.s/sum(wegt.s)
      
      ind.x.ses = sample(1:c.pool.n, size = batch, replace = FALSE, prob = ses.prob)
      sub.x = rbind(sub.x, pool.x[ind.x.ses,])
      sub.y = c(sub.y, pool.y[ind.x.ses])
      pool.x = pool.x[-ind.x.ses,]
      pool.y = pool.y[-ind.x.ses]
      
      numer = wegt.s[ind.x.ses]
      deno = rev(cumsum(c(sum(wegt.s[-ind.x.ses]), rev(numer)))[-1])
      sub.p = c(sub.p,  numer/deno)
      first = (t.n-length(sub.p))/(t.n-1:length(sub.p))
      second = (1/((t.n-(1:length(sub.p))+1)*sub.p) - 1)
      pi <- 1 + first * second
      c.beta = c.beta
   }   
   return(list(pool.x = pool.x, pool.y = pool.y, sub.x = sub.x, sub.y=sub.y, pi = pi,
               sub.p = sub.p, c.est = c.beta, i.model=i.model) )
}
##############active entropy####
entropy <- function(pool.x, pool.y, sub.x, sub.y, sub.p, batch, c.beta, setup = FALSE) 
{
   if(setup == TRUE)
   {
      pool.x = pool.x; pool.y = pool.y; sub.x = sub.x; sub.y=sub.y
      sub.p = sub.p; c.beta = c.beta; pi=NULL 
   }else{
      c.pool.n = length(pool.y); t.n = c.pool.n + length(sub.y)
      score = pool.x %*% c.beta
      p_e = c(1 - 1 / (1 + exp(score)))
      wegt.e = p_e* log(p_e) + (1-p_e) * log(1-p_e)
      if( sum(is.na(wegt.e)) != 0 )
      {
        wegt.e[is.na(wegt.e)] = min(wegt.e[!is.na(wegt.e)] ) 
       }
      en.prob = wegt.e/sum(wegt.e)
      
      ind.x.en = sample(1:c.pool.n, size = batch, replace = FALSE, prob = en.prob)
      sub.x = rbind(sub.x, pool.x[ind.x.en,])
      sub.y = c(sub.y, pool.y[ind.x.en])
      pool.x = pool.x[-ind.x.en,]
      pool.y = pool.y[-ind.x.en]
      
      numer = wegt.e[ind.x.en]
      deno = rev(cumsum(c(sum(wegt.e[-ind.x.en]), rev(numer)))[-1])
      sub.p = c(sub.p,  numer/deno)
      first = (t.n-length(sub.p))/(t.n-1:length(sub.p))
      second = (1/((t.n-(1:length(sub.p))+1)*sub.p) - 1)
      pi <- 1 + first * second
      c.beta = c.beta
   }   
   return(list(pool.x = pool.x, pool.y = pool.y, sub.x = sub.x, sub.y=sub.y, pi = pi,
               sub.p = sub.p, c.est = c.beta) )
}
##############active entropy####
LC <- function(pool.x, pool.y, sub.x, sub.y, sub.p, batch, c.beta, setup = FALSE) 
{
  if(setup == TRUE)
  {
    pool.x = pool.x; pool.y = pool.y; sub.x = sub.x; sub.y=sub.y
    sub.p = sub.p; c.beta = c.beta; pi=NULL 
  }else{
    c.pool.n = length(pool.y); t.n = c.pool.n + length(sub.y)
    score = pool.x %*% c.beta
    p_e = c(1 - 1 / (1 + exp(score)))

    wegt.e = 1-p_e
    ifelse(wegt.e ==0, min(wegt.e[wegt.e!=0])/10^2 ,wegt.e)
    if( sum(is.na(wegt.e)) != 0 )
    {
      wegt.e[is.na(wegt.e)] = min(wegt.e[!is.na(wegt.e)] ) 
    }
    en.prob = wegt.e/sum(wegt.e)
    
    ind.x.en = sample(1:c.pool.n, size = batch, replace = FALSE, prob = en.prob)
    sub.x = rbind(sub.x, pool.x[ind.x.en,])
    sub.y = c(sub.y, pool.y[ind.x.en])
    pool.x = pool.x[-ind.x.en,]
    pool.y = pool.y[-ind.x.en]
    
    numer = wegt.e[ind.x.en]
    deno = rev(cumsum(c(sum(wegt.e[-ind.x.en]), rev(numer)))[-1])
    sub.p = c(sub.p,  numer/deno)
    first = (t.n-length(sub.p))/(t.n-1:length(sub.p))
    second = (1/((t.n-(1:length(sub.p))+1)*sub.p) - 1)
    pi <- 1 + first * second
    c.beta = c.beta
    
  }   
  return(list(pool.x = pool.x, pool.y = pool.y, sub.x = sub.x, sub.y=sub.y, pi = pi,
              sub.p = sub.p, c.est = c.beta) )
}
###################
Uni <- function(pool.x, pool.y, sub.x, sub.y, sub.p, batch, c.est, setup = FALSE) {
  if(setup == TRUE)
  {
    pool.x = pool.x; pool.y = pool.y; sub.x = sub.x; sub.y=sub.y
    sub.p = sub.p; c.est = c.est; pi=NULL 
  }else{
    c.pool.n = length(pool.y); t.n = c.pool.n + length(sub.y)
    ind.x.u = sample(1:c.pool.n, size = batch, replace = FALSE)
    sub.x = rbind(sub.x, pool.x[ind.x.u,]);
    sub.y = c(sub.y, pool.y[ind.x.u])
    sub.p = c(sub.p, 1/(c.pool.n:1)[1:batch])
    pool.x = pool.x[-ind.x.u,]; pool.y = pool.y[-ind.x.u]
    
    first = (t.n-length(sub.p))/(t.n-1:length(sub.p))
    second = (1/((t.n-(1:length(sub.p))+1)*sub.p) - 1)
    pi <- 1 + first * second
    c.est = c.est;
  }  
  return(list(pool.x = pool.x, pool.y = pool.y, sub.x = sub.x, pi = pi,
              sub.y=sub.y, sub.p = sub.p, c.est = c.est))
}
###############
m.t <- function(sub.y, sub.x, res.x, pi, deg, method = 'IPW', design, link = 'mod2'){
   cur.n = length(sub.y); rest.n = dim(res.x)[1]; t.x = rbind(sub.x, res.x)
   if(method == 'No'){
     cur.est =  getest(sub.x, sub.y)$par 
     iw = pi; reg =NULL
   }
   if(method == 'IPW'){
      cur.est =  getest(sub.x, sub.y, w = pi)$par 
      iw = pi; reg =NULL
   }
   if(method == 'IMP'){
     x2 = sub.x[,2];x3 = sub.x[,3];x4 = sub.x[,4];x5 = sub.x[,5];x6 = sub.x[,6];x7 = sub.x[,7];
     x8 = sub.x[,8]
     reg = glm(sub.y~1+ns(x2,df=deg)+ns(x3,df=deg)+ns(x4,df=deg)+ns(x5,df=deg)+
               ns(x6,df=deg)+ns(x7,df=deg)+ns(x8,df=deg), family=binomial)
     Imp.y = predict(reg,newdata=data.frame(x2=t.x[,2],x3=t.x[,3],x4=t.x[,4],x5=t.x[,5],x6=t.x[,6],
                                            x7=t.x[,7],x8=t.x[,8]),type="response")
     cur.est = getest.aug.now(t.x, sub.y, Imp.y, w = pi)$par
     iw = pi
   }
   if(method == 'Single'){
     x2 = sub.x[,2];x3 = sub.x[,3];x4 = sub.x[,4];x5 = sub.x[,5];x6 = sub.x[,6];x7 = sub.x[,7];
     x8 = sub.x[,8]
     reg = glm(sub.y~1+x2+x3+ns(x4,df=deg)+x5+
                  x6+x7+x8, family=binomial)
     Imp.y = predict(reg,newdata=data.frame(x2=t.x[,2],x3=t.x[,3],x4=t.x[,4],x5=t.x[,5],x6=t.x[,6],
                                            x7=t.x[,7],x8=t.x[,8]),type="response")
     cur.est = getest.aug.now(t.x, sub.y, Imp.y, w = pi)$par
     iw = pi
   }
   if(method == 'Perfect'){
     if(link == 'mod2')
     {
       x2 = sub.x[,2];x3 = sub.x[,3];x4 = sub.x[,4];x5 = sub.x[,5];x6 = sub.x[,6];x7 = sub.x[,7];
       x8 = sub.x[,8]; x9= sin(0.5*sub.x[,2]); x10= sin(0.5*sub.x[,3]); x11= sin(0.2*sub.x[,8])
       reg = glm(sub.y~1+ns(x2,df=deg)+ns(x3,df=deg)+ns(x4,df=deg)+ns(x5,df=deg)+
                   ns(x6,df=deg)+ns(x7,df=deg)+ns(x8,df=deg)+ns(x9,df=deg)+ns(x10,df=deg)+ns(x11,df=deg), family=binomial)
       Imp.y = predict(reg,newdata=data.frame(x2=t.x[,2],x3=t.x[,3],x4=t.x[,4],x5=t.x[,5],x6=t.x[,6],
                                              x7=t.x[,7],x8=t.x[,8],x9=sin(0.5*t.x[,2]),
                                              x10=sin(0.5*t.x[,3]),x11=sin(0.2*t.x[,8])),type="response")
       cur.est = getest.aug.now(t.x, sub.y, Imp.y, w = pi)$par
       iw = pi
     }
     if(link == 'mod3')
     {
       x2 = sub.x[,2];x3 = sub.x[,3];x4 = sub.x[,4];x5 = sub.x[,5];x6 = sub.x[,6];x7 = sub.x[,7];
       x8 = sub.x[,8]; x9= sub.x[,6]^2; x10= sub.x[,8]^2; x11= exp( sub.x[,2]*0.3 ); x12= exp( sub.x[,3]*0.3 )
       reg = glm(sub.y~1+ns(x2,df=deg)+ns(x3,df=deg)+ns(x4,df=deg)+ns(x5,df=deg)+ns(x6,df=deg)+
                   ns(x7,df=deg)+ns(x8,df=deg)+ns(x9,df=deg)+ns(x10,df=deg)+ns(x11,df=deg)+ns(x12,df=deg), family=binomial)
       Imp.y = predict(reg,newdata=data.frame(x2=t.x[,2],x3=t.x[,3],x4=t.x[,4],x5=t.x[,5],x6=t.x[,6],
                                              x7=t.x[,7],x8=t.x[,8],x9=t.x[,6]^2,
                                              x10=t.x[,8]^2,x11=exp(t.x[,2]*0.3),x12=exp( t.x[,3]*0.3 )),type="response")
       cur.est = getest.aug.now(t.x, sub.y, Imp.y, w = pi)$par
       iw = pi
     }
   }
   return(list(cur.est=cur.est, iw = iw, i.model = reg))
}
###############
thres.inx <- function(e.proba, y, const = 0.3){
   thres = seq(0,1,0.0005)
   class.y = t(sapply(e.proba, function(x) ifelse(x>thres,1,0)))
   fpr.est.M = colMeans(class.y*!y)/(1-mean(y))
   inx = which(fpr.est.M <= const)[1]
   c.thres = thres[inx]
   return(c.thres)
}
#############################
getest <- function(x, y, w=1) {
  d = ncol(x)
  beta <- rep(0, d)
  dig <- diag(0.00001, d)
  loop  <- 1
  Loop  <- 100
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1 - 1 / (1 + exp(x %*% beta)))
    H <- t(x) %*% (pr * (1 - pr) * w * x)
    S <- colSums((y - pr) * w * x)
    tryCatch(
      {shs <- NA
      shs <- solve(H + dig, S) }, 
      error=function(e){
        cat("\n ERROR :", loop, conditionMessage(e), "\n")})
    if (is.na(shs[1])) {
      msg <- "Not converge"
      beta <- loop <- NA
      break
    }
    beta.new <- beta + shs
    tlr  <- sum((beta.new - beta)^2)
    beta  <- beta.new
    if(tlr < 0.000001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop) {
      warning("Maximum iteration reached")
      ## beta <- NA
    }
    loop  <- loop + 1
  }
  list(par=beta, message=msg, iter=loop)
}
#############################
getest.aug.now <- function(x, y, imp.y, w=1) {
   d = ncol(x)
   beta <- rep(0, d)
   dig <- diag(0.00001, d)
   loop  <- 1
   Loop  <- 100
   msg <- "NA"
   sub.size = length(y)
   while (loop <= Loop) {
      pr <- c(1 - 1 / (1 + exp(x %*% beta)))
      H <- t(x) %*% (pr * (1 - pr)   * x)
      S1 <- colSums((y - imp.y[1:sub.size])  * w * x[1:sub.size,]); 
      S2 <- colSums((imp.y - pr) * x); 
      S <- colSums(rbind(S1,S2))
      tryCatch(
         {shs <- NA
         shs <- solve(H + dig, S) }, 
         error=function(e){
            cat("\n ERROR :", loop, conditionMessage(e), "\n")})
      if (is.na(shs[1])) {
         msg <- "Not converge"
         beta <- loop <- NA
         break
      }
      beta.new <- beta + shs
      tlr  <- sum((beta.new - beta)^2)
      beta  <- beta.new
      if(tlr < 0.000001) {
         msg <- "Successful convergence"
         break
      }
      if (loop == Loop) {
         warning("Maximum iteration reached")
         ## beta <- NA
      }
      loop  <- loop + 1
   }
   list(par=beta, message=msg, iter=loop)
}



