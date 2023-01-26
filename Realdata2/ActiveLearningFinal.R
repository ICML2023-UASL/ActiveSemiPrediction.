library(mvtnorm)
Datageneration = function(n, beta, corr, Dist = 'mzNormal', link = 'mod1')
{
   d <- length(beta) - 2
   sigmax <- matrix(corr, d-1, d-1) + diag(2-corr, d-1)
   if( Dist == 'mzNormal' ){X  <- rmvnorm(n, rep(0, d-1), sigmax)}
   X =cbind(1, X,  runif(n, 0, 0.5), runif(n, 0, 0.5))
   
   #X = cbind(1, X)
   p <- d + 1
   if(link == 'mod1')
   {
      prob = 1 - 1 / (1 + exp(c(X %*% beta)))  
   }
   if(link == 'mod2')
   {
     prob = 1 - 1 / (1 + exp(c(X %*% beta + 0.5 * X[,2] * X[,3] + 0.5 * X[,5] * X[,7] - 0.5 * X[,4] * X[,8])))  
      #prob = 1 - 1 / (1 + exp(c(X %*% beta +  X[,4]^2 -  X[,7]^2 + X[,2] * X[,3] +  X[,5] * X[,7] - X[,4] * X[,8])))  
     
   }
   if(link == 'mod3')
   {
      prob = 1 - 1 / (1 + exp(c(X %*% beta) + 0.5* sin(0.5*X[,2]) - 0.5* sin(0.5*X[,3]) + 0.2* sin(0.2*X[,8])))  
   }
   if(link == 'mod4')
   {
     prob = 1 - 1 / (1 + exp(c(X %*% beta) + 0.5 * X[,6]^2 - 0.5 * X[,8]^2 + exp( X[,2:3]%*%rep(0.3,2) )))  
   }
   Y = rbinom(n,1,prob)
   list(Y = Y, X = X)
}
##############active entropy####
SES <- function(pool.x, pool.y, sub.x, sub.y, sub.p, batch, c.beta, i.model, setup = FALSE) 
{
   if(setup == TRUE)
   {
      pool.x = pool.x; pool.y = pool.y; sub.x = sub.x; sub.y=sub.y
      sub.p = sub.p; c.beta = c.beta; pi=NULL; i.model = i.model
   }else{
      full.x = rbind(sub.x, pool.x)
      #imp.v <- suppressMessages(predict(object=i.model, data.frame(pool.x), type='response'))
      imp.v = predict(i.model,newdata=data.frame(x2=pool.x[,2],x3=pool.x[,3],x4=pool.x[,4],x5=pool.x[,5],x6=pool.x[,6]),type="response")
      c.pool.n = length(pool.y); t.n = c.pool.n + length(sub.y)
      score = full.x %*% c.beta
      p_s = c(1 - 1 / (1 + exp(score)))
      w.s <- p_s * (1 - p_s)
      W.s <- solve(t(full.x) %*% (full.x * w.s))
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
m.t <- function(sub.y, sub.x, res.x, pi, deg, method = 'IPW', design){
   cur.n = length(sub.y); rest.n = dim(res.x)[1]; t.x = rbind(sub.x, res.x)
   if(method == 'IPW'){
      cur.est =  getest(sub.x, sub.y, w = pi)$par 
      iw = pi; reg =NULL
   }
   if(method == 'IMP'){
     x2 = sub.x[,2];x3 = sub.x[,3];x4 = sub.x[,4];x5 = sub.x[,5];x6 = sub.x[,6]
     reg = glm(sub.y~1+ns(x2,df=deg)+ns(x3,df=deg)+ns(x4,df=deg)+ns(x5,df=deg)+
                 ns(x6,df=deg), family=binomial)
     Imp.y = predict(reg,newdata=data.frame(x2=t.x[,2],x3=t.x[,3],x4=t.x[,4],x5=t.x[,5],x6=t.x[,6]),type="response")
     #data1 = data.frame(sub.y, sub.x)#,weights = pi
     #reg <- gbm(formula = sub.y ~ .,data=data1, distribution = "bernoulli")  
     #Imp.y <- suppressMessages(predict(object=reg, data.frame(t.x), type='response'))
     
     #cur.est =getest(t.x, Imp.y)$par 
     cur.est = getest.aug.now(t.x, sub.y, Imp.y, w = pi)$par
     iw = pi
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



