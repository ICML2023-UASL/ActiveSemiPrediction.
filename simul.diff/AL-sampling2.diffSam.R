library(splines);library(ggplot2);library(xtable);
source('ActiveLearningFinal.R')
set.seed(20230801)
beta = c(-3, rep(0.7,7));re=300; t.n = 10000; m1 = 10;deg = 2; batch = 100; pilot <- 150
t.beta <- pilot.est <- matrix(0,re, length(beta))
mod.uni.IMP <- mod.CC.IMP   <- mod.Ent.IMP <- mod.SES.IMP <- array(0, dim=c(m1,length(beta),re))
time.ent <- time.SES  <- matrix(0, re, m1)
for(i in 1:re)
{
  dat = Datageneration(t.n, beta = beta, Dist = 'mzNormal', corr = 2*0.4, link = 'mod2')   
  X = dat$X; Y = dat$Y; mean(Y)
  t.beta[i,]  = glm(Y~X-1, family = 'binomial')$coefficient
  ind.pilot = sample(1:t.n, size = pilot)
  pilot.X=  X[ind.pilot,]; pilot.Y = Y[ind.pilot]
  pilot.beta = glm(pilot.Y ~ pilot.X-1, family = "binomial")$coefficient
  sub.p.u = 1/(t.n:(t.n-pilot+1))
  pool.X = X[-ind.pilot,]; pool.Y = Y[-ind.pilot]
  imp.moe = m.t(pilot.Y, pilot.X, pool.X, sub.p.u, deg, method = 'IMP')$i.model
  #inital
  set.uni.IMP <-Uni(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, setup = TRUE) 
  set.CC.IMP <-LC(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, setup = TRUE)
  set.ent.IMP <-entropy(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, setup = TRUE) 
  set.ses.IMP <-SES(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, i.model = imp.moe, setup = TRUE) 
  for(j in 1:m1)
  { 
    ##Sampling
    set.uni.IMP = Uni(set.uni.IMP$pool.x, set.uni.IMP$pool.y, set.uni.IMP$sub.x, set.uni.IMP$sub.y, set.uni.IMP$sub.p, batch, set.uni.IMP$c.est) 
    set.CC.IMP = LC(set.CC.IMP$pool.x, set.CC.IMP$pool.y, set.CC.IMP$sub.x, set.CC.IMP$sub.y, set.CC.IMP$sub.p, batch, set.CC.IMP$c.est) 
    set.ent.IMP = entropy(set.ent.IMP$pool.x, set.ent.IMP$pool.y, set.ent.IMP$sub.x, set.ent.IMP$sub.y, set.ent.IMP$sub.p, batch, set.ent.IMP$c.est) 
    set.ses.IMP = SES(set.ses.IMP$pool.x, set.ses.IMP$pool.y, set.ses.IMP$sub.x, set.ses.IMP$sub.y, set.ses.IMP$sub.p, batch, set.ses.IMP$c.est, set.ses.IMP$i.model) 
    
    ##Model training
    mod.uni.IMP[j,,i] <- set.uni.IMP$c.est <- m.t(set.uni.IMP$sub.y, set.uni.IMP$sub.x, set.uni.IMP$pool.x, set.uni.IMP$pi, deg, method = 'IMP')$cur.est
    mod.CC.IMP[j,,i] <- set.CC.IMP$c.est <- m.t(set.CC.IMP$sub.y, set.CC.IMP$sub.x, set.CC.IMP$pool.x, set.CC.IMP$pi, deg, method = 'IMP')$cur.est
    mod.Ent.IMP[j,,i] <- set.ent.IMP$c.est <- m.t(set.ent.IMP$sub.y, set.ent.IMP$sub.x, set.ent.IMP$pool.x, set.ent.IMP$pi, deg, method = 'IMP')$cur.est
    
    fit.AI.IMP = m.t(set.ses.IMP$sub.y, set.ses.IMP$sub.x, set.ses.IMP$pool.x, set.ses.IMP$pi, deg, method = 'IMP')
    mod.SES.IMP[j,,i] <- set.ses.IMP$c.est <- fit.AI.IMP$cur.est
    set.ses.IMP$i.model <- fit.AI.IMP$i.model
    
    print(rbind(  t.beta[i,], set.ent.IMP$c.est, set.ses.IMP$c.est))
    print(j)
  }
  print(i)
}

save.image("Case2-diffSam.RData")

m.t.beta = colMeans(t.beta[1:re,])
s= 10
m.t.beta
rowMeans(mod.uni.IMP[s,,1:re])
rowMeans(mod.CC.IMP[s,,1:re])
rowMeans(mod.Ent.IMP[s,,1:re])
rowMeans(mod.SES.IMP[s,,1:re])

tab1 <- tab2 <- tab3 <-tab4<- c()
for(s in 1:m1)
{  
  tab1 = cbind(tab1, round(sqrt(mean(colSums((mod.uni.IMP[s,,1:re] -m.t.beta)^2))),4))
  tab2 = cbind(tab2, round(sqrt(mean(colSums((mod.CC.IMP[s,,1:re] -m.t.beta)^2))),4))
  tab3 = cbind(tab3, round(sqrt(mean(colSums((mod.Ent.IMP[s,,1:re] -m.t.beta)^2))),4))
  tab4 = cbind(tab4, round(sqrt(mean(colSums((mod.SES.IMP[s,,1:re] -m.t.beta)^2))),4))
}
tab1 = c(tab1);tab2=c(tab2);tab3=c(tab3);tab4 = c(tab4);

pdf(file = "~/Desktop/RPlot09.pdf", width = 8, height = 5)
par(mar = c(5.1, 6.1, 4.1, 1.1))
plot(tab1, ylim = c(min(tab4),max(tab2)) ,xaxt = "n",  cex.lab=2.5, lty=3,
     col = 'green3',cex.axis=2,lwd= 4,#,main = expression(F[1])
     cex.main=2, xlab='Batch',ylab = '', type = 'b', pch = 16,  cex = 2)
title(ylab=expression(sqrt('MSE')),cex.lab=3)
points( tab2, col = 'cyan2', type = 'b', lty=4,  pch = 3, lwd=4, cex = 1.5)
points( tab3, col = 'orange1', type = 'b', pch = 17, lwd= 4, cex = 2, lty=2)
points( tab4, col = 'red', type = 'b', pch = 8, lwd= 4, cex = 1.5, lty=6)
axis(1, at=1:10, labels=c(1:10), cex.axis=2)
legend("topright", legend=c('AI-Uni','AI-LC' , 'AI-Ent', 'AI-SBS'),
       col=c('green3','cyan2','orange1','red'),cex = 2,
       pch = c(16,3,17,8), lty = c(3,4,2,6), lwd = 4, y.intersp = 1.2)
dev.off()



tab1 <- tab2 <- tab3 <-tab4 <- c()
for(s in 1:m1)
{  
  tab1 = cbind(tab1, sum((rowMeans(mod.uni.IMP[s,,1:re]) - m.t.beta)^2))
  tab2 = cbind(tab2, sum((rowMeans(mod.CC.IMP[s,,1:re]) - m.t.beta)^2))
  tab3 = cbind(tab3, sum((rowMeans(mod.Ent.IMP[s,,1:re]) - m.t.beta)^2))
  tab4 = cbind(tab4, sum((rowMeans(mod.SES.IMP[s,,1:re]) - m.t.beta)^2))
}

tab1 = c(tab1);tab2=c(tab2);tab3=c(tab3);tab4 = c(tab4)

pdf(file = "~/Desktop/RPlot10.pdf", width = 8, height = 5)
par(mar = c(5.1, 6.1, 4.1, 1.1))
plot(tab1, ylim = c(-0.001,max(tab4)) ,xaxt = "n",  cex.lab=2.5, lty=3,
     col = 'green3',cex.axis=2,lwd= 4,#,main = expression(F[1])
     cex.main=2, xlab='Batch',ylab = '', type = 'b', pch = 16,  cex = 2)
title(ylab='Squared Bias',cex.lab=3)
points( tab2, col = 'cyan2', type = 'b', lty=4,  pch = 3, lwd=4, cex = 1.5)
points( tab3, col = 'orange1', type = 'b', pch = 17, lwd= 4, cex = 2, lty=2)
points( tab4, col = 'red', type = 'b', pch = 8, lwd= 4, cex = 1.5, lty=6)
abline(0,0, lwd  = 4, lty = 'dashed')
axis(1, at=1:10, labels=c(1:10), cex.axis=2)
legend("topright", legend=c('AI-Uni','AI-LC' , 'AI-Ent', 'AI-SBS'),
       col=c('green3','cyan2','orange1','red'),cex = 2,
       pch = c(16,3,17,8), lty = c(3,4,2,6), lwd = 4, y.intersp = 1.2)
dev.off()
