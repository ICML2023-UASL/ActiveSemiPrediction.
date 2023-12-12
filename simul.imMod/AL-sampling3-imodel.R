library(splines);library(ggplot2);library(xtable);
source('ActiveLearningFinal.R')
set.seed(20230810)
beta = c(-5.4, rep(0.7,7));re=300; t.n = 10000; m1 = 10;deg = 2; batch = 100; pilot <- 150
t.beta <- pilot.est <- matrix(0,re, length(beta))
mod.SES.CW  <- mod.SES.IMP <- 
   mod.SES.CW.one  <- mod.SES.IMP.one <-
   mod.SES.CW.all  <- mod.SES.IMP.all <-array(0, dim=c(m1,length(beta),re))
for(i in 1:re)
{
  dat = Datageneration(t.n, beta = beta, Dist = 'mzNormal', corr = 2*0.4, link = 'mod3')   
  X = dat$X; Y = dat$Y; mean(Y)
  t.beta[i,]  = glm(Y~X-1, family = 'binomial')$coefficient
  ind.pilot = sample(1:t.n, size = pilot)
  pilot.X=  X[ind.pilot,]; pilot.Y = Y[ind.pilot]
  pilot.beta = glm(pilot.Y ~ pilot.X-1, family = "binomial")$coefficient
  sub.p.u = 1/(t.n:(t.n-pilot+1))
  pool.X = X[-ind.pilot,]; pool.Y = Y[-ind.pilot]
  
  
  imp.moe = m.t(pilot.Y, pilot.X, pool.X, sub.p.u, deg, method = 'IMP')$i.model
  imp.moe.all = m.t(pilot.Y, pilot.X, pool.X, sub.p.u, deg, method = 'Perfect', link = 'mod3')$i.model
  imp.moe.one = m.t(pilot.Y, pilot.X, pool.X, sub.p.u, deg, method = 'Single')$i.model
  
  #inital
  set.ses.CW.one <-SES(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, i.model = imp.moe.one, setup = TRUE) 
  set.ses.IMP.one <-SES(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, i.model = imp.moe.one, setup = TRUE) 
  
  set.ses.CW.all <-SES(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, i.model = imp.moe.all, setup = TRUE) 
  set.ses.IMP.all <-SES(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, i.model = imp.moe.all, setup = TRUE) 
  
  set.ses.CW <-SES(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, i.model = imp.moe, setup = TRUE) 
  set.ses.IMP <-SES(pool.X, pool.Y, pilot.X, pilot.Y, sub.p.u, batch, pilot.beta, i.model = imp.moe, setup = TRUE) 
  
  for(j in 1:m1)
  { 
     ##Sampling
     set.ses.CW = SES(set.ses.CW$pool.x, set.ses.CW$pool.y, set.ses.CW$sub.x, set.ses.CW$sub.y, set.ses.CW$sub.p, batch, set.ses.CW$c.est, set.ses.CW$i.model)
     set.ses.IMP = SES(set.ses.IMP$pool.x, set.ses.IMP$pool.y, set.ses.IMP$sub.x, set.ses.IMP$sub.y, set.ses.IMP$sub.p, batch, set.ses.IMP$c.est, set.ses.IMP$i.model) 
     
     set.ses.CW.all = SES(set.ses.CW.all$pool.x, set.ses.CW.all$pool.y, set.ses.CW.all$sub.x, set.ses.CW.all$sub.y, set.ses.CW.all$sub.p, batch, set.ses.CW.all$c.est, set.ses.CW.all$i.model,link = 'mod3')
     set.ses.IMP.all = SES(set.ses.IMP.all$pool.x, set.ses.IMP.all$pool.y, set.ses.IMP.all$sub.x, set.ses.IMP.all$sub.y, set.ses.IMP.all$sub.p, batch, set.ses.IMP.all$c.est, set.ses.IMP.all$i.model,link = 'mod3') 
     
     set.ses.CW.one = SES(set.ses.CW.one$pool.x, set.ses.CW.one$pool.y, set.ses.CW.one$sub.x, set.ses.CW.one$sub.y, set.ses.CW.one$sub.p, batch, set.ses.CW.one$c.est, set.ses.CW.one$i.model,link = 'Single')
     set.ses.IMP.one = SES(set.ses.IMP.one$pool.x, set.ses.IMP.one$pool.y, set.ses.IMP.one$sub.x, set.ses.IMP.one$sub.y, set.ses.IMP.one$sub.p, batch, set.ses.IMP.one$c.est, set.ses.IMP.one$i.model,link = 'Single') 
     
     
     ##Model training
     mod.SES.CW[j,,i] <- set.ses.CW$c.est <- m.t(set.ses.CW$sub.y, set.ses.CW$sub.x, set.ses.CW$pool.x, set.ses.CW$pi, deg, method = 'IPW')$cur.est
     set.ses.CW$i.model <- m.t(set.ses.CW$sub.y, set.ses.CW$sub.x, set.ses.CW$pool.x, set.ses.CW$pi, deg, method = 'IMP')$i.model
     
     mod.SES.CW.all[j,,i] <- set.ses.CW.all$c.est <- m.t(set.ses.CW.all$sub.y, set.ses.CW.all$sub.x, set.ses.CW.all$pool.x, set.ses.CW.all$pi, deg, method = 'IPW')$cur.est
     set.ses.CW.all$i.model <- m.t(set.ses.CW.all$sub.y, set.ses.CW.all$sub.x, set.ses.CW.all$pool.x, set.ses.CW.all$pi, deg, method = 'Perfect', link = 'mod3')$i.model
     
     mod.SES.CW.one[j,,i] <- set.ses.CW.one$c.est <-m.t(set.ses.CW.one$sub.y, set.ses.CW.one$sub.x, set.ses.CW.one$pool.x, set.ses.CW.one$pi, deg, method = 'IPW')$cur.est
     set.ses.CW.one$i.model <- m.t(set.ses.CW.one$sub.y, set.ses.CW.one$sub.x, set.ses.CW.one$pool.x, set.ses.CW.one$pi, deg, method = 'Single')$i.model
     
     #SES mean
     fit.AI.IMP = m.t(set.ses.IMP$sub.y, set.ses.IMP$sub.x, set.ses.IMP$pool.x, set.ses.IMP$pi, deg, method = 'IMP')
     mod.SES.IMP[j,,i] <- set.ses.IMP$c.est <- fit.AI.IMP$cur.est
     set.ses.IMP$i.model <- fit.AI.IMP$i.model
     
     fit.AI.IMP.all = m.t(set.ses.IMP.all$sub.y, set.ses.IMP.all$sub.x, set.ses.IMP.all$pool.x, set.ses.IMP.all$pi, deg, method = 'Perfect', link = 'mod3')
     mod.SES.IMP.all[j,,i] <- set.ses.IMP.all$c.est <- fit.AI.IMP.all$cur.est
     set.ses.IMP.all$i.model <- fit.AI.IMP.all$i.model
     
     fit.AI.IMP.one = m.t(set.ses.IMP.one$sub.y, set.ses.IMP.one$sub.x, set.ses.IMP.one$pool.x, set.ses.IMP.one$pi, deg, method = 'Single')
     mod.SES.IMP.one[j,,i] <- set.ses.IMP.one$c.est <- fit.AI.IMP.one$cur.est
     set.ses.IMP.one$i.model <- fit.AI.IMP.one$i.model
     
     print(rbind(  t.beta[i,], set.ses.IMP.one$c.est, set.ses.IMP$c.est, set.ses.IMP.all$c.est))
     print(j)
  }
  print(i)
}
save.image("Case3-imodel.RData")

m.t.beta = colMeans(t.beta[1:re,])
s= 10
m.t.beta
rowMeans(mod.SES.CW.one[s,,1:re])
rowMeans(mod.SES.IMP.one[5,,1:re])
rowMeans(mod.SES.CW[s,,1:re])
rowMeans(mod.SES.IMP[s,,1:re])

tab1 <- tab2 <- tab3 <-tab4<- c()
for(s in 1:m1)
{  
  tab1 = cbind(tab1, round(sqrt(mean(colSums((mod.SES.CW.one[s,,1:re] -m.t.beta)^2))),4))
  tab2 = cbind(tab2, round(sqrt(mean(colSums((mod.SES.IMP.one[s,,1:re] -m.t.beta)^2))),4))
  tab3 = cbind(tab3, round(sqrt(mean(colSums((mod.SES.CW[s,,1:re] -m.t.beta)^2))),4))
  tab4 = cbind(tab4, round(sqrt(mean(colSums((mod.SES.IMP[s,,1:re] -m.t.beta)^2))),4))
}
tab1 = c(tab1);tab2=c(tab2);tab3=c(tab3);tab4 = c(tab4);


pdf(file = "~/Desktop/RPlot29.pdf", width = 8, height = 5)
par(mar = c(5.1, 6.1, 4.1, 1.1))
plot(tab3, ylim = c(min(tab4),max(tab1)+1) ,xaxt = "n",  cex.lab=2.5, lty = 6,
     col = 'blue',cex.axis=2,lwd= 4,#,main = expression(F[1])
     cex.main=2, xlab='Batch',ylab = '', type = 'b', pch = 4,  cex = 2)
title(ylab=expression(sqrt('MSE')),cex.lab=3)
points( tab1, col = 'deepskyblue', type = 'b', pch = 6, lwd= 4, cex = 2, lty=6)
points( tab4, col = 'red', type = 'b', pch = 8, lwd= 4,  lty=6, cex = 1.5)
points( tab2, col = 'pink', type = 'b', lty=6,  pch = 5, lwd=4, cex = 1.5)
axis(1, at=1:10, labels=c(1:10), cex.axis=2)
legend("topright", legend=c('CW-SBS-onlyX4','AI-SBS-onlyX4' , 'CW-SBS', 'AI-SBS'),
       col=c('deepskyblue','pink','blue','red'),cex = 2,
       pch = c(6,5,4,8), lty = c(6,6,6,6), lwd = 4, y.intersp = 1.2)
dev.off()



tab1 <- tab2 <- tab3 <-tab4 <- c()
for(s in 1:m1)
{  
  tab1 = cbind(tab1, sum((rowMeans(mod.SES.CW.one[s,,1:re]) - m.t.beta)^2))
  tab2 = cbind(tab2, sum((rowMeans(mod.SES.IMP.one[s,,1:re]) - m.t.beta)^2))
  tab3 = cbind(tab3, sum((rowMeans(mod.SES.CW[s,,1:re]) - m.t.beta)^2))
  tab4 = cbind(tab4, sum((rowMeans(mod.SES.IMP[s,,1:re]) - m.t.beta)^2))
}

tab1 = c(tab1);tab2=c(tab2);tab3=c(tab3);tab4 = c(tab4)

pdf(file = "~/Desktop/RPlot30.pdf", width = 8, height = 5)
par(mar = c(5.1, 6.1, 4.1, 1.1))
plot(tab3, ylim = c(-0.0001,max(tab1)) ,xaxt = "n",  cex.lab=2.5, lty = 6,
     col = 'blue',cex.axis=2,lwd= 4,#,main = expression(F[1])
     cex.main=2, xlab='Batch',ylab = '', type = 'b', pch = 4,  cex = 2)
title(ylab='Squared Bias',cex.lab=2.8)
points( tab1, col = 'deepskyblue', type = 'b', pch = 6, lwd= 4, cex = 2, lty=6)
points( tab4, col = 'red', type = 'b', pch = 8, lwd= 4,  lty=6, cex = 1.5)
points( tab2, col = 'pink', type = 'b', lty=6,  pch = 5, lwd=4, cex = 1.5)
abline(0,0, lwd  = 4, lty = 'dashed')
axis(1, at=1:10, labels=c(1:10), cex.axis=2)
legend("topright", legend=c('CW-SBS-onlyX4','AI-SBS-onlyX4' , 'CW-SBS', 'AI-SBS'),
       col=c('deepskyblue','pink','blue','red'),cex = 2,
       pch = c(6,5,4,8), lty = c(6,6,6,6), lwd = 4, y.intersp = 1.2)
dev.off()
