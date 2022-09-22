N <- 100
q <- c(10000,50000,100000)
s <- c(seq(from=10, to=100, by=10),seq(from=150, to=500, by=50))
beta0 <- rep(-3,p)
setting <- c(1,2,3)
prob <- runif(p,0.05,0.5)
nctrl <- 12000*2
ncase <- 6000*2
methods <- c("standard","lasso","ebay","fiqt")
all.pcc <- matrix(0,nrow=N,ncol=(length(s)+1)*4)
colnames(all.pcc) <- c(paste0(s,"_",methods[1]),"all_standard",
                       paste0(s,"_",methods[2]),"all_lasso",
                       paste0(s,"_",methods[3]),"all_ebay",
                       paste0(s,"_",methods[4]),"all_fiqt")

all.signifpv <- matrix(0,nrow=N,ncol=(length(s)+1))
colnames(all.signifpv) <- c(1:((length(s)+1)))
rownames(all.signifpv) <- c(1:N)

calculate.pcc <- function(beta,truebeta) {
  sum(beta*truebeta)/sqrt(sum(beta^2))
}

25

lasso.or <- function(n00,n01,n10,n11,lambda) {
  or <- (n00*n11)/(n01*n10)
  if(or > 1){
    estbetahat = ((n00-lambda) * (n11-lambda)) / ((n10+lambda) * (n01+lambda))
    if(estbetahat < 1) estbetahat = 1
  } else {
    estbetahat = ((n00+lambda) * (n11+lambda)) / ((n10-lambda) * (n01-lambda))
    if(estbetahat > 1) estbetahat = 1
  }
  log(estbetahat)
}

D1ss <- function(x, y, xout = x, spar.offset = 0.1384, spl.spar = NULL)
{
  sp <- if (is.null(spl.spar))
  {
    sp <- smooth.spline(x, y)
    smooth.spline(x, y, spar = sp$spar + spar.offset)
  } else smooth.spline(x, y, spar = spl.spar)
  list(sp=sp,pred=predict(sp, xout, deriv = 1)$y)
}

bias.corr.kernel <- function(zz, xout, bw = "nrd0")
{
  density.obj = density(zz, bw = bw)
  fz.func = splinefun(density.obj$x, density.obj$y)
  
  26
  
  Psi.z = log(fz.func(zz)/dnorm(zz))
  D1ss.list <- D1ss(x = zz, y = Psi.z, xout = xout)
  list(sp=D1ss.list$sp,EZ=D1ss.list$pred)
}

Ebay <- function(OR, SE) {
  z <- log(OR) / SE
  ebay.z <- bias.corr.kernel(z, z)$EZ
  c(exp(ebay.z * SE))
}

p.adjust.BH <- function(p.vector, total.length=length(p.vector)) {
  m <- length(p.vector)
  ord <- order(p.vector)
  p.ord <- p.vector[ord]
  p.ord <- p.ord * total.length / seq(1,m)
  p.ord <- pmin(p.ord,rep(1,m))
  p.ord <- sapply(1:m,function(x)min(p.ord[x:m]))
  out.pv <- p.vector
  out.pv[ord] <- p.ord
  out.pv
}

FIQT.adjust <- function(z, total.length=length(z), min.p=10^-30) {
  pvals <- 2*pnorm(abs(z), low=FALSE)
  pvals[pvals < min.p] <- min.p
  27
  
  adj.pvals <- p.adjust.BH(pvals, total.length=total.length)
  mu.z <- sign(z) * qnorm(adj.pvals/2, low=FALSE)
  mu.z[abs(z) > qnorm(min.p/2, low=FALSE)] <- z[abs(z) > qnorm(min.p/2,low=FALSE)]
  mu.z
}

all.bias.std <- matrix(nrow=N,ncol=100)
all.bias.lasso <- matrix(nrow=N,ncol=100)
all.bias.ebay <- matrix(nrow=N,ncol=100)
all.bias.fiqt <- matrix(nrow=N,ncol=100)
for (u in 1:3) {
  for(v in 1:3){
    p <- q[u]
    setting <- v
    if(setting==1){
      beta1 <- c(runif(1286,log(1),log(1.01))*replace(a<-rbinom(1286,1,0.5),a==0,-1),
                 runif(1920,log(1.01),log(1.02))*replace(a<-rbinom(1920,1,0.5),a==0,-1),
                 runif(825,log(1.02),log(1.03))*replace(a<-rbinom(825,1,0.5),a==0,-1),
                 runif(156,log(1.03),log(1.04))*replace(a<-rbinom(156,1,0.5),a==0,-1),
                 runif(17,log(1.04),log(1.05))*replace(a<-rbinom(17,1,0.5),a==0,-1),
                 runif(9,log(1.05),log(1.075))*replace(a<-rbinom(9,1,0.5),a==0,-1),
                 runif(8,log(1.075),log(1.1))*replace(a<-rbinom(8,1,0.5),a==0,-1),
                 runif(12,log(1.1),log(1.15))*replace(a<-rbinom(12,1,0.5),a==0,-1),
                 runif(5,log(1.15),log(1.2))*replace(a<-rbinom(5,1,0.5),a==0,-1),
                 runif(2,log(1.2),log(1.5))*replace(a<-rbinom(2,1,0.5),a==0,-1),
                 rep(0,p-4240))
      
      28
      
    }else if(setting==2){
      beta1 <- c(runif(81,log(1),log(1.01))*replace(a<-rbinom(81,1,0.5),a==0,-1),
                 runif(205,log(1.01),log(1.02))*replace(a<-rbinom(205,1,0.5),a==0,-1),
                 runif(245,log(1.02),log(1.03))*replace(a<-rbinom(245,1,0.5),a==0,-1),
                 runif(203,log(1.03),log(1.04))*replace(a<-rbinom(203,1,0.5),a==0,-1),
                 runif(135,log(1.04),log(1.05))*replace(a<-rbinom(135,1,0.5),a==0,-1),
                 runif(114,log(1.05),log(1.075))*replace(a<-rbinom(114,1,0.5),a==0,-1),
                 runif(11,log(1.075),log(1.1))*replace(a<-rbinom(11,1,0.5),a==0,-1),
                 runif(2,log(1.1),log(1.15))*replace(a<-rbinom(2,1,0.5),a==0,-1),
                 runif(2,log(1.15),log(1.2))*replace(a<-rbinom(2,1,0.5),a==0,-1),
                 runif(5,log(1.2),log(1.5))*replace(a<-rbinom(5,1,0.5),a==0,-1),
                 rep(0,p-1003))
    }else if(setting==3){
      beta1 <- c(rep(log(1.5),10),rep(0,p-10))
    }
    for (j in 1:N) {
      all.beta <- rep(0,p)
      all.betase <- rep(0,p)
      all.lasso.beta <- matrix(0,p,length(s)+1)
      all.pv <- rep(0,p)
      
      number <- t(sapply(1:p,function(ii){
        q0 <- (1+exp(beta0[ii]+beta1[ii]))*(1-prob[ii])/
          (1+prob[ii]*exp(beta0[ii])+(1-prob[ii])*exp(beta0[ii]+beta1[ii]))
        n00 <- rbinom(1,nctrl,q0)
        n10 <- nctrl-n00
        
        29
        
        q1 <- (1+exp(beta0[ii]+beta1[ii]))*(1-prob[ii])/
          (prob[ii]*exp(beta1[ii])+1-prob[ii]+exp(beta0[ii]+beta1[ii]))
        n01 <- rbinom(1,ncase,q1)
        n11 <- ncase-n01
        c(n00,n01,n10,n11)
      }))
      colnames(number) <- c("n00", "n01","n10","n11")
      
      res <- sapply(1:p,function(ii){
        betahat <- log((number[ii,1] * number[ii,4]) /
                         (number[ii,2] * number[ii,3]))
        betahatse <- sqrt((1/number[ii,1]) + (1/number[ii,2])
                          + (1/number[ii,3]) + (1/number[ii,4]))
        pv <- pnorm(abs(betahat / betahatse),lower=FALSE)*2
        c(betahat,betahatse,pv)
      })
      all.beta <- res[1,]
      all.betase <- res[2,]
      all.pv <- res[3,]
      
      pcc = calculate.pcc(all.beta*(prob*(1-prob)),beta1*(prob*(1-prob)))
      all.pcc[j,(length(s)+1)] <- pcc
      
      kk <- 1
      for (k in c(s)) {
        betatilde <- all.beta
        
        30
        
        betatilde[order(all.pv,decreasing=TRUE)[1:(p-k)]] <- 0
        all.pcc[j,kk] <- calculate.pcc(betatilde*(prob*(1-prob)),
                                       beta1*(prob*(1-prob)))
        all.signifpv[j,kk] <- sort(all.pv,decreasing=FALSE)[k]
        kk <- kk + 1
      }
      all.signifpv[j,(length(s)+1)] <- sort(all.pv,decreasing = FALSE)[p]
      
      all.pcc[j,((length(s)+1)*2)] <- pcc
      max.lam <- sapply(1:p,function(ii)(abs(number[ii,1]*number[ii,4]
                                             -number[ii,2]*number[ii,3])/sum(number[ii,]))
                        /sqrt(prob[ii]*(1-prob[ii])))
      sorted.lam <- sort(max.lam,decreasing=TRUE)
      all.lasso.beta <- t(sapply(1:p,function(ii){
        sapply(s,function(k) {
          l <- sorted.lam[k+1]
          lasso.or(number[ii,1],number[ii,2],number[ii,3],number[ii,4],
                   
                   l*sqrt(prob[ii]*(1-prob[ii])))
          
        })
      }))
      
      for (kk in 1:length(s)) {
        all.pcc[j,kk+(length(s)+1)]<-calculate.pcc(all.lasso.beta[,kk]*(prob*(1-prob)),
                                                   beta1*(prob*(1-prob)))
      }
      
      31
      
      all.ebay.beta <- log(Ebay(exp(all.beta),all.betase))
      ebaybeta <- all.ebay.beta
      all.pcc[j,(length(s)+1)*3] <- calculate.pcc(ebaybeta*prob*(1-prob),
                                                  beta1*prob*(1-prob))
      kk <- (length(s)+1)*2+1
      for (k in c(s)) {
        ebaybeta <- all.ebay.beta
        ebaybeta[order(all.pv,decreasing=TRUE)[1:(p-k)]] <- 0
        all.pcc[j,kk] <- calculate.pcc(ebaybeta*prob*(1-prob),
                                       beta1*prob*(1-prob))
        kk <- kk + 1
      }
      
      all.fiqt.beta <- FIQT.adjust(all.beta/all.betase)*all.betase
      fiqtbeta <- all.fiqt.beta
      all.pcc[j,(length(s)+1)*4] <- calculate.pcc(fiqtbeta*prob*(1-prob),
                                                  beta1*prob*(1-prob))
      kk <- (length(s)+1)*3+1
      for (k in c(s)) {
        fiqtbeta <- all.fiqt.beta
        fiqtbeta[order(all.pv,decreasing=TRUE)[1:(p-k)]] <- 0
        all.pcc[j,kk] <- calculate.pcc(fiqtbeta*prob*(1-prob),beta1*prob*(1-prob))
        kk <- kk + 1
      }
      
      sel <- order(all.pv)[1:100]
      32
      
      all.bias.std[j,] <- all.beta[sel]*sign(beta1[sel]+1e-10)
      - beta1[sel]*sign(beta1[sel]+1e-10)
      all.bias.lasso[j,] <- all.lasso.beta[sel,11]*sign(beta1[sel]+1e-10)
      - beta1[sel]*sign(beta1[sel]+1e-10)
      all.bias.ebay[j,] <- all.ebay.beta[sel]*sign(beta1[sel]+1e-10)
      - beta1[sel]*sign(beta1[sel]+1e-10)
      all.bias.fiqt[j,] <- all.fiqt.beta[sel]*sign(beta1[sel]+1e-10)
      - beta1[sel]*sign(beta1[sel]+1e-10)
      
    }
    
    d1 <- data.frame(x = c(s),
                     
                     std = apply(all.pcc[,1:length(s)],2,mean),
                     lasso = apply(all.pcc[,(length(s)+2):(length(s)*2+1)],2,mean),
                     ebay = apply(all.pcc[,(length(s)*2+3):(length(s)*3+2)],2,mean),
                     fiqt = apply(all.pcc[,(length(s)*3+4):(length(s)*4+3)],2,mean),
                     p = apply(all.signifpv[,1:18],2,mean),
                     logp = -log10(apply(all.signifpv[,1:18],2,mean)))
    
    par(mar = c(5,5,3,5))
    par(xpd=NA,oma=c(3,0,0,0))
    with(d1, plot(x,std,type="o",pch=20,col="#F8766D", xaxt = "n",
                  main=paste("p=",q[u]),xlab="Number of SNPs",ylab="PCC",
                  cex.axis=1.4,cex.lab=1.2,cex.main=1.4,
                  ylim=range(std,lasso,ebay,fiqt)))
    
    axis(1,xlim=c(0,500),cex.axis=1.4)
    33
    
    with(d1,lines(x,lasso,type="o",pch=20,col="#7CAE00"))
    with(d1,lines(x,ebay,type="o",pch=20,col="#00BFC4"))
    with(d1,lines(x,fiqt,type="o",pch=20,col="#C77CFF"))
    par(new = T)
    with(d1, plot(x, logp, type="o", col="black",
                  
                  pch=4,axes=FALSE, xlab=NA, ylab=NA,ylim=range(logp)))
    
    axis(4,ylim=range(d1$logp),cex.axis=1.4,cex.lab=1.1)
    mtext(side = 4, line = 3, expression(-log[10](italic(p))),cex.lab=1.1)
    
    legend("right",
           legend=c("Standard","Lasso","Ebay","FIQT",expression(-log[10](italic(p)))),
           text.col=c("#F8766D","#7CAE00","#00BFC4","#C77CFF","black"),
           pch=c(20,20,20,20,4),
           col=c("#F8766D","#7CAE00","#00BFC4","#C77CFF","black"),
           cex=0.8,pt.cex=0.8,ncol=3)
    
    d2 <- data.frame(x = c(1:100),
                     
                     Standard = apply(all.bias.std,2,mean),
                     Lasso = apply(all.bias.lasso,2,mean),
                     Ebay = apply(all.bias.ebay,2,mean),
                     FIQT = apply(all.bias.fiqt,2,mean))
    
    library(ggplot2)
    library(reshape2)
    df <- melt(d2, id.vars="x")
    
    34
    
    ggplot(df, aes(x,value, col=variable)) +
      ggtitle(paste("p=",q[u])) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_point() +
      geom_line() +
      geom_hline(yintercept=0,linetype="dashed") +
      ylab("bias") +
      xlab("SNP")+
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18))
    
  }
}