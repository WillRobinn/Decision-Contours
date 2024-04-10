
## illustrative dataset meta-analysis

RR <- c(0.157895,0.498016,0.324114,0.178804,0.333333) # relative risk
selnRR <- c(0.615225,0.701461,0.341625,0.403829,0.656167) # standard error for log relative risk
OR <- c(0.148361,0.491968,0.310397,0.164681,0.319149) # odds ratio
selnOR <- c(0.627499,0.71281,0.352263,0.415843,0.677454)
RD <- c(-0.0597,-0.012,-0.04148,-0.07765,-0.04167) # risk difference
seRD <- c(0.016944,0.011817,0.0118,0.015574,0.023421)
c <- c(61.8,20,41.8,0) # costs
u <- c(18.65652,20.58,18.65652,20.58) # utilities

install.packages("rmeta")
library("rmeta")

## all in one function for generating decision contours as shown in paper

decCont <- function(SS, seSS, method, p.usual, costs, utilities, wtp, outcome, sig.level=0.05, # data inputs: outcome is one of "lnRR", "lnOR", "RD". method is "fixed" or "random"
                    contour=TRUE, contour.points=200, # should contours be overlayed and how many NxN pixels for when numerical method is applied
                    uncertainty=FALSE, samples=1000, threshold=c(0.5,0.7,0.9),  greyscale=FALSE, # uncertainty or not. If greyscale is true then shades of 1-100 will be shown, if false then 50%, 70% and 90% thresholds are shown 
                    summ=FALSE, summ.pos=0, pred.interval=FALSE, plot.zero=FALSE, plot.summ=FALSE, ylim=NULL, xlim=NULL, legend=TRUE, # plotting options
                    expxticks=NULL, xticks=NULL, yticks=NULL, xlab=NULL, ylab=NULL, rand.load=10,
                    legendpos=c(xlim[2]+0.05*(xlim[2]-xlim[1]),ylim[2]), xpoints=NULL, ypoints=NULL, points=TRUE) {
  
  
  
  meta <- meta.summaries(SS, seSS, method=method, conf.level=(1-sig.level))
  tau2 <- meta$tau2
  ci <- qnorm(1-((sig.level)/2))
  
  length <- length(SS)
  
  c <- costs
  u <- utilities
  
  zero <- (p.usual*(wtp*(u[3]-u[4])-(c[3]-c[4])) + wtp*(u[4]-u[2])-c[4]+c[2])/(wtp*(u[1]-u[2])-(c[1]-c[2]))
  value1 = wtp*(u[1]-u[2])-(c[1]-c[2])
  value2 = p.usual*(wtp*(u[3]-u[4])-(c[3]-c[4]))+wtp*u[4]-c[4] - (wtp*u[2]-c[2])
  
  
  #####CURRENT WEIGHTINGS of studies
  
  if (method=="random") {
    df <- NROW(SS) - 1
    df2= df+1
    size <- 1/((seSS^2)+tau2)
  }
  
  else size <- 1/(seSS^2)
  
  
  sediff <- max(seSS) - min(seSS)
  
  if (!is.null(ylim)) if (ylim[1]<ylim[2]) ylim <- rev(ylim)
  
  if (is.null(ylim)) {
    ylim <- c(max(seSS) + 0.20*sediff, min(seSS) - 0.25*sediff)
    if (ylim[2]<0) ylim[2] <- 0
  }
  
  axisdiff <- ylim[2] - ylim[1]
  
  SSdiff <- max(SS) - min(SS)
  
  if (is.null(xlim)) {
    xlim <- c(min(SS) - 0.2*SSdiff, max(SS) + 0.2*SSdiff)
  }
  
  if (contour) {
    cSS <- seq(xlim[1], xlim[2], length.out=contour.points)
    csize <- seq(ylim[1], ylim[2], length.out=contour.points)
    csize[csize<=0] <- 0.0000001*min(seSS)
    for (k in 2:length(csize)) if (csize[k]==0 & csize[k-1]==0) csize[k] <- NA
    csize <- csize[!is.na(csize)]
  }
  
  
  if (contour) {
    
    if (method=="fixed" & uncertainty==FALSE) {
      
      vwt <- 1/(csize^2)
      
      if (outcome=="lnOR") {
        cSS <- (1/vwt)*((log(1/(p.usual*(1-zero))*zero*(1-p.usual)))*(sum(size) + vwt) - sum(size*SS))
      }
      
      if (outcome=="lnRR") {
        cSS <- (1/vwt)*(log((1/p.usual)*zero)*(sum(size) + vwt) - sum(size*SS))
      }
      
      if (outcome=="RD") {
        cSS <- (1/vwt)*((zero-p.usual)*(sum(size) + vwt) - sum(size*SS))
      }
      
    }
    
    if (method=="random" & uncertainty==FALSE)  {
      
      if (outcome=="RD") {
        
        matcont <- rep(0,times=length(csize))
        
        for (i in 1: length(csize))  {
          
          if (rand.load>0) {
            roundi<-i/rand.load
            flush.console()
            if (roundi==round(roundi,0)) {
              perc_complete <- (i/contour.points)*100
              cat(perc_complete, "%")
            }
            else cat(".")
          }
          
          if ( !is.na(csize[i]) ) {
            
            cSS1=min(cSS)
            
            while((meta.summaries(c(SS, cSS1), c(seSS, csize[i]), method=method)$summary + p.usual)*value1 - value2 > 0) {
              
              
              cSS1=cSS1+((max(cSS)-min(cSS))/contour.points)
              
            }
            
            matcont[i] <- cSS1
            
          }
          
        }
        
      }
      
      if (outcome=="lnRR") {
        
        matcont <- rep(0,times=length(csize)) 
        
        for (i in 1: length(csize))  {
          
          if (rand.load>0) {
            roundi<-i/rand.load
            flush.console()
            if (roundi==round(roundi,0)) {
              perc_complete <- (i/contour.points)*100
              cat(perc_complete, "%")
            }
            else cat(".")
          }
          
          if ( !is.na(csize[i]) ) {
            
            cSS1=min(cSS)
            
            while(p.usual*exp(meta.summaries(c(SS, cSS1), c(seSS, csize[i]), method=method)$summary)*value1 - value2 > 0) {
              
              
              cSS1=cSS1+((max(cSS)-min(cSS))/contour.points)
              
            }
            
            matcont[i] <- cSS1
            
          }
          
        }
        
      }
      
      if (outcome=="lnOR") {
        
        matcont <- rep(0,times=length(csize))
        
        for (i in 1: length(csize))  {
          
          if (rand.load>0) {
            roundi<-i/rand.load
            flush.console()
            if (roundi==round(roundi,0)) {
              perc_complete <- (i/contour.points)*100
              cat(perc_complete, "%")
            }
            else cat(".")
          }
          
          if ( !is.na(csize[i]) ) {
            
            cSS1=min(cSS)
            
            while(((exp(meta.summaries(c(SS, cSS1), c(seSS, csize[i]), method=method)$summary)*p.usual)/(1-p.usual+p.usual*exp(meta.summaries(c(SS, cSS1), c(seSS, csize[i]), method=method)$summary)))*value1 - value2 > 0) {
              
              
              cSS1=cSS1+((max(cSS)-min(cSS))/contour.points)
              
            }
            
            matcont[i] <- cSS1
            
          }
          
        }
        
      }
      
    }
    
    if (uncertainty==TRUE)  {
      
      matcont <- matrix(rep(NA,times=length(cSS)*length(csize)),nrow=length(csize))	
      overmax <- 0
      
      for (i in 1: length(csize))  {
        
        if (rand.load>0) {
          roundi<-i/rand.load
          flush.console()
          if (roundi==round(roundi,0)) {
            perc_complete <- (i/contour.points)*100
            cat(perc_complete, "%")
          }
          else cat(".")
        }
        
        if ( !is.na(csize[i]) ) {
          
          for (j in 1:length(cSS))  {
            
            metacont <- meta.summaries(c(SS, cSS[j]), c(seSS, csize[i]), method=method, conf.level=(1-sig.level))
            
            if (outcome=="RD") {
              s <- rnorm(samples,metacont$summary,metacont$se.summary)
              p.new = s + p.usual
            }
            
            if (outcome=="lnRR") {
              s <- rnorm(samples,metacont$summary,metacont$se.summary)
              p.new = exp(s)*p.usual
            }
            
            if (outcome=="lnOR") {
              s <- rnorm(samples,metacont$summary,metacont$se.summary)
              p.new = (exp(s)*p.usual)/(1-p.usual+exp(s)*p.usual)
            }
            
            if (greyscale==FALSE) {
              
              if ( sum(p.new*(wtp*(u[1]-u[2])-(c[1]-c[2]))+wtp*u[2]-c[2] > p.usual*(wtp*(u[3]-u[4])-(c[3]-c[4]))+wtp*u[4]-c[4]) < samples*threshold[1] ) matcont[i,j] <- 0
              if ( sum(p.new*(wtp*(u[1]-u[2])-(c[1]-c[2]))+wtp*u[2]-c[2] > p.usual*(wtp*(u[3]-u[4])-(c[3]-c[4]))+wtp*u[4]-c[4]) > samples*threshold[1]-1 ) matcont[i,j] <- 1
              if ( sum(p.new*(wtp*(u[1]-u[2])-(c[1]-c[2]))+wtp*u[2]-c[2] > p.usual*(wtp*(u[3]-u[4])-(c[3]-c[4]))+wtp*u[4]-c[4]) > samples*threshold[2]-1 ) matcont[i,j] <- 2
              if ( sum(p.new*(wtp*(u[1]-u[2])-(c[1]-c[2]))+wtp*u[2]-c[2] > p.usual*(wtp*(u[3]-u[4])-(c[3]-c[4]))+wtp*u[4]-c[4]) > samples*threshold[3]-1 ) matcont[i,j] <- 3
            }
            
            if (greyscale==TRUE) {
              
              x <- sum(p.new*(wtp*(u[1]-u[2])-(c[1]-c[2]))+wtp*u[2]-c[2] > p.usual*(wtp*(u[3]-u[4])-(c[3]-c[4]))+wtp*u[4]-c[4])
              
              x <- floor(200*(x/samples))
              
              if (x > 100) {
                x <- 200 - x                 
              }
              
              matcont[i,j] <- paste("gray",x,sep="")
              
            }
            
            
            
            
          }	   
          
        }
        
      }
      
    }
    
  }
  
  
  dev.new(width=11, height=11)
  par(mai = c(1, .8, .2, .2))
  par(plt=c(0.13, 0.95, 0.13, 0.95))
  
  
  if (is.null(expxticks) & is.null(xticks)) xaxis <- "s" else xaxis <- "n"
  
  if (is.null(yticks)) yaxis <- "s" else yaxis <- "n"
  
  if (!points) cexpoints <- 0
  else cexpoints <- 1
  
  plot(SS, seSS, ylim=ylim, xlim=xlim, xlab = xlab, ylab = ylab, pch=19, cex=cexpoints, cex.lab=1.35, cex.axis=1.35, col= "black", 
       xaxt=xaxis, yaxt=yaxis, xaxs="i", yaxs="i")
  
  if (contour & uncertainty==TRUE & greyscale==TRUE) {
    for (j in 1:(length(cSS)-1)) {
      for (i in 1:(length(csize)-1)) {
        
        polygon(c(cSS[j],cSS[j],cSS[j+1],cSS[j+1]),c(csize[i],csize[i+1],csize[i+1],csize[i]), 
                border=matcont[i,j], col = matcont[i,j])
        
      }
    }
  }
  
  if (contour & uncertainty==TRUE & greyscale==FALSE) {
    for (j in 1:(length(cSS)-1)) {
      for (i in 1:(length(csize)-1)) {
        
        if (matcont[i,j]==0)
          polygon(c(cSS[j],cSS[j],cSS[j+1],cSS[j+1]),c(csize[i],csize[i+1],csize[i+1],csize[i]), 
                  border="white", col = "white")
        if (matcont[i,j]==1)
          polygon(c(cSS[j],cSS[j],cSS[j+1],cSS[j+1]),c(csize[i],csize[i+1],csize[i+1],csize[i]), 
                  border="gray90", col = "gray90")
        if (matcont[i,j]==2)
          polygon(c(cSS[j],cSS[j],cSS[j+1],cSS[j+1]),c(csize[i],csize[i+1],csize[i+1],csize[i]), 
                  border="gray80", col = "gray80")
        if (matcont[i,j]==3)
          polygon(c(cSS[j],cSS[j],cSS[j+1],cSS[j+1]),c(csize[i],csize[i+1],csize[i+1],csize[i]), 
                  border="gray70", col = "gray70")
        
        
        
      }
    }
  }
  
  if (contour & uncertainty==FALSE & method=="random") {
    
    polygon(c(matcont,xlim[1],xlim[1]), c(csize,ylim[2],ylim[1]), border="gray72", col = "gray72")
  }
  
  if (contour & uncertainty==FALSE & method=="fixed") {
    
    polygon(c(cSS,xlim[1],xlim[1]), c(csize,ylim[2],ylim[1]), border="gray72", col = "gray72")
    
    
  }
  
  
  if (summ) {
    xsumm <- c(meta$summary - ci*meta$se.summary, meta$summary, meta$summary + ci*meta$se.summary, meta$summary)
    ysumm <- c(ylim[2]-0.10*axisdiff+summ.pos,ylim[2]-0.07*axisdiff+summ.pos,ylim[2]-0.10*axisdiff+summ.pos,ylim[2]-0.13*axisdiff+summ.pos)
    
    if (pred.interval) {	
      predint1 <- meta$summary - qt(p=sig.level,df=length-2)*(meta$tau2+(meta$se.summary^2))^0.5
      predint2 <- meta$summary + qt(p=sig.level,df=length-2)*(meta$tau2+(meta$se.summary^2))^0.5
    }
  }
  
  if (summ) {
    if (pred.interval & method=="random") {
      segments(x0=predint1, y0=ysumm[1], x1=predint2, y1=ysumm[1],col="black")
    }
    
    polygon(xsumm,ysumm, border="black", col = "lavenderblush4")
  }
  
  if (plot.summ) abline(v = meta$summary, col="slategrey")
  
  if (plot.zero) abline(v = 0, col="lightgrey", lty=1)
  
  box()
  
  if (greyscale==TRUE) {
    points(SS, seSS, ylim=ylim, xlim=xlim, xlab = xlab, ylab = ylab, pch=21, cex=cexpoints, col= "black",bg="white", xaxt=xaxis)
  }
  
  else {
    points(SS, seSS, ylim=ylim, xlim=xlim, xlab = xlab, ylab = ylab, pch=19, cex=cexpoints, col= "black", xaxt=xaxis)
  }
  
  if (outcome=="RD") {
    
    if ((meta.summaries(SS,seSS, method=method)$summary + p.usual)*value1 - value2 > 0) {
      title(main="Current Decision: NI is Cost-Effective",cex.main=1.5,line=0.6)
    }
    
    else {
      title(main="Current Decision: NI is not Cost-Effective",cex.main=1.5,line=0.6)
    }
    
    
  }
  
  if (outcome=="lnRR") {
    
    if (p.usual*exp(meta.summaries(SS,seSS, method=method)$summary)*value1 - value2 > 0) {
      title(main="Current Decision: NI is Cost-Effective",cex.main=1.5,line=0.6)
    }
    
    else {
      title(main="Current Decision: NI is not Cost-Effective",cex.main=1.5,line=0.6)
    }
    
    
  }
  
  if (outcome=="lnOR") {
    
    if (((exp(meta.summaries(SS, seSS, method=method)$summary)*p.usual)/(1-p.usual+p.usual*exp(meta.summaries(SS,seSS, method=method)$summary)))*value1 - value2 > 0) {
      title(main="Current Decision: NI is Cost-Effective",cex.main=1.5,line=0.6)
    }
    
    else {
      title(main="Current Decision: NI is not Cost-Effective",cex.main=1.5,line=0.6)
    }
    
    
  }
  
  if (uncertainty==FALSE) {
    legend("bottomleft",c(
      "Favour NI (Treatment A)","Favour Placebo (Treatment B)",
      "Pooled Effect","Null Effect"
    ),
    lty=c(
      0,0,
      1,1),pch=c(
        15,0,
        46,46),col=c(
          "gray72","black",
          "slategrey","lightgrey"),
    bg="white",cex=1)
  }
  
  if (uncertainty == TRUE & greyscale == FALSE) {
    legend("topright",c(paste(">",100*threshold[1],"% of samples favour placebo",sep=""),
                        paste(">",100*threshold[1],"% of samples favour NI",sep=""),
                        paste(">",100*threshold[2],"% of samples favour NI",sep=""),
                        paste(">",100*threshold[3],"% of samples favour NI",sep="")),
           bg="white",cex=1.15,
           fill=c("white","gray90","gray80","gray70"))
  }
  
  if (uncertainty & greyscale ==TRUE) {
    legend("bottomright",c("100%","90%","80%","70%","60%","50%"
    )
    ,bg="white",cex=1.15,
    fill=c("gray0","gray20","gray40","gray60","gray80","gray100"))
  } 
  
  
}


# An example call of this function would be:
decCont(SS=log(RR),seSS=selnRR,method="fixed",uncertainty=TRUE,greyscale=FALSE,
        summ=TRUE,pred.interval=TRUE,plot.zero=TRUE,plot.summ=TRUE,
        p.usual=0.05,costs=c,utilities=u,wtp=270,outcome="lnRR",
        xlab="Log Relative Risk",ylab="Standard Error",xlim=c(-2,2),
        ylim=c(0,1),contour=TRUE)


# function that overlays simulated studies onto the funnel plot (can be applied directly after contour generation)

sim <- function(SS,seSS,method,sig.level,samplesize,studies,p.usual,outcome) {
  
  meta <- meta.summaries(SS, seSS, method=method, conf.level=(1-sig.level))
  
  if (method=="fixed") {
    
    s <- rnorm(studies,meta$summary,meta$se.summary)
    
    if (outcome=="RD") {
      p.t = s + p.usual
    }
    
    if (outcome=="lnRR") {
      p.t = exp(s)*p.usual
    }
    
    if (outcome=="lnOR") {
      p.t = (exp(s)*p.usual)/(1-p.usual+exp(s)*p.usual)
    }
    
    p.c <- p.usual
    rc <- rep(NA,studies)
    rt <- rep(NA,studies)
    n <- samplesize
    
    
    for (i in 1:studies) {
      rc[i] <- rbinom(1,samplesize,p.c)
      rt[i] <- rbinom(1,samplesize,p.t[i])
    }
    
    
    if (outcome=="RD") {
      SS <- (rt/n)-(rc/n)
      seSS <- sqrt((rt*(n-rt)/(n^3))+rc*(n-rc)/(n^3))
    }
    
    if (outcome=="lnRR") {
      SS <- log(rt/rc)
      seSS <- sqrt((1/rt) + (1/rc) - (1/n) - (1/n))
    }
    
    if (outcome=="lnOR") {
      SS <- log((rt/n)/(1-(rt/n))*(1-(rc/n))/(rc/n))
      seSS <- sqrt((1/rt)+(1/(n-rt))+(1/rc)+(1/(n-rc)))
    }
    
  }
  
  if (method=="random") {
    
    t <- rt(studies,length(SS)-2)
    
    s <- meta$summary + t*sqrt(meta$tau2+meta$se.summary^2)
    
    if (outcome=="RD") {
      p.t = s + p.usual
    }
    
    if (outcome=="lnRR") {
      p.t = exp(s)*p.usual
    }
    
    if (outcome=="lnOR") {
      p.t = (exp(s)*p.usual)/(1-p.usual+exp(s)*p.usual)
    }
    
    p.c <- p.usual
    rc <- rep(NA,studies)
    rt <- rep(NA,studies)
    n <- samplesize
    p.t[p.t<0] <- 0
    
    for (i in 1:studies) {
      rc[i] <- rbinom(1,samplesize,p.c)
      rt[i] <- rbinom(1,samplesize,p.t[i])
    }
    
    
    if (outcome=="RD") {
      SS <- (rt/n)-(rc/n)
      seSS <- sqrt((rt*(n-rt)/(n^3))+rc*(n-rc)/(n^3))
    }
    
    if (outcome=="lnRR") {
      SS <- log(rt/rc)
      seSS <- sqrt((1/rt) + (1/rc) - (1/n) - (1/n))
    }
    
    if (outcome=="lnOR") {
      SS <- log((rt/n)/(1-(rt/n))*(1-(rc/n))/(rc/n))
      seSS <- sqrt((1/rt)+(1/(n-rt))+(1/rc)+(1/(n-rc)))
    }
    
  }
  
  points(jitter(SS,5),jitter(seSS,5),pch=3,col="Blue")
  
}


# An example call of this function (to be run after generating the contours) would be:
sim(RD,seRD,method="fixed",sig.level=0.05,samplesize=500,studies=200,p.usual=0.05,outcome="RD")

