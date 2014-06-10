qq.plot <- function(tab,lambda=F,stat="CHISQ",BA=F,plot=T,
                    pch.col="red",max.axis=10,scale.cex=T,dens=T,...) {

  ## plot obs:exp pvalues of assoc tests
  ## optionally, correct for lambda
  ## if plot=F, append to existing plot
  ## pch.col can be either a single color spec, or a vector of
  ## length nrow(tab) for differential display
  ## max.axis is the upper limit for xlim AND ylim
  ## scale.cex controls whether to scale plot character size with -log10(p) or not.

  ## Chris Cotsapas 2007.

  
  chi.ix <- match(stat, colnames(tab))

  if(lambda) {
    l.fac <- median(na.omit(tab[,chi.ix]))/0.456
    pvals <- -log10(1-pchisq(na.omit(tab[,chi.ix])/l.fac,df=1))
  }
  else {
    pvals <- -log10(1-pchisq(na.omit(tab[,chi.ix]),df=1))
  }

  ## If providing differential colors, sort color vector too!
  if(length(pch.col)>1) {
    s.obj <- sort(pvals, index.return=T)
    pvals <- pvals[s.obj$ix]
    pch.col <- pch.col[s.obj$ix]
  }
  else {
    pvals <- sort(pvals)
  }
  
  p.exp <- sort(-log10( c(1:length(pvals))/(length(pvals)+1) ))

  ## figure out scaling factor here
  if(scale.cex) {
    my.cex <- 0.25 * pvals
  }
  else{
    my.cex <- 1.5
  }
   


  
  if(BA) {
    bland.altman(pvals,p.exp,log=F,plot=plot,pch=23,
                 bg=pch.col,cex=my.cex,...)
  }

  else {
    if(plot) {
      par(mar=c(5,4,4,3)+0.1)
      plot(p.exp,pvals, xlab="Expected (-logP)", ylab="Observed (-logP)",
           xlim=c(0,max.axis), ylim=c(0,max.axis), type="n", xaxs="i", yaxs="i",bty="l")
      if(dens){    
        ## rescaled density plots in bg
        d.vals <- density(pvals)
        lines(d.vals$x, (max.axis/2)+((max.axis/2)*(d.vals$y/max(d.vals$y))),
              col="lightgray",lty="dotted",lwd=2)
        ## add density axis
        axis(side=4, at=c(max.axis/2,max.axis), labels=c(0,round(max(d.vals$y),1)))
        mtext("Data density",side=4,at=max.axis*0.75,line=1,adj=0.5)
      }
      lines(c(0,max.axis),c(0,max.axis), col="lightgray", lwd=2,lty=2)
      points(p.exp, pvals, pch=23, cex=my.cex, bg=pch.col)
    }
    else {
      points(p.exp,pvals, pch=23,cex=my.cex, bg=pch.col)
    }

  }
}
#Manhattan plot script for GWAMA
#Written by Joshua C Randall & Reedik Magi
for (e in commandArgs(trailingOnly=TRUE))
{
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2]))
  {
    assign(ta[[1]][1],ta[[1]][2])
  } else {
    assign(ta[[1]][1],TRUE)
  }
}
if(!exists("input"))
{
  input <- paste("gwama.out")
}
if(!exists("out")) {
  out <- paste(input,".manh.png",sep="")
}
data<-read.table(input,stringsAsFactors=FALSE,header=TRUE,sep = "\t",na.strings = "-9")
png(out,height=600,width=800)

obspval <- (data$p.value)
chr <- (data$chromosome)
pos <- (data$position)
obsmax <- trunc(max(-log10(obspval)))+1

sort.ind <- order(chr, pos) 
chr <- chr[sort.ind]
pos <- pos[sort.ind]
obspval <- obspval[sort.ind]

x <- 1:22
x2<- 1:22

for (i in 1:22)
{
	 curchr=which(chr==i)
	 x[i] <- trunc((max(pos[curchr]))/100) +100000
	 x2[i] <- trunc((min(pos[curchr]))/100) -100000
}

x[1]=x[1]-x2[1]
x2[1]=0-x2[1]

for (i in 2:24)
{
	x[i] <- x[i-1]-x2[i]+x[i]
	x2[i] <- x[i-1]-x2[i]

}
locX = trunc(pos/100) + x2[chr]
locY = -log10(obspval)
col1=rgb(0,0,108,maxColorValue=255)
col2=rgb(100,149,237,maxColorValue=255)
col3=rgb(0,205,102,maxColorValue=255)
col4 <- ifelse (chr%%2==0, col1, col2)
curcol <- ifelse (obspval<5e-8, col3, col4) 
plot(locX,locY,pch=20,col=curcol,axes=F,ylab="-log10 p-value",xlab="",bty="n",ylim=c(0,obsmax),cex=0.8)
axis(2,las=1)
for (i in 1:22)
{
	labpos = (x[i] + x2[i]) / 2
	mtext(i,1,at=labpos,cex=0.8,line=0)
}
mtext("Chromosome",1,at=x[22]/2,cex=1,line=1)
dev.off()
