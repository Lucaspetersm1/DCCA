#DCCA stands for Detrended cross-correlational analysis. It calculates a correlation coefficient 
#between two time series. In this case, speciation rate with temperature or sea level.
#DCCA does this using a 'detrended analysis'. This is because we need to remove long term trends
#that could distort the data and affect short term trends, which is what we are focusing on.
DCCA <- function(x,y,s){
  xx<-cumsum(x)
  yy<-cumsum(y)
  t<-1:length(xx)
  F2sj_xy<-runif(floor(length(xx)/s))
  F2sj_xx<-F2sj_xy
  F2sj_yy<-F2sj_xy
  for(ss in seq(1,(floor(length(xx)/s)*s),by=s)){
    F2sj_xy[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_xx[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_yy[(ss-1)/s+1]<-sum((summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
  }
  rho<-mean(F2sj_xy)/sqrt(mean(F2sj_xx)*mean(F2sj_yy))
  return(c(rho,1/sqrt(length(xx)),1-pnorm(abs(rho),mean=0,sd=1/sqrt(length(xx)))))
}