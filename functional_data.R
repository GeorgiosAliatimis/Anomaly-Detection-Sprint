data <- read.csv("Group1.csv")
x <- data$x

period = 1440 
dts = t(matrix(x,nrow=period))

plot_basic <- function(dt, median_index = 12,
                       outlier_indices = c() ,
                       other_indices = c(),
                       title="",ylabel="y"){
  plot(dt[median_index,],ylim=c( min(dt) , max(dt) ),
       main=title,ylab=ylabel,type="l", lwd=4)
  for (ind in outlier_indices){
    lines(dt[ind,],lty="dashed",col="red",lwd=3)
  } 
  for (ind in other_indices){
    lines(dt[ind,],lty="dotted")
  } 
}

plot_basic(dts,other_indices=c(1:50))


library(fdaoutlier)
source("my_functional_boxplot.R")



######## Seq_transform  ##########

# Transformations 
center_curves <- function(dt) dt - rowMeans(dt) 
normalize_curves <- function(dt) dt/sqrt(rowSums(dt^2))

### T0 
dts_T0 = dts 
out_T0 = my_functional_boxplot(dts_T0)

### T1 
dts_T1 = center_curves(dts_T0)
out_T1 = my_functional_boxplot(dts_T1)

### T2
dts_T2 = normalize_curves(dts_T1)
out_T2 = my_functional_boxplot(dts_T2)


plot_complete <- function(dts,out,title=""){
  plot_basic(dts, median_index = order(out$depth_values)[50],
             outlier_indices = out$outliers,
             other_indices = setdiff(c(1:50), out$outliers),
             title=title 
  )
  lines( out$lower, col="blue",lwd=2 )
  lines(out$upper, col="blue",lwd=2 )
  lines( out$inf, col="magenta",lwd=5)
  lines(out$sup, col="magenta",lwd=5 )
  time= c(1:period)
  polygon( c(time,rev(time)), c(out$sup, rev(out$inf)), col=rgb(1,0,1,0.5) )
}

plot_complete(dts_T0,out_T0,title="T0")
plot_complete(dts_T1,out_T1,title="T1")
plot_complete(dts_T2,out_T2,title="T2")

######## MUOD  ##########
N=19
plot( dts[N,] , colMeans(dts) , xlab=expression(Y[19]), ylab = expression(bar(Y)) )
abline( lm( colMeans(dts) ~ dts[N,] )$coefficients, col="red",lty="dashed", lwd=3)  

# Fast MUOD

fast_muod <- function(){
  N=50
  a <- rep(0,N)
  b <- rep(0,N)
  r <- rep(0,N)
  colmean = colMeans(dts)
  colmeanmean = mean(colmean)
  sy = sd(colmean)
  for(i in 1:N){
    sx = sd(dts[i,])
    sxy = cov(colmean,dts[i,])
    r[i] =  1 - sxy/(sx*sy)
    b[i] =  sxy/sx^2
    a[i] = abs ( colmeanmean - b[i] * mean(dts[i,]) )
    b[i] = abs(1 - b[i])
  }
  return(list( a=a,b=b,r=r ))
}
m = muod(dts)
m_fast = fast_muod()
a = m_fast$a
b = m_fast$b
r = m_fast$r

#### Muod histogram plots here
hist_plots <- function(vals1,vals2,xlab){
  p1 <- hist(vals1,breaks=9)                    
  p2 <- hist(vals2,breaks=9)                  
  plot( p1, col=rgb(0,0,1,1/4), xlab=xlab, main="" )  
  plot( p2, col=rgb(1,0,0,1/4), add=T)  
  legend("topright", legend=c('MUOD', 'Fast-MUOD'), 
         col=c('blue', 'pink'), lty=1, cex= .8)
}
# Shape
hist_plots(m$indices$shape,r,"Shape penalty")
outlier_cutoff_fast <- boxplot.stats( sort(r) )$stats[5]
outlier_cutoff <- boxplot.stats( sort(m$indices$shape) )$stats[5]
abline(v=outlier_cutoff,col=rgb(0,0,1,1), lty="dashed" )
abline(v=outlier_cutoff_fast,col=rgb(1,0,0,1),lty="dashed")
# Magnitude and amplitude
hist_plots(m$indices$magnitude,a,"Magnitude penalty")
hist_plots(m$indices$amplitude,b,"Amplitude penalty")

########## Timing ###########
N=100
time_seq_transform = system.time(for(i in 1:N) seq_transform(dts) )[3] / N 
time_boxplot = system.time(for(i in 1:N) functional_boxplot(dts) )[3] / N 
time_muod = system.time(for(i in 1:N) muod(dts) )[3] / N
time_fast = system.time(for(i in 1:N) fast_muod() )[3] / N
