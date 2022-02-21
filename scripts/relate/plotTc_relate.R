suppressMessages(library(ggplot2, quietly=TRUE))
suppressMessages(library(plotrix, quietly=TRUE))
suppressMessages(library(weights, quietly=TRUE))
suppressMessages(library(optparse, quietly=TRUE))
suppressMessages(library(LaplacesDemon, quietly=TRUE))

args_list <- list(make_option(c('-i', '--infile'), type='character',
                         help='bed file with coalescence times from argweaver'),
                  make_option(c('-N', '--Ne'), type='integer',
                              help='effective population size'))
args <- parse_args(OptionParser(option_list=args_list))

relatefile <- args$infile
Ne <- args$Ne

prefout <- sub('\\.txt$', '',relatefile)
cat(prefout)
relatetimes <- read.table(relatefile)
colnames(relatetimes) <- c('length', 'TrueTc', 1:(ncol(relatetimes)-2))

relatetimeslong <- data.frame(length=relatetimes$length,
                              TrueTc=relatetimes$TrueTc,
                              tmrca=rowMeans(relatetimes[,3:ncol(relatetimes)]))

# discretized theoretical
times <- c(0.000000,49.193481,122.586947,232.085215,395.449492,639.178343,1002.805899,1545.314509,2354.701987,3562.255340,5363.846221,8051.702368,12061.808515,18044.625462,26970.598323,40287.567936,60155.618452,89797.454603,134021.141756,200000.000000)/Ne
timebreaks <- numeric(length(times))
gettimes <- function(x, k=19 , delta=0.01, tmax=200000){(exp(x/k*log(1+delta*tmax))-1)/delta}
timebreaks <- gettimes(c(0,0:18+0.5,19))/Ne
f1 <- function(x, lambda=1){lambda*exp(-lambda*x)}
disc_exp1 <- numeric(length(timebreaks)-1)
for(ti in 2:length(timebreaks)-1){
    disc_exp1[ti] <- integrate(f1, timebreaks[ti], timebreaks[ti+1])[[1]]
}
#histogram/pmf
hall <- weighted.hist(relatetimeslong$tmrca, relatetimeslong$length, breaks=timebreaks, plot=F)
hall$counts <- hall$counts/sum(hall$counts)

mycolor='#DE4968FF' 
png(paste0(prefout, '_pmf_argwdisc_clean.png'), height=800, width=800, res=280)
par(mar=c(3,3,1,1), las=1, mgp=c(2,0.5,0))
plot(log(hall$mids), hall$counts, ylim=c(0,0.45), 
     xlog=T, tck=-0.03,
     main='', xlab='', ylab='')
points(log(hall$mids), hall$counts, , ylim=c(0,0.45),xlog=T, type='h')
legend(x=-8, y=.25, 'Relate\nsamples', col='black', pch=1, bty='n')
points(x=log(hall$mids), y=disc_exp1, col=mycolor, pch=8)
legend(x=-8, y=.4, 'Theoretical\nexpectation', col=mycolor, pch=8, bty='n')
dev.off()

png(paste0(prefout, '_pmf_argwdisc.png'), height=800, width=800, res=280)
par(mar=c(3,3,1,1), las=1, mgp=c(2,0.5,0))
plot(log(hall$mids), hall$counts, ylim=c(0,0.45), 
     xlog=T, tck=-0.03,
     xlab='log(Tcoal) (in 2Ne generations)',
     ylab='Frequency')
points(log(hall$mids), hall$counts, , ylim=c(0,0.45),xlog=T, type='h')
legend(x=-8, y=.25, 'Relate\nsamples', col='black', pch=1, bty='n')
points(x=log(hall$mids), y=disc_exp1, col=mycolor, pch=8)
legend(x=-8, y=.4, 'Theoretical\nexpectation', col=mycolor, pch=8, bty='n')
dev.off()

disc_exp1_mybreaks <- numeric(20)
mybreaks <- seq(0,max(relatetimeslong$tmrca), length.out=21)
f1 <- function(x, lambda=1){lambda*exp(-lambda*x)}
for(ti in 1:20){
    disc_exp1_mybreaks[ti] <- integrate(f1, mybreaks[ti], mybreaks[ti+1])[[1]]
}
##histogram/pdf
png(paste0(prefout, '_pmf.png'))
relhist <- wtd.hist(relatetimeslong$tmrca, weight=relatetimeslong$length, breaks=mybreaks)
relhist$counts <- relhist$counts/sum(relhist$counts)
plot(relhist,
              main=paste0('Theoretical vs. Relate PM'),
              xlab='Tcoal (in 2Ne generations)',
              ylim=c(0,0.35))
points(mybreaks[1:20]+(mybreaks[2:21]-mybreaks[1:20])/2, disc_exp1_mybreaks, type='h', col='red')
legend(x=3, y=0.22, 'Theoretical lambda=1 (discretized)', col='red', lty=1, bty='n')   
legend(x=3, y=0.18, 'Relate samples', col='black', pch=0, bty='n')
dev.off() 
