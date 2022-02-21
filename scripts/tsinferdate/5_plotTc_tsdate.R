suppressMessages(library(scales, quietly=TRUE))
suppressMessages(library(tidyr, quietly=TRUE))
suppressMessages(library(ggplot2, quietly=TRUE))
suppressMessages(library(plotrix, quietly=TRUE))
suppressMessages(library(weights, quietly=TRUE))
suppressMessages(library(spatstat, quietly=TRUE))
suppressMessages(library(optparse, quietly=TRUE))
suppressMessages(library(LaplacesDemon, quietly=TRUE))

args_list <- list(make_option(c('-p', '--pref'), type='character',
                        help='pref of simulation'),
                  make_option(c('-N', '--Ne'), type='integer',
                        help='effective population size (2Ne)'),
                  make_option(c('-g', '--grid'), type='character', 
                        help='prefix for grid used in tsdate'))
args <- parse_args(OptionParser(option_list=args_list))
pref <- args$pref
Ne <- args$Ne
grid <- args$grid

# Rscript 3_plotTc_tsdate.R -p l5mb -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p n4 -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p r2e9 -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p m2e9 -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p std -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p std -N 20000 -g prior100tp
# Rscript 3_plotTc_tsdate.R -p std -N 20000 -g priorGeomMax12
# Rscript 3_plotTc_tsdate.R -p r2e7 -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p m2e7 -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p n16 -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p n32 -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p n80 -N 20000 -g priorDef
# Rscript 3_plotTc_tsdate.R -p n200 -N 20000 -g priorDef

# discretized theoretical
times <- c(0.000000,49.193481,122.586947,232.085215,395.449492,639.178343,1002.805899,1545.314509,2354.701987,3562.255340,5363.846221,8051.702368,12061.808515,18044.625462,26970.598323,40287.567936,60155.618452,89797.454603,134021.141756,200000.000000)/Ne
#timebreaks <- numeric(length(times))
gettimes <- function(x, k=19 , delta=0.01, tmax=200000){(exp(x/k*log(1+delta*tmax))-1)/delta}
timebreaks <- gettimes(c(0,0:18+0.5,19))/Ne
f1 <- function(x, lambda=1){lambda*exp(-lambda*x)}
disc_exp1 <- numeric(length(timebreaks)-1)
for(ti in 2:length(timebreaks)-1){
    disc_exp1[ti] <- integrate(f1, timebreaks[ti], timebreaks[ti+1])[[1]]
}
sum(disc_exp1)

tsdpath=paste0(pref, '/')
filelist=dir(path=tsdpath, pattern=paste0('sim_',pref,'_',grid,'_spls[0-9]+-[0-9]+_inf.bed'))
prefout <- paste0(tsdpath, paste0("sim_",pref,"_",grid,"_inf"))
tsdfile=filelist[1]
tsdtimes <- read.table(paste0(tsdpath,tsdfile))
hall <- weighted.hist(tsdtimes$V4, tsdtimes$V3-tsdtimes$V2, breaks=timebreaks, plot=F)
hall_grid <- weighted.hist(tsdtimes$V4, tsdtimes$V3-tsdtimes$V2, plot=F)

for(tsdfile in filelist[-1]){
    print(tsdfile)
    tsdtimes <- read.table(paste0(tsdpath,tsdfile))
    tmp <- weighted.hist(tsdtimes$V4, tsdtimes$V3-tsdtimes$V2, breaks=timebreaks, plot=F)
    tmp2 <- weighted.hist(tsdtimes$V4, tsdtimes$V3-tsdtimes$V2, plot=F)
    hall$counts <- hall$counts+tmp$counts
    hall_grid$counts <- hall_grid$counts+tmp2$counts
    #print(all(round(hall_grid$breaks, 1) == round(tmp2$breaks, 1)))
}

hall$counts <- hall$counts/sum(hall$counts)
hall_grid$counts <- hall_grid$counts/sum(hall_grid$counts)
class(hall_grid)='histogram'

mycolor='#DE4968FF'
png(paste0(prefout, '_pmf_clean.png'),  height=800, width=800, res=280)
par(mar=c(3,3,1,1), las=1, mgp=c(2,0.5,0))
plot(log(hall$mids), hall$counts,
     ylim=c(0,0.45),
     xlog=T, tck=-0.03,
     main='', xlab='', ylab='')
points(log(hall$mids), hall$counts,ylim=c(0,0.45),xlog=T, type='h')
legend(x=-8, y=.25, 'tsdate\nestimates', col='black', pch=1, bty='n')
points(x=log(hall$mids), y=disc_exp1, col=mycolor, pch=8)
legend(x=-8, y=.4, 'Theoretical\nexpectation', col=mycolor, pch=8, bty='n')
dev.off()

png(paste0(prefout, '_pmf.png'),  height=800, width=800, res=280)
par(mar=c(3,3,1,1), las=1, mgp=c(2,0.5,0))
plot(log(hall$mids), hall$counts,
     ylim=c(0,0.45),
     xlog=T, tck=-0.03,
     xlab='log(Tcoal) (in 2Ne generations)',
     ylab='Frequency')
points(log(hall$mids), hall$counts,ylim=c(0,0.45),xlog=T, type='h')
legend(x=-8, y=.25, 'tsdate\nestimates', col='black', pch=1, bty='n')
points(x=log(hall$mids), y=disc_exp1, col=mycolor, pch=8)
legend(x=-8, y=.4, 'Theoretical\nexpectation', col=mycolor, pch=8, bty='n')
dev.off()

png(paste0(prefout, '_simplehist.png')) 
weighted.hist(tsdtimes$V4, w=tsdtimes$V3-tsdtimes$V2)
dev.off()

mybreaks <- seq(0,max(tsdtimes$V4), length.out=21)
disc_exp1_mybreaks <- numeric(20)
f1 <- function(x, lambda=1){lambda*exp(-lambda*x)}
for(ti in 1:20){
    disc_exp1_mybreaks[ti] <- integrate(f1, mybreaks[ti], mybreaks[ti+1])[[1]]
}

png(paste0(prefout, '_pmf_nodisc.png')) 
tsdhist <- wtd.hist(tsdtimes$V4, weight=tsdtimes$V3-tsdtimes$V2, breaks=mybreaks)
tsdhist$counts <- tsdhist$counts/sum(tsdhist$counts)
plot(tsdhist,
      main=paste0('Theoretical vs. tsdate PM'),
      xlab='Tcoal (in 2Ne generations)')
points(mybreaks[1:20]+(mybreaks[2:21]-mybreaks[1:20])/2, disc_exp1_mybreaks, type='h', col='red')
legend(x=0.25, y=0.2, 'Theoretical lambda=1 (discretized)', col='red', lty=1, bty='n')   
legend(x=0.25, y=0.16, 'tsdate', col='black', pch=0, bty='n')
dev.off() 

gridbreaks <- hall_grid$breaks
ngridbr <- length(gridbreaks)
disc_exp1_gridbreaks <- numeric(ngridbr-1)
f1 <- function(x, lambda=1){lambda*exp(-lambda*x)}
for(ti in 1:(ngridbr-1)){
    disc_exp1_gridbreaks[ti] <- integrate(f1, gridbreaks[ti], gridbreaks[ti+1])[[1]]
}

png(paste0(prefout, '_pmf_grid.png')) 
plot(hall_grid,
      main=paste0('Theoretical vs. tsdate PM'),
      xlab='Tcoal (in 2Ne generations)')
points(gridbreaks[1:(ngridbr-1)]+(gridbreaks[2:ngridbr]-gridbreaks[1:(ngridbr-1)])/2, disc_exp1_gridbreaks, type='h', col='red')
legend(x=0.25, y=0.2, 'Theoretical lambda=1 (discretized)', col='red', lty=1, bty='n')   
legend(x=0.25, y=0.16, 'tsdate', col='black', pch=0, bty='n')
dev.off() 
