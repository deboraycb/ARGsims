require(ggplot2)
require(dplyr)
#infile <- "/global/scratch/deboraycb/argbias/simulations_3/scripts/arg-sample_sim_std_0_0-1_Tcpair_posterior.csv"
suppressMessages(library(optparse, quietly=TRUE))

args_list <- list(make_option(c('-p', '--pref'), type='character',
                         help='prefix of simulation (eg std, r2e9 etc)'),
                  make_option(c('-d', '--dir'), type='character', 
                              help='directory for bed files'),
                  make_option(c('-n', '--nspl'), type='integer',
                              help='number of samples (haplotypes)'),
                  make_option(c('-s', '--mcspl'), type='integer',
                              help='MCMC samples taken (interval separated by - )'))
args <- parse_args(OptionParser(option_list=args_list))

nspl <- args$nspl
pref <- args$pref
dir <- args$dir
mcspl <- args$mcspl
outpref <- paste0(dir,'/arg-sample_sim_',pref,'_')

# 1) POINT ESTIMATES FROM MEAN OF MCMC SAMPLES

toplotall <- data.frame(Simulated=numeric(), 
                     PosteriorMean=numeric(),
                     x=numeric())
toplotallLOG <- data.frame(Simulated=numeric(), 
                     PosteriorMean=numeric(),
                     x=numeric())
for(s1 in 0:(nspl-2)){
    for(s2 in (s1+1):(nspl-1)){
    pair <- paste(s1,s2,sep='-')
    infile <- paste0(outpref,pair,'_Tcpair_spl',mcspl,'_posterior.bed')
    print(infile)
    df <- read.table(infile,
                     col.names=c('chr','start','end','Simulated','PosteriorMean','PosteriorMedian'))
    df$len <- df$end-df$start
    # for this pair of s1,s2,
    # aggregate total length for each combination of simulated vs. mean of estimated values, both rounded to .1
    toplot <- aggregate(df$len,
                    list(Simulated=round(df$Simulated, 1), PosteriorMean=round(df$PosteriorMean, 1)),
                    FUN=sum)
    # aggregate across pairs
    tmp <- rbind(toplot, toplotall)
    toplotall <- aggregate(tmp$x,
                           list(Simulated=tmp$Simulated, PosteriorMean=tmp$PosteriorMean),
                           FUN=sum)
    toplotLOG <- aggregate(df$len,
                     list(Simulated=round(log10(df$Simulated), 1),
                                  PosteriorMean=round(log10(df$PosteriorMean), 1)),
                     FUN=sum)
    tmpLOG <- rbind(toplotLOG,toplotallLOG)
    toplotallLOG <- aggregate(tmpLOG$x,
                           list(Simulated=tmpLOG$Simulated, PosteriorMean=tmpLOG$PosteriorMean),
                           FUN=sum)
}}

dim(toplotall)
dim(toplotallLOG)
ggplot(toplotall)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('Simulated Tcoal\n(2Ne generations)')+
ylab('ARGweaver Tcoal\n(2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref,'meanest_lin.png'), h=2, w=3)

ggplot(toplotallLOG)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('Simulated Tcoal\n(log10 2Ne generations)')+
ylab('ARGweaver Tcoal\n(log10 2Ne gen.)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma', 
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref,'meanest_log.png'), h=2, w=3)

ggplot(toplotall)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('')+ylab('')+ 
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10), guide=F)+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref,'meanest_lin_clean.png'), h=2, w=2)

ggplot(toplotallLOG)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('')+ylab('')+ 
scale_fill_viridis_c(name='nSites', trans='log10', option='magma', 
                     limits=c(1,1e+10), guide=F)+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref,'meanest_log_clean.png'), h=2, w=2)

# take means of all estimated values across all sites and pairs, for each bin of true value
meansall <- data.frame(SimRound=numeric(), 
                    WtdMean=numeric(),
                    WtdMed=numeric(),
                    len=numeric())
meansallLOG <- data.frame(SimlogRound=numeric(), 
                    WtdlogMean=numeric(),
                    WtdlogMed=numeric(),
                    len=numeric())
for(s1 in 0:(nspl-2)){
    for(s2 in (s1+1):(nspl-1)){
    pair <- paste(s1,s2,sep='-')
    infile <- paste0(outpref,pair,'_Tcpair_spl',mcspl,'_posterior.bed')
    print(infile)
    df <- read.table(infile,
                     col.names=c('chr','start','end','Simulated','PosteriorMean','PosteriorMedian'))
    df$len <- df$end-df$start
    df$SimRound <- round(df$Simulated, 1)
    df$WtdMean <- df$PosteriorMean*df$len
    df$WtdSqErr <- ((df$PosteriorMean-df$Simulated)^2)*df$len
    df$WtdMed <- df$PosteriorMedian*df$len
    means <- df %>% group_by(SimRound) %>% summarise_at(c('WtdMean','WtdSqErr','WtdMed','len'), sum)
    tmp <- rbind(means, meansall)
    meansall <- tmp %>% group_by(SimRound) %>%summarise_at(c('WtdMean','WtdSqErr','WtdMed','len'), sum)
}}

meansall$Mean <- meansall$WtdMean/meansall$len
meansall$MSE <- meansall$WtdSqErr/meansall$len
meansall$Median <- meansall$WtdMed/meansall$len
cat(sum(meansall$WtdSqErr)/sum(meansall$len),
    file=paste0(outpref, 'MSEall.txt'))

ggplot(toplotall)+
geom_line(data=meansall, aes(x=SimRound, y=Mean))+
geom_line(data=meansall, aes(x=SimRound, y=MSE/20),
           lty='dashed')+
xlab('Simulated Tcoal\n(2Ne generations)')+
ylab('ARGweaver Tcoal\n(2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8),
                   sec.axis = sec_axis(~.*20, name = "Mean square error"))+
geom_abline(slope=1,intercept=0, color='gray')+
theme_bw()
ggsave(paste0(outpref,'MeanMSE_lin.png'), h=2, w=3)
