require(ggplot2)
library(plyr)

suppressMessages(library(optparse, quietly=TRUE))

args_list <- list(make_option(c('-p', '--pref'), type='character',
                         help='prefix of simulation (eg std, r2e9 etc)'),
                  make_option(c('-d', '--dir'), type='character', 
                              help='directory for bed files'),
                  make_option(c('-n', '--nspl'), type='integer',
                              help='number of samples (haplotypes)'),
                  make_option(c('-s', '--skip'), type='integer',
                              help='number of samples to skip when extracting pairs'))
args <- parse_args(OptionParser(option_list=args_list))

nspl <- args$nspl
skip <- args$skip
dir <- args$dir
pref <- args$pref

#nspl <- 8
#skip <- 1
#dir <- '../std/'
#pref <- paste0(dir,'sim_std_priorDef_spls')
#pref <- paste0(dir,'sim_std_prior100tp_spls')
#pref <- paste0(dir,'sim_std_priorGeomMax12_spls')

toplotall <- data.frame(Simulated=numeric(), 
                     PosteriorMean=numeric(),
                     x=numeric())
toplotallLOG <- data.frame(Simulated=numeric(), 
                     PosteriorMean=numeric(),
                     x=numeric())
for(s1 in seq(0,nspl-2,skip)){
    for(s2 in seq(s1+1,nspl-1,skip)){
    pair <- paste(s1,s2,sep='-')
    infile <- paste0(pref,pair,'_post.bed.gz')
    print(infile)
    df <- read.table(infile, col.names=c('chr','start','end','Simulated','PosteriorMean'))
    df$len <- df$end-df$start
    toplot <- aggregate(df$len,
                    list(Simulated=round(df$Simulated, 1), PosteriorMean=round(df$PosteriorMean, 1)),
                    FUN=sum)
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
xlab('')+
ylab('')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+11), guide=FALSE)+
geom_abline(slope=1,intercept=0)+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
theme_bw()
ggsave(paste0(pref,'_ts_pointest_lin_clean.png'), h=2, w=2)

ggplot(toplotallLOG)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('')+
ylab('')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+11), guide=FALSE)+
geom_abline(slope=1,intercept=0)+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
theme_bw()
ggsave(paste0(pref,'_ts_pointest_log_clean.png'), h=2, w=2)

ggplot(toplotall)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
#geom_smooth(aes(x=Simulated, y=PosteriorMean, weight=x), method=loess)+
xlab('Simulated coalTime\n(2Ne generations)')+
ylab('tsdate coalTime\n(2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+11))+
geom_abline(slope=1,intercept=0)+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
theme_bw()
#ggsave(paste0(pref,'_ts_pointest_lin_loess.png'), h=3, w=4)
ggsave(paste0(pref,'_ts_pointest_lin.png'), h=2, w=3)

ggplot(toplotallLOG)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('Simulated coalTime\n(log10 2Ne generations)')+
ylab('tsdate coalTime\n(log10 2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+11))+
geom_abline(slope=1,intercept=0)+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
theme_bw()
ggsave(paste0(pref,'_ts_pointest_log.png'), h=2, w=3)

toplot2 <- ddply(toplotall, .(Simulated),
                 function(x) {
                     out <- data.frame(len=sum(x$x),
                        Mean=weighted.mean(x$PosteriorMean, x$x),
                        SqErr=sum(((x$PosteriorMean-x$Simulated)^2)*x$x))
                     out$MSE <- out$SqErr/out$len
                     return(out)
                 })
names(toplot2)[1] <- 'SimRound'

toplotLOG2 <- ddply(toplotallLOG, .(Simulated),
                 function(x) {
                     out <- data.frame(len=sum(x$x),
                        Mean=weighted.mean(x$PosteriorMean, x$x),
                        SqErr=sum(((x$PosteriorMean-x$Simulated)^2)*x$x))
                     out$MSE <- out$SqErr/out$len
                     return(out)
                 })
names(toplotLOG2)[1] <- 'SimRound'

write.table(toplot2[,c('SimRound', 'len', 'Mean', 'MSE')], 
            file=paste0(pref,'meanPerSim.txt'),
            sep='\t', quote=F, row.names=F)
(mseall <- sum(toplot2$SqErr)/sum(toplot2$len))
cat(mseall,
    file=paste0(pref, 'MSEall.txt'))

ggplot(toplot2)+
geom_line(aes(x=SimRound, y=Mean))+
geom_line(aes(x=SimRound, y=MSE/20),
           lty='dashed')+
xlab('Simulated Tcoal\n(2Ne generations)')+
ylab('tsdate Tcoal\n(2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+11))+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8),
                   sec.axis = sec_axis(~.*20, name = "Mean square error"))+
geom_abline(slope=1,intercept=0, color='gray')+
theme_bw()
ggsave(paste0(pref,'meanMSE_lin.png'), h=2, w=3)

ggplot(toplotall)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
geom_line(data=toplot2, aes(x=SimRound, y=Mean), color='gray')+
xlab('Simulated coalTime\n(2Ne generations)')+
ylab('tsdate coalTime\n(2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+11))+
geom_abline(slope=1,intercept=0)+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
theme_bw()
ggsave(paste0(pref,'_ts_meanest_mean_lin.png'), h=2, w=3)

ggplot(toplotallLOG)+
    geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
    geom_line(data=toplotLOG2, aes(x=SimRound, y=Mean), color='gray')+
    xlab('Simulated Tcoal\n(log10 2Ne generations)')+
    ylab('ARGweaver Tcoal\n(log10 2Ne gen.)')+
    scale_fill_viridis_c(name='nSites', trans='log10', option='magma', 
                                              limits=c(1,1e+11))+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(pref,'meanest_mean_log.png'), h=2, w=3)

