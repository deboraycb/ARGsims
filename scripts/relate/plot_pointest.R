library(ggplot2)
require(dplyr)
suppressMessages(library(optparse, quietly=TRUE))

args_list <- list(make_option(c('-i', '--infile'), type='character',
                         help='input file *pairs_thinned.txt'))
args <- parse_args(OptionParser(option_list=args_list))


infile <- args$infile
cat(paste0(infile, '\n'))    
outpref <- strsplit(infile, '.txt')
thinned <- read.table(infile)
decplaces <- 1

# LINEAR SCALE
df <- data.frame(len=thinned[,1],
                 Simulated=thinned[,2]/20000,
                 PosteriorMean=rowMeans(thinned[, 3:ncol(thinned)]),
                 PosteriorMedian=apply(thinned[, 3:ncol(thinned)], 1, median))
df$WtdSqErr <- ((df$PosteriorMean-df$Simulated)^2)*df$len
df$WtdMean <- df$PosteriorMean*df$len
df$WtdMed <- df$PosteriorMedian*df$len
df$SimRound <- round(df$Simulated, decplaces)
head(df)
sum(df$WtdSqErr)/sum(df$len)

toplot <- aggregate(df$len, 
                    list(Simulated=df$SimRound, 
                         PosteriorMean=round(df$PosteriorMean, decplaces)), 
                    FUN=sum)
toplotmeans <- df %>% group_by(SimRound) %>% summarise_at(c('WtdMean','WtdSqErr','WtdMed','len'), sum)
toplotmeans$Mean <- toplotmeans$WtdMean/toplotmeans$len
toplotmeans$MSE <- toplotmeans$WtdSqErr/toplotmeans$len
toplotmeans$Median <- toplotmeans$WtdMed/toplotmeans$len

write.table(toplotmeans[,c('SimRound', 'len', 'Mean', 'Median', 'MSE')], 
            file=paste0(outpref,'_meanPerSim.txt'),
            sep='\t', quote=F, row.names=F)
cat(sum(toplotmeans$WtdSqErr)/sum(toplotmeans$len),
    file=paste0(outpref, '_MSEall.txt'))

#mean+MSE
ggplot(toplotmeans)+
geom_line(aes(x=SimRound, y=Mean))+
geom_line(aes(x=SimRound, y=MSE/20),lty='dashed')+
xlab('Simulated Tcoal\n(2Ne generations)')+
ylab('Relate Tcoal\n(2Ne generations)')+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8),
                   sec.axis = sec_axis(~.*20, name = "Mean square error"))+
geom_abline(slope=1,intercept=0, color='gray')+
theme_bw()
ggsave(paste0(outpref,'_meanMSE_lin.png'), h=2, w=3)

# WITHOUT LABELS
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('')+
ylab('')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10), guide=FALSE)+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref, '_meanest_lin_dec', decplaces,'_clean.png'), h=2, w=2)

#WITH LABELS
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('Simulated Tcoal\n(2Ne generations)')+
ylab('Relate Tcoal\n(2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref, '_meanest_lin_dec', decplaces,'.png'), h=2, w=3)

# WITHOUT LABELS
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
geom_line(data=toplotmeans, aes(x=SimRound, y=Mean), color='gray')+
xlab('')+
ylab('')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10), guide=FALSE)+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref, '_meanest_mean_lin_dec', decplaces,'_clean.png'), h=2, w=2)

#WITH LABELS
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
geom_line(data=toplotmeans, aes(x=SimRound, y=Mean), color='gray')+
xlab('Simulated Tcoal\n(2Ne generations)')+
ylab('Relate Tcoal\n(2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref, '_meanest_mean_lin_dec', decplaces,'.png'), h=2, w=3)

#WITH LABELS-median
toplot <- aggregate(df$len, 
                    list(Simulated=df$SimRound, 
                         PosteriorMedian=round(df$PosteriorMedian, decplaces)), 
                    FUN=sum)
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMedian, fill=x))+
xlab('Simulated Tcoal\n(2Ne generations)')+
ylab('Relate Tcoal\n(2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(0, 16))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref, '_medest_lin_dec', decplaces,'.png'), h=2, w=3)

# LOG SCALE
df <- data.frame(len=thinned[,1],
                 Simulated=log10(thinned[,2]/20000),
                 PosteriorMean=log10(rowMeans(thinned[, 3:ncol(thinned)])),
                 PosteriorMedian=log10(apply(thinned[, 3:ncol(thinned)], 1, median)))
df$WtdSqErr <- ((df$PosteriorMean-df$Simulated)^2)*df$len
df$PosteriorMean2<-rowMeans(log10(thinned[, 3:ncol(thinned)]))
df$WtdSqErr2 <- ((df$PosteriorMean2-df$Simulated)^2)*df$len
df$WtdMean <- df$PosteriorMean*df$len
df$WtdMed <- df$PosteriorMedian*df$len
df$SimRound <- round(df$Simulated, decplaces)
sum(df$WtdSqErr)/sum(df$len)
sum(df$WtdSqErr2)/sum(df$len)

toplot <- aggregate(df$len, 
                    list(Simulated=df$SimRound, 
                         PosteriorMean=round(df$PosteriorMean, decplaces)), 
                    FUN=sum)
sum(((toplot$Simulated-toplot$PosteriorMean)^2)*toplot$x)/sum(toplot$x)
toplotmeans <- df %>% group_by(SimRound) %>% summarise_at(c('WtdMean','WtdMed','len'), sum)
toplotmeans$Mean <- toplotmeans$WtdMean/toplotmeans$len
toplotmeans$Median <- toplotmeans$WtdMed/toplotmeans$len

#WITHOUT LABELS
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('')+
ylab('')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10), guide=FALSE)+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref,'_meanest_log_dec', decplaces, '_clean.png'), h=2, w=2)

#WITH LABELS
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
geom_line(data=toplotmeans, aes(x=SimRound, y=Mean), color='gray')+
xlab('Simulated Tcoal\n(log10 2Ne generations)')+
ylab('Relate Tcoal\n(log10 2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref,'_meanest_mean_log_dec', decplaces, '.png'), h=2, w=3)
#WITH LABELS
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMean, fill=x))+
xlab('Simulated Tcoal\n(log10 2Ne generations)')+
ylab('Relate Tcoal\n(log10 2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref,'_meanest_log_dec', decplaces, '.png'), h=2, w=3)
#WITH LABELS-median
toplot <- aggregate(df$len, list(Simulated=df$Simulated, PosteriorMedian=df$PosteriorMed), FUN=sum)
ggplot(toplot)+
geom_tile(aes(x=Simulated, y=PosteriorMedian, fill=x))+
xlab('Simulated Tcoal\n(log10 2Ne generations)')+
ylab('Relate Tcoal\n(log10 2Ne generations)')+
scale_fill_viridis_c(name='nSites', trans='log10', option='magma',
                     limits=c(1,1e+10))+
scale_x_continuous(expand=c(0,0), limits = c(-4, 1.5))+
scale_y_continuous(expand=c(0,0), limits = c(-4, 1.5))+
geom_abline(slope=1,intercept=0)+
theme_bw()
ggsave(paste0(outpref,'_medest_log_dec', decplaces, '.png'), h=2, w=3)
