#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
 
option_list = list(
                   make_option(c("-n", "--nranks"), type="numeric", default=NULL,
                               help="number of levels in ranks (e.g. 101, 1001)", metavar="character"),
                   make_option(c("-i", "--quants"), type="character", default=NULL,
                               help="input: organized quantiles file", metavar="character"),
                   make_option(c("-s", "--nsamples"), type="numeric", default=8,
                               help="number of samples", metavar="numeric"),
                   make_option(c("-l", "--length"), type="numeric", default=100*10**6,
                               help="total length of sequence", metavar="numeric"),
                   make_option(c("-p", "--npairs"), type="numeric", default=28,
                               help="number of pairs of samples", metavar="numeric"),
                   make_option(c("-o", "--plot"), type="character", default="",
                               help="output: file name prefix for plot", metavar="character"))
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
infile<-opt$quants #read in the quantiles file
outpref <- opt$plot
if(outpref==""){
    tmp <- strsplit(infile, '\\.')[[1]]
    outpref <- paste0(tmp[1:(length(tmp)-1)], collapse=".")
}
toplot <- read.table(infile, col.names=c('Counts', 'mean_tcoal_msprime'))
nranks <- opt$nranks
nsamples <- opt$nsamples
npairs <- opt$npairs
seqlen <- opt$length

toplot$tcoal_rank <- 1:nrow(toplot)-1
mycolor='#DE4968FF' 

expect <- log10(npairs*seqlen/(nranks-1))
toround <- 0.5
tolim <- round(expect/toround)*toround

ggplot(toplot)+
    geom_point(aes(x=tcoal_rank, y=log10(Counts)))+
    xlab('Coalescence time rank')+
    ylab(expression('Number of sites (log'[10]*')'))+
    ylim(c(tolim-1,tolim+2))+
    geom_hline(yintercept=expect, color=mycolor)+
    theme_bw()
ggsave(paste0(outpref, '_compplot1_labels.png'), h=2, w=2)

ggplot(toplot)+
    geom_point(aes(x=tcoal_rank, y=log10(Counts)))+
    xlab('')+ylab('')+
    ylim(c(tolim-1,tolim+2))+
    geom_hline(yintercept=expect, color=mycolor)+
    theme_bw()
ggsave(paste0(outpref, '_compplot1.png'), h=2, w=2)

ggplot(toplot)+
    geom_point(aes(x=tcoal_rank, y=log10(Counts)))+
    xlab('')+ylab('')+
    geom_hline(yintercept=expect, color=mycolor)+
    theme_bw()
ggsave(paste0(outpref, '_compplot1_nolim.png'), h=2, w=2)
