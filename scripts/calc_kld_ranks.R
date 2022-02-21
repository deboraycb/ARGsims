suppressMessages(library(LaplacesDemon, quietly=TRUE))
suppressMessages(library(optparse, quietly=TRUE))

args_list <- list(make_option(c('-i', '--infile'), type='character',
                         help='bed file with coalescence times from argweaver'),
                  make_option(c('-r', '--nranks'), type='numeric',
                         help='number of ranks', default=101))
args <- parse_args(OptionParser(option_list=args_list))

suppressMessages(library(LaplacesDemon, quietly=TRUE))

file <- args$infile
nranks <- args$nranks

prefout <- paste0(sub('\\.txt$', '',file))
ranks <- read.table(file) 
nsites <- sum(ranks$V1)
cat('number of sites = ', nsites, '\n')
allranks <- rep(0:100, ranks$V1)
unif <- rep(1/nranks, nranks)
kld <- KLD(unif, ranks$V1/nsites)
cat('KLD = ', kld$sum.KLD.py.px, '\n')
round(kld$sum.KLD.py.px, 3)
write.table(kld$sum.KLD.py.px, paste0(prefout,'_KLD.txt'),
           quote=F, row.names=F, col.names=F) 
