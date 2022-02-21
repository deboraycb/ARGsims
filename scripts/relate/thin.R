suppressMessages(library(argparser))

parser <- arg_parser("Thin relate coalTimes")
parser <- add_argument(parser, 'infile', 'input file is the output of Aprils simulations_2/scripts/outputCoalTime: each line has [tree len] [trueTc] [sampledTc]*1000')
parser <- add_argument(parser, 'twoNe', 'Effective population size (haploid)')
argv <- parse_args(parser)

#infile <- 'coalTimeSample4Replicate1.txt'
infile <- argv$infile
pref <- sub('\\.txt$', '',(infile))
twoNe <- as.numeric(argv$twoNe)
data <- read.table(infile, header=F,
                   col.names=c('tree_len', 'trueTc', paste('sample', 1:1000, sep='_')))

thinned <- data[,c(1,2,seq(3,1002,10))]
write.table(thinned, paste0(pref, '_thinned.txt'),
            row.names=F, col.names=F)

outfile <- paste0(pref, '_thinned_ranks.txt')
datarank <- data.frame(rank=numeric(nrow(thinned)),
                       tree_len=thinned$tree_len,
                       trueTc=thinned$trueTc)
datarank$rank <- apply(thinned,1,function(x){
                       sum((x[2]/twoNe)<x[3:102])})
write.table(datarank,outfile,
            row.names=F, col.names=F)
