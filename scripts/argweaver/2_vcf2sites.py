#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser(description = "convert VCF file from msprime to sites file for argweaver")
parser.add_argument('vcf', help='vcf file from msprime')
parser.add_argument('nspl', help='number of diploid samples')
parser.add_argument('seqlen', help='length of sequence in bp')
parser.add_argument('outdir', help='dir to output sites file to')
args = parser.parse_args()
vcffile = args.vcf
nspl = int(args.nspl)
seqlen = args.seqlen
outdir = args.outdir

# Example to run in cluster with slurm:
# sbatch -A fc_popgen -p savio2_htc -t 48:00:30 --wrap "./scripts/argw/2_vcf2sites.py vcf/sim_m2e9_0.vcf 4 100000000 argweaver/input/"
# sbatch -A fc_popgen -p savio2_htc -t 48:00:30 --wrap "./scripts/argw/2_vcf2sites.py vcf/sim_r2e7_0.vcf 4 100000000 argweaver/input/"
# sbatch -A fc_popgen -p savio2_htc -t 48:00:30 --wrap "./scripts/argw/2_vcf2sites.py vcf/sim_l5mb_0.vcf 4 5000000 argweaver/input/"
# sbatch -A fc_popgen -p savio2_htc -t 48:00:30 --wrap "./scripts/argw/2_vcf2sites.py vcf/sim_l250kb_0.vcf 4 250000 argweaver/input/"

infh = open(vcffile, 'r')
pref = vcffile.split('/')[-1].split('.vcf')[0]
outfile = pref+'.sites'
outfh = open(outdir+'/'+outfile, 'w')
outfh.write('\t'.join(['NAMES']+['spl'+str(i)+'_1'+'\t'+'spl'+str(i)+'_2' for i in range(nspl)])+'\n')
outfh.write('REGION\t1\t1\t'+seqlen+'\n')
allpos=[]
for line in infh:
    if line.startswith('#'):
        continue
    fields = line.split()
    pos = int(fields[1])
    while True:
        if pos in allpos:
            pos = pos+1
        else:
            break
    allpos.append(pos)
    alleles_01 = []
    [alleles_01.extend(x.split('|')) for x in fields[9:]]
    alleles_tmp = [x.replace('0', 'A') for x in alleles_01]
    #two lines below added for JC finite sites simulation
    alleles_tmp = [x.replace('2', 'G') for x in alleles_tmp]
    alleles_tmp = [x.replace('3', 'T') for x in alleles_tmp]
    alleles_nuc = [x.replace('1', 'C') for x in alleles_tmp]
    outfh.write(str(pos)+'\t'+''.join(alleles_nuc)+'\n')

infh.close()
outfh.close()
