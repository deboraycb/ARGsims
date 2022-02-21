## import packages
#source activate tsk034 # for conda env with tskit 0.3.4
import msprime
import tskit
import tsinfer
import tsdate
import time
import cyvcf2
import numpy as np
import argparse
import gzip

parser = argparse.ArgumentParser(description = "run tsinfer and tsdate with vcf from msprime simulation as input")
parser.add_argument('pref', help='prefix from simulation')
parser.add_argument('nspl', help='number of samples in simulation', type=int)
parser.add_argument('mu', help='mutation rate in simulation', type=float)
args = parser.parse_args()

## How I ran this: 
#python3 ./2_tsinfer_date.py std 8 2e-8 > tsinfer_date_std.out2
#python3 ./2_tsinfer_date.py r2e9 8 2e-8 > tsinfer_date_r2e9.out2
#python3 ./2_tsinfer_date.py m2e9 8 2e-9 > tsinfer_date_m2e9.out2
#python3 ./2_tsinfer_date.py r2e7 8 2e-8 > tsinfer_date_r2e7.out2
#python3 ./2_tsinfer_date.py m2e7 8 2e-7 > tsinfer_date_m2e7.out2
#python3 ./2_tsinfer_date.py n4 4 2e-8 > tsinfer_date_n4.out2
#python3 ./2_tsinfer_date.py n16 16 2e-8 > tsinfer_date_n16.out2
#python3 ./2_tsinfer_date.py n32 32 2e-8 > tsinfer_date_n32.out2 
#python3 ./2_tsinfer_date.py n80 80 2e-8 > tsinfer_date_n80.out2
#python3 ./2_tsinfer_date.py n200 200 2e-8 > tsinfer_date_n200.out2
#python3 ./2_tsinfer_date.py l5mb 8 2e-8 > tsinfer_date_l5mb.out2
#python3 ./2_tsinfer_date.py l250kb 8 2e-8 > tsinfer_date_l250kb.out2

simpref=args.pref
nspl=args.nspl
mu=args.mu

def add_diploid_sites(vcf, samples):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
    pos = 0
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF] + variant.ALT
        ancestral = variant.INFO.get("AA", variant.REF)
        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        allele_index = {
            old_index: ordered_alleles.index(allele)
            for old_index, allele in enumerate(alleles)
        }
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [
            allele_index[old_index]
            for row in variant.genotypes
            for old_index in row[0:2]
        ]
        samples.add_site(pos, genotypes=genotypes, alleles=alleles)

def chromosome_length(vcf):
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]

def print_tc(dated_tree, n, pref, twoNe=20000):
    for s1 in range(0, n-1):
        for s2 in range(s1+1, n):
            coaltimefh = gzip.open(simpref+'/sim_'+pref+'_spls'+str(s1)+'-'+str(s2)+'_inf.bed.gz', 'wt')
            for tree in dated_tree.trees(): 
                treeInterval = tree.get_interval() # length of sequence that this tree spans
                coalescence_time = tree.tmrca(s1,s2)
                print('1',round(treeInterval[0]), round(treeInterval[1]), coalescence_time/twoNe, sep = "\t", file = coaltimefh)
            coaltimefh.close()

# read vcf file
vcfname = '../vcf/sim_'+simpref+'_0.vcf'
vcf = cyvcf2.VCF(vcfname)
with tsinfer.SampleData(
    path="../vcf/sim_"+simpref+"_0.samples", sequence_length=chromosome_length(vcf)
) as samples:
    add_diploid_sites(vcf, samples)

print(
    "Sample file created for {} samples ".format(samples.num_samples)
    + "({} individuals) ".format(samples.num_individuals)
    + "with {} variable sites.".format(samples.num_sites),
    flush=True,
)
#Sample file created for 8 samples (8 individuals) with 207247 variable sites.

# Do the inference with tsinfer
print('tsinfer.infer runtime')
start = time.time()
inferred_ts = tsinfer.infer(samples)
print(time.time() - start)
print(
        "Inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        inferred_ts.num_trees, inferred_ts.sequence_length / 1e6, inferred_ts.num_edges
    )
)
# Simplify tsinfer trees
inferred_ts_simpl = inferred_ts.simplify(keep_unary=False)
print(
        "Simplified inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        inferred_ts_simpl.num_trees, inferred_ts_simpl.sequence_length / 1e6, inferred_ts_simpl.num_edges
    )
)

# tsdate with  DEFAULT PRIOR
print('tsdate.build_prior_grid timepoints = 20 runtime:')
start = time.time()
priors_def_simpl = tsdate.build_prior_grid(inferred_ts_simpl, timepoints = 20)
print(time.time() - start)
print('Default time grid (simpl) = ', priors_def_simpl.timepoints)

pref=simpref+'_priorDef'
print('tsdate.date (simpl) runtime')
start = time.time()
dated_ts_simpl = tsdate.date(inferred_ts_simpl, Ne=10000, mutation_rate=mu, priors=priors_def_simpl, ignore_oldest_root=True)
print(time.time() - start)
print_tc(dated_ts_simpl, nspl, pref)

### PRIOR WITH FINER GRID
print('tsdate.build_prior_grid timepoints = 100 runtime:')
start = time.time()
priors_moretimes = tsdate.build_prior_grid(inferred_ts_simpl, timepoints=100)
print(time.time() - start)
print('100tp time grid = ', priors_moretimes.timepoints)

pref=simpref+'_prior100tp'
dated_ts = tsdate.date(inferred_ts_simpl, Ne=10000, mutation_rate=mu, priors=priors_moretimes, ignore_oldest_root=True)
print_tc(dated_ts, nspl, pref)

### PRIOR WITH MANUALLY DEFINED GRID
print('tsdate.build_prior_grid np.geomspace(1e-5, 12, 50) runtime:')
start = time.time()
priors_max12 = tsdate.build_prior_grid(inferred_ts_simpl, timepoints=np.geomspace(1e-5, 12, 50))
print(time.time() - start)
print('GeomMax12 time grid = ', priors_max12.timepoints)

pref=simpref+'_priorGeomMax12'
dated_ts = tsdate.date(inferred_ts_simpl, Ne=10000, mutation_rate=mu, priors=priors_max12, ignore_oldest_root=True)
print_tc(dated_ts, nspl, pref)
