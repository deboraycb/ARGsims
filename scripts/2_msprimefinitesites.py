#use source activate msp1 for msprime version1
import msprime
import tskit

vcfdir='../vcf/'
tcdir='../tcoalmsp/'
tsdir='../treeseq/'

def addmutations(outpref, ts, n, mutmodel='jc69', mu=2e-8, vcfdir=vcfdir, tsdir=tsdir):
    ts_mutated = msprime.sim_mutations(ts, rate=mu, model=mutmodel, keep=False)
    vcffh = open(vcfdir+'sim_'+outpref+str(i)+'.vcf', 'w')
    ts_mutated.write_vcf(vcffh, ploidy=2, position_transform='legacy', individual_names=['spl'+str(s) for s in range(int(n/2))])
    vcffh.close()
    tsfile = tsdir+'sim_'+outpref+str(i)+'.trees'
    ts_mutated.dump(tsfile)

i=0
n=8
mutmodel='jc69'

pref='std_'
pref='r2e9_'

tsfile = tsdir+'sim_'+pref+str(i)+'.trees'
ts_v0 = tskit.load(tsfile)
outpref=pref+mutmodel
addmutations(outpref, ts_v0, n, vcfdir=vcfdir, tsdir=tsdir)

pref='m2e7_'
tsfile = tsdir+'sim_'+pref+str(i)+'.trees'
ts_v0 = tskit.load(tsfile)
outpref=pref+mutmodel
addmutations(outpref, ts_v0, n, mu=2e-7, vcfdir=vcfdir, tsdir=tsdir)

pref='m4e8_'
tsfile = tsdir+'sim_'+pref+str(i)+'.trees'
ts_v0 = tskit.load(tsfile)
outpref=pref+mutmodel
addmutations(outpref, ts_v0, n, mu=4e-8, vcfdir=vcfdir, tsdir=tsdir)

pref='m8e8_'
tsfile = tsdir+'sim_'+pref+str(i)+'.trees'
ts_v0 = tskit.load(tsfile)
outpref=pref+mutmodel
addmutations(outpref, ts_v0, n, mu=8e-8, vcfdir=vcfdir, tsdir=tsdir)
