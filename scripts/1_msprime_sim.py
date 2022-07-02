## import packages
# use source activate msp0 to use msprime version 0.7.4
import msprime

vcfdir='../vcf/'
tcdir='../tcoalmsp/'
tsdir='../trees/'
nrep=1

def simulate(nrep, pref, n, mu=2e-8, rec=2e-8, Ne=10000, length=100*10**6, vcfdir=vcfdir, tcdir=tcdir, tsdir=tsdir, model='hudson'):
    for i in range(nrep):
        print('rep', i)
        tree_sequence = msprime.simulate(sample_size=n, Ne=Ne, recombination_rate=rec, mutation_rate=mu, length=length, model=model)
        vcffh = open(vcfdir+'sim_'+pref+str(i)+'.vcf', 'w')
        tree_sequence.write_vcf(vcffh, ploidy=2, position_transform='legacy', individual_names=['spl'+str(s) for s in range(int(n/2))])
        vcffh.close()
        tsfile = tsdir+'sim_'+pref+str(i)+'.trees'
        tree_sequence.dump(tsfile)
        for s1 in range(0, n-1):
            for s2 in range(s1+1, n):
                coaltimefh = open(tcdir+'rep'+str(i)+'/sim_'+pref+str(i)+'_spls'+str(s1)+'-'+str(s2)+'.tc', 'w')
                for tree in tree_sequence.trees(): 
                    treeInterval = tree.get_interval() # length of sequence that this tree spans
                    coalescence_time = tree.tmrca(s1,s2)
                    print(treeInterval[0], treeInterval[1], coalescence_time, sep = "\t", file = coaltimefh)
                coaltimefh.close()

## standard simulation
# coalescent model
# 4 diploids
# Ne 10000
# recombination 2e-8
# mutation 2e-8
# length 100mb
# replicates 10
simulate(nrep, pref='std_', n=8, rec=2e-8, Ne=10000)

## change number of samples
# coalescent model
# 40 diploids           *
# Ne 10000
# recombination 2e-8
# mutation 2e-8
# length 100mb
# replicates 10
simulate(20, pref='n200_', n=200, rec=2e-8, Ne=10000, length=5*10**6)
simulate(nrep, pref='n80_', n=80, rec=2e-8, Ne=10000)
simulate(nrep, pref='n32_', n=32, rec=2e-8, Ne=10000)
simulate(nrep, pref='n16_', n=16, rec=2e-8, Ne=10000)
simulate(nrep, pref='n4_', n=4, rec=2e-8, Ne=10000)

## /10 recombination
# coalescent model
# 4 diploids           
# Ne 10000
# recombination 2e-9    *
# mutation 2e-8
# length 100mb
# replicates 10
simulate(nrep, pref='r2e9_', n=8, rec=2e-9, Ne=10000)
## *10 recombination
simulate(nrep, pref='r2e7_', n=8, rec=2e-7, Ne=10000)

## *10 mutation
# coalescent model
# 4 diploids           
# Ne 10000
# recombination 2e-8    
# mutation 2e-7         *
# length 100mb
# replicates 10
simulate(nrep, pref='m2e7_', n=8, mu=2e-7, Ne=10000)
## /10 mutation
simulate(nrep, pref='m2e9_', n=8, mu=2e-9, Ne=10000)
    
# mut/rec ratio = 2,4
simulate(1, pref='m4e8_', n=8, mu=4e-8, Ne=10000)
simulate(1, pref='m8e8_', n=8, mu=8e-8, Ne=10000)
    
## standard simulation under smcprime
# coalescent model
# 4 diploids
# Ne 10000
# recombination 2e-8
# mutation 2e-8
# length 100mb
# replicates 10
simulate(nrep, pref='stdsmc_', n=8, rec=2e-8, Ne=10000, model='smc')
simulate(nrep, pref='stdsmcp_', n=8, rec=2e-8, Ne=10000, model='smc_prime')
## /10 rec simulation under smcprime
simulate(nrep, pref='recsmc_', n=8, rec=2e-9, Ne=10000, model='smc')
simulate(nrep, pref='recsmcp_', n=8, rec=2e-9, Ne=10000, model='smc_prime')

## shorter input sequences
# coalescent model
# 4 diploids
# Ne 10000
# recombination 2e-8
# mutation 2e-8
# length 100mb
simulate(nrep=1, pref='l5mb_', n=8, rec=2e-8, Ne=10000, length=5*10**6)
simulate(nrep=1, pref='l250kb_', n=8, rec=2e-8, Ne=10000, length=250*10**3)
