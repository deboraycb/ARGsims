#!/bin/bash
#
# Partition:
#SBATCH --partition=savio
#savio has 164 nodes, 20cores/node, 64GB/node, 0.75 SU/h
#savio2 has 144 (+35) nodes, 24 (or 28) cores/node, 64GB/node, 1.00 SU/h
#savio_bigmem has 4 nodes, 20cores/node, 512GB/node, 1.67 SU/h
#savio2_bigmem has 28nodes, 24 cores/node, 128GB/node, 1.20 SU/h
#savio2_htc can run single core, 20 nodes, 12 cores/node, 128GB/node, 1.2 SU/h
#
# Wall clock limit:
#SBATCH --time=12:00:00
#
# Mail type:
#SBATCH --mail-type=end
#
# Mail user:
#SBATCH --mail-user=deboraycb@berkeley.edu
#
# Account:
#SBATCH --account=ac_popgen
#
# Specify Faculty Computing Allowance:
#SBATCH --qos=savio_normal
#
# to pass bash variables:
#PREF=l5mb_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/l5mb_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_l5mb_0_spls
#MCSPL=200-1200
#NSPL=8
#sbatch --job-name=aw${PREF}AllPairsPointest --output=%A_aw${PREF}AllPairsPointest.out --export=NSPL=$NSPL,PREF=$PREF,DIR=$DIR,TRUEPREF=$TRUEPREF,MCSPL=$MCSPL 5_pointestimates.sh
#
## Command(s) to run:
module load gnu-parallel

for s1 in $(seq 0 $((NSPL-2)))
do for s2 in $(seq $((s1+1)) $((NSPL-1)))
    do 
        echo "$s1 $s2"  
    done
done | parallel --jobs $SLURM_CPUS_ON_NODE ./getPostMean_bedmap.sh $PREF $DIR $TRUEPREF $MCSPL
wait
Rscript plot_pointestimates.R -p $PREF -d $DIR -n $NSPL -s $MCSPL
