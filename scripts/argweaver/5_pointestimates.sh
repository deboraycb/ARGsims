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
#PREF=r2e9_jc690
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/r2e9hudjc69_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_r2e9_0_spls
#MCSPL=200-1200

#PREF=m4e8_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/m4e8hud_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_m4e8_0_spls
#MCSPL=1200-2200

#PREF=m8e8_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/m8e8hud_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_m8e8_0_spls
#MCSPL=1200-2200

#PREF=m4e8_jc690
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/m4e8hudjc69_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_m4e8_0_spls
#MCSPL=1200-2200

#PREF=m8e8_jc690
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/m8e8hudjc69_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_m8e8_0_spls
#MCSPL=1200-2200

#PREF=std_jc690
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/stdhudjc69_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_std_0_spls
#MCSPL=200-1200

#PREF=m2e7_jc690
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/m2e7hudjc69_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_m2e7_0_spls
#MCSPL=1200-2200

#PREF=m2e9_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/m2e9hud_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_m2e9_0_spls
#MCSPL=200-1200

#PREF=r2e7_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/r2e7hud_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_r2e7_0_spls
#MCSPL=200-1200

#PREF=l250kb_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/l250kb_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_l250kb_0_spls
#PREF=l5mb_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/l5mb_smcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_l5mb_0_spls
#MCSPL=200-1200
#
#PREF=stdsmc_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/stdsmc/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_stdsmc_0_spls
#MCSPL=200-1200
#
#PREF=stdsmcp_0
#DIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/stdsmcp/rep0
#TRUEPREF=/global/home/users/deboraycb/scratch/argbias/simulations_3/tcoalmsp/rep0/sim_stdsmcp_0_spls
#MCSPL=200-1200
#
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
Rscript ~/scratch/argbias/simulations_3/scripts/argw/plot_pointestimates.R -p $PREF -d $DIR -n $NSPL -s $MCSPL
