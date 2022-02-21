#!/bin/bash
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# Mail type:
#SBATCH --mail-type=end
#
# Mail user:
#SBATCH --mail-user=deboraycb@berkeley.edu
#
# Account:
#SBATCH --account=fc_popgen
#
# Specify Faculty Computing Allowance:
#SBATCH --qos=savio_normal
#
# to pass bash variables:
#
# cd /global/scratch/users/deboraycb/argbias/simulations_3/argweaver/
# REP=0; OUTPREF=m2e9hud_smcp; SITES=input/sim_m2e9_${REP}.sites; MU=2e-9; RE=2e-8; Ne=1e4; NITERS=1200; OUTDIR=${OUTPREF}/rep${REP}; LEN=100000000
# REP=0; OUTPREF=r2e7hud_smcp; SITES=input/sim_r2e7_${REP}.sites; MU=2e-8; RE=2e-7; Ne=1e4; NITERS=1200; OUTDIR=${OUTPREF}/rep${REP}; LEN=100000000
# REP=0; OUTPREF=l5mb_smcp; SITES=input/sim_l5mb_${REP}.sites; MU=2e-8; RE=2e-8; Ne=1e4; NITERS=1200; OUTDIR=${OUTPREF}/rep${REP}; LEN=5000000
# REP=0; OUTPREF=l250kb_smcp; SITES=input/sim_l250kb_${REP}.sites; MU=2e-8; RE=2e-8; Ne=1e4; NITERS=1200; OUTDIR=${OUTPREF}/rep${REP}; LEN=250000
#
# mkdir -p $OUTDIR; mkdir -p ${OUTDIR}_reps; mkdir -p ${OUTDIR}_long
# sbatch --partition=savio2_bigmem --nodes=1 --cpus-per-task=2 --job-name=AW${OUTPREF} --output=%A_AW${OUTPREF}.out --export=SITES=$SITES,MU=$MU,RE=$RE,Ne=$Ne,OUTDIR=$OUTDIR,NITERS=$NITERS,LEN=$LEN ../scripts/argw/3_run_argweaver_smcprime_hudsim.sh
#
## Command(s) to run:
module load gnu-parallel

export WDIR=/global/scratch/users/deboraycb/argbias/simulations_3/argweaver/
cd $WDIR

# run ARGweaver in 20 chunks of 5Mb
PREF=$(basename $SITES .sites)
COMMANDS=${OUTDIR}/commands_arg-sample_${SLURM_JOB_ID}.lst

S=1
CHUNKLEN=$((LEN/20))
E=$CHUNKLEN
OUTFILE=${OUTDIR}/arg-sample_${PREF}_${S}-${E}
if [ -f "${OUTFILE}.stats" ]; then
    echo "arg-sample --sites $SITES --region 1:$S-$E --mutrate $MU --recombrate $RE --popsize $Ne --output $OUTFILE --smc-prime --iters $NITERS --resume > ${OUTFILE}_${SLURM_JOB_ID}.out 2> ${OUTFILE}_${SLURM_JOB_ID}.err" > $COMMANDS
else
    echo "arg-sample --sites $SITES --region 1:$S-$E --mutrate $MU --recombrate $RE --popsize $Ne --output $OUTFILE --smc-prime --iters $NITERS > ${OUTFILE}_${SLURM_JOB_ID}.out 2> ${OUTFILE}_${SLURM_JOB_ID}.err" > $COMMANDS
fi

for i in {2..20}
do 
    S=$E
    E=$((i*CHUNKLEN))
    OUTFILE=${OUTDIR}/arg-sample_${PREF}_${S}-${E}
    if [ -f "${OUTFILE}.stats" ]; then
        echo "arg-sample --sites $SITES --region 1:$S-$E --mutrate $MU --recombrate $RE --popsize $Ne --output $OUTFILE --smc-prime --iters $NITERS --resume > ${OUTFILE}_${SLURM_JOB_ID}.out 2> ${OUTFILE}_${SLURM_JOB_ID}.err" >>  $COMMANDS
    else
        echo "arg-sample --sites $SITES --region 1:$S-$E --mutrate $MU --recombrate $RE --popsize $Ne --output $OUTFILE --smc-prime --iters $NITERS > ${OUTFILE}_${SLURM_JOB_ID}.out 2> ${OUTFILE}_${SLURM_JOB_ID}.err" >>  $COMMANDS
    fi
done

# set number of jobs based on number of cores available and number of threads per job
export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $SLURM_CPUS_PER_TASK ))

parallel --delay 1 -j $JOBS_PER_NODE --wd $WDIR --resume --joblog ${COMMANDS}.log < $COMMANDS
