#!/bin/bash
#
# Partition:
#SBATCH --partition=savio2_bigmem

# Processors:
#SBATCH --nodes=1
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
# Example - how to run this SLURM script:
#
#cd /global/scratch/users/deboraycb/argbias/simulations_3/relate/m2e9
#sbatch --job-name=RelM2e9 --output=%A_RelM2e9.out --export=PREF=sim_m2e9_0,NSPL=8,RECMAP=../input/sim_std.map,Ne=10000,MUT=2e-9,s0=0,LEN=100000000 /global/scratch/users/deboraycb/argbias/simulations_3/scripts/relate/2_run_relate2plots.sh
#
#cd /global/scratch/users/deboraycb/argbias/simulations_3/relate/r2e7
#sbatch --job-name=RelR2e7 --output=%A_RelR2e7.out --export=PREF=sim_r2e7_0,NSPL=8,RECMAP=../input/sim_r2e7.map,Ne=10000,MUT=2e-8,s0=0,LEN=100000000 /global/scratch/users/deboraycb/argbias/simulations_3/scripts/relate/2_run_relate2plots.sh
#
#cd /global/scratch/users/deboraycb/argbias/simulations_3/relate/l5mb
#sbatch --job-name=RelL5mb --output=%A_RelL5mb.out --export=PREF=sim_l5mb_0,NSPL=8,RECMAP=../input/sim_std.map,Ne=10000,MUT=2e-8,s0=0,LEN=5000000 /global/scratch/users/deboraycb/argbias/simulations_3/scripts/relate/2_run_relate2plots.sh
#
#cd /global/scratch/users/deboraycb/argbias/simulations_3/relate/l250kb
#sbatch --job-name=RelL250kb --output=%A_RelL250kb.out --export=PREF=sim_l250kb_0,NSPL=8,RECMAP=../input/sim_std.map,Ne=10000,MUT=2e-8,s0=0,LEN=250000 /global/scratch/users/deboraycb/argbias/simulations_3/scripts/relate/2_run_relate2plots.sh
#
## Command(s) to run:

RELDIR=/global/home/users/deboraycb/scratch/opt/relate_v1.1.2_x86_64_dynamic/
REP=0

date
$RELDIR/bin/RelateFileFormats --mode ConvertFromVcf -i ../../vcf/$PREF --haps ../input/$PREF.haps --sample ../input/$PREF.sample
date
$RELDIR/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps ../input/$PREF.haps --sample ../input/$PREF.sample --ancestor ../input/sim_hudson.anc -o ../input/${PREF}
date
gzip ../input/${PREF}.haps
gzip ../input/${PREF}.sample
date
$RELDIR/bin/Relate --mode All -m $MUT -N $((2*Ne)) --haps ../input/$PREF.haps.gz --sample ../input/$PREF.sample.gz --map $RECMAP -o $PREF
date

CHUNKLEN=$((LEN/20))

E=0
for i in {1..20}
do 
    S=$E
    E=$((i*CHUNKLEN))
    if [ -f "${PREF}_${S}-${E}_spl1000branch.dist" ]; then
        echo '---> skipping chunk already sampled: '${PREF}_${S}-${E}_spl1000branch
    else
        echo '---> sampling chunk: '${PREF}_${S}-${E}_spl1000branch
    $RELDIR/scripts/SampleBranchLengths/SampleBranchLengths.sh --mode All -m $MUT -i $PREF --map $RECMAP --num_samples 1000 --first_bp ${S} --last_bp ${E} --format a --coal ../input/sim_std.coal -o ${PREF}_${S}-${E}_spl1000branch &
    fi
done
wait
date

echo '---> concat chunks mut'
S=0
E=$CHUNKLEN
cat ${PREF}_${S}-${E}_spl1000branch.mut > ${PREF}_spl1000branch.mut
for i in {2..20}
do 
    S=$E
    E=$((i*$CHUNKLEN))
    tail -n +2 ${PREF}_${S}-${E}_spl1000branch.mut >> ${PREF}_spl1000branch.mut
done

echo '---> concat chunks anc'
totlines=0
for f in ${PREF}_*-*_spl1000branch.anc
do nlines=$(awk 'NR==2 {print $2;exit}' $f)
    totlines=$((totlines+nlines))
done
echo 'NUM_HAPLOTYPES '$NSPL > ${PREF}_spl1000branch.anc
echo 'NUM_TREES '$totlines >> ${PREF}_spl1000branch.anc
echo 'NUM_SAMPLES_PER_TREE 1000' >> ${PREF}_spl1000branch.anc
E=0
for i in {1..20}
do 
    S=$E
    E=$((i*CHUNKLEN))
    tail -n +4 ${PREF}_${S}-${E}_spl1000branch.anc >> ${PREF}_spl1000branch.anc
done 
date

for s1 in $(seq $s0 $((NSPL-2)))
do for s2 in $(seq $((s1+1)) $((NSPL-1)))
    do 
    if [ -s "relate_${PREF}_${s1}-${s2}_Tcpair_thinned.txt.gz" ]; then
        echo '---> skipping pair already done and thinned: '$s1 $s2
    else
        if [ -s "relate_${PREF}_${s1}-${s2}_Tcpair.txt" ]; then
            echo '---> thinning relate_'${PREF}'_'${s1}'-'${s2}'_Tcpair.txt'
            Rscript thin.R relate_${PREF}_${s1}-${s2}_Tcpair.txt $((2*Ne))
        else    
            echo '---> outputCoalTime samples '$s1 $s2
            ../outputCoalTime relate_${PREF}_${s1}-${s2}_Tcpair.txt ../../tcoalmsp/rep${REP}/${PREF}_spls${s1}-${s2}.tc ${PREF}_spl1000branch.anc ${PREF}_spl1000branch.mut $Ne $s1 $s2 
            echo '---> thinning relate_'${PREF}'_'${s1}'-'${s2}'_Tcpair.txt'
            Rscript thin.R relate_${PREF}_${s1}-${s2}_Tcpair.txt $((2*Ne))
        fi
        rm relate_${PREF}_${s1}-${s2}_Tcpair.txt
        gzip -f relate_${PREF}_*-*_Tcpair_thinned.txt
    fi
done
done
wait

rm relate_${PREF}_Tcallpair_thinned.txt
rm relate_${PREF}_Tcallpair_thinned_ranks.txt
for s1 in $(seq 0 $((NSPL-2)))
do for s2 in $(seq $((s1+1)) $((NSPL-1)))
    do 
    echo 'cat pairs thinned'$s1 $s2
    zcat relate_${PREF}_${s1}-${s2}_Tcpair_thinned.txt.gz >> relate_${PREF}_Tcallpair_thinned.txt
    cat relate_${PREF}_${s1}-${s2}_Tcpair_thinned_ranks.txt >> relate_${PREF}_Tcallpair_thinned_ranks.txt
done
done

date
echo 'organize all pairs thinned ranks'
../organizeQuantileOutput relate_${PREF}_Tcallpair_thinned_ranks.txt relate_${PREF}_Tcallpair_thinned_ranks_org.txt 101        
date
echo 'plot ranks'
Rscript ../plot_ranks.R -n 101 -i relate_${PREF}_Tcallpair_thinned_ranks_org.txt -s $NSPL
Rscript ../calc_kld_ranks.R -i relate_${PREF}_Tcallpair_thinned_ranks_org.txt
date
echo 'plot Tc'
Rscript plotTc_relate.R -i relate_${PREF}_Tcallpair_thinned.txt -N $((2*Ne))
date
echo 'plot point estimates heatmap'
Rscript plot_pointest.R -i relate_${PREF}_Tcallpair_thinned.txt
date
wait
