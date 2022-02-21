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
#SBATCH --account=ac_popgen
#
# Specify Faculty Computing Allowance:
#SBATCH --qos=savio_normal
#
# to pass bash variables:
# sbatch --job-name=$A.$b.run --output=$A.$b.out --export=A=$A,b=$b jobscript.sbatch
#
#DIR=m2e9hud_smcp; PREF=sim_m2e9_0; PREF1=sim_m2e9_0; REP=0; NSPL=8; BURN=200; NITER=1200; s0=0; LEN=100000000
#DIR=r2e7hud_smcp; PREF=sim_r2e7_0; PREF1=sim_r2e7_0; REP=0; NSPL=8; BURN=200; NITER=1200; s0=0; LEN=100000000
#DIR=l250kb_smcp; PREF=sim_l250kb_0; PREF1=sim_l250kb_0; REP=0; NSPL=8; BURN=200; NITER=1200; s0=0; LEN=250000
#DIR=l5mb_smcp; PREF=sim_l5mb_0; PREF1=sim_l5mb_0; REP=0; NSPL=8; BURN=200; NITER=1200; s0=0; LEN=5000000
#cd /global/scratch/users/deboraycb/argbias/simulations_3/argweaver/${DIR}/rep${REP}/
#sbatch --job-name=ARGw${PREF}ChunksPlt --output=%A_ARGw${PREF}ChunksPlt.out --export=NSPL=$NSPL,PREF=${PREF},PREF1=$PREF1,REP=$REP,Ne=10000,BURN=$BURN,NITER=$NITER,s0=$s0,LEN=$LEN /global/scratch/users/deboraycb/argbias/simulations_3/scripts/argw/4b_arg-summarize2plot_allchunks.sh
#
# Command(s) to run:
module load r-packages
module load r-spatial

### make list of sample names
rm allspls.txt
for i in $(seq 0 $(((NSPL-1)/2)))
do
    echo -e "spl${i}_1\nspl${i}_2" >> allspls.txt
done
#
### convert argweaver smc files to bed format and plots per rep
E=1
CHUNKLEN=$((LEN/20))
for i in {1..20}
do 
    S=$E
    E=$((i*CHUNKLEN))
    SMC=arg-sample_${PREF}_${S}-${E}
    if [ -f "${SMC}.bed.gz" ]; then
        echo '---> skipping rep already done smc2bed: '$SMC
    else
        echo '---> .smc 2 .bed: '$SMC
        smc2bed-all $SMC &
    fi
done
wait
for s1 in $(seq $s0 $((NSPL-2)))
do for s2 in $(seq $((s1+1)) $((NSPL-1)))
do
    awk -v s1=$s1 -v s2=$s2 'NR==s1+1; NR==s2+1' allspls.txt > spl.txt
    BED=arg-sample_${PREF}_${s1}-${s2}_Tcpair.bed
    ## summarize tmrca for all chunks in each pairs of samples
    echo '=============================================='
    if [ -s "$BED" ]; then
        echo '---> skipping pair already summarized: '${s1}'-'${s2}
    else
        echo '---> arg-summarize pair of samples: '${s1}'-'${s2}
        E=1
        for i in {1..20}
        do 
            S=$E
            E=$((i*CHUNKLEN))
            SMC=arg-sample_${PREF}_${S}-${E}
            arg-summarize -a ${SMC}.bed.gz --tmrca --subset spl.txt | grep -v '^#' >> $BED
        done
    fi
    wait
    # select iterations to plot        
    TOPLOT=arg-sample_${PREF}_${s1}-${s2}_Tcpair_spl${BURN}-${NITER}
    # get coal time ranks
    TCSIMFILE=../..//tcoalmsp/rep${REP}/${PREF1}_spls${s1}-${s2}
    if [ ! -f "${TCSIMFILE}.bed" ]; then
        python MSPcoalTime2bed.py ${TCSIMFILE}.tc $((2*Ne))
    fi
    if [ ! -s "${TOPLOT}.bed" ]; then
        awk -v BURN=$BURN -v NITER=$NITER '$4>BURN && $4<=NITER' $BED > ${TOPLOT}.bed
    fi
    NRANKS=$(((NITER-BURN)/10+1))
    if [ ! -s "${TOPLOT}_ranks.bed" ]; then
        echo "---> get ranks: python ~/scratch/argbias/simulations_3/scripts/argw/coalTimeRank.py ${TOPLOT}.bed ${TCSIMFILE}.bed $((NRANKS-1))" 
        python coalTimeRank.py ${TOPLOT}.bed ${TCSIMFILE}.bed $((NRANKS-1))
    fi
done
done

# concatenate all pw tmrca into one file and plot
echo '---> cat all pairs Tc'
echo '---> cat all pairs ranks'
ALLPAIRS=arg-sample_${PREF}_Tcallpair_spl${BURN}-${NITER}
rm ${ALLPAIRS}.bed
rm tmp1.bed
rm tmp2.bed
for s1 in $(seq 0 $((NSPL-2)))
do for s2 in $(seq $((s1+1)) $((NSPL-1)))
do
cat arg-sample_${PREF}_${s1}-${s2}_Tcpair_spl${BURN}-${NITER}.bed >> ${ALLPAIRS}.bed
cat arg-sample_${PREF}_${s1}-${s2}_Tcpair_spl${BURN}-${NITER}_ranks.bed >> tmp1.bed
cat arg-sample_${PREF}_${s1}-${s2}_Tcpair_spl${BURN}-${NITER}_ranks_newname.bed >> tmp2.bed
done
done
sort-bed tmp1.bed > ${ALLPAIRS}_ranks.bed
sort-bed tmp2.bed > ${ALLPAIRS}_ranks_newname.bed
rm tmp1.bed
rm tmp2.bed

echo '=============================================='
echo '---> plot Tc'
Rscript plotTc_argw.R -i ${ALLPAIRS}.bed -N $((2*Ne)) &
echo '---> organize all pairs ranks'
cut -f 4- ${ALLPAIRS}_ranks_newname.bed > ${ALLPAIRS}_ranks.txt
NRANKS=$(((NITER-BURN)/10+1))
../organizeQuantileOutput ${ALLPAIRS}_ranks.txt ${ALLPAIRS}_ranks_org.txt $NRANKS
echo '---> plot ranks'
Rscript ../plot_ranks.R -n $NRANKS -i ${ALLPAIRS}_ranks_org.txt -s $NSPL -l $LEN &
Rscript ../calc_kld_ranks.R -i ${ALLPAIRS}_ranks_org.txt
wait
