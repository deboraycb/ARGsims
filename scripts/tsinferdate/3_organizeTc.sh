NSPL=$1
sim=$2
pref=$3

#NSPL=8
#sim='std'
#pref='std_priorDef'
#pref='std_prior100tp'
#pref='std_priorGeomMax12'

cd $sim

for s1 in $(seq 0 $((NSPL-2)))
do for s2 in $(seq $((s1+1)) $((NSPL-1)))
do 
    bedops --element-of <(zcat sim_${sim}_0_spls${s1}-${s2}.bed.gz) <(zcat sim_${pref}_spls${s1}-${s2}_inf.bed.gz) | bedops --partition - <(zcat sim_${pref}_spls${s1}-${s2}_inf.bed.gz) | awk '$2!=$3{print}' | bedmap --echo --echo-map-id - <(zcat sim_${sim}_0_spls${s1}-${s2}.bed.gz) | bedmap --echo --echo-map-id - <( zcat sim_${pref}_spls${s1}-${s2}_inf.bed) | sed 's/|/\t/g' | gzip -c > sim_${pref}_spls${s1}-${s2}_post.bed.gz
done
done

zcat sim_${pref}_spls*_inf.bed.gz | gzip -c > sim_${pref}_allspls_inf.bed.gz
