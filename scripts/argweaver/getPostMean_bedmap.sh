#!/usr/bin/bash

PREF=$1
DIR=$2
TRUEPREF=$3
MCSPL=$4
s1=$(echo $5| cut -f 1 -d ' ')
s2=$(echo $5| cut -f 2 -d ' ')

EST=${DIR}'/arg-sample_sim_'${PREF}'_'$s1'-'$s2'_Tcpair_spl'${MCSPL}'.bed'
PART=${DIR}'/arg-sample_sim_'${PREF}'_'$s1'-'$s2'_Tcpair_spl'${MCSPL}'_partitions.bed'
OUT=${DIR}'/arg-sample_sim_'${PREF}'_'$s1'-'$s2'_Tcpair_spl'${MCSPL}'_posterior.bed'
TRUE=${TRUEPREF}$s1'-'$s2'.bed'
echo $s1'-'$s2
if [[ -s "$OUT" ]]; then
    printf "$OUT is not empty\n"
    ls -l $OUT
else
    bedmap --delim '\t' --echo --echo-map-id $PART $TRUE | bedmap --delim '\t' --echo --mean --median - $EST | awk -v OFS='\t' '{print $1,$2,$3,$4,$5/20000,$6/20000}' > $OUT
fi
