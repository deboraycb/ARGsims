import subprocess
import argparse
import numpy as np
import sys

parser = argparse.ArgumentParser(description = "Rank true tcoal in a set of ARGweaver samples tcoal")
parser.add_argument('tc1', help='bed file with ARGweaver TMRCAs')
parser.add_argument('tc2', help='bed file with true msprime TMRCAs')
parser.add_argument('r', help='number of samples/max rank')
args = parser.parse_args()
tc_argw = args.tc1
tc_msp_bed = args.tc2
maxrank=int(args.r)

# 1. partition all regions in msprime tcoal file (above) and in ARGweaver tcoal file
tc_argw_pref = '.'.join(tc_argw.split('.')[:-1])
partitionsfile = tc_argw_pref+'_partitions.bed'
partitionsfh = open(partitionsfile, 'w')
p1 = subprocess.Popen(['bedops', '--element-of', tc_msp_bed , tc_argw], stdout=subprocess.PIPE)
p2 = subprocess.Popen(['bedops', '--partition', '-', tc_argw], stdin=p1.stdout, stdout=partitionsfh)
p1.stdout.close()
p2.communicate()
partitionsfh.close()

# 2. paste back in the coal times for all partitions
tc_partitions = tc_argw_pref+'_ranks.bed'
#tc_part_check_fh = open(tc_argw_pref+'_CHECKranks.bed', 'w')
tc_part_fh = open(tc_partitions, 'w')
#tc_partitions_newname = 'partTc_newname'+tc_argw
tc_partitions_newname = tc_argw_pref+'_ranks_newname.bed'
tc_part_fh_newname = open(tc_partitions_newname, 'w')
res_bedmap_argw = subprocess.run(['bedmap', '--echo', '--echo-map-score', partitionsfile, tc_argw], stdout=subprocess.PIPE, text=True)
tc_part_argw = res_bedmap_argw.stdout.split()
res_bedmap_msp = subprocess.run(['bedmap', '--echo', '--echo-map-id', partitionsfile, tc_msp_bed], stdout=subprocess.PIPE, text=True)
tc_part_msp = res_bedmap_msp.stdout.split()
for i in range(2, len(tc_part_argw), 3):
    start = int(tc_part_argw[i-1])
    field = tc_part_argw[i].split('|')
    end = int(field[0])
    tc_list = field[1].split(';')
    if len(tc_list) != maxrank:
        print('wrong len '+str(len(tc_list))+ ' at line:'+ '\t'.join(tc_part_argw[i-2:i])+'\n'+tc_argw)
        tcs_argw = [float(x)/20000 for x in tc_list[0:maxrank]]
    else:
        tcs_argw = [float(x)/20000 for x in tc_list]
    tc_msp = float(tc_part_msp[i].split('|')[1])
    round_idxs = np.argmin(np.array([abs(x-tc_msp) for x in tcs_argw]))
    tc_msp_rounded = tcs_argw[round_idxs]
    rank = sum([x<tc_msp for x in tcs_argw])
    binnames = {0.000000000:0.0005534663, 0.002459674:0.0026082927, 0.006129347:0.0063502486, 0.011604261:0.0119319877, 0.019772475:0.0202573192, 0.031958917:0.0326731418, 0.050140295:0.0511855456, 0.077265725:0.0787799306, 0.117735099:0.1198935069, 0.178112767:0.1811089117, 0.268192311:0.2721639318, 0.402585118:0.4074025242, 0.603090426:0.6078185928, 0.902231273:0.9038469483, 1.348529916:1.3390060470, 2.014378397:1.9744168991, 3.007780923:2.8945134780, 4.489872730:4.2165786186, 6.701057088:6.1134094508, 10.000000000:8.728228632}
    # first, round argweaver coal times to nearest value in my dict
    oldbinnames = [] 
    for tcargw in tcs_argw:
        diffs = [abs(float(tcargw)-x) for x in binnames.keys()]
        keyindex = np.argmin(np.array(diffs))
        if min(diffs) > sorted(binnames.keys())[1]:
            raise Exception('dont know how to rename this interval: '+str(tcargw)+'. Closest one is '+str(binnames.keys()[keyindex]))
        oldbinnames.append(list(binnames.keys())[keyindex])
    newbinnames = []
    for k in oldbinnames:
        newbinnames.append(binnames[k])
    newrank = sum([x<tc_msp for x in newbinnames])
    length = end-start
    tc_part_fh.write('\t'.join(['1', str(start), str(end), str(rank), str(length), str(tc_msp)])+'\n')
    tc_part_fh_newname.write('\t'.join(['1', str(start), str(end), str(newrank), str(length), str(tc_msp)])+'\n')
tc_part_fh.close()
tc_part_fh_newname.close()
