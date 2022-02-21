import argparse

parser = argparse.ArgumentParser(description = "convert April's msprime coalTime into bed format, converting start and end positions to int")
parser.add_argument('tc', help='file with true msprime TMRCAs')
parser.add_argument('Ne', help='2Ne to divide coal time in generations')
args = parser.parse_args()
tc_msp_orig = args.tc
Ne=float(args.Ne)

# 1. convert msprime tcoal into bed format, converting start and end positions to int
tc_msp_bed = '.'.join(tc_msp_orig.split('.')[:-1])+'.bed'
infh = open(tc_msp_orig, 'r')
outfh = open(tc_msp_bed, 'w')
for line in infh:
    fields = line.strip().split()
    chrstart = int(round(float(fields[0])))
    chrend = int(round(float(fields[1])))
    coaltime = float(fields[2])/Ne
    outfh.write('1\t'+str(chrstart)+'\t'+str(chrend)+'\t'+str(coaltime)+'\n')

infh.close()
outfh.close()

