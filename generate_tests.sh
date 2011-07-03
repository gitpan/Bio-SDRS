#!/bin/sh -x

make

blib/script/sdrs.pl -multiple=1.05 -step=20 -ldose=0.4 -hdose=25000 -trim=0 -outdir=t t/OVCAR4_HCS_avg.txt
cd t
for f in sdrs.1.05.20.EC50.out sdrs.1.05.20.out sdrs.pval_FDR.out sdrs.sorted_probes.out
do
    mv ${f} ref.${f}
    sort ref.${f} >ref.${f}.srt
done

