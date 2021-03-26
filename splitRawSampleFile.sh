#!/bin/bash
for ((n=1; n<=14159; n++))
do
SPL=$n
let start=3*$n+1
let end=3*$n+3
cut -d "," -f 1-3,$start-$end /project/knathanslab/TECAC/CNV/tecac_all_cnv.csv > /project/knathanslab/TECAC/CNV/ByBID/sample${SPL}.txt &
if (( (n % 300) == 0 )); then
wait
fi
done
