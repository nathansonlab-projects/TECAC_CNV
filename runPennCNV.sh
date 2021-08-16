#!/bin/bash

# Run penncnv commands by specifying .pfb file, .list file and the output name. Please upload this script to the penncnv-1.0.5 folder, and run it from there.
#Example: ./runPennCNV.sh Infinium24HumanCore_KLN.pfb ListOfSamples.txt penncnv_run3
PFB=$1
LIST=$2
OUTNAME=$3


perl detect_cnv.pl -test -hmm lib/exome.hmm -pfb $PFB --list $LIST -log ${OUTNAME}.log -out ${OUTNAME}.rawcnv