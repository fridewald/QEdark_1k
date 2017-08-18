#!/bin/bash

# a script for creating kpoints folders
# addes the corresopnign k_value to the dm.in file in that folder

kpoints=243
input=ge.243
dm_in=dm.in
job_in=job.pbs
for k in $(seq -f %03g 1 $kpoints); do
    outdir="${input}_$k"
    if [ -d $outdir ]; then
        echo "dir $outdir exists"
    else
        echo "creating dir $outdir"
        mkdir $outdir
    fi

    sed -e "/\&dm_parameters/a\    k_value = $k" $dm_in  | cat   > $outdir/$dm_in
    cat < ${input}.in > $outdir/${input}.in
    cat < $job_in > $outdir/$job_in
done
