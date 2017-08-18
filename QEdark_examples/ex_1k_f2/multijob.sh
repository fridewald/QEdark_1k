#!/bin/bash

# start all the jobs in the created folder

kpoints=243
basedir=$PWD
input=ge.243
for k in $(seq -f %03g 1 $kpoints); do
    if [ -d "$basedir/${input}_$k" ]; then
        cd $basedir/${input}_$k
        echo `pwd`
        #qsub -cwd ./job.pbs -q all.q
    else
        echo "Folder $basedir/${input}_$k does not exist"
    fi
done
