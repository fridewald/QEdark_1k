#!/bin/bash
#

# set variables
. ./qedark_1k.sys
basedir=$PWD


function creation {
    echo "Update ${dm_in}"
    grep -q 'calculation_mode' $dm_in && sed -i "s/calculation_mode.*/calculation_mode='f2_1k'/" $dm_in || sed -i "/\&dm_parameters/a\    calculation_mode='f2_1k'" $dm_in
    grep -q 'nksf' $dm_in && sed -i "s/nksf.*/nksf=$kpoints/" $dm_in || sed -i "/\&dm_parameters/a\    nksf=$kpoints" $dm_in

    echo "Create/Update folders."
    for k in $(seq -f %03g 1 $kpoints); do
        outdir="ge_$k"
        if [ ! -d $outdir ]; then
            echo "creating dir $outdir"
            mkdir $outdir
        fi

        sed -e "/\&dm_parameters/a\    k_value = $k" $dm_in  | cat   > $outdir/$dm_in
        sed "s|pseudo_dir.*|pseudo_dir= '${pseudo_dir}'|" $GE_IN | cat > $outdir/$GE_IN
        #sed "s|GE_IN.*|GE_IN=${GE_IN}|" $job_f2  | cat > $outdir/$job_f2
        cat < $job_f2 > $outdir/$job_f2
    done
}

function start_jobs {
    for k in $(seq -f %03g 1 $kpoints); do
        if [ -d "$basedir/ge_$k" ]; then
            cd $basedir/ge_$k
            #echo `pwd`
            qsub -V -cwd ./$job_f2 -q all.q
        fi
    done
}

function sum {
    nice $EXEC_SUM > sum.out
}

function multimass {
    echo "Update ${dm_in}"
    grep -q 'calculation_mode' $dm_in && sed -i "s/calculation_mode.*/calculation_mode='multimass_f2'/" $dm_in || sed -i "/\&dm_parameters/a\    calculation_mode='multimass_f2'" $dm_in

    qsub -V -cwd ./$job_multi -q all.q
}

function clean {
    find . -name "germ.save" -exec rm -r {} \;
    find . -name "germ.igk" -exec rm {} \;
    find . -name "germ.mix" -exec rm {} \;
    find . -name "core.[0-9]*" -exec rm {} \;
    find . -name "job.pbs.*" -exec rm {} \;
}

function veryclean {
    echo ge_[0-9]*
    for remove_dir in "ge_[0-9]*" ; do
       rm -rf $remove_dir
    done
}

PS3='Please enter your choice: '
message=\
'##########################################################\n'\
'Control the calculation with qedark_1k over this script.\n'\
'You can set the path to executables, the names of the input files\n'\
'and the number of k-points in qedark_1k.sys.\n'\
'Other settings can be done in the dm.in and ge.##.in\n'\
'##########################################################\n'
options=("Create/Update tmp folders for formfactor calculation and copy ${dm_in}, ${GE_IN} and ${job_f2} inside."
         "Run pw_dark_1k.x. Start jobs parallel on Tesla."
         "Sum output of parallel jobs. The output-file is C.dat. It contains the formfactor."
         "Run qedark_1k_multimass.x"
         "Clean unused output of pw_darK_1k"
         "Remove the tmp folders used for formfactor calculation."
         )

echo -e $message
select opt in "${options[@]}" "Quit"
do
    case "$REPLY" in
        1 )
            creation
            ;;
        2 )
            echo "Start jobs."
            start_jobs
            ;;
        3 )
            echo "Start summation."
            sum
            ;;
        4 ) 
            echo "Start multimass."
            multimass
            ;;
        5 )
            echo "Clean"
            clean
            ;;
        6 )
            echo "Remove folders"
            veryclean
            ;;
        $(( ${#options[@]}+1 )) )
            break
            ;;
        *) echo invalid option;;
    esac
done
