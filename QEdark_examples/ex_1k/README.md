
This documentation aims to help with the reproduction or modification of the results of my thesis. The goal is to calculate the scattering rate of DM-electron scattering. <span style="background-color: transparent;">T</span><span style="background-color: transparent;">he main software used is <a href="http://www.quantum-espresso.org/" title="qe">QuantumESPRESSO</a> and <a href="https://github.com/adrian-soto/QEdark_repo" title="qedark">QEdark</a>:</span>

## HowTo QEdark_1k

This is a short incroduction, in how to use the code QEdark_1k.  
It is a modification to QEdark with the aim of faster execution on Tesla.

All executables can be found under

    /kalinka/home/lukwata/software/My_espresso-5.1.2/bin/

A example is in

    /kalinka/home/lukwata/software/QEdark_1k/QEdark_examples/ex_1k

Copy it to your home folder

    cp -r /kalinka/home/lukwata/software/QEdark_1k/QEdark_examples/ex_1k ~/

and change in it.

## There you find the file dm.in in which you can change the masses of the DM particales and the escape velocity.

There are also the two files `ge.243.in` and `ge.4.in`. In `ge.243.in` a 243 k-point mesh is defined, in `ge.4.in` a 4 k-point mesh is defined (it is really a reduced 27 k-point mesh, but that does not matter here).
A file with a description of the most important variables can be found here.
   * [[%ATTACHURL%/input-description.pdf][input-description.pdf]]: Description of variables in `dm.in` and `ge.##.in`
Let's use the 4 k-point mesh to generate the crystal formfaktor. The steps for `ge.243.in` are the same, but the execution time would be much longer.
For the sake of convenience I created a bash-skript `qedark_1k.sh`.
It lets you control all the following steps.

Before starting the script you need to set the path to the executables and the names of the input-files in `qedark_1k.sys`.
The content of `qedark_1k.sys` looks like this:
```bash
kpoints=243
GE_IN="ge.243.in"
dm_in="dm.in"
# use absolute path
pseudo_dir="/kalinka/home/lukwata/testQEdark/"

job_f2=job_f2.pbs
job_multi=job_multi.pbs

export EXEC_SUM="/kalinka/home/lukwata/software/test_compile/espresso-5.1.2/bin/qedark_1k_sum.x"
export EXEC_F2_1K="/kalinka/home/lukwata/software/test_compile/espresso-5.1.2/PW/src/pw_dark_1k.x"
export EXEC_MULTI="/kalinka/home/lukwata/software/test_compile/espresso-5.1.2/bin/qedark_1k_multimass.x"
export GE_IN
export KPOINTS="$kpoints"
```

There are some old paths for the executables. So change them to the right ones:

```bash
export EXEC_SUM="/kalinka/home/lukwata/software/My_espresso-5.1.2/bin/qedark_1k_sum.x"
export EXEC_F2_1K="/kalinka/home/lukwata/software/My_espresso-5.1.2/PW/src/pw_dark_1k.x"
export EXEC_MULTI="/kalinka/home/lukwata/software/My_espresso-5.1.2/bin/qedark_1k_multi.x"
```

Because we want to use the 4 k-point mesh change the corresponding lines to
```bash
GE_IN="ge.4.in"
kpoints=4
```

The path of the pseudo potential is also wrong. Change it to

    pseudo_dir="/<your homefolder>/ex_1k/"

Now we can start with the calculation of the crystal formfactor. Run
```bash
bash ./qedark_1k.sh
```

The output you will see is:
```
##########################################################
Control the calculation with qedark_1k over this script.
You can set the path to executables, the names of the input files
and the number of k-points in qedark_1k.sys.
Other settings can be done in the dm.in and ge.##.in
##########################################################

1) Create/Update tmp folders for formfactor calculation and copy dm.in, ge.##.in and job_f2.pbs inside.
2) Run pw_dark_1k.x. Start jobs parallel on Tesla.
3) Sum output of parallel jobs. The output-file is C.dat. It contains the formfactor.
4) Run qedark_1k_multimass.x
5) Clean unused output of pw_darK_1k
6) Remove the tmp folders used for formfactor calculation.
7) Quit
Please enter your choice:
```

Choose 1 and hit ENTER. This will create 4 folders `ge_###` , in each a partial sum of the squared formfaktor will be calculated.

Then choose 2. This will send 4 jobs to Tesla. The jobs put ther output in `ge_###`.  
Now you can hit 7 and take a coffee brack until all jobs are finished. To check if they are finished run

    qstat

If you get no output they are finished. If everything went right there should be a file K.###.dat in every ge_### folder.  
If something went wrong look in a folder ge_### for the job_f2.pbs.e######## file.

When changing anything in `dm.in`, `ge.#.in` or `qedark_1k.sys` rerun option 1 with `qedark_1k.sh`.

So finally everything went well, and there are 4 `K.###.dat` files.  
The content of these files need to be summed up.
Therefor run `qedark_1k.sh` and choose 3.
This creates the file `C.dat` which contains the squared formfactor.

You can use 5 or 6 to clean your folder a little bit befor the next step.
But be carefull 6 will delet all the partial sums so be shure you have your C.dat file

`qedark_1k_multimass.x` will use the DM-masses set in `dm.in` and `C.dat` to calculate the scattering for each mass.
To start it, again run `qedark_1k.sh` and choose 4.
It wil only take less then a minute to run then there is a file
`C.###.dat` for every DM-mass.

The `C.###.dat` files contain the scattering rate binned with the resolution set in dm.in.
A snippet of a `C.###.dat` file (less rows then the real one):
```
    01    01   0.139193079E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.390599274E-005
    02    01   0.103375018E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.231179746E-005
    03    01   0.245999056E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.171677492E-005
    01    02   0.139193079E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.390599274E-005
    02    02   0.103375018E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.231179746E-005
    03    02   0.245999056E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.171677492E-005
    01    03   0.139193079E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.390599274E-005
    02    03   0.103375018E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.231179746E-005
    03    03   0.245999056E-002   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.000000000E+000   0.171677492E-005</verbatim>
```

The first column stands for the DM-formfaktor ( 1, 1/q, 1/q^2).
The second column represents the month and so annual modulation of the velocity of the earth
1. December ve_mod=-15kmps
2. March ve_mod=0
3. June ve_mod=15kmps
The third column is the total rate divided by a prefactor the following bins represent the binned rate.
If you have not changed anything else then said above the files should have 503 columns.

## ToDo
- [ ] plot script
