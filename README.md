> edited by Julian Lukwata  
> All instructions concerning my code are labled QEdark_1k  
> 10.08.2017
<!---
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
--->
## 			CONTENTS
* QEdark basic description
* Important info, disclaimer and acknowledgements
* QEdark quick installation guide
* QEdark execution guidelines
<!---
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
--->
##		QEdark

 Code to evaluate and integrate the form factor as in eq (3.16)
 in Essig et.al. --  arXiV:1509.01598 and integrals of this
 quantity in order to calculate the DM-electron scattering rate.


 ### QEdark calculation modes:

 * `calculation_mode='f2'`
 
      gives eq. (4.4) as the output quantity.

 * `calculation_mode='multimass'` 
 
    computes the crystal form factor
    multiplied by the weight `eta(vmin(Ei, Ei'))`, giving an output
    quantity that is proportional to the scattering rate.
<!---
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
--->
## QEdark_1k

 QEdark calculation modes:

* `calculation_mode='f2_1k'`

    computes the fomrfaktor for one fix k-value devided by a prefactor

* `calculation_mode='multimass_f2'`

    computes the the scattering rate devided by a prefactor
    uses as input the file C.dat, which contains the formfaktor.

<!---
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
--->

## IMPORTANT

> QEdark is an extension to the PWscf code of the
> Quantum ESPRESSO package. It has been written to work with
> QE version v.5.1.2. It may work with other versions but this
> is not guranteed.  
> We greatly acknowledge the Quantum ESPRESSO team for making
> their code community accessible. See the following publication
> for details
> P. Giannozzi et.al. J. Phys. Condens. Matter 21, 395502 (2009)

<!---
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
--->

## QEdark quick installation guide

1. Download QuantumESPRESSO v.5.1.2 from http://www.qe-forge.org/gf/download/frsrelease/185/753/espresso-5.1.2.tar.gz

1. Modify `pwscf.f90` by including the code line

        call qedark()
        
    just above the line
    
        CALL stop_run( exit_status )
        
1. Modify `PW/src/Makefile` to include dependencies on two
    new files. Simply add
    
        qedark.o \
        qedark_routines.o \
        qedark_f2.o \
        qedark_multimass.o \
        qedark_read_input.o \

    directly below
    
        wyckoff.o \

 1. Copy the files `qedark.f90, qedark_routines.f90,
    qedark_f2.f90` and `qedark_multimass.f90`
    from QEdark_files into Quantum ESPRESSO's directory
    `PW/src/`


1. Run configure script in QE main directory as

        ./configure --disable-parallel

    For OpenMP multithreaded version, run configure as
    
        ./configure --disable-parallel --enable-openmp

    NOTE:QEdark v1.0 does not support OpenMPI parallelization

1. Compile by typing the following command in QE main directory

        make pw
<!---
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
--->

## QEdark_1k installation guide
> edit Julian Lukwata
> 10.08.2017


1.  Download QuantumESPRESSO v.5.1.2 from
    http://www.qe-forge.org/gf/download/frsrelease/185/753/espresso-5.1.2.tar.gz and unpack it
    
        tar -xvzf espresso-5.1.2.tar.gz

1.  Go in folder `espresso-5.1.2/PW/src/`

        cd espresso-5.1.2/PW/src/
        
    Create the a copie `pwscf_qedark_1k.f90` of `pwscf.f90`
    
        cp pwscf.f90 pwscf_dark_1k.f90

1.  Open the file and rename the program from `pwscf` to `pwscf_dark_1k`,
    you have to do it one time at the beginning and a second time at the end.

    insert the line
    
        call qedark_1k()
       
    directly above
    
        call stop_run( exit_status )

    if there is a line
    
        call qedark()
        
    remove it

1.  You can choose if
    1. just copy the file Makefile into the QuatumESPRESSO folder `PW/src/`.
    
            cp ../../../QEdark_files/Makefile ./
            
        It contains the instructions for QEdark and QEdark_1k
        So you need to do the instruction 2) & 4) of the QEdark installation guide above
    or if you only want QEdark_1k  
    1.  you can edit the Makefile manually.
        Insert in ./Makefile following lines.

        Below
        
            PWOBJS = \
            pwscf.o
            
        insert the lines
        
            PWDARKOBJS = \
            pwscf_dark_1k.o

            DARK1KSUMOBJS = \
            qedark_1k_sum.o

            DARK1KMULTIOBJS = \
            qedark_1k_multimass.o

        Below
        
            wyckoff.o \
            
        insert the lines
        
            qedark_read_input.o \
            qedark_routines.o \
            qedark_1k.o \
            qedark_1k_f2.o \
            qedark_1k_f2_3d.o \
            qedark_1k_f2_sum.o \
            qedark_1k_multimass_sub.o \ 

        Replace
        
            all : tldeps pw.x manypw.x generate_vdW_kernel_table.x generate_rVV10_kernel_table.x
            
        with
        
            all : tldeps pw.x pw_dark_1k.x qedark_1k_sum.x qedark_1k_multimass.x manypw.x generate_vdW_kernel_table.x generate_rVV10_kernel_table.x

        Below the expression
        
            pw.x:...
            
        insert the lines
        
            pw_dark_1k.x : $(PWDARKOBJS) libpw.a $(LIBOBJS) $(QEMODS)
                    $(LD) $(LDFLAGS) -o $@ \
                       $(PWDARKOBJS) libpw.a $(QEMODS) $(LIBOBJS) $(LIBS)
                    - ( cd ../../bin; ln -fs ../PW/src/$@ . )

            qedark_1k_sum.x : $(DARK1KSUMOBJS) libpw.a $(LIBOBJS) $(LIBS)
                    $(LD) $(LDFLAGS) -o $@ \
                            $(DARK1KSUMOBJS) libpw.a $(QEMODS) $(LIBOBJS) $(LIBS)
                    - ( cd ../../bin; ln -fs ../PW/src/$@ . )

            qedark_1k_multimass.x : $(DARK1KMULTIOBJS) libpw.a $(LIBOBJS) $(LIBS)
                    $(LD) $(LDFLAGS) -o $@ \
                            $(DARK1KMULTIOBJS) libpw.a $(QEMODS) $(LIBOBJS) $(LIBS)
                    - ( cd ../../bin; ln -fs ../PW/src/$@ . )

         Save the file.


1.  Copy the files `qedark_1k.f90, qedark_1k_f2.f90, qedark_1k_f2_3d.f90, qedark_f2_1k_sum.f90, qedark_1k_multimass_sub.f90` and `qedark_1k_multimass.f90`
    form QEdark_files into the QE folder `PW/src/`.
    
        cd ../../../QEdark_files
        cp qedark_1k* qedark_routines.f90 qedark_read_input.f90 ../espresso-5.1.2/PW/src/

1.  Run the configure script in QE main directory as

        cd ../espresso-5.1.2/
        ./configure --disable-parallel


1.  Compile by typing the following command in QE main directory

        make pw

<!---
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
--->

##         QEdark QUICK EXECUTION GUIDELINES

 1) The file dm.in must exist in the directory where the code
    is being run since it contains the input data to QEdark.

 2) QEdark is executed within Quantum ESPRESSO's pw.x program
    which is copied to the bin/ directory. Execute as
       QEdark_dir/bin/pw.x < QEinputfile > QEoutputfile

 3) In order to run the multithreaded version, the value of
    OMP_STACKSIZE will most likely need to be increased. For
    Si and Ge we have found that setting it to 512M has been
    enough.

    In bash this can be done by entering the commmand
       export OMP_STACKSIZE=512M

  4) Users are encouraged to get familiar with the code by
     going through the examples in the folder QEdark_examples/,
     which have been chosen to introduce the user into several
     features of QEdark.

  5) Previous knowledge of Quantum ESPRESSO
     is highly recommended. There are plenty of tutorials and
     examples at the link
     	 http://www.quantum-espresso.org/tutorials/

<!---
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
--->
## QEdark_1k execution

Look inside the folder [QEdark_examples/ex_1k/](QEdark_examples/ex_1k).

  
