program qedark_1k_sum
    !
    !
    !USE wavefunctions_module,           ONLY: evc                  ! For collinear, evc(npwx, nbnd) [look at allocate_wfc.f90]
    USE kinds,                          ONLY: DP

    !use omp_lib

    IMPLICIT NONE

    CHARACTER(20), PARAMETER :: infile="dm.in"                      ! Name of qedark input file


    CHARACTER(20) :: calculation_mode
    LOGICAL :: restart_mode
    LOGICAL :: f1exists, f2exists, f3exists

    REAL(DP) :: vearth_SI                                           ! Earth's average velocity around the galactic center in SI units
    REAL(DP) :: vesc_SI                                             ! Earth's escape velocity in SI units
    REAL(DP) :: deltav_SI(4)                                        ! Change in Earth's velocity over 4 seasons (spring, summer, fall winter) in SI
    REAL(DP) :: v0_SI                                               ! DM typical velocity in SI


    INTEGER, PARAMETER :: max_num_mx=99
    INTEGER :: num_mx
    REAL(DP) :: mx_NU(max_num_mx)                                   ! Dark matter particle mass in NU

    INTEGER :: er_bin_type                                          ! Type of bins for integrated form factors
    INTEGER :: num_er_bins                                          ! Number of bins for integrated form factors

    REAL(DP) :: ermax_NU                                            ! Max recoil energy cutoff in eV
    REAL(DP) :: er_binsize                                          ! Bin size (Ry) for linear bins with user-provided size

    INTEGER :: num_q_bins
    REAL(DP) :: deltaq

    INTEGER :: nksf                                                 ! Number of k-points in formfactor calculation.

    INTEGER :: nbndval, nbndcond                                    ! Number of occupied and unoccupied Kohn-Sham orbitals. TO BE SET BY USER!

    LOGICAL :: scissor_correction
    REAL(DP) :: scissorgap                                          ! Band gap (in eV) for scissor operator corrected band energies

    INTEGER :: k_value                                              ! k-value for which to calculate, only for ..._1k 
    !real(DP) :: kk_SI                                                   ! norm-factor for eta

    INTEGER :: ierr


    CALL start_clock( ' QEdark_1k_sum ')
    print *, " "
    print *, " "
    print *, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
    print *, "--- Entering QEdark_1k_sum"
    print *, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
    print *, " "
    print *, " "


    ! Read parameters from file
    WRITE(*,*), " "
    CALL qedark_read_input(&
       infile, calculation_mode, restart_mode, &
       nksf, nbndval, nbndcond, &
       vearth_SI, vesc_SI, v0_SI, deltav_SI, &
       num_mx, mx_NU, &
       Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
       num_q_bins, deltaq, &
       scissor_correction, scissorgap, &
       k_value)


    IF (num_mx > max_num_mx) &
         CALL errore('QEdark', 'Number of mass points exceeds max_num_mx.', ABS(ierr) )

    IF (calculation_mode == "f2_1k") THEN
        CALL qedark_1k_f2_sum( &
            nksf, &
            num_er_bins, &
            num_q_bins)

    ELSE
       print *, "Error"
       CALL errore('QEdark_1k_sum', 'calculation_mode not recognized. EXITING.', ABS(ierr) )
    END IF


    print *, " "
    print *, " "
    print *, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
    print *, "-- Exiting QEdark_1k_sum"
    print *, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
    print *, " "
    print *, " "


    CALL stop_clock( ' QEdark_1k_sum ')

end program qedark_1k_sum
