!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Julian Lukwata
!  04-04-2017
!  KIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
subroutine qedark_1k_multimass_sub( &
        vearth_SI, vesc_SI, v0_SI, deltav_SI, &
        num_mx, mx_NU, &
        Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
        numqbins, deltaq)

    !
    ! Evaluate and print the INTEGRATED form factors (scattering rate up to
    ! constant) squared weighted by the eta function for multiple DM masses.
    !
    !
    !USE wavefunctions_module,           ONLY: evc    ! For collinear, evc(npwx, nbnd) [look at allocate_wfc.f90]
    USE kinds,                          ONLY: DP

    !use omp_lib


    IMPLICIT NONE

    INTEGER, EXTERNAL :: find_bin
    REAL(DP), EXTERNAL :: vmin, eta, bzvol, eucl_norm_fast, find_ehomo, find_elumo, knorm


    ! Variables specific to int_formfact(..)
    REAL(DP), PARAMETER :: Ry2eV = 13.60569253_DP                   !
    REAL(DP), PARAMETER :: twobyalpha=274.07199814110116            ! ==2/alpha, which is the conversion factor for velocities from N.U. to Rydberg A.U.
    REAL(DP), PARAMETER :: speedoflight=3.0E08_DP                   ! Speed of light in SI units
    REAL(DP), PARAMETER :: RAU2kmps=(speedoflight/1000.0_DP)/twobyalpha
    real(DP), parameter :: alphame=3729.93



    !REAL(DP), PARAMETER :: vearth_SI = 2.40E5_DP                   ! Earth's average velocity around the galactic center in SI units
    REAL(DP) :: vearth_SI                                           ! Earth's average velocity around the galactic center in SI units
    REAL(DP) :: vearth_RAU !=(twobyalpha/speedoflight)*vearth_SI    ! Earth's average velocity in RAU
    REAL(DP) :: vearth_kmps !=vearth_SI/1000.0                      ! Earth's average velocity in km/s

    REAL(DP) :: vesc_SI                                             ! Earth's escape velocity in SI units
    REAL(DP) :: vesc_kmps != vesc_SI/1000.0                         ! Earth's escape velocity in km/s
    REAL(DP) :: vesc_RAU !=(twobyalpha/speedoflight)*vesc_SI        ! Earth's escape velocity in RAU

    REAL(DP) :: v0_SI                                               ! DM typical velocity in SI
    REAL(DP) :: v0_kmps !=v0_SI/1000.0                              ! DM typical velocity in km/s
    REAL(DP) :: v0_RAU !=(twobyalpha/speedoflight)*v0_SI            ! DM typical velocity in RAU

    INTEGER :: idmff, imonth
    INTEGER :: nmonths = 3
    REAL(DP) :: deltav_SI(3)
    REAL(DP) :: deltav_kmps(3) !=v0_SI/1000.0
    REAL(DP) :: deltav_RAU(3) !=(twobyalpha/speedoflight)*v0_SI

    REAL(DP) :: me_NU = 0.510998910E6_DP                            ! Electron mass in NU
    INTEGER, PARAMETER :: max_num_mx=99
    INTEGER :: imx, num_mx
    REAL(DP) :: mx_NU(max_num_mx) != 100.0E6_DP                     ! Dark matter particle mass in NU
    REAL(DP) :: mx_RAU(max_num_mx) != 0.5_DP*(mx_NU/me_NU)          ! Dark matter particle mass in RAU (me=0.5)

    REAL(DP), allocatable :: vmin_aux(:, :)                         ! auxiliary vmin for looping

    !REAL(DP) :: bb(6)                                              ! Dot products of reciprocal lattice basis vectors
    REAL(DP) :: bzv                                                 ! 1BZ volume

    REAL(DP) :: qnorm                                               ! Norm of q vector
    real(DP) :: q_RAU                                               ! abs of q in RAU
    !REAL(DP) :: wq                                                 ! weight of q-point for q-space integration

    INTEGER :: er_bin_type                                          ! Type of bins for integrated form factors
    INTEGER :: num_er_bins                                          ! Number of bins for integrated form factors
    REAL(DP) :: er_binsize                                          ! Bin size (Ry) for linear bins with user-provided size
    !real(DP) :: er_binsize_RAU
    integer :: iE                                                   ! iterate over Er bins

    REAL(DP) :: ermax_NU !=0.0_DP                                   ! Max recoil energy cutoff in eV
    REAL(DP) :: ermax_RAU ! = ercut_NU / Ry2eV                      ! Max recoil energy cutoff in Ry
    REAL(DP), ALLOCATABLE :: binedges(:)

    CHARACTER(2) :: numbins                                         ! String containing number of bins. It can't be larger than 99

    integer :: numqbins                                             ! Number of q bins for integrated form factor
    integer :: iq                                                   ! iterate over q bins
    real(DP) :: deltaq                                              ! q in RAU


    !REAL(DP) :: deltaE
    CHARACTER(26) :: FMT                                            ! Format specifier

    !real(DP) :: kk_SI                                                   ! norm-factor for eta
    real(DP) :: kk_kmps


    REAL(DP), ALLOCATABLE :: aux(:, :)                              ! tmp variable for output without dark matter formfactor
    REAL(DP), ALLOCATABLE :: aux_dmff(:, :)                         ! tmp variable for output with dark matter formfactor
    real(DP), allocatable :: eta_aux(:, :)                          ! tmp variable for eta


    REAL(DP), ALLOCATABLE :: ctot(:,:,:)
    REAL(DP), ALLOCATABLE :: cbinned(:,:,:,:)                       ! Binned integral. This is the final program output.
    REAL(DP), ALLOCATABLE :: formfactor(:,:)                        ! Integrated form factor. No band indices


    !REAL(DP) :: dk(3,nks,nks)                                       ! Table storing all k1-k2 values  TODO: this table is antisymmetric and can be reduced

    !INTEGER :: gsi(npwx, npwx)                                      ! G-vector sum index table for one pair of k-points

    !INTEGER :: nksf                                                 ! Number of k-points for formfactor calculation
    !integer :: ik                                                   ! Indices for k-point loops
    !integer :: ikinit


    INTEGER :: ierr                                                 ! Error index

    REAL(DP) :: tol = 1.0E-6                                        ! Tolerance for G-vector distances

    !DEBUGGING VARIABLES
    LOGICAL :: runff                                                ! For debugging purposes: set to false to skip the form factor evaluation


    character(20) :: infile = "C.dat"                               ! input file name
    character(20) :: outfile = "C.out.dat"                          ! output file name
    character(20) :: inprefix = "ge.243"                            ! prefix of input folders
    character(20) :: inrest


!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   SUBROUTINE INSTRUCTIONS
!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef __NORUN
      runff = .false.
#else
      runff = .true. ! Set to false to not skip ff calculation --for debugging purposes
#endif


    CALL start_clock( ' qedark_1k_multimass_sub ')


    print *, "           -------             "
    print *, " calculation_mode set to multimass"
    print *, " Calculating momentum-integrated, energy-dependent "
    print *, " DM-electron scattering rates for multiple DM masses."
    print *, "           -------             "


!    IF (nspin .ne. 1) THEN
!       CALL errore ('qedark_1k_multimass_sub', 'Form factor calculation works only for spin-unpolarized systems!', 1)
!    ENDIF


    !IF ( nksf > nks .or. nksf < 0 ) &
    !     CALL errore( 'qedark_1k_multimass_sub ',' nksf has a non-allowed value. Check input. ', ABS(ierr) )


    !CALL volume (tpiba, bg(:,1), bg(:,2), bg(:,3), bzv) ! could also do
    !bzv = bzvol(bg)

    !CALL qspace(bzv, .false.)

    IF (num_er_bins > 999) CALL errore( 'qedark_1k_multimass_sub ','Number of energy recoil bins cannot exceed 999 ', ABS(ierr) )

    vearth_RAU = (twobyalpha/speedoflight)*vearth_SI
    vearth_kmps = vearth_SI/1000.0
    vesc_kmps = vesc_SI/1000.0
    vesc_RAU = (twobyalpha/speedoflight)*vesc_SI
    v0_kmps = v0_SI/1000.00
    v0_RAU = (twobyalpha/speedoflight)*v0_SI
    deltav_kmps(:) = deltav_SI(:)/1000.0
    deltav_RAU(:) = (twobyalpha/speedoflight)*deltav_SI(:)
    mx_RAU(:) = 0.5_DP*(mx_NU(:)/me_NU)
    ermax_RAU = ermax_NU / Ry2eV
    !qunit = 1 / 137 * me_NU
    !er_binsize_RAU = er_binsize / Ry2eV
    kk_kmps = knorm(vesc_kmps, v0_kmps)


!!  !!!!!!!  CALL print_DM_data(vesc_kmps, vearth_kmps, v0_kmps, mx_NU, ermax_NU)


    ! Stuff for binning integral.
    ! Needs to go after ermax_RAU is calculated
    IF (num_er_bins > 0) THEN

        ALLOCATE (binedges(num_er_bins+1), STAT=ierr )
        IF( ierr /= 0 ) &
             CALL errore( 'qedark_1k_multimass_sub ',' error allocating alligk ', ABS(ierr) )

        CALL create_bins(er_bin_type, 0.0_DP, ermax_RAU, &
             num_er_bins, er_binsize, binedges)
        WRITE (*,*), " "
        WRITE (*,*), "Creating bins for C integral ..."
        WRITE (*,*), "bintype =", er_bin_type
        SELECT CASE(er_bin_type)
        CASE (1)
            WRITE (*,*), "Energy bin size=", ermax_NU/num_er_bins, "eV"
            er_binsize = ermax_NU/num_er_bins/Ry2eV
        CASE (2)
            WRITE (*,*), "Energy bin size=", Ry2eV*er_binsize, "eV"
        CASE (3)
            WRITE (*,*), "Energy bin are exponential"
        CASE default
            WRITE(*,*), "ERROR: wrong bin type selected. Creating linear bins."
        END SELECT


        ! This array contains the integrated formfactor. Second index labels month. Third array labels functional dependence of F_DM(q)
        ALLOCATE( ctot(3, nmonths, num_mx) , STAT=ierr )
        IF( ierr /= 0 ) &
             CALL errore( 'qedark_1k_multimass_sub',' cannot allocate ctot ', ABS(ierr) )
        ctot(:, :, :) = 0.0_DP


        ALLOCATE( cbinned(3, nmonths, num_er_bins+1, num_mx) , STAT=ierr )  ! +1 because last bin contains contributions with Er up to infinity
        IF( ierr /= 0 ) &
             CALL errore( 'qedark_1k_multimass_sub',' cannot allocate cbinned ', ABS(ierr) )
        cbinned(:, :, :, :) = 0.0_DP

        ! allocate arrays
        ! formfactor(iq, iE)
        allocate( formfactor(numqbins+1, num_er_bins+1) , STAT=ierr )
        if( ierr /= 0 ) &
             call errore('qedark_1k_multimass_sub',' cannot allocate formfactor ', ABS(ierr))
        formfactor(:, :) = 0.0_DP

        call C_from_file_onlyf2(infile, numqbins, num_er_bins, formfactor(:,:), ierr )

        allocate( vmin_aux(max_num_mx, num_er_bins), STAT=ierr )
        if(ierr /= 0) &
            call errore('qedark_1k_multimass_sub', 'cannot allocate vmin', ABS(ierr))
        vmin_aux(:, :) = 0.0_DP

        allocate( eta_aux(max_num_mx, num_er_bins), STAT=ierr )
        if(ierr /= 0) &
            call errore('qedark_1k_multimass_sub', 'cannot allocate eta_aux', ABS(ierr))
        eta_aux(:, :) = 0.0_DP

        allocate( aux(max_num_mx, num_er_bins), STAT=ierr )
        if(ierr /= 0) &
            call errore('qedark_1k_multimass_sub', 'cannot allocate aux', Abs(ierr))
        aux(:, :) = 0.0_DP

    ENDIF

    ! Initialize tables
    !CALL bdotb(bb)
    !CALL all_igk(nks, npwx, alligk)
    !CALL create_dk_table(nks, xk, dk)


    print *, ""
    print *, ""
    ! p = 10 RAU = 37299.2 eV
    !print *, "vmin in RAU: ", vmin(er_binsize, alphame*30, mx_RAU(2))
    !print *, "vmin in kmps: ", vmin(er_binsize, alphame*30, mx_RAU(2)) * RAU2kmps
    print *, "vmax in kmps: ", (vesc_SI + vearth_SI) / 1000_DP
    print *, "vesc in kmps: ", vesc_kmps
    print *, "vearth in kmps: ", vearth_kmps
    print *, "v0 in kmps: ", v0_kmps
    print *, "kk normalization factor in kmps: ", kk_kmps

    IF (runff) THEN
        qloop: do iq=1, numqbins
            !if (mod(iq, 100) == 1) print *, iq
            ! Make sure that G-vector exists
            !IF (alligk(ig2,ik2) < 1) CYCLE

            ! Components of q in cartesian coordinates
            !q(:) = xk(:,ik2) - xk(:,ik1) + g(:, alligk(ig2,ik2))

            ! Calculate |q| and convert to RAU multiplying by tpiba
            !qnorm = iq * deltaq * alphame
            q_RAU = iq * deltaq

            !IF (qnorm < tol) THEN
            !   ! |q|==0.0 and cannot divide by it --> leave this term out of the integral
            !   CYCLE
            !ENDIF

            ! Calculate vmin and convert to km/s
            DO imx=1, num_mx
                do iE=1, num_er_bins
                    !print *, iq*deltaq
                    vmin_aux(imx, iE) = vmin(iE*er_binsize, iq*deltaq, mx_RAU(imx) ) * RAU2kmps
                enddo
            ENDDO

            ! Loops over F_DM and month
            DO imonth=1, nmonths
                ! Current contribution to integral
                ! Weight by (eta/|q|) and integrate over k and q
                DO imx=1, num_mx
                    do iE=1, num_er_bins
                        eta_aux(imx, iE) = eta(vesc_kmps, vearth_kmps+deltav_kmps(imonth), v0_kmps, vmin_aux(imx, iE), kk_kmps)
                    end do
                    aux(imx, :) = formfactor(iq, :) * eta_aux(imx, :) * 0.25 / q_RAU
                ENDDO

                DO idmff=1, 3
                    ! Dark matter formfactor F_DM(q) ~ 1, 1/q or 1/q**2
                    ! Multiply by F_DM^2
                    IF (num_er_bins > 0) THEN
                    DO imx=1, num_mx
                        !aux_dmff(imx, :) = aux(imx, :) / (iq**( 2*(idmff-1)) )
                        cbinned(idmff, imonth, :, imx) = cbinned(idmff, imonth, :, imx) + aux(imx, :) / (q_RAU**(2 * (idmff-1)))
                        ctot(idmff, imonth, imx) = ctot(idmff, imonth, imx) + sum(aux(imx, :) / (q_RAU**(2 * (idmff-1))))
                    ENDDO
                    !if (mod(iq, 100) == 0) print *, cbinned(1, imonth, :, 3)
                    ENDIF

                    ! Add contribution to corresponding bin
                    ! Add contribution to corresponding bin
                    !DO imx=1, num_mx
                    !   cbinned(idmff, imonth, :, imx) = cbinned(idmff, imonth, :, imx) + aux_dmff(imx, :)
                    !ENDDO

                    ! Add contribution to total
                    !DO imx=1, num_mx
                    !    ctot(idmff, imonth, imx) = ctot(idmff, imonth, imx) + sum(aux_dmff(imx, :))
                    !ENDDO
                ENDDO
            enddo
        enddo qloop


        !print *, "vmin_aux: ", vmin_aux(3, :)
        !print *, "eta_aux: ", eta_aux(3, :)
        !print *, "aux: ", aux(3, :)
        !print *, "formfactor: ", formfactor(100, :)

        PRINT *, ""
        PRINT *, " Calculation finished"
        PRINT *, ""


        ! Print to files
        DO imx=1, num_mx
            WRITE (outfile, "(A2,I2.2,A4)") "C.", imx, ".dat"
            WRITE (*,*), "Writing output to ", outfile
            CALL C2file(outfile, num_er_bins, nmonths, ctot(:,:,imx), cbinned(:,:,:,imx) )
        ENDDO


    ENDIF  ! runff



    CALL stop_clock( ' qedark_1k_multimass_sub' )

end subroutine qedark_1k_multimass_sub

function muxe(mx)
    ! calculates reduced mass
    use kinds,      only: DP

    implicit none

    real(DP), intent(in) :: mx
    real(DP) :: muxe

    real(DP) :: me_NU = 0.510998910E6_DP    ! Electron mass in NU

    muxe = me_NU * mx / (me_NU + mx)
end function muxe
