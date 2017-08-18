subroutine qedark_1k_f2_sum(&
        nksf, &
        num_er_bins, &
        numqbins)
    ! read form factor from different files and sums it up.
    ! the files have to be stored in folder $input_prefix_$kvalue/

    ! import section
    use kinds,          only: DP

    ! variable defintion section
    implicit none

    integer :: ik                                           ! Indices for k-point loops
    integer :: ikinit
    integer :: nksf                                         ! Number of k-points for formfactor calculation

    integer :: num_er_bins                                  ! Number of Energy bins for integrated form factors
    integer :: numqbins                                     ! Number of q bins for integrated form factor
    REAL(DP), ALLOCATABLE :: ctot(:,:)                      ! Integrated form factor. No band indices
    real(DP), allocatable :: ctmp(:,:)                      ! tmp file for integrated form factor

    character(20) :: infile = "K.dat"                       ! input file name
    character(20) :: outfile = "C.dat"                      ! output file name
    character(20) :: inprefix = "ge"                        ! prefix of input folders
    character(20) :: inrest

    integer :: ierror                                       ! Error index
    character(100) :: ierror_message = ""                   ! error message

    ! code
    ! allocate arrays
    allocate( ctot(numqbins+1, num_Er_bins+1) , STAT=ierror )
    if( ierror /= 0 ) &
         call errore('qedark_f2_1k_sum',' cannot allocate ctot ', ABS(ierror))
    ctot(:,:) = 0.0_DP

    allocate( ctmp(numqbins+1, num_Er_bins+1) , STAT=ierror )
    if( ierror /= 0 ) &
         call errore('qedark_f2_1k_sum',' cannot allocate ctmp ', ABS(ierror))
    ctmp(:,:) = 0.0_DP

    ! read different files
    ikinit = 1
    readloop: do ik=ikinit, nksf
        write(inrest, "(A1,I3.3,A3,I3.3,A4)") "_", ik, "/K.", ik, ".dat"
        infile = trim(inprefix)//trim(inrest)
        call C_from_file_onlyf2(infile, numqbins, num_er_bins, ctmp(:,:), ierror )
        if (ierror == 10) then
            call errore('qedark_f2_1k_sum', 'number of bins must be the same as in file', abs(ierror))
        elseif (ierror /= 0) then
            write (ierror_message, *) ' error while reading file: ', infile
            call errore ('qedark_f2_1k_sum', ierror_message, abs(ierror))
        endif
        ctot(:,:) = ctot(:,:) + ctmp(:,:)
    enddo readloop

    ! save output
    call C2file_onlyf2(outfile, numqbins, num_er_bins, ctot(:,:))

end subroutine qedark_1k_f2_sum
