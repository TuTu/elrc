!This module provides functions to calculate the long-range corrections of LJ energies.
MODULE molecule_type
  TYPE molecule
     CHARACTER(LEN=128) :: name     
     INTEGER :: num_site
     REAL(KIND=8) :: num_density
     REAL(KIND=8), ALLOCATABLE :: lj(:,:) !(sigma:epsilon,num_site)
     LOGICAL :: is_solvent !only solvent has num_density, so far.
  END TYPE molecule
END MODULE molecule_type

PROGRAM elrc
  USE molecule_type
  IMPLICIT NONE
  INTERFACE
     ELEMENTAL REAL(KIND=8) FUNCTION geo_average(a, b)
       IMPLICIT NONE       
       REAL(KIND=8), INTENT(IN) :: a, b
     END FUNCTION geo_average
  END INTERFACE
  REAL(KIND=8), PARAMETER :: PI = 3.141592653589793238  
  REAL(KIND=8), PARAMETER :: JoulePerCal = 4.184
!  INTEGER, PARAMETER :: input_fileid = 10
  INTEGER, PARAMETER :: input_fileid = 5 !use stdin as input
  CHARACTER(LEN=128) :: input_filename = "input"

  !These variables are just for NAMELIST I/O for TYPE(molecule)
  !Because NAMELIST cannot deal with derived-type, yet
  CHARACTER(LEN=128) :: name
  INTEGER :: num_site
  REAL(KIND=8) :: num_density  
  REAL(KIND=8), ALLOCATABLE :: lj(:,:)
  LOGICAL :: is_solvent
  !------------------------------------
  
  REAL(KIND=8), ALLOCATABLE :: e_lrc(:,:) !This is our GOAL 
                              !e_lrc(num_slv_type:num_slt_type)
  INTEGER :: num_slv_type, num_slt_type
  REAL(KIND=8) :: r_switch, r_cutoff
  TYPE(molecule), ALLOCATABLE :: slv(:), slt(:)

  NAMELIST /general/ num_slv_type, num_slt_type, &
       r_switch, r_cutoff
  NAMELIST /solvent/ name, num_site, lj, num_density
  NAMELIST /solute/ name, num_site, lj
  
  INTEGER :: stat, i, j


  !----- Start reading parameters from input file -----!
!  open(input_fileid, FILE=input_filename, STATUS='OLD', IOSTAT=stat)
  open(input_fileid, IOSTAT=stat)
  if (stat /= 0) then
     write(*,*) 'Error: reading input file "', &
          TRIM(ADJUSTL(input_filename)), '" failed!'
     call EXIT(1)
  end if

  !read "general" parameters
  read(input_fileid, NML=general, IOSTAT=stat)
  if (stat /= 0) then
     write(*,*) 'Error: reading namelist "general" failed!'
     call EXIT(1)
  end if

  if (r_cutoff < r_switch) then
     write(*,*) "Error: r_cutoff < r_switch !"
     call EXIT(1)
  end if

  write(*,*) "num_slv_type=", num_slv_type
  write(*,*) "num_slt_type=", num_slt_type
  write(*,*) "r_switch=", r_switch
  write(*,*) "r_cutoff=", r_cutoff
  write(*,*)
  
  allocate(slv(num_slv_type))
  allocate(slt(num_slt_type))
  allocate(e_lrc(num_slv_type, num_slt_type))

  !read "solvent"
  do i = 1, num_slv_type
     !read name and num_site
     read(input_fileid, NML=solvent, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solvent" failed!'
        write(*,*) 'Trying to read "name" and "num_site"'
        call EXIT(1)
     end if

     !read lj     
     allocate(lj(2, num_site))
     read(input_fileid, NML=solvent, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solvent" failed!'
        write(*,*) 'Trying to read "lj"'
        call EXIT(1)
     end if     

     is_solvent = .TRUE.
     call set_molecule(slv(i))
     deallocate(lj)
     call output_molecule(slv(i))
  end do

  !read "solute"
  do i = 1, num_slt_type
     !read name and num_site
     read(input_fileid, NML=solute, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solute" failed!'
        write(*,*) 'Trying to read "name" and "num_site"'
        call EXIT(1)
     end if

     !read lj     
     allocate(lj(2, num_site))
     read(input_fileid, NML=solute, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solute" failed!'
        write(*,*) 'Trying to read "lj"'
        call EXIT(1)        
     end if     

     is_solvent = .FALSE.
     call set_molecule(slt(i))
     deallocate(lj)
     call output_molecule(slt(i))
  end do
  !----- End of reading input file -----!


  !----- Start e_lrc calculation -----!
  e_lrc = 0.
  do i = 1, num_slt_type
     do j = 1, num_slv_type
        call get_elrc(e_lrc(j,i), slv(j), slt(i))
     end do
  end do
  !----- End of e_lrc calculation -----!

  
  !----- Output results -----!
  write(*,*) "Unit: kJ/mol"
  write(*,"(10X)", ADVANCE='NO')
  do i = 1, num_slv_type
    if (i == num_slv_type) then
      write(*,"(1X,A15)") TRIM(ADJUSTL(slv(i)%name))
    else
      write(*,"(1X,A15)", ADVANCE='NO') TRIM(ADJUSTL(slv(i)%name))
    end if
  end do

  do i = 1, num_slt_type
     write(*,"(A10)", ADVANCE='NO') TRIM(ADJUSTL(slt(i)%name))
     do j = 1, num_slv_type
        if (j == num_slv_type) then
           write(*,"(1X,F15.6)") e_lrc(j, i)
        else
           write(*,"(1X,F15.6)", ADVANCE='NO') e_lrc(j, i)
        end if
     end do
  end do
  
  write(*,*)

  write(*,*) "Unit: kCal/mol"
  write(*,"(10X)", ADVANCE='NO')
  do i = 1, num_slv_type
    if (i == num_slv_type) then
      write(*,"(1X,A15)") TRIM(ADJUSTL(slv(i)%name))
    else
      write(*,"(1X,A15)", ADVANCE='NO') TRIM(ADJUSTL(slv(i)%name))
    end if
  end do

  do i = 1, num_slt_type
     write(*,"(A10)", ADVANCE='NO') TRIM(ADJUSTL(slt(i)%name))
     do j = 1, num_slv_type
        if (j == num_slv_type) then
           write(*,"(1X,F15.6)") e_lrc(j, i)/JoulePerCal
        else
           write(*,"(1X,F15.6)", ADVANCE='NO') e_lrc(j, i)/JoulePerCal
        end if
     end do
  end do
  !----- End of output -----!
  
  deallocate(slv)
  deallocate(slt)
  deallocate(e_lrc)

CONTAINS
  SUBROUTINE set_molecule(mol)
    !set molecule info from NAMELIST variables
    IMPLICIT NONE
    TYPE(molecule), INTENT(OUT) :: mol
    mol%name = name
    mol%num_site = num_site
    allocate(mol%lj(2, num_site))
    mol%lj = lj
    if (is_solvent) then
       mol%num_density = num_density
    else
       !so far, solute has no num_density
       mol%num_density = -1.
    end if
  END SUBROUTINE set_molecule

  SUBROUTINE output_molecule(mol)     
    IMPLICIT NONE
    TYPE(molecule), INTENT(INOUT) :: mol
    INTEGER :: i
    write(*,*) TRIM(ADJUSTL(mol%name)), mol%num_site
    do i = 1, num_site
       write(*,*) mol%lj(:,i)
    end do
    write(*,*)
  END SUBROUTINE output_molecule

  SUBROUTINE get_elrc(e_lrc, slv, slt)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: e_lrc
    TYPE(molecule), INTENT(IN) :: slv, slt
    REAL(KIND=8) :: eps, sig, const0, sig_12, sig_6
    REAL(KIND=8) :: term1, term2, term3, term4
    REAL(KIND=8) :: term3_1, term3_2, term3_3, term3_4
    REAL(KIND=8) :: term4_1, term4_2, term4_3, term4_4    
    REAL(KIND=8) :: const3_0, const3_2, const3_3, const3_4

    INTEGER :: i, j

    const0 = 16 * PI * slv%num_density

    if (.NOT. r_switch == r_cutoff) then
       const3_0 = (r_cutoff**2 - r_switch**2)**3
       const3_2 = 3 * r_switch**2 + 3 * r_cutoff**2
       const3_3 = 6 * r_cutoff**2 * r_switch**2
       const3_4 = r_cutoff**6 - 3 * r_cutoff**4 * r_switch**2
    end if
    e_lrc = 0.
    do i = 1, slv%num_site
       do j = 1, slt%num_site
          sig = geo_average(slv%lj(1,i), slt%lj(1,j))    
          eps = geo_average(slv%lj(2,i), slt%lj(2,j))
          sig_6 = sig**6
          sig_12 = sig_6 * sig_6

          term1 = sig_12 / (9 * r_switch**9)
          term2 = -sig_6 / (3 * r_switch**3)
          
          if (r_switch == r_cutoff) then
             term3 = 0.
             term4 = 0.
          else
             term3_1 = -2.0d0/3.0d0 * (1/r_cutoff**3 - 1/r_switch**3)
             term3_2 = const3_2 / 5 * (1/r_cutoff**5 - 1/r_switch**5)
             term3_3 = -const3_3 / 7 * (1/r_cutoff**7 - 1/r_switch**7)
             term3_4 = -const3_4 / 9 * (1/r_cutoff**9 - 1/r_switch**9)
             term3 = -sig_12 / const3_0 * &
                  &(term3_1 + term3_2 + term3_3 + term3_4)
             
             term4_1 = -2.0d0/3.0d0 * (r_cutoff**3 - r_switch**3)
             term4_2 = const3_2 * (r_cutoff - r_switch)
             term4_3 = const3_3 * (1/r_cutoff - 1/r_switch)
             term4_4 = const3_4 / 3 * (1/r_cutoff**3 - 1/r_switch**3)
             term4 = -sig_6 / const3_0 * &
                  &(term4_1 + term4_2 + term4_3 + term4_4)
          end if
          
          e_lrc = e_lrc + const0 * eps * (term1 + term2 + term3 + term4)
       end do
    end do
  END SUBROUTINE get_elrc
END PROGRAM elrc

ELEMENTAL REAL(KIND=8) FUNCTION geo_average(a, b)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: a, b
  geo_average = SQRT(a*b)
END FUNCTION geo_average
