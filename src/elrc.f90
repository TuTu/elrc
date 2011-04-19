!This program calculates the long-range corrections of LJ energy.
MODULE molecule_type
  TYPE molecule
     CHARACTER(LEN=128) :: name     
     INTEGER :: num_site
     REAL(KIND=8) :: num_density
     REAL(KIND=8), ALLOCATABLE :: lj(:,:) !(2,num_site)
     LOGICAL :: is_solvent
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
     
     REAL(KIND=8) FUNCTION get_elrc(slv, slt)
       USE molecule_type
       IMPLICIT NONE
       TYPE(molecule) :: slv, slt
     END FUNCTION get_elrc
  END INTERFACE

  REAL(KIND=8), PARAMETER :: PI = 3.141592653589793238
  INTEGER, PARAMETER :: input_fileid = 10
  CHARACTER(LEN=128) :: input_filename = "input"

  !These are just for NAMELIST I/O for TYPE(molecule)
  !Because NAMELIST cannot deal with derived-type, yet
  CHARACTER(LEN=128) :: name
  INTEGER :: num_site
  REAL(KIND=8) :: num_density  
  REAL(KIND=8), ALLOCATABLE :: lj(:,:)
  LOGICAL :: is_solvent
  !------------------------------------
  
  REAL(KIND=8), ALLOCATABLE :: e_lrc(:,:) !(num_slv_type:num_slt_type)
  INTEGER :: num_slv_type, num_slt_type
  REAL(KIND=8) :: r_switch, r_cutoff
  TYPE(molecule), ALLOCATABLE :: slv(:), slt(:)

  NAMELIST /general/ num_slv_type, num_slt_type, &
       r_switch, r_cutoff
  NAMELIST /solvent/ name, num_site, lj, num_density
  NAMELIST /solute/ name, num_site, lj
  
  INTEGER :: stat, i, j


  !----- Start reading parameters from input file -----!
  open(input_fileid, FILE=input_filename, STATUS='OLD', IOSTAT=stat)
  if (stat /= 0) then
     write(*,*) 'Error: reading input file "', &
          TRIM(ADJUSTL(input_filename)), '" failed!'
  end if
  
  read(input_fileid, NML=general, IOSTAT=stat)
  if (stat /= 0) then
     write(*,*) 'Error: reading namelist "general" failed!'
  end if
  
  allocate(slv(num_slv_type))
  allocate(slt(num_slt_type))

  do i = 1, num_slv_type
     !read name and num_site
     read(input_fileid, NML=solvent, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solvent" failed!'
        write(*,*) 'Trying to read "name" and "num_site"'
     end if

     !read lj     
     allocate(lj(2, num_site))
     read(input_fileid, NML=solvent, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solvent" failed!'
        write(*,*) 'Trying to read "lj"'
     end if     

     is_solvent = .TRUE.
     call set_molecule(slv(i))
     deallocate(lj)
     call output_molecule(slv(i))
  end do

  do i = 1, num_slt_type
     !read name and num_site
     read(input_fileid, NML=solute, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solute" failed!'
        write(*,*) 'Trying to read "name" and "num_site"'
     end if

     !read lj     
     allocate(lj(2, num_site))
     read(input_fileid, NML=solute, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solute" failed!'
        write(*,*) 'Trying to read "lj"'
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
        e_lrc(j,i) = get_elrc(slv(j), slt(i))
     end do
  end do
  !----- End of e_lrc calculation -----!

  
  deallocate(slv)
  deallocate(slt)

CONTAINS
  SUBROUTINE set_molecule(mol)
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
  END SUBROUTINE output_molecule
END PROGRAM elrc

ELEMENTAL REAL(KIND=8) FUNCTION geo_average(a, b)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: a, b
  geo_average = SQRT(a*b)
END FUNCTION geo_average

REAL(KIND=8) FUNCTION get_elrc(slv, slt)
  USE molecule_type
  IMPLICIT NONE
  TYPE(molecule) :: slv, slt
  
END FUNCTION get_elrc
