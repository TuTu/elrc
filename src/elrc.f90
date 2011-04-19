PROGRAM elrc
  IMPLICIT NONE
  TYPE molecule
     CHARACTER(LEN=128) :: name     
     INTEGER :: num_site
     REAL(KIND=8), ALLOCATABLE :: lj(:,:)
  END TYPE molecule
  
  INTEGER, PARAMETER :: input_fileid = 10
  CHARACTER(LEN=128) :: input_filename = "input"

  !These are just for NAMELIST I/O for TYPE(molecule)
  !Because NAMELIST cannot deal with derived-type, yet
  CHARACTER(LEN=128) :: name
  INTEGER :: num_site
  REAL(KIND=8), ALLOCATABLE :: lj(:,:)
  !------------------------------------
  
  INTEGER :: num_slv_type, num_slt_type
  TYPE(molecule), ALLOCATABLE :: slv(:), slt(:)

  NAMELIST /general/ num_slv_type, num_slt_type
  NAMELIST /solvent/ name, num_site, lj
  NAMELIST /solute/ name, num_site, lj
  
  INTEGER :: stat, i, j

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
     !read slv(i)'s name and num_site
     read(input_fileid, NML=solvent, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solvent" failed!'
        write(*,*) 'Trying to read "name" and "num_site"'
     end if

     deallocate(lj)
     allocate(lj(num_site))
     
     !read slv(i)'s lj
     read(input_fileid, NML=solvent, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solvent" failed!'
        write(*,*) 'Trying to read "name" and "num_site"'
     end if     

     allocate(slv(i)%lj(2, slv(i)%num_site))
     read(input_fileid, NML=solvent, IOSTAT=stat)
     if (stat /= 0) then
        write(*,*) 'Error: reading namelist "solvent" failed!'
        write(*,*) 'Trying to read "lj"'               
     end if
  end do
  
  deallocate(slv)
  deallocate(slt)

CONTAINS
  SUBROUTINE set_molecule(mol)
    IMPLICIT NONE
    TYPE(molecule), INTENT(OUT) :: mol
    mol%name = name
    mol%num_site = num_site
    allocate(mol%lj(num_site))
    mol%lj = lj
  END SUBROUTINE set_molecule
END PROGRAM elrc
