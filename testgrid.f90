 
! This simple program tests a type definition
! for having several grids with different sizes
! for a set of arrays that have the same name
! on all grids.



! This module defines the type 'gridfuncs'
! with all the arrays the code will use.
! They have to be define as pointers since
! since they will be the target of other
! pointers and the 'target' attribute is
! not allowed on  a type definition.

  module mytypes

  type gridfuncs
     integer, pointer :: u(:)
     integer, pointer :: v(:)
  end type gridfuncs

  end module mytypes



! This module again contains all the
! arrays defined above as pointers.
! They have the same names on purpose,
! but in principle are different things.
! The ones here are generic pointers that
! will later point to the corresponding
! elements of a given grid.

! It also has a pointer array for
! printing to screen.

  module arrays

  integer, pointer :: printvar(:)

  integer, pointer :: u(:)
  integer, pointer :: v(:)

  end module arrays





! Main program.

  program testgrid

! Load modules.

  use mytypes
  use arrays

! Extra varibles.

  implicit none

  integer i,j  ! Counters
  integer Nl   ! Number of different grid levels

  integer, allocatable :: N(:)  ! Number of ints on each grid level.

! Define an allocatable array of grids of type 'gridfuncs'.

  type(gridfuncs), allocatable :: grid(:)

! Number of grid levels.

  Nl = 4

! Allocate N and grid.

  allocate(N(1:Nl),grid(1:Nl))

! Now allocate the different grid arrays
! for each grid level.

  do i=1,Nl
     N(i) = 2*i
     allocate(grid(i)%u(1:N(i)),grid(i)%v(1:N(i)))
  end do

! Give values to the different grid levels
! and print them to screen.

  do i=1,Nl

     call currentgrid(grid(i))

     do j=1,N(i)
       u(j) = 10*i + j
       v(j) = 100*i + j
     end do

     printvar => u
     print *, printvar

     printvar => v
     print *, printvar

  end do

  end program testgrid





! Subroutine for pointing the generic
! pointer arrays to the corresponding
! elements of a given grid.

  subroutine currentgrid(grid)

  use mytypes
  use arrays

  implicit none

  type(gridfuncs) :: grid

  u => grid%u
  v => grid%v

  end subroutine currentgrid
