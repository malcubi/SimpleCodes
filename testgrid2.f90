 
  module mytypes

  type grid
     integer, pointer :: u(:)
     integer, pointer :: v(:)
  end type grid

  end module mytypes




  module arrays

  integer, pointer :: printvar(:)
  integer, pointer :: u(:)
  integer, pointer :: v(:)

  end module arrays





  program testgrid

  use mytypes
  use arrays

  implicit none

  integer i,Nl

  integer, allocatable :: Nb(:)

  type(grid), allocatable :: box(:)

  Nl = 4

  allocate(Nb(1:Nl),box(1:Nl))

  do i=1,Nl
     Nb(i) = 2*i
     allocate(box(i)%u(1:Nb(i)),box(i)%v(1:Nb(i)))
  end do

  do i=1,Nl

     call currentgrid(box(i))

     u = +i
     v = -i

     printvar => u
     print *, pvar

     printvar => v
     print *, pvar

  end do

  end program testgrid





  subroutine currentgrid(box)

  use mytypes
  use arrays

  implicit none

  type(grid) :: box

  u => box%u
  v => box%v

  end subroutine currentgrid
