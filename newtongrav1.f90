
 program newtongrav1

! *********************************************
! ***   NEWTONIAN GRAVITATIONAL POTENTIAL   ***
! *********************************************

! This program solves the Poisson equation for a Newtonian
! gravitational potential in spherical symmetry, given a
! distribution of mass.
!
! The equation to solve is (in units such that G=1):
!
!  2
! d phi + (2/r) d phi = 4 pi rho
!  r             r
!
! with rho the mass density.

! Declare variables.

  implicit none

  integer i,j          ! Counters.
  integer Nr           ! Total number of grid points.

  real(8) dr           ! Grid spacing.
  real(8) mass         ! Total integrated mass.
  real(8) pi,aux

  real(8), allocatable, dimension(:) :: r       ! Radial coordinate.
  real(8), allocatable, dimension(:) :: phi     ! Gravitational potential.
  real(8), allocatable, dimension(:) :: psi     ! dphi/dx.
  real(8), allocatable, dimension(:) :: rho     ! Mass density.
  real(8), allocatable, dimension(:) :: phiasym ! Asymptotic phi.
  real(8), allocatable, dimension(:) :: res     ! Residual.

  character(20) density ! Density profile.


! *******************
! ***   NUMBERS   ***
! *******************

  pi = acos(-1.0d0)


! **************************
! ***   GET PARAMETERS   ***
! **************************

! Grid spacing.

  print *
  print *, 'Give grid spacing dr'
  read(*,*) dr

! Number of grid points.

  print *
  print *, 'Give total number of grid points Nr'
  read(*,*) Nr

! Density profile.

  print *
  print *, 'Give the mass density profile (constant,gaussian)'
  read(*,*) density

  print *


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

! Allocate arrays.

  allocate(r(0:Nr))
  allocate(phi(0:Nr),psi(0:Nr),rho(0:Nr),phiasym(0:Nr),res(0:Nr))


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do i=0,Nr
     r(i) = dble(i)*dr
  end do


! ************************
! ***   MASS DENSITY   ***
! ************************

! Give mass density.

  if (density=="constant") then

!    Constant density star.

     do i=0,Nr

        if (r(i)<=1.d0) then
           rho(i) = 1.d0
        else
           rho(i) = 0.d0
        end if

     end do

  else if (density=="gaussian") then

!    Gaussian density star.

     rho = exp(-r**2)

  else

!    Unknown density profile.

     print *, 'Unknown density profile, aborting ...'
     print *
     stop

  end if


! ***************************
! ***   FIND TOTAL MASS   ***
! ***************************

! Integrate the density to find total mass
! using trapecium rule.

  mass = 0.d0

  do i=1,Nr

     if ((density/='constant').or. &
        ((density=='constant').and.(r(i)<=1.d0))) then
        mass = mass + 2.d0*pi*dr &
             *(rho(i-1)*r(i-1)**2 + rho(i)*r(i)**2)
     end if

  end do

! Print value of total mass to screen.

  print *, 'Total integrated mass =',mass
  print *

! Find asymtptotic phi.  We only do it from r=1 out
! to avoid having very large values close to r=0.

  do i=0,Nr
     if (r(i)>=1.d0) then
        phiasym(i) = - mass/r(i)
     else
        phiasym(i) = 0.d0
     end if
  end do


! ****************************************
! ***   FIND GRAVITATIONAL POTENTIAL   ***
! ****************************************

! Call subroutine for solving Poisson equation.

  call poisson(Nr,dr,r,phi,psi,rho)


! *************************************************************
! ***   ADD CONSTANT TO MATCH EXTERIOR BOUNDARY CONDITION   ***
! *************************************************************

! We now need to add a constant to the solution to make
! sure that far away the potential behaves as phi ~ 1/r.
!
! Notice that this condition implies psi = -1/r**2 = phi/r,
! so we use this to determine the constant.

  phi = phi - (phi(Nr) + psi(Nr)*r(Nr))


! ******************************
! ***   CALCULATE RESIDUAL   ***
! ******************************

  res = 0.d0

  do i=1,Nr-1
     res(i) = (phi(i+1) - 2.d0*phi(i) + phi(i-1))/dr**2 &
            + 1.d0/r(i)*(phi(i+1) - phi(i-1))/dr &
            - 4.d0*pi*rho(i)
  end do


! *********************
! ***   SAVE DATA   ***
! *********************

! Newtonian potential.

  open(1,file='newton_phi.rl',form='formatted',status='replace')

  do i=0,Nr
     if (dabs(phi(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") r(i),phi(i)
     else
        write(1,"(2ES16.8)") r(i),0.0D0
     end if
  end do

  close(1)

! Derivative of potential.

  open(1,file='newton_psi.rl',form='formatted',status='replace')

  do i=0,Nr
     if (dabs(psi(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") r(i),psi(i)
     else
        write(1,"(2ES16.8)") r(i),0.0D0
     end if
  end do

  close(1)

! Density

  open(1,file='newton_rho.rl',form='formatted',status='replace')

  do i=0,Nr
     if (dabs(rho(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") r(i),rho(i)
     else
        write(1,"(2ES16.8)") r(i),0.0D0
     end if
  end do

  close(1)

! Asymptotic phi.

  open(1,file='newton_phiasym.rl',form='formatted',status='replace')

  do i=0,Nr
     if (dabs(phiasym(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") r(i),phiasym(i)
     else
        write(1,"(2ES16.8)") r(i),0.0D0
     end if
  end do

  close(1)

! Residual.

  open(1,file='newton_res.rl',form='formatted',status='replace')

  do i=0,Nr
     if (dabs(res(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") r(i),res(i)
     else
        write(1,"(2ES16.8)") r(i),0.0D0
     end if
  end do

  close(1)


! ***************
! ***   END   ***
! ***************

  end program newtongrav1







  subroutine poisson(Nr,dr,r,phi,psi,rho)

! Here we solve the Poisson equation by writing
! it in first order form and using Runge-Kutta.
!
! In order to rewrite it in first order form we
! introduce the auxiliary variable:  psi = dphi/dr.

  implicit none

  integer i,Nr

  real(8) dr
  real(8) sphi,spsi
  real(8) phih,psih
  real(8) pi

  real(8), dimension(0:Nr) :: r
  real(8), dimension(0:Nr) :: phi
  real(8), dimension(0:Nr) :: psi
  real(8), dimension(0:Nr) :: rho


! *******************
! ***   NUMBERS   ***
! *******************

  pi = acos(-1.d0)


! *****************************
! ***   INITIALIZE ARRAYS   ***
! *****************************

  phi = 0.0D0
  psi = 0.0D0


! **************************
! ***   SOLVE EQUATION   ***
! **************************

! For the moment we use only second order Runge-Kutta.

! Give data at the origin.

  phi(0) = 0.d0
  psi(0) = 0.d0

! Note: The first point must be treated differently
! since we have a division by zero at the origin.
! If we take  psi~a*r  for r<<1, with "a" constant,
! one can show that a=(4/3)*pi*rho(r=0).
! (But notice that this is only second order)

  do i=1,Nr

!    Calculate sources at left point.

     sphi = psi(i-1)

     if (i==1) then
        spsi = (4.d0/3.d0)*pi*rho(0)
     else
        spsi = - 2.d0*psi(i-1)/r(i-1) + 4.d0*pi*rho(i-1)
     end if

!    Advance half a step.

     phih = phi(i-1) + 0.5d0*dr*sphi
     psih = psi(i-1) + 0.5d0*dr*spsi

!    Calculate sources at intermediate point.

     sphi = psih
     spsi = - 4.d0*psih/(r(i-1) + r(i)) &
          + 2.d0*pi*(rho(i-1) + rho(i))

!    Advanced full step.

     phi(i) = phi(i-1) + dr*sphi
     psi(i) = psi(i-1) + dr*spsi

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine poisson

