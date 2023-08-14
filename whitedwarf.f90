
 program whitedwarf

! *****************************
! ***   WHITE DWARF STARS   ***
! *****************************

! This program solves the structure equation for a white
! dwarf supported by electron degeneracy pressure.
!
! The equation to solve is a generalized Lane-Emden
! equation of the form:
!
!             
!  d    /     xi^2  theta       d theta  \
! ––––  | ––––––––––––––––––––  –––––––  |  + xi^2 theta^3 = 0
! d xi  \ (theta^2 + 1/k^2)^(1/2) d xi     /
!
! with k a constant defined below.
!
! In the above equation "xi" is a radial coordinate related
! to the physical radius "r" as  r = a xi, with:
!
! a = (3 pi)^(1/2) (Le  MP / 4 k mH)
!
! with:
!
! Le = hbar/ me c       Compton wavelength of the electron (me = electron mass).
! MP = sqrt(hbar c/G)   Planck's mass
! mH =                  mass of the hydrogen atom.
!
! On the other hand, "theta" is a variable related to the
! density "rho" as:
!
! rho = rho0 theta^3
!
! with "rho" the central value of the density. Finally, "k" is a
! constant related to "rho0" as:
!
! k^3 = ( 3 pi^2 Le^3 / 2 mH ) rho0

! Declare variables.

  implicit none

  integer i           ! Counter.
  integer N           ! Total number of grid points.
  integer type        ! Type of run.

  real(8) MP          ! Plack mass.
  real(8) mH          ! Mass of hydrogen atom.
  real(8) me          ! Electron mass.
  real(8) MS          ! Solar mass.
  real(8) RE          ! Earth radius.
  real(8) c           ! speed of light.
  real(8) hbar        ! Planck constant (reduced).
  real(8) Le          ! Compton wavelength of electron.

  real(8) k,rho0      ! Central densities.
  real(8) Mstar       ! Total mass of star.
  real(8) Rstar       ! Radius of star.           
  real(8) a           ! Scale factor.
  real(8) dxi,dr      ! Grid spacings.
  real(8) ximax       ! Boundary of star in xi coordinate.
  real(8) pi

  real(8), allocatable, dimension(:) :: r,xi        ! Radial coordinates.
  real(8), allocatable, dimension(:) :: rho,theta   ! Density variables.


! *******************
! ***   NUMBERS   ***
! *******************

  pi = acos(-1.0d0)


! *********************
! ***   CONSTANTS   ***
! *********************

! Planck mass.

  MP = 2.18D-8

! Mass of hydrogen atom.

  mH = 1.67D-27

! Mass of electron.

  me = 9.1D-31

! Speed of light.

  c = 2.998D+8

! Plack constant.

  hbar = 6.626D-34/(2.d0*pi)

! Solar mass

  MS = 1.99D+30

! Earth radius.

  RE = 6.37D+6

! Compton wavelength of electron.

  Le = hbar/(me*c)


! **************************
! ***   GET PARAMETERS   ***
! **************************

! Grid spacing.

  print *
  print *, 'Give grid spacing dxi'
  read(*,*) dxi

! Number of grid points.

  print *
  print *, 'Give total number of grid points N'
  read(*,*) N

! Single density plot or full mass plot?

  print *
  print *, 'Choose type of run:'
  print *, '  1  Single density plot.'
  print *, '  2  Mass versus radius plot.'
  read(*,*) type

! For single density plot give value of parameter k.'

  if (type==1) then
     print *
     print *, 'Give value of parameter k'
     read(*,*) k
  end if

  print *


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

! Allocate arrays.

  allocate(r(0:N),xi(0:N))
  allocate(rho(0:N),theta(0:N))


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do i=0,N
     xi(i) = dble(i)*dxi
  end do


! *******************************
! ***   SINGLE DENSITY PLOT   ***
! *******************************

  if (type==1) then


!    ******************************
!    ***   INTEGRATE EQUATION   ***
!    ******************************

     call integrate(N,k,dxi,ximax,xi,theta)


!    *******************************
!    ***   PHYSICAL QUANTITIES   ***
!    *******************************

!    Scale factor.

     a = dsqrt(3.d0*pi)*Le*MP/(4.d0*k*mH)

!    Find physical radial coordinate r and increment dr.

     r  = a*xi
     dr = a*dxi

!    Find physical density rho.

     rho0 = 2.d0*mH*k**3/(3.d0*pi**2*Le**3)
     rho  = rho0*theta**3
     print *, 'Central density (SI units kg/m^3): ',rho0

!    Find radius of star.

     Rstar = a*ximax
     print *, 'Star radius (meters,Earth radius): ',Rstar,Rstar/RE

!    Integrate mass of star.

     Mstar = 0.d0

     do i=1,N
        Mstar = Mstar + 2.d0*pi*dr &
              *(rho(i-1)*r(i-1)**2 + rho(i)*r(i)**2)
     end do

     print *, 'Star mass   (kg,Solar masses):     ',Mstar,Mstar/MS
     print *


!    *********************
!    ***   SAVE DATA   ***
!    *********************

!    Density variable theta in terms of xi.

     open(1,file='whitedwarf_xi.rl',form='formatted',status='replace')

     do i=0,N
        if (dabs(theta(i))>1.d-50) then
           write(1,"(2ES16.8)") xi(i),theta(i)
        else
           write(1,"(2ES16.8)") xi(i),0.0D0
        end if
     end do

     close(1)

!    Density rho in terms of r.

     open(1,file='whitedwarf_rho.rl',form='formatted',status='replace')

     do i=0,N
        if (dabs(rho(i))>1.d-50) then
           write(1,"(2ES16.8)") r(i)/RE,rho(i)
        else
           write(1,"(2ES16.8)") r(i)/RE,0.0D0
        end if
     end do

     close(1)


! ***********************************
! ***   MASS VERSUS RADIUS PLOT   ***
! ***********************************

  else


!    **************************
!    ***   OPEN DATA FILE   ***
!    **************************

     open(1,file='whitedwarf_MvR.rl',form='formatted',status='replace')


!    *********************************
!    ***   LOOP OVER VALUES OF k   ***
!    *********************************

     k = 0.1d0

     do while (k<=100.0)

!       ******************************
!       ***   INTEGRATE EQUATION   ***
!       ******************************

        call integrate(N,k,dxi,ximax,xi,theta)


!       *******************************
!       ***   PHYSICAL QUANTITIES   ***
!       *******************************

!       Scale factor.

        a = dsqrt(3.d0*pi)*Le*MP/(4.d0*k*mH)

!       Find physical radial coordinate r and increment dr.

        r  = a*xi
        dr = a*dxi

!       Find physical density rho.

        rho0 = 2.d0*mH*k**3/(3.d0*pi**2*Le**3)
        rho  = rho0*theta**3
        !print *, 'Central density (SI units kg/m^3): ',rho0

!       Find radius of star.

        Rstar = a*ximax
        !print *, 'Star radius (meters,Earth radius): ',Rstar,Rstar/RE

!       Integrate mass of star.

        Mstar = 0.d0

        do i=1,N
           Mstar = Mstar + 2.d0*pi*dr &
                 *(rho(i-1)*r(i-1)**2 + rho(i)*r(i)**2)
        end do

        !print *, 'Star mass   (kg,Solar masses):     ',Mstar,Mstar/MS

!       Save data to file.

        write(1,"(2ES16.8)") Mstar/MS,Rstar/RE


!       ***********************
!       ***   INCREMENT k   ***
!       ***********************

        k = k + 0.01d0

     end do


!    ***************************
!    ***   CLOSE DATA FILE   ***
!    ***************************

     close(1)

  end if


! ***************
! ***   END   ***
! ***************

  end program whitedwarf










  subroutine integrate(N,k,dxi,ximax,xi,theta)

! Here we integrate the structure equation in
! first order form. We first define:
!
! u = theta^2
!
! F = xi^2 / (u + 1/k^2)^(1/2)  du/dxi
!
! The equation to solve then becomes the following
! system:
!
! du/dxi  =  [ (u + 1/k^2)^(1/2) / xi^2 ] F
!
! dF/dxi  =  - 2 xi^2 u^(3/2)
!
! This is integrated from the origin out taking:
!
! u(0) = 1,    F(0) = 0
!
! We integrate until we find u=0, which marks the surface
! of the star.

  implicit none

  integer i,N

  real(8) k,dxi,ximax
  real(8) su,sF,uh,Fh,xih

  real(8), dimension(0:N) :: xi,theta
  real(8), dimension(0:N) :: u,F


! *****************************
! ***   INITIALIZE ARRAYS   ***
! *****************************

  theta = 0.d0

  u = 0.d0
  F = 0.d0


! **************************
! ***   SOLVE EQUATION   ***
! **************************

! For the moment we use only second order Runge-Kutta,
! might improve this later.

! Give data at the origin.

  u(0) = 1.d0
  F(0) = 0.d0

! Note: The first point must be treated differently
! since we have a division by zero at the origin.
!
! We first assume F = F0*xi^3. Substituting in the
! equation fof dF/dxi we find:  F0 = - (2/3).
!

! Next we assume u = 1 + u0 xi^2. Substituting in the
! equation for du/dxi we now find:  u0 = - (1 + k^2)^(1/2) / 3

  u(1) = 1.d0 - (dsqrt(1.d0 + 1.d0/k**2)/3.d0)*dxi**2
  F(1) = - 2.d0/3.d0*dxi**3

! Loop over grid points.

  do i=1,N

!    Calculate sources at left point.

     su = dsqrt(u(i) + 1.d0/k**2)/xi(i)**2*F(i)
     sF = - 2.d0*xi(i)**2*u(i)**1.5d0

!    Advance half a step.

     xih = xi(i) + 0.5d0*dxi

     uh = u(i) + 0.5d0*dxi*su
     Fh = F(i) + 0.5d0*dxi*sF

!    Check if uh is negative. In that
!    case we reached the surface.
 
     if (uh<=0.d0) then
        u(i+1) = 0.d0
        goto 100
     end if

!    Calculate sources at intermediate point.

     su = dsqrt(uh + 1.d0/k**2)/xih**2*Fh
     sF = - 2.d0*xih**2*uh**1.5d0

!    Advance full step.

     u(i+1) = u(i) + dxi*su
     F(i+1) = F(i) + dxi*sF

!    Check if u(i+1) is negative. In that
!    case we reached the surface.
 
     if (u(i+1)<=0.d0) then
        u(i+1) = 0.d0
        goto 100
     end if

  end do

! Find theta from u.

  100 continue

  ximax = xi(i)

  theta = sqrt(u)


! ***************
! ***   END   ***
! ***************

  end subroutine integrate


