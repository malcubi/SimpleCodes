program trumpet

! ********************************
! ***   TRUMPET INITIAL DATA   ***
! ********************************

! This program solves for "trumpet" black hole
! initial data.
!
! Then original equation for the conformal factor
! has the form:
!
!  2                          8          7  6
! d psi  +  2 d psi / 2  +  k1  / ( 4 psi  r  )  =  0
!  r           r
!
! with k1 some constant. Notice that this equation has
! an exact solution of the form:
!
! psi = k1/sqrt(r)
!
! Unfortunately, this solution does not satisfy the
! boundary asymptotic boundary, which requires that
! psi goes to 1 far away.
!
! In order to solve the equation with the correct asymptotic
! condition, we first define:
!
! u  :=  sqrt(r) psi
!
! The equation to solve now becomes:
!
!  2                        2      8   7
! d u  +  d u / r  +  ( 1/4r ) [ k1 / u  -  u ]  =  0
!  r       r
!
! Notice that u=k1 solves the equation.
!
! One could now try to solve this equation integrating
! outwards from r=0, and using the boundarty conditions:
! u(0) = k1, u'(0) = 0. Unfortunately, one then only recovers
! the exact solution above, which does not satisfy the asymptotic
! condition.
!
! One can instead ask that, for small r, u behaves as:
!
! u ~ k1 + k2 r^n
!
! When we substitute this in the equation we find that the
! constant k2 is arbitrary, and the power "n" must be equal
! to sqrt(2).
!
! This has the disadvantage that for small r the second derivative
! of u diverges (u is only C1 at r=0). But we can use the expansion
! above to start the integration at some small value of r.  Here
! we in fact stagger the origin, and set initial data at r=dr/2.
!
! The value of the constant k2 can then be used in a shooting method
! to guarantee that psi goes to 1 at infinity.
!
! Notice that once we fix the value of K2 with the shooting,
! we will also get a specific value for the mass M of the
! black hole:
!
! psi  ~  1 + M/2r
!
! The value of the mass will depend on the value of the constant k1.
! So if we want to solve for a specific value of the mass, we will
! need a double shooting: One to find k2 for a given value of k1,
! and a second one to find k1 for a given value of M.

! Declare variables.

  implicit none

  integer i            ! Counters.
  integer Nr           ! Total number of grid points.
  integer maxiter      ! Maximum number of iterations.

  real(8) dr,hdr       ! Grid spacing.
  real(8) k1,k2        ! Coefficients of u expansion for small r:  u~k1+k2*sqrt(r)
  real(8) u_rk,su1,su2 ! Runge-Kutta coefficients for u.
  real(8) v_rk,sv1,sv2 ! Runge-Kutta coefficients for v.
  real(8) psifar       ! Asymptotic value of psi.
  real(8) mass         ! Black hole mass.
  real(8) epsilon      ! Tolerance.
  real(8) rm,aux       ! Auxiliary quantities.

  real(8), allocatable, dimension(:) :: r    ! Radial coordinate.
  real(8), allocatable, dimension(:) :: u    ! Solution of equation.
  real(8), allocatable, dimension(:) :: v    ! du/dr.

  real(8), allocatable, dimension(:) :: psi  ! Conformal factor.


! **************************
! ***   GET PARAMETERS   ***
! **************************

! Grid spacing.

  print *
  print *, 'Give grid spacing dr'
  read(*,*) dr

  hdr = 0.5d0*dr

! Number of grid points.

  print *
  print *, 'Give total number of grid points Nr'
  read(*,*) Nr

! Coefficient of sqrt(r) in expansion of u.

  print *
  print *, 'Coefficient of sqrt(r) in expansion of u'
  read(*,*) k2

  print *


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

! Allocate arrays.

  allocate(r(0:Nr))
  allocate(u(0:Nr))
  allocate(v(0:Nr))

  allocate(psi(0:Nr))


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

! Notice that we stagger the origin.

  do i=0,Nr
     r(i) = (dble(i)-0.5d0)*dr
  end do


! ************************
! ***   INITIAL DATA   ***
! ************************

  k1 = 1.d0

! Initialize u and v.

  u = k1
  v = 0.d0

! First grid point.  For this we use the fact that
! for small r the solution goes as:
!
! u  ~  K1 + k2 r^sqrt(2)
!
! so that:
!
! v  ~ k2 sqrt(2) r^(sqrt(2)-1)
!
! with "a" a constant that will be used to find the
! correct asymptotic behaviour.
!
! Also, remember that the first grid point is at r=dr/2.

  aux = dsqrt(2.d0)

  u(1) = k1 + k2*hdr**aux
  v(1) = k2*aux*hdr**(aux-1.d0)


! *****************************
! ***   START INTEGRATION   ***
! *****************************

! We start integrating from point r=dr/2 outwards.
! For the moment we use second order Runge-Kutta.

  do i=2,Nr

!    Advance half a time step.

     su1 = v(i-1)
     sv1 = - v(i-1)/r(i-1) + 0.25d0*(u(i-1) - k1/u(i-1)**7)/r(i-1)**2

     u_rk = u(i-1) + hdr*su1
     v_rk = v(i-1) + hdr*sv1

!    Advance full time step.

     rm = 0.5d0*(r(i-1)+r(i))

     su2 = v_rk
     sv2 = - v_rk/rm + 0.25d0*(u_rk - k1/u_rk**7)/rm**2

     u(i) = u(i-1) + dr*su2
     v(i) = v(i-1) + dr*sv2

  end do

! Symmetries for ghost point at r=-dr/2.
! This is in fact not used, but makes
! for nicer plots.

  u(0) = u(1)
  v(0) = v(1)

! Conformal factor.

  psi = u/sqrt(abs(r))

! Find asymptotic value of psi.

  psifar = psi(Nr) + r(Nr)*(psi(Nr) - psi(Nr-1))/dr

  print *, 'Asymptotic value of psi: ',psifar

! Find mass term.  Remember that for Schwarzschild
! asymptotically we must have psi ~ 1 + M/2r, so
! that the mass is:
!
! M  ~ 2 r (psi-1)
!
! I use the asymptotic value of psi, which should be equal
! to 1 within the tolerance if we did the shooting correctly.

  mass = 2.d0*r(Nr)*(psi(Nr) - psifar)

  print *, 'Mass parameter:          ',mass
  print *


! *********************
! ***   SAVE DATA   ***
! *********************

! Save solution for u.

  open(1,file='trumpet_u.rl',form='formatted',status='replace')

  do i=0,Nr
     if (dabs(u(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") r(i),u(i)
     else
        write(1,"(2ES16.8)") r(i),0.0D0
     end if
  end do

  close(1)

! Save v=du/dr.

  open(1,file='trumpet_v.rl',form='formatted',status='replace')

  do i=0,Nr
     if (dabs(v(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") r(i),v(i)
     else
        write(1,"(2ES16.8)") r(i),0.0D0
     end if
  end do

  close(1)

! Save conformal factor psi.

  open(1,file='trumpet_psi.rl',form='formatted',status='replace')

  do i=0,Nr
     if (dabs(psi(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") r(i),psi(i)
     else
        write(1,"(2ES16.8)") r(i),0.0D0
     end if
  end do

  close(1)


! ***************
! ***   END   ***
! ***************

  print *, 'PROGRAM TRUMPET HAS FINISHED'
  print *
  print *, 'Have a nice day!'
  print *
  print *

  end program trumpet


