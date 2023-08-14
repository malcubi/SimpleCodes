
 program tortoise

! ****************************
! ***   PROGRAM TORTOISE   ***
! ****************************

! This program solves the stationary Klein-Gordon
! equation for a scalar field around a black hole
! in tortoise coordinates.  We take c=G=hbar=1.
!
! We assume the the wave function is of the form:
!
! psi = exp(iwt) u(r)
!
! The equation to solve is then
!
!  2                 2
! d  u  =  ( V(r) - w  ) u
!  r*
!
! with r* the tortoise coordinate which is given
! in terms of the Schwarzschild radius as:
!
! r*  =  r + 2 M ln (r/2M - 1)
!
! and with the effective potential V(r) given by:
!
!                              2          3      2
! V(r) =  (1 - 2M/r) [ l(l+1)/r  +  2M / r  +  mu  ]
!
!
! Notice that this potential is given in terms of the
! Schwarzschild radius.

! Declare variables.

  implicit none

  logical minimum      ! Find local minimum of external amplitude?
  logical map          ! Map the Energy space?
  logical bracket      ! Did we bracket the minimum?

  integer i,j          ! Counters.
  integer Nr           ! Total number of grid points.
  integer l            ! Angular momentum quantum number.

  real(8) dr           ! Grid spacing.
  real(8) rmin         ! Minimum value of rstar.
  real(8) w            ! Frequency.
  real(8) dw           ! Step in frequency.
  real(8) M            ! Black hole mass.
  real(8) mu           ! Mass parameter.
  real(8) epsilon      ! Tolerance.

  real(8) wl,wc,wr
  real(8) A,Al,Ac,Ar,A1,A2

  real(8) zero,one,two

  real(8), allocatable, dimension(:) :: rstar   ! Tortoise radius.
  real(8), allocatable, dimension(:) :: rschw   ! Schwarzschild radius.
  real(8), allocatable, dimension(:) :: V       ! Potential.
  real(8), allocatable, dimension(:) :: u       ! Wave function.
  real(8), allocatable, dimension(:) :: D1_u    ! du/dr.


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.0D0
  one  = 1.0D0
  two  = 2.0D0


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

! Minimum value of rstar.

  print *
  print *, 'Give the minimum value of rstar'
  read(*,*) rmin

! Black hole mass.

  print *
  print *, 'Give the black hole mass M'
  read(*,*) M

! Mass parameter.

  print *
  print *, 'Give the mass parameter mu'
  read(*,*) mu

! Angular momentum quantum number.

  print *
  print *, 'Give the angular momentum quantum number l'
  read(*,*) l

! Frequency.

  print *
  print *, 'Give the frequency of the solution w'
  read(*,*) w

! Are we trying to look for a local minimum of external amplitude?

  print *
  print *, 'Do we look for minimum of external amplitude?'
  read(*,'(L1)') minimum

! Do we want to map the frequency space?

  print *
  print *, 'Do we map the frequency space?'
  read(*,'(L1)') map

! Step in frequency.

  print *
  print *, 'Give the step in frequency dw'
  read(*,*) dw

  print *


! ************************
! ***   SANITY CHECK   ***
! ************************

! Check value of frequency.

  if (w>=mu) then
     print *
     print *, 'The frequency w must be smaller than mu'
     print *
     stop
  end if

! Check what we want to do.

  if ((minimum).and.(map)) then
     print *
     print *, 'It is incompatible to have both minimum=.true. and map=.true.'
     print *
     stop
  end if


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

! Allocate arrays.

  allocate(rstar(0:Nr),rschw(0:Nr))
  allocate(u(0:Nr),D1_u(0:Nr),V(0:Nr))


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do i=0,Nr
     rstar(i) = rmin + dble(i)*dr
  end do


! *************************************
! ***   FIND SCHWARZSCHILD RADIUS   ***
! *************************************

! Find the Schwarzschild radius.

  call schwarz(Nr,M,rstar,rschw)


! **************************
! ***   FIND POTENTIAL   ***
! **************************

  do i=0,Nr
     V(i) = (1.0D0 - 2.0D0*M/rschw(i))*(dble(l*(l+1))/rschw(i)**2 &
          + 2.0D0*M/rschw(i)**3 + mu**2)
  end do


! ***************************************
! ***   SOLVE KLEIN-GORDON EQUATION   ***
! ***************************************

! Just a simple solution.

  if ((.not.minimum).and.(.not.map)) then

     call kleingordon(Nr,dr,w,mu,A,rstar,u,D1_u,V)

! Here we try to find a minimum of the external amplitude.

  else if (minimum) then

!    First we need to bracket the minimum.

     epsilon = 1.0D-10

     bracket = .false.

!    Find first 3 values of A.

     wl = w
     wc = wl + dw
     wr = wc + dw

     call kleingordon(Nr,dr,wl,mu,Al,rstar,u,D1_u,V)
     call kleingordon(Nr,dr,wc,mu,Ac,rstar,u,D1_u,V)
     call kleingordon(Nr,dr,wr,mu,Ar,rstar,u,D1_u,V)

!    Now locate minimum.

     j = 0

     do while ((dw>epsilon).and.(w<mu))

        j = j+1
!       print *, j,wc,Al,Ac,Ar,dw

!       We have bracketed the minimum

        if ((Ac<Al).and.(Ac<Ar)) then

           bracket = .true.

           dw = 0.5D0*dw

           w = wc - dw
           call kleingordon(Nr,dr,w,mu,A1,rstar,u,D1_u,V)

           w = wc + dw
           call kleingordon(Nr,dr,w,mu,A2,rstar,u,D1_u,V)

           if ((A1<Al).and.(A1<Ac)) then
              Ar = Ac
              Ac = A1
              wc = wc - dw
           else if ((A2<Ac).and.(A2<Ar)) then
              Al = Ac
              Ac = A2
              wc = wc + dw
           else
              Al = A1
              Ar = A2
           end if

!       We have not yet bracketed the minimum.

        else

           Al = Ac
           Ac = Ar

           wl = wc
           wc = wr
           wr = wr + dw

           call kleingordon(Nr,dr,wr,mu,Ar,rstar,u,D1_u,V)

        end if

     end do

!    Find solution minimum.

     call kleingordon(Nr,dr,wc,mu,A,rstar,u,D1_u,V)

     print *
     print *, 'Minimum found for frequency  w =',wc,',  w**2 = ',wc**2, &
              ', (w/mu) = ',W/mu
     print *, 'External wave amplitude   A =',A
     print *

! Here we map the energy space.

  else if (map) then

     j = 0

     open(1,file='As.tl',form='formatted',status='replace')

     do while (w<mu)

        j = j+1

        call kleingordon(Nr,dr,w,mu,A,rstar,u,D1_u,V)
        write(1,"(2ES16.8)") w,A

        w = w+dw

     end do

     close(1)

  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  if (.not.map) then

!    Schwarzschild radius.

     open(1,file='rschw.rl',form='formatted',status='replace')

     do i=0,Nr
        if (dabs(rschw(i)) > 1.0D-50) then
           write(1,"(2ES16.8)") rstar(i),rschw(i)
        else
           write(1,"(2ES16.8)") rstar(i),0.0D0
        end if
     end do

     close(1)

!    Potential.

     open(1,file='Vstar.rl',form='formatted',status='replace')

     do i=0,Nr
        if (dabs(V(i)) > 1.0D-50) then
           write(1,"(2ES16.8)") rstar(i),V(i)
        else
           write(1,"(2ES16.8)") rstar(i),0.0D0
        end if
     end do

     close(1)

!    Wave function.

     open(1,file='ustar.rl',form='formatted',status='replace')

     do i=0,Nr
        if (dabs(u(i)) > 1.0D-50) then
           write(1,"(2ES16.8)") rstar(i),u(i)
        else
           write(1,"(2ES16.8)") rstar(i),0.0D0
        end if
     end do

     close(1)

!    rstar from rscwh.

     open(1,file='rstar.rl',form='formatted',status='replace')

     do i=0,Nr
        write(1,"(2ES16.8)") rstar(i),rschw(i)+2.0D0*M*log(0.5D0*rschw(i)/M-1.0D0)
     end do

     close(1)

  end if


! ***************
! ***   END   ***
! ***************

  end program tortoise







  subroutine kleingordon(Nr,dr,w,mu,A,rstar,u,D1_u,V)

! ***************************************
! ***   SOLVE KLEIN-GORDON EQUATION   ***
! ***************************************

! Here we solve the Klein-Gordon equation by writing
! it in first order form and using Runge-Kutta.

  implicit none

  integer i,Nr

  real(8) dr,w,mu
  real(8) A,Anorm,q
  real(8) su,sup,uh,uph

  real(8), dimension(0:Nr) :: rstar
  real(8), dimension(0:Nr) :: u,D1_u
  real(8), dimension(0:Nr) :: V


! *****************************
! ***   INITIALIZE ARRAYS   ***
! *****************************

  u = 0.0D0
  D1_u = 0.0D0


! **************************
! ***   SOLVE EQUATION   ***
! **************************

! Find parameter q for exponential solution for large r.

  if (w<mu) then
     q = sqrt(mu**2-w**2)
  else
     print *
     print *, 'w>=mu.  This should not happen.'
     print *
     stop
  end if

! At the outer boundary we set the wave function
! and its derivative to decaying exponential.

  u(Nr) = exp(-q*rstar(Nr))
  D1_u(Nr) = -q*u(Nr)

! Second order Runge-Kutta.  Notice that we integrate
! from Nr to 0 (backward).

  A = 0.0D0
  Anorm = 0.0D0

  do i=Nr-1,0,-1

!    Calculate sources.

     su  = D1_u(i+1)
     sup = (V(i+1) - w**2)*u(i+1)

!    Advance half a step.

     uh  = u(i+1)  - 0.5D0*dr*su
     uph = D1_u(i+1) - 0.5D0*dr*sup

!    Calculate sources again.

     su  = uph
     sup = (0.5D0*(V(i)+V(i+1)) - w**2)*uh

!    Advanced full step.

     u(i) = u(i+1) - dr*su
     D1_u(i) = D1_u(i+1) - dr*sup

!    Find maximum overall amplitude.

     if (rstar(i)>10.0D0) then
        Anorm = max(Anorm,abs(u(i)))
     end if

!    Find maximum amplitude for r<10 (outside of the barrier.

     if (rstar(i)<-10.0D0) then
        A = max(A,abs(u(i)))
     end if

  end do

! Normalize solution.

  u = u/Anorm
  A = A/Anorm


! ***************
! ***   END   ***
! ***************

  end subroutine kleingordon







  subroutine schwarz(Nr,M,rstar,rschw)

! ************************************
! ***   FIND SCWARZSCHILD RADIUS   ***
! ************************************

! Here we invert the relation:
!
! r*  =  r + 2 M ln (r/2M - 1)
!
! In order to solve for the Schwarzschild radius
! in terms of the tortoise radius.

  implicit none

  integer i,j,jmax
  integer Nr

  real(8) M
  real(8) dr,epsilon
  real(8) rs,resl,resr,resc

  real(8), dimension(0:Nr) :: rstar
  real(8), dimension(0:Nr) :: rschw

! Set tolerance.

  epsilon = 1.0D-10

! Loop over values of rstar.

  do i=0,Nr

     dr = 1.0D-7

!    Initial guess.

     if (i==0) then
        rs = 2.0D0*M + epsilon
     else
        rs = rschw(i-1)
     end if

!    Bracket the minimum.

     resl = rs + 2.0D0*M*dlog(0.5D0*rs/M - 1.0D0) - rstar(i)
     resr = (rs+dr) + 2.0D0*M*dlog(0.5D0*(rs+dr)/M - 1.0D0) - rstar(i)

     j = 0
     jmax = 100000

     do while((dr>epsilon).and.(j<jmax))

        j=j+1

!       Did we bracket the minimum?

        if (resl*resr>0.0D0) then

!          We have not bracketed the minimum.
!          Increase dr and move to the right.

           dr = 1.005D0*dr
           rs = rs + dr

           resl = resr
           resr = (rs+dr) + 2.0D0*M*dlog(0.5D0*(rs+dr)/M - 1.0D0) - rstar(i)

        else

!          We have bracketed the minimum.

           dr = 0.5D0*dr

           resc = (rs+dr) + 2.0D0*M*dlog(0.5D0*(rs+dr)/M - 1.0D0) - rstar(i)

           if (resl*resc>0.0D0) then
              rs = rs + 3.0D0*dr
              resl=resc
              resc = rs + 2.0D0*M*dlog(0.5D0*rs/M - 1.0D0) - rstar(i)
           else
              rs = rs + dr
              resr=resc
              resc = rs + 2.0D0*M*dlog(0.5D0*rs/M - 1.0D0) - rstar(i)
           end if

        end if

     end do

!    Check if we exceded maximum number of iterations.

     if (j>=jmax) then
        print *,'Maximum number of iterations exceded at point: ',i,j,dr
        stop
     end if

!    Find Schwarzschild radius.

     rschw(i) = rs
     resc =  rschw(i) + 2.0D0*M*log(0.5D0*rschw(i)/M - 1.0D0) - rstar(i)
!    print *, i,rschw(i),resc,dr

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine schwarz
