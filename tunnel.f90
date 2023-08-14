
  program tunnel

! *************************
! ***   TUNNEL EFFECT   ***
! *************************

! This program solves the Schroedinger equation
! for a simple tunnel effect.  We take both the mass
! of the particle and hbar equal to 1.
!
! The equation to solve is:
!
!     2
!    d phi  =   2 (V - E) phi
!     x
!
! with E the energy of the mode, and with the
! potential V given by:
!
!    x  <= 0              V = infinity
!    0  < x < xb          V = Vmin
!    xb < x < xb + d      V = Vmax
!    xb + d < x           V = 0

! Declare variables.

  implicit none

  logical minimum      ! Find local minimum of external amplitude?
  logical map          ! Map the Energy space?
  logical bracket      ! Did we bracket the minimum?

  integer i,j          ! Counters.
  integer Nx           ! Total number of grid points.

  real(8) dx           ! Grid spacing.
  real(8) E            ! Energy of the solution.
  real(8) dE           ! Step in energy.
  real(8) xb           ! Barrier position.
  real(8) d            ! Barrier width.
  real(8) Vmin         ! Minimum of potential well.
  real(8) Vmax         ! Height of potential barrier.
  real(8) epsilon      ! Tolerance.

  real(8) El,Ec,Er
  real(8) A,Al,Ac,Ar,A1,A2

  real(8) zero,half,one,two,smallpi

  real(8), allocatable, dimension(:) :: x       ! Spatial coordinate.
  real(8), allocatable, dimension(:) :: phi     ! Wave function.
  real(8), allocatable, dimension(:) :: psi     ! dphi/dx.
  real(8), allocatable, dimension(:) :: V       ! Potential.


! *******************
! ***   NUMBERS   ***
! *******************

  zero = 0.0D0
  half = 0.5D0
  one  = 1.0D0
  two  = 2.0D0

  smallpi = acos(-one)


! **************************
! ***   GET PARAMETERS   ***
! **************************

! Grid spacing.

  print *
  print *, 'Give grid spacing dx'
  read(*,*) dx

! Number of grid points.

  print *
  print *, 'Give total number of grid points Nx'
  read(*,*) Nx

! Energy.

  print *
  print *, 'Give the energy of the solution (E)'
  read(*,*) E

! Minimum of potential well.

  print *
  print *, 'Give minimum of potential well (Vmin)'
  read(*,*) Vmin

! Barrier height.

  print *
  print *, 'Give the barrier height (Vmax)'
  read(*,*) Vmax

! Barrier position.

  print *
  print *, 'Give the barrier position (xb)'
  read(*,*) xb

! Barrier width.

  print *
  print *, 'Give the barrier width (d)'
  read(*,*) d

! Are we trying to look for a local minimum of external amplitude?

  print *
  print *, 'Do we look for minimum of external amplitude?'
  read(*,'(L1)') minimum

! Do we want to map the energy space?

  print *
  print *, 'Do we map the energy space?'
  read(*,'(L1)') map

! Step in energy.

  print *
  print *, 'Give the step in energy (dE)'
  read(*,*) dE

  print *


! ************************
! ***   SANITY CHECK   ***
! ************************

! Check values of potential.

  if (Vmax<=Vmin) then
     print *
     print *, 'Vmax must be larger than Vmin'
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

  allocate(x(0:Nx))
  allocate(phi(0:Nx),psi(0:Nx),V(0:Nx))


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do i=0,Nx
     x(i) = dble(i)*dx
  end do


! **************************
! ***   FIND POTENTIAL   ***
! **************************

  do i=0,Nx

     if (x(i)<xb) then
        V(i)=Vmin
     else if ((x(i)>=xb).and.(x(i)<xb+d)) then
        V(i) = Vmax
     else if (x(i)>=xb+d) then
        v(i) = zero
     end if

  end do


! ***************************************
! ***   SOLVE SCHROEDINGER EQUATION   ***
! ***************************************

! Just a simple solution.

  if ((.not.minimum).and.(.not.map)) then

     call schroedinger(Nx,dx,E,xb,d,A,x,phi,psi,V)

! Here we try to find a minimum of the external amplitude.

  else if (minimum) then

!    First we need to bracket the minimum.

     epsilon = 1.0D-10

     bracket = .false.

!    Find first 3 values of A.

     El = E
     Ec = El + dE
     Er = Ec + dE

     call schroedinger(Nx,dx,El,xb,d,Al,x,phi,psi,V)
     call schroedinger(Nx,dx,Ec,xb,d,Ac,x,phi,psi,V)
     call schroedinger(Nx,dx,Er,xb,d,Ar,x,phi,psi,V)

!    Now locate minimum.

     j = 0

     do while ((dE>epsilon).and.(E<Vmax))

        j = j+1
        print *, j,Ec,Al,Ac,Ar,dE

!       We have bracketed the minimum

        if ((Ac<Al).and.(Ac<Ar)) then

           bracket = .true.

           dE = 0.5D0*dE

           E = Ec - dE
           call schroedinger(Nx,dx,E,xb,d,A1,x,phi,psi,V)

           E = Ec + dE
           call schroedinger(Nx,dx,E,xb,d,A2,x,phi,psi,V)

           if ((A1<Al).and.(A1<Ac)) then
              Ar = Ac
              Ac = A1
              Ec = Ec - dE
           else if ((A2<Ac).and.(A2<Ar)) then
              Al = Ac
              Ac = A2
              Ec = Ec + dE
           else
              Al = A1
              Ar = A2
           end if

!       We have not yet bracketed the minimum.

        else

           Al = Ac
           Ac = Ar

           El = Ec
           Ec = Er
           Er = Er + dE

           call schroedinger(Nx,dx,Er,xb,d,Ar,x,phi,psi,V)

        end if

     end do

!    Find solution minimum.

     call schroedinger(Nx,dx,Ec,xb,d,A,x,phi,psi,V)

     print *
     print *, 'Minimum found for energy  E =',Ec
     print *, 'External wave amplitude   A =',A
     print *

! Here we map the energy space.

  else if (map) then

     j = 0

     open(1,file='As.tl',form='formatted',status='replace')

     do while (E<Vmax)

        j = j+1

        call schroedinger(Nx,dx,E,xb,d,A,x,phi,psi,V)
        write(1,"(2ES16.8)") E,A
!       print *, j,E,A

        E = E+dE

     end do

     close(1)

  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  if (.not.map) then

!    Potential

     open(1,file='Vs.xl',form='formatted',status='replace')

     do i=0,Nx
        if (dabs(V(i)) > 1.0D-50) then
           write(1,"(2ES16.8)") x(i),V(i)
        else
           write(1,"(2ES16.8)") x(i),0.0D0
        end if
     end do

     close(1)

!    Wave function.

     open(1,file='phis.xl',form='formatted',status='replace')

     do i=0,Nx
        if (dabs(phi(i)) > 1.0D-50) then
           write(1,"(2ES16.8)") x(i),phi(i)
        else
           write(1,"(2ES16.8)") x(i),0.0D0
        end if
     end do

     close(1)

!    Derivative of wave function.

     open(1,file='psis.xl',form='formatted',status='replace')

     do i=0,Nx
        if (dabs(psi(i)) > 1.0D-50) then
           write(1,"(2ES16.8)") x(i),psi(i)
        else
           write(1,"(2ES16.8)") x(i),0.0D0
       end if
     end do

     close(1)

  end if


! ***************
! ***   END   ***
! ***************

  end program tunnel










  subroutine schroedinger(Nx,dx,E,xb,d,A,x,phi,psi,V)

! ***************************************
! ***   SOLVE SCHROEDINGER EQUATION   ***
! ***************************************

! Here we solve the Schroedinger equation by writing
! it in first order form and using Runge-Kutta.
! In order to rewrite it in first order form we
! introduce the auxiliary variable:  psi = dphi/dr.

  implicit none

  integer i,Nx

  real(8) dx
  real(8) k,E
  real(8) sphi,spsi
  real(8) phih,psih
  real(8) xb,d,A

  real(8), dimension(0:Nx) :: x
  real(8), dimension(0:Nx) :: phi
  real(8), dimension(0:Nx) :: psi
  real(8), dimension(0:Nx) :: V


! *****************************
! ***   INITIALIZE ARRAYS   ***
! *****************************

  phi = 0.0D0
  psi = 0.0D0


! **************************
! ***   SOLVE EQUATION   ***
! **************************

! Initialize wave function and derivative at x=0.

  k = sqrt(2.0D0*E)

  phi(0) = 0.0d0
  psi(0) = k

! Second order Runge-Kutta

  A = 0.0D0

  do i=1,Nx

!    Calculate sources.

     sphi = psi(i-1)
     spsi = 2.0D0*(V(i-1) - E)*phi(i-1)

!    Advance half a step.

     phih = phi(i-1) + 0.5D0*dx*sphi
     psih = psi(i-1) + 0.5D0*dx*spsi
 
!    Calculate sources again.

     sphi = psih
     spsi = 2.0D0*(0.5D0*(V(i-1)+V(i)) - E)*phih

!    Advanced full step.

     phi(i) = phi(i-1) + dx*sphi
     psi(i) = psi(i-1) + dx*spsi

!    Find maximum external amplitude.

     if (x(i)>xb+d) then
        A = max(A,abs(phi(i)))
     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine schroedinger
