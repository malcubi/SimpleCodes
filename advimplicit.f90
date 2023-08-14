!$Header: /usr/local/ollincvs/Codes/SimpleCodes/advimplicit.f90,v 1.22 2017/03/21 19:44:10 malcubi Exp $

  program advimplicit

! This is a simple code for the advection equation
! that uses implicit methods.

! ******************************
! ***   ADVECTION EQUATION   ***
! ******************************

! Declare variables.

  implicit none

  integer i,l           ! Counters
  integer Nt            ! Total number of time steps.
  integer Nx            ! Total number of grid points.
  integer Noutput       ! Frequency of output.

  real(8) dx            ! Grid spacing.
  real(8) dt            ! Time step.
  real(8) rho,rhoh      ! Courant parameter rho=dt/dx, rhoh=rho/2.
  real(8) t             ! Time.
  real(8) v             ! Wave speed.
  real(8) theta         ! Implicit parameter.

  real(8) x0            ! Center of initial gaussian.
  real(8) s0            ! Width of initial gaussian.
  real(8) a0            ! Amplitude of initial gaussian.

  real(8) sloper,slopel ! Right and left slopes for TVD schemes.

  real(8) aux

  character(20) method  ! Integration method.

  real(8), allocatable, dimension(:) :: x       ! Position.
  
  real(8), allocatable, dimension(:) :: phi     ! Wave function at current time step.
  real(8), allocatable, dimension(:) :: phi_p   ! Wave function at previous time step.

  real(8), allocatable, dimension (:,:) :: M    ! Matrix elements defining the implicit system.
  real(8), allocatable, dimension (:) :: S      ! Source vector.

  real(8), allocatable, dimension (:) :: flux   ! Fluxes for TVD schemes.
  real(8), allocatable, dimension (:,:) :: C    ! Coefficients for calculating fluxes.


! **************************
! ***   GET PARAMETERS   ***
! **************************

  print *
  print *, 'Give grid spacing dx'
  read(*,*) dx

  print *
  print *, 'Give time step dt'
  read(*,*) dt

  print *
  print *, 'Give wave speed'
  read(*,*) v
  
  print *
  print *, 'Give total number of grid points Nx'
  read(*,*) Nx

  print *
  print *, 'Give total number of time steps'
  read(*,*) Nt

  print *
  print *, 'Give amplitude of initial gaussian'
  read(*,*) a0

  print *
  print *, 'Give width of initial gaussian'
  read(*,*) s0

  print *
  print *, 'Give position of initial gaussian'
  read(*,*) x0

  print *
  print *, 'Give frequency of output'
  read(*,*) Noutput

  print *
  print *, 'Give integration method (UP=upwind,AUP=averagedupwind,CN=cranknicolson,MM=minmod,IE-CN,IE-UP)'
  read(*,*) method

  print *
  print *, 'Give implicit parameter theta [0,1] (theta=0 explicit, theta=1 full implicit)'
  read(*,*) theta

  print *

  if ((theta<0.d0).or.(theta>1.d0)) then
     print *, 'The implicit parameter theta must be between 0 and 1.'
     print *, 'Aborting ...'
     print *
     stop
  end if


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

  allocate(x(0:Nx),phi(0:Nx),phi_p(0:Nx))
  allocate(M(0:Nx,-2:2),S(0:Nx))

  allocate(C(0:Nx,-2:2),flux(0:Nx))


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do i=0,Nx
     x(i) = dble(i)*dx
  end do


! ****************************
! ***   OUTPUT TO SCREEN   ***
! ****************************

  print *,'------------------------------'
  print *,'|  Time step  |     Time     |'
  print *,'------------------------------'


! ************************
! ***   INITIAL DATA   ***
! ************************

! Initialize time.

  t = 0.0D0

! Calculate Courant parameter.

  rho  = dt/dx
  rhoh = 0.5d0*rho

! Initial data (gaussian).

  do i=0,Nx
     phi(i) = a0*exp(-(x(i)-x0)**2/s0**2)
  end do

! Output to screen.

  write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t,'  | '


! *****************************
! ***   OPEN OUTPUT FILES   ***
! *****************************

  open(1,file='phi_advimp.xl',form='formatted',status='replace')
  open(2,file='phi_tv.tl',form='formatted',status='replace')


! *****************************
! ***   SAVE INITIAL DATA   ***
! *****************************

! Save wave function.

  write(1,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(phi(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") x(i),phi(i)
     else
        write(1,"(2ES16.8)") x(i),0.0D0
     end if
  end do

  write(1,*)

! Save total variation.

  aux = 0.d0

  do i=1,Nx
     aux = aux + abs(phi(i)-phi(i-1))
  end do

  write(2,*) '"phi_tv.tl'

  if (abs(aux)>1.d-50) then
     write(2,"(2ES16.8)") t,aux
  else
     write(2,"(2ES16.8)") t,0.d0
  end if


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do l=1,Nt

!    Time.

     t = t + dt

!    Save old time step.

     phi_p = phi


!    ******************
!    ***   UPWIND   ***
!    ******************

!    Here we use upwinded differences in space,
!    averaged between the old and new time level with
!    weight theta.
!
!    The computational molecule has the form:
!
!    n+1   *-----*               *-----*
!                |               |
!                |  for v>0      |       for v<0
!                |               |
!    n     *-----*               *-----*
!         i-1    i               i    i+1
!
!    Special cases are:
!
!    theta = 0:      Standard explicit upwind (first order in both time and space).
!    theta = 1/2:    Averaged uwpind (second order in time, first order in space).
!    theta = 1       Fully implicit upwind (first order in both time and space).

     if (method=='UP') then

!       Initialize (M,S) to zero.

        M = 0.d0
        S = 0.d0

!       Positive speed.

        if (v>=0.d0) then

!          Source term, all the way to right boundary.

           do i=1,Nx
              S(i) = phi_p(i) - v*rho*(1.d0-theta)*(phi_p(i) - phi_p(i-1))
           end do

!          Matrix components, all the way to right boundary.

           do i=1,Nx
              M(i,-1) = - v*rho*theta
              M(i, 0) = 1.d0 + v*rho*theta
              M(i,+1) = 0.d0
           end do

!          Set phi at left boundary to zero.

           M(0,0) = 1.d0

!       Negative speed.

        else

!          Source term, all the way to left boundary.

           do i=0,Nx-1
              S(i) = phi_p(i) - v*rho*(1.d0-theta)*(phi_p(i+1) - phi_p(i))
           end do

!          Matrix components, all the way to left boundary.

           do i=0,Nx-1
              M(i,-1) = 0.d0
              M(i, 0) = 1.d0 - v*rho*theta
              M(i,+1) = + v*rho*theta
           end do

!          Set phi at right boundary to zero.

           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S


!    ***************************
!    ***   AVERAGED UPWIND   ***
!    ***************************

!    Here we use upwinded differences in space 
!    averaged between the old and new time level with
!    weight theta.  We also use averaged differences
!    in time along the upwinded direcion.
!
!!    The computational molecule has the form:
!
!    n+1   *-----*               *-----*
!          |     |               |     |
!          |     |  for v>0      |     |  for v<0
!          |     |               |     |
!    n     *-----*               *-----*
!         i-1    i               i    i+1
!
!    The special case theta=0.5 is the method used at the
!    boundaries, and is second order in both time and space.

     else if (method=='AUP') then

!       Initialize (M,S) to zero.

        M = 0.d0
        S = 0.d0

!       Positive speed.

        if (v>=0.d0) then

!          Source term, all the way to right boundary.

           do i=1,Nx
              S(i) = phi_p(i  )*(1.d0 - 2.d0*v*rho*(1.d0-theta)) &
                   + phi_p(i-1)*(1.d0 + 2.d0*v*rho*(1.d0-theta))
           end do

!          Matrix components, all the way to right boundary.

           do i=1,Nx
              M(i,-1) = 1.d0 - 2.d0*v*rho*theta
              M(i, 0) = 1.d0 + 2.d0*v*rho*theta
              M(i,+1) = 0.d0
           end do

!          Set phi at left boundary to zero.

           M(0,0) = 1.d0

!       Negative speed

        else

!          Source term, all the way to left boundary.

           do i=0,Nx-1
              S(i) = phi_p(i+1)*(1.d0 - 2.d0*v*rho*(1.d0-theta)) &
                   + phi_p(i  )*(1.d0 + 2.d0*v*rho*(1.d0-theta))
           end do

!          Matrix components, all the way to left boundary.

           do i=0,Nx-1
              M(i,-1) = 0.d0
              M(i, 0) = 1.d0 - 2.d0*v*rho*theta
              M(i,+1) = 1.d0 + 2.d0*v*rho*theta
           end do

!          Set phi at right boundary to zero.

           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S


!    **************************
!    ***   CRANK-NICOLSON   ***
!    **************************

!    Here we used centered differences in space,
!    averaged between the old and new time level with
!    weight theta.
!
!    The computational molecule has the form:
!
!    n+1   *-----*-----*
!                |
!                |
!                | 
!    n     *-----*–––––*
!         i-1    i    i+1
!
!    Special cases are:
!
!    theta = 0:      Forward Euler (first order in time, second in space).
!    theta = 1/2:    Crank-Nicolson (second order in both time and space).
!    theta = 1:      Backward Euler (first order in time, second in space).

     else if (method=='CN') then

!       Initialize (M,S) to zero.

        M = 0.d0
        S = 0.d0

!       Source term.

        do i=1,Nx-1
           S(i) = phi_p(i) - 0.5d0*v*rho*(1.d0-theta)*(phi_p(i+1) - phi_p(i-1))
        end do

!       Matrix components.

        do i=1,Nx-1
           M(i,-1) = - 0.5d0*v*rho*theta
           M(i, 0) = 1.d0
           M(i,+1) = + 0.5d0*v*rho*theta
        end do

!       Boundaries using upwind.

        if (v>=0.d0) then

!          Left boundary.

           M(0,0) = 1.d0

!          Right boundary.

           S(Nx) = phi_p(Nx) - v*rho*(1.d0-theta)*(phi_p(Nx) - phi_p(Nx-1))

           M(Nx,-1) = - v*rho*theta
           M(Nx, 0) = 1.d0 + v*rho*theta
           M(Nx,+1) = 0.d0

        else

!          Left boundary.

           S(0) = phi_p(0) - v*rho*(1.d0-theta)*(phi_p(1) - phi_p(0))

           M(0,-1) = 0.d0
           M(0, 0) = 1.d0 - v*rho*theta
           M(0,+1) = + v*rho*theta

!          Right boundary.

           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S


!    ******************
!    ***   MINMOD   ***
!    ******************

!    This is implicit minmod.
!
!    Notice that explicit minmod (theta=0) requires a stronger
!    Courant condition:  dt < dx/2.

     else if (method=='MM') then

!       Initialize (M,S) to zero.

        M = 0.d0
        S = 0.d0
        C = 0.d0

!       Positive speed.

        if (v>=0.d0) then

!          Calculate fluxes.

           flux(0) = 0.5d0*v*(phi_p(0) + phi_p(1))

           do i=1,Nx-1

!             Calculate right and left slopes.

              sloper = phi_p(i+1) - phi_p(i  )
              slopel = phi_p(i  ) - phi_p(i-1)

!             Find coefficients for flux.

              if (sloper*slopel>0.d0) then
                 if (dabs(sloper)<dabs(slopel)) then
                    C(i,-1) = 0.d0
                    C(i, 0) = 0.5d0
                    C(i,+1) = 0.5d0
                 else
                    C(i,-1) = - 0.5d0
                    C(i, 0) = + 1.5d0
                    C(i,+1) = 0.d0
                 end if
              else
                 C(i,-1) = 0.d0
                 C(i, 0) = 1.d0
                 C(i,+1) = 0.d0
              end if

!             Calculate explicit fluxes.

              flux(i) = v*(C(i,-1)*phi_p(i-1) + C(i,0)*phi_p(i) + C(i,1)*phi_p(i+1))

           end do

!       Negative speed.

        else

!          Calculate fluxes.

           flux(Nx-1) = 0.5d0*v*(phi_p(Nx) + phi_p(Nx-1))

           do i=0,Nx-2

!             Calculate right and left slopes.

              sloper = phi_p(i+2) - phi_p(i+1)
              slopel = phi_p(i+1) - phi_p(i  )

              if (sloper*slopel>0.d0) then
                 if (dabs(sloper)<dabs(slopel)) then
                    C(i,0) = 0.d0
                    C(i,1) = 1.5d0
                    C(i,2) = -0.5d0
                 else
                    C(i,0) = 0.5d0
                    C(i,1) = 0.5d0
                    C(i,2) = 0.d0
                 end if
              else
                 C(i,0) = 0.d0
                 C(i,1) = 1.d0
                 C(i,2) = 0.d0
              end if

!             Calculate explicit fluxes.

              flux(i) = v*(C(i,0)*phi_p(i) + C(i,1)*phi_p(i+1) + C(i,2)*phi_p(i+2))

           end do

        end if

!       Find source term.

        do i=1,Nx-1
           S(i) = phi_p(i) - rho*(1.d0-theta)*(flux(i) - flux(i-1))
        end do

!       Matrix components.

        aux = v*rho*theta

        do i=1,Nx-1
           M(i,-2) = - aux*C(i-1,-1)
           M(i,-1) = + aux*(C(i,-1) - C(i-1,0))
           M(i, 0) = + aux*(C(i, 0) - C(i-1,1)) + 1.d0
           M(i,+1) = + aux*(C(i, 1) - C(i-1,2))
           M(i,+2) = + aux*C(i,2)
        end do

!       Boundaries using upwind.

        if (v>=0.d0) then

!          Left boundary.

           S(0) = 0.d0
           M(0,0) = 1.d0

!          Right boundary.

           S(Nx) = phi_p(Nx) - v*rho*(1.d0-theta)*(phi_p(Nx) - phi_p(Nx-1))

           M(Nx,-1) = - v*rho*theta
           M(Nx, 0) = 1.d0 + v*rho*theta
           M(Nx,+1) = 0.d0

        else

!          Left boundary.

           S(0) = phi_p(0) - v*rho*(1.d0-theta)*(phi_p(1) - phi_p(0))

           M(0,-1) = 0.d0
           M(0, 0) = 1.d0 - v*rho*theta
           M(0,+1) = + v*rho*theta

!          Right boundary.

           S(Nx) = 0.d0
           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S


!    ************************************
!    ***   IMPLICIT-EXPLICIT UPWIND   ***
!    ************************************

!    Here we first do half an implicit time step, and then
!    a full explicit time step calculating the sources with
!    the result of the implicit half time step.

     else if  (method=='IE-UP') then

!       HALF A TIME STEP WITH IMPLICIT METHOD.

        M = 0.d0
        S = 0.d0

!       Positive speed.

        if (v>=0.d0) then

!          Source term, all the way to right boundary.

           do i=1,Nx
              S(i) = phi_p(i) - v*rhoh*(1.d0-theta)*(phi_p(i) - phi_p(i-1))
           end do

!          Matrix components, all the way to right boundary.

           do i=1,Nx
              M(i,-1) = - v*rhoh*theta
              M(i, 0) = 1.d0 + v*rhoh*theta
              M(i,+1) = 0.d0
           end do

!          Set phi at left boundary to zero.

           M(0,0) = 1.d0

!       Negative speed.

        else

!          Source term, all the way to left boundary.

           do i=0,Nx-1
              S(i) = phi_p(i) - v*rhoh*(1.d0-theta)*(phi_p(i+1) - phi_p(i))
           end do

!          Matrix components, all the way to left boundary.

           do i=0,Nx-1
              M(i,-1) = 0.d0
              M(i, 0) = 1.d0 - v*rhoh*theta
              M(i,+1) = + v*rhoh*theta
           end do

!          Set phi at right boundary to zero.

           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S

!       FULL TIME STEP WITH EXPLICIT METHOD.

        M = 0.d0
        S = 0.d0

!       Positive speed.

        if (v>=0.d0) then

!          Source term, all the way to right boundary.

           do i=1,Nx
              S(i) = phi_p(i) - v*rho*(phi(i) - phi(i-1))
           end do

!          Matrix components, all the way to right boundary.

           do i=1,Nx
              M(i,-1) = 0.d0
              M(i, 0) = 1.d0
              M(i,+1) = 0.d0
           end do

!          Set phi at left boundary to zero.

           M(0,0) = 1.d0

!       Negative speed.

        else

!          Source term, all the way to left boundary.

           do i=0,Nx-1
              S(i) = phi_p(i) - v*rho*(phi(i+1) - phi(i))
           end do

!          Matrix components, all the way to left boundary.

           do i=0,Nx-1
              M(i,-1) = 0.d0
              M(i, 0) = 1.d0
              M(i,+1) = 0.d0
           end do

!          Set phi at right boundary to zero.

           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S


!    ********************************************
!    ***   IMPLICIT-EXPLICIT CRANK-NICOLSON   ***
!    ********************************************

!    Here we first do half an implicit time step, and then
!    a full explicit time step calculating the sources with
!    the result of the implicit half time step.
!
!    For theta=1, this seems to inherit the unconditional
!    stability from the implict method while becoming second
!    order in time.

     else if  (method=='IE-CN') then

!       HALF A TIME STEP WITH IMPLICIT METHOD.

        M = 0.d0
        S = 0.d0

!       Source term.

        do i=1,Nx-1
           S(i) = phi_p(i) - 0.5d0*v*rhoh*(1.d0-theta)*(phi_p(i+1) - phi_p(i-1))
        end do

!       Matrix components.

        do i=1,Nx-1
           M(i,-1) = - 0.5d0*v*rhoh*theta
           M(i, 0) = 1.d0
           M(i,+1) = + 0.5d0*v*rhoh*theta
        end do

!       Boundaries using upwind.

        if (v>=0.d0) then

!          Left boundary.

           M(0,0) = 1.d0

!          Right boundary.

           S(Nx) = phi_p(Nx) - v*rhoh*(1.d0-theta)*(phi_p(Nx) - phi_p(Nx-1))

           M(Nx,-1) = - v*rhoh*theta
           M(Nx, 0) = 1.d0 + v*rhoh*theta
           M(Nx,+1) = 0.d0

        else

!          Left boundary.

           S(0) = phi_p(0) - v*rhoh*(1.d0-theta)*(phi_p(1) - phi_p(0))

           M(0,-1) = 0.d0
           M(0, 0) = 1.d0 - v*rhoh*theta
           M(0,+1) = + v*rhoh*theta

!          Right boundary.

           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S

!       FULL TIME STEP WITH EXPLICIT METHOD.

        M = 0.d0
        S = 0.d0

!       Source term.

        do i=1,Nx-1
           S(i) = phi_p(i) - 0.5d0*v*rho*(phi(i+1) - phi(i-1))
        end do

!       Matrix components.

        do i=1,Nx-1
           M(i,-1) = 0.d0
           M(i, 0) = 1.d0
           M(i,+1) = 0.d0
        end do

!       Boundaries using upwind.

        if (v>=0.d0) then

!          Left boundary.

           M(0,0) = 1.d0

!          Right boundary.

           S(Nx) = phi_p(Nx) - v*rho*(phi(Nx) - phi(Nx-1))

           M(Nx,-1) = 0.d0
           M(Nx, 0) = 1.d0
           M(Nx,+1) = 0.d0

        else

!          Left boundary.

           S(0) = phi_p(0) - v*rho*(phi(1) - phi(0))

           M(0,-1) = 0.d0
           M(0, 0) = 1.d0
           M(0,+1) = 0.d0

!          Right boundary.

           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S


!    ************************************
!    ***   IMPLICIT-EXPLICIT MINMOD   ***
!    ************************************

!    Here we first do half an implicit time step, and then
!    a full explicit time step calculating the sources with
!    the result of the implicit half time step.

     else if  (method=='IE-MM') then

!       HALF A TIME STEP WITH IMPLICIT METHOD.

        M = 0.d0
        S = 0.d0
        C = 0.d0

!       Positive speed.

        if (v>=0.d0) then

!          Calculate fluxes.

           flux(0) = 0.5d0*v*(phi_p(0) + phi_p(1))

           do i=1,Nx-1

!             Calculate right and left slopes.

              sloper = phi_p(i+1) - phi_p(i  )
              slopel = phi_p(i  ) - phi_p(i-1)

!             Find coefficients for flux.

              if (sloper*slopel>0.d0) then
                 if (dabs(sloper)<dabs(slopel)) then
                    C(i,-1) = 0.d0
                    C(i, 0) = 0.5d0
                    C(i,+1) = 0.5d0
                 else
                    C(i,-1) = - 0.5d0
                    C(i, 0) = + 1.5d0
                    C(i,+1) = 0.d0
                 end if
              else
                 C(i,-1) = 0.d0
                 C(i, 0) = 1.d0
                 C(i,+1) = 0.d0
              end if

!             Calculate explicit fluxes.

              flux(i) = v*(C(i,-1)*phi_p(i-1) + C(i,0)*phi_p(i) + C(i,1)*phi_p(i+1))

           end do

!       Negative speed.

        else

!          Calculate fluxes.

           flux(Nx-1) = 0.5d0*v*(phi_p(Nx) + phi_p(Nx-1))

           do i=0,Nx-2

!             Calculate right and left slopes.

              sloper = phi_p(i+2) - phi_p(i+1)
              slopel = phi_p(i+1) - phi_p(i  )

              if (sloper*slopel>0.d0) then
                 if (dabs(sloper)<dabs(slopel)) then
                    C(i,0) = 0.d0
                    C(i,1) = 1.5d0
                    C(i,2) = -0.5d0
                 else
                    C(i,0) = 0.5d0
                    C(i,1) = 0.5d0
                    C(i,2) = 0.d0
                 end if
              else
                 C(i,0) = 0.d0
                 C(i,1) = 1.d0
                 C(i,2) = 0.d0
              end if

!             Calculate explicit fluxes.

              flux(i) = v*(C(i,0)*phi_p(i) + C(i,1)*phi_p(i+1) + C(i,2)*phi_p(i+2))

           end do

        end if

!       Find source term.

        do i=1,Nx-1
           S(i) = phi_p(i) - rhoh*(1.d0-theta)*(flux(i) - flux(i-1))
        end do

!       Matrix components.

        aux = v*rhoh*theta

        do i=1,Nx-1
           M(i,-2) = - aux*C(i-1,-1)
           M(i,-1) = + aux*(C(i,-1) - C(i-1,0))
           M(i, 0) = + aux*(C(i, 0) - C(i-1,1)) + 1.d0
           M(i,+1) = + aux*(C(i, 1) - C(i-1,2))
           M(i,+2) = + aux*C(i,2)
        end do

!       Boundaries using upwind.

        if (v>=0.d0) then

!          Left boundary.

           S(0) = 0.d0
           M(0,0) = 1.d0

!          Right boundary.

           S(Nx) = phi_p(Nx) - v*rhoh*(1.d0-theta)*(phi_p(Nx) - phi_p(Nx-1))

           M(Nx,-1) = - v*rhoh*theta
           M(Nx, 0) = 1.d0 + v*rhoh*theta
           M(Nx,+1) = 0.d0

        else

!          Left boundary.

           S(0) = phi_p(0) - v*rhoh*(1.d0-theta)*(phi_p(1) - phi_p(0))

           M(0,-1) = 0.d0
           M(0, 0) = 1.d0 - v*rhoh*theta
           M(0,+1) = + v*rhoh*theta

!          Right boundary.

           S(Nx) = 0.d0
           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S

!       FULL TIME STEP WITH EXPLICIT METHOD.

        M = 0.d0
        S = 0.d0
        C = 0.d0

!       Positive speed.

        if (v>=0.d0) then

!          Calculate fluxes.

           flux(0) = 0.5d0*v*(phi(0) + phi(1))

           do i=1,Nx-1

!             Calculate right and left slopes.

              sloper = phi(i+1) - phi(i  )
              slopel = phi(i  ) - phi(i-1)

!             Find coefficients for flux.

              if (sloper*slopel>0.d0) then
                 if (dabs(sloper)<dabs(slopel)) then
                    C(i,-1) = 0.d0
                    C(i, 0) = 0.5d0
                    C(i,+1) = 0.5d0
                 else
                    C(i,-1) = - 0.5d0
                    C(i, 0) = + 1.5d0
                    C(i,+1) = 0.d0
                 end if
              else
                 C(i,-1) = 0.d0
                 C(i, 0) = 1.d0
                 C(i,+1) = 0.d0
              end if

!             Calculate explicit fluxes.

              flux(i) = v*(C(i,-1)*phi(i-1) + C(i,0)*phi(i) + C(i,1)*phi(i+1))

           end do

!       Negative speed.

        else

!          Calculate fluxes.

           flux(Nx-1) = 0.5d0*v*(phi(Nx) + phi(Nx-1))

           do i=0,Nx-2

!             Calculate right and left slopes.

              sloper = phi(i+2) - phi(i+1)
              slopel = phi(i+1) - phi(i  )

              if (sloper*slopel>0.d0) then
                 if (dabs(sloper)<dabs(slopel)) then
                    C(i,0) = 0.d0
                    C(i,1) = 1.5d0
                    C(i,2) = -0.5d0
                 else
                    C(i,0) = 0.5d0
                    C(i,1) = 0.5d0
                    C(i,2) = 0.d0
                 end if
              else
                 C(i,0) = 0.d0
                 C(i,1) = 1.d0
                 C(i,2) = 0.d0
              end if

!             Calculate explicit fluxes.

              flux(i) = v*(C(i,0)*phi(i) + C(i,1)*phi(i+1) + C(i,2)*phi(i+2))

           end do

        end if

!       Find source term.

        do i=1,Nx-1
           S(i) = phi_p(i) - rho*(flux(i) - flux(i-1))
        end do

!       Matrix components.

        aux = v*rho*theta

        do i=1,Nx-1
           M(i,-2) = 0.d0
           M(i,-1) = 0.d0
           M(i, 0) = 1.d0
           M(i,+1) = 0.d0
           M(i,+2) = 0.d0
        end do

!       Boundaries using upwind.

        if (v>=0.d0) then

!          Left boundary.

           M(0,0) = 1.d0

!          Right boundary.

           S(Nx) = phi_p(Nx) - v*rho*(phi(Nx) - phi(Nx-1))

           M(Nx,-1) = 0.d0
           M(Nx, 0) = 1.d0
           M(Nx,+1) = 0.d0

        else

!          Left boundary.

           S(0) = phi_p(0) - v*rho*(phi(1) - phi(0))

           M(0,-1) = 0.d0
           M(0, 0) = 1.d0
           M(0,+1) = 0.d0

!          Right boundary.

           M(Nx,0) = 1.d0

        end if

!       Invert matrix.

        call invert(Nx,M,S)

!       Copy solution.

        phi = S


!    **************************
!    ***   UNKNOWN METHOD   ***
!    **************************

     else

        print *
        print *, 'Unknown integration method.'
        print *, 'Aborting ...'
        print *
        stop

     end if


!    *****************************
!    ***   SAVE DATA TO FILE   ***
!    *****************************

     if (mod(l,Noutput).eq.0) then

!       Wave function.

        write(1,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(phi(i)) > 1.0D-50) then
              write(1,"(2ES16.8)") x(i),phi(i)
           else
              write(1,"(2ES16.8)") x(i),0.0D0
           end if 
        end do

        write(1,*)

!       Save total variation.

        aux = 0.d0

        do i=1,Nx
           aux = aux + abs(phi(i)-phi(i-1))
        end do

        if (abs(aux)>1.d-50) then
           write(2,"(2ES16.8)") t,aux
        else
           write(2,"(2ES16.8)") t,0.d0
        end if

     end if


!    ***********************************
!    ***   END MAIN EVOLUTION LOOP   ***
!    ***********************************

!    Time step information to screen.

     if (mod(l,Noutput).eq.0) then
        write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',l,'   | ',t,'  | '
     end if

  end do

  print *,'------------------------------'

  
! ******************************
! ***   CLOSE OUTPUT FILES   ***
! ******************************

  close(1)
  close(2)

  
! ***************
! ***   END   ***
! ***************

  print *
  print *, 'PROGRAM HAS FINISHED'
  print *
  print *, 'Have a nice day!'
  print *
  print *
  print *

  end program advimplicit








  subroutine invert(Nx,M,S)

! ******************
! ***   INVERT   ***
! ******************

! This routine just calls the matrix inversion routines
! routines from Numerical Recipes.

  implicit none

  integer Nx,NP

  real(8), dimension(0:Nx) :: S
  real(8), dimension (0:Nx,-2:2) :: M,Ml

  integer, dimension (0:Nx) :: indx

! Call matrix inversion routines.

  ML = 0.d0
  NP = Nx+1

  call bandec(M,NP,2,2,ML,indx)
  call banbks(M,NP,2,2,ML,indx,S)


! ***************
! ***   END   ***
! ***************

  end subroutine invert






  subroutine bandec(a,n,m1,m2,al,indx)

! This routine performs an LU decomposition of a
! band diagonal matrix.
!
! On input, a(n,m1+m2+1) contains the matrix elements
! defining the system.  Here "n" is the number of
! equations, and (m1,m2) indicate how many subdiagonals
! and superdiagonals are different from zero (e.g. for
! a tridiagonal matrix m1=m2=1, and for a pentadiagonal
! (m1=m2=2).
!
! On output, "a" now contains the upper triangular matrix
! and "al" the lower tringular matrix.  The integer vector "indx"
! records the row permutation effected by the partial pivoting
! and is required by subroutine "banbks".

  implicit none

  integer n,m1,m2
  integer i,j,k,l,mm

  integer indx(n)

  real(8) a(n,m1+m2+1),al(n,m1+m2+1)
  real(8) d,dum,TINY

  parameter (TINY=1.e-20)

  mm=m1+m2+1
  l=m1

  do i=1,m1
     do j=m1+2-i,mm
        a(i,j-l)=a(i,j)
     end do
     l=l-1
     do j=mm-l,mm
        a(i,j)=0.d0
     end do
  end do
 
  d=1.d0
  l=m1

  do k=1,n

    dum=a(k,1)
    i=k
 
    if (l.lt.n) l=l+1

    do j=k+1,l
       if (abs(a(j,1)).gt.abs(dum)) then
          dum=a(j,1)
          i=j
       end if
    end do

    indx(k)=i

    if (dum.eq.0.) a(k,1)=TINY

    if (i.ne.k) then
       d=-d
       do j=1,mm
          dum=a(k,j)
          a(k,j)=a(i,j)
          a(i,j)=dum
       end do
    end if

    do i=k+1,l
       dum=a(i,1)/a(k,1)
       al(k,i-k)=dum
       do j=2,mm
          a(i,j-1)=a(i,j)-dum*a(k,j)
       end do
       a(i,mm)=0.d0
    end do
 
  end do

  end subroutine bandec








  subroutine banbks(a,n,m1,m2,al,indx,b)

! This routine solves a band-diagonal system of equations.
! The subroutine must be called immediately after "bandec".
!
! On input, "a" contains the upper triangular matrix,
! and "al" the lower triangular.  Here "n" is the number of
! equations, and (m1,m2) indicate how many subdiagonals
! and superdiagonals are different from zero (e.g. for
! a tridiagonal matrix m1=m2=1, and for a pentadiagonal
! (m1=m2=2).
!
! Also on input "b" contains the right hand side source vector,
! and "indx" records the row permutation effected by the partial
! pivoting as obtained in "bandec".
!
! On output, the vector "b" contains the solution.

  implicit none

  integer n,m1,m2
  integer i,k,l,mm
  
  integer indx(n)

  real(8) b(n),a(n,m1+m2+1),al(n,m1+m2+1)
  real(8) dum  

  mm=m1+m2+1
  l=m1

  do k=1,n
     i=indx(k)
     if (i.ne.k) then
        dum=b(k)
        b(k)=b(i)
        b(i)=dum 
     end if
     if (l.lt.n) l=l+1
     do i=k+1,l
        b(i)=b(i)-al(k,i-k)*b(k)
     end do
  end do

  l=1

  do i=n,1,-1
     dum=b(i)
     do k=2,l
        dum=dum-a(i,k)*b(k+i-1)
     end do
     b(i)=dum/a(i,1)
     if (l.lt.mm) l=l+1
  end do

  end subroutine banbks




