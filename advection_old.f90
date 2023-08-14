!$Header: /usr/local/ollincvs/Codes/SimpleCodes/advection.f90,v 1.9 2016/10/04 19:32:00 malcubi Exp $

  program advection

! This is a simple code for the advection equation
! that can use several different numerical methods:
! euler, upwind, lax-wendroff and ICN.


! ******************************
! ***   ADVECTION EQUATION   ***
! ******************************

! Declare variables.

  implicit none

  integer i,j,l         ! Counters
  integer Nt            ! Total number of time steps.
  integer Nx            ! Total number of grid points.
  integer Noutput       ! Frequency of output.

  real(8) dx            ! Grid spacing.
  real(8) dt            ! Time step.
  real(8) t             ! Time.

  real(8) v             ! Wave speed.

  real(8) x0            ! Center of initial gaussian.
  real(8) s0            ! Width of initial gaussian.
  real(8) a0            ! Amplitude of initial gaussian.

  character(20) method  ! Integration method.

  real(8), allocatable, dimension(:) :: x       ! Position.

  real(8), allocatable, dimension(:) :: phi     ! Wave function.
  real(8), allocatable, dimension(:) :: phi_p   ! Old wave function.
  real(8), allocatable, dimension(:) :: sphi    ! Source for wave function.

  real(8), allocatable, dimension(:) :: res     ! Residual evaluator.


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
  print *, 'Give frequency of output'
  read(*,*) Noutput

  print *
  print *, 'Give integration method (euler,upwind,laxwendroff,icn)'
  read(*,*) method

  print *

  
! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

  allocate(x(0:Nx),phi(0:Nx),phi_p(0:Nx),sphi(0:Nx),res(0:Nx))


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

! Initial data (gaussian).

  a0 = 1.0D0
  x0 = dble(Nx/2)*dx
  s0 = 1.0D0

  do i=0,Nx
     phi(i) = a0*exp(-(x(i)-x0)**2/s0**2)
  end do

! Initialize residual.

  res = 0.0D0

! Output to screen.

  write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t,'  | '


! *****************************
! ***   OPEN OUTPUT FILES   ***
! *****************************

  open(1,file='phi.xl',form='formatted',status='replace')
  open(2,file='res.xl',form='formatted',status='replace')


! *********************************
! ***   SAVE THE INITIAL DATA   ***
! *********************************

! Wave function.

  write(1,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(phi(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") x(i),phi(i)
     else
        write(1,"(2ES16.8)") x(i),0.0D0
     end if
  end do

! Residual.
 
  write(2,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     write(2,"(2ES16.8)") x(i),0.0D0
  end do

! Leave blank space before next time level.

  write(1,*)
  write(2,*)


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do l=1,Nt

!    Time.

     t = t + dt

!    Save old time step.

     phi_p = phi

!    Euler method.  This uses forward differencing in time
!    and central differences in space.  It is unstable.

      if (method=='euler') then

!       Centered differences for interior points.

        do i=1,Nx-1
           sphi(i) = - 0.5D0*v*(phi(i+1) - phi(i-1))/dx
        end do

!       Boundaries with one-sided differences.

        if (v>0.0D0) then
           sphi(0 ) = 0.0D0
           sphi(Nx) = - v*(phi(Nx) - phi(Nx-1))/dx
        else
           sphi(0 ) = - v*(phi(1 ) - phi(0   ))/dx
           sphi(Nx) = 0.0D0
        end if

!       Update wave function.

        phi = phi_p + dt*sphi

!    Upwind.  This is first order in both time and space.
!    It is stable as long as the one-sided spatial
!    differences are done in the upwind direction.

     else if (method=='upwind') then

!       Interior points for positive speed.

        if (v>0.0D0) then

           do i=1,Nx-1
              sphi(i) = - v*(phi(i) - phi(i-1))/dx
           end do

!       Interior points for negative speed.

        else

           do i=1,Nx-1
              sphi(i) = - v*(phi(i+1) - phi(i))/dx
           end do

        end if

!       Boundaries with one-sided differences.

        if (v>0.0D0) then
           sphi(0 ) = 0.0D0
           sphi(Nx) = - v*(phi(Nx) - phi(Nx-1))/dx  
        else
           sphi(0 ) = - v*(phi(1 ) - phi(0   ))/dx
           sphi(Nx) = 0.0D0
        end if

!       Update wave function.

        phi = phi_p + dt*sphi

!    Lax-Wendroff.  This method is second order in both space
!    and time and is also stable.  It is based on doing a
!    Taylor expansion in time to second order and then
!    substituting the time derivatives for spatial derivatives
!    using the original differential equation.

     else if (method=='laxwendroff') then

!       Interior points using centered differences plus correction.

        do i=1,Nx-1
           sphi(i) = - 0.5D0*v*(phi(i+1) - phi(i-1))/dx &
              + 0.5D0*dt*v**2*(phi(i+1) - 2.0D0*phi(i) + phi(i-1))/dx**2
        end do

!       Boundaries with one-sided differences.

        if (v>0.0D0) then
           sphi(0 ) = 0.0D0
           sphi(Nx) = - v*(phi(Nx) - phi(Nx-1))/dx  
        else
           sphi(0 ) = - v*(phi(1 ) - phi(0   ))/dx
           sphi(Nx) = 0.0D0
        end if

!       Update wave function.

        phi = phi_p + dt*sphi

!    Iterative Crank-Nicholson (3 steps).  This uses centered
!    differences in space, but iterates three times in order
!    to get second order accuracy in time and also stability.
!    The first two iterations go only half a time step forward,
!    and the third one goes the full time step.

     else if (method=='icn') then

!       Three ICN iterations.

        do j=1,3

!          Centered differences.

           do i=1,Nx-1
              sphi(i) = - 0.5D0*v*(phi(i+1) - phi(i-1))/dx
           end do

!          Boundaries with one-sided differences.

           if (v>0.0D0) then
              sphi(0 ) = 0.0D0
              sphi(Nx) = - v*(phi(Nx) - phi(Nx-1))/dx  
           else
              sphi(0 ) = - v*(phi(1 ) - phi(0   ))/dx
              sphi(Nx) = 0.0D0
           end if

!          Update wave function.  First two iterations we only
!          advance half a time step.  The third time we advance
!          the full step.

           if (j<3) then
              phi = phi_p + 0.5D0*dt*sphi
           else
              phi = phi_p + dt*sphi
           end if

        end do

!    Unknown method.

     else

        print *
        print *, 'Unknown integration method.'
        print *
        stop

     end if


!    ******************************
!    ***   CALCULATE RESIDUAL   ***
!    ******************************

!    Independent residual evaluation.  Here we calculate the value
!    of a different finite difference approximation to the advection
!    equation for the purpose of testing convergence.  This is done
!    to second order.

     do i=1,Nx-1
        res(i) = 0.5D0*((phi(i+1) - phi_p(i+1)) + (phi(i-1) - phi_p(i-1)))/dt &
               + 0.25D0*v*((phi(i+1) - phi(i-1)) + (phi_p(i+1) - phi_p(i-1)))/dx
     end do


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

!       Residual.

        write(2,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(res(i)) > 1.0D-50) then
              write(2,"(2ES16.8)") x(i),res(i)
           else
              write(2,"(2ES16.8)") x(i),0.0D0
           end if 
        end do

!       Leave blank space before next time level.

        write(1,*)
        write(2,*)

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

  end program advection

