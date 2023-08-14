!$Header: /usr/local/ollincvs/Codes/SimpleCodes/toygammadriver.f90,v 1.5 2020/02/17 16:56:15 malcubi Exp $

  program toygammadriver

! This is a simple code to test a toy model
! for the Gammadriver shift in 1+1.

! Declare variables.

  implicit none

  integer i,j,l         ! Counters.
  integer Nt            ! Total number of time steps.
  integer Nx            ! Total number of grid points.
  integer Noutput       ! Frequency of output.

  real(8) dx            ! Grid spacing.
  real(8) dt            ! Time step.
  real(8) t             ! Time.

  real(8) csi           ! Gammadriver parameter.
  real(8) sigma         ! Advection term coefficient.

  real(8) x0            ! Center of initial gaussian.
  real(8) s0            ! Width of initial gaussian.
  real(8) a0            ! Amplitude of initial gaussian.

  character(20) method  ! Integration method.

  real(8), allocatable, dimension(:) :: x         ! Position.

  real(8), allocatable, dimension(:) :: g         ! Metric coefficient.
  real(8), allocatable, dimension(:) :: g_p       ! Old g.
  real(8), allocatable, dimension(:) :: sg        ! Source for g.

  real(8), allocatable, dimension(:) :: beta      ! Shift.
  real(8), allocatable, dimension(:) :: beta_p    ! Old beta.
  real(8), allocatable, dimension(:) :: sbeta     ! Source for beta.

  real(8), allocatable, dimension(:) :: lambdap   ! + characteristic speed.
  real(8), allocatable, dimension(:) :: lambdam   ! - characteristic speed.

  real(8), allocatable, dimension(:) :: wp        ! Right propagatung eigenfield.
  real(8), allocatable, dimension(:) :: wm        ! Left  propagatung eigenfield.


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
  print *, 'Give total number of grid points Nx'
  read(*,*) Nx

  print *
  print *, 'Give total number of time steps'
  read(*,*) Nt

  print *
  print *, 'Give frequency of output'
  read(*,*) Noutput

  print *
  print *, 'Give Gammadriver parameter'
  read(*,*) csi

  print *
  print *, 'Give advection term coefficient'
  read(*,*) sigma

  print *


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

  allocate(x(0:Nx))

  allocate(g(0:Nx),g_p(0:Nx),sg(0:Nx))
  allocate(beta(0:Nx),beta_p(0:Nx),sbeta(0:Nx))

  allocate(lambdap(0:Nx),lambdam(0:Nx))
  allocate(wp(0:Nx),wm(0:Nx))


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

  t = 0.d0

! Parameters for initial data.

  a0 = 1.d0
  x0 = dble(Nx/2)*dx
  s0 = 1.d0

! Initial data (gaussian in g).

  beta = 0.d0

  do i=0,Nx
     g(i) = 1.d0 + a0*exp(-(x(i)-x0)**2/s0**2)
  end do

! Initial characteristic speeds.

  lambdap = - beta + sqrt(csi/g)
  lambdam = - beta - sqrt(csi/g)

! Output to screen.

  write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t,'  | '


! *****************************
! ***   OPEN OUTPUT FILES   ***
! *****************************

  open(1,file='g.xl',form='formatted',status='replace')
  open(2,file='beta.xl',form='formatted',status='replace')
  open(3,file='lambdap.xl',form='formatted',status='replace')
  open(4,file='lambdam.xl',form='formatted',status='replace')


! *********************************
! ***   SAVE THE INITIAL DATA   ***
! *********************************

! Metric function.

  write(1,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(g(i)) > 1.d-50) then
        write(1,"(2ES16.8)") x(i),g(i)
     else
        write(1,"(2ES16.8)") x(i),0.d0
     end if
  end do

! Shift.

  write(2,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(beta(i)) > 1.d-50) then
        write(2,"(2ES16.8)") x(i),beta(i)
     else
        write(2,"(2ES16.8)") x(i),0.d0
     end if
  end do

! Positive characteristic speed.

  write(3,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(lambdap(i)) > 1.d-50) then
        write(3,"(2ES16.8)") x(i),lambdap(i)
     else
        write(3,"(2ES16.8)") x(i),0.d0
     end if
  end do

! Negative characteristic speed.

  write(4,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(lambdam(i)) > 1.d-50) then
        write(4,"(2ES16.8)") x(i),lambdam(i)
     else
        write(4,"(2ES16.8)") x(i),0.d0
     end if
  end do


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

! For the moment always use ICN method.

  method = 'icn'

! Start iterations.

  do l=1,Nt

!    Time.

     t = t + dt

!    Save old time step.

     g_p    = g
     beta_p = beta

!    Iterative Crank-Nicholson (3 steps).

     if (method=='icn') then

!       ICN iterations.

        do j=1,3

!          Sources for g:
!
!          dg/dt  =  beta dg/dx  +  2 g dbeta/dx

           do i=1,Nx-1
              sg(i) = 0.5d0*beta(i)*(g(i+1) - g(i-1))/dx &
                    + g(i)*(beta(i+1) - beta(i-1))/dx 
           end do

!          Sources for beta:
!
!          dbeta/dt  =  (csi / 2 g^2) dg/dx  +  sigma beta dbeta/dx

           do i=1,Nx-1
              sbeta(i) = 0.25d0*csi/g(i)**2*(g(i+1) - g(i-1))/dx &
                       + 0.5d0*sigma*beta(i)*(beta(i+1) - beta(i-1))/dx
           end do

!          Boundaries: simple first order outgoing waves.

           sg(0 ) = (beta(0 ) + sqrt(csi/g(0 )))*(g(1 ) - g(0   ))/dx
           sg(Nx) = (beta(Nx) - sqrt(csi/g(Nx)))*(g(Nx) - g(Nx-1))/dx

           sbeta(0 ) = (beta(0 ) + sqrt(csi/g(0 )))*(beta(1 ) - beta(0   ))/dx
           sbeta(Nx) = (beta(Nx) - sqrt(csi/g(Nx)))*(beta(Nx) - beta(Nx-1))/dx

!          Update functions.  First two iterations we only
!          advance half a time step.  The third time we advance
!          the full step.

           if (j<3) then
              g = g_p + 0.5d0*dt*sg
              beta = beta_p + 0.5d0*dt*sbeta
           else
              g = g_p + dt*sg
              beta = beta_p + dt*sbeta
           end if

        end do

!    Unknown.

     else

        print *, 'Unknown integration method.'
        print *
        stop

     end if


!    **************************************
!    ***   FIND CHARACTERISTIC SPEEDS   ***
!    **************************************

     lambdap = - beta + sqrt(csi/g)
     lambdam = - beta - sqrt(csi/g)


!    *****************************
!    ***   SAVE DATA TO FILE   ***
!    *****************************

     if (mod(l,Noutput).eq.0) then

!       Metric function.

        write(1,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(g(i)) > 1.d-50) then
              write(1,"(2ES16.8)") x(i),g(i)
           else
              write(1,"(2ES16.8)") x(i),0.d0
           end if
        end do

!       Shift.

        write(2,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(beta(i)) > 1.d-50) then
              write(2,"(2ES16.8)") x(i),beta(i)
           else
              write(2,"(2ES16.8)") x(i),0.d0
           end if
        end do

!       Positive characteristic speed.

        write(3,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(lambdap(i)) > 1.d-50) then
              write(3,"(2ES16.8)") x(i),lambdap(i)
           else
              write(3,"(2ES16.8)") x(i),0.d0
           end if
        end do

!       Negative characteristic speed.

        write(4,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(lambdam(i)) > 1.d-50) then
              write(4,"(2ES16.8)") x(i),lambdam(i)
           else
              write(4,"(2ES16.8)") x(i),0.d0
           end if
        end do

!       Leave blank spaces before next time level.

        write(1,*)
        write(2,*)
        write(3,*)
        write(4,*)

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
  close(3)
  close(4)


! ***************
! ***   END   ***
! ***************

  print *
  print *, 'PROGRAM HAS FINISHED'
  print *
  print *, 'Have a nice day!'
  print *
  print *

  end program toygammadriver


