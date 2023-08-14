!$Header: /usr/local/ollincvs/Codes/SimpleCodes/wave.f90,v 1.8 2023/06/01 18:33:08 malcubi Exp $

  program wave

! This is a simple code for the wave equation in first order form
! using ICN (Iterative Crank-Nicholson).  For this we define:
!
! pi  = dphi/dt
! psi = dphi/dx
!
! The equations to solve take the simple form:
!
! dphi/dt = pi
!
! dpsi/dt = dpi/dx
!
! dpi/dt  = v^2 dpsi/dx  
!
! where v is the wave speed.


! *************************
! ***   WAVE EQUATION   ***
! *************************

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
  real(8), allocatable, dimension(:) :: phi_p   ! Old phi.
  real(8), allocatable, dimension(:) :: sphi    ! Source for phi.

  real(8), allocatable, dimension(:) :: psi     ! dphi/dx.
  real(8), allocatable, dimension(:) :: psi_p   ! Old psi.
  real(8), allocatable, dimension(:) :: spsi    ! Source for psi.

  real(8), allocatable, dimension(:) :: pi      ! dphi/dt
  real(8), allocatable, dimension(:) :: pi_p    ! Old pi.
  real(8), allocatable, dimension(:) :: spi     ! Source for pi.


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
  print *, 'Give integration method (only ICN for now)'
  read(*,*) method


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

  allocate(x(0:Nx))

  allocate(phi(0:Nx),phi_p(0:Nx),sphi(0:Nx))
  allocate(psi(0:Nx),psi_p(0:Nx),spsi(0:Nx))
  allocate(pi(0:Nx),pi_p(0:Nx),spi(0:Nx))


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

! Parameters for initial data.

  a0 = 1.0D0
  x0 = dble(Nx/2)*dx
  s0 = 1.0D0

! Initial data (stationary gaussian).

  pi = 0.0D0

  do i=0,Nx
     phi(i) = a0*exp(-(x(i)-x0)**2/s0**2)
     psi(i) = - 2.0D0*(x(i)-x0)/s0**2*phi(i)
  end do

! Output to screen.

  write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t,'  | '


! *****************************
! ***   OPEN OUTPUT FILES   ***
! *****************************

  open(1,file='phi_w.xl',form='formatted',status='replace')
  open(2,file='psi_w.xl',form='formatted',status='replace')
  open(3,file='pi_w.xl',form='formatted',status='replace')
  open(4,file='wp_w.xl',form='formatted',status='replace')
  open(5,file='wm_w.xl',form='formatted',status='replace')


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

! Spatial derivative.
 
  write(2,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(psi(i)) > 1.0D-50) then
        write(2,"(2ES16.8)") x(i),psi(i)
     else
        write(2,"(2ES16.8)") x(i),0.0D0
     end if
  end do

! Time derivative.

  write(3,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(pi(i)) > 1.0D-50) then
        write(3,"(2ES16.8)") x(i),pi(i)
     else
        write(3,"(2ES16.8)") x(i),0.0D0
     end if
  end do

! Eigenfields.

  write(4,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(pi(i)) > 1.0D-50) then
        write(4,"(2ES16.8)") x(i),pi(i) - psi(i)
     else
        write(4,"(2ES16.8)") x(i),0.0D0
     end if
  end do

  write(5,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(pi(i)) > 1.0D-50) then
        write(5,"(2ES16.8)") x(i),pi(i) + psi(i)
     else
        write(5,"(2ES16.8)") x(i),0.0D0
     end if
  end do

! Leave blank spaces before next time level.

  write(1,*)
  write(2,*)
  write(3,*)
  write(4,*)
  write(5,*)


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do l=1,Nt

!    Time.

     t = t + dt

!    Save old time step.

     phi_p = phi
     psi_p = psi
     pi_p = pi

!    Iterative Crank-Nicholson (3 steps).  This uses centered
!    differences in space, but iterates three times in order
!    to get second order accuracy in time and also stability.
!    The first two iterations go only half a time step forward,
!    and the third one goes the full time step.

     if (method=='icn') then

!       Three ICN iterations.

        do j=1,3

!          Source for phi.

           sphi = pi

!          Sources for psi and pi using centered differences.

           do i=1,Nx-1
              spsi(i) = 0.5D0*(pi(i+1) - pi(i-1))/dx
              spi(i)  = v**2*0.5D0*(psi(i+1) - psi(i-1))/dx
           end do

!          Boundary.  For the boundary we use one-sided differences
!          for the source of psi (they can be first order or second order).
!          The source for pi is then obtained by asking for the incoming
!          wave to be zero.

!          spsi(0) = + v*(psi(1) - psi(0))/dx
           spsi(0) = + v*(2.0D0*(psi(1) - psi(0)) - 0.5D0*(psi(2) - psi(0)))/dx
           spi(0 ) = + spsi(0)

!          spsi(Nx) = - v*(psi(Nx) - psi(Nx-1))/dx
           spsi(Nx) = - v*(2.0D0*(psi(Nx) - psi(Nx-1)) - 0.5D0*(psi(Nx) - psi(Nx-2)))/dx
           spi(Nx ) = - spsi(Nx)

!          Update functions.  First two iterations we only
!          advance half a time step.  The third time we advance
!          the full step.

           if (j<3) then
              phi = phi_p + 0.5D0*dt*sphi
              psi = psi_p + 0.5D0*dt*spsi
              pi  = pi_p  + 0.5D0*dt*spi
           else
              phi = phi_p + dt*sphi
              psi = psi_p + dt*spsi
              pi  = pi_p  + dt*spi
           end if

        end do

!    Unknown.

     else

        print *, 'Unknown integration method.'
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

!       Spatial derivative.
 
        write(2,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(psi(i)) > 1.0D-50) then
              write(2,"(2ES16.8)") x(i),psi(i)
           else
              write(2,"(2ES16.8)") x(i),0.0D0
           end if
        end do

!       Time derivative.

        write(3,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(pi(i)) > 1.0D-50) then
              write(3,"(2ES16.8)") x(i),pi(i)
           else
              write(3,"(2ES16.8)") x(i),0.0D0
           end if
        end do

!       Eigenfields.

        write(4,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(pi(i)) > 1.0D-50) then
              write(4,"(2ES16.8)") x(i),pi(i) - psi(i)
           else
              write(4,"(2ES16.8)") x(i),0.0D0
           end if
        end do

        write(5,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(pi(i)) > 1.0D-50) then
              write(5,"(2ES16.8)") x(i),pi(i) + psi(i)
           else
              write(5,"(2ES16.8)") x(i),0.0D0
           end if
        end do

!       Leave blank spaces before next time level.

        write(1,*)
        write(2,*)
        write(3,*)
        write(4,*)
        write(5,*)

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
  close(5)


! ***************
! ***   END   ***
! ***************

  print *
  print *, 'PROGRAM HAS FINISHED'
  print *
  print *, 'Have a nice day!'
  print *
  print *

  end program wave

