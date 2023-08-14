!$Header: /usr/local/ollincvs/Codes/SimpleCodes/fluid.f90,v 1.15 2010/04/10 05:14:53 malcubi Exp $

  module arrays

! *****************************
! ***   MODULE FOR ARRAYS   ***
! *****************************

! This module defines the different arrays.
!
! The primitive fluid variables are:
!
! fluid_rho0      Rest-mass energy density in comoving frame.
! fluid_p         Pressure in comoving frame.
! fluid_v         Velocity of the fluid elements in laboratory frame.
! fluid_e         Specific internal energy in comoving frame.
! fluid_h         Specific enthalpy:  h = 1 + e + p/rho0
!
! The conserved variables are:
! 
! fluid_cD        Rest-mass energy density in laboratory frame:  D = rho0 W
! fluid_cE        Energy density (minus rest mass) in laboratory frame:  rho0 h W**2 - p - D
! fluid_cS        Momentum density in laboratory frame:  rho0 h W**2 v
!
! where W is the Lorentz factor given by:
!
! W  =  1 / sqrt( 1 - v**2 )
!
! Finally, when we use artificial viscosity we have:
!
! fluid_q         Contribution to the pressure from artificial viscosity.

  implicit none

  real(8), allocatable, dimension(:) :: x
  real(8), allocatable, dimension(:) :: fluid_rho0
  real(8), allocatable, dimension(:) :: fluid_e
  real(8), allocatable, dimension(:) :: fluid_v
  real(8), allocatable, dimension(:) :: fluid_p
  real(8), allocatable, dimension(:) :: fluid_h
  real(8), allocatable, dimension(:) :: fluid_cD,fluid_cD_p,sfluid_cD
  real(8), allocatable, dimension(:) :: fluid_cE,fluid_cE_p,sfluid_cE
  real(8), allocatable, dimension(:) :: fluid_cS,fluid_cS_p,sfluid_cS
  real(8), allocatable, dimension(:) :: fluid_q
  real(8), allocatable, dimension(:) :: fluid_vs
  real(8), allocatable, dimension(:) :: fluid_rho0l,fluid_rho0r
  real(8), allocatable, dimension(:) :: fluid_el,fluid_er
  real(8), allocatable, dimension(:) :: fluid_vl,fluid_vr
  real(8), allocatable, dimension(:) :: fluid_pl,fluid_pr


! ***************
! ***   END   ***
! ***************

  end module arrays








  program fluid

! This is a simple code for the relativistic Euler equations
! in Minkowski spacetime.


! ****************************************
! ***   RELATIVISTIC EULER EQUATIONS   ***
! ****************************************

! Declare variables.

  use arrays

  implicit none

  integer i,j,l         ! Counters
  integer Nt            ! Total number of time steps.
  integer Nx            ! Total number of grid points.
  integer Noutput       ! Frequency of output.

  real(8) dx            ! Grid spacing.
  real(8) dt            ! Time step.
  real(8) dtfac         ! Courant parameter.
  real(8) t             ! Time.

  real(8) gamma         ! Adiabatic index.
  real(8) q1,q2         ! Artificial viscosity parameters.

  real(8) rho0_left,rho0_right,p_left,p_right  ! Initial data parameters.

  character(20) initial ! Type of initial data.
  character(20) method  ! Integration method.


! **************************
! ***   GET PARAMETERS   ***
! **************************

! Grid parameters.

  print *
  print *, 'Give grid spacing dx'
  read(*,*) dx

  print *
  print *, 'Give Courant parameter dtfac'
  read(*,*) dtfac

  dt = dtfac*dx

  print *
  print *, 'Give total number of grid points Nx'
  read(*,*) Nx

  print *
  print *, 'Give total number of time steps'
  read(*,*) Nt

! Output.

  print *
  print *, 'Give frequency of output'
  read(*,*) Noutput

! Integration method.

  print *
  print *, 'Give integration method'
  read(*,*) method

! Artificial viscosity parameters.

  print *
  print *, 'Give artificial viscosity parameters: q1,q2'
  read(*,*) q1,q2

! Type of initial data.

  print *
  print *, 'Give initial data type'
  read(*,*) initial

! Shocktube parameters.

  if (initial=="shocktube") then

     print *
     print *, 'Give left density and pressure (rho0,p)'
     read(*,*) rho0_left,p_left

     print *
     print *, 'Give right density and pressure (rho0,p)'
     read(*,*) rho0_right,p_right

! Unknown initial data.

  else

     print *
     print *, 'Unknown type of initial data, aborting ...'
     print *
     stop

  end if 

  print *


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

  allocate(x(0:Nx))

  allocate(fluid_rho0(0:Nx),fluid_p(0:Nx),fluid_v(0:Nx))
  allocate(fluid_e(0:Nx),fluid_h(0:Nx))

  allocate(fluid_cD(0:Nx),fluid_cD_p(0:Nx),sfluid_cD(0:Nx))
  allocate(fluid_cE(0:Nx),fluid_cE_p(0:Nx),sfluid_cE(0:Nx))
  allocate(fluid_cS(0:Nx),fluid_cS_p(0:Nx),sfluid_cS(0:Nx))

  allocate(fluid_q(0:Nx),fluid_vs(0:Nx))

  allocate(fluid_rho0l(0:Nx),fluid_rho0r(0:Nx))
  allocate(fluid_el(0:Nx),fluid_er(0:Nx))
  allocate(fluid_vl(0:Nx),fluid_vr(0:Nx))
  allocate(fluid_pl(0:Nx),fluid_pr(0:Nx))


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

! Choose adiabatic index.

  gamma = 5.0/3.0

! Initialize time.

  t = 0.0D0

! Shock tube initial data.

  if (initial=="shocktube") then

!    Set fluid velocity to zero.

     fluid_v = 0.0

!    Left-side initial data.

     do i=0,Nx/2
        fluid_rho0(i) = rho0_left
        fluid_p(i) = p_left
     end do

!    Right-side initial data.

     do i=Nx/2,Nx
        fluid_rho0(i) = rho0_right
        fluid_p(i) = p_right
     end do

!    Calculate internal energy and enthalpy.

     fluid_e = fluid_p/(gamma-1.0)/fluid_rho0
     fluid_h = 1.0 + fluid_e + fluid_p/fluid_rho0

!    Find speed of sound.

     fluid_vs = sqrt(gamma*fluid_p/(fluid_rho0*fluid_h))

!    Find conserved quantities.

     fluid_cD = fluid_rho0
     fluid_cE = fluid_rho0*fluid_e
     fluid_cS = 0.0

  end if

! Initialize artificial viscosity to 0.

  fluid_q = 0.0


! *****************************
! ***   OPEN OUTPUT FILES   ***
! *****************************

   open(1,file='fluid_rho0.xl',form='formatted',status='replace')
   open(2,file='fluid_p.xl',form='formatted',status='replace')
   open(3,file='fluid_v.xl',form='formatted',status='replace')


! *********************************
! ***   SAVE THE INITIAL DATA   ***
! *********************************

! Rest-mass density.

  write(1,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(fluid_rho0(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") x(i),fluid_rho0(i)
     else
        write(1,"(2ES16.8)") x(i),0.0D0
     end if
  end do

! Pressure.

  write(2,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(fluid_p(i)) > 1.0D-50) then
        write(2,"(2ES16.8)") x(i),fluid_p(i)
     else
        write(2,"(2ES16.8)") x(i),0.0D0
     end if
  end do

! Fluid velocity.

  write(3,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(fluid_v(i)) > 1.0D-50) then
        write(3,"(2ES16.8)") x(i),fluid_v(i)
     else
        write(3,"(2ES16.8)") x(i),0.0D0
     end if
  end do

! Leave blank spaces before next time level.

  write(1,*)
  write(2,*)
  write(3,*)


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do l=1,Nt

!    Time.

     t = t + dt

!    Save old time step.

     fluid_cD_p = fluid_cD
     fluid_cE_p = fluid_cE
     fluid_cS_p = fluid_cS

!    Start Iterative Crank-Nicholson (ICN).  This iterates
!    three times in order to get second order accuracy in time
!    and also stability. The first two iterations go only half
!    a time step forward, and the third one goes the full time step.

     do j=1,3

!       Calculate sources.

        call sources(method,Nx,dx,dt,gamma)

!       Update functions.  First two iterations we only
!       advance half a time step.  The third time we advance
!       the full step.

        if (j<3) then
           fluid_cD = fluid_cD_p + 0.5D0*dt*sfluid_cD
           fluid_cE = fluid_cE_p + 0.5D0*dt*sfluid_cE
           fluid_cS = fluid_cS_p + 0.5D0*dt*sfluid_cS
        else
           fluid_cD = fluid_cD_p + dt*sfluid_cD
           fluid_cE = fluid_cE_p + dt*sfluid_cE
           fluid_cS = fluid_cS_p + dt*sfluid_cS
        end if

     end do


!    ***************************************
!    ***   RECOVER PRIMITIVE VARIABLES   ***
!    ***************************************

     call primitive(Nx,dx,gamma,q1,q2)


!    *****************************
!    ***   SAVE DATA TO FILE   ***
!    *****************************

     if (mod(l,Noutput).eq.0) then

!       Rest-mass density.

        write(1,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(fluid_rho0(i)) > 1.0D-50) then
              write(1,"(2ES16.8)") x(i),fluid_rho0(i)
           else
              write(1,"(2ES16.8)") x(i),0.0D0
           end if
        end do

!       Pressure.

        write(2,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(fluid_p(i)) > 1.0D-50) then
              write(2,"(2ES16.8)") x(i),fluid_p(i)
           else
              write(2,"(2ES16.8)") x(i),0.0D0
           end if
        end do

!       Fluid velocity.

        write(3,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(fluid_v(i)) > 1.0D-50) then
              write(3,"(2ES16.8)") x(i),fluid_v(i)
           else
              write(3,"(2ES16.8)") x(i),0.0D0
           end if
        end do

!       Leave blank spaces before next time level.

        write(1,*)
        write(2,*)
        write(3,*)

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


! *****************************
! ***   DEALLOCATE ARRAYS   ***
! *****************************

  deallocate(x)
  deallocate(fluid_rho0,fluid_p,fluid_v,fluid_e,fluid_h)
  deallocate(fluid_cD,fluid_cD_p,sfluid_cD)
  deallocate(fluid_cE,fluid_cE_p,sfluid_cE)
  deallocate(fluid_cS,fluid_cS_p,sfluid_cS)
  deallocate(fluid_q,fluid_vs)
  deallocate(fluid_rho0l,fluid_rho0r,fluid_el,fluid_er,fluid_vl,fluid_vr,fluid_pl,fluid_pr)


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

  end program fluid








  subroutine sources(method,Nx,dx,dt,gamma)

! ***************************************
! ***   SOURCES FOR EULER EQUATIONS   ***
! ***************************************

! This routine calculates the sources for the relativistic
! Euler equations written in conservative form. These equations
! take the general form:
!
! d D  =  - d  [ v D ]
!  t         x
!
! d E  =  - d  [ v (E + p) ]
!  t         x
!
! d S  =  - d  [ v S + p ]
!  t         x
!
! Notice that when we use artificial viscosity we need
! to add an extra contribution to the pressure:
!
! p  ->  p + q
!
! The routine can use several different differencing methods
! for the flux terms.
!
! We use a flux conservative formulation, so that the
! sources are calculated as:
!
! source = ( F      -  F      ) / dx
!             i+1/2     i-1/2
!
! Notice that for notation we in fact use:
!
! F      =  F(i-1)
!  i-1/2
!
! F      =  F(i)
!  i+1/2

! Declare variables.

  use arrays

  implicit none

  character(20) method

  integer i,Nx

  real(8) dx,dt,idx,gamma
  real(8) slope1,slope2,slopelim
  real(8) Wl,Wr,hl,hr
  real(8) vsl,vsr,vpl,vpr,vml,vmr,vpp,vmm
  real(8) Dl,Dr,El,Er,Sl,Sr
  real(8) aux,auxl,auxr

  real(8) :: flux_D(-1:Nx),flux_E(-1:Nx),flux_S(-1:Nx)


! *******************
! ***   NUMBERS   ***
! *******************

  idx = 1.0/dx


! ***************************
! ***   METHOD = CENTER   ***
! ***************************

! The fluxes at cell interfaces at just averaged over
! adjacent grid points.

  if (method=="center") then

     do i=0,Nx-1

!       Source for D.

        flux_D(i) = 0.5*(fluid_v(i)*fluid_cD(i) + fluid_v(i+1)*fluid_cD(i+1))

!       Source for E.

        flux_E(i) = 0.5*(fluid_v(i  )*(fluid_cE(i  ) + fluid_p(i  ) + fluid_q(i  )) &
                       + fluid_v(i+1)*(fluid_cE(i+1) + fluid_p(i+1) + fluid_q(i+1))) 

!       Source for S.

        flux_S(i) = 0.5*(fluid_v(i  )*fluid_cS(i  ) + fluid_p(i  ) + fluid_q(i  ) &
                       + fluid_v(i+1)*fluid_cS(i+1) + fluid_p(i+1) + fluid_q(i+1))

     end do


! ****************************
! ***   METHOD = UPWIND1   ***
! ****************************

! The fluxes at cell interfaces are just copied
! from the adjacent grid point depending on the
! sign of the fluid speed.

  else if (method=="upwind1") then

     do i=0,Nx-1

        aux = 0.5*(fluid_v(i) + fluid_v(i+1))

!       Source for D.

        flux_D(i) = (fluid_cD(i  )*(0.5 + sign(0.5D0,aux)) &
                  +  fluid_cD(i+1)*(0.5 - sign(0.5D0,aux)))*aux

!       Source for E.

        flux_E(i) = (fluid_cE(i  )*(0.5 + sign(0.5D0,aux)) &
                  +  fluid_cE(i+1)*(0.5 - sign(0.5D0,aux)))*aux &
                  + 0.5*(fluid_v(i)*(fluid_p(i) + fluid_q(i)) &
                  + fluid_v(i+1)*(fluid_p(i+1) + fluid_q(i+1)))

!       Source for S.

        flux_S(i) = (fluid_cS(i  )*(0.5 + sign(0.5D0,aux)) &
                  +  fluid_cS(i+1)*(0.5 - sign(0.5D0,aux)))*aux &
                  + 0.5*(fluid_p(i) + fluid_q(i) + fluid_p(i+1) + fluid_q(i+1))

     end do


! ****************************
! ***   METHOD = UPWIND2   ***
! ****************************

! The fluxes at cell interfaces are either interpolated
! (averaged) or extrapolated from one side depending on the
! sign of the fluid speed.

  else if (method=="upwind2") then

     do i=0,Nx-1

        aux = 0.5*(fluid_v(i) + fluid_v(i+1))

!       Source for D.

        if (aux>=0.0) then
           flux_D(i) = aux*(1.5*fluid_cD(i  ) - 0.5*fluid_cD(i-1))
        else
           flux_D(i) = aux*(1.5*fluid_cD(i+1) - 0.5*fluid_cD(i+2))
        end if

!       Source for E.

        if (aux>=0.0) then
           flux_E(i) = aux*(1.5*fluid_cE(i  ) - 0.5*fluid_cE(i-1))
        else
           flux_E(i) = aux*(1.5*fluid_cE(i+1) - 0.5*fluid_cE(i+2))
        end if

        flux_E(i) = flux_E(i) + 0.5*(fluid_v(i)*(fluid_p(i) + fluid_q(i)) &
                  + fluid_v(i+1)*(fluid_p(i+1) + fluid_q(i+1)))

!       Source for S.

        if (aux>=0.0) then
           flux_S(i) = aux*(1.5*fluid_cS(i  ) - 0.5*fluid_cS(i-1))
        else
           flux_S(i) = aux*(1.5*fluid_cS(i+1) - 0.5*fluid_cS(i+2))
        end if

        flux_S(i) = flux_S(i) + 0.5*(fluid_p(i) + fluid_q(i) + fluid_p(i+1) + fluid_q(i+1))

     end do


! *****************************
! ***   METHOD = LIMITER1   ***
! *****************************

! Mimmod limited second-order differences depending
! on the sign of the fluid velocity.

  else if (method=="limiter1") then

     do i=0,Nx-1

        aux = 0.5*(fluid_v(i) + fluid_v(i+1))

!       Source for D.

        if (aux>=0.0) then

           slope1 = fluid_cD(i+1) - fluid_cD(i  )
           slope2 = fluid_cD(i  ) - fluid_cD(i-1)

           if (slope1*slope2>0.0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.0
           end if

           flux_D(i) = aux*(fluid_cD(i  ) + 0.5*(1.0-aux*dt/dx)*slopelim)

        else

           slope1 = fluid_cD(i+1) - fluid_cD(i  )
           slope2 = fluid_cD(i+2) - fluid_cD(i+1)

           if (slope1*slope2>0.0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.0
           end if

           flux_D(i) = aux*(fluid_cD(i+1) - 0.5*(1.0+aux*dt/dx)*slopelim)

        end if

!    Source for E.

        if (aux>=0.0) then

           slope1 = fluid_cE(i+1) + fluid_p(i+1) + fluid_q(i+1) &
                  - fluid_cE(i  ) - fluid_p(i  ) - fluid_q(i  )
           slope2 = fluid_cE(i  ) + fluid_p(i  ) + fluid_q(i  ) &
                  - fluid_cE(i-1) - fluid_p(i-1) - fluid_q(i-1)

           if (slope1*slope2>0.0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.0
           end if

           flux_E(i) = aux*(fluid_cE(i  ) + fluid_p(i  ) + fluid_q(i  ) + 0.5*(1.0-aux*dt/dx)*slopelim)

        else

           slope1 = fluid_cE(i+1) + fluid_p(i+1) + fluid_q(i+1) &
                  - fluid_cE(i  ) - fluid_p(i  ) - fluid_q(i  )
           slope2 = fluid_cE(i+2) + fluid_p(i+2) + fluid_q(i+2) &
                  - fluid_cE(i+1) - fluid_p(i+1) - fluid_q(i+1)

           if (slope1*slope2>0.0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.0
           end if

           flux_E(i) = aux*(fluid_cE(i+1) + fluid_p(i+1) + fluid_q(i+1) - 0.5*(1.0+aux*dt/dx)*slopelim)

        end if

!       Source for S.

        if (aux>=0.0) then

           slope1 = fluid_cS(i+1) - fluid_cS(i  )
           slope2 = fluid_cS(i  ) - fluid_cS(i-1)

           if (slope1*slope2>0.0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.0
           end if

           flux_S(i) = aux*(fluid_cS(i  ) + 0.5*(1.0-aux*dt/dx)*slopelim)

        else

           slope1 = fluid_cS(i+1) - fluid_cS(i  )
           slope2 = fluid_cS(i+2) - fluid_cS(i+1)

           if (slope1*slope2>0.0) then
              if (abs(slope1)<abs(slope2)) then
                 slopelim = slope1
              else
                 slopelim = slope2
              end if
           else
              slopelim = 0.0
           end if

           flux_S(i) = aux*(fluid_cS(i+1) - 0.5*(1.0+aux*dt/dx)*slopelim)

        end if

        flux_S(i) = flux_S(i) + 0.5*(fluid_p(i) + fluid_q(i) + fluid_p(i+1) + fluid_q(i+1))

     end do


! *****************************
! ***   METHOD = LIMITER2   ***
! *****************************

! VanLeer limited second-order differences depending
! on the sign of the fluid velocity.

  else if (method=="limiter2") then

     do i=0,Nx-1

        aux = 0.5*(fluid_v(i) + fluid_v(i+1))

!       Source for D.

        if (aux>=0.0) then

           slope1 = fluid_cD(i+1) - fluid_cD(i  )
           slope2 = fluid_cD(i  ) - fluid_cD(i-1)

           if (slope1*slope2>0.0) then
              slopelim = 2.0*slope1*slope2/(slope1 + slope2)
           else
              slopelim = 0.0
           end if

           flux_D(i) = aux*(fluid_cD(i  ) + 0.5*(1.0-aux*dt/dx)*slopelim)

        else

           slope1 = fluid_cD(i+1) - fluid_cD(i  )
           slope2 = fluid_cD(i+2) - fluid_cD(i+1)

           if (slope1*slope2>0.0) then
              slopelim = 2.0*slope1*slope2/(slope1 + slope2)
           else
              slopelim = 0.0
           end if

           flux_D(i) = aux*(fluid_cD(i+1) - 0.5*(1.0+aux*dt/dx)*slopelim)

        end if

!       Source for E.

        if (aux>=0.0) then

           slope1 = fluid_cE(i+1) + fluid_p(i+1) + fluid_q(i+1) &
                  - fluid_cE(i  ) - fluid_p(i  ) - fluid_q(i  )
           slope2 = fluid_cE(i  ) + fluid_p(i  ) + fluid_q(i  ) &
                  - fluid_cE(i-1) - fluid_p(i-1) - fluid_q(i-1)

           if (slope1*slope2>0.0) then
              slopelim = 2.0*slope1*slope2/(slope1 + slope2)
           else
              slopelim = 0.0
           end if

           flux_E(i) = aux*(fluid_cE(i  ) + fluid_p(i  ) + fluid_q(i  ) + 0.5*(1.0-aux*dt/dx)*slopelim)

        else

           slope1 = fluid_cE(i+1) + fluid_p(i+1) + fluid_q(i+1) &
                  - fluid_cE(i  ) - fluid_p(i  ) - fluid_q(i  )
           slope2 = fluid_cE(i+2) + fluid_p(i+2) + fluid_q(i+2) &
                  - fluid_cE(i+1) - fluid_p(i+1) - fluid_q(i+1)

           if (slope1*slope2>0.0) then
              slopelim = 2.0*slope1*slope2/(slope1 + slope2)
           else
              slopelim = 0.0
           end if

           flux_E(i) = aux*(fluid_cE(i+1) + fluid_p(i+1) + fluid_q(i+1) - 0.5*(1.0+aux*dt/dx)*slopelim)

        end if

!       Source for S.

        if (aux>=0.0) then

           slope1 = fluid_cS(i+1) - fluid_cS(i  )
           slope2 = fluid_cS(i  ) - fluid_cS(i-1)

           if (slope1*slope2>0.0) then
              slopelim = 2.0*slope1*slope2/(slope1 + slope2)
           else
              slopelim = 0.0
           end if

           flux_S(i) = aux*(fluid_cS(i  ) + 0.5*(1.0-aux*dt/dx)*slopelim)

        else

           slope1 = fluid_cS(i+1) - fluid_cS(i  )
           slope2 = fluid_cS(i+2) - fluid_cS(i+1)

           if (slope1*slope2>0.0) then
              slopelim = 2.0*slope1*slope2/(slope1 + slope2)
           else
              slopelim = 0.0
           end if

           flux_S(i) = aux*(fluid_cS(i+1) - 0.5*(1.0+aux*dt/dx)*slopelim)

        end if

        flux_S(i) = flux_S(i) + 0.5*(fluid_p(i) + fluid_q(i) + fluid_p(i+1) + fluid_q(i+1))

     end do


! **************************
! ***   METHOD = HLLE1   ***
! **************************

! The HLLE method requires knowledge of the
! speed of sound, as well as the left and right
! characteritic speeds.

  else if (method=="hlle1") then

!    Calculate left and right characteristic speeds
!    using the relativistic composition of v and vs.

     fluid_vl = (fluid_v - fluid_vs)/(1.0 - fluid_v*fluid_vs)
     fluid_vr = (fluid_v + fluid_vs)/(1.0 + fluid_v*fluid_vs)

     do i=0,Nx-1

        auxl = min(0.0,fluid_vl(i),fluid_vl(i+1))
        auxr = max(0.0,fluid_vr(i),fluid_vr(i+1))

!       Source for D.

        flux_D(i) = (auxr*fluid_v(i  )*fluid_cD(i  ) &
                  -  auxl*fluid_v(i+1)*fluid_cD(i+1) &
                  +  auxr*auxl*(fluid_cD(i+1) - fluid_cD(i))) &
                  / (auxr - auxl)

!       Source for E.

        flux_E(i) = (auxr*fluid_v(i  )*(fluid_cE(i  ) + fluid_p(i  ) + fluid_q(i  )) &
                  -  auxl*fluid_v(i+1)*(fluid_cE(i+1) + fluid_p(i+1) + fluid_q(i+1)) &
                  +  auxr*auxl*(fluid_cE(i+1) - fluid_cE(i))) &
                  / (auxr - auxl)

!       Source for S.

        flux_S(i) = (auxr*(fluid_v(i  )*fluid_cS(i  ) + fluid_p(i  ) + fluid_q(i  )) &
                  -  auxl*(fluid_v(i+1)*fluid_cS(i+1) + fluid_p(i+1) + fluid_q(i+1)) &
                  +  auxr*auxl*(fluid_cS(i+1) - fluid_cS(i))) &
                  / (auxr - auxl)

     end do


! **************************
! ***   METHOD = HLLE2   ***
! **************************

  else if (method=="hlle2") then

!    Reconstruct primitive variables at grid interfaces.

     call reconstruct(Nx,dx)

!    Find fluxes.

     do i=0,Nx-1

!       Calculate left and right Lorentz factor.

        Wl = 1.0/sqrt(1.0 - fluid_vl(i)**2)
        Wr = 1.0/sqrt(1.0 - fluid_vr(i)**2)

!       Calculate left and right enthalpy.

        hl = 1.0 + fluid_el(i) + fluid_pl(i)/fluid_rho0l(i)
        hr = 1.0 + fluid_er(i) + fluid_pr(i)/fluid_rho0r(i)

!       Calculate left and right speed of sound.

        vsl = sqrt(gamma*fluid_pl(i)/(hl*fluid_rho0l(i)))
        vsr = sqrt(gamma*fluid_pr(i)/(hr*fluid_rho0r(i)))

!       Calculate left and right characteristic speeds.

        vml = (fluid_vl(i) - vsl)/(1.0 - fluid_vl(i)*vsl)
        vmr = (fluid_vr(i) - vsr)/(1.0 - fluid_vr(i)*vsr)

        vpl = (fluid_vl(i) + vsl)/(1.0 + fluid_vl(i)*vsl)
        vpr = (fluid_vr(i) + vsr)/(1.0 + fluid_vr(i)*vsr)

        vmm = min(0.0,vml,vmr)
        vpp = max(0.0,vpl,vpr)

!       Calculate left and right states of conserved quantities.

        Dl = Wl*fluid_rho0l(i)
        Dr = Wr*fluid_rho0r(i)

        El = Wl**2*hl*fluid_rho0l(i) - fluid_pl(i) - Dl
        Er = Wr**2*hr*fluid_rho0r(i) - fluid_pr(i) - Dr

        Sl = Wl**2*hl*fluid_rho0l(i)*fluid_vl(i)
        Sr = Wr**2*hr*fluid_rho0r(i)*fluid_vr(i)

!       Source for D.

        flux_D(i) = (vpp*fluid_vl(i)*Dl &
                   - vmm*fluid_vr(i)*Dr &
                   + vpp*vmm*(Dr - Dl))/(vpp - vmm)

!       Source for E.

        flux_E(i) = (vpp*fluid_vl(i)*(El + fluid_pl(i)) &
                   - vmm*fluid_vr(i)*(Er + fluid_pr(i)) &
                   + vpp*vmm*(Er - El))/(vpp - vmm)

!       Source for S.

        flux_S(i) = (vpp*(fluid_vl(i)*Sl + fluid_pl(i)) &
                   - vmm*(fluid_vr(i)*Sr + fluid_pr(i)) &
                   + vpp*vmm*(Sr - Sl))/(vpp - vmm)

  end do


! **************************
! ***   UNKNOWN METHOD   ***
! **************************

  else

     print *
     print *, 'Unknown method, aborting ...'
     print *
     stop

  end if


! ****************************
! ***   OUTER BOUNDARIES   ***
! ****************************

! At the moment just copied.

  flux_D(-1) = flux_D(0)
  flux_E(-1) = flux_E(0)
  flux_S(-1) = flux_S(0)

  flux_D(Nx) = flux_D(Nx-1)
  flux_E(Nx) = flux_E(Nx-1)
  flux_S(Nx) = flux_S(Nx-1)


! *****************************************
! ***   CALCULATE SOURCES FROM FLUXES   ***
! *****************************************

  do i=0,Nx
     sfluid_cD(i) = - (flux_D(i) - flux_D(i-1))*idx
     sfluid_cE(i) = - (flux_E(i) - flux_E(i-1))*idx
     sfluid_cS(i) = - (flux_S(i) - flux_S(i-1))*idx
  end do


! ***************
! ***   END   ***
! ***************

  end subroutine sources








  subroutine primitive(Nx,dx,gamma,q1,q2)

! ***************************************
! ***   RECOVER PRIMITIVE VARIABLES   ***
! ***************************************

! This routine recovers the fluid primitive variables (rho0,p,v)
! starting from the conserved quantities (D,E,S).
!
! Since the system of equations we need to invert is coupled and
! non-linear, we do the inversion using Newton-Raphson's method.
!
! The idea is to choose a trial value of the pressure

! Declare variables.

  use arrays

  implicit none

  integer i,l,Nx
  integer maxiter

  real(8) :: p1,p2,f1,f2
  real(8) :: dx,gamma,q1,q2
  real(8) :: W,res,small
  real(8) :: aux


! *********************************
! ***   LOOP OVER GRID POINTS   ***
! *********************************

  maxiter = 100

! Loop over grid points.

  do i=0,Nx

!    Initialize counter and residue.

     l = 0
     res = 1.0

!    If the density is too low we run the risk of dividing by a
!    very small quantity to recover the fluid speed which will
!    introduce large round-off errors.  In order to avoid this
!    here we just set the density to a small value.

     small = 1.0D-10

     if (fluid_cD(i)<=small) then
        fluid_cD(i) = small
     end if

!    Start iterations for Newton-Rapson's method.

     p1 = fluid_p(i)

     do while ((abs(res)>1.0D-12).and.(l<maxiter))

!       Increment counter.

        l = l + 1

!       Find value of v using old value of pressure.

        fluid_v(i) = fluid_cS(i)/(fluid_cE(i) + fluid_cD(i) + p1 + fluid_q(i))

!       Find Lorentz factor.

        aux = 1.0 - fluid_v(i)**2

        if (aux>0.0) then
           W = 1.0/sqrt(aux)
        else
           print *
           print *, 'Negative root in Lorentz factor at point: i = ',i,', x = ',x(i)
           print *, fluid_v(i),fluid_cS(i),fluid_cE(i),fluid_cD(i),p1,fluid_q(i)
           stop
        end if

!       Find rest-mass density rho0.

        fluid_rho0(i) = fluid_cD(i)/W

!       Find specific internal energy.

        fluid_e(i) = (fluid_cE(i) + fluid_cD(i)*(1.0 - W) &
                   + (p1 + fluid_q(i))*(1.0 - W**2))/(fluid_cD(i)*W)

        if (fluid_e(i)<0.0) fluid_e(i) = 0.0

!       Find enthalpy.

        fluid_h(i) = 1.0 + fluid_e(i) + (p1 + fluid_q(i))/fluid_rho0(i)

!       Find speed of sound.

        aux = gamma*(p1+ fluid_q(i))/(fluid_rho0(i)*fluid_h(i))

        if (aux>=0.0) then
           fluid_vs(i) = sqrt(aux)
        else
           print *
           print *, 'Negative root in speed of sound at point: i = ',i,', x = ',x(i)
           print *, p1,fluid_rho0(i),fluid_h(i)
           print *
           stop
        end if

!       Update artificial viscosity.

        aux = 0.5*(fluid_cS(i+1) - fluid_cS(i-1))

        if (aux>=0.0) then
            fluid_q(i) = 0.0
        else
            fluid_q(i) = abs(aux)*(q1*fluid_vs(i) + q2*abs(aux)/(fluid_cE(i) + fluid_cD(i)))
        end if

!       Calculate diference between trial value of pressure
!       and predicted value.

        f1 = (gamma - 1.0)*fluid_rho0(i)*fluid_e(i) - p1

!       Update the value of the pressure for the next iteration.

        if (l==1) then
           aux = (gamma - 1.0)*fluid_rho0(i)*fluid_e(i)
        else
           aux = p1 - f1*(p2-p1)/(f2-f1)
        end if

        p2 = p1
        f2 = f1

        p1 = aux

        if (p1<0.0) p1=0.0

!       Calculate residual.

        res = abs(p2 - p1)

!       print *, i,res,(fluid_cD(i)-fluid_rho0(i)*W)

     end do

!    Save last value of p1 as the new pressure.

     fluid_p(i) = p1

!    If we exceeded the maximum number of iterations send message to screen.

     if (l==maxiter) then
        print *
        print *, 'Maximum iteration number reached in subroutine primitive at point: i = ',i,', x = ',x(i)
        print *, 'Residue = ',res
        print *
        stop
     end if

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine primitive







  subroutine reconstruct(Nx,dx)

! This subroutine makes a linear reconstruction of
! the primitive variables at the cell interfaces
! (i.e. in between grid points) using a limiter
! with both left and right sided extrapolations.
!
! Having such reconstructions available is important
! for several of the methods used in the sources.

! Declare variables.

  use arrays

  implicit none

  integer i,Nx

  real(8) :: dx
  real(8) :: slope1,slope2,slope3,slopelim


! *********************************
! ***   LOOP OVER GRID POINTS   ***
! *********************************

  do i=0,Nx

!    Reconstruct rho0.

     slope1 = fluid_rho0(i  ) - fluid_rho0(i-1)
     slope2 = fluid_rho0(i+1) - fluid_rho0(i  )
     slope3 = fluid_rho0(i+2) - fluid_rho0(i+1)

     if (slope1*slope2>0.0) then
        if (abs(slope1)<abs(slope2)) then
           slopelim = slope1
        else
           slopelim = slope2
        end if
!       slopelim = 2.0*slope1*slope2/(slope1 + slope2)
     else
        slopelim = 0.0
     end if

     fluid_rho0l(i) = fluid_rho0(i  ) + 0.5*slopelim

     if (slope2*slope3>0.0) then
        if (abs(slope2)<abs(slope3)) then
           slopelim = slope2
        else
           slopelim = slope3
        end if
!       slopelim = 2.0*slope2*slope3/(slope2 + slope3)
     else
        slopelim = 0.0
     end if

     fluid_rho0r(i) = fluid_rho0(i+1) - 0.5*slopelim

!    Reconstruct e.

     slope1 = fluid_e(i  ) - fluid_e(i-1)
     slope2 = fluid_e(i+1) - fluid_e(i  )
     slope3 = fluid_e(i+2) - fluid_e(i+1)

     if (slope1*slope2>0.0) then
        if (abs(slope1)<abs(slope2)) then
           slopelim = slope1
        else
           slopelim = slope2
        end if
!       slopelim = 2.0*slope1*slope2/(slope1 + slope2)
     else
        slopelim = 0.0
     end if

     fluid_el(i) = fluid_e(i  ) + 0.5*slopelim

     if (slope2*slope3>0.0) then
        if (abs(slope2)<abs(slope3)) then
           slopelim = slope2
        else
           slopelim = slope3
        end if
!       slopelim = 2.0*slope2*slope3/(slope2 + slope3)
     else
        slopelim = 0.0
     end if

     fluid_er(i) = fluid_e(i+1) - 0.5*slopelim

!    Reconstruct v.

     slope1 = fluid_v(i  ) - fluid_v(i-1)
     slope2 = fluid_v(i+1) - fluid_v(i  )
     slope3 = fluid_v(i+2) - fluid_v(i+1)

     if (slope1*slope2>0.0) then
        if (abs(slope1)<abs(slope2)) then
           slopelim = slope1
        else
           slopelim = slope2
        end if
        slopelim = 2.0*slope1*slope2/(slope1 + slope2)
     else
        slopelim = 0.0
     end if

     fluid_vl(i) = fluid_v(i  ) + 0.5*slopelim

     if (slope2*slope3>0.0) then
        if (abs(slope2)<abs(slope3)) then
           slopelim = slope2
        else
           slopelim = slope3
        end if
        slopelim = 2.0*slope2*slope3/(slope2 + slope3)
     else
        slopelim = 0.0
     end if

     fluid_vr(i) = fluid_v(i+1) - 0.5*slopelim

!    Reconstruct p.

     slope1 = (fluid_p(i  ) + fluid_q(i  )) - (fluid_p(i-1) + fluid_q(i-1))
     slope2 = (fluid_p(i+1) + fluid_q(i+1)) - (fluid_p(i  ) + fluid_q(i  ))
     slope3 = (fluid_p(i+2) + fluid_q(i+2)) - (fluid_p(i+1) + fluid_q(i+1))

     if (slope1*slope2>0.0) then
        if (abs(slope1)<abs(slope2)) then
           slopelim = slope1
        else
           slopelim = slope2
        end if
 !      slopelim = 2.0*slope1*slope2/(slope1 + slope2)
     else
        slopelim = 0.0
     end if

     fluid_pl(i) = fluid_p(i  ) + fluid_q(i  ) + 0.5*slopelim

     if (slope2*slope3>0.0) then
        if (abs(slope2)<abs(slope3)) then
           slopelim = slope2
        else
           slopelim = slope3
        end if
!       slopelim = 2.0*slope2*slope3/(slope2 + slope3)
     else
        slopelim = 0.0
     end if

     fluid_pr(i) = fluid_p(i+1) + fluid_q(i+1) - 0.5*slopelim

  end do


! ***************
! ***   END   ***
! ***************

  end subroutine reconstruct
