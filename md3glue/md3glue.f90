module particles

!*****************************************************************************80
!
!! PARTICLES is a module for atom property data.
!
!  Modified:
!
!    01 November 2005
!
!  Author:
!
!    Furio Ercolessi
!
!  Parameters:
!
!    integer :: DIM, the spatial dimension.  Reasonable values are 2 or 3.
!
!    real ( kind = 8 ) :: ACC(N), the particle accelerations.
!
!    real ( kind = 8 ) :: BOX_SIZE(DIM), the length of the sides of the box.
!
!    real ( kind = 8 ) :: COORD(N), the particle coordinations.
!
!    real ( kind = 8 ) :: DELTAR(N), the particle displacements.
!
!    real ( kind = 8 ) :: DENSITY, the density of the box.
!
!    real ( kind = 8 ) :: DERU(N), U'(n_i).
!
!    real ( kind = 8 ) :: ENE_KIN(N), the particle kinetic energies.
!
!    real ( kind = 8 ) :: ENE_POT(N), the particle potential energies.
!
!    integer :: N, the  number of particles.
!
!    real ( kind = 8 ) :: POS(N), the particle positions.
!
!    real ( kind = 8 ) :: VEL(N), the particle velocities.
!
!    logical :: VEL_ACC, is TRUE if the input file includes velocities
!    and accelerations.
!
!    real ( kind = 8 ) :: VIRIAL, the virial term.
!
!    real ( kind = 8 ) :: VOLUME, the volume of the box.
!
  integer, parameter :: DIM = 3

  real ( kind = 8 ), dimension(:,:), allocatable :: acc
  real ( kind = 8 ), dimension(DIM) :: box_size
  real ( kind = 8 ), dimension(:), allocatable :: coord
  real ( kind = 8 ), dimension(:,:), allocatable :: deltar
  real ( kind = 8 ) :: density
  real ( kind = 8 ), dimension(:), allocatable :: deru
  real ( kind = 8 ), dimension(:), allocatable :: ene_kin
  real ( kind = 8 ), dimension(:), allocatable :: ene_pot
  integer :: N = 0
  real ( kind = 8 ), dimension(:,:), allocatable :: pos
  real ( kind = 8 ), dimension(:,:), allocatable :: vel
  logical :: vel_acc = .false.
  real ( kind = 8 ) :: virial
  real ( kind = 8 ) :: volume

end
module simulation_control

!*****************************************************************************80
!
!! SIMULATION_CONTROL is a module for simulation control data.
!
!  Discussion:
!
!    This data is supplied from user input. 
!
!  Modified:
!
!    31 October 2005
!
!  Author:
!
!    Furio Ercolessi
!
!  Parameters:
!
!    character ( len = 80 ) :: TITLE, an arbitrary title string.
!
!    character ( len = 80 ) :: SAMPIN, the name of the input file.
!
!    character ( len = 80 ) :: SAMPOUT, the name of the output file.
!
!    integer :: STEP_NUM, the number of time steps to do.
!
!    real ( kind = 8 ) :: DELTAT, the time steps, in reduced units.
!
!    real ( kind = 8 ) :: RHO_REQUESTED, the desired density.  A value
!    of 0 leaves it unchanged.
!
!    real ( kind = 8 ) :: TEMP_REQUESTED, the desired temperature.  Any value
!    less than zero leaves it constant.
!
!    logical :: RHO_CHANGE, is TRUE if the user is changing the density.
!
!    logical :: TEMP_CONSTANT, is TRUE if TEMP_REQUESTED is nonnegative.
!
!    real ( kind = 8 ) :: SKIN, extra range for the neighbor list.
!    (used by MD3 and MD3GLUE only).
!
!    real ( kind = 8 ) :: k_B, the Boltzmann constant.
!    (used by MD3GLUE only).
!
  character ( len = 80 ) :: title
  character ( len = 80 ) :: sampin
  character ( len = 80 ) :: sampout
  integer :: step_num
  real ( kind = 8 ) :: deltat
  real ( kind = 8 ) :: rho_requested
  real ( kind = 8 ) :: temp_requested
  logical :: rho_change
  logical :: temp_constant
  real ( kind = 8 ) :: skin
  real ( kind = 8 ), parameter :: k_B = 8.617385D-05

end
module potential_glue

!*****************************************************************************80
!
!! POTENTIAL_GLUE is a module containaing glue potential data.
!
!  Discussion:
!
!    This module contains parameters and tables of the ab-initio derived 
!    glue potential for Aluminum [F. Ercolessi and James B. Adams, 
!    Europhysics Letters 26, 583 (1994) ].
!
!    The tables are filled by calling externally defined subroutines
!    V2, RH and UU which return the actual potential.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  integer, parameter :: TableSize = 20001

  real ( kind = 8 ), parameter :: CoordMax = 1.2d0
  real ( kind = 8 ), dimension(TableSize) :: DPhiTab ! 1/r dphi/dr
  real ( kind = 8 ), dimension(TableSize) :: DRhoTab ! 1/r drho/dr
  real ( kind = 8 ), dimension(TableSize) :: DUTab      !du/dn
  real ( kind = 8 ), dimension(TableSize) :: PhiTab !phi(r)
  real ( kind = 8 ), parameter :: r_cutoff = 0.555805441821810D+01 ! cutoff
  real ( kind = 8 ), dimension(TableSize) :: RhoTab !rho(r)
  real ( kind = 8 ), parameter :: r_min = 0.055805441821810D+01 ! cutoff-5 A
  real ( kind = 8 ), dimension(TableSize) :: UTab !U(n)

  real ( kind = 8 ), parameter :: DeltaCoord = CoordMax / ( TableSize - 1 )

  real ( kind = 8 ), parameter :: InvDeltaCoord = 1.0D+00 / DeltaCoord

  real ( kind = 8 ), parameter :: RsqMin = r_min**2

  real ( kind = 8 ), parameter :: DeltaRsq = &
                                  ( r_cutoff**2 - RsqMin ) / ( TableSize - 1 )

  real ( kind = 8 ), parameter :: InvDeltaRsq = 1.0D+00 / DeltaRsq

end
module statistics

!*****************************************************************************80
!
!! STATISTICS is a module storing statistical quantities.
!
!  Discussion:
!
!    This module contains statistical quantities accumulated during the run.
!    All quantities are set to zero at the beginning of the run.
!    At the end of each time step, the new value is added to the running sum.
!    At the end of all the steps, the running sums are used to compute
!    a mean value.
!
!  Modified:
!
!    31 October 2005
!
!  Author:
!
!    Furio Ercolessi
!
!  Parameters:
!
!    real ( kind = 8 ) :: TEMPERATURE_SUM, the running sum of the temperature.
!
!    real ( kind = 8 ) :: ENE_KIN_SUM, the running sum of the kinetic energy.
!
!    real ( kind = 8 ) :: ENE_POT_SUM, the running sum of the potential energy.
!
!    real ( kind = 8 ) :: PRESSURE_SUM, the running sum of the pressure.
!
  real ( kind = 8 ) :: temperature_sum = 0.0D+00
  real ( kind = 8 ) :: ene_kin_sum = 0.0D+00
  real ( kind = 8 ) :: ene_pot_sum = 0.0D+00
  real ( kind = 8 ) :: pressure_sum = 0.0D+00

end
module neighbor_list

!*****************************************************************************80
!
!! NEIGHBOR_LIST is a module containing data about the neighbor list.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  integer, parameter :: MaxPairsPerAtom = 100
  integer, dimension (:), allocatable :: Advance
  integer, dimension (:), allocatable :: List
  integer, dimension (:), allocatable :: Marker1
  integer, dimension (:), allocatable :: Marker2
  integer :: MaxListLength, ListLength
  real ( kind = 8 ), dimension (:,:), allocatable :: displace_list

end
program main

!*****************************************************************************80
!
!! MAIN is the main program for MD3GLUE.
!
!  Discussion:
!
!    MD3GLUE is a simple, minimal molecular dynamics program.
!
!    MD3GLUE is the same as MD3, but uses an ab-initio derived glue 
!    many-body potential for Aluminum, instead of a Lennard-Jones potential.
!
!    Now units are:
!      Distances in A, 
!      energies in eV, 
!      temperatures in K.
!
!    Equilibrium structure at T=0: fcc, a=4.032 A, E=-3.36 eV/atom 
!
!    This program must be linked with the supplied file 'aluminum.f'.
!
!    Files used by this program:
!
!    Unit  I/O  Meaning
!    ----  ---  ----------------------------------------------------------------
!      1    I   Input sample (coordinates, and perhaps also velocities and
!               accelerations) read at the beginning of the run
!      2    O   Output sample (coordinates, velocities, accelerations) written
!               at the end of the run
!      *    I   Standard input for the simulation control
!      *    O   Standard output containing various informations on the run
!
!    Output suitable to be directly fed to 'gnuplot' to produce plots as a
!    function of time.
!
!  Modified:
!
!    09 January 2006
!
!  Author:
!
!    Furio Ercolessi
!
!  Reference:
!
!    M P Allen and D J Tildsley,
!    Computer Simulation of Liquids,
!    Oxford University Press, 1987.
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) '#  MD3GLUE'
  write ( *, '(a)' ) '#    FORTRAN90 version'
  write ( *, '(a)' ) '#    A Molecular Dynamics simulation.'
  write ( *, '(a)' ) '#    Last modified on 09 January 2006.'

  call initialize ( )

  call evolve_sample ( )

  call terminate ( )

  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) '#  MD3_GLUE:'
  write ( *, '(a)' ) '#    Normal end of execution.'

  write ( *, '(a)' ) '#'
  call timestamp ( )

  stop
end
subroutine initialize ( )

!*****************************************************************************80
!
!! INITIALIZE is the initialization procedure,
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  use Particles

  implicit none
!
!  Read the user directives controlling the simulation.
!
  call read_input ( )
!
!  Read the input sample containing the initial particle coordinates
!  (and perhaps also velocities and accelerations)
!
  call read_sample ( )
!
!  Print informations on the run on the standard output file.
!
  call initial_printout ( )
!
!  Define the potential tables.
!
  call define_potential_tables ( )

  return
end
subroutine read_sample ( )

!*****************************************************************************80
!
!! READ_SAMPLE reads the initial sample from file unit 1.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  use Particles
  use Neighbor_List
  use Simulation_Control

  implicit none

  real ( kind = 8 ), dimension(DIM) :: AccAtomReal
  real ( kind = 8 ), dimension(DIM) :: PosAtomReal
  real ( kind = 8 ), dimension(DIM) :: VelAtomReal
  real ( kind = 8 ), dimension(DIM) :: mass_center
  real ( kind = 8 ) :: scale
  integer :: i,k,lis

  lis = len_trim(SampIn)
  open(unit=1,file=SampIn(1:lis),status='old', &
     action='read',err=700)   
  read(1,'(1X,L2,I7,3E23.15)',end=600,err=900) vel_acc, N, box_size(1:DIM)

  if ( N <= 0 ) then
    print*,'Read_Sample: FATAL: N is',N
    stop
  end if
!
!  Compute volume and density once for all (they do not change in the run)
!
  volume  = product ( box_size(1:DIM) )
  density = N / volume
!
!  Unless the user wants to change the density, in this case we do
!  it here:
!
  if ( rho_change ) then
    scale = ( density / rho_requested ) ** ( 1.d0 / DIM )
    box_size(1:DIM) = scale * box_size(1:DIM)
    volume  = product ( box_size(1:DIM) )
    density = N / volume
  end if
!
!  Now that we know the system size, we can dynamically allocate the
!  arrays containing atomic informations
!
  allocate ( deltar(DIM,N) )
  allocate ( pos(DIM,N) )
  allocate ( vel(DIM,N) )
  allocate ( acc(DIM,N) )
  allocate ( ene_pot(N) )
  allocate ( ene_kin(N) )
  allocate ( coord(N) )
  allocate ( deru(N) )
!
!  Also those related to the neighbor list:
!
  MaxListLength = MaxPairsPerAtom*N
  allocate ( List(MaxListLength) )
  allocate ( Marker1(N) )
  allocate ( Marker2(N) )
  allocate ( Advance(N) )
  allocate ( displace_list(DIM,N) )
!
!  Read the coordinates from the file (one line per atom), normalize
!  them to the box size along each direction and store them.
!  Energies are set initially to zero.
!
  do i=1,N
    read(1,*,end=800,err=900) PosAtomReal(1:DIM)
    pos(1:DIM,i) = PosAtomReal(1:DIM) / box_size(1:DIM)
  end do

  ene_pot(1:N) = 0.0D+00
  ene_kin(1:N) = 0.0D+00
!
!  For "new" samples (that is, samples just created by defining the atomic
!  coordinates and not the result of previous simulations), we have now
!  read everything, and velocities and accelerations are set to zero.
!  For sample which have been produced by previous simulations, we also
!  have to read velocities and accelerations.
!  The logical variable VEL_ACC distinguishes between these two cases.
!
  if ( vel_acc ) then

    do i=1,N
      read(1,'(1X,3E23.15)',end=800,err=900) VelAtomReal(1:DIM)
      vel(1:DIM,i) = VelAtomReal(1:DIM) / box_size(1:DIM)
    enddo

    do i=1,N
      read(1,'(1X,3E23.15)',end=800,err=900) AccAtomReal(1:DIM)
      acc(1:DIM,i) = AccAtomReal(1:DIM) / box_size(1:DIM)
    enddo

  else

    vel(1:DIM,1:N) = 0.0D+00
    acc(1:DIM,1:N) = 0.0D+00

  endif
!
!  Compute center of mass coordinates.
!
  mass_Center = sum( pos , dim=2 ) / N
!
!  Translate the atoms so that center of mass is at the origin.
!
  if ( .false. ) then

    write ( *, * ) '#'
    write ( *, * ) '#READ_SAMPLE: Warning!'
    write ( *, * ) '#  The coordinate data is translated to have'
    write ( *, * ) '#  zero center of mass.'

    do k = 1, dim
      pos(k,1:n) = pos(k,1:n) - mass_center(k)
    end do

  else

    write ( *, * ) '#'
    write ( *, * ) '#READ_SAMPLE: Warning!'
    write ( *, * ) '#  The coordinate data is NOT translated to have'
    write ( *, * ) '#  zero center of mass.'

  end if
!
!  All coordinates read successfully if we get to this point.
!
  close ( unit = 1 )

  return
!
!  handling of various kinds of errors
!
600 continue
   print*,'Read_Sample: FATAL: ',SampIn(1:lis),' is empty?'
   stop
700 continue
   print*,'Read_Sample: FATAL: ',SampIn(1:lis),' not found.'
   stop
800 continue
   print*,'Read_Sample: FATAL: premature end-of-file at atom ',i
   close(unit=1)
   stop   
900 continue
   print*,'Read_Sample: FATAL: read error in ',SampIn(1:lis)
   close(unit=1)
   stop
end
subroutine read_input ( )

!*****************************************************************************80
!
!! READ_INPUT reads simulation parameters.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  use Simulation_Control

  implicit none
!
!  Read the input parameters controlling the simulation from the standard
!  input.
!
  read(*,'(a)',end=200,err=800) title
  read(*,'(a)',end=200,err=800) SampIn
  read(*,'(a)',end=200,err=800) SampOut
  read(*,  *  ,end=200,err=800) step_num
  read(*,  *  ,end=200,err=800) deltat
  read(*,  *  ,end=200,err=800) rho_requested
  read(*,  *  ,end=200,err=800) temp_requested
  read(*,  *  ,end=210,err=800) Skin
  temp_requested = k_B * temp_requested              ! user sees K, we use eV internally

  rho_change = ( rho_requested > 0 )
  temp_constant = ( temp_requested >= 0 )
  return   !  successful exit
210 continue
   print*,'The last input line (containing the Skin) is missing. This is md3!'
200 continue
   print*,'Read_Input: FATAL: premature end-of-file in standard input'
   stop
800 continue
   print*,'Read_Input: FATAL: read error in standard input'
   stop
end
subroutine initial_printout ( )

!*****************************************************************************80
!
!! INITIAL_PRINTOUT prints simulation parameter information.
!
!  Discussion:
!
!    Leading # are to directly use the output as a gnuplot input.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  use Particles
  use Neighbor_List
  use Simulation_Control

  implicit none

  print '(a,/,2a,/,a)','#','# MD3glue: a minimal molecular dynamics program,', &
                 ' using a glue model for Al','#'
  print '(2a)','# ',title(1:len_trim(title))
  print '(2a)','#  Input sample: ',SampIn(1:len_trim(SampIn))
  if ( vel_acc ) then
    print '(a)','#          (positions, velocities, accelerations read from file)'
  else
    print '(a)','#                                (only positions read from file)'
  end if

  print '(2a)','# Output sample: ',SampOut(1:len_trim(SampOut))
  print '(a,i8,a,f7.4,a,f12.4)', &
   '# Number of steps:', step_num,', time step:',deltat, &
   ', total time:', step_num * deltat
  print '(a,i6)','# Number of atoms:',N
  print '(a,3f12.6,a,f15.3)', &
   '# Box size:',box_size(1:DIM),', Volume:',volume

  if ( rho_change ) then
    print '(a,f12.6,a)','# Density:',density,' (changed)'
  else
    print '(a,f12.6,a)','# Density:',density,' (unchanged)'
  end if

  if (temp_constant) then
    print '(a,f12.2)','# Constant T run with T =',temp_requested/k_B
  else
    print '(a)','# Free evolution run.'
  end if

  print '(a,f8.4,a,i9)', &
      '# Skin:',Skin,' ,  maximum neighbor list length:',MaxListLength
!
!  Now print headers of columns
!
  print '(a,/,a,/,a)', '#', &
'#  Step   Temperature     Kinetic      Potential   Total Energy    Pressure',&
'# -----  ------------  ------------  ------------  ------------  ------------'

  return
end
subroutine define_potential_tables ( )

!*****************************************************************************80
!
!! DEFINE_POTENTIAL_TABLES computes tables for the glue potential for Al.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  use potential_glue

  implicit none

  real ( kind = 8 ) :: c
  real ( kind = 8 ) :: d2glue
  real ( kind = 8 ) :: d2phi
  real ( kind = 8 ) :: d2rho
  real ( kind = 8 ) :: dglue
  real ( kind = 8 ) :: dphi
  real ( kind = 8 ) :: drho
  real ( kind = 8 ) :: glue
  integer :: k
  real ( kind = 8 ) :: phi
  real ( kind = 8 ) :: r
  real ( kind = 8 ) :: rho
  real ( kind = 8 ) :: Rsq

  do k = 1, TableSize
    Rsq = RsqMin + dble ( k - 1 ) * DeltaRsq
    r = sqrt ( Rsq )
    call v2 ( r, phi, dphi, d2phi )
    call rh ( r, rho, drho, d2rho )
    PhiTab(k)  = phi
    DPhiTab(k) = - dphi / r
    RhoTab(k)  = rho
    DRhoTab(k) = - drho / r
  end do

  do k = 1, TableSize
    c = dble ( k - 1 ) * DeltaCoord
    call uu ( c, glue, dglue, d2glue )
    UTab(k) = glue
    DUTab(k) = dglue
  end do

  return
end
subroutine evolve_sample ( )

!*****************************************************************************80
!
!! EVOLVE_SAMPLE controls the time evolution of the system.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  use Particles
  use Neighbor_List
  use Simulation_Control
  use potential_glue
  use Statistics

  implicit none

  integer step,i
  real ( kind = 8 ) :: ene_kin_aver
  real ( kind = 8 ) :: ene_pot_aver
  real ( kind = 8 ) :: ene_tot_aver
  real ( kind = 8 ) :: temperature,pressure
  real ( kind = 8 ) :: chi
  logical, external :: MovedTooMuch
  logical ::  ListUpdateRequested = .TRUE. ! always need an update when starting
!
!  We need to have the initial temperature ready in case we are going
!  at constant T:
!
  call compute_temperature ( ene_kin_aver, temperature ) 
!
!  "Velocity Verlet" integrator (see e.g. Allen and Tildesley book, p. 81).
!  Simple velocity scaling (done on velocities at the previous step)
!  applied when TEMP_CONSTANT is enabled.
!
time: do step=1, step_num

   call refold_positions ( )

   deltar = deltat*vel + 0.5d0*(deltat**2)*acc         ! dr = r(t+dt) - r(t)
!
!  deltar has now defined to allow updating of displace_list
!
   pos = pos + deltar                                  ! r(t+dt)
   if (temp_constant .and. (temperature > 0) ) then	  !  veloc rescale for const T
      call Compute_Temperature(ene_kin_aver,temperature) ! T(t)
      chi = sqrt( temp_requested / temperature )
      vel = chi*vel + 0.5d0*deltat*acc                 ! v(t+dt/2)
   else                          !  regular constant E dynamics
      vel = vel + 0.5d0*deltat*acc                     ! v(t+dt/2)
   endif   
   if (ListUpdateRequested) then                       ! list update required !
      call Update_List(r_cutoff+Skin)                   ! so we update the list
      ListUpdateRequested = .FALSE.                    ! request fulfilled
   endif
   call Compute_Forces ( )                                 ! a(t+dt),ene_pot,virial
   vel = vel + 0.5d0*deltat*acc                        ! v(t+dt)
   call Compute_Temperature(ene_kin_aver,temperature)  ! at t+dt, also ene_kin
   ene_pot_aver = sum( ene_pot ) / N
   ene_tot_aver = ene_kin_aver + ene_pot_aver
!
!  For the pressure calculation, see the Allen and Tildesley book, section 2.4
!
   pressure = density*temperature + virial/volume
!
!  Update DISPLACE_LIST, which contains the displacements of atoms
!  since the last list update.  Note that this is in real space units,
!  therefore the scaled displacement deltar is multiplied by box_size.
!
   do i=1,N
      displace_list(:,i) = displace_list(:,i) + box_size(1:DIM) * deltar(:,i)
   enddo
!
!  Now we perform the list deterioration test. If particles 'moved too
!  much' relative to Skin, a list update should be scheduled to occur
!  at the next step.  See MovedTooMuch for the definition of 'moved too much'.
!
   ListUpdateRequested = MovedTooMuch(Skin)
!
!  Print step report with instantaneous values
!
   print '(1x,i6,f14.2,4f14.6)',step,temperature/k_B, &
                          ene_kin_aver,ene_pot_aver,ene_tot_aver,pressure
!
!  Accumulate statistics:
!
   temperature_sum = temperature_sum + temperature
   ene_kin_sum     = ene_kin_sum     + ene_kin_aver
   ene_pot_sum     = ene_pot_sum     + ene_pot_aver
   pressure_sum    = pressure_sum    + pressure
end do time

  return
end
subroutine refold_positions ( )

!*****************************************************************************80
!
!! REFOLD_POSITIONS folds exiting particles back into the box.
!
!  Discussion:
!
!    Periodic boundary conditions are used.
!
!  Modified:
!
!    28 October 2005
!
!  Author:
!
!    Furio Ercolessi
!
  use particles

  implicit none

  where ( 0.5D+00 < pos ) 
    pos = pos - 1.0D+00
  end where

  where ( pos < -0.50D+00 ) 
    pos = pos + 1.0D+00
  end where

  return
end
subroutine compute_forces ( )

!*****************************************************************************80
!
!! COMPUTE_FORCES computes forces on atoms.
!
!  Discussion:
!
!    Compute forces on atoms from the positions, using the potential stored
!    in the potential tables.  Note how this routine is valid for any
!    two-body potential.
!
!    This version makes use of an existing neighbor list, and runs in O(N)
!    time.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  use Particles
  use Neighbor_List
  use potential_glue

  implicit none

  real ( kind = 8 ), dimension(DIM) :: Sij,Rij
  real ( kind = 8 ) :: Rsqij,rho,drho,phi,dphi,rk,weight
  integer :: i,j,k,L
!
!  Reset to zero coordinations, potential energies
!
  coord = 0.d0
  ene_pot = 0.d0
! 
!  Loop over particles in the neighbor list a first time, to construct the
!  coordinations.
!
do i = 1,N
   do L = Marker1(i),Marker2(i)                   ! part of List useful for i
      j = List(L)                                 ! take neigh index out of List
      Sij = pos(:,i) - pos(:,j)                   ! distance vector between i j
      where ( abs(Sij) > 0.5d0 )                  ! (in box scaled units)
         Sij = Sij - sign(1.d0,Sij)               ! periodic boundary conditions
      end where                                   ! applied where needed.
      Rij(1:DIM) = box_size(1:DIM) * Sij(1:DIM)                           ! go to real space units
      Rsqij = dot_product(Rij,Rij)                ! compute square distance
      if ( Rsqij < r_cutoff**2 ) then              ! particles are interacting?
         rk = (Rsqij - RsqMin)*InvDeltaRsq + 1.d0 ! "continuous" index in table
         k = int(rk)                              ! discretized index
         if (k < 1) k = 1                         ! unlikely but just to protect
         weight = rk - k                          ! fractional part, in [0,1]
         rho  = weight* RhoTab(k+1) + (1.d0-weight)* RhoTab(k)  ! do linear int.
         coord(i) = coord(i) + rho                ! accumulate atomic density
         coord(j) = coord(j) + rho                ! (i and j share it)
      endif
   enddo
enddo
!
!  Now that we know n_i (coord), calculate glue part of energies U(N_i),
!  and derivatives U'(n_i) (deru):
!
  do i = 1,N
    rk = coord(i)*InvDeltaCoord + 1.d0             ! "continuous" index in table
    k = int(rk)                                    ! discretized index
    if ( k >= TableSize ) k = TableSize - 1        ! unlikely but just to protect
    weight = rk - k                                ! fractional part, in [0,1]
    ene_pot(i) = weight*UTab(k+1) + (1.d0-weight)*UTab(k) ! linear interpolation
    deru(i) = weight*DUTab(k+1) + (1.d0-weight)*DUTab(k)  ! for U(n_i)/U'(n_i)
  end do
!
!  Reset to zero forces, virial term
!
  acc = 0.d0
  virial = 0.d0
!
!  Loop over particles in the neighbor list a second time, to compute
!  the forces:
!
do i = 1,N
   do L = Marker1(i),Marker2(i)                   ! part of List useful for i
      j = List(L)                                 ! take neigh index out of List
      Sij = pos(:,i) - pos(:,j)                   ! distance vector between i j
      where ( abs(Sij) > 0.5d0 )                  ! (in box scaled units)
         Sij = Sij - sign(1.d0,Sij)               ! periodic boundary conditions
      end where                                   ! applied where needed.
      Rij(1:DIM) = box_size(1:DIM) * Sij(1:DIM)                           ! go to real space units
      Rsqij = dot_product(Rij,Rij)                ! compute square distance
      if ( Rsqij < r_cutoff**2 ) then              ! particles are interacting?
         rk = (Rsqij - RsqMin)*InvDeltaRsq + 1.d0 ! "continuous" index in table
         k = int(rk)                              ! discretized index
         if (k < 1) k = 1                         ! unlikely but just to protect
         weight = rk - k                          ! fractional part, in [0,1]
         phi  = weight* PhiTab(k+1) + (1.d0-weight)* PhiTab(k)  ! do linear
         ene_pot(i) = ene_pot(i) + 0.5d0*phi      ! accumulate energy
         ene_pot(j) = ene_pot(j) + 0.5d0*phi      ! (i and j share it)
         dphi = weight*DPhiTab(k+1) + (1.d0-weight)*DPhiTab(k)  ! interpolation
         drho = weight*DRhoTab(k+1) + (1.d0-weight)*DRhoTab(k)  ! interpolation
         dphi = dphi + ( deru(i) + deru(j) )*drho ! effective potential
         virial = virial - dphi*Rsqij             ! accumul. virial=sum r(dV/dr)
         acc(:,i) = acc(:,i) + dphi*Sij           ! accumulate forces
         acc(:,j) = acc(:,j) - dphi*Sij           !    (Fji = -Fij)
      endif
   enddo
enddo
virial = - virial/DIM                             ! definition of virial term\

  return
end
subroutine compute_temperature ( ene_kin_aver, temperature )

!*****************************************************************************80
!
!! COMPUTE_TEMPERATURE updates the kinetic energy and temperature.
!
!  Discussion:
!
!    Starting from the velocities currently stored in VEL, update
!    the kinetic energy array ENE_KIN, and compute ENE_KIN_AVER,
!    the average kinetic energy per particle, and the
!    instantaneous temperature.
!
!  Modified:
!
!    28 October 2005
!
!  Author:
!
!    Furio Ercolessi
!
  use particles

  implicit none

  real ( kind = 8 ) :: ene_kin_aver
  integer :: i
  real ( kind = 8 ), dimension(DIM) :: real_vel
  real ( kind = 8 ) :: temperature

  do i = 1, N
    real_vel(1:DIM) = box_size(1:DIM) * vel(1:DIM,i)
    ene_kin(i) = 0.5D+00 * dot_product ( real_vel(1:DIM), real_vel(1:DIM) )
  end do

  ene_kin_aver = sum ( ene_kin(1:N) ) / N

  temperature = 2.0D+00 * ene_kin_aver / DIM

  return
end
subroutine update_list ( range )

!*****************************************************************************80
!
!! UPDATE_LIST updates the neighbor list.
!
!  Discussion:
!
!    Update the neighbor list, including all atom pairs that are within a
!    distance Range to each other.  Range may be equal to the interaction
!    range but is usually a made somewhat bigger to have the list remain
!    valid for a while.
!
!    The Fincham-Ralston method is used, which allows vectorization of the
!    first of the two inner loop on vector computers and therefore gives
!    good performance.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
!  Reference:
!
!    David Fincham, BJ Ralston,
!    Molecular Dynamics Simulation Using the Cray-1 Vector Processing Computer,
!    Computer Physics Communications,
!    Volume 23, pages 127-134, 1981,
!
  use Particles
  use Neighbor_List

  implicit none

  real ( kind = 8 ) :: Range, RangeSq
  real ( kind = 8 ), dimension(DIM) :: Sij,Rij
  real ( kind = 8 ) :: Rsqij
  integer :: i,j,L

  write(6,'(a)',advance='NO') '# Neighbor list update ..'   ! beginning of report
                        ! line, # is to have gnuplot forgetting about this line
  RangeSq = Range**2

!
!  The following is the Fincham-Ralston loop for list updating
!
  L = 1                                             ! initialize index into List
  do i = 1,N                                        ! looping an all pairs
    do j = i+1,N
      Sij = pos(:,i) - pos(:,j)                   ! distance vector between i j
      where ( abs(Sij) > 0.5d0 )                  ! (in box scaled units)
         Sij = Sij - sign(1.d0,Sij)               ! periodic boundary conditions
      end where                                   ! applied where needed.
      Rij(1:DIM) = box_size(1:DIM) * Sij(1:DIM)                           ! go to real space units
      Rsqij = dot_product(Rij,Rij)                ! compute square distance
      if ( Rsqij < RangeSq ) then                 ! is j to be a neighbor of i?
         Advance(j) = 1                           !    if yes, put a 1,
      else                                        !    otherwise
         Advance(j) = 0                           !    put a 0, in Advance
      endif
    enddo

    Marker1(i) = L                                 ! list for i starts here

    do j = i+1,N
      if ( L > MaxListLength) goto 9999           ! panic exit if no more room
      List(L) = j                                 ! j will be really included ..
      L = L + Advance(j)                          ! only of Advance(j) is 1
    enddo

    Marker2(i) = L - 1                             ! list for i terminates here
  enddo

  ListLength = L - 1                                ! final used length of List
  write(6,'(i8,a)') ListLength,' indexes in list'   ! complete report line
  displace_list = 0.d0                               ! reset displacement vectors
  return                                            ! regular exit

9999 continue                                     ! here if no more room
   print '(/,a)','Update_List: FATAL: List too small for such a large Skin.'
   print '(a)','If you really want to use this Skin, you must increase the'
   print '(a)','value of parameter MaxPairsPerAtom in module Neighbor_List from'
   print '(i4,a)',MaxPairsPerAtom,' to something bigger, recompile and rerun.'
   stop
end
logical function MovedTooMuch ( Skin )

!*****************************************************************************80
!
!! MOVEDTOOMUCH determines if too much motion has taken place.
!
!  Discussion:
!
!    This routine scans the DISPLACE_LIST vector, containing the displacement
!    vectors of particles since the last time the list was updated.  It 
!    looks for the magnitude of the two largest displacements in the system.
!    These are summed together.  If the sum is larger than Skin, the function
!    returns TRUE, otherwise FALSE.
!
!  Modified:
!
!    25 February 2007
!
!  Author:
!
!    Furio Ercolessi
!
  use Particles
  use Neighbor_List

  implicit none

  real ( kind = 8 ) :: Skin
  real ( kind = 8 ) :: Displ,Displ1,Displ2
  integer :: i

  Displ1 = 0.0D+00
  Displ2 = 0.0D+00

  do i = 1, N

    Displ = sqrt( dot_product(displace_list(:,i),displace_list(:,i)) )

    if (Displ >= Displ1) then        !  1st position taken
      Displ2 = Displ1               !  push current 1st into 2nd place
      Displ1 = Displ                !  and put this one into current 1st
    else if (Displ >= Displ2) then   !  2nd position taken
      Displ2 = Displ
    end if

  end do

  MovedTooMuch = ( Displ1 + Displ2 > Skin )

  return
end
subroutine terminate ( )

!*****************************************************************************80
!
!! TERMINATE carries out the termination procedures.
!
!  Modified:
!
!    06 November 2005
!
!  Author:
!
!    Furio Ercolessi
!
  use particles
  use neighbor_list

  implicit none
!
!  Print a line with averages.
!
  call print_statistics ( )
!
!  Write the final sample.
!
  call write_sample ( )
!
!  Deallocate dynamic arrays.
!
  deallocate ( acc )
  deallocate ( advance )
  deallocate ( coord )
  deallocate ( deltar )
  deallocate ( deru )
  deallocate ( displace_list )
  deallocate ( ene_kin )
  deallocate ( ene_pot )
  deallocate ( list )
  deallocate ( marker1 )
  deallocate ( marker2 )
  deallocate ( pos )
  deallocate ( vel )

  return
end
subroutine print_statistics ( )

!*****************************************************************************80
!
!! PRINT_STATISTICS prints statistics from the calculation.
!
!  Discussion:
!
!    This routine prints the mean value, averaged during the run, of the
!    statistical quantities which have been accumulated.
!  
!  Modified:
!
!    04 November 2005
!
!  Author:
!
!    Furio Ercolessi
!
  use simulation_control
  use statistics

  implicit none

  if ( step_num <= 0 ) then
    return               
  end if

  write ( *, '(a)' ) '#'
  write ( *, '(a,f14.6,f14.6,f14.6,f14.6,f14.6)' ) &
    '# Means', &
    temperature_sum / step_num, &
    ene_kin_sum / step_num, &
    ene_pot_sum / step_num, &
    ( ene_kin_sum + ene_pot_sum ) / step_num, &
    pressure_sum / step_num         

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine write_sample

!*****************************************************************************80
!
!! WRITE_SAMPLE writes the final sample to file unit 2.
!
!  Discussion:
!
!    This sample file may be read by the BALLROOM program, for
!    visualization.  Since it only wants to read point coordinates,
!    the other data is "hidden" by "%" signs, which BALLROOM takes
!    as comment indicators.
!
!    This sample file may also be read as the input file for another
!    run of this program.  In that case, the routine which reads the
!    input data will skip over the comment character and "see"
!    the velocity and acceleration data.
!
!  Modified:
!
!    28 October 2005
!
!  Author:
!
!    Furio Ercolessi
!
  use particles
  use simulation_control

  implicit none

  real ( kind = 8 ), dimension(DIM) :: AccAtomReal
  integer :: i
  integer :: ios
  real ( kind = 8 ), dimension(DIM) :: PosAtomReal
  real ( kind = 8 ), dimension(DIM) :: VelAtomReal

  open ( unit = 2, file = sampout, status = 'replace', &
    action = 'write', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if
!
!  We are going to write velocities and accelerations to the file.
!
  vel_acc = .true.
!
!  The '%' signals to ignore this line to the BallRoom program producing
!  an image from the simulation
!
  write ( 2, '(a1,l2,i7,3e23.15)' ) '%', vel_acc, N, box_size(1:DIM)
!
!  Multiply coordinates (which are scaled by the box size) by box_size in
!  order to have them in real, unscaled units, then write them in the
!  file (one line per atom).
!
  do i = 1, N
    PosAtomReal(1:DIM) = pos(1:DIM,i) * box_size(1:DIM)
    write ( 2, '(1x,3e23.15)' ) PosAtomReal(1:DIM)
  end do
!
!  Do the same for velocities and accelerations.
!
  do i = 1, N
    VelAtomReal(1:DIM) = vel(1:DIM,i) * box_size(1:DIM)
    write ( 2, '(a1,3e23.15)' ) '%', VelAtomReal(1:DIM)
  end do

  do i = 1, N
    AccAtomReal(1:DIM) = acc(1:DIM,i) * box_size(1:DIM)
    write ( 2, '(a1,3e23.15)' ) '%', AccAtomReal(1:DIM)
  end do

  close ( unit = 2 )

  return
end
