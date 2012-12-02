module particles

!*****************************************************************************80
!
!! PARTICLES is a module for atom property data.
!
!  Modified:
!
!    01 November 2005
!
!  Parameters:
!
!    integer :: DIM, the spatial dimension.  Reasonable values are 2 or 3.
!
!    real ( kind = 8 ) :: ACC(N), the particle accelerations.
!
!    real ( kind = 8 ) :: BOX_SIZE(DIM), the length scale factor.
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
!    logical :: VEL_ACC, is TRUE if the file includes velocities
!    and accelerations.
!
!    real ( kind = 8 ) :: VIRIAL, the virial term.
!
!    real ( kind = 8 ) :: VOLUME, the volume of the box.
!
  integer, parameter :: dim = 3

  real ( kind = 8 ), dimension(:,:), allocatable :: acc
  real ( kind = 8 ), dimension (dim) :: box_size
  real ( kind = 8 ), dimension(:), allocatable :: coord
  real ( kind = 8 ), dimension(:,:), allocatable :: deltar
  real ( kind = 8 ) :: density
  real ( kind = 8 ), dimension(:), allocatable :: deru
  real ( kind = 8 ), dimension(:), allocatable :: ene_kin
  real ( kind = 8 ), dimension(:), allocatable :: ene_pot
  integer :: n = 0
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
!  Parameters:
!
!    real ( kind = 8 ) :: DELTAT, the time steps, in reduced units.
!
!    character ( len = 80 ) :: SAMPIN, the name of the input file.
!
!    character ( len = 80 ) :: SAMPOUT, the name of the output file.
!
!    real ( kind = 8 ) :: SKIN, extra range for the neighbor list.
!    (used by MD3 and MD3GLUE only).
!
!    integer :: STEP_NUM, the number of time steps to do.
!
!    logical :: TEMP_CONSTANT, is TRUE if TEMP_REQUESTED is nonnegative.
!
!    real ( kind = 8 ) :: TEMP_REQUESTED, the desired temperature.  Any value
!    less than zero leaves it constant.
!
!    character ( len = 80 ) :: TITLE, an arbitrary title string.
!
!    real ( kind = 8 ) :: RHO_REQUESTED, the desired density.  A value
!    of 0 leaves it unchanged.
!
!    logical :: RHO_CHANGE, is TRUE if the user is changing the density.
!
!    real ( kind = 8 ) :: k_B, the Boltzmann constant.
!    (used by MD3GLUE only).
!
  real ( kind = 8 ) :: deltat
  character ( len = 80 ) :: sampin
  character ( len = 80 ) :: sampout
  real ( kind = 8 ) :: skin
  integer :: step_num
  logical :: temp_constant
  real ( kind = 8 ) :: temp_requested
  character ( len = 80 ) :: title
  real ( kind = 8 ) :: rho_requested
  logical :: rho_change
  real ( kind = 8 ), parameter :: k_B = 8.617385D-05

end
module potential

!*****************************************************************************80
!
!! POTENTIAL is a module with the parameters of the Lennard-Jones potential.
!
!  Modified:
!
!    02 November 2005
!
!  Parameters:
!
!    real ( kind = 8 ), R_CUTOFF, the cutoff distance.
!
!    real ( kind = 8 ), phicutoff, the value of the potential at the
!    cutoff distance.
!
!    integer, TableSize, the size of the table.
!
!    real ( kind = 8 ), R_MIN
!
!    real ( kind = 8 ), rsqMin = R_MIN**2.
!
!    real ( kind = 8 ), DeltaRsq, the spacing in the table data, as a function
!    of radius squared. 
!
!    real ( kind = 8 ), InvDeltaRsq = 1 / DeltaRsq.
!
!    real ( kind = 8 ), PhiTab, the tabulated value of the potential.
!
!    real ( kind = 8 ), DPhiTab, the value 1/r dphi/dr.
!
  real ( kind = 8 ), parameter :: r_cutoff = 2.5D+00
  real ( kind = 8 ), parameter :: r_min = 0.5D+00
  integer, parameter :: TableSize = 20001
  real ( kind = 8 ), dimension(TableSize) :: PhiTab
  real ( kind = 8 ), dimension(TableSize) :: DPhiTab

  real ( kind = 8 ), parameter :: phicutoff = &
                               4.0D+00/(r_cutoff**12) - 4.0D+00/(r_cutoff**6)
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
!  Parameters:
!
!    real ( kind = 8 ) :: ENE_KIN_SUM, the running sum of the kinetic energy.
!
!    real ( kind = 8 ) :: ENE_POT_SUM, the running sum of the potential energy.
!
!    real ( kind = 8 ) :: PRESSURE_SUM, the running sum of the pressure.
!
!    real ( kind = 8 ) :: TEMPERATURE_SUM, the running sum of the temperature.
!
  real ( kind = 8 ) :: ene_kin_sum = 0.0D+00
  real ( kind = 8 ) :: ene_pot_sum = 0.0D+00
  real ( kind = 8 ) :: pressure_sum  = 0.0D+00
  real ( kind = 8 ) :: temperature_sum = 0.0D+00

end
program main

!*****************************************************************************80
!
!! MAIN is the main program for MD1.
!
!  Discussion:
!
!    MD1 is a simple, minimal molecular dynamics program in Fortran90
!
!    MD1 uses the Lennard-Jones potential and the 'velocity' Verlet 
!    time integration algorithm.
!
!    It computes kinetic and potential energy, density and pressure.
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
!    The standard output of this program is suitable to be directly fed to 
!    'gnuplot' to produce plots as a function of time.  This is because
!    everything but data lines is commented out with an initial '#' character.
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
!    MP Allen, DJ Tildsley,
!    Computer Simulation of Liquids,
!    Oxford University Press, 1987.
!
!    Furio Ercolessi,
!    A Molecular Dynamics Primer.
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) '#  MD1'
  write ( *, '(a)' ) '#    FORTRAN90 version'
  write ( *, '(a)' ) '#    A Molecular Dynamics simulation.'
  write ( *, '(a)' ) '#    Last modified on 09 January 2006.'

  call initialize ( )

  call evolve_sample ( )

  call terminate ( )

  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) '#  MD1:'
  write ( *, '(a)' ) '#    Normal end of execution.'

  write ( *, '(a)' ) '#'
  call timestamp ( )

  stop
end
subroutine compute_forces ( )

!*****************************************************************************80
!
!! COMPUTE_FORCES computes the forces on atoms.
!
!  Discussion:
!
!    The forces are computed based on the atomic positions, using the
!    Lennard-Jones potential.
!    Note double nested loop, giving O(N^2) time: this is a SLOW ROUTINE,
!    unsuitable for large systems.
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
  use potential

  implicit none

  real ( kind = 8 ) :: dphi
  integer :: i
  integer :: j
  real ( kind = 8 ) :: phi
  real ( kind = 8 ), dimension(dim) :: rij
  real ( kind = 8 ) :: rm12
  real ( kind = 8 ) :: rm2
  real ( kind = 8 ) :: rm6
  real ( kind = 8 ) :: rsqij
  real ( kind = 8 ), dimension(dim) :: sij
!
!  Reset to zero the potential energies, forces, virial term.
!
  ene_pot = 0.0D+00
  acc = 0.0D+00
  virial = 0.0D+00
!
!  Loop over all pairs of particles.
!
  do i = 1, n - 1
    do j = i + 1, n

      sij(1:dim) = pos(1:dim,i) - pos(1:dim,j)

      where ( 0.5D+00 < abs ( sij ) )
        sij = sij - sign ( 1.0D+00, sij )
      end where 
 
      rij(1:dim) = box_size(1:dim) * sij(1:dim)

      rsqij = dot_product ( rij(1:dim), rij(1:dim) )
!
!  If the particles are close enough, compute the interaction.
!
      if ( rsqij < r_cutoff**2 ) then

        rm2 = 1.0D+00 / rsqij
        rm6 = rm2**3
        rm12 = rm6**2

        phi = 4.0D+00 * ( rm12 - rm6 ) - phicutoff

        dphi = 24.0D+00 * rm2 * ( 2.0D+00 * rm12 - rm6 )
!
!  Accumulate energy, virial, and forces.
!
        ene_pot(i) = ene_pot(i) + 0.5D+00 * phi
        ene_pot(j) = ene_pot(j) + 0.5D+00 * phi
        virial = virial - dphi * rsqij
        acc(1:dim,i) = acc(1:dim,i) + dphi * sij(1:dim)
        acc(1:dim,j) = acc(1:dim,j) - dphi * sij(1:dim)

      end if

    end do
  end do
!
!  Definition of the virial term.
!
  virial = - virial / real ( dim, kind = 8 )                       

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
!  Parameters:
!
!    Output, real ( kind = 8 ) ENE_KIN_AVER, the average kinetic
!    energy per particle.
!
!    Output, real ( kind = 8 ) TEMPERATURE, the instantaneous
!    temperature of the particle system.
!
  use particles

  implicit none

  real ( kind = 8 ) :: ene_kin_aver
  integer :: i
  real ( kind = 8 ), dimension(dim) :: real_vel
  real ( kind = 8 ) :: temperature

  do i = 1, n
    real_vel(1:dim) = box_size(1:dim) * vel(1:dim,i)
    ene_kin(i) = 0.5D+00 * dot_product ( real_vel(1:dim), real_vel(1:dim) )
  end do

  ene_kin_aver = sum ( ene_kin(1:n) ) / real ( n, kind = 8 )

  temperature = 2.0D+00 * ene_kin_aver / real ( dim, kind = 8 )

  return
end
subroutine evolve_sample ( )

!*****************************************************************************80
!
!! EVOLVE_SAMPLE controls the time evolution of the system.
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
  use statistics

  implicit none

  real ( kind = 8 ) :: chi
  real ( kind = 8 ) :: ene_kin_aver
  real ( kind = 8 ) :: ene_pot_aver
  real ( kind = 8 ) :: ene_tot_aver
  real ( kind = 8 ) :: pressure
  integer step
  real ( kind = 8 ) :: temperature
!
!  We need to have the initial temperature ready in case we are going
!  at constant T:
!
  call compute_temperature ( ene_kin_aver, temperature ) 
!
!  "Velocity Verlet" integrator (see e.g. Allen and Tildesley book, p. 81).
!  Simple velocity scaling (done on velocities at the previous step)
!  applied when temp_constant is enabled.
!
  do step = 1, step_num

    call refold_positions ( )

    pos = pos + deltat * vel + 0.5D+00 * deltat * deltat * acc
!
!  For constant T, rescale the velocity.
!
    if ( temp_constant .and. 0.0D+00 < temperature ) then
      call compute_temperature ( ene_kin_aver, temperature )
      chi = sqrt ( temp_requested / temperature )
      vel = chi * vel + 0.5D+00 * deltat * acc
    else
      vel = vel + 0.5D+00 * deltat * acc
    end if  
 
    call compute_forces ( )

    vel = vel + 0.5D+00 * deltat * acc

    call compute_temperature ( ene_kin_aver, temperature )

    ene_pot_aver = sum( ene_pot ) / real ( n, kind = 8 )
    ene_tot_aver = ene_kin_aver + ene_pot_aver
!
!  For the pressure calculation, see the Allen and Tildesley book, section 2.4
!
    pressure = density * temperature + virial / volume

    write ( *, '(i6,f14.6,f14.6,f14.6,f14.6,f14.6)' ) step, temperature, &
      ene_kin_aver, ene_pot_aver, ene_tot_aver, pressure
!
!  Accumulate statistics.
!
    temperature_sum = temperature_sum + temperature
    ene_kin_sum = ene_kin_sum + ene_kin_aver
    ene_pot_sum = ene_pot_sum + ene_pot_aver
    pressure_sum = pressure_sum + pressure

  end do

  return
end
subroutine initial_printout ( )

!*****************************************************************************80
!
!! INITIAL_PRINTOUT prints information on the run parameters.
!
!  Discussion:
!
!    Leading '#' characters are to directly use the output as a gnuplot input.
!
!  Modified:
!
!    01 November 2005
!
!  Author:
!
!    Furio Ercolessi
!
  use particles
  use simulation_control

  implicit none

  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) '# MD1: a minimal molecular dynamics program'
  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) '#' // trim ( title )
  write ( *, '(a)' ) '#  Input sample: "' // trim ( sampin ) // '".'

  if ( vel_acc ) then
    write ( *, '(a)' ) &
      '#          (positions, velocities, accelerations read from file)'
  else
    write ( *, '(a)' ) &
       '#                                (only positions read from file)'
  end if

  write ( *, '(a)' ) '# Output sample: "' // trim ( sampout ) // '".'

  write ( *, '(a,i8,a,f7.4,a,f12.4)' ) &
   '# Number of steps:',step_num,', time step:', deltat, &
   ', total time:', step_num * deltat
  write ( *, '(a,i6)' ) '# Number of atoms:', n
  write ( *, '(a,3f12.6,a,f15.3)' ) &
   '# Box size:', box_size(1:dim), ', Volume:', volume

  if ( rho_change ) then
    write ( *, '(a,f12.6,a)' ) '# Density:', density, ' (changed)'
  else
    write ( *, '(a,f12.6,a)' ) '# Density:', density, ' (unchanged)'
  end if

  if ( temp_constant ) then
    write ( *, '(a,f12.6)' ) '# Constant T run with T =', temp_requested
  else
    write ( *, '(a)' ) '# Free evolution run.'
  end if
!
!  Print headers of columns.
!
  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) &
    '#  Step   Temperature     Kinetic      Potential   ' // &
    '  Total Energy    Pressure'
  write ( *, '(a)' ) &
    '# -----  ------------  ------------  ------------  ' // &
    '------------  ------------'
  write ( *, '(a)' ) '#'

  return
end
subroutine initialize ( )

!*****************************************************************************80
!
!! INITIALIZE controls the initialization procedure.
!
!  Discussion:
!
!    This routine is called once at the beginning, before the
!    time evolution.
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
!
!  Read the user directives controlling the simulation.
!
  call read_input ( )
!
!  Read the input sample containing the initial particle coordinates
!  and perhaps also velocities and accelerations.
!
  call read_sample ( )
!
!  Print information.
!
  call initial_printout ( )

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
    temperature_sum               / real ( step_num, kind = 8 ), &
    ene_kin_sum                   / real ( step_num, kind = 8 ), &
    ene_pot_sum                   / real ( step_num, kind = 8 ), &
    ( ene_kin_sum + ene_pot_sum ) / real ( step_num, kind = 8 ), &
    pressure_sum                  / real ( step_num, kind = 8 )

  return
end
subroutine read_input ( )

!*****************************************************************************80
!
!! READ_INPUT reads the parameters controlling the simulation.
!
!  Discussion:
!
!    This routine reads an input file from the user containing
!    the values of certain program parameters.
!
!  Modified:
!
!    28 October 2005
!
!  Author:
!
!    Furio Ercolessi
!
  use simulation_control

  implicit none

  integer io_status

  read ( *, '(a)', iostat = io_status ) title

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_INPUT - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading TITLE.'
    stop
  end if

  read ( *, '(a)', iostat = io_status ) sampin

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_INPUT - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading SAMPIN.'
    stop
  end if

  read ( *, '(a)', iostat = io_status ) sampout

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_INPUT - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading SAMPOUT.'
    stop
  end if

  read ( *, *, iostat = io_status ) step_num

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_INPUT - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading STEP_NUM.'
    stop
  end if

  read ( *, *, iostat = io_status ) deltat

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_INPUT - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading DELTAT.'
    stop
  end if

  read ( *, *, iostat = io_status ) rho_requested

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_INPUT - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading RHO_REQUESTED.'
    stop
  end if

  read ( *, *, iostat = io_status ) temp_requested

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_INPUT - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading TEMP_REQUESTED.'
    stop
  end if

  rho_change = ( 0.0D+00 < rho_requested )
  temp_constant = ( 0.0D+00 <= temp_requested )

  return
end
subroutine read_sample ( )

!*****************************************************************************80
!
!! READ_SAMPLE reads the initial sample from file unit 1.
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

  real ( kind = 8 ), dimension(dim) :: AccAtomReal
  integer :: i
  integer :: io_status
  integer, parameter :: io_unit = 1
  integer :: k
  real ( kind = 8 ), dimension(dim) :: mass_center
  real ( kind = 8 ), dimension(dim) :: PosAtomReal
  real ( kind = 8 ) :: scale
  real ( kind = 8 ), dimension(dim) :: VelAtomReal

  open ( unit = io_unit, file = sampin, status = 'old', &
    action = 'read', iostat = io_status )   

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    stop
  end if

  read ( io_unit, '(1x,l2,i7,3e23.15)', iostat = io_status ) &
    vel_acc, n, box_size(1:dim)

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Error reading line 1 of the input file.'
    stop
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_SAMPLE - Fatal error!'
    write ( *, '(a,i6)' ) '  N is ', n
    stop
  end if
!
!  Compute the volume and density.  They do not change in the run.
!
  volume = product ( box_size(1:dim) )
  density = real ( n, kind = 8 ) / volume
!
!  Unless the user wants to change the density, in this case we do
!  it here:
!
  if ( rho_change ) then
    scale = ( density / rho_requested ) ** ( 1.0D+00 / dim )
    box_size(1:dim) = scale * box_size(1:dim)
    volume  = product ( box_size(1:dim) )
    density = real ( n, kind = 8 ) / volume
  end if
!
!  Now that we know the system size, we can dynamically allocate the
!  arrays containing atomic information.
!
  allocate ( acc(dim,n) )
  allocate ( ene_kin(n) )
  allocate ( ene_pot(n) )
  allocate ( pos(dim,n) )
  allocate ( vel(dim,n) )
!
!  Read the coordinates from the file (one line per atom), normalize
!  them to the box size along each direction and store them.
!  Energies are set initially to zero.
!
  do i = 1, n

    read ( io_unit, *, iostat = io_status ) PosAtomReal(1:dim)

    if ( io_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'READ_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  Error reading PosAtomReal array.'
      stop
    end if

    pos(1:dim,i) = PosAtomReal(1:dim) / box_size(1:dim)

  end do

  ene_pot(1:n) = 0.0D+00
  ene_kin(1:n) = 0.0D+00
!
!  For "new" samples (that is, samples just created by defining the atomic
!  coordinates and not the result of previous simulations), we have now
!  read everything, and velocities and accelerations are set to zero.
!
!  For samples which have been produced by previous simulations, we also
!  have to read velocities and accelerations.
!
!  The logical variable VEL_ACC distinguishes between these two cases.
!
  if ( vel_acc ) then

    do i = 1, n

      read ( io_unit, '(1x,3e23.15)', iostat = io_status ) VelAtomReal

      if ( io_status /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'READ_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Error reading VelAtomReal array.'
        stop
      end if

      vel(1:dim,i) = VelAtomReal / box_size(1:dim)

    end do

    do i = 1, n

      read ( io_unit, '(1x,3e23.15)', iostat = io_status ) AccAtomReal

      if ( io_status /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'READ_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Error reading AccAtomReal array.'
        stop
      end if

      acc(1:dim,i) = AccAtomReal / box_size(1:dim)

    end do

  else

    vel(1:dim,1:n) = 0.0D+00

    acc(1:dim,1:n) = 0.0D+00

  end if
!
!  Compute the center of mass.
!
  mass_center(1:dim) = sum ( pos(1:dim,1:n), dim = 2 ) / real ( n, kind = 8 )
!
!  Translate the atoms so that center of mass is at the origin.
!
  if ( .false. ) then

    write ( *, '(a)' ) '#'
    write ( *, '(a)' ) '#READ_SAMPLE: Warning!'
    write ( *, '(a)' ) '#  The coordinate data is translated to have'
    write ( *, '(a)' ) '#  zero center of mass.'

    do k = 1, dim
      pos(k,1:n) = pos(k,1:n) - mass_center(k)
    end do

  else

    write ( *, '(a)' ) '#'
    write ( *, '(a)' ) '#READ_SAMPLE: Warning!'
    write ( *, '(a)' ) '#  The coordinate data is NOT translated to have'
    write ( *, '(a)' ) '#  zero center of mass.'

  end if

  close ( unit = io_unit )

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
subroutine terminate ( )

!*****************************************************************************80
!
!! TERMINATE carries out the termination procedures.
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
!
!  Print the mean statistics.
!
  call print_statistics ( )
!
!  Write the final sample to a file.
!
  call write_sample ( )
!
!  Deallocate all the dynamic arrays.
!
  deallocate ( acc )
  deallocate ( ene_kin )
  deallocate ( ene_pot )
  deallocate ( pos )
  deallocate ( vel )

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

  write ( *, '(a,2x,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    '#', d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine write_sample ( )

!*****************************************************************************80
!
!! WRITE_SAMPLE writes the final sample.
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

  real ( kind = 8 ), dimension(dim) :: AccAtomReal
  integer :: i
  integer :: io_status
  integer, parameter :: io_unit = 2
  real ( kind = 8 ), dimension(dim) :: PosAtomReal
  real ( kind = 8 ), dimension(dim) :: VelAtomReal

  open ( unit = io_unit, file = sampout, status = 'replace', &
    action = 'write', iostat = io_status )

  if ( io_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if
!
!  VEL_ACC = TRUE because we will write velocities and accelerations.
!
  vel_acc = .true.
!
!  The '%' is a comment character, and allows information to be inserted
!  into the file that can be used by some postprocessing programs,
!  and ignored by others.
!
  write ( io_unit, '(a1,l2,i7,3e23.15)' ) '%', vel_acc, n, box_size(1:dim)
!
!  Multiply coordinates (which are scaled by the box size) by box_size in
!  order to have them in real, unscaled units, then write them in the
!  file (one line per atom).
!
  do i = 1, n
    PosAtomReal(1:dim) = pos(1:dim,i) * box_size(1:dim)
    write ( io_unit, '(1x,3e23.15)' ) PosAtomReal(1:dim)
  end do
!
!  Do the same for velocities and accelerations.
!
  do i = 1, n
    VelAtomReal(1:dim) = vel(1:dim,i) * box_size(1:dim)
    write ( io_unit, '(a1,3e23.15)' ) '%', VelAtomReal(1:dim)
  end do

  do i = 1, n
    AccAtomReal(1:dim) = acc(1:dim,i) * box_size(1:dim)
    write ( io_unit, '(a1,3e23.15)' ) '%', AccAtomReal(1:dim)
  end do

  close ( unit = io_unit )

  return
end
