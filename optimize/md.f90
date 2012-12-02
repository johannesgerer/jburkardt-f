program main

!*****************************************************************************80
!
!! MAIN is the main program for MD.
!
!  Discussion:
!
!    MD implements a simple molecular dynamics simulation.
!
!    The velocity Verlet time integration scheme is used. 
!
!    The particles interact with a central pair potential.
!
!  Usage:
!
!    md nd np step_num
!    where
!    * nd is the spatial dimension (2 or 3);
!    * np is the number of particles (500, for instance);
!    * step_num is the number of time steps (500, for instance).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2009
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This version by John Burkardt
!
  real ( kind = 8 ), allocatable :: acc(:,:)
  integer ( kind = 4 ) arg_num
  real ( kind = 8 ), allocatable :: box(:)
  real ( kind = 8 ) ctime
  real ( kind = 8 ) ctime1
  real ( kind = 8 ) ctime2
  real ( kind = 8 ), parameter :: dt = 0.0001D+00
  real ( kind = 8 ) e0
  real ( kind = 8 ), allocatable :: force(:,:)
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) kinetic
  integer ( kind = 4 ) last
  real ( kind = 8 ), parameter :: mass = 1.0D+00
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) np
  real ( kind = 8 ), allocatable :: pos(:,:)
  real ( kind = 8 ) potential
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) step
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) step_print
  integer ( kind = 4 ) step_print_index
  integer ( kind = 4 ) step_print_num
  character ( len = 255 ) string
  real ( kind = 8 ), allocatable :: vel(:,:)
  real ( kind = 8 ) wtime

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MD'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A molecular dynamics program.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get ND, the number of spatial dimensions.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, nd, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter ND, the spatial dimension (2 or 3 ):'
    read ( *, * ) nd
  end if
!
!  Get NP, the number of particles.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, np, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter NP, the number of particles (500, for instance):'
    read ( *, * ) np
  end if
!
!  Get STEP_NUM, the number of time steps.
!
  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, step_num, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter STEP_NUM, the number of time steps (500, for instance):'
    read ( *, * ) step_num
  end if
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ND, the spatial dimension, is ', nd
  write ( *, '(a,i8)' ) '  NP, the number of particles in the simulation is ', np
  write ( *, '(a,i8)' ) '  STEP_NUM, the number of time steps, is ', step_num
  write ( *, '(a,g14.6)' ) '  DT, the size of each time step, is ', dt
!
!  Allocate memory.
!
  allocate ( acc(nd,np) )
  allocate ( box(nd) )
  allocate ( force(nd,np) )
  allocate ( pos(nd,np) )
  allocate ( vel(nd,np) )
!
!  Set the dimensions of the box.
!
  box(1:nd) = 10.0D+00
!
!  Set initial positions, velocities, and accelerations.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initializing positions, velocities, and accelerations.'

  seed = 123456789
  call initialize ( np, nd, box, seed, pos, vel, acc )
!
!  Compute the forces and energies.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing initial forces and energies.'

  call compute ( np, nd, pos, vel, mass, force, potential, kinetic )
!
!  Save the initial total energy for use in the accuracy check.
!
  e0 = potential + kinetic

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  At each step, we report the potential and kinetic energies.'
  write ( *, '(a)' ) '  The sum of these energies should be a constant.'
  write ( *, '(a)' ) '  As an accuracy check, we also print the relative error'
  write ( *, '(a)' ) '  in the total energy.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Step      Potential       Kinetic        (P+K-E0)/E0'
  write ( *, '(a)' ) '                Energy P        Energy K       Relative Energy Error'
  write ( *, '(a)' ) ' '
!
!  This is the main time stepping loop:
!    Compute forces and energies,
!    Update positions, velocities, accelerations.
!
  step_print = 0
  step_print_index = 0
  step_print_num = 10
  
  step = 0
  write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    step, potential, kinetic, ( potential + kinetic - e0 ) / e0
  step_print_index = step_print_index + 1
  step_print = ( step_print_index * step_num ) / step_print_num

  call cpu_time ( ctime1 )

  do step = 1, step_num

    call compute ( np, nd, pos, vel, mass, force, potential, kinetic )

    if ( step == step_print ) then

      write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        step, potential, kinetic, ( potential + kinetic - e0 ) / e0

      step_print_index = step_print_index + 1
      step_print = ( step_print_index * step_num ) / step_print_num

    end if

    call update ( np, nd, pos, vel, force, acc, mass, dt )

  end do

  call cpu_time ( ctime2 )

  ctime = ctime2 - ctime1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Elapsed cpu time for main computation:'
  write ( *, '(2x,g14.6,a)' ) ctime, ' seconds'
!
!  Deallocate memory.
!
  deallocate ( acc )
  deallocate ( box )
  deallocate ( force )
  deallocate ( pos )
  deallocate ( vel )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MD:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine compute ( np, nd, pos, vel, mass, f, pot, kin )

!*****************************************************************************80
!
!! COMPUTE computes the forces and energies.
!
!  Discussion:
!
!    The computation of forces and energies is fully parallel.
!
!    The potential function V(X) is a harmonic well which smoothly
!    saturates to a maximum value at PI/2:
!
!      v(x) = ( sin ( min ( x, PI2 ) ) )**2
!
!    The derivative of the potential is:
!
!      dv(x) = 2.0D+00 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
!            = sin ( 2.0 * min ( x, PI2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2008
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NP, the number of particles.
!
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!
!    Input, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!
!    Input, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!
!    Input, real ( kind = 8 ) MASS, the mass of each particle.
!
!    Output, real ( kind = 8 ) F(ND,NP), the forces.
!
!    Output, real ( kind = 8 ) POT, the total potential energy.
!
!    Output, real ( kind = 8 ) KIN, the total kinetic energy.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) nd

  real ( kind = 8 ) d
  real ( kind = 8 ) d2
  real ( kind = 8 ) f(nd,np)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) kin
  real ( kind = 8 ) mass
  real ( kind = 8 ), parameter :: PI2 = 3.141592653589793D+00 / 2.0D+00
  real ( kind = 8 ) pos(nd,np)
  real ( kind = 8 ) pot
  real ( kind = 8 ) rij(nd)
  real ( kind = 8 ) vel(nd,np)

  pot = 0.0D+00

  do i = 1, np
!
!  Compute the potential energy and forces.
!
    f(1:nd,i) = 0.0D+00

    do j = 1, np

      if ( i /= j ) then

        rij(1:nd) = pos(1:nd,i) - pos(1:nd,j)

        d = sqrt ( sum ( rij(1:nd)**2 ) )
!
!  Truncate the distance.
!
        d2 = min ( d, PI2 )
!
!  Attribute half of the total potential energy to particle J.
!
        pot = pot + 0.5D+00 * sin ( d2 ) * sin ( d2 )
!
!  Add particle J's contribution to the force on particle I.
!
        f(1:nd,i) = f(1:nd,i) - rij(1:nd) * sin ( 2.0D+00 * d2 ) / d

      end if

    end do

  end do
!
!  Compute the total kinetic energy.
!
  kin = 0.5D+00 * mass * sum ( vel(1:nd,1:np)**2 )
  
  return
end
subroutine dist ( nd, r1, r2, dr, d )

!*****************************************************************************80
!
!! DIST computes the displacement and distance between two particles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2002
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!
!    Input, real ( kind = 8 ) R1(ND), R2(ND), the positions of the particles.
!
!    Output, real ( kind = 8 ) DR(ND), the displacement vector.
!
!    Output, real ( kind = 8 ) D, the Euclidean norm of the displacement,
!    in other words, the distance between the two particles.
!
  implicit none

  integer ( kind = 4 ) nd

  real ( kind = 8 ) d
  real ( kind = 8 ) dr(nd)
  real ( kind = 8 ) r1(nd)
  real ( kind = 8 ) r2(nd)

  dr(1:nd) = r1(1:nd) - r2(1:nd)

  d = sqrt ( sum ( dr(1:nd)**2 ) )

  return
end
subroutine initialize ( np, nd, box, seed, pos, vel, acc )

!*****************************************************************************80
!
!! INITIALIZE initializes the positions, velocities, and accelerations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 2007
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NP, the number of particles.
!
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!
!    Input, real ( kind = 8 ) BOX(ND), specifies the maximum position
!    of particles in each dimension.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!
!    Output, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!
!    Output, real ( kind = 8 ) ACC(ND,NP), the acceleration of each particle.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) nd

  real ( kind = 8 ) acc(nd,np)
  real ( kind = 8 ) box(nd)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) pos(nd,np)
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) vel(nd,np)
!
!  Start by setting the positions to random numbers between 0 and 1.
!
  call random_number ( harvest = pos(1:nd,1:np) )
!
!  Use these random values as scale factors to pick random locations
!  inside the box.
!
  do i = 1, nd
    pos(i,1:np) = box(i) * pos(i,1:np)
  end do
!
!  Velocities and accelerations begin at 0.
!
  vel(1:nd,1:np) = 0.0D+00
  acc(1:nd,1:np) = 0.0D+00

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters 
!    of S used to make the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * )  s
  integer ( kind = 4 ) state
  character :: TAB = achar ( 9 )
  integer ( kind = 4 ) value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
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
  integer ( kind = 4 ) d
  character ( len = 8 )  date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine update ( np, nd, pos, vel, f, acc, mass, dt )

!*****************************************************************************80
!
!! UPDATE updates positions, velocities and accelerations.
!
!  Discussion:
!
!    The time integration is fully parallel.
!
!    A velocity Verlet algorithm is used for the updating.
!
!    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
!    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
!    a(t+dt) = f(t) / m
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 2007
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NP, the number of particles.
!
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!
!    Input/output, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!
!    Input/output, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!
!    Input, real ( kind = 8 ) F(ND,NP), the force on each particle.
!
!    Input/output, real ( kind = 8 ) ACC(ND,NP), the acceleration of each
!    particle.
!
!    Input, real ( kind = 8 ) MASS, the mass of each particle.
!
!    Input, real ( kind = 8 ) DT, the time step.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) nd

  real ( kind = 8 ) acc(nd,np)
  real ( kind = 8 ) dt
  real ( kind = 8 ) f(nd,np)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mass
  real ( kind = 8 ) pos(nd,np)
  real ( kind = 8 ) rmass
  real ( kind = 8 ) vel(nd,np)

  rmass = 1.0D+00 / mass

  do j = 1, np
    do i = 1, nd
      pos(i,j) = pos(i,j) + vel(i,j) * dt + 0.5D+00 * acc(i,j) * dt * dt
      vel(i,j) = vel(i,j) + 0.5D+00 * dt * ( f(i,j) * rmass + acc(i,j) )
      acc(i,j) = f(i,j) * rmass
    end do
  end do

  return
end
