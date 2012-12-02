program main

!*****************************************************************************80
!
!! MAIN is the main program for MD_OPENMP.
!
!  Discussion:
!
!    MD_OPENMP implements a simple molecular dynamics simulation.
!
!    The program uses OpenMP directives to allow parallel computation.
!
!    The velocity Verlet time integration scheme is used. 
!
!    The particles interact with a central pair potential.
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
!    This FORTRAN90 version by John Burkardt.
!
  use omp_lib

  implicit none

  integer ( kind = 4 ), parameter :: nd = 3
  integer ( kind = 4 ), parameter :: np = 1000

  real ( kind = 8 ) acc(nd,np)
  real ( kind = 8 ) box(nd)
  real ( kind = 8 ), parameter :: dt = 0.0001D+00
  real ( kind = 8 ) e0
  real ( kind = 8 ) force(nd,np)
  integer ( kind = 4 ) id
  real ( kind = 8 ) kinetic
  real ( kind = 8 ), parameter :: mass = 1.0D+00
  real ( kind = 8 ) pos(nd,np)
  real ( kind = 8 ) potential
  integer ( kind = 4 ) proc_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_num = 400
  integer ( kind = 4 ) step_print
  integer ( kind = 4 ) step_print_index
  integer ( kind = 4 ) step_print_num
  integer ( kind = 4 ) thread_num
  real ( kind = 8 ) vel(nd,np)
  real ( kind = 8 ) wtime

  call timestamp ( )

  proc_num = omp_get_num_procs ( )
  thread_num = omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MD_OPENMP'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A molecular dynamics program.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NP, the number of particles in the simulation is ', np
  write ( *, '(a,i8)' ) '  STEP_NUM, the number of time steps, is ', step_num
  write ( *, '(a,g14.6)' ) '  DT, the size of each time step, is ', dt
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of processors available is: ', proc_num
  write ( *, '(a,i8)' ) '  The number of threads available is:    ', thread_num
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

  wtime = omp_get_wtime ( )

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

  wtime = omp_get_wtime ( ) - wtime
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Elapsed time for main computation:'
  write ( *, '(2x,g14.6,a)' ) wtime, ' seconds'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MD_OPENMP'
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
!    This FORTRAN90 version by John Burkardt.
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
  kin = 0.0D+00

!$omp parallel &
!$omp shared ( f, nd, np, pos, vel ) &
!$omp private ( d, d2, i, j, rij )

!$omp do reduction ( + : pot, kin )

  do i = 1, np
!
!  Compute the potential energy and forces.
!
    f(1:nd,i) = 0.0D+00

    do j = 1, np

      if ( i /= j ) then

        call dist ( nd, pos(1,i), pos(1,j), rij, d )
!
!  Attribute half of the potential energy to particle J.
!
        d2 = min ( d, PI2 )

        pot = pot + 0.5D+00 * ( sin ( d2 ) )**2

        f(1:nd,i) = f(1:nd,i) - rij(1:nd) * sin ( 2.0D+00 * d2 ) / d

      end if

    end do
!
!  Compute the kinetic energy.
!
    kin = kin + sum ( vel(1:nd,i)**2 )

  end do
!$omp end do

!$omp end parallel

  kin = kin * 0.5D+00 * mass
  
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
!    This FORTRAN90 version by John Burkardt.
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
!    This FORTRAN90 version by John Burkardt.
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
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
!$omp parallel &
!$omp shared ( box, nd, np, pos ) &
!$omp private ( j )

!$omp do

  do j = 1, np
    pos(1:nd,j) = box(1:nd) * pos(1:nd,j)
  end do

!$omp end do
!$omp end parallel
!
!  Velocities and accelerations begin at 0.
!
!$omp parallel &
!$omp shared ( acc, nd, np, vel )

!$omp workshare

  vel(1:nd,1:np) = 0.0D+00
  acc(1:nd,1:np) = 0.0D+00

!$omp end workshare

!$omp end parallel

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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
  character ( len = 5 )  zone

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
!    This FORTRAN90 version by John Burkardt.
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

!$omp parallel &
!$omp shared ( acc, dt, f, nd, np, pos, rmass, vel ) &
!$omp private ( i, j )

!$omp do
  do j = 1, np
    do i = 1, nd
      pos(i,j) = pos(i,j) + vel(i,j) * dt + 0.5D+00 * acc(i,j) * dt * dt
      vel(i,j) = vel(i,j) + 0.5D+00 * dt * ( f(i,j) * rmass + acc(i,j) )
      acc(i,j) = f(i,j) * rmass
    end do
  end do
!$omp end do

!$omp end parallel

  return
end
