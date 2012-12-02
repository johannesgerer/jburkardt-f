program main

!*****************************************************************************80
!
!! MAIN is the main program for the reactor shielding simulation.
!
!  Discussion:
!
!    This is a Monte Carlo simulation, using
!    uniform random numbers, which investigates the
!    effectiveness of a shield intended to absorb the
!    neutrons emitted from a nuclear reactor.
!   
!    The reactor is modeled as a point source,
!    located at (0,0,0).
!   
!    A particle emitted from the reactor has a random
!    initial direction, and an energy selected from
!    [Emin,Emax] with a 1/Sqrt(E) distribution.
!   
!    The shield is modeled as a wall of thickness THICK,
!    extending from 0 to THICK in the X direction, and
!    extending forever in the Y and Z directions.
!   
!    Based on the particle energy, a distance D is computed
!    which measures how far the particle could travel through
!    the shield before colliding.
!   
!    Based on the particle direction, the position is updated
!    by D units.
!   
!    If the particle is now to the left of the shield, it is
!    counted as being REFLECTED.
!   
!    If the particle is to the right of the shield, it is 
!    counted as being ABSORBED.
!   
!    If the particle is inside the shield, it has COLLIDED.
!    A particle that collides is either absorbed (end of story)
!    or SCATTERED with a new random direction and a new (lower)
!    energy.
!   
!    Every particle is followed from origin to its final fate,
!    which is reflection, transmission, or absorption.
!    At the end, a summary is printed, giving the number of
!    particles with each fate, and the average energy of each
!    group of particles.
!   
!    Increasing NTOT, the number of particles used, will improve the
!    expected reliability of the results.
!   
!    Increasing THICK, the thickness of the shield, should 
!    result in more absorptions and reflections.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) AZM, the azimuthal angle of the particle's
!    direction.
!
!    Local, real ( kind = 8 ) D, the distance that the particle can
!    travel through the slab, given its current energy.
!
!    Local, real ( kind = 8 ) E, the energy of the particle.
!
!    Local, real ( kind = 8 ) EA, energy absorbed by the slab.
!
!    Local, real ( kind = 8 ) ER, energy reflected by the slab.
!
!    Local, real ( kind = 8 ) ET, energy transmitted through the slab.
!
!    Local, real ( kind = 8 ) MU, the cosine of the angle between the
!    particle's direction and the X axis.
!
!    Local, integer ( kind = 4 ) NA, number of particles absorbed by the slab.
!
!    Local, integer ( kind = 4 ) NPART, the index of the current particle.
!
!    Local, integer ( kind = 4 ) NR, number of particles reflected by the slab.
!
!    Local, integer ( kind = 4 ) NT, number of particles transmitted by the slab.
!
!    Local, integer ( kind = 4 ) NTOT, the total number of particles to be
!    emitted from the neutron source.
!
!    Local, real ( kind = 8 ) SA, standard deviation of absorbed energy.
!
!    Local, real ( kind = 8 ) SR, standard deviation of reflected energy.
!
!    Local, real ( kind = 8 ) ST, standard deviation of transmitted energy.
!
!    Local, real ( kind = 8 ) THICK, the thickness of the slab that is
!    intended to absorb most of the particles.
!
!    Local, real ( kind = 8 ) X, Y, Z, the current position of the particle.
!
  implicit none

  logical absorb
  real ( kind = 8 ) azm
  real ( kind = 8 ) d
  real ( kind = 8 ) dist2c
  real ( kind = 8 ) e
  real ( kind = 8 ) ea
  real ( kind = 8 ) er
  real ( kind = 8 ) et
  integer ( kind = 4 ) i
  real ( kind = 8 ) mu
  integer ( kind = 4 ) na
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt
  integer ( kind = 4 ), parameter :: ntot = 100000
  real ( kind = 8 ) sa
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sr
  real ( kind = 8 ) st
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 8 ), parameter :: thick = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REACTOR_SIMULATION'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  The reactor shielding simulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Shield thickness is THICK = ', thick
  write ( *, '(a,i8)' ) '  Number of simulated particles is NTOT = ', ntot
  write ( *, '(a,i8)' ) '  Number of tests TEST_NUM = ', test_num

  seed = 123456789

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i2)' ) '  Test # ', test
    write ( *, '(a,i12)' ) '  SEED = ', seed
!
!  Initialize.
!
    ea = 0.0D+00
    er = 0.0D+00
    et = 0.0D+00
    na = 0
    nr = 0
    nt = 0
    sa = 0.0D+00
    sr = 0.0D+00
    st = 0.0D+00
!
!  Loop over the particles.
!
    do npart = 1, ntot
!
!  Generate a new particle.
!
      call source ( seed, e, mu, azm, x, y, z )

      do
!
!  Compute the distance that the particle can travel through the slab,
!  based on its current energy.
!
        d = dist2c ( e, seed )
!
!  Update the particle's position by D units.
!
        call update ( mu, azm, d, x, y, z )
!
!  The particle was reflected by the shield, and this path is complete.
!
        if ( x < 0.0D+00 ) then

          nr = nr + 1
          er = er + e
          sr = sr + e * e
          exit
!
!  The particle was transmitted through the shield, and this path is complete.
!
        else if ( thick < x ) then

          nt = nt + 1
          et = et + e
          st = st + e * e
          exit
!
!  The particle collided with the shield, and was absorbed.  This path is done.
!
        else if ( absorb ( seed ) ) then

          na = na + 1
          ea = ea + e
          sa = sa + e * e
          exit
!
!  The particle collided with the shield and was scattered.
!  Find the scattering angle and energy, and continue along the new path.
!
        else

          call scatter ( seed, e, mu, azm )

        end if

      end do

    end do
!
!  Print the results of the simulation.
!
    call output ( na, ea, sa, nr, er, sr, nt, et, st, ntot )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REACTOR_SIMULATION:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function absorb ( seed )

!*****************************************************************************80
!
!! ABSORB determines if a colliding particle is absorbed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, logical ABSORB, is TRUE if the particle is absorbed.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) PA, the probability of absorption.
!
  implicit none

  logical absorb
  real ( kind = 8 ), parameter :: pa = 0.1D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u

  u = r8_uniform_01 ( seed )

  if ( u <= pa ) then
    absorb = .true.
  else
    absorb = .false.
  end if

  return
end
function cross ( e )

!*****************************************************************************80
!
!! CROSS returns the "cross section" of a particle based on its energy.
!
!  Discussion:
!
!    The particle's cross section is a measure of its likelihood to collide
!    with the material of the slab.  This quantity typically depends on both
!    the particle's energy and the kind of medium through which it is traveling.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) E, the energy of the particle.
!
!    Output, real ( kind = 8 ) CROSS, the cross section.
!
  implicit none

  real ( kind = 8 ) cross
  real ( kind = 8 ) e
  real ( kind = 8 ) s
  real ( kind = 8 ) y

  s = abs ( sin ( 100.0D+00 * ( exp ( e ) - 1.0D+00 ) ) &
    + sin ( 18.81D+00 * ( exp ( e ) - 1.0D+00 ) ) )

  y = max ( 0.02D+00, s )

  cross = 10.0D+00 * exp ( -0.1D+00 / y )

  return
end
function dist2c ( e, seed )

!*****************************************************************************80
!
!! DIST2C returns the distance to collision.
!
!  Discussion:
!
!    Assuming the particle has a given energy, and assuming it is currently
!    somewhere inside the shield, it is possible to determine a typical distance
!    which the particle can travel before it collides with the material of
!    the shield.
!
!    The computation of the collision distance is made by estimating a
!    "cross section" (as though having more energy made the particle "bigger"
!    and hence more likely to collide) and then randomly selecting a distance
!    that is logarithmically distributed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) E, the energy of the particle.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) DIST2C, the distance the particle can travel
!    through the slab before colliding.
!
  implicit none

  real ( kind = 8 ) cross
  real ( kind = 8 ) dist2c
  real ( kind = 8 ) e
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u

  u = r8_uniform_01 ( seed )

  dist2c = - log ( u ) / cross ( e )

  return
end
function energy ( seed )

!*****************************************************************************80
!
!! ENERGY assigns an energy to an emitted particle.
!
!  Discussion:
!
!    The energy E is in the range [EMIN,EMAX], with distribution
!    const/sqrt(energy).
!
!    An inverse function approach is used to compute this.
!
!    The energies are measured in MeV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) ENERGY, a randomly chosen energy that is
!    distributed as described above.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) EMIN, EMAX, the minimum and maximum
!    energies.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: emax = 2.5D+00
  real ( kind = 8 ), parameter :: emin = 1.0D-03
  real ( kind = 8 ) energy
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u

  u = r8_uniform_01 ( seed )

  c = 1.0D+00 / ( 2.0D+00 * ( sqrt ( emax ) - sqrt ( emin ) ) )

  energy = ( u / ( 2.0D+00 * c ) + sqrt ( emin ) )
  energy = energy * energy

  return
end
subroutine output ( na, ea, sa, nr, er, sr, nt, et, st, ntot )

!*****************************************************************************80
!
!! OUTPUT prints the results of the reactor shielding simulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, number of particles absorbed by the slab.
!
!    Input, real ( kind = 8 ) EA, energy absorbed by the slab.
!
!    Input, real ( kind = 8 ) SA, the sum of the squares of the 
!    absorbed energies.
!
!    Input, integer ( kind = 4 ) NR, number of particles reflected by the slab.
!
!    Input, real ( kind = 8 ) ER, energy reflected by the slab.
!
!    Input, real ( kind = 8 ) SR, the sum of the squares of the 
!    reflected energies.
!
!    Input, integer ( kind = 4 ) NT, number of particles transmitted by the slab.
!
!    Input, real ( kind = 8 ) ET, energy transmitted through the slab.
!
!    Input, real ( kind = 8 ) ST, the sum of the squares of the 
!    transmitted energies.
!
!    Input, integer ( kind = 4 ) NTOT, the total number of particles.
!
  implicit none

  real ( kind = 8 ) ea
  real ( kind = 8 ) ea_ave
  real ( kind = 8 ) er
  real ( kind = 8 ) er_ave
  real ( kind = 8 ) et
  real ( kind = 8 ) et_ave
  real ( kind = 8 ) etot
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) ntot
  real ( kind = 8 ) pa
  real ( kind = 8 ) pr
  real ( kind = 8 ) pt
  real ( kind = 8 ) ptot
  real ( kind = 8 ) sa
  real ( kind = 8 ) sr
  real ( kind = 8 ) st

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Reactor Shielding Problem:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '                           Total                   Average'
  write ( *, '(a)' ) '                    #      Energy      ' &
    // 'Percent     Energy         StDev'
  write ( *, '(a)' ) ' '

  etot = ea + er + et

  if ( 0 < na ) then
    ea_ave = ea / real ( na, kind = 8 )
    sa = sqrt ( sa / real ( na, kind = 8 ) - ea_ave * ea_ave )
  else
    ea_ave = 0.0D+00
  end if

  pa = real ( na * 100, kind = 8 ) / real ( ntot, kind = 8 )

  write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' ) &
    'Absorbed   ', na, ea, pa, ea_ave, sa

  if ( 0 < nr ) then
    er_ave = er / real ( nr, kind = 8 )
    sr = sqrt ( sr / real ( nr, kind = 8 ) - er_ave * er_ave )
  else
    er_ave = 0.0D+00
  end if

  pr = real ( nr * 100, kind = 8 ) / real ( ntot, kind = 8 )

  write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  &
    'Reflected  ', nr, er, pr, er_ave, sr

  if ( 0 < nt ) then
    et_ave = et / real ( nt, kind = 8 )
    st = sqrt ( st / real ( nt, kind = 8 ) - et_ave * et_ave )
  else
    et_ave = 0.0D+00
  end if

  pt = real ( nt * 100, kind = 8 ) / real ( ntot, kind = 8 )

  write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  &
    'Transmitted', nt, et, pt, et_ave, st

  ptot = 100.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  &
    'Total      ', ntot, etot, ptot

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine scatter ( seed, e, mu, azm )

!*****************************************************************************80
!
!! SCATTER returns the new direction and energy of a particle that is scattered.
!
!  Discussion:
!
!    The scattering direction is chosen uniformly on the sphere.
!
!    The energy of the scattered particle is chosen uniformly in
!    [ 0.3*E, E ].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Input/output, real ( kind = 8 ) E.  On input, the particle energy
!    before collision.  On output, the particle energy after collision
!    and scattering.
!
!    Output, real ( kind = 8 ) MU, the cosine of the angle between the
!    particle's direction and the X axis.
!
!    Output, real ( kind = 8 ) AZM, the azimuthal angle of the particle's
!    direction.
!
  implicit none

  real ( kind = 8 ) azm
  real ( kind = 8 ) e
  real ( kind = 8 ) mu
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u

  u = r8_uniform_01 ( seed )
  mu = - 1.0D+00 + 2.0D+00 * u

  u = r8_uniform_01 ( seed )
  azm = u * 2.0D+00 * pi

  u = r8_uniform_01 ( seed )
  e = ( u * 0.7D+00 + 0.3D+00 ) * e

  return
end
subroutine source ( seed, e, mu, azm, x, y, z )

!*****************************************************************************80
!
!! SOURCE generates a new particle from the neutron source.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) E, the initial energy of the particle.
!
!    Output, real ( kind = 8 ) MU, the cosine of the angle between the
!    particle's direction and the X axis.
!
!    Output, real ( kind = 8 ) AZM, the azimuthal angle of the particle's
!    direction.
!
!    Output, real ( kind = 8 ) X, Y, Z, the initial coordinates of the particle.
!
  implicit none

  real ( kind = 8 ) azm
  real ( kind = 8 ) e
  real ( kind = 8 ) energy
  real ( kind = 8 ) mu
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  u = r8_uniform_01 ( seed )
  mu = u

  u = r8_uniform_01 ( seed )
  azm = u * 2.0D+00 * pi

  x = 0.0D+00
  y = 0.0D+00
  z = 0.0D+00

  e = energy ( seed )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
subroutine update ( mu, azm, d, x, y, z )

!*****************************************************************************80
!
!! UPDATE determines the position of the particle after it has traveled D units.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) MU, the cosine of the angle between the
!    particle's direction and the X axis.
!
!    Input, real ( kind = 8 ) AZM, the azimuthal angle of the particle's
!    direction.
!
!    Input, real ( kind = 8 ) D, the distance the particle traveled.
!
!    Input/output, real ( kind = 8 ) X, Y, Z.  On input, the previous
!    coordinates of the particle.  On output, the updated coordinates of the
!    particle.
!
  implicit none

  real ( kind = 8 ) azm
  real ( kind = 8 ) d
  real ( kind = 8 ) mu
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  s = sqrt ( 1.0D+00 - mu * mu )

  x = x + d * mu
  y = y + d * s * cos ( azm )
  z = z + d * s * sin ( azm )

  return
end
