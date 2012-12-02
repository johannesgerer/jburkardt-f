subroutine asset_path ( s0, mu, sigma, t1, n, seed, s )

!*****************************************************************************80
!
!! ASSET_PATH simulates the behavior of an asset price over time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2012
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) S0, the asset price at time 0.
!
!    Input, real ( kind = 8 ) MU, the expected growth rate.
!
!    Input, real ( kind = 8 ) SIGMA, the volatility of the asset.
!
!    Input, real ( kind = 8 ) T1, the expiry date.
!
!    Input, integer ( kind = 4 ) N, the number of steps to take between 0 and T1.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) S(0:N), the option values from time 0 to T1 
!    in equal steps.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dt
  integer ( kind = 4 ) i
  real ( kind = 8 ) mu
  real ( kind = 8 ) p
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) s(0:n)
  real ( kind = 8 ) s0
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t1

  dt = t1 / real ( n, kind = 8 )

  call r8vec_normal_01 ( n, seed, r )

  s(0) = s0
  p = s0
  do i = 1, n
    p = p * exp ( ( mu - sigma * sigma ) * dt + sigma * sqrt ( dt ) * r(i) )
    s(i) = p
  end do

  return
end
subroutine binomial ( s0, e, r, sigma, t1, m, c )

!*****************************************************************************80
!
!! BINOMIAL uses the binomial method for a European call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2012
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) S0, the asset price at time 0.
!
!    Input, real ( kind = 8 ) E, the exercise price.
!
!    Input, real ( kind = 8 ) R, the interest rate.
!
!    Input, real ( kind = 8 ) SIGMA, the volatility of the asset.
!
!    Input, real ( kind = 8 ) T1, the expiry date.
!
!    Input, integer ( kind = 4 ) M, the number of steps to take 
!    between 0 and T1.
!
!    Output, real ( kind = 8 ) C, the option value at time 0.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dt
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) s0
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t1
  real ( kind = 8 ) u
  real ( kind = 8 ) w(1:m+1)
!
!  Time stepsize.
!
  dt = t1 / real ( m, kind = 8 )

  a = 0.5D+00 * ( exp ( - r * dt ) + exp ( ( r + sigma**2 ) * dt ) )

  d = a - sqrt ( a * a - 1.0D+00 )
  u = a + sqrt ( a * a - 1.0D+00 )

  p = ( exp ( r * dt ) - d ) / ( u - d )

  do i = 1, m + 1
    w(i) = max ( s0 * d**(m+1-i) * u**(i-1) - e, 0.0D+00 )
  end do
!
!  Trace backwards to get the option value at time 0.
!
  do n = m, 1, -1
    do i = 1, n
      w(i) = ( 1.0D+00 - p ) * w(i) + p * w(i+1)
    end do
  end do

  w(1:m+1) = exp ( - r * t1 ) * w(1:m+1)

  c = w(1)

  return
end
subroutine bsf ( s0, t0, e, r, sigma, t1, c )

!*****************************************************************************80
!
!! BSF evaluates the Black-Scholes formula for a European call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2012
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) S0, the asset price at time T0.
!
!    Input, real ( kind = 8 ) T0, the time at which the asset price is known.
!
!    Input, real ( kind = 8 ) E, the exercise price.
!
!    Input, real ( kind = 8 ) R, the interest rate.
!
!    Input, real ( kind = 8 ) SIGMA, the volatility of the asset.
!
!    Input, real ( kind = 8 ) T1, the expiry date.
!
!    Output, real ( kind = 8 ) C, the value of the call option.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) e
  real ( kind = 8 ) n1
  real ( kind = 8 ) n2
  real ( kind = 8 ) r
  real ( kind = 8 ) s0
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ) tau

  tau = t1 - t0

  if ( 0.0D+00 < tau ) then

    d1 = ( log ( s0 / e ) + ( r + 0.5D+00 * sigma * sigma ) * tau ) &
      / ( sigma * sqrt ( tau ) )

    d2 = d1 - sigma * sqrt ( tau )

    n1 = 0.5D+00 * ( 1.0D+00 + erf ( d1 / sqrt ( 2.0D+00 ) ) )
    n2 = 0.5D+00 * ( 1.0D+00 + erf ( d2 / sqrt ( 2.0D+00 ) ) )

    c = s0 * n1 - e * exp ( - r * tau ) * n2

  else

    c = max ( s0 - e, 0.0D+00 )

  end if

  return
end
subroutine forward ( e, r, sigma, t1, nx, nt, smax, u )

!*****************************************************************************80
!
!! FORWARD uses the forward difference method to value a European call option.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2012
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) E, the exercise price.
!
!    Input, real ( kind = 8 ) R, the interest rate.
!
!    Input, real ( kind = 8 ) SIGMA, the volatility of the asset.
!
!    Input, real ( kind = 8 ) T1, the expiry date.
!
!    Input, integer ( kind = 4 ) NX, the number of "space" steps used to 
!    divide the interval [0,L].
!
!    Input, integer ( kind = 4 ) NT, the number of time steps.
!
!    Input, real ( kind = 8 ) SMAX, the maximum value of S to consider.
!
!    Output, real ( kind = 8 ) U(NX-1,NT+1), the value of the European 
!    call option.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nx

  real ( kind = 8 ) a(2:nx-1)
  real ( kind = 8 ) b(1:nx-1)
  real ( kind = 8 ) c(1:nx-2)
  real ( kind = 8 ) dt
  real ( kind = 8 ) dx
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) sigma
  real ( kind = 8 ) smax
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) u(nx-1,nt+1)
  real ( kind = 8 ) u0

  dt = t1 / real ( nt, kind = 8 )
  dx = smax / real ( nx, kind = 8 )

  do i = 1, nx - 1
    b(i) = 1.0D+00 - r * dt - dt * ( sigma * i )**2
  end do

  do i = 1, nx - 2
    c(i) = 0.5D+00 * dt * ( sigma * i )**2 + 0.5D+00 * dt * r * i
  end do

  do i = 2, nx - 1
    a(i) = 0.5D+00 * dt * ( sigma * i )**2 - 0.5D+00 * dt * r * i
  end do

  u0 = 0.0D+00
  do i = 1, nx - 1
    u0 = u0 + dx
    u(i,1) = max ( u0 - e, 0.0D+00 )
  end do
  
  do j = 1, nt

    t = real ( j - 1, kind = 8 ) * t1 / real ( nt, kind = 8 )

    p = 0.5D+00 * dt * ( nx - 1 ) * ( sigma * sigma * ( nx - 1 ) + r ) &
      * ( smax - e * exp ( - r * t ) )

    u(1:nx-1,j+1) =                 b(1:nx-1) * u(1:nx-1,j)
    u(1:nx-2,j+1) = u(1:nx-2,j+1) + c(1:nx-2) * u(2:nx-1,j)
    u(2:nx-1,j+1) = u(2:nx-1,j+1) + a(2:nx-1) * u(1:nx-2,j)

    u(nx-1,j+1) = u(nx-1,j+1) + p

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine mc ( s0, e, r, sigma, t1, m, seed, conf )

!*****************************************************************************80
!
!! MC uses Monte Carlo valuation on a European call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) S0, the asset price at time 0.
!
!    Input, real ( kind = 8 ) E, the exercise price.
!
!    Input, real ( kind = 8 ) R, the interest rate.
!
!    Input, real ( kind = 8 ) SIGMA, the volatility of the asset.
!
!    Input, real ( kind = 8 ) T1, the expiry date.
!
!    Input, integer ( kind = 4 ) M, the number of simulations.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) CONF(2), the estimated range of the valuation.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) conf(2)
  real ( kind = 8 ) e
  real ( kind = 8 ) pmean
  real ( kind = 8 ) pvals(m)
  real ( kind = 8 ) r
  real ( kind = 8 ) s0
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sigma
  real ( kind = 8 ) std
  real ( kind = 8 ) svals(m)
  real ( kind = 8 ) t1
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) width

  call r8vec_normal_01 ( m, seed, u )

  svals(1:m) = s0 * exp ( ( r - 0.5D+00 * sigma * sigma ) * t1 &
    + sigma * sqrt ( t1 ) * u(1:m) )

  pvals(1:m) = exp ( - r * t1 ) * max ( svals(1:m) - e, 0.0D+00 )

  pmean = sum ( pvals(1:m) ) / real ( m, kind = 8 )

  std = sqrt ( sum ( ( pvals(1:m) - pmean )**2 ) / real ( m - 1, kind = 8 ) )

  width = 1.96D+00 * std / sqrt ( real ( m, kind = 8 ) )

  conf(1) = pmean - width
  conf(2) = pmean + width

  return
end
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!    This means that on every other call, essentially, the input value of
!    SEED is ignored, since the code saves the second normal random value.
!
!    If you didn't know this, you might be confused since, usually, the
!    output of a random number generator can be completely controlled by
!    the input value of the SEED.  If I were more careful, I could rewrite
!    this routine so that it would distinguish between cases where the input
!    value of SEED is the output value from the previous call (all is well)
!    and those cases where it is not (the user has decided to do something
!    new.  Restart the uniform random number sequence.)  But I'll leave
!    that for later.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a sample of the standard
!    normal PDF.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, 2 ) == 0 ) then

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r8_uniform_01 ( seed2 )

    x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( - 2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r8_normal_01 = x

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
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
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
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) SAVED, is 0 or 1 depending on whether there
!    is a single saved value left over from the previous call.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute.  This starts off as 1:N, but
!    is adjusted if we have a saved value that can be immediately stored
!    in X(1), and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: made = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( - 2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine r8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_PART prints "part" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_write ( output_filename, n, x )

!*****************************************************************************80
!
!! R8VEC_WRITE writes an R8VEC file.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the data.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) x(n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if

  if ( 0 < n ) then
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, '(2x,g24.16)' ) x(j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

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
