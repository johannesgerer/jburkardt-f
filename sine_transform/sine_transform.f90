subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

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
subroutine sine_transform_data ( n, d, s )

!*****************************************************************************80
!
!! SINE_TRANSFORM_DATA does a sine transform on a vector of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, real ( kind = 8 ) D(N), the vector of data.
!
!    Output, real ( kind = 8 ) S(N), the sine transform coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) s(n)

  do i = 1, n
    s(i) = 0.0D+00
    do j = 1, n
      angle = pi * real ( i * j, kind = 8 ) / real ( n + 1, kind = 8 )
      s(i) = s(i) + sin ( angle ) * d(j)
    end do
    s(i) = s(i) * sqrt ( 2.0D+00 / real ( n + 1, kind = 8 ) )
  end do

  return
end
subroutine sine_transform_function ( n, a, b, f, s )

!*****************************************************************************80
!
!! SINE_TRANSFORM_FUNCTION does a sine transform on functional data.
!
!  Discussion:
!
!    The interval [A,B] is divided into N+1 intervals using N+2 points,
!    which are indexed by 0 through N+1.
!
!    The original function F(X) is regarded as the sum of a linear function 
!    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
!    which is 0 at A and B.
!
!    The sine transform coefficients for F2 are then computed.
!
!    To recover the interpolant of F(X), it is necessary to combine the
!    linear part F1 with the sine transform interpolant:
!
!      Interp(F)(X) = F1(X) + F2(X)
!
!    This can be done by calling SINE_TRANSFORM_INTERPOLANT().
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
!    Input, external, real ( kind = 8 ) F, a pointer to the function.
!
!    Output, real ( kind = 8 ) S(N), the sine transform coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) angle
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) f2(n)
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) x(n)
!
!  Evenly spaced points between A and B, but omitting
!  A and B themselves.
!
  do i = 1, n
    x(i) = ( real ( n - i + 1, kind = 8 ) * a   &
           + real (     i,     kind = 8 ) * b ) &
           / real ( n     + 1, kind = 8 )
  end do
!
!  Subtract F1(X) from F(X) to get F2(X).
!
  fa = f ( a )
  fb = f ( b )

  do i = 1, n
    f2(i) = f ( x(i) )                &
          - ( ( b - x(i)     ) * fa   &
          +   (     x(i) - a ) * fb ) &
          /   ( b        - a )
  end do
!
!  Compute the sine transform of F2(X).
!
  do i = 1, n
    s(i) = 0.0D+00
    do j = 1, n
      angle = pi * real ( i * j, kind = 8 ) / real ( n + 1, kind = 8 )
      s(i) = s(i) + sin ( angle ) * f2(j)
    end do
    s(i) = s(i) * sqrt ( 2.0D+00 / real ( n + 1, kind = 8 ) )
  end do

  return
end
subroutine sine_transform_interpolant ( n, a, b, fa, fb, s, nx, x, value )

!*****************************************************************************80
!
!! SINE_TRANSFORM_INTERPOLANT evaluates the sine transform interpolant.
!
!  Discussion:
!
!    The interval [A,B] is divided into N+1 intervals using N+2 points,
!    which are indexed by 0 through N+1.
!
!    The original function F(X) is regarded as the sum of a linear function 
!    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
!    which is 0 at A and B.
!
!    The function F2 has been approximated using the sine transform,
!    and the interpolant is then evaluated as:
!
!      Interp(F)(X) = F1(X) + F2(X)
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of terms in the approximation.
!
!    Input, real ( kind = 8 ) A, B, the interval over which the approximant 
!    was defined.
!
!    Input, real ( kind = 8 ) FA, FB, the function values at A and B.
!
!    Input, real ( kind = 8 ) S(N), the approximant coefficients.
!
!    Input, integer ( kind = 4 ) NX, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(NX), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(NX), the value of the interpolant.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nx

  real ( kind = 8 ) a
  real ( kind = 8 ) angle
  real ( kind = 8 ) b
  real ( kind = 8 ) f1(nx)
  real ( kind = 8 ) f2(nx)
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) value(nx)
  real ( kind = 8 ) x(nx)
!
!  Compute linear function F1(X).
!
  f1(1:nx) = ( ( b - x(1:nx)     ) * fa   &
             + (     x(1:nx) - a ) * fb ) &
             / ( b           - a )
!
!  Compute sine interpolant F2(X).
!
  f2(1:nx) = 0.0D+00

  do i = 1, nx
    do j = 1, n
      angle = real ( j, kind = 8 ) * ( x(i) - a ) * pi / ( b - a )
      f2(i) = f2(i) + s(j) * sin ( angle )
    end do
  end do

  f2(1:nx) = f2(1:nx) * sqrt ( 2.0D+00 / real ( n + 1, kind = 8 ) )
!
!  Interpolant = F1 + F2.
! 
  value(1:nx) = f1(1:nx) + f2(1:nx)

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
