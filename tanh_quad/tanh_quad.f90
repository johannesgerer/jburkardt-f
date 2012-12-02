subroutine midpoint_rule ( n, x, w )

!*****************************************************************************80
!
!! MIDPOINT_RULE computes a variant of the midpoint quadrature rule.
!
!  Discussion:
!
!    2*N+1 equally spaced points in [-1,1] are defined, but the endpoints
!    are not used in the rule.
!
!    This rule is useful in cases where there is a singularity at one
!    or both endpoints of the interval.
!
!    This rule will nest, with N = 0, 1, 3, 7, 15, 31, ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the quadrature order.
!
!    Output, real ( kind = 8 ) X(-N:N), the abscissas.
!
!    Output, real ( kind = 8 ) W(-N:N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real    ( kind = 8 ) w(-n:n)
  real    ( kind = 8 ) x(-n:n)

  do i = -n, n
    x(i) = real ( i, kind = 8 ) / real ( n + 1, kind = 8 )
  end do

  w(-n:n) =         2.0D+00 / real ( 2 * ( n + 1 ), kind = 8 )
!
!  If I correct W(-N) and W(N) in this slightly awkward way,
!  I correctly include the special case where N = 0.
!
  w(-n) = w(-n) + 1.0D+00 / real ( 2 * ( n + 1 ), kind = 8 )
  w( n) = w( n) + 1.0D+00 / real ( 2 * ( n + 1 ), kind = 8 )

  return
end
subroutine rule_adjust ( a, b, c, d, order, x, w )

!*****************************************************************************80
!
!! RULE_ADJUST maps a quadrature rule from [A,B] to [C,D].
!
!  Discussion:
!
!    Most quadrature rules are defined on a special interval, like
!    [-1,1] or [0,1].  To integrate over an interval, the abscissas
!    and weights must be adjusted.  This can be done on the fly,
!    or by calling this routine.
!
!    If the weight function W(X) is not 1, then the weight vector W will
!    require further adjustment by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the definition interval.
!
!    Input, real ( kind = 8 ) C, D, the endpoints of the integration interval.
!
!    Input, integer ( kind = 4 ) ORDER, the number of abscissas and weights.
!
!    Input/output, real ( kind = 8 ) X(ORDER), W(ORDER), the abscissas
!    and weights.
!
  implicit none

  integer ( kind = 4 ) order

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) w(order)
  real    ( kind = 8 ) x(order)

  x(1:order) = ( ( b - x(1:order)     ) * c   &
               + (     x(1:order) - a ) * d ) &
               / ( b              - a )

  w(1:order) = ( ( d - c ) / ( b - a ) ) * w(1:order)

  return
end
subroutine tanh_h_to_n ( h, tol, n )

!*****************************************************************************80
!
!! TANH_H_TO_N computes N as a function of H and TOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H, the spacing.
!
!    Input, real ( kind = 8 ) TOL, the tolerance.
!
!    Output, integer ( kind = 4 ) N, the corresponding quadrature order.
!
  implicit none

  real    ( kind = 8 ) ct
  real    ( kind = 8 ) h
  integer ( kind = 4 ) n
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) w

  n = 0

  do

    t = real ( n, kind = 8 ) * h / 2.0D+00

    ct = cosh ( t )

    w = 0.5D+00 * h / ct / ct

    if ( w <= tol ) then
      exit
    end if

    n = n + 1

  end do

  return
end
subroutine tanh_m_to_h ( m, h )

!*****************************************************************************80
!
!! TANH_M_TO_H computes H as a function of M.
!
!  Discussion:
!
!    H = 2^(-M).
!
!    This is simply an orderly way to index a family of decreasing values of H.
!
!  Example:
!
!     M      H
!    --  -----
!    -2      4
!    -1      2
!     0      1
!     1     1/2
!     2     1/4
!     3     1/8
!     4     1/16
!     5     1/32
!   ...     ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the level.
!
!    Output, real ( kind = 8 ) H, the spacing.
!
  implicit none

  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m

  h = 1.0D+00

  do i = -1, m, -1
    h = h * 2.0D+00
  end do

  do i = 1, m
    h = h / 2.0D+00
  end do

  return
end
subroutine tanh_n_to_h ( n, h  )

!*****************************************************************************80
!
!! TANH_N_TO_H computes N as a function of H for the tanh rule.
!
!  Discussion:
!
!    This formula for N(H) is suggested in Kahaner, Moler and Nash.
!    Note, however, that using this formula means that it is not possible
!    to make families of nested rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) N, the quadrature order.
!
!    Output, real ( kind = 8 ) H, the spacing.
!
  implicit none

  real    ( kind = 8 ) h
  integer ( kind = 4 ) n
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  h = pi * sqrt ( 2.0D+00 / real ( n, kind = 8 ) ) &
    - 1.0D+00 / real ( n, kind = 8 )

  return
end
subroutine tanh_rule ( n, h, x, w )

!*****************************************************************************80
!
!! TANH_RULE computes a tanh quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the quadrature order.
!
!    Input, real ( kind = 8 ) H, the spacing.
!
!    Output, real ( kind = 8 ) X(-N:N), the abscissas.
!
!    Output, real ( kind = 8 ) W(-N:N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) ct
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) t
  real    ( kind = 8 ) w(-n:n)
  real    ( kind = 8 ) x(-n:n)

  do i = -n, n

    t = real ( i, kind = 8 ) * h / 2.0D+00

    ct = cosh ( t )

    x(i) = tanh ( t )
    w(i) = 0.5D+00 * h / ct / ct

  end do

  return
end
subroutine tanh_sinh_h_to_n ( h, tol, n )

!*****************************************************************************80
!
!! TANH_SINH_H_TO_N computes N as a function of H and TOL for tanh-sinh rule.
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
!    Input, real ( kind = 8 ) H, the spacing.
!
!    Input, real ( kind = 8 ) TOL, the tolerance.
!
!    Output, integer ( kind = 4 ) N, the corresponding quadrature order.
!
  implicit none

  real    ( kind = 8 ) ct
  real    ( kind = 8 ) ct2
  real    ( kind = 8 ) h
  integer ( kind = 4 ) n
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) st
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) w

  n = 0

  do

    t = real ( n, kind = 8 ) * h

    ct = cosh ( t )
    st = sinh ( t )
    ct2 = cosh ( 0.5D+00 * pi * st )

    w = 0.5D+00 * pi * h * ct / ct2 / ct2

    if ( w <= tol ) then
      exit
    end if

    n = n + 1

  end do

  return
end
subroutine tanh_sinh_rule ( n, h, x, w )

!*****************************************************************************80
!
!! TANH_SINH_RULE computes a tanh-sinh quadrature rule.
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
!    Input, integer ( kind = 4 ) N, the quadrature order.
!
!    Input, real ( kind = 8 ) H, the spacing.
!
!    Output, real ( kind = 8 ) X(-N:N), the abscissas.
!
!    Output, real ( kind = 8 ) W(-N:N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) ct
  real    ( kind = 8 ) ct2
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) st
  real    ( kind = 8 ) t
  real    ( kind = 8 ) w(-n:n)
  real    ( kind = 8 ) x(-n:n)

  do i = -n, n

    t = real ( i, kind = 8 ) * h

    ct = cosh ( t )
    st = sinh ( t )
    ct2 = cosh ( 0.5D+00 * pi * st )

    x(i) = tanh ( 0.5D+00 * pi * st )

    w(i) = 0.5D+00 * pi * h * ct / ct2 / ct2

  end do

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
subroutine trap_rule ( n, x, w )

!*****************************************************************************80
!
!! TRAP_RULE computes a trapezoid quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the quadrature order.
!
!    Output, real ( kind = 8 ) X(-N:N), the abscissas.
!
!    Output, real ( kind = 8 ) W(-N:N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real    ( kind = 8 ) w(-n:n)
  real    ( kind = 8 ) x(-n:n)

  do i = -n, n
    x(i) = real ( i, kind = 8 ) / real ( n, kind = 8 )
  end do

  w(-n)       = 1.0D+00
  w(-n+1:n-1) = 2.0D+00
  w(+n)       = 1.0D+00

  w(-n:n) = w(-n:n) / real ( 2 * n, kind = 8 )

  return
end
