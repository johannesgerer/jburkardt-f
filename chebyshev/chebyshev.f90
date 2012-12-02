subroutine chebyshev_coefficients ( a, b, n, f, c )

!*****************************************************************************80
!
!! CHEBYSHEV_COEFFICIENTS determines Chebyshev interpolation coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Second Edition,
!    Cambridge University Press, 1992,
!    ISBN: 0-521-43064-X,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the domain of definition.
!
!    Input, integer ( kind = 4 ) N, the order of the interpolant.
!
!    Input, real ( kind = 8 ), external :: F ( X ), an external function.
!
!    Output, real ( kind = 8 ) C(N), the Chebyshev coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) angle
  real ( kind = 8 ) b
  real ( kind = 8 ) c(n)
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  do i = 1, n
    angle = real ( 2 * i - 1, kind = 8 ) * pi / real ( 2 * n, kind = 8 )
    x = cos ( angle )
    x = 0.5D+00 * ( a + b ) + x * 0.5D+00 * ( b - a )
    fx(i) = f ( x );
  end do

  do i = 1, n
    c(i) = 0.0D+00
    do j = 1, n
      angle = real ( ( i - 1 ) * ( 2 * j - 1 ), kind = 8 ) * pi &
        / real ( 2 * n, kind = 8 )
      c(i) = c(i) + fx(j) * cos ( angle )
    end do
  end do

  c(1:n) = 2.0D+00 * c(1:n) / real ( n, kind = 8 )

  return
end
subroutine chebyshev_interpolant ( a, b, n, c, m, x, cf )

!*****************************************************************************80
!
!! CHEBYSHEV_INTERPOLANT evaluates a Chebyshev interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Second Edition,
!    Cambridge University Press, 1992,
!    ISBN: 0-521-43064-X,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the domain of definition.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) C(N), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) M, the number of points.
!
!    Input, real ( kind = 8 ) X(M), the point at which the polynomial is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CF(M), the value of the Chebyshev
!    polynomial at X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) cf(m)
  real ( kind = 8 ) di
  real ( kind = 8 ) dip1
  real ( kind = 8 ) dip2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) y

  do j = 1, m

    dip1 = 0.0D+00
    di = 0.0D+00
    y = ( 2.0D+00 * x(j) - a  - b ) / ( b - a )

    do i = n, 2, -1
      dip2 = dip1
      dip1 = di
      di = 2.0D+00 * y * dip1 - dip2 + c(i)
    end do

    cf(j) = y * di - dip1 + 0.5D+00 * c(1)

  end do

  return
end
subroutine chebyshev_zeros ( n, x )

!*****************************************************************************80
!
!! CHEBYSHEV_ZEROS returns zeroes of the Chebyshev polynomial T(N,X).
!
!  Discussion:
!
!    We produce the Chebyshev zeros in ascending order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X(N), the zeroes of T(N)(X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  do i = 1, n
    angle = real ( 2 * ( n - i ) + 1, kind = 8 ) * pi / real ( 2 * n, kind = 8 )
    x(i) = cos ( angle )
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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
