subroutine lagcheby1_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! LAGCHEBY1_INTERP_1D evaluates the Lagrange Chebyshev 1 interpolant.
!
!  Discussion:
!
!    The weight vector WD computed below is only valid if the data points
!    XD are, as expected, the Chebyshev Type 1 points for [-1,+1], or a linearly 
!    mapped version for [A,B].  The XD values may be computed by:
!
!      xd = r8vec_cheby1space ( nd, a, b )
!
!    for instance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Paul Berrut, Lloyd Trefethen,
!    Barycentric Lagrange Interpolation,
!    SIAM Review,
!    Volume 46, Number 3, September 2004, pages 501-517.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) denom(ni)
  integer ( kind = 4 ) exact(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) numer(ni)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) t
  real ( kind = 8 ) theta
  real ( kind = 8 ) wd
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  exact(1:ni) = 0

  do j = 1, nd

    theta = real ( 2 * j - 1, kind = 8 ) * pi / real ( 2 * nd, kind = 8 )
    wd = r8_mop ( j + 1 ) * sin ( theta )

    do i = 1, ni

      if ( xi(i) == xd(j) ) then
        exact(i) = j
        numer(i) = yd(j)
        denom(i) = 1.0D+00
      end if

      if ( exact(i) == 0 ) then
        t = wd / ( xi(i) - xd(j) )
        numer(i) = numer(i) + t * yd(j)
        denom(i) = denom(i) + t
      end if

    end do
  end do

  yi(1:ni) = numer(1:ni) / denom(1:ni)

  return
end
subroutine lagcheby2_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! LAGCHEBY2_INTERP_1D evaluates the Lagrange Chebyshev 2 interpolant.
!
!  Discussion:
!
!    The weight vector WD computed below is only valid if the data points
!    XD are, as expected, the Chebyshev Type 2 points for [-1,+1], or a linearly 
!    mapped version for [A,B].  The XD values may be computed by:
!
!      xd = r8vec_cheby2space ( nd, a, b )
!
!    for instance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Paul Berrut, Lloyd Trefethen,
!    Barycentric Lagrange Interpolation,
!    SIAM Review,
!    Volume 46, Number 3, September 2004, pages 501-517.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) denom(ni)
  integer ( kind = 4 ) exact(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) numer(ni)
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) t
  real ( kind = 8 ) wd
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  exact(1:ni) = 0

  do j = 1, nd

    wd = r8_mop ( j + 1 )
    if ( j == 1 .or. j == nd ) then
      wd = 0.5D+00 * wd
    end if

    do i = 1, ni

      if ( xi(i) == xd(j) ) then
        exact(i) = j
        numer(i) = yd(j)
        denom(i) = 1.0D+00
      end if

      if ( exact(i) == 0 ) then
        t = wd / ( xi(i) - xd(j) )
        numer(i) = numer(i) + t * yd(j)
        denom(i) = denom(i) + t
      end if

    end do
  end do

  yi(1:ni) = numer(1:ni) / denom(1:ni)

  return
end
subroutine lageven_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! LAGEVEN_VALUE_1D evaluates the Lagrange evenly-spaced interpolant.
!
!  Discussion:
!
!    The weight vector WD computed below is only valid if the data points
!    XD are, as expected, evenly spaced in an interval [A,B] with
!    spacing (B-A)/N.  The XD values might be computed by:
!
!      xd(i) = ( ( 2 * nd - 2 * i + 1 ) * a 
!              + (          2 * i - 1 ) * b ) 
!              / ( 2 * nd             )
!
!    for instance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Paul Berrut, Lloyd Trefethen,
!    Barycentric Lagrange Interpolation,
!    SIAM Review,
!    Volume 46, Number 3, September 2004, pages 501-517.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) denom(ni)
  integer ( kind = 4 ) exact(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) numer(ni)
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) t
  real ( kind = 8 ) wd
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  exact(1:ni) = 0

  do j = 1, nd

    wd = r8_mop ( j ) * r8_choose ( nd, j )

    do i = 1, ni

      if ( xi(i) == xd(j) ) then
        exact(i) = j
        numer(i) = yd(j)
        denom(i) = 1.0D+00
      end if

      if ( exact(i) == 0 ) then
        t = wd / ( xi(i) - xd(j) )
        numer(i) = numer(i) + t * yd(j)
        denom(i) = denom(i) + t
      end if

    end do
  end do

  yi(1:ni) = numer(1:ni) / denom(1:ni)

  return
end
