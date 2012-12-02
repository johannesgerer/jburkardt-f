subroutine chebyshev_coef_1d ( nd, xd, yd, c, xmin, xmax )

!*****************************************************************************80
!
!! CHEBYSHEV_COEF_1D determines the Chebyshev interpolant coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data locations.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Output, real ( kind = 8 ) C(ND), the Chebyshev coefficients.
!
!    Output, real ( kind = 8 ) XMIN, XMAX, the interpolation interval.
!
  implicit none

  integer ( kind = 4 ) nd

  real ( kind = 8 ) a(nd,nd)
  real ( kind = 8 ) c(nd)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) yd(nd)

  if ( nd == 1 ) then
    xmin = xd(1)
    xmax = xd(1)
    c(1) = 1.0D+00
    return
  end if

  xmin = minval ( xd(1:nd) )
  xmax = maxval ( xd(1:nd) )
!
!  Map XD to [-1,+1].
!
  x(1:nd) = ( 2.0D+00 * xd(1:nd) - xmin - xmax ) / ( xmax - xmin )
!
!  Form the Chebyshev Vandermonde matrix.
!
  do j = 1, nd
    do i = 1, nd
      a(i,j) = cos ( acos ( x(i) ) * real ( j - 1, kind = 8 ) )
    end do
  end do 
!
!  Solve for the expansion coefficients.
!
  call qr_solve ( nd, nd, a, yd, c )

  return
end
subroutine chebyshev_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! CHEBYSHEV_INTERP_1D determines and evaluates the Chebyshev interpolant.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data locations.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points, which
!    must be each be in the interval [ min(XD), max(XD)].
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) c(nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  call chebyshev_coef_1d ( nd, xd, yd, c, xmin, xmax )

  call chebyshev_value_1d ( nd, c, xmin, xmax, ni, xi, yi )

  return
end
subroutine chebyshev_value_1d ( nd, c, xmin, xmax, ni, xi, yi )

!*****************************************************************************80
!
!! CHEBYSHEV_VALUE_1D evaluates a Chebyshev interpolant, given its coefficients.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) C(ND), the Chebyshev coefficients.
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the interpolation interval.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points, which
!    must be each be in the interval [XMIN,XMAX].
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a(ni,nd)
  real ( kind = 8 ) c(nd)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(ni)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) yi(ni)

  if ( nd == 1 ) then
    yi(1) = c(1)
    return
  end if
!
!  Map XI to [-1,+1].
!
  x(1:ni) = ( 2.0D+00 * xi(1:ni) - xmin - xmax ) / ( xmax - xmin )

  do j = 1, nd
    do i = 1, ni
      a(i,j) = cos ( acos ( x(i) ) * real ( j - 1, kind = 8 ) )
    end do
  end do 

  yi(1:ni) = matmul ( a(1:ni,1:nd), c(1:nd) )

  return
end
