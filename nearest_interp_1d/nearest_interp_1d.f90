subroutine nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! NEAREST_INTERP_1D evaluates the nearest neighbor interpolant.
!
!  Discussion:
!
!    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
!    constant function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2012
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

  real ( kind = 8 ) d
  real ( kind = 8 ) d2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  do i = 1, ni

    k = 1
    d = abs ( xi(i) - xd(k) )

    do j = 2, nd

      d2 = abs ( xi(i) - xd(j) )

      if ( d2 < d ) then
        k = j
        d = d2
      end if

    end do

    yi(i) = yd(k)

  end do

  return
end
