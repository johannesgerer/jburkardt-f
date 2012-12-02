subroutine shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi, zi )

!*****************************************************************************80
!
!! SHEPARD_INTERP_2D evaluates a 2D Shepard interpolant.
!
!  Discussion:
!
!    This code should be vectorized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Shepard,
!    A two-dimensional interpolation function for irregularly spaced data,
!    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
!    ACM, pages 517-524, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), YD(ND), the data points.
!
!    Input, real ( kind = 8 ) ZD(ND), the data values.
!
!    Input, real ( kind = 8 ) P, the power.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), YI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ) s
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)
  integer ( kind = 4 ) z
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do i = 1, ni

    if ( p == 0.0D+00 ) then

      w(1:nd) = 1.0D+00 / real ( nd, kind = 8 )

    else

      z = -1
      do j = 1, nd
        w(j) = sqrt ( ( xi(i) - xd(j) ) ** 2 + ( yi(i) - yd(j) ) ** 2 )
        if ( w(j) == 0.0D+00 ) then
          z = j
          exit
        end if
      end do

      if ( z /= -1 ) then
        w(1:nd) = 0.0D+00
        w(z) = 1.0D+00
      else
        w(1:nd) = 1.0D+00 / w(1:nd) ** p
        s = sum ( w )
        w(1:nd) = w(1:nd) / s
      end if

    end if

    zi(i) = dot_product ( w, zd )

  end do

  return
end
