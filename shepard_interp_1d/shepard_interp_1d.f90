subroutine shepard_basis_1d ( nd, xd, k, p, ni, xi, bk )

!*****************************************************************************80
!
!! SHEPARD_BASIS_1D evaluates a 1D Shepard basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2012
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
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, integer ( kind = 4 ) K, the index of the desired basis function,
!    1 <= K <= ND.
!
!    Input, real ( kind = 8 ) P, the power.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) BK(NI), the basis function at the interpolation 
!    points.
! 
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) bk(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) p
  real ( kind = 8 ) s
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  integer ( kind = 4 ) z

  do i = 1, ni

    if ( p == 0.0D+00 ) then

      w(1:nd) = 1.0D+00 / real ( nd, kind = 8 )

    else

      z = -1
      do j = 1, nd
        w(j) = abs ( xi(i) - xd(j) )
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
        s = sum ( w(1:nd) )
        w(1:nd) = w(1:nd) / s
      end if

    end if

    bk(i) = w(k)

  end do

  return
end
subroutine shepard_interp_1d ( nd, xd, yd, p, ni, xi, yi )

!*****************************************************************************80
!
!! SHEPARD_INTERP_1D evaluates a 1D Shepard interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2012
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
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, real ( kind = 8 ) P, the power.
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) p
  real ( kind = 8 ) s
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)
  integer ( kind = 4 ) z

  do i = 1, ni

    if ( p == 0.0D+00 ) then

      w(1:nd) = 1.0D+00 / real ( nd, kind = 8 )

    else

      z = 0
      do j = 1, nd
        w(j) = abs ( xi(i) - xd(j) )
        if ( w(j) == 0.0D+00 ) then
          z = j
          exit
        end if
      end do

      if ( z /= 0 ) then
        w(1:nd) = 0.0D+00
        w(z) = 1.0D+00
      else
        w(1:nd) = 1.0D+00 / w(1:nd) ** p
        s = sum ( w )
        w(1:nd) = w(1:nd) / s
      end if

    end if

    yi(i) = dot_product ( w(1:nd), yd(1:nd) )

  end do

  return
end
