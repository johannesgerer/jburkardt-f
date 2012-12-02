subroutine lagrange_approx_1d ( m, nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! LAGRANGE_APPROX_1D evaluates the Lagrange approximant of degree M.
!
!  Discussion:
!
!    The Lagrange approximant L(M,ND,XD,YD)(X) is a polynomial of
!    degree M which approximates the data (XD(I),YD(I)) for I = 1 to ND.
!
!    We can represent any polynomial of degree M+1 as the sum of the Lagrange 
!    basis functions at the M+1 Chebyshev points.
!
!      L(M)(X) = sum ( 1 <= I <= M+1 ) C(I) LB(M,XC)(X)
!
!    Given our data, we can seek the M+1 unknown coefficients C which minimize
!    the norm of || L(M)(XD(1:ND)) - YD(1:ND) ||.
!
!    Given the coefficients, we can then evaluate the polynomial at the
!    points XI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the polynomial degree.
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

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ld(nd,m+1)
  real ( kind = 8 ) li(ni,m+1)
  integer ( kind = 4 ) nc
  real ( kind = 8 ) xc(m+1)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yc(m+1)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  nc = m + 1
!
!  Evaluate the Chebyshev points.
!
  a = -1.0D+00
  b = +1.0D+00
  call r8vec_chebyspace ( nc, a, b, xc )
!
!  Evaluate the Lagrange basis functions for the Chebyshev points 
!  at the data points.
!
  call lagrange_basis_1d ( nc, xc, nd, xd, ld )
!
!  The value of the Lagrange approximant at each data point should
!  approximate the data value: LD * YC = YD, where YC are the unknown
!  coefficients.
!
  call qr_solve ( nd, nc, ld, yd, yc )
!
!  Now we want to evaluate the Lagrange approximant at the "interpolant
!  points": LI * YC = YI
!
  call lagrange_basis_1d ( nc, xc, ni, xi, li )
  yi(1:ni) = matmul ( li(1:ni,1:nc), yc(1:nc) )

  return
end
subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 

!*****************************************************************************80
!
!! LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the interpolation nodes.
!
!    Input, integer ( kind = 4 ) NI, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XI(NI), the evaluation points.
!
!    Output, real ( kind = 8 ) LB(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lb(ni,nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  
  do i = 1, ni
    do j = 1, nd
      lb(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end
