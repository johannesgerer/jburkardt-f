subroutine vandermonde_interp_1d_coef ( n, x, y, c )

!*****************************************************************************80
!
!! VANDERMONDE_INTERP_1D_COEF computes a 1D polynomial interpolant.
!
!  Discussion:
!
!    We assume the interpolant has the form
!
!      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
!
!    We have n data values (x(i),y(i)) which must be interpolated:
!
!      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
!
!    This can be cast as an NxN linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
!      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
!      [ ...................... ] [ ... ] = [ ... ]
!      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
!
!    and if the x values are distinct, the system is theoretically
!    invertible, so we can retrieve the coefficient vector c and
!    evaluate the interpolant.
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
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the data values.
!
!    Output, real ( kind = 8 ) C(N), the coefficients of the interpolating
!    polynomial.  C(1) is the constant term, and C(N) multiplies X^(N-1).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  call vandermonde_interp_1d_matrix ( n, x, a )

  call qr_solve ( n, n, a, y, c )

  return
end
subroutine vandermonde_interp_1d_matrix ( n, x, a )

!*****************************************************************************80
!
!! VANDERMONDE_INTERP_1D_MATRIX computes a Vandermonde 1D interpolation matrix.
!
!  Discussion:
!
!    We assume the interpolant has the form
!
!      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
!
!    We have n data values (x(i),y(i)) which must be interpolated:
!
!      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
!
!    This can be cast as an NxN linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
!      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
!      [ ...................... ] [ ... ] = [ ... ]
!      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
!
!    and if the x values are distinct, the matrix A is theoretically
!    invertible (though in fact, generally badly conditioned).
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
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, real ( kind = 8 ) X(N), the data values.
!
!    Output, real ( kind = 8 ) A(N,N), the Vandermonde matrix for X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  a(1:n,1) = 1.0D+00
  do j = 2, n
    a(1:n,j) = a(1:n,j-1) * x(1:n)
  end do

  return
end
