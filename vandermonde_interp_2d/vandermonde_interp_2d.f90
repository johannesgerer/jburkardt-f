function triangle_num ( n )

!*****************************************************************************80
!
!! TRIANGLE_NUM returns the N-th triangular number.
!
!  Discussion:
!
!    The N-th triangular number T(N) is formed by the sum of the first
!    N integers:
!
!      T(N) = sum ( 1 <= I <= N ) I
!
!    By convention, T(0) = 0.
!
!    T(N) can be computed quickly by the formula:
!
!      T(N) = ( N * ( N + 1 ) ) / 2
!
!  First Values:
!
!     0
!     1
!     3
!     6
!    10
!    15
!    21
!    28
!    36
!    45
!    55
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired number, 
!    which must be at least 0.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NUM, the N-th triangular number.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num

  triangle_num = ( n * ( n + 1 ) ) / 2

  return
end
subroutine vandermonde_interp_2d_matrix ( n, m, x, y, a )

!*****************************************************************************80
!
!! VANDERMONDE_INTERP_2D_MATRIX computes a Vandermonde 2D interpolation matrix.
!
!  Discussion:
!
!    We assume the approximating function has the form of a polynomial
!    in X and Y of total degree M.
!
!      p(x,y) = c00 
!             + c10 * x                + c01 * y
!             + c20 * x^2   + c11 * xy + c02 * y^2
!             + ...
!             + cm0 * x^(m) + ...      + c0m * y^m.
!
!    If we let T(K) = the K-th triangular number 
!            = sum ( 1 <= I <= K ) I
!    then the number of coefficients in the above polynomial is T(M+1).
!
!    We have n data locations (x(i),y(i)) and values z(i) to approximate:
!
!      p(x(i),y(i)) = z(i)
!
!    and we assume that N = T(M+1).
!
!    This can be cast as an NxN linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
!      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
!      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
!      [ ...................... ] [ ... ] = [ ... ]
!      [ 1 xn yn  xn^2 ... yn^m ] [ c0n ] = [  zn ]
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  It is necessary 
!    that N = T(M+1), where T(K) is the K-th triangular number.
!
!    Input, integer ( kind = 4 ) M, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the data locations.
!
!    Output, real ( kind = 8 ) A(N,N), the Vandermonde matrix for X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) ex
  integer ( kind = 4 ) ey
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) s
  integer ( kind = 4 ) tmp1
  integer ( kind = 4 ) triangle_num
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  tmp1 = triangle_num ( m + 1 )

  if ( n /= tmp1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VANDERMONDE_INTERP_2D_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  For interpolation, we need N = T(M+1).'
    write ( *, '(a,i6)' ) '  But we have N = ', n
    write ( *, '(a,i6)' ) '  M = ', m
    write ( *, '(a,i6)' ) '  and T(M+1) = ', tmp1
    stop
  end if

  j = 0

  do s = 0, m
    do ex = s, 0, -1
      ey = s - ex
      j = j + 1
      a(1:n,j) = x(1:n) ** ex * y(1:n) ** ey
    end do
  end do
 
  return
end
