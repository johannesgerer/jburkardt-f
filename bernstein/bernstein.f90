subroutine bernstein_matrix ( n, a )

!*****************************************************************************80
!
!! BERNSTEIN_MATRIX returns the Bernstein matrix.
!
!  Discussion:
!
!    The Bernstein matrix of order N is an NxN matrix A which can be used to
!    transform a vector of power basis coefficients C representing a polynomial 
!    P(X) to a corresponding Bernstein basis coefficient vector B:
!
!      B = A * C
!
!    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
!    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the Bernstein matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) n0
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop

  a(1:n,1:n) = 0.0D+00

  n0 = n - 1

  do j0 = 0, n0
    do i0 = 0, j0
      a(i0+1,j0+1) = r8_mop ( j0 - i0 ) * r8_choose ( n0 - i0, j0 - i0 ) &
        * r8_choose ( n0, i0 )
    end do
  end do

  return
end
subroutine bernstein_matrix_inverse ( n, a )

!*****************************************************************************80
!
!! BERNSTEIN_MATRIX_INVERSE returns the inverse Bernstein matrix.
!
!  Discussion:
!
!    The inverse Bernstein matrix of order N is an NxN matrix A which can 
!    be used to transform a vector of Bernstein basis coefficients B
!    representing a polynomial P(X) to a corresponding power basis 
!    coefficient vector C:
!
!      C = A * B
!
!    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
!    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the inverse Bernstein matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) n0
  real ( kind = 8 ) r8_choose

  a(1:n,1:n) = 0.0D+00

  n0 = n - 1

  do j0 = 0, n0
    do i0 = 0, j0
      a(i0+1,j0+1) = r8_choose ( j0, i0 ) / r8_choose ( n0, i0 )
    end do
  end do

  return
end
subroutine bernstein_poly ( n, x, bern )

!*****************************************************************************80
!
!! BERNSTEIN_POLY evaluates the Bernstein polynomials at a point X.
!
!  Discussion:
!
!    The Bernstein polynomials are assumed to be based on [0,1].
!
!    The formula is:
!
!      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
!
!  First values:
!
!    B(0,0)(X) = 1
!
!    B(1,0)(X) =      1-X
!    B(1,1)(X) =                X
!
!    B(2,0)(X) =     (1-X)^2
!    B(2,1)(X) = 2 * (1-X)    * X
!    B(2,2)(X) =                X^2
!
!    B(3,0)(X) =     (1-X)^3
!    B(3,1)(X) = 3 * (1-X)^2 * X
!    B(3,2)(X) = 3 * (1-X)   * X^2
!    B(3,3)(X) =               X^3
!
!    B(4,0)(X) =     (1-X)^4
!    B(4,1)(X) = 4 * (1-X)^3 * X
!    B(4,2)(X) = 6 * (1-X)^2 * X^2
!    B(4,3)(X) = 4 * (1-X)   * X^3
!    B(4,4)(X) =               X^4
!
!  Special values:
!
!    B(N,I)(X) has a unique maximum value at X = I/N.
!
!    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
!
!    B(N,I)(1/2) = C(N,K) / 2^N
!
!    For a fixed X and N, the polynomials add up to 1:
!
!      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bernstein polynomials 
!    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
!    each of degree N, which form a basis for polynomials on [0,1].
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) BERN(0:N), the values of the N+1 
!    Bernstein polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bern(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  if ( n == 0 ) then
 
    bern(0) = 1.0D+00
 
  else if ( 0 < n ) then
 
    bern(0) = 1.0D+00 - x
    bern(1) = x
 
    do i = 2, n
      bern(i) = x * bern(i-1)
      do j = i - 1, 1, -1
        bern(j) =             x   * bern(j-1) &
                + ( 1.0D+00 - x ) * bern(j)
      end do
      bern(0) = ( 1.0D+00 - x ) * bern(0)
    end do
 
  end if
 
  return
end
subroutine bernstein_poly_values ( n_data, n, k, x, b )

!*****************************************************************************80
!
!! BERNSTEIN_POLY_VALUES returns some values of the Bernstein polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and 
!    N_DATA is set to the index of the test data.  On each subsequent 
!    call, N_DATA is incremented and that test data is returned.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Output, integer ( kind = 4 ) K, the index of the polynomial.
! 
!    Output, real ( kind = 8 ) X, the argument of the polynomial.
!
!    Output, real ( kind = 8 ) B, the value of the polynomial B(N,K)(X).
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 15

  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( nmax ) :: b_vec = (/ &
    1.0D+00, &
    0.75D+00,       0.25D+00, &
    0.5625D+00,     0.3750D+00,   0.0625D+00, &
    0.421875D+00,   0.421875D+00, 0.140625D+00,  0.015625D+00, &
    0.31640625D+00, 0.421875D+00, 0.2109375D+00, 0.046875D+00, 0.00390625D+00 /)
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save, dimension ( nmax ) :: k_vec = (/ &
    0, &
    0, 1, &
    0, 1, 2, &
    0, 1, 2, 3, &
    0, 1, 2, 3, 4 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
    0, &
    1, 1, &
    2, 2, 2, &
    3, 3, 3, 3, &
    4, 4, 4, 4, 4 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    k = 0
    x = 0.0D+00
    b = 0.0D+00
  else
    n = n_vec(n_data)
    k = k_vec(n_data)
    x = x_vec(n_data)
    b = b_vec(n_data)
  end if

  return
end
subroutine bpab ( n, a, b, x, bern )

!*****************************************************************************80
!
!! BPAB evaluates at X the Bernstein polynomials based in [A,B].
!
!  Discussion:
!
!    The formula is:
!
!      BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
!
!  First values:
!
!    B(0,0)(X) =   1
!
!    B(1,0)(X) = (      B-X                ) / (B-A)
!    B(1,1)(X) = (                 X-A     ) / (B-A)
!
!    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
!    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
!    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
!
!    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
!    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
!    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
!    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
!
!    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
!    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
!    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
!    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
!    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bernstein polynomials 
!    to be used.  For any N, there is a set of N+1 Bernstein polynomials, 
!    each of degree N, which form a basis for polynomials on [A,B].
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval on which the
!    polynomials are to be based.  A and B should not be equal.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials 
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) BERN(0:N), the values of the N+1
!    Bernstein polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bern(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  if ( b == a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BPAB - Fatal error!'
    write ( *, '(a,g14.6)' ) '  A = B = ', a
    stop
  end if

  if ( n == 0 ) then
 
    bern(0) = 1.0D+00
 
  else if ( 0 < n ) then
 
    bern(0) = ( b - x ) / ( b - a )
    bern(1) = ( x - a ) / ( b - a )
 
    do i = 2, n
      bern(i) = ( x - a ) * bern(i-1) / ( b - a )
      do j = i - 1, 1, -1
        bern(j) = ( ( b - x     ) * bern(j)     &
                  + (     x - a ) * bern(j-1) ) &
                  / ( b     - a )
      end do
      bern(0) = ( b - x ) * bern(0) / ( b - a )
    end do
 
  end if
 
  return
end
subroutine bpab_approx ( n, a, b, ydata, nval, xval, yval )

!*****************************************************************************80
!
!! BPAB_APPROX evaluates the Bernstein polynomial approximant to F(X) on [A,B].
!
!  Formula:
!
!    BPAB(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
!
!    where
!
!      X(I) = ( ( N - I ) * A + I * B ) / N
!      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
!
!  Discussion:
!
!    The Bernstein polynomial BPAB(F) for F(X) over [A,B] is an approximant, 
!    not an interpolant; in other words, its value is not guaranteed to equal
!    that of F at any particular point.  However, for a fixed interval
!    [A,B], if we let N increase, the Bernstein polynomial converges
!    uniformly to F everywhere in [A,B], provided only that F is continuous.
!    Even if F is not continuous, but is bounded, the polynomial converges
!    pointwise to F(X) at all points of continuity.  On the other hand,
!    the convergence is quite slow compared to other interpolation
!    and approximation schemes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2004
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
!    Input, integer ( kind = 4 ) N, the degree of the Bernstein polynomial
!    to be used.  N must be at least 0.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval on which the
!    approximant is based.  A and B should not be equal.
!
!    Input, real ( kind = 8 ) YDATA(0:N), the data values at N+1 equally
!    spaced points in [A,B].  If N = 0, then the evaluation point should
!    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
!    ( (N-I)*A + I*B ) / N ).
!
!    Input, integer ( kind = 4 ) NVAL, the number of points at which the
!    approximant is to be evaluated.
!
!    Input, real ( kind = 8 ) XVAL(NVAL), the point at which the Bernstein 
!    polynomial approximant is to be evaluated.  The entries of XVAL do not 
!    have to lie in the interval [A,B].
!
!    Output, real ( kind = 8 ) YVAL(NVAL), the values of the Bernstein 
!    polynomial approximant for F, based in [A,B], evaluated at XVAL.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nval

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bvec(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) xval(nval)
  real ( kind = 8 ) ydata(0:n)
  real ( kind = 8 ) yval(nval)

  do i = 1, nval
!
!  Evaluate the Bernstein basis polynomials at XVAL.
!
    call bpab ( n, a, b, xval(i), bvec )
!
!  Now compute the sum of YDATA(I) * BVEC(I).
!
    yval(i) = dot_product ( ydata(0:n), bvec(0:n) )

  end do

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, real ( kind = 8 ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = 8 )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = 8 ) ) / real ( i, kind = 8 )
    end do

  end if

  r8_choose = value

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, real ( kind = 8 ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_is_identity ( n, a, error_frobenius )

!*****************************************************************************80
!
!! R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!    The routine returns the Frobenius norm of A - I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A - I, which would be exactly zero
!    if A were the identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) value

  error_frobenius = 0.0D+00

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        error_frobenius = error_frobenius + ( a(i,j) - 1.0D+00 )**2
      else
        error_frobenius = error_frobenius + a(i,j)**2
      end if
    end do 
  end do

  error_frobenius = sqrt ( error_frobenius )

  return
end
function r8mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO returns the Frobenius norm of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) Sum ( 1 <= J <= N ) A(I,J)^2 )
!
!    The matrix Frobenius-norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A*x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_norm_fro

  r8mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer( kind = 4 ) m
  integer( kind = 4 ) n

  real   ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * a_first &
             + real (     i - 1, kind = 8 ) * a_last ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

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

  character ( len = 8 ) ampm
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
