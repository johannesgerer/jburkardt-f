subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to 
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine. 
!
!    It has been modified to produce the product Q' * Z, where Z is an input 
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!    The changes consist (essentially) of applying the orthogonal 
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the 
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine p_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! P_EXPONENTIAL_PRODUCT: exponential products for P(n,x).
!
!  Discussion:
!
!    Let P(n,x) represent the Legendre polynomial of degree n.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -1.0 <= X <= +1.0 ) exp(B*X) * P(I,X) * P(J,X) dx
!
!    We will estimate these integrals using Gauss-Legendre quadrature.
!    Because of the exponential factor exp(B*X), the quadrature will not 
!    be exact.
!
!    However, when B = 0, the quadrature is exact.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, real ( kind = 8 ) B, the coefficient of X in the exponential factor.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = ( 3 * p + 4 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call p_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call p_polynomial ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    do j = 0, p
      do i = 0, p
        table(i,j) = table(i,j) &
          + w_table(k) * exp ( b * x ) * h_table(i) * h_table(j)
      end do
    end do

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine p_integral ( n, value )

!*****************************************************************************80
!
!! P_INTEGRAL evaluates a monomial integral associated with P(n,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x < +1 ) x^n dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the exponent.
!    0 <= N.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) exact
  integer ( kind = 4 ) n
  real ( kind = 8 ) value

  if ( mod ( n, 2 ) == 1 ) then
    value = 0.0D+00
  else
    value = 2.0D+00 / real ( n + 1, kind = 8 )
  end if

  return
end
subroutine p_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! P_POLYNOMIAL evaluates the Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(n,1) = 1.
!    P(n,-1) = (-1)^N.
!    | P(n,x) | <= 1 in [-1,1].
!
!    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!    The Legendre polynomials are orthogonal under the inner product defined
!    as integration from -1 to 1:
!
!      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
!        = 0 if I =/= J
!        = 2 / ( 2*I+1 ) if I = J.
!
!    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
!    where
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
!
!    The formula is:
!
!      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Differential equation:
!
!    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,x) =      1
!    P( 1,x) =      1 X
!    P( 2,x) = (    3 X^2 -       1)/2
!    P( 3,x) = (    5 X^3 -     3 X)/2
!    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
!
!  Recursion:
!
!    P(0,x) = 1
!    P(1,x) = x
!    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
!
!    P'(0,x) = 0
!    P'(1,x) = 1
!    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials 
!    of order 0 through N at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)

 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )
 
  end do
 
  return
end
subroutine p_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).
!
!  First terms:
!
!     1
!     0     1
!    -1/2   0      3/2
!     0    -3/2    0     5/2
!     3/8   0    -30/8   0     35/8
!     0    15/8    0   -70/8    0     63/8
!    -5/16  0    105/16  0   -315/16   0    231/16
!     0   -35/16   0   315/16   0   -693/16   0    429/16
!
!     1.00000
!     0.00000  1.00000
!    -0.50000  0.00000  1.50000
!     0.00000 -1.50000  0.00000  2.5000
!     0.37500  0.00000 -3.75000  0.00000  4.37500
!     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
!    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
!     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the 
!    Legendre polynomials of degree 0 through N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n <= 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0:i-2) =          real (   - i + 1, kind = 8 ) * c(i-2,0:i-2) &
                        / real (     i,     kind = 8 )
    c(i,1:i) = c(i,1:i) + real ( i + i - 1, kind = 8 ) * c(i-1,0:i-1) &
                        / real (     i,     kind = 8 )
  end do
 
  return
end
subroutine p_polynomial_prime ( m, n, x, vp )

!*****************************************************************************80
!
!! P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P'(n,x).
!
!  Discussion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) VP(M,0:N), the values of the derivatives of the
!    Legendre polynomials of order 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) vp(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00
  vp(1:m,0) = 0.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
  vp(1:m,1) = 1.0D+00
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )
 
    vp(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * ( v(1:m,i-1) + x(1:m) * vp(1:m,i-1) ) &
                - real (     i - 1, kind = 8 ) *   vp(1:m,i-2)               ) &
                / real (     i,     kind = 8 )
 
  end do
 
  return
end
subroutine p_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! P_POLYNOMIAL_VALUES returns values of the Legendre polynomials P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.2500000000000000D+00, &
    -0.4062500000000000D+00, &
    -0.3359375000000000D+00, &
     0.1577148437500000D+00, &
     0.3397216796875000D+00, &
     0.2427673339843750D-01, &
    -0.2799186706542969D+00, &
    -0.1524540185928345D+00, &
     0.1768244206905365D+00, &
     0.2212002165615559D+00, &
     0.0000000000000000D+00, &
    -0.1475000000000000D+00, &
    -0.2800000000000000D+00, &
    -0.3825000000000000D+00, &
    -0.4400000000000000D+00, &
    -0.4375000000000000D+00, &
    -0.3600000000000000D+00, &
    -0.1925000000000000D+00, &
     0.8000000000000000D-01, &
     0.4725000000000000D+00, &
     0.1000000000000000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine p_polynomial_zeros ( nt, t )

!*****************************************************************************80
!
!! P_POLYNOMIAL_ZEROS: zeros of Legendre function P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), the zeros.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  
  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = real ( i * i, kind = 8 ) / real ( 4 * i * i - 1, kind = 8 )
  end do
  bj(1:nt) = sqrt ( bj(1:nt) )

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( 2.0D+00 )

  call imtqlx ( nt, t, bj, wts )

  return
end
subroutine p_power_product ( p, e, table )

!*****************************************************************************80
!
!! P_POWER_PRODUCT: power products for Legendre polynomial P(n,x).
!
!  Discussion:
!
!    Let P(n,x) represent the Legendre polynomial of degree n.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -1.0 <= X <= +1.0 ) X^E * P(i,x) * P(j,x) dx
!
!    We will estimate these integrals using Gauss-Legendre quadrature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, integer ( kind = 4 ) E, the exponent of X in the integrand.
!    0 <= E.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  integer ( kind = 4 ) e
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1 + ( ( e + 1 ) / 2 )

  allocate ( x_table(order) )
  allocate ( w_table(order) )

  call p_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call p_polynomial ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    if ( e == 0 ) then
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) + w_table(k) * h_table(i) * h_table(j)
        end do
      end do
    else
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) &
            + w_table(k) * x ** e * h_table(i) * h_table(j)
        end do
      end do
    end if

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine p_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! P_QUADRATURE_RULE: quadrature for Legendre function P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), WTS(NT), the points and weights
!    of the rule.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  
  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = real ( i * i, kind = 8 ) / real ( 4 * i * i - 1, kind = 8 )
  end do
  bj(1:nt) = sqrt ( bj(1:nt) )

  wts(1) = sqrt ( 2.0D+00 )
  wts(2:nt) = 0.0D+00

  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt) ** 2

  return
end
subroutine pm_polynomial ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! PM_POLYNOMIAL evaluates the Legendre polynomials Pm(n,m,x).
!
!  Differential equation:
!
!    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N,X) )
!
!    Pm(0,0,x) =    1
!    Pm(1,0,x) =    1 X
!    Pm(2,0,x) = (  3 X^2 -   1)/2
!    Pm(3,0,x) = (  5 X^3 -   3 X)/2
!    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
!    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
!    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
!    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
!
!    M = 1
!
!    Pm(0,1,x) =   0
!    Pm(1,1,x) =   1 * SQRT(1-X^2)
!    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
!    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
!    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
!
!    M = 2
!
!    Pm(0,2,x) =   0
!    Pm(1,2,x) =   0
!    Pm(2,2,x) =   3 * (1-X^2)
!    Pm(3,2,x) =  15 * (1-X^2) * X
!    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)
!
!    M = 3
!
!    Pm(0,3,x) =   0
!    Pm(1,3,x) =   0
!    Pm(2,3,x) =   0
!    Pm(3,3,x) =  15 * (1-X^2)^1.5
!    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X
!
!    M = 4
!
!    Pm(0,4,x) =   0
!    Pm(1,4,x) =   0
!    Pm(2,4,x) =   0
!    Pm(3,4,x) =   0
!    Pm(4,4,x) = 105 * (1-X^2)^2
!
!  Recursion:
!
!    if N < M:
!      Pm(N,M,x) = 0
!    if N = M:
!      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
!    if M+1 < N:
!      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X(MM), the point at which the function is to be
!    evaluated.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(mm)

  cx(1:mm,0:n) = 0.0D+00
!
!  J = M is the first nonzero function.
!
  if ( m <= n ) then
    cx(1:mm,m) = 1.0D+00

    fact = 1.0D+00
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * fact * sqrt ( 1.0D+00 - x(1:mm)**2 )
      fact = fact + 2.0D+00
    end do

  end if
!
!  J = M + 1 is the second nonzero function.
!
  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
  end if
!
!  Now we use a three term recurrence.
!
  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &  
                 / real (     j - m,     kind = 8 )
  end do

  return
end
subroutine pm_polynomial_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! PM_POLYNOMIAL_VALUES returns values of Legendre polynomials Pm(n,m,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
    -0.5000000000000000D+00, &
     0.0000000000000000D+00, &
     0.3750000000000000D+00, &
     0.0000000000000000D+00, &
    -0.8660254037844386D+00, &
    -0.1299038105676658D+01, &
    -0.3247595264191645D+00, &
     0.1353164693413185D+01, &
    -0.2800000000000000D+00, &
     0.1175755076535925D+01, &
     0.2880000000000000D+01, &
    -0.1410906091843111D+02, &
    -0.3955078125000000D+01, &
    -0.9997558593750000D+01, &
     0.8265311444100484D+02, &
     0.2024442836815152D+02, &
    -0.4237997531890869D+03, &
     0.1638320624828339D+04, &
    -0.2025687389227225D+05 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 0, 0, &
    0, 1, 1, 1, &
    1, 0, 1, 2, &
    3, 2, 2, 3, &
    3, 4, 4, 5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    1,  2,  3,  4, &
    5,  1,  2,  3, &
    4,  3,  3,  3, &
    3,  4,  5,  6, &
    7,  8,  9, 10 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine pmn_polynomial ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! PMN_POLYNOMIAL: normalized Legendre polynomial Pmn(n,m,x).
!
!  Discussion:
!
!    The unnormalized associated Legendre functions P_N^M(X) have
!    the property that
!
!      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX 
!      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
!
!    By dividing the function by the square root of this term,
!    the normalized associated Legendre functions have norm 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) x(mm)

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMN_POLYNOMIAL - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop
  end if
 
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMN_POLYNOMIAL - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    stop
  end if

  cx(1:mm,0:n) = 0.0D+00

  if ( m <= n ) then
    cx(1:mm,m) = 1.0D+00
    factor = 1.0D+00
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * factor * sqrt ( 1.0D+00 - x(1:mm)**2 )
      factor = factor + 2.0D+00
    end do
  end if

  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
  end if

  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &  
                 / real (     j - m,     kind = 8 )
  end do
!
!  Normalization.
!
  do j = m, n
    factor = sqrt ( ( real ( 2 * j + 1, kind = 8 ) * r8_factorial ( j - m ) ) &
      / ( 2.0D+00 * r8_factorial ( j + m ) ) )
    cx(1:mm,j) = cx(1:mm,j) * factor
  end do

  return
end
subroutine pmn_polynomial_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! PMN_POLYNOMIAL_VALUES: normalized Legendre polynomial Pmn(n,m,x).
!
!  Discussion:
!
!    In Mathematica, the unnormalized function can be evaluated by:
!
!      LegendreP [ n, m, x ]
!
!    The function is normalized by dividing by the factor:
!
!      sqrt ( 2 * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7071067811865475D+00, &
    0.6123724356957945D+00, &
   -0.7500000000000000D+00, &
   -0.1976423537605237D+00, &
   -0.8385254915624211D+00, &
    0.7261843774138907D+00, &
   -0.8184875533567997D+00, &
   -0.1753901900050285D+00, &
    0.9606516343087123D+00, &
   -0.6792832849776299D+00, &
   -0.6131941618102092D+00, &
    0.6418623720763665D+00, &
    0.4716705890038619D+00, &
   -0.1018924927466445D+01, &
    0.6239615396237876D+00, &
    0.2107022704608181D+00, &
    0.8256314721961969D+00, &
   -0.3982651281554632D+00, &
   -0.7040399320721435D+00, &
    0.1034723155272289D+01, &
   -0.5667412129155530D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 1, 0, &
    1, 2, 0, 1, &
    2, 3, 0, 1, &
    2, 3, 4, 0, &
    1, 2, 3, 4, &
    5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  1,  2, &
    2,  2,  3,  3, &
    3,  3,  4,  4, &
    4,  4,  4,  5, &
    5,  5,  5,  5, &
    5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine pmns_polynomial ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! PMNS_POLYNOMIAL: sphere-normalized Legendre polynomial Pmns(n,m,x).
!
!  Discussion:
!
!    The unnormalized associated Legendre functions P_N^M(X) have
!    the property that
!
!      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX 
!      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
!
!    By dividing the function by the square root of this term,
!    the normalized associated Legendre functions have norm 1.
!
!    However, we plan to use these functions to build spherical
!    harmonics, so we use a slightly different normalization factor of
!
!      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) x(mm)

  cx(1:m,0:n) = 0.0D+00

  if ( m <= n ) then
    cx(1:mm,m) = 1.0D+00
    factor = 1.0D+00
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * factor * sqrt ( 1.0D+00 - x(1:mm)**2 )
      factor = factor + 2.0D+00
    end do
  end if

  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
  end if

  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &  
                 / real (     j - m,     kind = 8 )
  end do
!
!  Normalization.
!
  do j = m, n
    factor = sqrt ( ( real ( 2 * j + 1, kind = 8 ) * r8_factorial ( j - m ) ) &
      / ( 4.0D+00 * pi * r8_factorial ( j + m ) ) )
    cx(1:mm,j) = cx(1:mm,j) * factor
  end do

  return
end
subroutine pmns_polynomial_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! PMNS_POLYNOMIAL_VALUES: sphere-normalized Legendre polynomial Pmns(n,m,x).
!
!  Discussion:
!
!    In Mathematica, the unnormalized function can be evaluated by:
!
!      LegendreP [ n, m, x ]
!
!    The function is normalized by dividing by the factor:
!
!      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.2820947917738781D+00, &
     0.2443012559514600D+00, &
    -0.2992067103010745D+00, &
    -0.07884789131313000D+00, &
    -0.3345232717786446D+00, &
     0.2897056515173922D+00, &
    -0.3265292910163510D+00, &
    -0.06997056236064664D+00, &
     0.3832445536624809D+00, &
    -0.2709948227475519D+00, &
    -0.2446290772414100D+00, &
     0.2560660384200185D+00, &
     0.1881693403754876D+00, &
    -0.4064922341213279D+00, &
     0.2489246395003027D+00, &
     0.08405804426339821D+00, &
     0.3293793022891428D+00, &
    -0.1588847984307093D+00, &
    -0.2808712959945307D+00, &
     0.4127948151484925D+00, &
    -0.2260970318780046D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 1, 0, &
    1, 2, 0, 1, &
    2, 3, 0, 1, &
    2, 3, 4, 0, &
    1, 2, 3, 4, &
    5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  1,  2, &
    2,  2,  3,  3, &
    3,  3,  4,  4, &
    4,  4,  4,  5, &
    5,  5,  5,  5, &
    5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine pn_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! PN_POLYNOMIAL evaluates the normalized Legendre polynomials Pn(n,x).
!
!  Discussion:
!
!    The normalized Legendre polynomials are orthonormal under the inner product 
!    defined as integration from -1 to 1:
!
!      Integral ( -1 <= x <= +1 ) Pn(i,x) * Pn(j,x) dx = delta(i,j)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials 
!    of order 0 through N at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) norm
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  call p_polynomial ( m, n, x, v )

  do j = 0, n
    norm = sqrt ( 2.0D+00 / real ( 2 * j + 1, kind = 8 ) )
    v(1:m,j) = v(1:m,j) / norm
  end do
 
  return
end
subroutine pn_pair_product ( p, table )

!*****************************************************************************80
!
!! PN_PAIR_PRODUCT: pair products for normalized Legendre polynomial Pn(n,x).
!
!  Discussion:
!
!    Let Pn(n,x) represent the normalized Legendre polynomial of degree n.  
!
!    To check orthonormality, we compute
!
!      Tij = Integral ( -1.0 <= X <= +1.0 ) Pn(i,x) * Pn(j,x) dx
!
!    We will estimate these integrals using Gauss-Legendre quadrature.
!
!    The computed table should be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1

  allocate ( x_table(order) )
  allocate ( w_table(order) )

  call p_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call pn_polynomial ( 1, p, x, h_table )

    do i = 0, p
      do j = 0, p
        table(i,j) = table(i,j) + w_table(k) * h_table(i) * h_table(j)
      end do
    end do

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
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

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

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
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
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
