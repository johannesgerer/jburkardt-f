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
function j_double_product_integral ( i, j, a, b )

!*****************************************************************************80
!
!! J_DOUBLE_PRODUCT_INTEGRAL: integral of J(i,x)*J(j,x)*(1-x)^a*(1+x)^b.
!
!  Discussion:
!
!    VALUE = integral ( -1 <= x <= +1 ) J(i,x)*J(j,x)*(1-x)^a*(1+x)^b dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the polynomial indices.
!
!    Input, real ( kind = 8 ) A, B, the parameters.
!    -1 < A, B.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  integer ( kind = 4 ) j
  real ( kind = 8 ) j_double_product_integral
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) value

  if ( i /= j ) then

    value = 0.0D+00

  else

    i_r8 = real ( i, kind = 8 )

    value = 2.0D+00 ** ( a + b + 1.0D+00 ) &
      / ( 2.0D+00 * i_r8 + a + b + 1.0D+00 ) &
      * gamma ( i_r8 + a + 1.0D+00 ) &
      * gamma ( i_r8 + b + 1.0D+00 ) &
      / r8_factorial ( i ) &
      / gamma ( i_r8 + a + b + 1.0D+00 )

  end if

  j_double_product_integral = value

  return
end
function j_integral ( n )

!*****************************************************************************80
!
!! J_INTEGRAL evaluates a monomial integral associated with J(n,a,b,x).
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
!    19 April 2012
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
!    Output, real ( kind = 8 ) J_INTEGRAL, the value of the integral.
!
  implicit none

  real ( kind = 8 ) j_integral
  integer ( kind = 4 ) n
  real ( kind = 8 ) value

  if ( mod ( n, 2 ) == 1 ) then
    value = 0.0D+00
  else
    value = 2.0D+00 / real ( n + 1, kind = 8 )
  end if

  j_integral = value

  return
end
subroutine j_polynomial ( m, n, alpha, beta, x, cx )

!*****************************************************************************80
!
!! J_POLYNOMIAL evaluates the Jacobi polynomials J(n,a,b,x).
!
!  Differential equation:
!
!    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
!
!  Recursion:
!
!    P(0,ALPHA,BETA,X) = 1,
!
!    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
!
!    P(N,ALPHA,BETA,X)  = 
!      ( 
!        (2*N+ALPHA+BETA-1) 
!        * ((ALPHA^2-BETA^2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X) 
!        * P(N-1,ALPHA,BETA,X)
!        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
!      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
!
!  Restrictions:
!
!    -1 < ALPHA
!    -1 < BETA
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA 
!      * P(N,ALPHA,BETA,X)^2 dX 
!    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
!      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
!
!  Special values:
!
!    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2012
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
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.  
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) ALPHA, one of the parameters defining the Jacobi
!    polynomials, ALPHA must be greater than -1.
!
!    Input, real ( kind = 8 ) BETA, the second parameter defining the Jacobi
!    polynomials, BETA must be greater than -1.
!
!    Input, real ( kind = 8 ) X(M), the point at which the polynomials are 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CX(M,0:N), the values of the first N+1 Jacobi
!    polynomials at the point X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) cx(1:m,0:n)
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) c4
  integer ( kind = 4 ) i
  real ( kind = 8 ) r_i
  real ( kind = 8 ) x(m)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'J_POLYNOMIAL - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of ALPHA = ', alpha
    write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
    stop
  end if
 
  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'J_POLYNOMIAL - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of BETA = ', beta
    write ( *, '(a)' ) '  But BETA must be greater than -1.'
    stop
  end if
  
  if ( n < 0 ) then
    return
  end if

  cx(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1:m,1) = ( 1.0D+00 + 0.5D+00 * ( alpha + beta ) ) * x(1:m) &
    + 0.5D+00 * ( alpha - beta )
 
  do i = 2, n

    r_i = real ( i, kind = 8 ) 

    c1 = 2.0D+00 * r_i * ( r_i + alpha + beta ) &
       * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

    c2 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) &
       * ( 2.0D+00 * r_i  + alpha + beta ) &
       * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

    c3 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) &
      * ( alpha + beta ) * ( alpha - beta )

    c4 = - 2.0D+00 * ( r_i - 1.0D+00 + alpha ) &
      * ( r_i - 1.0D+00 + beta )  * ( 2.0D+00 * r_i + alpha + beta )

    cx(1:m,i) = ( ( c3 + c2 * x(1:m) ) * cx(1:m,i-1) + c4 * cx(1:m,i-2) ) / c1

  end do

  return
end
subroutine j_polynomial_values ( n_data, n, a, b, x, fx )

!*****************************************************************************80
!
!! J_POLYNOMIAL_VALUES returns some values of the Jacobi polynomial.
!
!  Discussion:
!
!    In Mathematica, the function
!
!      JacobiP[ n, a, b, x ]
!
!    returns the value of the Jacobi polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2004
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
!    Output, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Output, real ( kind = 8 ) A, B, parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 26

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, &
     3.0D+00, 4.0D+00, 5.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00, 5.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
      1.000000000000000D+00, &
      0.250000000000000D+00, &
     -0.375000000000000D+00, &
     -0.484375000000000D+00, &
     -0.132812500000000D+00, &
      0.275390625000000D+00, &
     -0.164062500000000D+00, &
     -1.174804687500000D+00, &
     -2.361328125000000D+00, &
     -2.616210937500000D+00, &
      0.117187500000000D+00, &
      0.421875000000000D+00, &
      0.504882812500000D+00, &
      0.509765625000000D+00, &
      0.430664062500000D+00, &
     -6.000000000000000D+00, &
      0.038620000000000D+00, &
      0.811840000000000D+00, &
      0.036660000000000D+00, &
     -0.485120000000000D+00, &
     -0.312500000000000D+00, &
      0.189120000000000D+00, &
      0.402340000000000D+00, &
      0.012160000000000D+00, &
     -0.439620000000000D+00, &
      1.000000000000000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, 1, 2, 3, &
     4, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00, -1.0D+00, &
    -0.8D+00, -0.6D+00, -0.4D+00, -0.2D+00, &
     0.0D+00,  0.2D+00,  0.4D+00,  0.6D+00, &
     0.8D+00,  1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine j_polynomial_zeros ( n, alpha, beta, x )

!*****************************************************************************80
!
!! J_POLYNOMIAL_ZEROS: zeros of Jacobi polynomial J(n,a,b,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ), N, the order.
!
!    Input, real ( kind = 8 ), ALPHA, BETA, the parameters.
!    -1 < ALPHA, BETA.
!
!    Output, real ( kind = 8 ) X(N), the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a2b2
  real ( kind = 8 ) ab
  real ( kind = 8 ) abi
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu

  ab = alpha + beta
  abi = 2.0D+00 + ab
!
!  Define the zero-th moment.
!
  zemu = 2.0D+00 ** ( ab + 1.0D+00 ) * gamma ( alpha + 1.0D+00 ) &
    * gamma ( beta + 1.0D+00 ) / gamma ( abi )
!
!  Define the Jacobi matrix.
!
  x(1) = ( beta - alpha ) / abi
  x(2:n) = 0.0D+00

  bj(1) = 4.0D+00 * ( 1.0D+00 + alpha ) * ( 1.0D+00 + beta ) &
    / ( ( abi + 1.0D+00 ) * abi * abi )
  bj(2:n) = 0.0D+00

  a2b2 = beta * beta - alpha * alpha

  do i = 2, n
    i_r8 = real ( i, kind = 8 )
    abi = 2.0D+00 * i_r8 + ab
    x(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
    abi = abi ** 2
    bj(i) = 4.0D+00 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) &
      * ( i_r8 + ab ) / ( ( abi - 1.0D+00 ) * abi )
  end do
  bj(1:n) =  sqrt ( bj(1:n) )

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  return
end
subroutine j_quadrature_rule ( n, alpha, beta, x, w )

!*****************************************************************************80
!
!! J_QUADRATURE_RULE: Gauss-Jacobi quadrature based on J(n,a,b,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ), N, the order.
!
!    Input, real ( kind = 8 ), ALPHA, BETA, the parameters.
!    -1 < ALPHA, BETA.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a2b2
  real ( kind = 8 ) ab
  real ( kind = 8 ) abi
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu

  ab = alpha + beta
  abi = 2.0D+00 + ab
!
!  Define the zero-th moment.
!
  zemu = 2.0D+00 ** ( ab + 1.0D+00 ) * gamma ( alpha + 1.0D+00 ) &
    * gamma ( beta + 1.0D+00 ) / gamma ( abi )
!
!  Define the Jacobi matrix.
!
  x(1) = ( beta - alpha ) / abi
  x(2:n) = 0.0D+00

  bj(1) = 4.0D+00 * ( 1.0D+00 + alpha ) * ( 1.0D+00 + beta ) &
    / ( ( abi + 1.0D+00 ) * abi * abi )
  bj(2:n) = 0.0D+00

  a2b2 = beta * beta - alpha * alpha

  do i = 2, n
    i_r8 = real ( i, kind = 8 )
    abi = 2.0D+00 * i_r8 + ab
    x(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
    abi = abi ** 2
    bj(i) = 4.0D+00 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) &
      * ( i_r8 + ab ) / ( ( abi - 1.0D+00 ) * abi )
  end do
  bj(1:n) =  sqrt ( bj(1:n) )

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n) ** 2

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
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value = +1 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 8 ) R8_SIGN, the sign of X:
!
  implicit none

  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

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
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
