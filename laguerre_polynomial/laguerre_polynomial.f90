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
subroutine l_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! L_EXPONENTIAL_PRODUCT: exponential product table for L(n,x).
!
!  Discussion:
!
!    Let L(n,x) represent the Laguerre polynomial of degree n.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( 0 <= X < +oo ) exp(b*x) * L(i,x) * L(j,x) * exp (-x) dx
!
!    Because of the exponential factor, the quadrature will not be exact.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2012
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
!    TABLE(I,J) represents the weighted integral of exp(B*X) * L(i,x) * L(j,x).
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), allocatable :: l_table(:)
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = ( 3 * p + 4 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call l_quadrature_rule ( order, x_table, w_table )

  allocate ( l_table(0:p) )

  do k = 1, order

    x = x_table(k)
    call l_polynomial ( 1, p, x, l_table )

    do j = 0, p
      do i = 0, p
        table(i,j) = table(i,j) &
          + w_table(k) * exp ( b * x ) * l_table(i) * l_table(j)
      end do
    end do

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine l_integral ( n, exact )

!*****************************************************************************80
!
!! L_INTEGRAL evaluates a monomial integral associated with L(n,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial

  exact = r8_factorial ( n )

  return
end
subroutine l_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! L_POLYNOMIAL evaluates the Laguerre polynomial L(n,x).
!
!  First terms:
!
!      1
!     -X     +  1
!   (  X^2 -  4 X      +  2 ) / 2
!   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
!   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
!   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -   600 X    +  120 ) / 120
!   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 +  5400 X^2 -  4320 X     + 720 ) 
!     / 720
!   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 + 52920 X^2 - 35280 X 
!     + 5040 ) / 5040
!
!  Recursion:
!
!    L(0,X) = 1
!    L(1,X) = 1 - X
!    L(N,X) = (2*N-1-X)/N * L(N-1,X) - (N-1)/N * L(N-2,X)
!
!  Orthogonality:
!
!    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX = delta ( M, N )
!
!  Relations:
!
!    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
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
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  v(1:m,1) = 1.0D+00 - x(1:m)
 
  do j = 2, n

    v(1:m,j) = ( ( real ( 2 * j - 1, kind = 8 ) - x(1:m) ) * v(1:m,j-1)   &
                 - real (     j - 1, kind = 8 )            * v(1:m,j-2) ) &
                 / real (     j,     kind = 8 )

  end do

  return
end
subroutine l_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! L_POLYNOMIAL_COEFFICIENTS: coefficients of the Laguerre polynomial L(n,x).
!
!  First terms:
!
!    0: 1
!    1: 1  -1
!    2: 1  -2  1/2
!    3: 1  -3  3/2  1/6
!    4: 1  -4  4   -2/3  1/24
!    5: 1  -5  5   -5/3  5/24  -1/120
!
!  Recursion:
!
!    L(0) = ( 1,  0, 0, ..., 0 )
!    L(1) = ( 1, -1, 0, ..., 0 )
!    L(N) = (2*N-1-X) * L(N-1) - (N-1) * L(N-2) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2003
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the
!    Laguerre polynomials of degree 0 through N.  Each polynomial 
!    is stored as a row.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0:n,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = -1.0D+00
 
  do i = 2, n

    c(i,1:n) = ( &
        real ( 2 * i - 1, kind = 8 ) * c(i-1,1:n)     &
      + real (   - i + 1, kind = 8 ) * c(i-2,1:n)     &
      -                                c(i-1,0:n-1) ) &
      / real (     i,     kind = 8 )

  end do

  return
end
subroutine l_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! L_POLYNOMIAL_VALUES: some values of the Laguerre polynomial L(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,x]
!
!  Differential equation:
!
!    X * Y'' + (1-X) * Y' + N * Y = 0
!
!  First terms:
!
!      1
!     -X    +  1
!   (  X^2 -  4 X     +  2 ) / 2
!   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
!   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
!   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
!   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
!   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 + 52920 X^2 - 35280 X 
!     + 5040 ) / 5040
!
!  Recursion:
!
!    L(0,X) = 1,
!    L(1,X) = 1-X,
!    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
!
!  Orthogonality:
!
!    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX
!    = 0 if N /= M
!    = 1 if N == M
!
!  Special values:
!
!    L(N,0) = 1.
!
!  Relations:
!
!    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2004
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
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.0000000000000000D+00, &
    -0.5000000000000000D+00, &
    -0.6666666666666667D+00, &
    -0.6250000000000000D+00, &
    -0.4666666666666667D+00, &
    -0.2569444444444444D+00, &
    -0.4047619047619048D-01, &
     0.1539930555555556D+00, &
     0.3097442680776014D+00, &
     0.4189459325396825D+00, &
     0.4801341790925124D+00, &
     0.4962122235082305D+00, &
    -0.4455729166666667D+00, &
     0.8500000000000000D+00, &
    -0.3166666666666667D+01, &
     0.3433333333333333D+02 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    0.5D+00, &
    3.0D+00, &
    5.0D+00, &
    1.0D+01 /)

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
subroutine l_polynomial_zeros ( n, x )

!*****************************************************************************80
!
!! L_POLYNOMIAL_ZEROS: zeros of the Laguerre polynomial L(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X(N), the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = 1.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    x(i) = real ( 2 * i - 1, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  return
end
subroutine l_power_product ( p, e, table )

!*****************************************************************************80
!
!! L_POWER_PRODUCT: power product table for L(n,x).
!
!  Discussion:
!
!    Let L(n,x) represent the Laguerre polynomial of degree n.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X^E with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( 0 <= X < +oo ) x^e * L(i,x) * L(j,x) * exp (-x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2012
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
!    Input, integer ( kind = 4 ) E, the exponent of X.
!    0 <= E.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!    TABLE(I,J) represents the weighted integral of x^E * L(i,x) * L(j,x).
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), allocatable :: l_table(:)
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1 + ( e + 1 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call l_quadrature_rule ( order, x_table, w_table )

  allocate ( l_table(0:p) )

  do k = 1, order

    x = x_table(k)
    call l_polynomial ( 1, p, x, l_table )

    if ( e == 0 ) then

      do j = 0, p
        do i = 0, p
          table(i,j) = table(i,j) + w_table(k) * l_table(i) * l_table(j)
        end do
      end do

    else

      do j = 0, p
        do i = 0, p
          table(i,j) = table(i,j) &
            + w_table(k) * ( x ** e ) * l_table(i) * l_table(j)
        end do
      end do

    end if

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine l_quadrature_rule ( n, x, w )

!*****************************************************************************80
!
!! L_QUADRATURE_RULE: Gauss-Laguerre quadrature based on L(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = 1.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    x(i) = real ( 2 * i - 1, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine lf_integral ( n, alpha, exact )

!*****************************************************************************80
!
!! LF_INTEGRAL evaluates a monomial integral associated with Lf(n,alpha,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^n * x^alpha * exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) exact
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_gamma

  arg = alpha + real ( n + 1, kind = 8 )

  exact = r8_gamma ( arg )

  return
end
subroutine lf_function ( m, n, alpha, x, cx )

!*****************************************************************************80
!
!! LF_FUNCTION evaluates the Laguerre function Lf(n,alpha,x).
!
!  Recursion:
!
!    Lf(0,ALPHA,X) = 1
!    Lf(1,ALPHA,X) = 1+ALPHA-X
!
!    Lf(N,ALPHA,X) = (2*N-1+ALPHA-X)/N * Lf(N-1,ALPHA,X) 
!                      - (N-1+ALPHA)/N * Lf(N-2,ALPHA,X)
!
!  Restrictions:
!
!    -1 < ALPHA
!
!  Special values:
!
!    Lf(N,0,X) = L(N,X).
!    Lf(N,ALPHA,X) = LM(N,ALPHA,X) for ALPHA integral.
!
!  Norm:
!
!    Integral ( 0 <= X < +oo ) exp ( - X ) * Lf(N,ALPHA,X)^2 dX
!    = Gamma ( N + ALPHA + 1 ) / N!
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order function to compute.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.  -1 < ALPHA is required.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(1:M,0:N), the functions of 
!    degrees 0 through N evaluated at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cx(1:m,0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(1:m)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LF_FUNCTION - Fatal error!'
    write ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
    write ( *, '(a)' ) '  but ALPHA must be greater than -1.'
    stop
  end if
 
  if ( n < 0 ) then
    return
  end if

  cx(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1:m,1) = 1.0D+00 + alpha - x(1:m)

  do i = 2, n
    cx(1:m,i) = ( ( real ( 2 * i - 1, kind = 8 ) + alpha - x(1:m) ) * cx(1:m,i-1)   &
                + ( real (   - i + 1, kind = 8 ) - alpha          ) * cx(1:m,i-2) ) &
                  / real (     i,     kind = 8 )
  end do

  return
end
subroutine lf_function_values ( n_data, n, a, x, fx )

!*****************************************************************************80
!
!! LF_FUNCTION_VALUES returns values of the Laguerre function Lf(n,alpha,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,a,x]
!
!    The functions satisfy the following differential equation:
!
!      X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0
!
!    Function values can be generated by the recursion:
!
!      Lf(0,ALPHA,X) = 1
!      Lf(1,ALPHA,X) = 1+ALPHA-X
!
!      Lf(N,ALPHA,X) = ( (2*N-1+ALPHA-X) * Lf(N-1,ALPHA,X)
!                          - (N-1+ALPHA) * Lf(N-2,ALPHA,X) ) / N
!
!    The parameter ALPHA is required to be greater than -1.
!
!    For ALPHA = 0, the generalized Laguerre function Lf(N,ALPHA,X)
!    is equal to the Laguerre polynomial L(N,X).
!
!    For ALPHA integral, the generalized Laguerre function
!    Lf(N,ALPHA,X) equals the associated Laguerre polynomial Lm(N,ALPHA,X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
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
!    Output, real ( kind = 8 ) A, the parameter.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.00D+00, &
    0.25D+00, &
    0.50D+00, &
    0.75D+00, &
    1.50D+00, &
    2.50D+00, &
    5.00D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    1.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00, &
    5.20D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.3726399739583333D-01, &
     0.3494791666666667D+00, &
     0.8710042317708333D+00, &
     0.1672395833333333D+01, &
     0.6657625325520833D+01, &
     0.2395726725260417D+02, &
     0.2031344319661458D+03, &
     0.1284193996800000D+02, &
     0.5359924801587302D+01, &
     0.9204589064126984D+00, &
    -0.1341585114857143D+01, &
    -0.2119726307555556D+01, &
    -0.1959193658349206D+01, &
     0.1000000000000000D+01, &
     0.5450000000000000D+01, &
     0.1720125000000000D+02, &
     0.4110393750000000D+02, &
     0.8239745859375000D+02, &
     0.1460179186171875D+03, &
     0.2359204608298828D+03 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     8, &
     8, &
     8, &
     8, &
     8, &
     8, &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     6 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    1.00D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00, &
    0.75D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine lf_function_zeros ( n, alpha, x )

!*****************************************************************************80
!
!! LF_FUNCTION_ZEROS returns the zeros of Lf(n,alpha,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_gamma ( alpha + 1.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    bj(i) = i_r8 * ( i_r8 + alpha )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    x(i) = 2.0D+00 * i_r8 - 1.0D+00 + alpha
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  return
end
subroutine lf_quadrature_rule ( n, alpha, x, w )

!*****************************************************************************80
!
!! LF_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lf(n,alpha,x);
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( - x ) * x^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2011
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_gamma ( alpha + 1.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    bj(i) = i_r8 * ( i_r8 + alpha )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    x(i) = 2.0D+00 * i_r8 - 1.0D+00 + alpha
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine lm_integral ( n, m, exact )

!*****************************************************************************80
!
!! LM_INTEGRAL evaluates a monomial integral associated with Lm(n,m,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^n * x^m * exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
!    Input, integer ( kind = 4 ) M, the parameter.
!    0 <= M.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial

  exact = r8_factorial ( n + m )

  return
end
subroutine lm_polynomial ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
!
!  First terms:
!
!    M = 0
!
!    Lm(0,0,X) =   1
!    Lm(1,0,X) =  -X   +  1
!    Lm(2,0,X) =   X^2 -  4 X   +  2
!    Lm(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
!    Lm(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
!    Lm(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
!    Lm(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
!
!    M = 1
!
!    Lm(0,1,X) =    0
!    Lm(1,1,X) =   -1,
!    Lm(2,1,X) =    2 X - 4,
!    Lm(3,1,X) =   -3 X^2 + 18 X - 18,
!    Lm(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
!
!    M = 2
!
!    Lm(0,2,X) =    0
!    Lm(1,2,X) =    0,
!    Lm(2,2,X) =    2,
!    Lm(3,2,X) =   -6 X + 18,
!    Lm(4,2,X) =   12 X^2 - 96 X + 144
!
!    M = 3
!
!    Lm(0,3,X) =    0
!    Lm(1,3,X) =    0,
!    Lm(2,3,X) =    0,
!    Lm(3,3,X) =   -6,
!    Lm(4,3,X) =   24 X - 96
!
!    M = 4
!
!    Lm(0,4,X) =    0
!    Lm(1,4,X) =    0
!    Lm(2,4,X) =    0
!    Lm(3,4,X) =    0
!    Lm(4,4,X) =   24
!
!  Recursion:
!
!    Lm(0,M,X)   = 1 
!    Lm(1,M,X)   = (M+1-X)
!
!    if 2 <= N:
!
!      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X) 
!                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
!
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal 
!    to the Laguerre polynomials L(N,X).
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.  M must be nonnegative.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the associated Laguerre polynomials 
!    of degrees 0 through N evaluated at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(mm)

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LM_POLYNOMIAL - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

  cx(1:mm,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1:mm,1) = real ( m + 1, kind = 8 ) - x(1:mm)

  do i = 2, n
    cx(1:mm,i) = &
      ( ( real (   m + 2 * i - 1, kind = 8 ) - x(1:mm) ) * cx(1:mm,i-1)   &
        + real ( - m     - i + 1, kind = 8 )             * cx(1:mm,i-2) ) &
        / real (           i,     kind = 8 )
  end do

  return
end
subroutine lm_polynomial_coefficients ( n, m, c )

!*****************************************************************************80
!
!! LM_POLYNOMIAL_COEFFICIENTS: coefficients of Laguerre polynomial Lm(n,m,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the
!    Laguerre polynomials of degree 0 through N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,0) = real ( m + 1, kind = 8 )
  c(1,1) = -1.0D+00
 
  do i = 2, n

    c(i,0:i) = ( real (   m + 2 * i - 1, kind = 8 ) * c(i-1,0:i)   &
               + real ( - m     - i + 1, kind = 8 ) * c(i-2,0:i) ) &
               / real (           i,     kind = 8 )

    c(i,1:i) = c(i,1:i) - c(i-1,0:i) / real ( i, kind = 8 )

  end do

  return
end
subroutine lm_polynomial_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! LM_POLYNOMIAL_VALUES returns values of Laguerre polynomials Lm(n,m,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,m,x]
!
!    The associated Laguerre polynomials may be generalized so that the
!    parameter M is allowed to take on arbitrary noninteger values.
!    The resulting function is known as the generalized Laguerre function.
!
!    The polynomials satisfy the differential equation:
!
!      X * Y'' + (M+1-X) * Y' + (N-M) * Y = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
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
!    Output, integer ( kind = 4 ) M, the parameter.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.1625000000000000D+01, &
    0.1479166666666667D+01, &
    0.1148437500000000D+01, &
    0.4586666666666667D+00, &
    0.2878666666666667D+01, &
    0.8098666666666667D+01, &
    0.1711866666666667D+02, &
    0.1045328776041667D+02, &
    0.1329019368489583D+02, &
    0.5622453647189670D+02, &
    0.7484729341779436D+02, &
    0.3238912982762806D+03, &
    0.4426100000097533D+03, &
    0.1936876572288250D+04 /)
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
subroutine lm_polynomial_zeros ( n, m, x )

!*****************************************************************************80
!
!! LM_POLYNOMIAL_ZEROS returns the zeros for Lm(n,m,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) M, the parameter.
!    0 <= M.
!
!    Output, real ( kind = 8 ) X(N), the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_factorial ( m )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i * ( i + m ), kind = 8 )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    x(i) = real ( 2 * i - 1 + m, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  return
end
subroutine lm_quadrature_rule ( n, m, x, w )

!*****************************************************************************80
!
!! LM_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lm(n,m,x);
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( - x ) * x^m * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) M, the parameter.
!    0 <= M.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_factorial ( m )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i * ( i + m ), kind = 8 )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    x(i) = real ( 2 * i - 1 + m, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( n ) = product ( 1 <= i <= n ) i
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
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none
!
!  Coefficients for minimax approximation over (12, INF).
!
  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

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
