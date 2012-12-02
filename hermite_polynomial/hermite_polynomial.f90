subroutine h_integral ( n, value )

!*****************************************************************************80
!
!! H_INTEGRAL evaluates a monomial physicist's Hermite integral for H(n,x).
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
!    Input, integer ( kind = 4 ) N, the order of the integral.  
!    0 <= N.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_factorial2
  real ( kind = 8 ) value

  if ( mod ( n, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = r8_factorial2 ( n - 1 ) * sqrt ( pi ) / 2.0D+00**( n / 2 )

  end if

  return
end
subroutine h_polynomial ( m, n, x, p )

!*****************************************************************************80
!
!! H_POLYNOMIAL evaluates the physicist's Hermite polynomial H(n,x).
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!    H(N,X) = (-1)^N * exp ( X^2 ) * dn/dXn ( exp(-X^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2002
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
!    Larry Andrews,
!    Special Functions of Mathematics for Engineers,
!    Second Edition, 
!    Oxford University Press, 1998.
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
!    Output, real ( kind = 8 ) P(M,0:N), the values of the first N+1 Hermite
!    polynomials at the point X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) p(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = 2.0D+00 * x(1:m)
 
  do j = 2, n
    p(1:m,j) = 2.0D+00 * x(1:m) * p(1:m,j-1) &
             - 2.0D+00 * real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
 
  return
end
subroutine h_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! H_POLYNOMIAL_COEFFICIENTS: coeffs of physicist's Hermite polynomial H(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     2
!     2     -2     0      4
!     3      0   -12      0      8
!     4     12     0    -48      0      16
!     5      0   120      0   -160       0    32
!     6   -120     0    720      0    -480     0     64
!     7      0 -1680      0   3360       0 -1344      0   128
!     8   1680     0 -13440      0   13440     0  -3584     0    256
!     9      0 30240      0 -80640       0 48384      0 -9216      0 512
!    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024 
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
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the polynomials
!    of orders 0 through N.
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

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 2.0D+00
 
  do i = 2, n
    c(i,0)     =  -2.0D+00 * real ( i - 1, kind = 8 ) * c(i-2,0)
    c(i,1:i-2) =   2.0D+00                            * c(i-1,0:i-3)  &
                  -2.0D+00 * real ( i - 1, kind = 8 ) * c(i-2,1:i-2)
    c(i,  i-1) =   2.0D+00                            * c(i-1,  i-2)
    c(i,  i  ) =   2.0D+00                            * c(i-1,  i-1)
  end do
 
  return
end
subroutine h_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! H_POLYNOMIAL_VALUES: values of the physicist's Hermite polynomial H(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      HermiteH[n,x]
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
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

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
      0.1000000000000000D+01, &
      0.1000000000000000D+02, &
      0.9800000000000000D+02, &
      0.9400000000000000D+03, &
      0.8812000000000000D+04, &
      0.8060000000000000D+05, &
      0.7178800000000000D+06, &
      0.6211600000000000D+07, &
      0.5206568000000000D+08, &
      0.4212712000000000D+09, &
      0.3275529760000000D+10, &
      0.2432987360000000D+11, &
      0.1712370812800000D+12, &
      0.0000000000000000D+00, &
      0.4100000000000000D+02, &
     -0.8000000000000000D+01, &
      0.3816000000000000D+04, &
      0.3041200000000000D+07 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5, 5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    3.0D+00, &
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
subroutine h_polynomial_zeros ( nt, z )

!*****************************************************************************80
!
!! H_POLYNOMIAL_ZEROS: zeros of the physicist's Hermite polynomial H(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the degree of the polynomial.
!
!    Output, real ( kind = 8 ) Z(NT), the zeros of the polynomial.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) z(nt)

  z(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( pi ) )

  call imtqlx ( nt, z, bj, wts )

  return
end
subroutine h_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! H_QUADRATURE_RULE: quadrature for physicist's Hermite polynomial H(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
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
  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( pi ) )

  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt)**2

  return
end
function he_double_product_integral ( i, j )

!*****************************************************************************80
!
!! HE_DOUBLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*e^(-x^2/2).
!
!  Discussion:
!
!    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
!    Princeton, 2010,
!    ISBN13: 978-0-691-14212-8,
!    LC: QA274.23.X58.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the polynomial indices.
!
!    Output, real ( kind = 8 ) HE_DOUBLE_PRODUCT_INTEGRAL, the value of 
!    the integral.
!
  implicit none

  real ( kind = 8 ) he_double_product_integral
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) value

  if ( i == j ) then
    value = r8_factorial ( i )
  else
    value = 0.0D+00
  end if

  he_double_product_integral = value

  return
end
subroutine he_integral ( n, value )

!*****************************************************************************80
!
!! HE_INTEGRAL evaluates a monomial probabilist's Hermite integral for He(n,x).
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
!    Input, integer ( kind = 4 ) N, the order of the integral.  
!    0 <= N.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_factorial2
  real ( kind = 8 ) value

  if ( mod ( n, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = r8_factorial2 ( n - 1 ) * sqrt ( 2.0D+00 * pi )

  end if

  return
end
subroutine he_polynomial ( m, n, x, p )

!*****************************************************************************80
!
!! HE_POLYNOMIAL evaluates the probabilist's Hermite polynomial He(n,x).
!
!  Differential equation:
!
!    ( exp ( - 0.5 * x^2 ) * He(n,x)' )' + n * exp ( - 0.5 * x^2 ) * He(n,x) = 0
!
!  First terms:
!
!   1
!   X
!   X^2  -  1
!   X^3  -  3 X
!   X^4  -  6 X^2 +   3
!   X^5  - 10 X^3 +  15 X
!   X^6  - 15 X^4 +  45 X^2 -   15
!   X^7  - 21 X^5 + 105 X^3 -  105 X
!   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
!   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
!   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
!
!  Recursion:
!
!    He(0,X) = 1,
!    He(1,X) = X,
!    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
!
!  Orthogonality:
!
!    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
!    = sqrt ( 2 * pi ) * N! * delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2012
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
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real X(M), the evaluation points.
!
!    Output, real P(M,0:N), the values of the probabilist's Hermite polynomials 
!    of index 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) p(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = x(1:m)
 
  do j = 2, n
    p(1:m,j) = x(1:m) * p(1:m,j-1) - real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
 
  return
end
subroutine he_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! HE_POLYNOMIAL_COEFFICIENTS: coeffs probabilist's Hermite polynomial He(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     1
!     2     -1     0      1
!     3      0    -3      0      1
!     4      3     0     -6      0       1
!     5      0    15      0    -10       0     1
!     6    -15     0     45      0     -15     0      1
!     7      0  -105      0    105       0   -21      0     1
!     8    105     0   -420      0     210     0    -28     0      1
!     9      0   945      0  -1260       0   378      0   -36      0   1
!    10   -945     0   4725      0   -3150     0    630     0    -45   0    1
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
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the polynomials
!    of orders 0 through N.
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

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0)     =              - real ( i - 1, kind = 8 ) * c(i-2,0)
    c(i,1:i-2) = c(i-1,0:i-3) - real ( i - 1, kind = 8 ) * c(i-2,1:i-2)
    c(i,  i-1) = c(i-1,  i-2)
    c(i,  i  ) = c(i-1,  i-1)
  end do
 
  return
end
subroutine he_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! HE_POLYNOMIAL_VALUES: values of the probabilist's Hermite polynomial He(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      He(n,x) = HermiteH[n,x/Sqrt[2]] / Sqrt [ 2^n ] 
!
!  First terms:
!
!   1
!   X
!   X^2  -  1
!   X^3  -  3 X
!   X^4  -  6 X^2 +   3
!   X^5  - 10 X^3 +  15 X
!   X^6  - 15 X^4 +  45 X^2 -   15
!   X^7  - 21 X^5 + 105 X^3 -  105 X
!   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
!   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
!   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
!
!  Recursion:
!
!    He(0,X) = 1,
!    He(1,X) = X,
!    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
!    = sqrt ( 2 * pi ) * N! * delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
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

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.000000000000000D+00, &
    5.000000000000000D+00, &
    24.00000000000000D+00, &
    110.0000000000000D+00, &
    478.0000000000000D+00, &
    1950.000000000000D+00, &
    7360.000000000000D+00, &
    25100.00000000000D+00, &
    73980.00000000000D+00, &
    169100.0000000000D+00, &
    179680.0000000000D+00, &
   -792600.0000000000D+00, &
   -5939480.000000000D+00, &
    0.000000000000000D+00, &
    6.281250000000000D+00, &
    6.000000000000000D+00, &
    18.00000000000000D+00, &
    90150.00000000000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    3.0D+00, &
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
subroutine he_polynomial_zeros ( nt, z )

!*****************************************************************************80
!
!! HE_POLYNOMIAL_ZEROS: zeros of the probabilist's Hermite polynomial He(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the degree of the polynomial.
!
!    Output, real ( kind = 8 ) Z(NT), the zeros of the polynomial.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) z(nt)

  z(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( pi ) )

  call imtqlx ( nt, z, bj, wts )

  z(1:nt) = z(1:nt) * sqrt ( 2.0D+00 )

  return
end
subroutine he_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! HE_QUADRATURE_RULE: quadrature for probabilist's Hermite polynomial He(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
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
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( pi ) )

  call imtqlx ( nt, t, bj, wts )

  t(1:nt) = t(1:nt) * sqrt ( 2.0D+00 )
  wts(1:nt) = wts(1:nt)**2 * sqrt ( 2.0D+00 )

  return
end
function he_triple_product_integral ( i, j, k )

!*****************************************************************************80
!
!! HE_TRIPLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
!
!  Discussion:
!
!    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
!    Princeton, 2010,
!    ISBN13: 978-0-691-14212-8,
!    LC: QA274.23.X58.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, K, the polynomial indices.
!
!    Output, real ( kind = 8 ) HE_TRIPLE_PRODUCT_INTEGRAL, the value of the integral.
!
  implicit none

  real ( kind = 8 ) he_triple_product_integral
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) s
  real ( kind = 8 ) value

  s = ( i + j + k ) / 2

  if ( s < max ( i, j, k ) ) then
    value = 0.0D+00
  else if ( mod ( i + j + k, 2 ) /= 0 ) then
    value = 0.0D+00
  else
    value = r8_factorial ( i ) / r8_factorial ( s - i ) &
          * r8_factorial ( j ) / r8_factorial ( s - j ) &
          * r8_factorial ( k ) / r8_factorial ( s - k )
  end if

  he_triple_product_integral = value

  return
end
subroutine hen_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! HEN_EXPONENTIAL_PRODUCT: probabilist Hermite exponential product, Hen(n,x).
!
!  Discussion:
!
!    Let Hen(I,X) represent the normalized probabilist's Hermite polynomial 
!    of degree I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) 
!        exp(B*X) * Hen(I,X) * Hen(J,X) exp(-0.5*X*X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!    Because of the exponential factor exp(B*X), the quadrature will not 
!    be exact.
!
!    However, when B = 0, the quadrature is exact, and moreoever, the
!    table will be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the 
!    polyonomial factors.  0 <= P.
!
!    Input, real ( kind = 8 ) B, the coefficient of X in the exponential factor.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!    TABLE(I,J) represents the weighted integral of 
!    exp(B*X) * Hen(I,X) * Hen(J,X).
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

  call he_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k);
    call hen_polynomial ( 1, p, x, h_table )
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
subroutine hen_polynomial ( m, n, x, p )

!*****************************************************************************80
!
!! HEN_POLYNOMIAL: evaluate normalized probabilist's Hermite poly Hen(n,x).
!
!  Discussion:
!
!    These polynomials satisfy the orthonormality condition:
!
!      Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * Hen(M,X) Hen(N,X) dX 
!      = delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
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
!    Output, real ( kind = 8 ) P(M,0:N), the values of the polynomials of 
!    index 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fact
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(m,0:n)
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) x(m)

  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = x(1:m)
 
  do j = 2, n
    p(1:m,j) = x(1:m) * p(1:m,j-1) - real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
!
!  Normalize.
!
  fact = 1.0D+00
  do j = 0, n
    p(1:m,j) = p(1:m,j) / sqrt ( fact * sqrt ( 2.0D+00 * pi ) )
    fact = fact * real ( j + 1, kind = 8 )
  end do

  return
end
subroutine hen_power_product ( p, e, table )

!*****************************************************************************80
!
!! HEN_POWER_PRODUCT: power products, normalized probabilist's Hermite Hen(n,x).
!
!  Discussion:
!
!    Let Hen(I,X) represent the normalized probabilist's Hermite  polynomial 
!    of degree I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) 
!        X^E * Hen(I,X) * Hen(J,X) exp(-0.5*X*X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!
!    When E is 0, the computed table should be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
!    TABLE(I,J) represents the weighted integral of 
!    X^E * Hen(I,X) * Hen(J,X).
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

  call he_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hen_polynomial ( 1, p, x, h_table )
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
subroutine hf_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! HF_EXPONENTIAL_PRODUCT: exponential products, Hermite function Hf(n,x).
!
!  Discussion:
!
!    Let Hf(I,X) represent the Hermite function of "degree" I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) exp(B*X) * Hf(I,X) * Hf(J,X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!    Because of the exponential factor exp(B*X), the quadrature will not 
!    be exact.
!
!    However, when B = 0, the quadrature is exact, and moreoever, the
!    table will be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
!    TABLE(I,J) represents the integral of exp(B*X) * Hf(I,X) * Hf(J,X).
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

  call hf_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hf_function ( 1, p, x, h_table )
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
subroutine hf_function ( m, n, x, f )

!*****************************************************************************80
!
!! HF_FUNCTION evaluates the Hermite function Hf(n,x).
!
!  Discussion:
!
!    The Hermite function of degree n is related to the physicist's
!    Hermite polynomial H(n,x):
!
!      Hf(n,x) = H(n,x) * exp ( - 0.5 * x^2 ) / sqrt ( 2^n n! sqrt ( pi ) )
!
!    The Hermite functions are orthonormal:
!
!      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2012
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
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
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
!    Output, real ( kind = 8 ) F(M,0:N), the values of the Hermite functions 
!    of index 0 through N at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(m,0:n)
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) x(m)

  f(1:m,0) = exp ( - 0.5D+00 * x(1:m)**2 ) / sqrt ( sqrt ( pi ) )

  if ( n == 0 ) then
    return
  end if

  f(1:m,1) = 2.0D+00 * exp ( - 0.5D+00 * x(1:m)**2 ) * x(1:m) &
    / sqrt ( 2.0D+00 * sqrt ( pi ) )

  do j = 2, n
    f(1:m,j) = ( sqrt ( 2.0D+00 ) * x(1:m) * f(1:m,j-1) &
      - sqrt ( real ( j - 1, kind = 8 ) ) * f(1:m,j-2) ) &
      / sqrt ( real ( j, kind = 8 ) )
  end do
 
  return
end
subroutine hf_function_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! HF_FUNCTION_VALUES: values of the Hermite function Hf(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Hf(n,x) = HermiteH[n,x] 
!        * Exp [ -1/2 * x^2] / Sqrt [ 2^n * n! * Sqrt[Pi] ]
!
!    The Hermite functions are orthonormal:
!
!      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
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

  integer ( kind = 4 ), parameter :: n_max = 23

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7511255444649425D+00,  0.0000000000000000D+00, -0.5311259660135985D+00, &
    0.0000000000000000D+00,  0.4599685791773266D+00,  0.0000000000000000D+00, &
    0.4555806720113325D+00,  0.6442883651134752D+00,  0.3221441825567376D+00, &
   -0.2630296236233334D+00, -0.4649750762925110D+00, -0.5881521185179581D-01, &
    0.3905052515434106D+00,  0.2631861423064045D+00, -0.2336911435996523D+00, &
   -0.3582973361472840D+00,  0.6146344487883041D-01,  0.3678312067984882D+00, &
    0.9131969309166278D-01,  0.4385750950032321D+00, -0.2624689527931006D-01, &
    0.5138426125477819D+00,  0.9355563118061758D-01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  2,  &
    3,  4,  5,  &
    0,  1,  2,  &
    3,  4,  5,  &
    6,  7,  8,  &
    9, 10, 11,  &
   12,  5,  5,  &
    5,  5  /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 0.5D+00, 2.0D+00, &
    3.0D+00, 4.0D+00 /)

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
subroutine hf_power_product ( p, e, table )

!*****************************************************************************80
!
!! HF_POWER_PRODUCT: power products for Hermite function Hf(n,x).
!
!  Discussion:
!
!    Let Hf(I,X) represent the Hermite function of "degree" I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) X^E * Hf(I,X) * Hf(J,X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!
!    When E is 0, the computed table should be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
!    TABLE(I,J) represents the integral of X^E * Hf(I,X) * Hf(J,X).
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

  call hf_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hf_function ( 1, p, x, h_table )
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
subroutine hf_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! HF_QUADRATURE_RULE: quadrature for Hermite function Hf(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( pi ) )

  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt) ** 2 * exp ( t(1:nt) ** 2 )

  return
end
subroutine hn_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! HN_EXPONENTIAL_PRODUCT: exponential products for Hn(n,x).
!
!  Discussion:
!
!    Let Hn(I,X) represent the normalized physicist's Hermite polynomial 
!    of degree I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) 
!        exp(B*X) * Hn(I,X) * Hn(J,X) exp(-X*X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!    Because of the exponential factor exp(B*X), the quadrature will not 
!    be exact.
!
!    However, when B = 0, the quadrature is exact, and moreoever, the
!    table will be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
!    TABLE(I,J) represents the weighted integral of 
!    exp(B*X) * Hn(I,X) * Hn(J,X).
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

  call h_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hn_polynomial ( 1, p, x, h_table )
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
subroutine hn_polynomial ( m, n, x, p )

!*****************************************************************************80
!
!! HN_POLYNOMIAL evaluates normalized physicist's Hermite polynomials Hn(n,x).
!
!  Discussion:
!
!    These polynomials satisfy the orthonormality condition:
!
!      Integral ( -oo < X < +oo ) 
!        exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
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
!    Output, real ( kind = 8 ) P(M,0:N), the values of the polynomials of 
!    index 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fact
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(m,0:n)
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) two
  real ( kind = 8 ) x(m)

  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = 2.0D+00 * x(1:m)
 
  do j = 2, n
    p(1:m,j) = 2.0D+00 * x(1:m) * p(1:m,j-1) &
      - 2.0D+00 * real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
!
!  Normalize.
!
  fact = 1.0D+00
  two = 1.0D+00
  do j = 0, n
    p(1:m,j) = p(1:m,j) / sqrt ( fact * two * sqrt ( pi ) )
    fact = fact * real ( j + 1, kind = 8 )
    two = two * 2.0D+00
  end do

  return
end
subroutine hn_power_product ( p, e, table )

!*****************************************************************************80
!
!! HN_POWER_PRODUCT: power products for normalized physicist's Hermite Hn(n,x).
!
!  Discussion:
!
!    Let Hn(I,X) represent the normalized physicist's Hermite polynomial 
!    of degree I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) X^E * Hn(I,X) * Hn(J,X) exp(-X*X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!
!    When E is 0, the computed table should be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
!    TABLE(I,J) represents the weighted integral of X^E * Hn(I,X) * Hn(J,X).
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

  call h_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hn_polynomial ( 1, p, x, h_table )
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
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    N!!
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial
!    function.  If N is less than 1, the value is returned as 1.0.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL2, the value.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial2
  real ( kind = 8 ) r8_n

  if ( n < 1 ) then
    r8_factorial2 = 1.0D+00
    return
  end if

  r8_n = real ( n, kind = 8 )
  r8_factorial2 = 1.0D+00

  do while ( 1.0D+00 < r8_n )
    r8_factorial2 = r8_factorial2 * r8_n
    r8_n = r8_n - 2.0D+00
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
!    value =  0 if X => 0.
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
