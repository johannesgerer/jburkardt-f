subroutine hermite_cubic_integral ( x1, f1, d1, x2, f2, d2, q )

!*****************************************************************************80
!
!! HERMITE_CUBIC_INTEGRAL returns the integral of a Hermite cubic polynomial.
!
!  Discussion:
!
!    The integral is taken over the definition interval [X1,X2].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, F1, D1, the left endpoint, function value
!    and derivative.
!
!    Input, real ( kind = 8 ) X2, F2, D2, the right endpoint, function value
!    and derivative.
!
!    Output, real ( kind = 8 ) Q, the integral of the Hermite cubic polynomial
!    over the interval X1 <= X <= X2.
!
  implicit none

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) h
  real ( kind = 8 ) q
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h = x2 - x1

  q = 0.5D+00 * h * ( f1 + f2 + h * ( d1 - d2 ) / 6.0D+00 )

  return
end
subroutine hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, q )

!*****************************************************************************80
!
!! HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic polynomial from A to B.
!
!  Discussion:
!
!    A and B may be scalars, or one may be a vector, or both
!    may be vectors of the same size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, F1, D1, the left endpoint, function value
!    and derivative.
!
!    Input, real ( kind = 8 ) X2, F2, D2, the right endpoint, function value
!    and derivative.
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints of the interval
!    of integration.
!
!    Output, real ( kind = 8 ) Q, the integral of the Hermite cubic polynomial
!    over the interval A <= X <= B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dterm
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fterm
  real ( kind = 8 ) h
  real ( kind = 8 ) phia1
  real ( kind = 8 ) phia2
  real ( kind = 8 ) phib1
  real ( kind = 8 ) phib2
  real ( kind = 8 ) psia1
  real ( kind = 8 ) psia2
  real ( kind = 8 ) psib1
  real ( kind = 8 ) psib2
  real ( kind = 8 ) q
  real ( kind = 8 ) ta1
  real ( kind = 8 ) ta2
  real ( kind = 8 ) tb1
  real ( kind = 8 ) tb2
  real ( kind = 8 ) ua1
  real ( kind = 8 ) ua2
  real ( kind = 8 ) ub1
  real ( kind = 8 ) ub2
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h = x2 - x1

  ta1 = ( a - x1 ) / h
  ta2 = ( x2 - a ) / h
  tb1 = ( b - x1 ) / h
  tb2 = ( x2 - b ) / h

  ua1 = ta1 * ta1 * ta1
  phia1 = ua1 * ( 2.0D+00 - ta1 )
  psia1 = ua1 * ( 3.0D+00 * ta1 - 4.0D+00 )

  ua2 = ta2 * ta2 * ta2
  phia2 =  ua2 * ( 2.0D+00 - ta2 )
  psia2 = -ua2 * ( 3.0D+00 * ta2 - 4.0D+00 )

  ub1 = tb1 * tb1 * tb1
  phib1 = ub1 * ( 2.0D+00 - tb1 )
  psib1 = ub1 * ( 3.0D+00 * tb1 - 4.0D+00 )

  ub2 = tb2 * tb2 * tb2
  phib2 =  ub2 * ( 2.0D+00 - tb2 )
  psib2 = -ub2 * ( 3.0D+00 * tb2 - 4.0D+00 )

  fterm =   f1 * ( phia2 - phib2 ) + f2 * ( phib1 - phia1 )
  dterm = ( d1 * ( psia2 - psib2 ) + d2 * ( psib1 - psia1 ) ) * ( h / 6.0D+00 )

  q = 0.5D+00 * h * ( fterm + dterm )

  return
end
subroutine hermite_cubic_lagrange_integral ( x1, x2, q )

!*****************************************************************************80
!
!! HERMITE_CUBIC_LAGRANGE_INTEGRAL: Hermite cubic Lagrange integrals.
!
!  Discussion:
!
!    The Hermite cubic polynomial P(X) for interval (X1,X2) and data
!    (F1,D1,F2,D2) satisfies:
!
!      P(X1) = F1,
!      P'(X1) = D1,
!      P(X2) = F2,
!      P'(X2) = D2.
!
!    We can determine four Lagrange polynomials L1(X) through L4(X) so that
!
!      P(X) = F1 * L1(X) + D1 * L2(X) + F2 * L3(X) + D2 * L4(X).
!
!    This function returns the integrals of these four polynomials over
!    the domain of definition [X1,X2].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints.
!
!    Output, real ( kind = 8 ) Q(4), the integrals of the Hermite cubic
!    Lagrange polynomials from X1 to X2.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) q(4)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h = x2 - x1

  q(1) =   h     / 2.0D+00
  q(2) =   h * h / 12.0D+00
  q(3) =   h     / 2.0D+00
  q(4) = - h * h / 12.0D+00

  return
end
subroutine hermite_cubic_lagrange_integrate ( x1, x2, a, b, q )

!*****************************************************************************80
!
!! HERMITE_CUBIC_LAGRANGE_INTEGRATE: integrate Hermite cubic Lagrange polys.
!
!  Discussion:
!
!    A and B may be scalars, or one may be a vector, or both
!    may be vectors of the same size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real X1, X2, the endpoints of the interval of definition.
!
!    Input, real A, B, the left and right endpoints of the interval
!    of integration.
!
!    Output, real Q(4), the integrals of the Hermite cubic Lagrange polynomials
!    over the interval A <= X <= B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) h
  real ( kind = 8 ) phia1
  real ( kind = 8 ) phia2
  real ( kind = 8 ) phib1
  real ( kind = 8 ) phib2
  real ( kind = 8 ) psia1
  real ( kind = 8 ) psia2
  real ( kind = 8 ) psib1
  real ( kind = 8 ) psib2
  real ( kind = 8 ) q(4)
  real ( kind = 8 ) ta1
  real ( kind = 8 ) ta2
  real ( kind = 8 ) tb1
  real ( kind = 8 ) tb2
  real ( kind = 8 ) ua1
  real ( kind = 8 ) ua2
  real ( kind = 8 ) ub1
  real ( kind = 8 ) ub2
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h = x2 - x1
  ta1 = ( a - x1 ) / h
  ta2 = ( x2 - a ) / h
  tb1 = ( b - x1 ) / h
  tb2 = ( x2 - b ) / h

  ua1 = ta1 * ta1 * ta1
  phia1 = ua1 * ( 2.0D+00 - ta1 )
  psia1 = ua1 * ( 3.0D+00 * ta1 - 4.0D+00 )

  ua2 = ta2 * ta2 * ta2
  phia2 =  ua2 * ( 2.0D+00 - ta2 )
  psia2 = -ua2 * ( 3.0D+00 * ta2 - 4.0D+00 )

  ub1 = tb1 * tb1 * tb1
  phib1 = ub1 * ( 2.0D+00 - tb1 )
  psib1 = ub1 * ( 3.0D+00 * tb1 - 4.0D+00 )

  ub2 = tb2 * tb2 * tb2
  phib2 =  ub2 * ( 2.0D+00 - tb2 )
  psib2 = -ub2 * ( 3.0D+00 * tb2 - 4.0D+00 )

  q(1) = 0.5D+00 * h * ( phia2 - phib2 )
  q(2) = 0.5D+00 * h * ( psia2 - psib2 ) * ( h / 6.0D+00 )
  q(3) = 0.5D+00 * h * ( phib1 - phia1 )
  q(4) = 0.5D+00 * h * ( psib1 - psia1 ) * ( h / 6.0D+00 )

  return
end
subroutine hermite_cubic_lagrange_value ( x1, x2, n, x, f, d, s, t )

!*****************************************************************************80
!
!! HERMITE_CUBIC_LAGRANGE_VALUE: evaluate Hermite cubic Lagrange polynomials.
!
!  Discussion:
!
!    The Hermite cubic polynomial P(X) for interval (X1,X2) and data
!    (F1,D1,F2,D2) satisfies:
!
!      P(X1) = F1,
!      P'(X1) = D1,
!      P(X2) = F2,
!      P'(X2) = D2.
!
!    We can determine four Lagrange polynomials L1(X) through L4(X) so that
!
!      P(X) = F1 * L1(X) + D1 * L2(X) + F2 * L3(X) + D2 * L4(X).
!
!    This function returns the values and derivatives of these four
!    polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints.
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!
!    Input, real ( kind = 8 ) X(N), the sample points.
!
!    Output, real ( kind = 8 ) F(4,N), D(4,N), S(4,N), T(4,N), the value
!    and first three derivatives of the Hermite cubic Lagrange polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d(4,n)
  real ( kind = 8 ) dx(n)
  real ( kind = 8 ) f(4,n)
  real ( kind = 8 ) h
  real ( kind = 8 ) s(4,n)
  real ( kind = 8 ) t(4,n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h = x2 - x1
  dx(1:n) = x(1:n) - x1
!
!  F1.
!
  f(1,1:n) = 1.0D+00 &
           + ( dx(1:n)**2 / h**2 ) * ( - 3.0D+00 + ( dx(1:n) / h ) *  2.0D+00 )
  d(1,1:n) = ( dx(1:n)    / h**2 ) * ( - 6.0D+00 + ( dx(1:n) / h ) *  6.0D+00 )
  s(1,1:n) = ( 1.0D+00    / h**2 ) * ( - 6.0D+00 + ( dx(1:n) / h ) * 12.0D+00 )
  t(1,1:n) = ( 1.0D+00    / h**3 )                                 * 12.0D+00
!
!  D1
!
  f(2,1:n) = dx(1:n) + ( dx(1:n)**2   / h    ) * ( - 2.0D+00 + ( dx(1:n) / h )           )
  d(2,1:n) = 1.0D+00 + ( dx(1:n)      / h    ) * ( - 4.0D+00 + ( dx(1:n) / h ) * 3.0D+00 )
  s(2,1:n) =           ( 1.0D+00      / h    ) * ( - 4.0D+00 + ( dx(1:n) / h ) * 6.0D+00 )
  t(2,1:n) =           ( 1.0D+00      / h**2 )                                 * 6.0D+00
!
!  F2
!
  f(3,1:n) = ( dx(1:n)**2 / h**2 ) * ( 3.0D+00 -  2.0D+00 * ( dx(1:n) / h ) )
  d(3,1:n) = ( dx(1:n)    / h**2 ) * ( 6.0D+00 -  6.0D+00 * ( dx(1:n) / h ) )
  s(3,1:n) = ( 1.0D+00    / h**2 ) * ( 6.0D+00 - 12.0D+00 * ( dx(1:n) / h ) )
  t(3,1:n) = ( 1.0D+00    / h**3 ) * (         - 12.0D+00              )
!
!  D2
!
  f(4,1:n) = ( dx(1:n)**2 / h ) * ( - 1.0D+00 + ( dx(1:n) / h )           )
  d(4,1:n) = ( dx(1:n)    / h ) * ( - 2.0D+00 + ( dx(1:n) / h ) * 3.0D+00 )
  s(4,1:n) = ( 1.0D+00    / h ) * ( - 2.0D+00 + ( dx(1:n) / h ) * 6.0D+00 )
  t(4,1:n) = ( 1.0D+00    / h )                                 * 6.0D+00

  return
end
subroutine hermite_cubic_spline_integral ( nn, xn, fn, dn, q )

!*****************************************************************************80
!
!! HERMITE_CUBIC_SPLINE_INTEGRAL: Hermite cubic spline integral.
!
!  Discussion:
!
!    The integral is taken over the definition interval [X(1),X(NN)].
!
!    Note that if the intervals are equal in size, then the derivative
!    information in DN has no effect on the integral value,
!    except for the first and last entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NN, the number of data points.
!
!    Input, real ( kind = 8 ) XN(NN), the coordinates of the data points.
!    The entries in XN must be in strictly ascending order.
!
!    Input, real ( kind = 8 ) FN(NN), the function values.
!
!    Input, real ( kind = 8 ) DN(NN), the derivative values.
!
!    Output, real ( kind = 8 ) Q, the integral of the Hermite cubic spline
!    over the interval X(1) <= X <= X(NN).
!
  implicit none

  integer ( kind = 4 ) nn

  real ( kind = 8 ) dn(nn)
  real ( kind = 8 ) fn(nn)
  real ( kind = 8 ) q
  real ( kind = 8 ) xn(nn)

  q = sum ( &
    0.5D+00 * ( xn(2:nn) - xn(1:nn-1) ) * ( fn(1:nn-1) + fn(2:nn) &
    + ( xn(2:nn) - xn(1:nn-1) ) * ( dn(1:nn-1) - dn(2:nn) ) / 6.0D+00 ) )

  return
end
subroutine hermite_cubic_spline_quad_rule ( nn, xn, w )

!*****************************************************************************80
!
!! HERMITE_CUBIC_SPLINE_QUAD_RULE: Hermite cubic spline quadrature rule.
!
!  Discussion:
!
!    The integral is taken over the definition interval [X(1),X(NN)].
!
!    Note that if the intervals are equal in size, then the derivative
!    information in DN has no effect on the integral value,
!    except for the first and last entries.
!
!    The quadrature rule is
!
!      Integral ( XN(1) <= X <= XN(NN) ) F(X) dX is approximately
!
!      Sum ( 1 <= I <= NN ) W(1,I) * F(X(I)) + W(2,I) * F'(X(I))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NN, the number of data points.
!
!    Input, real ( kind = 8 ) XN(NN), the coordinates of the data points.
!    The entries in XN must be in strictly ascending order.
!
!    Output, real ( kind = 8 ) W(2,NN), the quadrature weights for F(1:NN)
!    and DN(1:NN).
!
  implicit none

  integer ( kind = 4 ) nn

  real ( kind = 8 ) w(2,nn)
  real ( kind = 8 ) xn(nn)

  w(1,1)      = 0.5D+00 * ( xn(2)    - xn(1)      )
  w(1,2:nn-1) = 0.5D+00 * ( xn(3:nn) - xn(1:nn-2) )
  w(1,nn)     = 0.5D+00 * ( xn(nn)   - xn(nn-1)   )

  w(2,1)      =   ( xn(2) - xn(1) )**2 / 12.0D+00
  w(2,2:nn-1) =   ( xn(3:nn) - xn(1:nn-2) ) &
                * ( xn(3:nn) - 2.0D+00 * xn(2:nn-1) + xn(1:nn-2) ) / 12.0D+00
  w(2,nn)     = - ( xn(nn-1) - xn(nn) )**2 / 12.0D+00

  return
end
subroutine hermite_cubic_spline_integrate ( nn, xn, fn, dn, n, a, b, q )

!*****************************************************************************80
!
!! HERMITE_CUBIC_SPLINE_INTEGRATE integrates a Hermite cubic spline over [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NN, the number of data points.
!
!    Input, real ( kind = 8 ) XN(NN), the coordinates of the data points.
!    The entries in XN must be in strictly ascending order.
!
!    Input, real ( kind = 8 ) FN(NN), the function values.
!
!    Input, real ( kind = 8 ) DN(NN), the derivative values.
!
!    Input, integer ( kind = 4 ) N, the number of integration intervals.
!
!    Input, real ( kind = 8 ) A(N), B(N), the integration endpoints.
!
!    Output, real ( kind = 8 ) Q(N), the integral over the interval [A,B].
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nn

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aa
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) bb
  real ( kind = 8 ) dn(nn)
  real ( kind = 8 ) fn(nn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) q(n)
  real ( kind = 8 ) qq
  real ( kind = 8 ) s
  real ( kind = 8 ) xn(nn)

  q(1:n) = 0.0D+00

  do ii = 1, n

    if ( a(ii) <= b(ii) ) then
      aa = a(ii)
      bb = b(ii)
      s = + 1.0D+00
    else
      aa = b(ii)
      bb = a(ii)
      s = - 1.0D+00
    end if

    call r8vec_bracket3 ( nn, xn, aa, i )
    call r8vec_bracket3 ( nn, xn, bb, j )
!
!  Evaluate the polynomial with the appropriate data.
!
    if ( i == j ) then

      call hermite_cubic_integrate ( xn(i), fn(i), dn(i), &
        xn(i+1), fn(i+1), dn(i+1), aa, bb, q(ii) )

    else

      call hermite_cubic_integrate ( xn(i), fn(i), dn(i), &
        xn(i+1), fn(i+1), dn(i+1), aa, xn(i+1), qq )

      q(ii) = qq

      do k = i + 1, j - 1

        call hermite_cubic_integral ( xn(k), fn(k), dn(k), &
          xn(k+1), fn(k+1), dn(k+1), qq )

        q(ii) = q(ii) + qq

      end do

      call hermite_cubic_integrate ( xn(j), fn(j), dn(j), &
        xn(j+1), fn(j+1), dn(j+1), xn(j), bb, qq )

      q(ii) = q(ii) + qq

    end if

    q(ii) = s * q(ii)

  end do

  return
end
subroutine hermite_cubic_spline_value ( nn, xn, fn, dn, n, x, f, d, s, t )

!*****************************************************************************80
!
!! HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NN, the number of data points.
!
!    Input, real ( kind = 8 ) XN(NN), the coordinates of the data points.
!    The entries in XN must be in strictly ascending order.
!
!    Input, real ( kind = 8 ) FN(NN), the function values.
!
!    Input, real ( kind = 8 ) DN(NN), the derivative values.
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!
!    Input, real ( kind = 8 ) X(N), the coordinates of the sample points.
!
!    Output, real ( kind = 8 ) F(N), the function value at the sample points.
!
!    Output, real ( kind = 8 ) D(N), the derivative value at the sample points.
!
!    Output, real ( kind = 8 ) S(N), the second derivative value at the
!    sample points.
!
!    Output, real ( kind = 8 ) T(N), the third derivative value at the
!    sample points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nn

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dn(nn)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fn(nn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xn(nn)

  do i = 1, n

    call r8vec_bracket3 ( nn, xn, x(i), left )

    call hermite_cubic_value ( xn(left), fn(left), dn(left), xn(left+1), &
      fn(left+1), dn(left+1), 1, x(i), f(i), d(i), s(i), t(i) )

  end do

  return
end
subroutine hermite_cubic_to_power_cubic ( x1, f1, d1, x2, f2, d2, c0, c1, &
  c2, c3 )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TO_POWER_CUBIC converts a Hermite cubic to power form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, F1, D1, the left endpoint, function value
!    and derivative.
!
!    Input, real ( kind = 8 ) X2, F2, D2, the right endpoint, function value
!    and derivative.
!
!    Output, real ( kind = 8 ) C0, C1, C2, C3, the power form of the polynomial.
!
  implicit none

  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) df
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) h
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h =    x2 - x1
  df = ( f2 - f1 ) / h
!
!  Polynomial in terms of X - X1:
!
  c0 = f1
  c1 = d1
  c2 = - ( 2.0D+00 * d1 - 3.0D+00 * df + d2 ) / h
  c3 =   (           d1 - 2.0D+00 * df + d2 ) / h / h
!
!  Shift polynomial to X.
!
  c2 = c2 - x1 * c3
  c1 = c1 - x1 * c2
  c0 = c0 - x1 * c1
  c2 = c2 - x1 * c3
  c1 = c1 - x1 * c2
  c2 = c2 - x1 * c3

  return
end
subroutine hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, d, s, t )

!*****************************************************************************80
!
!! HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.
!
!  Discussion:
!
!    The input arguments can be vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, F1, D1, the left endpoint, function value
!    and derivative.
!
!    Input, real ( kind = 8 ) X2, F2, D2, the right endpoint, function value
!    and derivative.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the points at which the Hermite cubic
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), D(N), S(N), T(N), the value and first
!    three derivatives of the Hermite cubic at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) df
  real ( kind = 8 ) dx(n)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) h
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h =    x2 - x1
  df = ( f2 - f1 ) / h

  c2 = - ( 2.0 * d1 - 3.0 * df + d2 ) / h
  c3 =   (       d1 - 2.0 * df + d2 ) / h / h

  dx(1:n) = x(1:n) - x1
  f(1:n) = f1 + dx(1:n) &
             * ( d1 + dx(1:n) * (           c2 + dx(1:n) *           c3 ) )
  d(1:n) =       d1 + dx(1:n) * ( 2.0D+00 * c2 + dx(1:n) * 3.0D+00 * c3 )
  s(1:n) =                        2.0D+00 * c2 + dx(1:n) * 6.0D+00 * c3
  t(1:n) =                                                 6.0D+00 * c3

  return
end
subroutine power_cubic_to_hermite_cubic ( c0, c1, c2, c3, x1, x2, f1, d1, &
  f2, d2 )

!*****************************************************************************80
!
!! POWER_CUBIC_TO_HERMITE_CUBIC converts a power cubic to Hermite form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C0, C1, C2, C3, the power form of the
!    polynomial.
!
!    Input, real ( kind = 8 ) X1, X2, the left and right endpoints of
!    the Hermite form.
!
!    Output, real ( kind = 8 ) F1, D1, the function and derivative values at X1.
!
!    Output, real ( kind = 8 ) F2, D2, the function and derivative values at X2.
!
  implicit none

  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  f1 = c0 + x1 * ( c1 + x1 * (           c2 + x1           * c3 ) )
  d1 =             c1 + x1 * ( 2.0D+00 * c2 + x1 * 3.0D+00 * c3 )

  f2 = c0 + x2 * ( c1 + x2 * (           c2 + x2           * c3 ) )
  d2 =             c1 + x2 * ( 2.0D+00 * c2 + x2 * 3.0D+00 * c3 )

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
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
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
  real    ( kind = 8 ) r8_uniform_01
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
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8vec_bracket3 ( n, t, tval, left )

!*****************************************************************************80
!
!! R8VEC_BRACKET3 finds the interval containing or nearest a given value.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of the input array.
!
!    Input, real ( kind = 8 ) T(N), an array that has been sorted
!    into ascending order.
!
!    Input, real ( kind = 8 ) TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer ( kind = 4 ) LEFT.
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
!
!  Check the input data.
!
  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  if ( left < 1 .or. n - 1 < left ) then
    left = ( n + 1 ) / 2
  end if
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( t(left-1) <= tval ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  else if ( t(left+1) < tval ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( t(n-1) <= tval ) then
      left = n - 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  else

  end if

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) ahi
  real    ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
