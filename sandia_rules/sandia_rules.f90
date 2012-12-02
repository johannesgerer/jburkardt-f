subroutine binary_vector_next ( n, bvec )

!*****************************************************************************80
!
!! BINARY_VECTOR_NEXT generates the next binary vector.
!
!  Discussion:
!
!    A binary vector is a vector whose entries are 0 or 1.
!
!    The user inputs an initial zero vector to start.  The program returns
!    the "next" vector.
!
!    The vectors are produced in the order:
!
!    ( 0, 0, 0, ..., 0 )
!    ( 1, 0, 0, ..., 0 )
!    ( 0, 1, 0, ..., 0 )
!    ( 1, 1, 0, ..., 0 )
!    ( 0, 0, 1, ..., 0 )
!    ( 1, 0, 1, ..., 0 )
!               ...
!    ( 1, 1, 1, ..., 1)
!
!    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
!    we allow wrap around.
!
!  Example:
!
!    N = 3
!
!    Input      Output
!    -----      ------
!    0 0 0  =>  1 0 0
!    1 0 0  =>  0 1 0
!    0 1 0  =>  1 1 0
!    1 1 0  =>  0 0 1
!    0 0 1  =>  1 0 1
!    1 0 1  =>  0 1 1
!    0 1 1  =>  1 1 1
!    1 1 1  =>  0 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input/output, integer ( kind = 4 ) BVEC(N), on output, the successor
!    to the input vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i

  do i = 1, n

    if ( bvec(i) == 1 ) then
      bvec(i) = 0
    else
      bvec(i) = 1
      exit
    end if

  end do

  return
end
subroutine ccn_compute ( n, x, w )

!*****************************************************************************80
!
!! CCN_COMPUTE computes a nested Clenshaw Curtis quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call ccn_compute_points ( n, x )
  call ccn_compute_weights ( n, w )

  return
end
subroutine ccn_compute_np ( n, np, p, x, w )

!*****************************************************************************80
!
!! CCN_COMPUTE_NP computes a nested Clenshaw Curtis rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the points.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call ccn_compute_points ( n, x )
  call ccn_compute_weights ( n, w )

  return
end
subroutine ccn_compute_points ( n, x )

!*****************************************************************************80
!
!! CCN_COMPUTE_POINTS: compute nested Clenshaw Curtis points.
!
!  Discussion:
!
!    We want to compute the following sequence:
!
!    1/2,
!    0, 1
!    1/4, 3/4
!    1/8, 3/8, 5/8, 7/8,
!    1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.
!
!    But we'd prefer that the numbers in each row be regrouped in pairs
!    that are symmetric about 1/2, with the number above 1/2 coming first.
!    Thus, the last row might become:
!    (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).
!
!    Once we have our sequence, we apply the Chebyshev transformation
!    which maps [0,1] to [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements to compute.
!
!    Output, real ( kind = 8 ) X(N), the elements of the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) td
  integer ( kind = 4 ) tu
  real ( kind = 8 ) x(n)
!
!  Handle first three entries specially.
!
  if ( 1 <= n ) then
    x(1) = 0.5D+00
  end if

  if ( 2 <= n ) then
    x(2) = 1.0D+00
  end if

  if ( 3 <= n ) then
    x(3) = 0.0D+00
  end if

  m = 3
  d = 2

  do while ( m < n )

    tu = d + 1
    td = d - 1

    k = min ( d, n - m )

    do i = 1, k
      if ( mod ( i, 2 ) == 1 ) then
        x(m+i) = real ( tu, kind = 8 ) / 2.0D+00 / real ( k, kind = 8 )
        tu = tu + 2
      else
        x(m+i) = real ( td, kind = 8 ) / 2.0D+00 / real ( k, kind = 8 )
        td = td - 2
      end if
    end do

    m = m + k
    d = d * 2

  end do
!
!  Apply the Chebyshev transformation.
!
  x(1:n) = cos ( x(1:n) * pi )

  x(1) = 0.0D+00

  if ( 2 <= n ) then
    x(2) = -1.0D+00
  end if

  if ( 3 <= n ) then
    x(3) = +1.0D+00
  end if

  return
end
subroutine ccn_compute_points_np ( n, np, p, x )

!*****************************************************************************80
!
!! CCN_COMPUTE_POINTS_NP: abscissas of a nested Clenshaw Curtis rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) x(n)

  call ccn_compute_points ( n, x )

  return
end
subroutine ccn_compute_weights ( n, w )

!*****************************************************************************80
!
!! CCN_COMPUTE_WEIGHTS: weights for nested Clenshaw Curtis rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min

  call ccn_compute_points ( n, x )
!
!  Get the weights.
!
  x_min = -1.0D+00
  x_max = +1.0D+00

  call nc_compute ( n, x_min, x_max, x, w )

  return
end
subroutine ccn_compute_weights_np ( n, np, p, w )

!*****************************************************************************80
!
!! CCN_COMPUTE_WEIGHTS_NP computes nested Clenshaw Curtis weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) w(n)

  call ccn_compute_weights ( n, w )

  return
end
subroutine chebyshev1_compute ( n, x, w )

!*****************************************************************************80
!
!! CHEBYSHEV1_COMPUTE computes a Chebyshev type 1 quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0 / sqrt ( 1 - x^2 ).
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) / sqrt ( 1 - x^2 ) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEBYSHEV1_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop
  end if

  w(1:n) = pi / real ( n, kind = 8 )

  do i = 1, n
    x(i) = cos ( pi * real ( 2 * n + 1 - 2 * i, kind = 8 ) &
                    / real ( 2 * n, kind = 8 ) )
  end do

  return
end
subroutine chebyshev1_integral ( expon, exact )

!*****************************************************************************80
!
!! CHEBYSHEV1_INTEGRAL evaluates a monomial Chebyshev type 1 integral.
!
!  Discussion:
!
!    To test a Chebyshev type 1 quadrature rule, we use it to approximate the
!    integral of a monomial:
!
!      integral ( -1 <= x <= +1 ) x^n / sqrt ( 1 - x^2 ) dx
!
!    This routine is given the value of the exponent, and returns the
!    exact value of the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Output, real ( kind = 8 ) EXACT, the value of the exact integral.
!
  implicit none

  real ( kind = 8 ) bot
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) top
!
!  Get the exact value of the integral.
!
  if ( mod ( expon, 2 ) == 0 ) then

    top = 1
    bot = 1
    do i = 2, expon, 2
      top = top * ( i - 1 )
      bot = bot *   i
    end do

    exact = pi * real ( top, kind = 8 ) / real ( bot, kind = 8 )

  else

    exact = 0.0D+00

  end if

  return
end
subroutine chebyshev2_compute ( n, x, w )

!*****************************************************************************80
!
!! CHEBYSHEV2_COMPUTE computes a Chebyshev type 2 quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = sqrt ( 1 - x^2 ).
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X)  sqrt ( 1 - x^2 )  dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) w(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEBYSHEV2_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop
  end if

  do i = 1, n
    angle = pi * real ( n + 1 - i, kind = 8 ) / real ( n + 1, kind = 8 )
    w(i) = pi / real ( n + 1, kind = 8 ) * ( sin ( angle ) )**2
    x(i) = cos ( angle )
  end do

  return
end
subroutine chebyshev2_integral ( expon, exact )

!*****************************************************************************80
!
!! CHEBYSHEV2_INTEGRAL evaluates a monomial Chebyshev type 2 integral.
!
!  Discussion:
!
!    To test a Chebyshev type 2 quadrature rule, we use it to approximate the
!    integral of a monomial:
!
!      integral ( -1 <= x <= +1 ) x^n * sqrt ( 1 - x^2 ) dx
!
!    This routine is given the value of the exponent, and returns the
!    exact value of the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Output, real ( kind = 8 ) EXACT, the value of the exact integral.
!
  implicit none

  real ( kind = 8 ) bot
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) top
!
!  Get the exact value of the integral.
!
  if ( mod ( expon, 2 ) == 0 ) then

    top = 1
    bot = 1
    do i = 2, expon, 2
      top = top * ( i - 1 )
      bot = bot *   i
    end do

    bot = bot * real ( expon + 2, kind = 8 )

    exact = pi * real ( top, kind = 8 ) / real ( bot, kind = 8 )

  else

    exact = 0.0D+00

  end if

  return
end
subroutine clenshaw_curtis_compute ( n, x, w )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop
  end if

  if ( n == 1 ) then
    x(1) = 0.0D+00
    w(1) = 2.0D+00
    return
  end if

  do i = 1, n
    x(i) = cos ( real ( n - i, kind = 8 ) * pi &
               / real ( n - 1, kind = 8 ) )
  end do

  x(1) = -1.0D+00
  if ( mod ( n, 2 ) == 1 ) then
    x((n+1)/2) = 0.0D+00
  end if
  x(n) = +1.0D+00

  do i = 1, n

    theta = real ( i - 1, kind = 8 ) * pi &
          / real ( n - 1, kind = 8 )

    w(i) = 1.0D+00

    do j = 1, ( n - 1 ) / 2

      if ( 2 * j == ( n - 1 ) ) then
        b = 1.0D+00
      else
        b = 2.0D+00
      end if

      w(i) = w(i) - b * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta ) &
           / real ( 4 * j * j - 1, kind = 8 )

    end do

  end do

  w(1)     =           w(1)     / real ( n - 1, kind = 8 )
  w(2:n-1) = 2.0D+00 * w(2:n-1) / real ( n - 1, kind = 8 )
  w(n)     =           w(n)     / real ( n - 1, kind = 8 )

  return
end
subroutine clenshaw_curtis_compute_points ( n, points )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE_POINTS: abscissas of a Clenshaw Curtis rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) POINTS(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) points(n)

  if ( n < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  else if ( n == 1 ) then

    points(1) = 0.0D+00

  else

    do i = 1, n
      points(i) = cos ( real ( n - i, kind = 8 ) * pi &
                      / real ( n - 1, kind = 8 ) )
    end do

    points(1) = -1.0D+00
    if ( mod ( n, 2 ) == 1 ) then
      points((n+1)/2) = 0.0D+00
    end if
    points(n) = +1.0D+00

  end if

  return
end
subroutine clenshaw_curtis_compute_points_np ( n, np, p, x )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE_POINTS_NP: abscissas of a Clenshaw Curtis rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) x(n)

  call clenshaw_curtis_compute_points ( n, x )

  return
end
subroutine clenshaw_curtis_compute_weights ( n, w )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE_WEIGHTS computes Clenshaw Curtis weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) w(n)

  if ( n < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  else if ( n == 1 ) then

    w(1) = 2.0D+00

  else

    do i = 1, n

      theta = real ( i - 1, kind = 8 ) * pi &
            / real ( n - 1, kind = 8 )

      w(i) = 1.0D+00

      do j = 1, ( n - 1 ) / 2

        if ( 2 * j == ( n - 1 ) ) then
          b = 1.0D+00
        else
          b = 2.0D+00
        end if

        w(i) = w(i) - b * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta ) &
             / real ( 4 * j * j - 1, kind = 8 )

      end do

    end do

    w(1)     =           w(1)     / real ( n - 1, kind = 8 )
    w(2:n-1) = 2.0D+00 * w(2:n-1) / real ( n - 1, kind = 8 )
    w(n)     =           w(n)     / real ( n - 1, kind = 8 )

  end if

  return
end
subroutine clenshaw_curtis_compute_weights_np ( n, np, p, w )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE_WEIGHTS_NP computes Clenshaw Curtis weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) w(n)

  call clenshaw_curtis_compute_weights ( n, w )

  return
end
subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose compositions are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer ( kind = 4 )  H, T, two internal parameters needed
!    for the computation.  The user should allocate space for these in the
!    calling program, include them in the calling sequence, but never alter
!    them!
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) h
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else
!
!  If the first entry A(1) is positive, then set H to zero,
!  so that when we increment H, it points to A(1); we will decrement A(1) by 1
!  and increment A(2).
!
    if ( 1 < t ) then
      h = 0
    end if
!
!  Otherwise, A(1) is 0.  Then by H + 1 is the entry we incremented last time.
!  Set H = H + 1, zero A(H), adding all but one of its value to A(1),
!  and incrementing A(H+1) by 1.
!
    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

  return
end
subroutine dif_deriv ( nd, xd, yd, ndp, xdp, ydp )

!*****************************************************************************80
!
!! DIF_DERIV computes the derivative of a polynomial in divided difference form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the size of the input table.
!
!    Input, real ( kind = 8 ) XD(ND), the abscissas for the divided
!    difference table.
!
!    Input, real ( kind = 8 ) YD(ND), the divided difference table.
!
!    Output, integer ( kind = 4 ) NDP, the size of the output table, 
!    which is ND-1.
!
!    Input, real ( kind = 8 ) XDP(NDP), the abscissas for the divided
!    difference table for the derivative.
!
!    Output, real ( kind = 8 ) YDP(NDP), the divided difference 
!    table for the derivative.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ndp
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xd_temp(nd)
  real ( kind = 8 ) xdp(nd-1)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yd_temp(nd)
  real ( kind = 8 ) ydp(nd-1)
!
!  Using a temporary copy of the difference table, shift the abscissas to zero.
!
  xd_temp(1:nd) = xd(1:nd)
  yd_temp(1:nd) = yd(1:nd)

  call dif_shift_zero ( nd, xd_temp, yd_temp )
!
!  Now construct the derivative.
!
  ndp = nd - 1

  xdp(1:ndp) = 0.0D+00

  do i = 1, ndp
    ydp(i) = real ( i, kind = 8 ) * yd_temp(i+1)
  end do

  return
end
subroutine dif_shift_x ( nd, xd, yd, xv )

!*****************************************************************************80
!
!! DIF_SHIFT_X replaces one abscissa of a divided difference table.
!
!  Discussion:
!
!    The routine shifts the representation of a divided difference polynomial by
!    dropping the last X value in XD, and adding a new value XV to the
!    beginning of the XD array, suitably modifying the coefficients stored
!    in YD.
!
!    The representation of the polynomial is changed, but the polynomial itself
!    should be identical.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of divided difference 
!    coefficients, and the number of entries in XD.
!
!    Input/output, real ( kind = 8 ) XD(ND), the X values used in 
!    the representation of the divided difference polynomial.
!    After a call to this routine, the last entry of XD has been dropped,
!    the other entries have shifted up one index, and XV has been inserted
!    at the beginning of the array.
!
!    Input/output, real ( kind = 8 ) YD(ND), the divided difference
!    coefficients corresponding to the XD array.  On output, this 
!    array has been adjusted.
!
!    Input, real ( kind = 8 ) XV, a new X value which is to be used in 
!    the representation of the polynomial.  On output, XD(1) equals 
!    XV and the representation of the polynomial has been suitably changed.
!    Note that XV does not have to be distinct from any of the original XD
!    values.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) i
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xv
  real ( kind = 8 ) yd(nd)
!
!  Recompute the divided difference coefficients.
!
  do i = nd - 1, 1, -1
    yd(i) = yd(i) + ( xv - xd(i) ) * yd(i+1)
  end do
!
!  Shift the XD values up one position, and insert XV at the beginning.
!
  xd(2:nd) = xd(1:nd-1)

  xd(1) = xv

  return
end
subroutine dif_shift_zero ( nd, xd, yd )

!*****************************************************************************80
!
!! DIF_SHIFT_ZERO shifts a divided difference table so all abscissas are zero.
!
!  Discussion:
!
!    When the abscissas are changed, the coefficients naturally
!    must also be changed.
!
!    The resulting pair (XD, YD) still represents the
!    same polynomial, but the entries in YD are now the
!    standard polynomial coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the length of the XD and YD arrays.
!
!    Input/output, real ( kind = 8 ) XD(ND), the X values that 
!    correspond to the divided difference table.  On output, XD 
!    contains only zeroes.
!
!    Input/output, real ( kind = 8 ) YD(ND), the divided difference table
!    for the polynomial.  On output, YD is also the coefficient array for 
!    the standard representation of the polynomial.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) i
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xv
  real ( kind = 8 ) yd(nd)

  xv = 0.0D+00

  do i = 1, nd
    call dif_shift_x ( nd, xd, yd, xv )
  end do

  return
end
subroutine dif_to_r8poly ( nd, xd, yd, c )

!*****************************************************************************80
!
!! DIF_TO_R8POLY converts a divided difference table to a standard polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of coefficients, and abscissas.
!
!    Input, real ( kind = 8 ) XD(ND), the X values used in the divided
!    difference representation of the polynomial.
!
!    Input, real ( kind = 8 ) YD(ND) the divided difference table.
!
!    Output, real ( kind = 8 ) C(ND), the polyomial coefficients.
!    C(1) is the constant term, and C(ND) is the coefficient of X^(ND-1).
!
  implicit none

  integer ( kind = 4 ) nd

  real ( kind = 8 ) c(nd)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)

  c(1:nd) = yd(1:nd)
!
!  Recompute the divided difference coefficients.
!
  do j = 1, nd - 1
    do i = 1, nd - j
      c(nd-i) = c(nd-i) - xd(nd-i-j+1) * c(nd-i+1)
    end do
  end do

  return
end
subroutine fejer2_compute ( order, x, w )

!*****************************************************************************80
!
!! FEJER2_COMPUTE computes a Fejer Type 2 quadrature rule.
!
!  Discussion:
!
!    Our convention is that the points are numbered from left to right.
!
!    The rule is defined on [-1,+1].
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEJER2_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if

  if ( order == 1 ) then
    x(1) = 0.0D+00
    w(1) = 2.0D+00
    return
  end if

  do i = 1, order
    x(i) = cos ( real ( order + 1 - i, kind = 8 ) * pi &
               / real ( order + 1,        kind = 8 ) )
  end do

  if ( mod ( order, 2 ) == 1 ) then
    x((order+1)/2) = 0.0D+00
  end if

  if ( order == 2 ) then

    w(1:2) = 1.0D+00

  else

    do i = 1, order

      theta = real ( order + 1 - i, kind = 8 ) * pi &
            / real ( order + 1, kind = 8 )

      w(i) = 1.0D+00

      do j = 1, ( ( order - 1 ) / 2 )
        w(i) = w(i) - 2.0D+00 &
          * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta ) &
          / real ( 4 * j * j - 1, kind = 8 )
      end do

      p = 2.0D+00 * real ( ( ( order + 1 ) / 2 ), kind = 8 ) - 1.0D+00
      w(i) = w(i) - cos ( ( p + 1.0D+00 ) * theta ) / p

    end do

    w(1:order) = 2.0D+00 * w(1:order) / real ( order + 1, kind = 8 )

  end if

  return
end
subroutine fejer2_compute_points ( order, points )

!*****************************************************************************80
!
!! FEJER2_COMPUTE_POINTS returns the abscissas of a Fejer type 2 rule.
!
!  Discussion:
!
!    Our convention is that the points are numbered from left to right.
!
!    The rule is defined on [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) indx
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) points(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEJER2_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if

  if ( order == 1 ) then

    points(1) = 0.0D+00

  else

    do indx = 1, order
      points(indx) = cos ( real ( order + 1 - indx, kind = 8 ) * pi &
                         / real ( order + 1,        kind = 8 ) )
    end do

    if ( mod ( order, 2 ) == 1 ) then
      points((order+1)/2) = 0.0D+00
    end if

  end if

  return
end
subroutine fejer2_compute_points_np ( order, np, p, points )

!*****************************************************************************80
!
!! FEJER2_COMPUTE_POINTS_NP returns the abscissas of a Fejer type 2 rule.
!
!  Discussion:
!
!    Our convention is that the points are numbered from left to right.
!
!    The rule is defined on [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) points(order)

  call fejer2_compute_points ( order, points )

  return
end
subroutine fejer2_compute_weights ( order, w )

!*****************************************************************************80
!
!! FEJER2_COMPUTE_WEIGHTS computes weights for a Fejer type 2 quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) w(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEJER2_COMPUTE_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if

  if ( order == 1 ) then

    w(1) = 2.0D+00

  else if ( order == 2 ) then

    w(1:2) = 1.0D+00

  else

    do i = 1, order

      theta = real ( order + 1 - i, kind = 8 ) * pi &
            / real ( order + 1,     kind = 8 )

      w(i) = 1.0D+00

      do j = 1, ( ( order - 1 ) / 2 )
        w(i) = w(i) - 2.0D+00 &
          * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta ) &
          / real ( 4 * j * j - 1, kind = 8 )
      end do

      p = 2.0D+00 * real ( ( ( order + 1 ) / 2 ), kind = 8 ) - 1.0D+00
      w(i) = w(i) - cos ( ( p + 1.0D+00 ) * theta ) / p

    end do

    w(1:order) = 2.0D+00 * w(1:order) / real ( order + 1, kind = 8 )

  end if

  return
end
subroutine fejer2_compute_weights_np ( order, np, p, w )

!*****************************************************************************80
!
!! FEJER2_COMPUTE_WEIGHTS_NP: weights for a Fejer type 2 quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) w(order)

  call fejer2_compute_weights ( order, w )

  return
end
subroutine gegenbauer_compute ( order, alpha, x, w )

!*****************************************************************************80
!
!! GEGENBAUER_COMPUTE computes a Gegenbauer quadrature rule.
!
!  Discussion:
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!    Thanks to Janiki Raman for pointing out a problem in an earlier
!    version of the code that occurred when ALPHA was -0.5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X^2) in the weight.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) an
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) cc
  real ( kind = 8 ) delta
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(order)
  real ( kind = 8 ) x0
!
!  Check ORDER.
!
  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEGENBAUER_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  1 <= ORDER is required.'
    stop
  end if
!
!  Check ALPHA.
!
  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEGENBAUER_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < ALPHA is required.'
    stop
  end if
!
!  Set the recursion coefficients.
!
  c(1) = 0.0D+00

  if ( 2 <= order ) then
    c(2) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
  end if

  do i = 3, order

    c(i) = real ( i - 1, kind = 8 ) &
          * ( alpha + alpha + real ( i - 1, kind = 8 ) ) / &
          ( ( alpha + alpha + real ( 2 * i - 1, kind = 8 ) ) &
          * ( alpha + alpha + real ( 2 * i - 3, kind = 8 ) ) )

  end do

  delta = r8_gamma ( alpha         + 1.0D+00 ) &
        * r8_gamma (         alpha + 1.0D+00 ) &
        / r8_gamma ( alpha + alpha + 2.0D+00 )

  cc = delta * 2.0D+00**( 2.0D+00 * alpha + 1.0D+00 ) * product ( c(2:order) )

  do i = 1, order

    if ( i == 1 ) then

      an = alpha / real ( order, kind = 8 )

      r1 = ( 1.0D+00 + alpha ) &
        * ( 2.78D+00 / ( 4.0D+00 + real ( order * order, kind = 8 ) ) &
        + 0.768D+00 * an / real ( order, kind = 8 ) )

      r2 = 1.0D+00 + 2.44D+00 * an + 1.282D+00 * an * an

      x0 = ( r2 - r1 ) / r2

    else if ( i == 2 ) then

      r1 = ( 4.1D+00 + alpha ) / &
        ( ( 1.0D+00 + alpha ) * ( 1.0D+00 + 0.156D+00 * alpha ) )

      r2 = 1.0D+00 + 0.06D+00 * ( real ( order, kind = 8 ) - 8.0D+00 ) * &
        ( 1.0D+00 + 0.12D+00 * alpha ) / real ( order, kind = 8 )

      r3 = 1.0D+00 + 0.012D+00 * alpha * &
        ( 1.0D+00 + 0.25D+00 * abs ( alpha ) ) / real ( order, kind = 8 )

      x0 = x0 - r1 * r2 * r3 * ( 1.0D+00 - x0 )

    else if ( i == 3 ) then

      r1 = ( 1.67D+00 + 0.28D+00 * alpha ) / ( 1.0D+00 + 0.37D+00 * alpha )

      r2 = 1.0D+00 + 0.22D+00 * ( real ( order, kind = 8 ) - 8.0D+00 ) &
        / real ( order, kind = 8 )

      r3 = 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( order * order, kind = 8 ) )

      x0 = x0 - r1 * r2 * r3 * ( x(1) - x0 )

    else if ( i < order - 1 ) then

      x0 = 3.0D+00 * x(i-1) - 3.0D+00 * x(i-2) + x(i-3)

    else if ( i == order - 1 ) then

      r1 = ( 1.0D+00 + 0.235D+00 * alpha ) / ( 0.766D+00 + 0.119D+00 * alpha )

      r2 = 1.0D+00 / ( 1.0D+00 + 0.639D+00 &
        * ( real ( order, kind = 8 ) - 4.0D+00 ) &
        / ( 1.0D+00 + 0.71D+00 * ( real ( order, kind = 8 ) - 4.0D+00 ) ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 20.0D+00 * alpha / ( ( 7.5D+00 + alpha ) * &
        real ( order**2, kind = 8 ) ) )

      x0 = x0 + r1 * r2 * r3 * ( x0 - x(i-2) )

    else if ( i == order ) then

      r1 = ( 1.0D+00 + 0.37D+00 * alpha ) / ( 1.67D+00 + 0.28D+00 * alpha )

      r2 = 1.0D+00 / &
        ( 1.0D+00 + 0.22D+00 * ( real ( order, kind = 8 ) - 8.0D+00 ) &
        / real ( order, kind = 8 ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( order * order, kind = 8 ) ) )

      x0 = x0 + r1 * r2 * r3 * ( x0 - x(i-2) )

    end if

    call gegenbauer_root ( x0, order, alpha, dp2, p1, c )

    x(i) = x0
    w(i) = cc / ( dp2 * p1 )

  end do
!
!  Reverse the order of the data.
!
  x(1:order) = x(order:1:-1)
  w(1:order) = w(order:1:-1)

  return
end
subroutine gegenbauer_integral ( expon, alpha, value )

!*****************************************************************************80
!
!! GEGENBAUER_INTEGRAL integrates a monomial with Gegenbauer weight.
!
!  Discussion:
!
!    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x^2)^ALPHA dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X^2)
!    in the weight factor.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) c
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) value
  real ( kind = 8 ) value1

  if ( mod ( expon, 2 ) == 1 ) then
    value = 0.0D+00
    return
  end if

  c = real ( expon, kind = 8 )

  arg1 = - alpha
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + alpha + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

  value = 2.0D+00 * r8_gamma ( 1.0D+00 + c ) * r8_gamma ( 1.0D+00 + alpha ) &
    * value1 / r8_gamma ( 2.0D+00 + alpha  + c )

  return
end
subroutine gegenbauer_recur ( p2, dp2, p1, x, order, alpha, c )

!*****************************************************************************80
!
!! GEGENBAUER_RECUR finds the value and derivative of a Gegenbauer polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of J(ORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of J'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X^2).
!
!    Input, real ( kind = 8 ) C(ORDER), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = x * p1 - c(i) * p0
    dp2 = x * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine gegenbauer_root ( x, order, alpha, dp2, p1, c )

!*****************************************************************************80
!
!! GEGENBAUER_ROOT improves an approximate root of a Gegenbauer polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X^2).
!
!    Output, real ( kind = 8 ) DP2, the value of J'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) C(ORDER), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = r8_epsilon ( )

  do step = 1, step_max

    call gegenbauer_recur ( p2, dp2, p1, x, order, alpha, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine gen_hermite_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The code uses an algorithm by Elhay and Kautsky.
!
!    The integral:
!
!      integral ( -oo < x < +oo ) |x|^alpha exp(-x^2) f(x) dx
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
!    30 April 2011
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
!    Input, integer ( kind = 4 ) N, the number of abscissas.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
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
  zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    if ( mod ( i, 2 ) == 1 ) then
      bj(i) = ( i_r8 + alpha ) / 2.0D+00
    else
      bj(i) = i_r8 / 2.0D+00
    end if
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine gen_hermite_compute_points ( order, alpha, points )

!*****************************************************************************80
!
!! GEN_HERMITE_COMPUTE_POINTS: abscissas of a Generalized Hermite rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  call gen_hermite_compute ( order, alpha, points, weight )

  return
end
subroutine gen_hermite_compute_points_np ( order, np, p, points )

!*****************************************************************************80
!
!! GEN_HERMITE_COMPUTE_POINTS_NP: abscissas of a Generalized Hermite rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) p(*)
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  alpha = p(1)
  call gen_hermite_compute ( order, alpha, points, weight )

  return
end
subroutine gen_hermite_compute_weights ( order, alpha, weight )

!*****************************************************************************80
!
!! GEN_HERMITE_COMPUTE_WEIGHTS: weights of a Generalized Hermite rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    Set ALPHA = 0.0 for the simplest rule.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  call gen_hermite_compute ( order, alpha, points, weight )

  return
end
subroutine gen_hermite_compute_weights_np ( order, np, p, weight )

!*****************************************************************************80
!
!! GEN_HERMITE_COMPUTE_WEIGHTS_NP: weights of a Generalized Hermite rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) p(*)
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  alpha = p(1)
  call gen_hermite_compute ( order, alpha, points, weight )

  return
end
subroutine gen_hermite_dr_compute ( order, alpha, x, w )

!*****************************************************************************80
!
!! GEN_HERMITE_DR_COMPUTE computes a Generalized Hermite rule.
!
!  Discussion:
!
!    The integral to be approximated has the form:
!
!      Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!
!    Output, real ( kind = 8 ) X(ORDER), W(ORDER), the abscissas and weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_laguerre
  real ( kind = 8 ) arg
  integer ( kind = 4 ) order_laguerre
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(order)
  real ( kind = 8 ), allocatable, dimension ( : ) :: w_laguerre
  real ( kind = 8 ) x(order)
  real ( kind = 8 ), allocatable, dimension ( : ) :: x_laguerre

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEN_HERMITE_DR_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if
!
!  Generate the related Generalized Laguerre rule.
!
  if ( order == 1 ) then
    arg = ( alpha + 1.0D+00 ) / 2.0D+00
    x(1) = 0.0D+00
    w(1) = r8_gamma ( arg )
    return
  end if

  if ( mod ( order, 2 ) == 0 ) then
    order_laguerre = order / 2
    alpha_laguerre = ( alpha - 1.0D+00 ) / 2.0D+00
  else if ( mod ( order, 2 ) == 1 ) then
    order_laguerre = ( order - 1 ) / 2
    alpha_laguerre = ( alpha + 1.0D+00 ) / 2.0D+00
  end if

  allocate ( w_laguerre(order_laguerre) )
  allocate ( x_laguerre(order_laguerre) )

  call gen_laguerre_ss_compute ( order_laguerre, alpha_laguerre, x_laguerre, &
    w_laguerre )

  if ( mod ( order, 2 ) == 0 ) then

    x(1:order_laguerre) = - sqrt ( x_laguerre(order_laguerre:1:-1) )
    x(order_laguerre+1:order_laguerre+order_laguerre) &
      = sqrt ( x_laguerre(1:order_laguerre) )

    w(1:order_laguerre) = 0.5D+00 * w_laguerre(order_laguerre:1:-1)
    w(order_laguerre+1:order_laguerre+order_laguerre) &
      = 0.5D+00 * w_laguerre(1:order_laguerre)

  else if ( mod ( order, 2 ) == 1 ) then

    x(1:order_laguerre) = - sqrt ( x_laguerre(order_laguerre:1:-1) )
    x(order_laguerre+1) = 0.0D+00
    x(order_laguerre+2:order_laguerre+order_laguerre+1) &
      = sqrt ( x_laguerre(1:order_laguerre) )

    w(1:order_laguerre) = 0.5D+00 * w_laguerre(order_laguerre:1:-1) &
                                  / x_laguerre(order_laguerre:1:-1)

    arg = ( alpha + 1.0D+00 ) / 2.0D+00
    w(order_laguerre+1) &
      = r8_gamma ( arg ) &
      - sum ( w_laguerre(1:order_laguerre) / x_laguerre(1:order_laguerre) )

    w(order_laguerre+2:order_laguerre+order_laguerre+1) &
      = 0.5D+00 * w_laguerre(1:order_laguerre) &
                / x_laguerre(1:order_laguerre)

  end if

  deallocate ( w_laguerre )
  deallocate ( x_laguerre )

  return
end
subroutine gen_hermite_integral ( expon, alpha, value )

!*****************************************************************************80
!
!! GEN_HERMITE_INTEGRAL evaluates a monomial Generalized Hermite integral.
!
!  Discussion:
!
!    H(n,alpha) = Integral ( -oo < x < +oo ) x^n |x|^alpha exp(-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent of the monomial.
!    0 <= EXPON.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of |X| in the
!    weight function.  -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) value

  if ( mod ( expon, 2 ) == 1 ) then

    value = 0.0D+00

  else

    a = alpha + real ( expon, kind = 8 )

    if ( a <= -1.0D+00 ) then

      value = - huge ( value )

    else

      value = r8_gamma ( ( a + 1.0D+00 ) / 2.0D+00 )

    end if

  end if

  return
end
subroutine gen_laguerre_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! GEN_LAGUERRE_COMPUTE: generalized Gauss-Laguerre quadrature rule.
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
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * exp ( x(i) ) * f ( x(i) )
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
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = gamma ( alpha + 1.0D+00 )
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
subroutine gen_laguerre_compute_points ( order, alpha, points )

!*****************************************************************************80
!
!! GEN_LAGUERRE_COMPUTE_POINTS: points of a Generalized Laguerre rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    Set ALPHA = 0.0 for the simplest rule.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  call gen_laguerre_compute ( order, alpha, points, weight )

  return
end
subroutine gen_laguerre_compute_points_np ( order, np, p, points )

!*****************************************************************************80
!
!! GEN_LAGUERRE_COMPUTE_POINTS_NP: points of a Generalized Laguerre rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) p(*)
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  alpha = p(1)

  call gen_laguerre_compute ( order, alpha, points, weight )

  return
end
subroutine gen_laguerre_compute_weights ( order, alpha, weight )

!*****************************************************************************80
!
!! GEN_LAGUERRE_COMPUTE_WEIGHTS: weights of a Generalized Laguerre rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    Set ALPHA = 0.0 for the simplest rule.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  call gen_laguerre_compute ( order, alpha, points, weight )

  return
end
subroutine gen_laguerre_compute_weights_np ( order, np, p, weight )

!*****************************************************************************80
!
!! GEN_LAGUERRE_COMPUTE_WEIGHTS_NP: weights of a Generalized Laguerre rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) p(*)
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  alpha = p(1)

  call gen_laguerre_compute ( order, alpha, points, weight )

  return
end
subroutine gen_laguerre_integral ( expon, alpha, exact )

!*****************************************************************************80
!
!! GEN_LAGUERRE_INTEGRAL evaluates a monomial genearlized Laguerre integral.
!
!  Discussion:
!
!    The integral being computed is
!
!      integral ( 0 <= x < +oo ) x^n * x^alpha * exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma

  arg = alpha + real ( expon + 1, kind = 8 )

  exact = r8_gamma ( arg )

  return
end
subroutine gen_laguerre_ss_compute ( order, alpha, x, w )

!*****************************************************************************80
!
!! GEN_LAGUERRE_SS_COMPUTE computes a Generalized Laguerre quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ 0, +oo ).
!
!    The weight function is w(x) = exp ( -x ) * x**alpha.
!
!
!    If the integral to approximate is:
!
!        Integral ( 0 <= X < +oo ) exp ( - X ) * X^ALPHA * F(X) dX
!
!    then the quadrature rule is:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!
!    If the integral to approximate is:
!
!        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX
!
!    then the quadrature rule is:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * exp ( X(I) ) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    Original FORTRAN77 by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) ratio
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(order)
  real ( kind = 8 ) x0

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEN_LAGUERRE_SS_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if
!
!  Set the recursion coefficients.
!
  do i = 1, order
    b(i) = ( alpha + real ( 2 * i - 1, kind = 8 ) )
  end do

  do i = 1, order
    c(i) = real ( i - 1, kind = 8 ) * ( alpha + real ( i - 1, kind = 8 ) )
  end do

  cc = r8_gamma ( alpha + 1.0D+00 ) * product ( c(2:order) )

  do i = 1, order
!
!  Compute an estimate for the root.
!
    if ( i == 1 ) then

      x0 = ( 1.0D+00 + alpha ) * ( 3.0D+00 + 0.92 * alpha ) / &
        ( 1.0D+00 + 2.4D+00 * real ( order, kind = 8 ) + 1.8D+00 * alpha )

    else if ( i == 2 ) then

      x0 = x0 + ( 15.0D+00 + 6.25D+00 * alpha ) / &
        ( 1.0D+00 + 0.9D+00 * alpha + 2.5D+00 * real ( order, kind = 8 ) )

    else

      r1 = ( 1.0D+00 + 2.55D+00 * real ( i - 2, kind = 8 ) ) &
        / ( 1.9D+00 * real ( i - 2, kind = 8 ) )

      r2 = 1.26D+00 * real ( i - 2, kind = 8 ) * alpha / &
        ( 1.0D+00 + 3.5D+00 * real ( i - 2, kind = 8 ) )

      ratio = ( r1 + r2 ) / ( 1.0D+00 + 0.3D+00 * alpha )

      x0 = x0 + ratio * ( x0 - x(i-2) )

    end if
!
!  Use iteration to find the root.
!
    call gen_laguerre_ss_root ( x0, order, alpha, dp2, p1, b, c )
!
!  Set the abscissa and weight.
!
    x(i) = x0
    w(i) = ( cc / dp2 ) / p1

  end do

  return
end
subroutine gen_laguerre_ss_recur ( p2, dp2, p1, x, order, alpha, b, c )

!*****************************************************************************80
!
!! GEN_LAGUERRE_SS_RECUR evaluates a Generalized Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of L(ORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of L'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor in the
!    integrand.
!
!    Input, real ( kind = 8 ) B(ORDER), C(ORDER), the recursion
!    coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x - alpha - 1.0D+00
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine gen_laguerre_ss_root ( x, order, alpha, dp2, p1, b, c )

!*****************************************************************************80
!
!! GEN_LAGUERRE_SS_ROOT seeks roots of a Generalized Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!
!    Output, real ( kind = 8 ) DP2, the value of L'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) B(ORDER), C(ORDER), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = r8_epsilon ( )

  do step = 1, step_max

    call gen_laguerre_ss_recur ( p2, dp2, p1, x, order, alpha, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine hc_compute_weights_from_points ( m, x, w )

!*****************************************************************************80
!
!! HC_COMPUTE_WEIGHTS_FROM_POINTS: Hermite-Cubic weights, user-supplied points.
!
!  Discussion:
!
!    An interval [A,B] has been divided by NHALF points X; at each
!    point both function and derivative information is available.
!
!    The piecewise cubic Hermite interpolant is constructed for this data.
!
!    A quadrature rule is determined for the interpolant.
!
!    There will be N=2*M weights.  If the quadrature rule is to be written 
!    out, one would normally list each point twice, so that the number of points
!    and weights are equal.  The listing of the same point value twice is an
!    implicit indication that both function and derivative values should be
!    used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of points, not counting 
!    repetitions.
!
!    Input, real ( kind = 8 ) X(M), the points, without repetition.
!
!    Output, real ( kind = 8 ) W(2*M), the weights.  The first two weights are 
!    associated with the first point, and so on.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) j
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m)

  w(1) = 0.5D+00 * ( x(2) - x(1) )
  w(2) = ( x(2) - x(1) )**2 / 12.0D+00

  do j = 2, m - 1
    w(1+(j-1)*2) = 0.5D+00 * ( x(j+1) - x(j-1) )
    w(2+(j-1)*2) = ( x(j+1) - x(j-1) ) &
      * ( x(j+1) - 2.0D+00 * x(j) + x(j-1) ) / 12.0D+00
  end do

  w(1+(m-1)*2) = 0.5D+00 * ( x(m) - x(m-1) )
  w(2+(m-1)*2) = - ( x(m-1) - x(m) )**2 / 12.0D+00

  return
end
subroutine hcc_compute ( n, x, w )

!*****************************************************************************80
!
!! HCC_COMPUTE computes a Hermite-Cubic-Chebyshev-Spacing quadrature rule.
!
!  Discussion:
!
!    For the HCC rule, we assume that an interval has been divided by
!    M nodes X into Chebyshev-spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hcc_compute_points ( n, x )
  call hcc_compute_weights ( n, w )

  return
end
subroutine hcc_compute_np ( n, np, p, x, w )

!*****************************************************************************80
!
!! HCC_COMPUTE_NP computes a Hermite-Cubic-Chebyshev-Spacing quadrature rule.
!
!  Discussion:
!
!    For the HCC rule, we assume that an interval has been divided by
!    M nodes X into Chebyshev-spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hcc_compute ( n, x, w )

  return
end
subroutine hcc_compute_points ( n, x )

!*****************************************************************************80
!
!! HCC_COMPUTE_POINTS: abscissas of a Hermite-Cubic-Chebyshev-Spacing rule.
!
!  Discussion:
!
!    For the HCC rule, we assume that an interval has been divided by
!    M nodes X into Chebyshev-spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_value

  if ( mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HCC_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Order of rule N is not even.'
    stop
  end if

  m = n / 2

  if ( m == 1 ) then
    x(1:2) = 0.0D+00
    return
  end if

  j = 0
  do i = 1, m

    if ( i == 1 ) then
      x_value = - 1.0D+00
    else if ( 2 * i - 1 == m ) then
      x_value = 0.0D+00
    else if ( i == m ) then
      x_value = 0.0D+00
    else
      x_value = cos ( real ( m - i, kind = 8 ) * pi &
                    / real ( m - 1, kind = 8 ) )
    end if

    j = j + 1
    x(j)   = x_value
    j = j + 1
    x(j) = x_value

  end do

  return
end
subroutine hcc_compute_points_np ( n, np, p, x )

!*****************************************************************************80
!
!! HCC_COMPUTE_POINTS_NP: abscissas of a Hermite-Cubic-Chebyshev-Spacing rule.
!
!  Discussion:
!
!    For the HCC rule, we assume that an interval has been divided by
!    M nodes X into Chebyshev-spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) x(n)

  call hcc_compute_points ( n, x )

  return
end
subroutine hcc_compute_weights ( n, w )

!*****************************************************************************80
!
!! HCC_COMPUTE_WEIGHTS computes Hermite-Cubic-Chebyshev-Spacing weights.
!
!  Discussion:
!
!    For the HCC rule, we assume that an interval has been divided by
!    M nodes X into Chebyshev-spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) w(n)
  real ( kind = 8 ), allocatable :: x(:)

  if ( mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HCC_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Order of rule N is not even.'
    stop
  end if

  m = n / 2

  allocate ( x(1:m) )

  call clenshaw_curtis_compute_points ( m, x )

  call hc_compute_weights_from_points ( m, x, w )

  deallocate ( x )

  return
end
subroutine hcc_compute_weights_np ( n, np, p, w )

!*****************************************************************************80
!
!! HCC_COMPUTE_WEIGHTS_NP computes Hermite-Cubic-Chebyshev-Spacing weights.
!
!  Discussion:
!
!    For the HCC rule, we assume that an interval has been divided by
!    M nodes X into Chebyshev-spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) w(n)

  call hcc_compute_weights ( n, w )

  return
end
subroutine hce_compute ( n, x, w )

!*****************************************************************************80
!
!! HCE_COMPUTE computes a Hermite-Cubic-Equal-Spacing quadrature rule.
!
!  Discussion:
!
!    For the HCE rule, we assume that an interval has been divided by
!    M nodes X into equally spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hce_compute_points ( n, x )
  call hce_compute_weights ( n, w )

  return
end
subroutine hce_compute_np ( n, np, p, x, w )

!*****************************************************************************80
!
!! HCE_COMPUTE_NP computes a Hermite-Cubic-Equal-Spacing quadrature rule.
!
!  Discussion:
!
!    For the HCE rule, we assume that an interval has been divided by
!    M nodes X into equally spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hce_compute ( n, x, w )

  return
end
subroutine hce_compute_points ( n, x )

!*****************************************************************************80
!
!! HCE_COMPUTE_POINTS: abscissas of a Hermite-Cubic-Equal-Spacing rule.
!
!  Discussion:
!
!    For the HCE rule, we assume that an interval has been divided by
!    M nodes X into equally spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_value

  if ( mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HCE_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Order of rule N is not even.'
    stop
  end if

  m = n / 2

  do j = 1, m

    x_value = real ( 2 * j - 1 - m, kind = 8 ) / real ( m - 1, kind = 8  )

    do i = 1, 2
      x(i+(j-1)*2) = x_value
    end do

  end do

  return
end
subroutine hce_compute_points_np ( n, np, p, x )

!*****************************************************************************80
!
!! HCE_COMPUTE_POINTS_NP: abscissas of a Hermite-Cubic-Equal-Spacing rule.
!
!  Discussion:
!
!    For the HCE rule, we assume that an interval has been divided by
!    M nodes X into equally spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) x(n)

  call hce_compute_points ( n, x )

  return
end
subroutine hce_compute_weights ( n, w )

!*****************************************************************************80
!
!! HCE_COMPUTE_WEIGHTS computes Hermite-Cubic-Equal-Spacing weights.
!
!  Discussion:
!
!    For the HCE rule, we assume that an interval has been divided by
!    M nodes X into equally spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) w(n)
  real ( kind = 8 ), allocatable :: x(:)

  if ( mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HCE_COMPUTE_WEIGHTS - Fatal error!'
    write ( *, '(a)' ) '  Order of rule N is not even.'
    stop
  end if

  m = n / 2

  allocate ( x(1:m) )

  do j = 1, m
    x(j) = real ( 2 * j - 1 - m, kind = 8 ) / real ( m - 1, kind = 8  )
  end do

  call hc_compute_weights_from_points ( m, x, w )

  deallocate ( x )

  return
end
subroutine hce_compute_weights_np ( n, np, p, w )

!*****************************************************************************80
!
!! HCE_COMPUTE_WEIGHTS_NP computes Hermite-Cubic-Equal-Spacing weights.
!
!  Discussion:
!
!    For the HCE rule, we assume that an interval has been divided by
!    M nodes X into equally spaced subintervals, and that at each
!    abscissa both function and derivative information is available.
!    The piecewise cubic Hermite interpolant is constructed for this data.
!    The quadrature rule uses N = 2 * M abscissas, where each node is
!    listed twice, and the weights occur in pairs, with the first multiplying
!    the function value and the second the derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(np)
  real ( kind = 8 ) w(n)

  call hce_compute_weights ( n, w )

  return
end
subroutine hermite_compute ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The code uses an algorithm by Elhay and Kautsky.
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!    The integral:
!
!      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
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
!    09 May 2012
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
!    Input, integer ( kind = 4 ) N, the number of abscissas.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_gamma ( 1.0D+00 / 2.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = 8 ) / 2.0D+00
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )
!
!  If N is odd, force the middle X to be zero.
!
  if ( mod ( n, 2 ) == 1 ) then
    x((n+1)/2) = 0.0D+00
  end if

  w(1:n) = w(1:n)**2

  return
end
subroutine hermite_compute_points ( n, x )

!*****************************************************************************80
!
!! HERMITE_COMPUTE_POINTS computes points of a Hermite quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hermite_compute ( n, x, w )

  return
end
subroutine hermite_compute_points_np ( n, np, p, x )

!*****************************************************************************80
!
!! HERMITE_COMPUTE_POINTS_NP computes points of a Hermite quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hermite_compute ( n, x, w )

  return
end
subroutine hermite_compute_weights ( n, w )

!*****************************************************************************80
!
!! HERMITE_COMPUTE_WEIGHTS computes weights of a Hermite quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hermite_compute ( n, x, w )

  return
end
subroutine hermite_compute_weights_np ( n, np, p, w )

!*****************************************************************************80
!
!! HERMITE_COMPUTE_WEIGHTS_NP computes weights of a Hermite quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hermite_compute ( n, x, w )

  return
end
subroutine hermite_genz_keister_lookup ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_GENZ_KEISTER_LOOKUP returns a Genz-Keister rule for Hermite problems.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+?, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and a final rule of order
!    35, 37, 41 or 43.
!
!    The precisions of these rules are P = 1, 5, 15, 29, 
!    with the final rule of precision 51, 55, 63 or 67.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, 35, 37, 41 and 43.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call hermite_genz_keister_lookup_points ( n, x )
  call hermite_genz_keister_lookup_weights ( n, w );

  return
end
subroutine hermite_genz_keister_lookup_points ( n, x )

!*****************************************************************************80
!
!! HERMITE_GENZ_KEISTER_LOOKUP_POINTS: abscissas of a Genz-Keister Hermite rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+?, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and a final rule of order
!    35, 37, 41 or 43.
!
!    The precisions of these rules are P = 1, 5, 15, 29, 
!    with the final rule of precision 51, 55, 63 or 67.
!
!    Three related families begin the same way, but end with a different final
!    rule.  As a convenience, this function includes these final rules as well:
!
!    Designation  Orders       Precisions
!
!    1+2+6+10+16  1,3,9,19,35  1,5,15,29,51
!    1+2+6+10+18  1,3,9,19,37  1,5,15,29,55
!    1+2+6+10+22  1,3,9,19,41  1,5,15,29,63
!    1+2+6+10+24  1,3,9,19,43  1,5,15,29,67
!
!    Some of the data in this function was kindly supplied directly by
!    Alan Genz on 24 April 2011.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, 35, 37, 41, or 43.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x( 1) =   0.0000000000000000D+00

  else if ( n == 3 ) then

    x( 1) =  -1.2247448713915889D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   1.2247448713915889D+00

  else if ( n == 9 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -2.0232301911005157D+00
    x( 3) =  -1.2247448713915889D+00
    x( 4) =  -5.2403354748695763D-01
    x( 5) =   0.0000000000000000D+00
    x( 6) =   5.2403354748695763D-01
    x( 7) =   1.2247448713915889D+00
    x( 8) =   2.0232301911005157D+00
    x( 9) =   2.9592107790638380D+00

  else if ( n == 19 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.2665132620567876D+00
    x( 5) =  -2.0232301911005157D+00
    x( 6) =  -1.8357079751751868D+00
    x( 7) =  -1.2247448713915889D+00
    x( 8) =  -8.7004089535290285D-01
    x( 9) =  -5.2403354748695763D-01
    x(10) =   0.0000000000000000D+00
    x(11) =   5.2403354748695763D-01
    x(12) =   8.7004089535290285D-01
    x(13) =   1.2247448713915889D+00
    x(14) =   1.8357079751751868D+00
    x(15) =   2.0232301911005157D+00
    x(16) =   2.2665132620567876D+00
    x(17) =   2.9592107790638380D+00
    x(18) =   3.6677742159463378D+00
    x(19) =   4.4995993983103881D+00

  else if ( n == 35 ) then

    x( 1) =  -6.3759392709822356D+00
    x( 2) =  -5.6432578578857449D+00
    x( 3) =  -5.0360899444730940D+00
    x( 4) =  -4.4995993983103881D+00
    x( 5) =  -4.0292201405043713D+00
    x( 6) =  -3.6677742159463378D+00
    x( 7) =  -3.3491639537131945D+00
    x( 8) =  -2.9592107790638380D+00
    x( 9) =  -2.5705583765842968D+00
    x(10) =  -2.2665132620567876D+00
    x(11) =  -2.0232301911005157D+00
    x(12) =  -1.8357079751751868D+00
    x(13) =  -1.5794121348467671D+00
    x(14) =  -1.2247448713915889D+00
    x(15) =  -8.7004089535290285D-01
    x(16) =  -5.2403354748695763D-01
    x(17) =  -1.7606414208200893D-01
    x(18) =   0.0000000000000000D+00
    x(19) =   1.7606414208200893D-01
    x(20) =   5.2403354748695763D-01
    x(21) =   8.7004089535290285D-01
    x(22) =   1.2247448713915889D+00
    x(23) =   1.5794121348467671D+00
    x(24) =   1.8357079751751868D+00
    x(25) =   2.0232301911005157D+00
    x(26) =   2.2665132620567876D+00
    x(27) =   2.5705583765842968D+00
    x(28) =   2.9592107790638380D+00
    x(29) =   3.3491639537131945D+00
    x(30) =   3.6677742159463378D+00
    x(31) =   4.0292201405043713D+00
    x(32) =   4.4995993983103881D+00
    x(33) =   5.0360899444730940D+00
    x(34) =   5.6432578578857449D+00
    x(35) =   6.3759392709822356D+00

  else if ( n == 37 ) then

    x( 1) =  -6.853200069757519D+00
    x( 2) =  -6.124527854622158D+00
    x( 3) =  -5.521865209868350D+00
    x( 4) =  -4.986551454150765D+00
    x( 5) =  -4.499599398310388D+00
    x( 6) =  -4.057956316089741D+00
    x( 7) =  -3.667774215946338D+00
    x( 8) =  -3.315584617593290D+00
    x( 9) =  -2.959210779063838D+00
    x(10) =  -2.597288631188366D+00
    x(11) =  -2.266513262056788D+00
    x(12) =  -2.023230191100516D+00
    x(13) =  -1.835707975175187D+00
    x(14) =  -1.561553427651873D+00
    x(15) =  -1.224744871391589D+00
    x(16) =  -0.870040895352903D+00
    x(17) =  -0.524033547486958D+00
    x(18) =  -0.214618180588171D+00
    x(19) =   0.000000000000000D+00
    x(20) =   0.214618180588171D+00
    x(21) =   0.524033547486958D+00
    x(22) =   0.870040895352903D+00
    x(23) =   1.224744871391589D+00
    x(24) =   1.561553427651873D+00
    x(25) =   1.835707975175187D+00
    x(26) =   2.023230191100516D+00
    x(27) =   2.266513262056788D+00
    x(28) =   2.597288631188366D+00
    x(29) =   2.959210779063838D+00
    x(30) =   3.315584617593290D+00
    x(31) =   3.667774215946338D+00
    x(32) =   4.057956316089741D+00
    x(33) =   4.499599398310388D+00
    x(34) =   4.986551454150765D+00
    x(35) =   5.521865209868350D+00
    x(36) =   6.124527854622158D+00
    x(37) =   6.853200069757519D+00

  else if ( n == 41 ) then

    x( 1) =  -7.251792998192644D+00
    x( 2) =  -6.547083258397540D+00
    x( 3) =  -5.961461043404500D+00
    x( 4) =  -5.437443360177798D+00
    x( 5) =  -4.953574342912980D+00
    x( 6) =  -4.4995993983103881D+00
    x( 7) =  -4.070919267883068D+00
    x( 8) =  -3.6677742159463378D+00
    x( 9) =  -3.296114596212218D+00
    x(10) =  -2.9592107790638380D+00
    x(11) =  -2.630415236459871D+00
    x(12) =  -2.2665132620567876D+00
    x(13) =  -2.043834754429505D+00
    x(14) =  -2.0232301911005157D+00
    x(15) =  -1.8357079751751868D+00
    x(16) =  -1.585873011819188D+00
    x(17) =  -1.2247448713915889D+00
    x(18) =  -0.87004089535290285D+00
    x(19) =  -0.52403354748695763D+00
    x(20) =  -0.195324784415805D+00
    x(21) =   0.0000000000000000D+00
    x(22) =   0.195324784415805D+00
    x(23) =   0.52403354748695763D+00
    x(24) =   0.87004089535290285D+00
    x(25) =   1.2247448713915889D+00
    x(26) =   1.585873011819188D+00
    x(27) =   1.8357079751751868D+00
    x(28) =   2.0232301911005157D+00
    x(29) =   2.043834754429505D+00
    x(30) =   2.2665132620567876D+00
    x(31) =   2.630415236459871D+00
    x(32) =   2.9592107790638380D+00
    x(33) =   3.296114596212218D+00
    x(34) =   3.6677742159463378D+00
    x(35) =   4.070919267883068D+00
    x(36) =   4.4995993983103881D+00
    x(37) =   4.953574342912980D+00
    x(38) =   5.437443360177798D+00
    x(39) =   5.961461043404500D+00
    x(40) =   6.547083258397540D+00
    x(41) =   7.251792998192644D+00

  else if ( n == 43 ) then

    x( 1) = -10.167574994881873D+00
    x( 2) =  -7.231746029072501D+00
    x( 3) =  -6.535398426382995D+00
    x( 4) =  -5.954781975039809D+00
    x( 5) =  -5.434053000365068D+00
    x( 6) =  -4.952329763008589D+00
    x( 7) =  -4.4995993983103881D+00
    x( 8) =  -4.071335874253583D+00
    x( 9) =  -3.6677742159463378D+00
    x(10) =  -3.295265921534226D+00
    x(11) =  -2.9592107790638380D+00
    x(12) =  -2.633356763661946D+00
    x(13) =  -2.2665132620567876D+00
    x(14) =  -2.089340389294661D+00
    x(15) =  -2.0232301911005157D+00
    x(16) =  -1.8357079751751868D+00
    x(17) =  -1.583643465293944D+00
    x(18) =  -1.2247448713915889D+00
    x(19) =  -0.87004089535290285D+00
    x(20) =  -0.52403354748695763D+00
    x(21) =  -0.196029453662011D+00
    x(22) =   0.0000000000000000D+00
    x(23) =   0.196029453662011D+00
    x(24) =   0.52403354748695763D+00
    x(25) =   0.87004089535290285D+00
    x(26) =   1.2247448713915889D+00
    x(27) =   1.583643465293944D+00
    x(28) =   1.8357079751751868D+00
    x(29) =   2.0232301911005157D+00
    x(30) =   2.089340389294661D+00
    x(31) =   2.2665132620567876D+00
    x(32) =   2.633356763661946D+00
    x(33) =   2.9592107790638380D+00
    x(34) =   3.295265921534226D+00
    x(35) =   3.6677742159463378D+00
    x(36) =   4.071335874253583D+00
    x(37) =   4.4995993983103881D+00
    x(38) =   4.952329763008589D+00
    x(39) =   5.434053000365068D+00
    x(40) =   5.954781975039809D+00
    x(41) =   6.535398426382995D+00
    x(42) =   7.231746029072501D+00
    x(43) =  10.167574994881873D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GENZ_KEISTER_LOOKUP_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 9, 19, 35, 37, 41, or 43.'
    stop

  end if

  return
end
subroutine hermite_genz_keister_lookup_points_np ( n, np, p, x )

!*****************************************************************************80
!
!! HERMITE_GENZ_KEISTER_LOOKUP_POINTS_NP: Genz-Keister Hermite abscissas.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+?, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and a final rule of order
!    35, 37, 41 or 43.
!
!    The precisions of these rules are P = 1, 5, 15, 29, 
!    with the final rule of precision 51, 55, 63 or 67.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, 35, 37, 41 and 43.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) x(n)

  call hermite_genz_keister_lookup_points ( n, x )

  return
end
subroutine hermite_genz_keister_lookup_weights ( n, w )

!*****************************************************************************80
!
!! HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS: weights for Genz-Keister Hermite rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+?, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and a final rule of order
!    35, 37, 41 or 43.
!
!    The precisions of these rules are P = 1, 5, 15, 29, 
!    with the final rule of precision 51, 55, 63 or 67.
!
!    Three related families begin the same way, but end with a different final
!    rule.  As a convenience, this function includes these final rules as well:
!
!    Designation  Orders       Precisions
!
!    1+2+6+10+16, 1,3,9,19,35  1,5,15,29,51
!    1+2+6+10+18  1,3,9,19,37  1,5,15,29,55
!    1+2+6+10+22  1,3,9,19,41  1,5,15,29,63
!    1+2+6+10+24  1,3,9,19,43  1,5,15,29,67
!
!    Some of the data in this function was kindly supplied directly by
!    Alan Genz on 24 April 2011.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, 35, 37, 41 and 43.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)

  if ( n == 1 ) then

    w( 1) =   1.7724538509055159D+00

  else if ( n == 3 ) then

    w( 1) =   2.9540897515091930D-01
    w( 2) =   1.1816359006036772D+00
    w( 3) =   2.9540897515091930D-01

  else if ( n == 9 ) then

    w( 1) =   1.6708826306882348D-04
    w( 2) =   1.4173117873979098D-02
    w( 3) =   1.6811892894767771D-01
    w( 4) =   4.7869428549114124D-01
    w( 5) =   4.5014700975378197D-01
    w( 6) =   4.7869428549114124D-01
    w( 7) =   1.6811892894767771D-01
    w( 8) =   1.4173117873979098D-02
    w( 9) =   1.6708826306882348D-04

  else if ( n == 19 ) then

    w( 1) =   1.5295717705322357D-09
    w( 2) =   1.0802767206624762D-06
    w( 3) =   1.0656589772852267D-04
    w( 4) =   5.1133174390883855D-03
    w( 5) =  -1.1232438489069229D-02
    w( 6) =   3.2055243099445879D-02
    w( 7) =   1.1360729895748269D-01
    w( 8) =   1.0838861955003017D-01
    w( 9) =   3.6924643368920851D-01
    w(10) =   5.3788160700510168D-01
    w(11) =   3.6924643368920851D-01
    w(12) =   1.0838861955003017D-01
    w(13) =   1.1360729895748269D-01
    w(14) =   3.2055243099445879D-02
    w(15) =  -1.1232438489069229D-02
    w(16) =   5.1133174390883855D-03
    w(17) =   1.0656589772852267D-04
    w(18) =   1.0802767206624762D-06
    w(19) =   1.5295717705322357D-09

  else if ( n == 35 ) then

    w( 1) =   1.8684014894510604D-18
    w( 2) =   9.6599466278563243D-15
    w( 3) =   5.4896836948499462D-12
    w( 4) =   8.1553721816916897D-10
    w( 5) =   3.7920222392319532D-08
    w( 6) =   4.3737818040926989D-07
    w( 7) =   4.8462799737020461D-06
    w( 8) =   6.3328620805617891D-05
    w( 9) =   4.8785399304443770D-04
    w(10) =   1.4515580425155904D-03
    w(11) =   4.0967527720344047D-03
    w(12) =   5.5928828911469180D-03
    w(13) =   2.7780508908535097D-02
    w(14) =   8.0245518147390893D-02
    w(15) =   1.6371221555735804D-01
    w(16) =   2.6244871488784277D-01
    w(17) =   3.3988595585585218D-01
    w(18) =   9.1262675363737921D-04
    w(19) =   3.3988595585585218D-01
    w(20) =   2.6244871488784277D-01
    w(21) =   1.6371221555735804D-01
    w(22) =   8.0245518147390893D-02
    w(23) =   2.7780508908535097D-02
    w(24) =   5.5928828911469180D-03
    w(25) =   4.0967527720344047D-03
    w(26) =   1.4515580425155904D-03
    w(27) =   4.8785399304443770D-04
    w(28) =   6.3328620805617891D-05
    w(29) =   4.8462799737020461D-06
    w(30) =   4.3737818040926989D-07
    w(31) =   3.7920222392319532D-08
    w(32) =   8.1553721816916897D-10
    w(33) =   5.4896836948499462D-12
    w(34) =   9.6599466278563243D-15
    w(35) =   1.8684014894510604D-18

  else if ( n == 37 ) then

    w( 1) = 0.337304188079177058D-20
    w( 2) = 0.332834739632930463D-16
    w( 3) = 0.323016866782871498D-13
    w( 4) = 0.809333688669950037D-11
    w( 5) = 0.748907559239519284D-09
    w( 6) = 0.294146671497083432D-07
    w( 7) = 0.524482423744884136D-06
    w( 8) = 0.586639457073896277D-05
    w( 9) = 0.571885531470621903D-04
    w(10) = 0.41642095727577091D-03
    w(11) = 0.174733389581099482D-02
    w(12) = 0.313373786000304381D-02
    w(13) = 0.768092665770660459D-02
    w(14) = 0.274962713372148476D-01
    w(15) = 0.783630990508037449D-01
    w(16) = 0.16611584261479281D+00
    w(17) = 0.253636910481387185D+00
    w(18) = 0.261712932511430884D+00
    w(19) = 0.171719680968980257D+00
    w(20) = 0.261712932511430884D+00
    w(21) = 0.253636910481387185D+00
    w(22) = 0.16611584261479281D+00
    w(23) = 0.783630990508037449D-01
    w(24) = 0.274962713372148476D-01
    w(25) = 0.768092665770660459D-02
    w(26) = 0.313373786000304381D-02
    w(27) = 0.174733389581099482D-02
    w(28) = 0.41642095727577091D-03
    w(29) = 0.571885531470621903D-04
    w(30) = 0.586639457073896277D-05
    w(31) = 0.524482423744884136D-06
    w(32) = 0.294146671497083432D-07
    w(33) = 0.748907559239519284D-09
    w(34) = 0.809333688669950037D-11
    w(35) = 0.323016866782871498D-13
    w(36) = 0.332834739632930463D-16
    w(37) = 0.337304188079177058D-20

  else if ( n == 41 ) then

    w( 1) =   0.117725656974405367D-22
    w( 2) =   0.152506745534300636D-18
    w( 3) =   0.202183949965101288D-15
    w( 4) =   0.724614869051195508D-13
    w( 5) =   0.103121966469463034D-10
    w( 6) =   0.710371395169350952D-09
    w( 7) =   0.264376044449260516D-07
    w( 8) =   0.558982787078644997D-06
    w( 9) =   0.675628907134744976D-05
    w(10) =   0.512198007019776873D-04
    w(11) =   0.335013114947200879D-03
    w(12) =   0.249379691096933139D-02
    w(13) = - 0.25616995850607458D-01
    w(14) =   0.317007878644325588D-01
    w(15) =   0.125041498584003435D-02
    w(16) =   0.293244560924894295D-01
    w(17) =   0.799536390803302298D-01
    w(18) =   0.164543666806555251D+00
    w(19) =   0.258718519718241095D+00
    w(20) =   0.293588795735908566D+00
    w(21) =   0.997525375254611951D-01
    w(22) =   0.293588795735908566D+00
    w(23) =   0.258718519718241095D+00
    w(24) =   0.164543666806555251D+00
    w(25) =   0.799536390803302298D-01
    w(26) =   0.293244560924894295D-01
    w(27) =   0.125041498584003435D-02
    w(28) =   0.317007878644325588D-01
    w(29) = - 0.25616995850607458D-01
    w(30) =   0.249379691096933139D-02
    w(31) =   0.335013114947200879D-03
    w(32) =   0.512198007019776873D-04
    w(33) =   0.675628907134744976D-05
    w(34) =   0.558982787078644997D-06
    w(35) =   0.264376044449260516D-07
    w(36) =   0.710371395169350952D-09
    w(37) =   0.103121966469463034D-10
    w(38) =   0.724614869051195508D-13
    w(39) =   0.202183949965101288D-15
    w(40) =   0.152506745534300636D-18
    w(41) =   0.117725656974405367D-22

  else if ( n == 43 ) then

    w( 1) =   0.968100020641528185D-37
    w( 2) =   0.15516931262860431D-22
    w( 3) =   0.175937309107750992D-18
    w( 4) =   0.217337608710893738D-15
    w( 5) =   0.747837010380540069D-13
    w( 6) =   0.104028132097205732D-10
    w( 7) =   0.70903573389336778D-09
    w( 8) =   0.263481722999966618D-07
    w( 9) =   0.560127964848432175D-06
    w(10) =   0.680410934802210232D-05
    w(11) =   0.508343873102544037D-04
    w(12) =   0.32753080006610181D-03
    w(13) =   0.267479828788552937D-02
    w(14) = - 0.687704270963253854D-02
    w(15) =   0.119383201790913588D-01
    w(16) =   0.248083722871002796D-02
    w(17) =   0.29000335749726387D-01
    w(18) =   0.798689557875757008D-01
    w(19) =   0.164609842422580606D+00
    w(20) =   0.258535954731607738D+00
    w(21) =   0.292243810406117141D+00
    w(22) =   0.102730713753441829D+00
    w(23) =   0.292243810406117141D+00
    w(24) =   0.258535954731607738D+00
    w(25) =   0.164609842422580606D+00
    w(26) =   0.798689557875757008D-01
    w(27) =   0.29000335749726387D-01
    w(28) =   0.248083722871002796D-02
    w(29) =   0.119383201790913588D-01
    w(30) = - 0.687704270963253854D-02
    w(31) =   0.267479828788552937D-02
    w(32) =   0.32753080006610181D-03
    w(33) =   0.508343873102544037D-04
    w(34) =   0.680410934802210232D-05
    w(35) =   0.560127964848432175D-06
    w(36) =   0.263481722999966618D-07
    w(37) =   0.70903573389336778D-09
    w(38) =   0.104028132097205732D-10
    w(39) =   0.747837010380540069D-13
    w(40) =   0.217337608710893738D-15
    w(41) =   0.175937309107750992D-18
    w(42) =   0.15516931262860431D-22
    w(43) =   0.968100020641528185D-37

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 9, 19, 35, 37, 41 or 43.'
    stop

  end if

  return
end
subroutine hermite_genz_keister_lookup_weights_np ( n, np, p, w )

!*****************************************************************************80
!
!! HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS_NP sets weights for a Patterson rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+?, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and a final rule of order
!    35, 37, 41 or 43.
!
!    The precisions of these rules are P = 1, 5, 15, 29, 
!    with the final rule of precision 51, 55, 63 or 67.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, 35, 37, 41 and 43.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) np

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) w(n)

  call hermite_genz_keister_lookup_weights ( n, w )

  return
end
subroutine hermite_gk18_lookup_points ( n, x )

!*****************************************************************************80
!
!! HERMITE_GK18_LOOKUP_POINTS: abscissas of a Hermite Genz-Keister 18 rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+18, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and 37.
!
!    The precisions of these rules are P = 1, 5, 15, 29, and 55.
!
!    Some of the data in this function was kindly supplied directly by
!    Alan Genz on 24 April 2011.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be 1, 3, 9, 19, or 37.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x( 1) =   0.0000000000000000D+00

  else if ( n == 3 ) then

    x( 1) =  -1.2247448713915889D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   1.2247448713915889D+00

  else if ( n == 9 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -2.0232301911005157D+00
    x( 3) =  -1.2247448713915889D+00
    x( 4) =  -5.2403354748695763D-01
    x( 5) =   0.0000000000000000D+00
    x( 6) =   5.2403354748695763D-01
    x( 7) =   1.2247448713915889D+00
    x( 8) =   2.0232301911005157D+00
    x( 9) =   2.9592107790638380D+00

  else if ( n == 19 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.2665132620567876D+00
    x( 5) =  -2.0232301911005157D+00
    x( 6) =  -1.8357079751751868D+00
    x( 7) =  -1.2247448713915889D+00
    x( 8) =  -8.7004089535290285D-01
    x( 9) =  -5.2403354748695763D-01
    x(10) =   0.0000000000000000D+00
    x(11) =   5.2403354748695763D-01
    x(12) =   8.7004089535290285D-01
    x(13) =   1.2247448713915889D+00
    x(14) =   1.8357079751751868D+00
    x(15) =   2.0232301911005157D+00
    x(16) =   2.2665132620567876D+00
    x(17) =   2.9592107790638380D+00
    x(18) =   3.6677742159463378D+00
    x(19) =   4.4995993983103881D+00

  else if ( n == 37 ) then

    x( 1) =  -6.853200069757519D+00
    x( 2) =  -6.124527854622158D+00
    x( 3) =  -5.521865209868350D+00
    x( 4) =  -4.986551454150765D+00
    x( 5) =  -4.499599398310388D+00
    x( 6) =  -4.057956316089741D+00
    x( 7) =  -3.667774215946338D+00
    x( 8) =  -3.315584617593290D+00
    x( 9) =  -2.959210779063838D+00
    x(10) =  -2.597288631188366D+00
    x(11) =  -2.266513262056788D+00
    x(12) =  -2.023230191100516D+00
    x(13) =  -1.835707975175187D+00
    x(14) =  -1.561553427651873D+00
    x(15) =  -1.224744871391589D+00
    x(16) =  -0.870040895352903D+00
    x(17) =  -0.524033547486958D+00
    x(18) =  -0.214618180588171D+00
    x(19) =   0.000000000000000D+00
    x(20) =   0.214618180588171D+00
    x(21) =   0.524033547486958D+00
    x(22) =   0.870040895352903D+00
    x(23) =   1.224744871391589D+00
    x(24) =   1.561553427651873D+00
    x(25) =   1.835707975175187D+00
    x(26) =   2.023230191100516D+00
    x(27) =   2.266513262056788D+00
    x(28) =   2.597288631188366D+00
    x(29) =   2.959210779063838D+00
    x(30) =   3.315584617593290D+00
    x(31) =   3.667774215946338D+00
    x(32) =   4.057956316089741D+00
    x(33) =   4.499599398310388D+00
    x(34) =   4.986551454150765D+00
    x(35) =   5.521865209868350D+00
    x(36) =   6.124527854622158D+00
    x(37) =   6.853200069757519D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GK18_LOOKUP_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 9, 19, or 37.'
    stop

  end if

  return
end
subroutine hermite_gk22_lookup_points ( n, x )

!*****************************************************************************80
!
!! HERMITE_GK22_LOOKUP_POINTS: abscissas of a Genz-Keister 22 Hermite rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+22, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and 41.
!
!    The precisions of these rules are P = 1, 5, 15, 29, and 63.
!
!    Some of the data in this function was kindly supplied directly by
!    Alan Genz on 24 April 2011.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, and 41.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x( 1) =   0.0000000000000000D+00

  else if ( n == 3 ) then

    x( 1) =  -1.2247448713915889D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   1.2247448713915889D+00

  else if ( n == 9 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -2.0232301911005157D+00
    x( 3) =  -1.2247448713915889D+00
    x( 4) =  -0.52403354748695763D+00
    x( 5) =   0.0000000000000000D+00
    x( 6) =   0.52403354748695763D+00
    x( 7) =   1.2247448713915889D+00
    x( 8) =   2.0232301911005157D+00
    x( 9) =   2.9592107790638380D+00

  else if ( n == 19 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.2665132620567876D+00
    x( 5) =  -2.0232301911005157D+00
    x( 6) =  -1.8357079751751868D+00
    x( 7) =  -1.2247448713915889D+00
    x( 8) =  -0.87004089535290285D+00
    x( 9) =  -0.52403354748695763D+00
    x(10) =   0.0000000000000000D+00
    x(11) =   0.52403354748695763D+00
    x(12) =   0.87004089535290285D+00
    x(13) =   1.2247448713915889D+00
    x(14) =   1.8357079751751868D+00
    x(15) =   2.0232301911005157D+00
    x(16) =   2.2665132620567876D+00
    x(17) =   2.9592107790638380D+00
    x(18) =   3.6677742159463378D+00
    x(19) =   4.4995993983103881D+00

  else if ( n == 41 ) then

    x( 1) =  -7.251792998192644D+00
    x( 2) =  -6.547083258397540D+00
    x( 3) =  -5.961461043404500D+00
    x( 4) =  -5.437443360177798D+00
    x( 5) =  -4.953574342912980D+00
    x( 6) =  -4.4995993983103881D+00
    x( 7) =  -4.070919267883068D+00
    x( 8) =  -3.6677742159463378D+00
    x( 9) =  -3.296114596212218D+00
    x(10) =  -2.9592107790638380D+00
    x(11) =  -2.630415236459871D+00
    x(12) =  -2.2665132620567876D+00
    x(13) =  -2.043834754429505D+00
    x(14) =  -2.0232301911005157D+00
    x(15) =  -1.8357079751751868D+00
    x(16) =  -1.585873011819188D+00
    x(17) =  -1.2247448713915889D+00
    x(18) =  -0.87004089535290285D+00
    x(19) =  -0.52403354748695763D+00
    x(20) =  -0.195324784415805D+00
    x(21) =   0.0000000000000000D+00
    x(22) =   0.195324784415805D+00
    x(23) =   0.52403354748695763D+00
    x(24) =   0.87004089535290285D+00
    x(25) =   1.2247448713915889D+00
    x(26) =   1.585873011819188D+00
    x(27) =   1.8357079751751868D+00
    x(28) =   2.0232301911005157D+00
    x(29) =   2.043834754429505D+00
    x(30) =   2.2665132620567876D+00
    x(31) =   2.630415236459871D+00
    x(32) =   2.9592107790638380D+00
    x(33) =   3.296114596212218D+00
    x(34) =   3.6677742159463378D+00
    x(35) =   4.070919267883068D+00
    x(36) =   4.4995993983103881D+00
    x(37) =   4.953574342912980D+00
    x(38) =   5.437443360177798D+00
    x(39) =   5.961461043404500D+00
    x(40) =   6.547083258397540D+00
    x(41) =   7.251792998192644D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GK22_LOOKUP_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 9, 19, or 41.'
    stop

  end if

  return
end
subroutine hermite_gk24_lookup_points ( n, x )

!*****************************************************************************80
!
!! HERMITE_GK24_LOOKUP_POINTS: abscissas of a Genz-Keister 24 Hermite rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+24, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and 43.
!
!    The precisions of these rules are P = 1, 5, 15, 29, and 67.
!
!    Some of the data in this function was kindly supplied directly by
!    Alan Genz on 24 April 2011.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, and 43.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x( 1) =   0.0000000000000000D+00

  else if ( n == 3 ) then

    x( 1) =  -1.2247448713915889D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   1.2247448713915889D+00

  else if ( n == 9 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -2.0232301911005157D+00
    x( 3) =  -1.2247448713915889D+00
    x( 4) =  -0.52403354748695763D+00
    x( 5) =   0.0000000000000000D+00
    x( 6) =   0.52403354748695763D+00
    x( 7) =   1.2247448713915889D+00
    x( 8) =   2.0232301911005157D+00
    x( 9) =   2.9592107790638380D+00

  else if ( n == 19 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.2665132620567876D+00
    x( 5) =  -2.0232301911005157D+00
    x( 6) =  -1.8357079751751868D+00
    x( 7) =  -1.2247448713915889D+00
    x( 8) =  -0.87004089535290285D+00
    x( 9) =  -0.52403354748695763D+00
    x(10) =   0.0000000000000000D+00
    x(11) =   0.52403354748695763D+00
    x(12) =   0.87004089535290285D+00
    x(13) =   1.2247448713915889D+00
    x(14) =   1.8357079751751868D+00
    x(15) =   2.0232301911005157D+00
    x(16) =   2.2665132620567876D+00
    x(17) =   2.9592107790638380D+00
    x(18) =   3.6677742159463378D+00
    x(19) =   4.4995993983103881D+00

  else if ( n == 43 ) then

    x( 1) = -10.167574994881873D+00
    x( 2) =  -7.231746029072501D+00
    x( 3) =  -6.535398426382995D+00
    x( 4) =  -5.954781975039809D+00
    x( 5) =  -5.434053000365068D+00
    x( 6) =  -4.952329763008589D+00
    x( 7) =  -4.4995993983103881D+00
    x( 8) =  -4.071335874253583D+00
    x( 9) =  -3.6677742159463378D+00
    x(10) =  -3.295265921534226D+00
    x(11) =  -2.9592107790638380D+00
    x(12) =  -2.633356763661946D+00
    x(13) =  -2.2665132620567876D+00
    x(14) =  -2.089340389294661D+00
    x(15) =  -2.0232301911005157D+00
    x(16) =  -1.8357079751751868D+00
    x(17) =  -1.583643465293944D+00
    x(18) =  -1.2247448713915889D+00
    x(19) =  -0.87004089535290285D+00
    x(20) =  -0.52403354748695763D+00
    x(21) =  -0.196029453662011D+00
    x(22) =   0.0000000000000000D+00
    x(23) =   0.196029453662011D+00
    x(24) =   0.52403354748695763D+00
    x(25) =   0.87004089535290285D+00
    x(26) =   1.2247448713915889D+00
    x(27) =   1.583643465293944D+00
    x(28) =   1.8357079751751868D+00
    x(29) =   2.0232301911005157D+00
    x(30) =   2.089340389294661D+00
    x(31) =   2.2665132620567876D+00
    x(32) =   2.633356763661946D+00
    x(33) =   2.9592107790638380D+00
    x(34) =   3.295265921534226D+00
    x(35) =   3.6677742159463378D+00
    x(36) =   4.071335874253583D+00
    x(37) =   4.4995993983103881D+00
    x(38) =   4.952329763008589D+00
    x(39) =   5.434053000365068D+00
    x(40) =   5.954781975039809D+00
    x(41) =   6.535398426382995D+00
    x(42) =   7.231746029072501D+00
    x(43) =  10.167574994881873D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GK24_LOOKUP_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 9, 19, or 43.'
    stop

  end if

  return
end
subroutine hermite_integral ( n, value )

!*****************************************************************************80
!
!! HERMITE_INTEGRAL evaluates a monomial Hermite integral.
!
!  Discussion:
!
!    H(n) = Integral ( -oo < x < +oo ) x^n exp(-x^2) dx
!
!    H(n) is 0 for n odd.
!
!    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
!
!    Note that it is difficult to correctly evaluate the double factorial
!    quantity.  It "blows up", though somewhat less rapidly than the
!    standard factorial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
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

  if ( n < 0 ) then

    value = - huge ( value )

  else if ( mod ( n, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = r8_factorial2 ( n - 1 ) * sqrt ( pi ) / 2.0D+00**( n / 2 )

  end if

  return
end
subroutine hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT sets up a divided difference table from Hermite data.
!
!  Discussion:
!
!    The polynomial represented by the divided difference table can be
!    evaluated by calling HERMITE_INTERPOLANT_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data 
!    ( X(I), Y(I), YP(I) ).
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!    These values must be distinct.
!
!    Input, real ( kind = 8 ) Y(N), YP(N), the function and 
!    derivative values.
!
!    Output, real ( kind = 8 ) XD(2*N), YD(2*N), the divided difference table
!    for the interpolant value.
!
!    Output, real ( kind = 8 ) XDP(2*N-1), YDP(2*N-1), the divided difference 
!    table for the interpolant derivative.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ndp
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xd(2*n)
  real ( kind = 8 ) xdp(2*n-1)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yd(2*n)
  real ( kind = 8 ) ydp(2*n-1)
  real ( kind = 8 ) yp(n)
!
!  Copy the data:
!
  nd = 2 * n
  xd(1:nd-1:2) = x(1:n)
  xd(2:nd  :2) = x(1:n)
!
!  Carry out the first step of differencing.
!
  yd(1) = y(1)
  yd(3:nd-1:2) = ( y(2:n) - y(1:n-1) ) / ( x(2:n) - x(1:n-1) )
  yd(2:nd  :2) = yp(1:n)  
!
!  Carry out the remaining steps in the usual way.
!
  do i = 3, nd
    do j = nd, i, -1

      yd(j) = ( yd(j) - yd(j-1) ) / ( xd(j) - xd(j+1-i) )

    end do

  end do
!
!  Compute the difference table for the derivative.
!
  call dif_deriv ( nd, xd, yd, ndp, xdp, ydp )

  return
end
subroutine hermite_interpolant_rule ( n, a, b, x, w )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT_RULE: quadrature rule for a Hermite interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of abscissas.
!
!    Input, real ( kind = 8 ) A, B, the integration limits.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(2*N), the quadrature coefficients, given as
!    pairs for function and derivative values at each abscissa.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) a_value
  real ( kind = 8 ) b
  real ( kind = 8 ) b_value
  real ( kind = 8 ) c(0:2*n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nd
  real ( kind = 8 ) w(2*n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xd(2*n)
  real ( kind = 8 ) xdp(2*n-1)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yd(2*n)
  real ( kind = 8 ) ydp(2*n-1)
  real ( kind = 8 ) yp(n)

  nd = 2 * n

  k = 0

  do i = 1, n

    k = k + 1
    y(1:n) = 0.0D+00
    y(i) = 1.0D+00
    yp(1:n) = 0.0D+00
    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
    call dif_to_r8poly ( nd, xd, yd, c )
    call r8poly_ant_val ( n, c, a, a_value )
    call r8poly_ant_val ( n, c, b, b_value )
    w(k) = b_value - a_value

    k = k + 1
    y(1:n) = 0.0D+00
    yp(1:n) = 0.0D+00
    yp(i) = 1.0D+00
    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
    call dif_to_r8poly ( nd, xd, yd, c )
    call r8poly_ant_val ( n, c, a, a_value )
    call r8poly_ant_val ( n, c, b, b_value )
    w(k) = b_value - a_value

  end do

  return
end
subroutine hermite_interpolant_value ( nd, xd, yd, xdp, ydp, nv, xv, yv, yvp )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT_VALUE evaluates the Hermite interpolant polynomial.
!
!  Discussion:
!
!    In fact, this function will evaluate an arbitrary polynomial that is
!    represented by a difference table.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the order of the difference table.
!
!    Input, real ( kind = 8 ) XD(ND), YD(ND), the difference table for the
!    interpolant value.
!
!    Input, real ( kind = 8 ) XDP(ND-1), YDP(ND-1), the difference table for
!    the interpolant derivative.
!
!    Input, integer ( kind = 4 ) NV, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XV(NV), the evaluation points.
!
!    Output, real ( kind = 8 ) YV(NV), YVP(NV), the value of the interpolant and
!    its derivative at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nv

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ndp
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xdp(nd-1)
  real ( kind = 8 ) xv(nv)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) ydp(nd-1)
  real ( kind = 8 ) yv(nv)
  real ( kind = 8 ) yvp(nv)

  yv(1:nv) = yd(nd)
  do i = nd - 1, 1, -1
    yv(1:nv) = yd(i) + ( xv(1:nv) - xd(i) ) * yv(1:nv)
  end do

  ndp = nd - 1

  yvp(1:nv) = ydp(ndp)
  do i = ndp - 1, 1, -1
    yvp(1:nv) = ydp(i) + ( xv(1:nv) - xdp(i) ) * yvp(1:nv)
  end do

  return
end
subroutine hermite_lookup_points ( order, points )

!*****************************************************************************80
!
!! HERMITE_LOOKUP_POINTS returns the abscissas of a Hermite rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on (-oo,+oo).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    Legal values are 1, 3, 7, 15, 31, 63 and 127.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) points(order)
  real ( kind = 8 ), save, dimension ( 1 ) :: x_001 = (/ &
    0.0D+00 /)
  real ( kind = 8 ), save, dimension ( 3 ) :: x_003 = (/ &
   -0.122474487139158904909864203735D+01, &
    0.0D+00, &
    0.122474487139158904909864203735D+01 /)
  real ( kind = 8 ), save, dimension ( 7 ) :: x_007 = (/ &
   -0.265196135683523349244708200652D+01, &
   -0.167355162876747144503180139830D+01, &
   -0.816287882858964663038710959027D+00, &
    0.0D+00, &
    0.816287882858964663038710959027D+00, &
    0.167355162876747144503180139830D+01, &
    0.265196135683523349244708200652D+01 /)
  real ( kind = 8 ), save, dimension ( 15 ) :: x_015 = (/ &
   -0.449999070730939155366438053053D+01, &
   -0.366995037340445253472922383312D+01, &
   -0.296716692790560324848896036355D+01, &
   -0.232573248617385774545404479449D+01, &
   -0.171999257518648893241583152515D+01, &
   -0.113611558521092066631913490556D+01, &
   -0.565069583255575748526020337198D+00, &
    0.0D+00, &
    0.565069583255575748526020337198D+00, &
    0.113611558521092066631913490556D+01, &
    0.171999257518648893241583152515D+01, &
    0.232573248617385774545404479449D+01, &
    0.296716692790560324848896036355D+01, &
    0.366995037340445253472922383312D+01, &
    0.449999070730939155366438053053D+01 /)
  real ( kind = 8 ), save, dimension ( 31 ) :: x_031 = (/ &
   -6.9956801237185402753248521473232D+00, &
   -6.2750787049428601427036567812530D+00, &
   -5.6739614446185883296332558789276D+00, &
   -5.1335955771123807045862968913996D+00, &
   -4.6315595063128599420667997654336D+00, &
   -4.1562717558181451724831352315314D+00, &
   -3.7007434032314694224497164589673D+00, &
   -3.2603207323135408104645401509648D+00, &
   -2.8316804533902054557015640151425D+00, &
   -2.4123177054804201051740184582119D+00, &
   -2.0002585489356389657975562598571D+00, &
   -1.5938858604721398261388419455550D+00, &
   -1.1918269983500464260821358649242D+00, &
   -0.79287697691530893968593032998830D+00, &
   -0.39594273647142311094670041663436D+00, &
    0.0000000000000000000000000000000D+00, &
    0.39594273647142311094670041663436D+00, &
    0.79287697691530893968593032998830D+00, &
    1.1918269983500464260821358649242D+00, &
    1.5938858604721398261388419455550D+00, &
    2.0002585489356389657975562598571D+00, &
    2.4123177054804201051740184582119D+00, &
    2.8316804533902054557015640151425D+00, &
    3.2603207323135408104645401509648D+00, &
    3.7007434032314694224497164589673D+00, &
    4.1562717558181451724831352315314D+00, &
    4.6315595063128599420667997654336D+00, &
    5.1335955771123807045862968913996D+00, &
    5.6739614446185883296332558789276D+00, &
    6.2750787049428601427036567812530D+00, &
    6.9956801237185402753248521473232D+00 /)
  real ( kind = 8 ), save, dimension ( 63 ) :: x_063 = (/ &
   -10.435499877854168053468115427285D+00, &
   -9.8028759912974963635223935286507D+00, &
   -9.2792019543050391319404745506496D+00, &
   -8.8118581437284546442526628275570D+00, &
   -8.3807683451863219343010651043788D+00, &
   -7.9755950801420373181541806298501D+00, &
   -7.5901395198641066762479783194468D+00, &
   -7.2203167078889678461161324222529D+00, &
   -6.8632544331795368527353285876066D+00, &
   -6.5168348106821160605273395854042D+00, &
   -6.1794379922705969862418461787263D+00, &
   -5.8497884000810673462526582961482D+00, &
   -5.5268572526403031425047575122840D+00, &
   -5.2097979830408354861575136416263D+00, &
   -4.8979018644975742350745099214868D+00, &
   -4.5905665744435190229271294569091D+00, &
   -4.2872733352824404031727616199454D+00, &
   -3.9875699104197157485227052068068D+00, &
   -3.6910577000963465117322810559754D+00, &
   -3.3973817713303911852755941806287D+00, &
   -3.1062230279282566329138616746036D+00, &
   -2.8172919672837977750747135657355D+00, &
   -2.5303236304712010926855221718499D+00, &
   -2.2450734604812066298995918179330D+00, &
   -1.9613138583081485293922008411321D+00, &
   -1.6788312791720137520802800622638D+00, &
   -1.3974237486049625107570752063702D+00, &
   -1.1168987050996462690510970277840D+00, &
   -0.83707109558947615977737795461293D+00, &
   -0.55776166427908221668763665253822D+00, &
   -0.27879538567115223986687628627202D+00, &
    0.00000000000000000000000000000000D+00, &
    0.27879538567115223986687628627202D+00, &
    0.55776166427908221668763665253822D+00, &
    0.83707109558947615977737795461293D+00, &
    1.1168987050996462690510970277840D+00, &
    1.3974237486049625107570752063702D+00, &
    1.6788312791720137520802800622638D+00, &
    1.9613138583081485293922008411321D+00, &
    2.2450734604812066298995918179330D+00, &
    2.5303236304712010926855221718499D+00, &
    2.8172919672837977750747135657355D+00, &
    3.1062230279282566329138616746036D+00, &
    3.3973817713303911852755941806287D+00, &
    3.6910577000963465117322810559754D+00, &
    3.9875699104197157485227052068068D+00, &
    4.2872733352824404031727616199454D+00, &
    4.5905665744435190229271294569091D+00, &
    4.8979018644975742350745099214868D+00, &
    5.2097979830408354861575136416263D+00, &
    5.5268572526403031425047575122840D+00, &
    5.8497884000810673462526582961482D+00, &
    6.1794379922705969862418461787263D+00, &
    6.5168348106821160605273395854042D+00, &
    6.8632544331795368527353285876066D+00, &
    7.2203167078889678461161324222529D+00, &
    7.5901395198641066762479783194468D+00, &
    7.9755950801420373181541806298501D+00, &
    8.3807683451863219343010651043788D+00, &
    8.8118581437284546442526628275570D+00, &
    9.2792019543050391319404745506496D+00, &
    9.8028759912974963635223935286507D+00, &
    10.435499877854168053468115427285D+00 /)
  real ( kind = 8 ), save, dimension ( 127 ) :: x_127 = (/ &
   -15.228338148167350978246954433464D+00, &
   -14.669595158833972632746354112896D+00, &
   -14.209085995284870755168244250887D+00, &
   -13.799722290211676634645246746673D+00, &
   -13.423518590070950062438258321855D+00, &
   -13.071208660474601901583995439649D+00, &
   -12.737235652415686338138003924072D+00, &
   -12.417939378869715805445879624069D+00, &
   -12.110749020947747600132123508132D+00, &
   -11.813772198267727195134584136191D+00, &
   -11.525565112572696599167888588564D+00, &
   -11.244994583785543445194384194300D+00, &
   -10.971150569840247423423040263881D+00, &
   -10.703288201027481347670940744690D+00, &
   -10.440787957772772867742591798027D+00, &
   -10.183127473450343888624126450357D+00, &
   -9.9298610495114250736847004273684D+00, &
   -9.6806044412474728038150712732737D+00, &
   -9.4350233389881650135019598506287D+00, &
   -9.1928244988460305715774195052527D+00, &
   -8.9537488108565404323807890169970D+00, &
   -8.7175658087076307363833999548548D+00, &
   -8.4840692689832473326097180339984D+00, &
   -8.2530736454457156579694124243888D+00, &
   -8.0244111514703375578594739796798D+00, &
   -7.7979293513870105420829120455591D+00, &
   -7.5734891556083454022834960763301D+00, &
   -7.3509631392269052701961258043733D+00, &
   -7.1302341220350710668064025713431D+00, &
   -6.9111939615465713197465633109366D+00, &
   -6.6937425208758294190074417381666D+00, &
   -6.4777867811645365448144903821487D+00, &
   -6.2632400742737354345609723857092D+00, &
   -6.0500214161419845694465474482388D+00, &
   -5.8380549248774187386601690807757D+00, &
   -5.6272693105464816659423455794909D+00, &
   -5.4175974259243240722848425872924D+00, &
   -5.2089758693153983587570258372239D+00, &
   -5.0013446320386360038520809107373D+00, &
   -4.7946467843764925009748509930857D+00, &
   -4.5888281947698372951606485031212D+00, &
   -4.3838372778464736294253744407459D+00, &
   -4.1796247675352031349421189892408D+00, &
   -3.9761435120673355916035814195920D+00, &
   -3.7733482881250526721004678400057D+00, &
   -3.5711956317782180447199756485249D+00, &
   -3.3696436841717397896643629240035D+00, &
   -3.1686520501953630191857798261495D+00, &
   -2.9681816685955910267761649521505D+00, &
   -2.7681946921824058801226545958892D+00, &
   -2.5686543769473501723144013022363D+00, &
   -2.3695249790490401080012474645702D+00, &
   -2.1707716587411506879498498083695D+00, &
   -1.9723603904195020079324743227565D+00, &
   -1.7742578780516791584676442103681D+00, &
   -1.5764314753267801315519597621879D+00, &
   -1.3788491099261778091441557053728D+00, &
   -1.1814792113700685848678583598423D+00, &
   -0.98429064194027277726568984213773D+00, &
   -0.78725263021825034151596831878971D+00, &
   -0.59033470680942102142230439346102D+00, &
   -0.39350664185130136568037826200185D+00, &
   -0.19673838392423251964272239737078D+00, &
    0.0000000000000000000000000000000D+00, &
    0.19673838392423251964272239737078D+00, &
    0.39350664185130136568037826200185D+00, &
    0.59033470680942102142230439346102D+00, &
    0.78725263021825034151596831878971D+00, &
    0.98429064194027277726568984213773D+00, &
    1.1814792113700685848678583598423D+00, &
    1.3788491099261778091441557053728D+00, &
    1.5764314753267801315519597621879D+00, &
    1.7742578780516791584676442103681D+00, &
    1.9723603904195020079324743227565D+00, &
    2.1707716587411506879498498083695D+00, &
    2.3695249790490401080012474645702D+00, &
    2.5686543769473501723144013022363D+00, &
    2.7681946921824058801226545958892D+00, &
    2.9681816685955910267761649521505D+00, &
    3.1686520501953630191857798261495D+00, &
    3.3696436841717397896643629240035D+00, &
    3.5711956317782180447199756485249D+00, &
    3.7733482881250526721004678400057D+00, &
    3.9761435120673355916035814195920D+00, &
    4.1796247675352031349421189892408D+00, &
    4.3838372778464736294253744407459D+00, &
    4.5888281947698372951606485031212D+00, &
    4.7946467843764925009748509930857D+00, &
    5.0013446320386360038520809107373D+00, &
    5.2089758693153983587570258372239D+00, &
    5.4175974259243240722848425872924D+00, &
    5.6272693105464816659423455794909D+00, &
    5.8380549248774187386601690807757D+00, &
    6.0500214161419845694465474482388D+00, &
    6.2632400742737354345609723857092D+00, &
    6.4777867811645365448144903821487D+00, &
    6.6937425208758294190074417381666D+00, &
    6.9111939615465713197465633109366D+00, &
    7.1302341220350710668064025713431D+00, &
    7.3509631392269052701961258043733D+00, &
    7.5734891556083454022834960763301D+00, &
    7.7979293513870105420829120455591D+00, &
    8.0244111514703375578594739796798D+00, &
    8.2530736454457156579694124243888D+00, &
    8.4840692689832473326097180339984D+00, &
    8.7175658087076307363833999548548D+00, &
    8.9537488108565404323807890169970D+00, &
    9.1928244988460305715774195052527D+00, &
    9.4350233389881650135019598506287D+00, &
    9.6806044412474728038150712732737D+00, &
    9.9298610495114250736847004273684D+00, &
    10.183127473450343888624126450357D+00, &
    10.440787957772772867742591798027D+00, &
    10.703288201027481347670940744690D+00, &
    10.971150569840247423423040263881D+00, &
    11.244994583785543445194384194300D+00, &
    11.525565112572696599167888588564D+00, &
    11.813772198267727195134584136191D+00, &
    12.110749020947747600132123508132D+00, &
    12.417939378869715805445879624069D+00, &
    12.737235652415686338138003924072D+00, &
    13.071208660474601901583995439649D+00, &
    13.423518590070950062438258321855D+00, &
    13.799722290211676634645246746673D+00, &
    14.209085995284870755168244250887D+00, &
    14.669595158833972632746354112896D+00, &
    15.228338148167350978246954433464D+00 /)

  if ( order == 1 ) then

    points(1:order) = x_001(1:order)

  else if ( order == 3 ) then

    points(1:order) = x_003(1:order)

  else if ( order == 7 ) then

    points(1:order) = x_007(1:order)

  else if ( order == 15 ) then

    points(1:order) = x_015(1:order)

  else if ( order == 31 ) then

    points(1:order) = x_031(1:order)

  else if ( order == 63 ) then

    points(1:order) = x_063(1:order)

  else if ( order == 127 ) then

    points(1:order) = x_127(1:order)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_LOOKUP_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Unexpected value of ORDER = ', order
    stop

  end if

  return
end
subroutine hermite_lookup_weights ( n, w )

!*****************************************************************************80
!
!! HERMITE_LOOKUP_WEIGHTS returns weights for Hermite quadrature rules.
!
!  Discussion:
!
!    The allowed orders are 1, 3, 7, 15, 31, 63, and 127.
!
!    The weights are positive, symmetric and should sum to SQRT(PI).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2011
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
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be 1, 3, 7, 15, 31, 63 or 127.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)

  if ( n == 1 ) then

    w(1) = 1.77245385090551602729816748334D+00

  else if ( n == 3 ) then

    w(1) = 0.295408975150919337883027913890D+00
    w(2) = 0.118163590060367735153211165556D+01
    w(3) = 0.295408975150919337883027913890D+00


  else if ( n == 7 ) then

    w(1) = 0.971781245099519154149424255939D-03
    w(2) = 0.545155828191270305921785688417D-01
    w(3) = 0.425607252610127800520317466666D+00
    w(4) = 0.810264617556807326764876563813D+00
    w(5) = 0.425607252610127800520317466666D+00
    w(6) = 0.545155828191270305921785688417D-01
    w(7) = 0.971781245099519154149424255939D-03

  else if ( n == 15 ) then

    w(1) =  0.152247580425351702016062666965D-08
    w(2) =  0.105911554771106663577520791055D-05
    w(3) =  0.100004441232499868127296736177D-03
    w(4) =  0.277806884291277589607887049229D-02
    w(5) =  0.307800338725460822286814158758D-01
    w(6) =  0.158488915795935746883839384960D+00
    w(7) =  0.412028687498898627025891079568D+00
    w(8) =  0.564100308726417532852625797340D+00
    w(9) =  0.412028687498898627025891079568D+00
    w(10) = 0.158488915795935746883839384960D+00
    w(11) = 0.307800338725460822286814158758D-01
    w(12) = 0.277806884291277589607887049229D-02
    w(13) = 0.100004441232499868127296736177D-03
    w(14) = 0.105911554771106663577520791055D-05
    w(15) = 0.152247580425351702016062666965D-08

  else if ( n == 31 ) then
 
    w(1) = 4.61896839446420502132944426974D-22
    w(2) = 5.11060900792715640739422641166D-18
    w(3) = 5.89955649875387299038431589378D-15
    w(4) = 1.86037352145214652437380892603D-12
    w(5) = 2.35249200320864163398597795323D-10
    w(6) = 1.46119883449105307352780323055D-08
    w(7) = 5.04371255893979974253745671633D-07
    w(8) = 0.0000104986027576756063228123279208D+00
    w(9) = 0.000139520903950470433823653754396D+00
    w(10) = 0.00123368330730688826551750402319D+00
    w(11) = 0.00748279991403519848345678003016D+00
    w(12) = 0.0318472307313003327772087235339D+00
    w(13) = 0.0967179481608704535580338478886D+00
    w(14) = 0.212132788668764779877735637343D+00
    w(15) = 0.338772657894107724675701919878D+00
    w(16) = 0.395778556098609545141783810611D+00
    w(17) = 0.338772657894107724675701919878D+00
    w(18) = 0.212132788668764779877735637343D+00
    w(19) = 0.0967179481608704535580338478886D+00
    w(20) = 0.0318472307313003327772087235339D+00
    w(21) = 0.00748279991403519848345678003016D+00
    w(22) = 0.00123368330730688826551750402319D+00
    w(23) = 0.000139520903950470433823653754396D+00
    w(24) = 0.0000104986027576756063228123279208D+00
    w(25) = 5.04371255893979974253745671633D-07
    w(26) = 1.46119883449105307352780323055D-08
    w(27) = 2.35249200320864163398597795323D-10
    w(28) = 1.86037352145214652437380892603D-12
    w(29) = 5.89955649875387299038431589378D-15
    w(30) = 5.11060900792715640739422641166D-18
    w(31) = 4.61896839446420502132944426974D-22

  else if ( n == 63 ) then

    w(1) = 3.70992064349030055823376157823D-48
    w(2) = 1.04007786152246672212559599908D-42
    w(3) = 1.97968047083199197900260998813D-38
    w(4) = 8.46874781919035663281042885251D-35
    w(5) = 1.30713059308206243904769877879D-31
    w(6) = 9.34378371756582396450246862195D-29
    w(7) = 3.60274266352851638202340658522D-26
    w(8) = 8.29638631162099766157527065317D-24
    w(9) = 1.22666299091434557721622529775D-21
    w(10) = 1.22884356288353036990240371039D-19
    w(11) = 8.69255369584585252225619256428D-18
    w(12) = 4.48570586893158184069444097978D-16
    w(13) = 1.73358179557891044383064226749D-14
    w(14) = 5.1265062385197846998384009333D-13
    w(15) = 1.18089218445696923817995132237D-11
    w(16) = 2.15086982978749617679069862879D-10
    w(17) = 3.13719295353830786449435629291D-09
    w(18) = 3.70416259848969809883356560995D-08
    w(19) = 3.57347329499908777461505032558D-07
    w(20) = 2.83931144984692884712301165567D-06
    w(21) = 0.0000187091130037887216027832755405D+00
    w(22) = 0.000102848808006856425543062213642D+00
    w(23) = 0.000474117026103206754395975199216D+00
    w(24) = 0.0018409222622442103760124297917D+00
    w(25) = 0.00604360445513757113209247151533D+00
    w(26) = 0.0168292991996521044559098701555D+00
    w(27) = 0.0398582640278170328649908688578D+00
    w(28) = 0.0804670879942008323850873860195D+00
    w(29) = 0.138719508176584635072239096351D+00
    w(30) = 0.204486953468973988225911656103D+00
    w(31) = 0.25799889943138332612723393346D+00
    w(32) = 0.278766948849251654365527505911D+00
    w(33) = 0.25799889943138332612723393346D+00
    w(34) = 0.204486953468973988225911656103D+00
    w(35) = 0.138719508176584635072239096351D+00
    w(36) = 0.0804670879942008323850873860195D+00
    w(37) = 0.0398582640278170328649908688578D+00
    w(38) = 0.0168292991996521044559098701555D+00
    w(39) = 0.00604360445513757113209247151533D+00
    w(40) = 0.0018409222622442103760124297917D+00
    w(41) = 0.000474117026103206754395975199216D+00
    w(42) = 0.000102848808006856425543062213642D+00
    w(43) = 0.0000187091130037887216027832755405D+00
    w(44) = 2.83931144984692884712301165567D-06
    w(45) = 3.57347329499908777461505032558D-07
    w(46) = 3.70416259848969809883356560995D-08
    w(47) = 3.13719295353830786449435629291D-09
    w(48) = 2.15086982978749617679069862879D-10
    w(49) = 1.18089218445696923817995132237D-11
    w(50) = 5.1265062385197846998384009333D-13
    w(51) = 1.73358179557891044383064226749D-14
    w(52) = 4.48570586893158184069444097978D-16
    w(53) = 8.69255369584585252225619256428D-18
    w(54) = 1.22884356288353036990240371039D-19
    w(55) = 1.22666299091434557721622529775D-21
    w(56) = 8.29638631162099766157527065317D-24
    w(57) = 3.60274266352851638202340658522D-26
    w(58) = 9.34378371756582396450246862195D-29
    w(59) = 1.30713059308206243904769877879D-31
    w(60) = 8.46874781919035663281042885251D-35
    w(61) = 1.97968047083199197900260998813D-38
    w(62) = 1.04007786152246672212559599908D-42
    w(63) = 3.70992064349030055823376157823D-48

  else if ( n == 127 ) then  

    w(1) = 1.25044975770895101066558695394D-101
    w(2) = 1.72727980594728851329952877284D-94
    w(3) = 8.93216815722645216635320162557D-89
    w(4) = 7.7306185241134158744827181222D-84
    w(5) = 2.01439576527109443920782513994D-79
    w(6) = 2.15037147336771602203551878273D-75
    w(7) = 1.13419242086298913875376620343D-71
    w(8) = 3.34891390118992716444169809114D-68
    w(9) = 6.04865489642049179016214753843D-65
    w(10) = 7.13750929465743002965122123123D-62
    w(11) = 5.78845633750656959788340019085D-59
    w(12) = 3.3581166223962736386929935773D-56
    w(13) = 1.4394641949298720336141068619D-53
    w(14) = 4.68218083833618292793410025836D-51
    w(15) = 1.18170544407210392716367665268D-48
    w(16) = 2.35816591560823143778744566357D-46
    w(17) = 3.78144279409152203964384313149D-44
    w(18) = 4.9411031115925407477456893331D-42
    w(19) = 5.32553037755907921458489847863D-40
    w(20) = 4.78543906802804099967221020647D-38
    w(21) = 3.61918834460649868835433546523D-36
    w(22) = 2.3232083386415854084664074623D-34
    w(23) = 1.27533314110484056196532640642D-32
    w(24) = 6.02777538509463291699314327193D-31
    w(25) = 2.4679773241854004762148469348D-29
    w(26) = 8.8019567691972403392314252914D-28
    w(27) = 2.74824892121260880467531987939D-26
    w(28) = 7.54682189033203465872349657723D-25
    w(29) = 1.83031346363374264415878982576D-23
    w(30) = 3.93559908609832906838466602268D-22
    w(31) = 7.52931616388155067444192947319D-21
    w(32) = 1.28579977867628696999762170542D-19
    w(33) = 1.96593268885070384943390296306D-18
    w(34) = 2.69865119072980851232572568063D-17
    w(35) = 3.33444143033026256341061235315D-16
    w(36) = 3.71733031252663248624409938613D-15
    w(37) = 3.74739544729563577089986076081D-14
    w(38) = 3.42300944935037851188976963928D-13
    w(39) = 2.83853037250817094975750489262D-12
    w(40) = 2.14069202905212884993201956606D-11
    w(41) = 1.47063312734774830028408333227D-10
    w(42) = 9.21739409677215086782446989876D-10
    w(43) = 5.27816639371369729333040255118D-09
    w(44) = 2.76504970450371674155194812923D-08
    w(45) = 1.32678558425807549298485884004D-07
    w(46) = 5.83809442762947462901022315301D-07
    w(47) = 2.35815617248490159838145978859D-06
    w(48) = 8.75244680345528247507614056972D-06
    w(49) = 0.0000298767905360019901790649251988D+00
    w(50) = 0.0000938744357203646866361259710004D+00
    w(51) = 0.000271707626280157286781639661883D+00
    w(52) = 0.000724939297427239633212185817821D+00
    w(53) = 0.0017841208326818955520088211458D+00
    w(54) = 0.00405248551861722466559241860023D+00
    w(55) = 0.00850002630418086349941683729112D+00
    w(56) = 0.0164711422416609467530350356258D+00
    w(57) = 0.0294992962483054353948393364098D+00
    w(58) = 0.0488473871144520262535428484316D+00
    w(59) = 0.074807989768816537216026182806D+00
    w(60) = 0.10598520508123912472195529192D+00
    w(61) = 0.138939453090947794093360848265D+00
    w(62) = 0.168562360742603870987330592834D+00
    w(63) = 0.189278495801793364889704841035D+00
    w(64) = 0.196733406888845140995323677102D+00
    w(65) = 0.189278495801793364889704841035D+00
    w(66) = 0.168562360742603870987330592834D+00
    w(67) = 0.138939453090947794093360848265D+00
    w(68) = 0.10598520508123912472195529192D+00
    w(69) = 0.074807989768816537216026182806D+00
    w(70) = 0.0488473871144520262535428484316D+00
    w(71) = 0.0294992962483054353948393364098D+00
    w(72) = 0.0164711422416609467530350356258D+00
    w(73) = 0.00850002630418086349941683729112D+00
    w(74) = 0.00405248551861722466559241860023D+00
    w(75) = 0.0017841208326818955520088211458D+00
    w(76) = 0.000724939297427239633212185817821D+00
    w(77) = 0.000271707626280157286781639661883D+00
    w(78) = 0.0000938744357203646866361259710004D+00
    w(79) = 0.0000298767905360019901790649251988D+00
    w(80) = 8.75244680345528247507614056972D-06
    w(81) = 2.35815617248490159838145978859D-06
    w(82) = 5.83809442762947462901022315301D-07
    w(83) = 1.32678558425807549298485884004D-07
    w(84) = 2.76504970450371674155194812923D-08
    w(85) = 5.27816639371369729333040255118D-09
    w(86) = 9.21739409677215086782446989876D-10
    w(87) = 1.47063312734774830028408333227D-10
    w(88) = 2.14069202905212884993201956606D-11
    w(89) = 2.83853037250817094975750489262D-12
    w(90) = 3.42300944935037851188976963928D-13
    w(91) = 3.74739544729563577089986076081D-14
    w(92) = 3.71733031252663248624409938613D-15
    w(93) = 3.33444143033026256341061235315D-16
    w(94) = 2.69865119072980851232572568063D-17
    w(95) = 1.96593268885070384943390296306D-18
    w(96) = 1.28579977867628696999762170542D-19
    w(97) = 7.52931616388155067444192947319D-21
    w(98) = 3.93559908609832906838466602268D-22
    w(99) = 1.83031346363374264415878982576D-23
    w(100) = 7.54682189033203465872349657723D-25
    w(101) = 2.74824892121260880467531987939D-26
    w(102) = 8.8019567691972403392314252914D-28
    w(103) = 2.4679773241854004762148469348D-29
    w(104) = 6.02777538509463291699314327193D-31
    w(105) = 1.27533314110484056196532640642D-32
    w(106) = 2.3232083386415854084664074623D-34
    w(107) = 3.61918834460649868835433546523D-36
    w(108) = 4.78543906802804099967221020647D-38
    w(109) = 5.32553037755907921458489847863D-40
    w(110) = 4.9411031115925407477456893331D-42
    w(111) = 3.78144279409152203964384313149D-44
    w(112) = 2.35816591560823143778744566357D-46
    w(113) = 1.18170544407210392716367665268D-48
    w(114) = 4.68218083833618292793410025836D-51
    w(115) = 1.4394641949298720336141068619D-53
    w(116) = 3.3581166223962736386929935773D-56
    w(117) = 5.78845633750656959788340019085D-59
    w(118) = 7.13750929465743002965122123123D-62
    w(119) = 6.04865489642049179016214753843D-65
    w(120) = 3.34891390118992716444169809114D-68
    w(121) = 1.13419242086298913875376620343D-71
    w(122) = 2.15037147336771602203551878273D-75
    w(123) = 2.01439576527109443920782513994D-79
    w(124) = 7.7306185241134158744827181222D-84
    w(125) = 8.93216815722645216635320162557D-89
    w(126) = 1.72727980594728851329952877284D-94
    w(127) = 1.25044975770895101066558695394D-101

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_LOOKUP_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1, 3, 7, 15, 31, 63 and 127.'
    stop

  end if

  return
end
subroutine hermite_ss_compute ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_SS_COMPUTE computes a Hermite quadrature rule.
!
!  Discussion:
!
!    The code uses an algorithm by Stroud and Secrest.
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!    The integration interval is ( -oo, +oo ).
!
!    The weight function is w(x) = exp ( - x^2 ).
!
!    The integral to approximate:
!
!      Integral ( -oo < X < +oo ) exp ( - X^2 ) * F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) s
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_SS_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', n
    stop
  end if

  cc = 1.7724538509D+00 * r8_gamma ( real ( n, kind = 8 ) ) &
    / ( 2.0D+00**( n - 1 ) )

  s = ( 2.0D+00 * real ( n, kind = 8 ) + 1.0D+00 )**( 1.0D+00 / 6.0D+00 )

  do i = 1, ( n + 1 ) / 2

    if ( i == 1 ) then

      x0 = s * s * s - 1.85575D+00 / s

    else if ( i == 2 ) then

      x0 = x0 - 1.14D+00 * ( ( real ( n, kind = 8 ) )**0.426D+00 ) / x0

    else if ( i == 3 ) then

      x0 = 1.86D+00 * x0 - 0.86D+00 * x(1)

    else if ( i == 4 ) then

      x0 = 1.91D+00 * x0 - 0.91D+00 * x(2)

    else

      x0 = 2.0D+00 * x0 - x(i-2)

    end if

    call hermite_ss_root ( x0, n, dp2, p1 )

    x(i) = x0
    w(i) = ( cc / dp2 ) / p1

    x(n-i+1) = - x0
    w(n-i+1) = w(i)

  end do
!
!  Reverse the order of the abscissas.
!  Because of symmetry, the weights are unchanged,
!  and the abscissas simply change sign.
!
  x(1:n) = - x(1:n)

  return
end
subroutine hermite_ss_recur ( p2, dp2, p1, x, order )

!*****************************************************************************80
!
!! HERMITE_SS_RECUR finds the value and derivative of a Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of H(ORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of H'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of H(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) order
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2  = x * p1 - 0.5D+00 * ( real ( i, kind = 8 ) - 1.0D+00 ) * p0
    dp2 = x * dp1 + p1 - 0.5D+00 * ( real ( i, kind = 8 ) - 1.0D+00 ) * dp0

  end do

  return
end
subroutine hermite_ss_root ( x, order, dp2, p1 )

!*****************************************************************************80
!
!! HERMITE_SS_ROOT improves an approximate root of a Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Output, real ( kind = 8 ) DP2, the value of H'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of H(ORDER-1)(X).
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  integer ( kind = 4 ) order
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = r8_epsilon ( )

  do step = 1, step_max

    call hermite_ss_recur ( p2, dp2, p1, x, order )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = factorial ( N ) / ( factorial ( K ) * factorial ( N - K ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
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
!    Output, integer ( kind = 4 ) I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
!
!  Discussion:
!
!    For positive I4_LOG_2(I), it should be true that
!      2^I4_LOG_2(X) <= |I| < 2^(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 2
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_2, the integer part of the
!    logarithm base 2 of the absolute value of I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647

  if ( i == 0 ) then

    i4_log_2 = - i4_huge

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  character ( len = * )  title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
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
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 10
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  character ( len = 8 )  ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_WRITE writes an I4MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  integer   ( kind = 4 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine i4vec_min_mv ( m, n, u, v, w )

!*****************************************************************************80
!
!! I4VEC_MIN_MV determines U(1:N) /\ V for vectors U and a single vector V.
!
!  Discussion:
!
!    For two vectors U and V, each of length M, we define
!
!      ( U /\ V ) (I) = min ( U(I), V(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the vectors.
!
!    Input, integer ( kind = 4 ) N, the number of vectors in U.
!
!    Input, integer ( kind = 4 ) U(M,N), N vectors, each of length M.
!
!    Input, integer ( kind = 4 ) V(M), a vector of length M.
!
!    Output, integer ( kind = 4 ) W(M,N), the value of U /\ W.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) u(m,n)
  integer ( kind = 4 ) v(m)
  integer ( kind = 4 ) w(m,n)

  do j = 1, n
    do i = 1, m
      w(i,j) = min ( u(i,j), v(i) )
    end do
  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

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
subroutine jacobi_compute ( n, alpha, beta, x, w )

!*****************************************************************************80
!
!! JACOBI_COMPUTE: Elhay-Kautsky method for Gauss-Jacobi quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) (1-x)^alpha * (1+x)^beta * f(x) dx
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
!    30 April 2011
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
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
!    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) aj(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) abi
  real ( kind = 8 ) beta
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
  zemu = 2.0D+00**( alpha + beta + 1.0D+00 ) &
    * r8_gamma ( alpha + 1.0D+00 ) &
    * r8_gamma ( beta + 1.0D+00 ) &
    / r8_gamma ( 2.0D+00 + alpha + beta )
!
!  Define the Jacobi matrix.
!
  x(1) = ( beta - alpha ) / ( 2.0D+00 + alpha + beta )

  bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
    / ( ( 3.0D+00 + alpha + beta ) * ( 2.0D+00 + alpha + beta )**2 )

  do i = 2, n
    i_r8 = real ( i, kind = 8 )
    abi = 2.0D+00 * i_r8 + alpha + beta
    x(i) = ( beta + alpha ) * ( beta - alpha ) / ( ( abi - 2.0D+00 ) * abi )
    bj(i) = 4.0D+00 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) &
      * ( i_r8 + alpha + beta ) &
      / ( ( abi - 1.0D+00 ) * ( abi + 1.0D+00 ) * abi * abi )
  end do

  bj(1:n) = sqrt ( bj(1:n) )

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine jacobi_compute_points ( order, alpha, beta, points )

!*****************************************************************************80
!
!! JACOBI_COMPUTE_POINTS returns the points of a Jacobi rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of the (1-X) and
!    (1+X) factors in the weight.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  call jacobi_compute ( order, alpha, beta, points, weight )

  return
end
subroutine jacobi_compute_points_np ( order, np, p, points )

!*****************************************************************************80
!
!! JACOBI_COMPUTE_POINTS_NP returns the points of a Jacobi rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) p(*)
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  alpha = p(1)
  beta = p(2)

  call jacobi_compute ( order, alpha, beta, points, weight )

  return
end
subroutine jacobi_compute_weights ( order, alpha, beta, weight )

!*****************************************************************************80
!
!! JACOBI_COMPUTE_WEIGHTS returns the weights of a Jacobi rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of the (1-X) and
!    (1+X) factors in the weight.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  call jacobi_compute ( order, alpha, beta, points, weight )

  return
end
subroutine jacobi_compute_weights_np ( order, np, p, weight )

!*****************************************************************************80
!
!! JACOBI_COMPUTE_WEIGHTS_NP returns the weights of a Jacobi rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) p(*)
  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weight(order)

  alpha = p(1)
  beta = p(2)

  call jacobi_compute ( order, alpha, beta, points, weight )

  return
end
subroutine jacobi_integral ( expon, alpha, beta, value )

!*****************************************************************************80
!
!! JACOBI_INTEGRAL evaluates the integral of a monomial with Jacobi weight.
!
!  Discussion:
!
!    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x)^ALPHA (1+x)^BETA dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X) in the weight factor.
!
!    Input, real ( kind = 8 ) BETA, the exponent of (1+X) in the weight factor.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) s
  real ( kind = 8 ) value
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2

  c = real ( expon, kind = 8 )

  if ( mod ( expon, 2 ) == 0 ) then
    s = +1.0D+00
  else
    s = -1.0D+00
  end if

  arg1 = - alpha
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + beta + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

  arg1 = - beta
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + alpha + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value2 )

  value = r8_gamma ( 1.0D+00 + c ) * ( &
      s * r8_gamma ( 1.0D+00 + beta  ) * value1 &
    / r8_gamma ( 2.0D+00 + beta  + c ) &
    +     r8_gamma ( 1.0D+00 + alpha ) * value2 &
    / r8_gamma ( 2.0D+00 + alpha + c ) )

  return
end
subroutine jacobi_ss_compute ( order, alpha, beta, x, w )

!*****************************************************************************80
!
!! JACOBI_SS_COMPUTE computes a Jacobi quadrature rule.
!
!  Discussion:
!
!    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) (1-X)^ALPHA * (1+X)^BETA * F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!    Thanks to Xu Xiang of Fudan University for pointing out that
!    an earlier implementation of this routine was incorrect!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2007
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.  For simple Legendre quadrature,
!    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) an
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) beta
  real ( kind = 8 ) bn
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) cc
  real ( kind = 8 ) delta
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(order)
  real ( kind = 8 ) x0

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_SS_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if
!
!  Check ALPHA and BETA.
!
  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_SS_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < ALPHA is required.'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_SS_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < BETA is required.'
    stop
  end if
!
!  Set the recursion coefficients.
!
  do i = 1, order

    if ( alpha + beta == 0.0D+00 .or. beta - alpha == 0.0D+00 ) then

      b(i) = 0.0D+00

    else

      b(i) = ( alpha + beta ) * ( beta - alpha ) / &
            ( ( alpha + beta + real ( 2 * i, kind = 8 ) ) &
            * ( alpha + beta + real ( 2 * i - 2, kind = 8 ) ) )

    end if

    if ( i == 1 ) then

      c(i) = 0.0D+00

    else

      c(i) = 4.0D+00 * real ( i - 1, kind = 8 ) &
            * ( alpha + real ( i - 1, kind = 8 ) ) &
            * ( beta + real ( i - 1, kind = 8 ) ) &
            * ( alpha + beta + real ( i - 1, kind = 8 ) ) / &
            ( ( alpha + beta + real ( 2 * i - 1, kind = 8 ) ) &
            * ( alpha + beta + real ( 2 * i - 2, kind = 8 ) )**2 &
            * ( alpha + beta + real ( 2 * i - 3, kind = 8 ) ) )

    end if

  end do

  delta = r8_gamma ( alpha        + 1.0D+00 ) &
        * r8_gamma (         beta + 1.0D+00 ) &
        / r8_gamma ( alpha + beta + 2.0D+00 )

  cc = delta * 2.0D+00**( alpha + beta + 1.0D+00 ) * product ( c(2:order) )

  do i = 1, order

    if ( i == 1 ) then

      an = alpha / real ( order, kind = 8 )
      bn = beta / real ( order, kind = 8 )

      r1 = ( 1.0D+00 + alpha ) &
        * ( 2.78D+00 / ( 4.0D+00 + real ( order**2, kind = 8 ) ) &
        + 0.768D+00 * an / real ( order, kind = 8 ) )

      r2 = 1.0D+00 + 1.48D+00 * an + 0.96D+00 * bn &
        + 0.452D+00 * an**2 + 0.83D+00 * an * bn

      x0 = ( r2 - r1 ) / r2

    else if ( i == 2 ) then

      r1 = ( 4.1D+00 + alpha ) / &
        ( ( 1.0D+00 + alpha ) * ( 1.0D+00 + 0.156D+00 * alpha ) )

      r2 = 1.0D+00 + 0.06D+00 * ( real ( order, kind = 8 ) - 8.0D+00 ) * &
        ( 1.0D+00 + 0.12D+00 * alpha ) / real ( order, kind = 8 )

      r3 = 1.0D+00 + 0.012D+00 * beta * &
        ( 1.0D+00 + 0.25D+00 * abs ( alpha ) ) / real ( order, kind = 8 )

      x0 = x0 - r1 * r2 * r3 * ( 1.0D+00 - x0 )

    else if ( i == 3 ) then

      r1 = ( 1.67D+00 + 0.28D+00 * alpha ) / ( 1.0D+00 + 0.37D+00 * alpha )

      r2 = 1.0D+00 + 0.22D+00 * ( real ( order, kind = 8 ) - 8.0D+00 ) &
        / real ( order, kind = 8 )

      r3 = 1.0D+00 + 8.0D+00 * beta / &
        ( ( 6.28D+00 + beta ) * real ( order**2, kind = 8 ) )

      x0 = x0 - r1 * r2 * r3 * ( x(1) - x0 )

    else if ( i < order - 1 ) then

      x0 = 3.0D+00 * x(i-1) - 3.0D+00 * x(i-2) + x(i-3)

    else if ( i == order - 1 ) then

      r1 = ( 1.0D+00 + 0.235D+00 * beta ) / ( 0.766D+00 + 0.119D+00 * beta )

      r2 = 1.0D+00 / ( 1.0D+00 + 0.639D+00 &
        * ( real ( order, kind = 8 ) - 4.0D+00 ) &
        / ( 1.0D+00 + 0.71D+00 * ( real ( order, kind = 8 ) - 4.0D+00 ) ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 20.0D+00 * alpha / ( ( 7.5D+00 + alpha ) * &
        real ( order**2, kind = 8 ) ) )

      x0 = x0 + r1 * r2 * r3 * ( x0 - x(i-2) )

    else if ( i == order ) then

      r1 = ( 1.0D+00 + 0.37D+00 * beta ) / ( 1.67D+00 + 0.28D+00 * beta )

      r2 = 1.0D+00 / &
        ( 1.0D+00 + 0.22D+00 * ( real ( order, kind = 8 ) - 8.0D+00 ) &
        / real ( order, kind = 8 ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( order**2, kind = 8 ) ) )

      x0 = x0 + r1 * r2 * r3 * ( x0 - x(i-2) )

    end if

    call jacobi_ss_root ( x0, order, alpha, beta, dp2, p1, b, c )

    x(i) = x0
    w(i) = cc / ( dp2 * p1 )

  end do
!
!  Reverse the order of the data.
!
  x(1:order) = x(order:1:-1)
  w(1:order) = w(order:1:-1)

  return
end
subroutine jacobi_ss_recur ( p2, dp2, p1, x, order, alpha, beta, b, c )

!*****************************************************************************80
!
!! JACOBI_SS_RECUR finds the value and derivative of a Jacobi polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of J(ORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of J'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.
!
!    Input, real ( kind = 8 ) B(ORDER), C(ORDER), the recursion
!    coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0D+00 )
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine jacobi_ss_root ( x, order, alpha, beta, dp2, p1, b, c )

!*****************************************************************************80
!
!! JACOBI_SS_ROOT improves an approximate root of a Jacobi polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2000
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.
!
!    Output, real ( kind = 8 ) DP2, the value of J'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) B(ORDER), C(ORDER), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = r8_epsilon ( )

  do step = 1, step_max

    call jacobi_ss_recur ( p2, dp2, p1, x, order, alpha, beta, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine laguerre_compute ( n, x, w )

!*****************************************************************************80
!
!! LAGUERRE_COMPUTE: Laguerre quadrature rule by the Elhay-Kautsky method.
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
subroutine laguerre_compute_points ( n, x )

!*****************************************************************************80
!
!! LAGUERRE_COMPUTE_POINTS computes points of a Laguerre quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2009
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call laguerre_compute ( n, x, w )

  return
end
subroutine laguerre_compute_points_np ( n, np, p, x )

!*****************************************************************************80
!
!! LAGUERRE_COMPUTE_POINTS_NP computes points of a Laguerre quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) n

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call laguerre_compute ( n, x, w )

  return
end
subroutine laguerre_compute_weights ( n, w )

!*****************************************************************************80
!
!! LAGUERRE_COMPUTE_WEIGHTS computes weights of a Laguerre quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2009
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) WN), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call laguerre_compute ( n, x, w )

  return
end
subroutine laguerre_compute_weights_np ( n, np, p, w )

!*****************************************************************************80
!
!! LAGUERRE_COMPUTE_WEIGHTS_NP computes weights of a Laguerre quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) n

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call laguerre_compute ( n, x, w )

  return
end
subroutine laguerre_integral ( expon, exact )

!*****************************************************************************80
!
!! LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.
!
!  Discussion:
!
!    The integral being computed is
!
!      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_factorial

  exact = r8_factorial ( expon )

  return
end
subroutine laguerre_lookup_points ( n, points )

!*****************************************************************************80
!
!! LAGUERRE_LOOKUP_POINTS returns the abscissas of a Laguerre rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [0,+oo).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 7, 15, 31, 63 and 127.
!
!    Output, real ( kind = 8 ) POINTS(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) points(n)
  real ( kind = 8 ), save, dimension ( 1 ) :: x_001 = (/ &
    1.0D+00 /)
  real ( kind = 8 ), save, dimension ( 3 ) :: x_003 = (/ &
    0.415774556783479083311533873128D+00, &
    0.229428036027904171982205036136D+01, &
    0.628994508293747919686641576551D+01 /)
  real ( kind = 8 ), save, dimension ( 7 ) :: x_007 = (/ &
    0.193043676560362413838247885004D+00, &
    0.102666489533919195034519944317D+01, &
    0.256787674495074620690778622666D+01, &
    0.490035308452648456810171437810D+01, &
    0.818215344456286079108182755123D+01, &
    0.127341802917978137580126424582D+02, &
    0.193957278622625403117125820576D+02 /)
  real ( kind = 8 ), save, dimension ( 15 ) :: x_015 = (/ &
    0.933078120172818047629030383672D-01, &
    0.492691740301883908960101791412D+00, &
    0.121559541207094946372992716488D+01, &
    0.226994952620374320247421741375D+01, &
    0.366762272175143727724905959436D+01, &
    0.542533662741355316534358132596D+01, &
    0.756591622661306786049739555812D+01, &
    0.101202285680191127347927394568D+02, &
    0.131302824821757235640991204176D+02, &
    0.166544077083299578225202408430D+02, &
    0.207764788994487667729157175676D+02, &
    0.256238942267287801445868285977D+02, &
    0.314075191697539385152432196202D+02, &
    0.385306833064860094162515167595D+02, &
    0.480260855726857943465734308508D+02 /)
  real ( kind = 8 ), save, dimension ( 31 ) :: x_031 = (/ &
    0.45901947621108290743496080275224D-01, &
    0.24198016382477204890408974151714D+00, &
    0.59525389422235073707330165005414D+00, &
    1.1066894995329987162111308789792D+00, &
    1.7775956928747727211593727482675D+00, &
    2.6097034152566806503893375925315D+00, &
    3.6051968023400442698805817554243D+00, &
    4.7667470844717611313629127271123D+00, &
    6.0975545671817409269925429328463D+00, &
    7.6014009492331374229360106942867D+00, &
    9.2827143134708894182536695297710D+00, &
    11.146649755619291358993815629587D+00, &
    13.199189576244998522464925028637D+00, &
    15.447268315549310075809325891801D+00, &
    17.898929826644757646725793817752D+00, &
    20.563526336715822170743048968779D+00, &
    23.451973482011858591050255575933D+00, &
    26.577081352118260459975876986478D+00, &
    29.953990872346445506951917840024D+00, &
    33.600759532902202735410313885784D+00, &
    37.539164407330440882887902558001D+00, &
    41.795830870182219981347945853330D+00, &
    46.403866806411123136029227604386D+00, &
    51.405314476797755161861461088395D+00, &
    56.854992868715843620511922055660D+00, &
    62.826855908786321453677523304806D+00, &
    69.425277191080345623322251656443D+00, &
    76.807047763862732837609972285484D+00, &
    85.230358607545669169387065607043D+00, &
    95.188939891525629981308606853957D+00, &
    107.95224382757871475002440117666D+00 /)
  real ( kind = 8 ), save, dimension ( 63 ) :: x_063 = (/ &
    0.22768893732576153785994330248562D-01, &
    0.11998325242727824715771416426383D+00, &
    0.29494185444770149577427738517405D+00, &
    0.54779087896237725363865073775856D+00, &
    0.87869061179931901673895567052285D+00, &
    1.2878464335919706302309207788611D+00, &
    1.7755123815388553763979463268728D+00, &
    2.3419925567085989256055628337716D+00, &
    2.9876423223246473939976731053629D+00, &
    3.7128695992018000346299637413422D+00, &
    4.5181363349503584391105568561550D+00, &
    5.4039601781825946286902599782736D+00, &
    6.3709163787865330220392250891777D+00, &
    7.4196399339311711154888493199004D+00, &
    8.5508280008403328312589048722235D+00, &
    9.7652425999245366807004592977996D+00, &
    11.063713635140661736220550410604D+00, &
    12.447142262356492749798687569289D+00, &
    13.916504641057818562912967008183D+00, &
    15.472856110036296424777143607779D+00, &
    17.117335833863588753116900303886D+00, &
    18.851171974154856850873483787506D+00, &
    20.675687448056515660377265667433D+00, &
    22.592306346311528381292277759986D+00, &
    24.602561094972638883700642760037D+00, &
    26.708100458737343969779087998829D+00, &
    28.910698500451382640177718103234D+00, &
    31.212264631175912885477773820802D+00, &
    33.614854909101154836598842888345D+00, &
    36.120684774484823056306328740825D+00, &
    38.732143442933582145626041607663D+00, &
    41.451810222318741191114726181363D+00, &
    44.282473071479233839358857134636D+00, &
    47.227149784295686898935095231536D+00, &
    50.289112264240695761749021839419D+00, &
    53.471914456788652808348280619542D+00, &
    56.779424636342062213099781057119D+00, &
    60.215862909019862886417550114424D+00, &
    63.785845004235974631701139601836D+00, &
    67.494433702293885830374325695045D+00, &
    71.347199604295266286654803376075D+00, &
    75.350293425653234254290504744279D+00, &
    79.510532629986309149555391354778D+00, &
    83.835506080872257843339817658508D+00, &
    88.333701570354369086112766326498D+00, &
    93.014662728558547405303399037100D+00, &
    97.889184147578140043386727677112D+00, &
    102.96955690741381650783952746778D+00, &
    108.26988161961595392226350967206D+00, &
    113.80647350287462738934485955901D+00, &
    119.59839538830458666962452963285D+00, &
    125.66817255856119431291196303280D+00, &
    132.04277272091165746585590583045D+00, &
    138.75498418103789078167590567526D+00, &
    145.84541318313540358283994248439D+00, &
    153.36548459497863623710815962660D+00, &
    161.38215194813761243562172669592D+00, &
    169.98570600665839438795175301156D+00, &
    179.30366247401580910251827858515D+00, &
    189.52789596532475473668721332981D+00, &
    200.97521159924656741628671841018D+00, &
    214.25368536638788642698056296400D+00, &
    230.93465747089703971246562985079D+00 /)
  real ( kind = 8 ), save, dimension ( 127 ) :: x_127 = (/ &
    0.11339635298518611691893169631306D-01, &
    0.59749753435726620281348237057387D-01, &
    0.14685098690746167612388223687431D+00, &
    0.27267590735859553131378008278900D+00, &
    0.43724600644192665554577035869932D+00, &
    0.64058688222566929533576416399983D+00, &
    0.88272968639058364481487653650042D+00, &
    1.1637114160166537661560584700951D+00, &
    1.4835750152834613891313584861012D+00, &
    1.8423694351613565380686320809853D+00, &
    2.2401496839579024244513315656522D+00, &
    2.6769768780141303692167869961238D+00, &
    3.1529182957082825565771508308846D+00, &
    3.6680474360304752540226339926515D+00, &
    4.2224440823301888455977876667425D+00, &
    4.8161943715870502475665535087286D+00, &
    5.4493908694559416755862178908416D+00, &
    6.1221326512997254193944584763155D+00, &
    6.8345253894122668112237994973336D+00, &
    7.5866814466367472174205986836847D+00, &
    8.3787199765932725254842120659452D+00, &
    9.2107670307426558777922506102445D+00, &
    10.082955672528643809166439353647D+00, &
    10.995426098858125429803147358780D+00, &
    11.948325769197725997610605127857D+00, &
    12.941809542585531053723381098192D+00, &
    13.976039822878506520014405668679D+00, &
    15.051186712579523631574796365435D+00, &
    16.167428175612852922977395051768D+00, &
    17.324950209443673446561163712616D+00, &
    18.523947026965688560811711309349D+00, &
    19.764621248611504104071669386884D+00, &
    21.047184105173183606877044020054D+00, &
    22.371855651855542817648123918101D+00, &
    23.738864994122497183652313788712D+00, &
    25.148450525937368234077278385644D+00, &
    26.600860181041749607253384279755D+00, &
    28.096351697964619201753961292129D+00, &
    29.635192899504178910610227138642D+00, &
    31.217661987479759144214467152615D+00, &
    32.844047853610430460522951341338D+00, &
    34.514650407441149149105635947422D+00, &
    36.229780922306804019615388508885D+00, &
    37.989762400399956435968780140278D+00, &
    39.794929958089961778396437141707D+00, &
    41.645631232730180705153990897484D+00, &
    43.542226812286859549950892993822D+00, &
    45.485090689228791137996151336673D+00, &
    47.474610740231964719468766599146D+00, &
    49.511189233379087716728884584381D+00, &
    51.595243364671244443182771266934D+00, &
    53.727205825819316758288140069145D+00, &
    55.907525405447553305830605991732D+00, &
    58.136667626022439197077526025660D+00, &
    60.415115419018590295707192053805D+00, &
    62.743369841051809700207126742685D+00, &
    65.121950833949996311956025417139D+00, &
    67.551398031997886314411872443149D+00, &
    70.032271619884584511229871192030D+00, &
    72.565153245206849090888669416801D+00, &
    75.150646989739935299354362325096D+00, &
    77.789380404085816000647405462136D+00, &
    80.482005610750729205803962926758D+00, &
    83.229200481195914886796120019048D+00, &
    86.031669892953582966798238732643D+00, &
    88.890147073512051099652518544282D+00, &
    91.805395038358177994971250170499D+00, &
    94.778208131331583205387031034825D+00, &
    97.809413676305116411054110115424D+00, &
    100.89987375017285940371939762172D+00, &
    104.05048708821598934704076845022D+00, &
    107.26219113414600428423116401414D+00, &
    110.53596424851500530602771351277D+00, &
    113.87282809075839485348376187652D+00, &
    117.27385019192517774095477886379D+00, &
    120.74014673718880106173978002719D+00, &
    124.27288557955698354259506446928D+00, &
    127.87328950885942645093841745425D+00, &
    131.54263980314366921809377742137D+00, &
    135.28228009311836970132738106369D+00, &
    139.09362057432970013964422086977D+00, &
    142.97814260643601776808227753574D+00, &
    146.93740374437366549441080969072D+00, &
    150.97304325252187127492511437460D+00, &
    155.08678816034612572229641420609D+00, &
    159.28045992663288235401956989889D+00, &
    163.55598178957571104015967182053D+00, &
    167.91538689194360134245547184721D+00, &
    172.36082728473812536838156191681D+00, &
    176.89458392960192176311674993508D+00, &
    181.51907784036813069227528834025D+00, &
    186.23688252828112373861202530357D+00, &
    191.05073794450929196790836610789D+00, &
    195.96356614879879837839002542988D+00, &
    200.97848897600025153696475526130D+00, &
    206.09884802468871112127283042753D+00, &
    211.32822735671655260572377256981D+00, &
    216.67047937658230323477089465777D+00, &
    222.12975445929687246267304963754D+00, &
    227.71053502072232419089132431317D+00, &
    233.41767488282602453367775322563D+00, &
    239.25644498830308620018749667089D+00, &
    245.23258677871567172531254018984D+00, &
    251.35237488718128030005500991754D+00, &
    257.62269123792061413076191882313D+00, &
    264.05111322908240551754377241831D+00, &
    270.64601945722796749299111718606D+00, &
    277.41671750163651071798388218104D+00, &
    284.37359974220870326674402873120D+00, &
    291.52833521346495719581282021650D+00, &
    298.89410837028248600878895615414D+00, &
    306.48591978262611320418112423947D+00, &
    314.32096986471177487400007507615D+00, &
    322.41915589128679683349440361344D+00, &
    330.80372663802405651933847334878D+00, &
    339.50216127832433747735367595958D+00, &
    348.54737559472697355480761787441D+00, &
    357.97942028029845454049007443090D+00, &
    367.84794520076004578858341422871D+00, &
    378.21590623135532818332979188889D+00, &
    389.16539141251004101579475325153D+00, &
    400.80729331451702589996361286427D+00, &
    413.29853681779384418008260081859D+00, &
    426.87579153663675538288509017051D+00, &
    441.93085485310841412460309271842D+00, &
    459.21804639888429981971267313224D+00, &
    480.69378263388373859704269229304D+00 /)

  if ( n == 1 ) then

    points(1:n) = x_001(1:n)

  else if ( n == 3 ) then

    points(1:n) = x_003(1:n)

  else if ( n == 7 ) then

    points(1:n) = x_007(1:n)

  else if ( n == 15 ) then

    points(1:n) = x_015(1:n)

  else if ( n == 31 ) then

    points(1:n) = x_031(1:n)

  else if ( n == 63 ) then

    points(1:n) = x_063(1:n)

  else if ( n == 127 ) then

    points(1:n) = x_127(1:n)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_LOOKUP_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Unexpected value of N = ', n
    stop

  end if

  return
end
subroutine laguerre_lookup_weights ( order, weight )

!*****************************************************************************80
!
!! LAGUERRE_LOOKUP_WEIGHTS returns weights for Laguerre quadrature rules.
!
!  Discussion:
!
!    The allowed orders are 1, 3, 7, 15, 31, 63, and 127.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2007
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
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    Legal values are 1, 3, 7, 15, 31, 63 or 127.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric and should sum to 1.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) weight(order)

  if ( order == 1 ) then

    weight(1) = 1.0D+00

  else if ( order == 3 ) then

    weight(1) =  0.711093009929173015449590191143D+00
    weight(2) =  0.278517733569240848801444888457D+00
    weight(3) =  0.103892565015861357489649204007D-01

  else if ( order == 7 ) then

    weight(1) =  0.409318951701273902130432880018D+00
    weight(2) =  0.421831277861719779929281005417D+00
    weight(3) =  0.147126348657505278395374184637D+00
    weight(4) =  0.206335144687169398657056149642D-01
    weight(5) =  0.107401014328074552213195962843D-02
    weight(6) =  0.158654643485642012687326223234D-04
    weight(7) =  0.317031547899558056227132215385D-07

  else if ( order == 15 ) then

    weight(1) =  0.218234885940086889856413236448D+00
    weight(2) =  0.342210177922883329638948956807D+00
    weight(3) =  0.263027577941680097414812275022D+00
    weight(4) =  0.126425818105930535843030549378D+00
    weight(5) =  0.402068649210009148415854789871D-01
    weight(6) =  0.856387780361183836391575987649D-02
    weight(7) =  0.121243614721425207621920522467D-02
    weight(8) =  0.111674392344251941992578595518D-03
    weight(9) =  0.645992676202290092465319025312D-05
    weight(10) = 0.222631690709627263033182809179D-06
    weight(11) = 0.422743038497936500735127949331D-08
    weight(12) = 0.392189726704108929038460981949D-10
    weight(13) = 0.145651526407312640633273963455D-12
    weight(14) = 0.148302705111330133546164737187D-15
    weight(15) = 0.160059490621113323104997812370D-19

  else if ( order == 31 ) then

    weight(  1) =   0.11252789550372583820847728082801D+00
    weight(  2) =   0.21552760818089123795222505285045D+00
    weight(  3) =   0.23830825164569654731905788089234D+00
    weight(  4) =   0.19538830929790229249915303390711D+00
    weight(  5) =   0.12698283289306190143635272904602D+00
    weight(  6) =   0.67186168923899300670929441993508D-01
    weight(  7) =   0.29303224993879487404888669311974D-01
    weight(  8) =   0.10597569915295736089529380314433D-01
    weight(  9) =   0.31851272582386980320974842433019D-02
    weight( 10) =   0.79549548307940382922092149012477D-03
    weight( 11) =   0.16480052126636687317862967116412D-03
    weight( 12) =   0.28229237864310816393860971468993D-04
    weight( 13) =   0.39802902551008580387116174900106D-05
    weight( 14) =   0.45931839841801061673729694510289D-06
    weight( 15) =   0.43075545187731100930131457465897D-07
    weight( 16) =   0.32551249938271570855175749257884D-08
    weight( 17) =   0.19620246675410594996247151593142D-09
    weight( 18) =   0.93190499086617587129534716431331D-11
    weight( 19) =   0.34377541819411620520312597898311D-12
    weight( 20) =   0.96795247130446716997405035776206D-14
    weight( 21) =   0.20368066110115247398010624219291D-15
    weight( 22) =   0.31212687280713526831765358632585D-17
    weight( 23) =   0.33729581704161052453395678308350D-19
    weight( 24) =   0.24672796386616696011038363242541D-21
    weight( 25) =   0.11582201904525643634834564576593D-23
    weight( 26) =   0.32472922591425422434798022809020D-26
    weight( 27) =   0.49143017308057432740820076259666D-29
    weight( 28) =   0.34500071104808394132223135953806D-32
    weight( 29) =   0.87663710117162041472932760732881D-36
    weight( 30) =   0.50363643921161490411297172316582D-40
    weight( 31) =   0.19909984582531456482439549080330D-45

  else if ( order == 63 ) then

    weight(  1) =   0.57118633213868979811587283390476D-01
    weight(  2) =   0.12067476090640395283319932036351D+00
    weight(  3) =   0.15925001096581873723870561096472D+00
    weight(  4) =   0.16875178327560799234596192963585D+00
    weight(  5) =   0.15366641977668956696193711310131D+00
    weight(  6) =   0.12368770614716481641086652261948D+00
    weight(  7) =   0.89275098854848671545279150057422D-01
    weight(  8) =   0.58258485446105944957571825725160D-01
    weight(  9) =   0.34546657545992580874717085812508D-01
    weight( 10) =   0.18675685985714656798286552591203D-01
    weight( 11) =   0.92233449044093536528490075241649D-02
    weight( 12) =   0.41671250684839592762582663470209D-02
    weight( 13) =   0.17238120299900582715386728541955D-02
    weight( 14) =   0.65320845029716311169340559359043D-03
    weight( 15) =   0.22677644670909586952405173207471D-03
    weight( 16) =   0.72127674154810668410750270234861D-04
    weight( 17) =   0.21011261180466484598811536851241D-04
    weight( 18) =   0.56035500893357212749181536071292D-05
    weight( 19) =   0.13673642785604888017836641282292D-05
    weight( 20) =   0.30507263930195817240736097189550D-06
    weight( 21) =   0.62180061839309763559981775409241D-07
    weight( 22) =   0.11566529551931711260022448996296D-07
    weight( 23) =   0.19614588267565478081534781863335D-08
    weight( 24) =   0.30286171195709411244334756404054D-09
    weight( 25) =   0.42521344539400686769012963452599D-10
    weight( 26) =   0.54202220578073819334698791381873D-11
    weight( 27) =   0.62627306838597672554166850420603D-12
    weight( 28) =   0.65474443156573322992307089591924D-13
    weight( 29) =   0.61815575808729181846302500000047D-14
    weight( 30) =   0.52592721363507381404263991342633D-15
    weight( 31) =   0.40230920092646484015391506025408D-16
    weight( 32) =   0.27600740511819536505013824207729D-17
    weight( 33) =   0.16936946756968296053322009855265D-18
    weight( 34) =   0.92689146872177087314963772462726D-20
    weight( 35) =   0.45093739060365632939780140603959D-21
    weight( 36) =   0.19435162876132376573629962695374D-22
    weight( 37) =   0.73926270895169207037999639194513D-24
    weight( 38) =   0.24714364154434632615980126000066D-25
    weight( 39) =   0.72288649446741597655145390616476D-27
    weight( 40) =   0.18407617292614039362985209905608D-28
    weight( 41) =   0.40583498566841960105759537058880D-30
    weight( 42) =   0.77000496416438368114463925286343D-32
    weight( 43) =   0.12488505764999334328843314866038D-33
    weight( 44) =   0.17185000226767010697663950619912D-35
    weight( 45) =   0.19896372636672396938013975755522D-37
    weight( 46) =   0.19199671378804058267713164416870D-39
    weight( 47) =   0.15278588285522166920459714708240D-41
    weight( 48) =   0.99054752688842142955854138884590D-44
    weight( 49) =   0.51597523673029211884228858692990D-46
    weight( 50) =   0.21249846664084111245693912887783D-48
    weight( 51) =   0.67903852766852910591172042494884D-51
    weight( 52) =   0.16466654148296177467908300517887D-53
    weight( 53) =   0.29509065402691055027053659375033D-56
    weight( 54) =   0.37838420647571051984882241014675D-59
    weight( 55) =   0.33358130068542431878174667995217D-62
    weight( 56) =   0.19223461022273880981363303073329D-65
    weight( 57) =   0.67812696961083016872779388922288D-69
    weight( 58) =   0.13404752802440604607620468935693D-72
    weight( 59) =   0.13109745101805029757648048223928D-76
    weight( 60) =   0.52624863881401787388694579143866D-81
    weight( 61) =   0.63780013856587414257760666006511D-86
    weight( 62) =   0.12997078942372924566347473916943D-91
    weight( 63) =   0.10008511496968754063443740168421D-98

  else if ( order == 127 ) then

    weight(  1) =   0.28773246692000124355770010301506D-01
    weight(  2) =   0.63817468175134649363480949265236D-01
    weight(  3) =   0.91919669721570571389864194652717D-01
    weight(  4) =   0.11054167914413766381245463002967D+00
    weight(  5) =   0.11879771633375850188328329422643D+00
    weight(  6) =   0.11737818530052695148804451630074D+00
    weight(  7) =   0.10819305984180551488335145581193D+00
    weight(  8) =   0.93827075290489628080377261401107D-01
    weight(  9) =   0.76966450960588843995822485928431D-01
    weight( 10) =   0.59934903912939714332570730063476D-01
    weight( 11) =   0.44417742073889001371708316272923D-01
    weight( 12) =   0.31385080966252320983009372215062D-01
    weight( 13) =   0.21172316041924506411370709025015D-01
    weight( 14) =   0.13650145364230541652171185564626D-01
    weight( 15) =   0.84172852710599172279366657385445D-02
    weight( 16) =   0.49674990059882760515912858620175D-02
    weight( 17) =   0.28069903895001884631961957446400D-02
    weight( 18) =   0.15192951003941952460445341057817D-02
    weight( 19) =   0.78789028751796084086217287140548D-03
    weight( 20) =   0.39156751064868450584507324648999D-03
    weight( 21) =   0.18652434268825860550093566260060D-03
    weight( 22) =   0.85173160415576621908809828160247D-04
    weight( 23) =   0.37285639197853037712145321577724D-04
    weight( 24) =   0.15648416791712993947447805296768D-04
    weight( 25) =   0.62964340695224829035692735524979D-05
    weight( 26) =   0.24288929711328724574541379938222D-05
    weight( 27) =   0.89824607890051007201922871545035D-06
    weight( 28) =   0.31844174740760353710742966328091D-06
    weight( 29) =   0.10821272905566839211861807542741D-06
    weight( 30) =   0.35245076750635536015902779085340D-07
    weight( 31) =   0.11001224365719347407063839761738D-07
    weight( 32) =   0.32904079616717932125329343003261D-08
    weight( 33) =   0.94289145237889976419772700772988D-09
    weight( 34) =   0.25882578904668318184050195309296D-09
    weight( 35) =   0.68047437103370762630942259017560D-10
    weight( 36) =   0.17131398805120837835399564475632D-10
    weight( 37) =   0.41291744524052865469443922304935D-11
    weight( 38) =   0.95264189718807273220707664873469D-12
    weight( 39) =   0.21032604432442425932962942047474D-12
    weight( 40) =   0.44427151938729352860940434285789D-13
    weight( 41) =   0.89760500362833703323319846405449D-14
    weight( 42) =   0.17341511407769287074627948346848D-14
    weight( 43) =   0.32028099548988356631494379835210D-15
    weight( 44) =   0.56531388950793682022660742095189D-16
    weight( 45) =   0.95329672799026591234588044025896D-17
    weight( 46) =   0.15353453477310142565288509437552D-17
    weight( 47) =   0.23608962179467365686057842132176D-18
    weight( 48) =   0.34648742794456611332193876653230D-19
    weight( 49) =   0.48515241897086461320126957663545D-20
    weight( 50) =   0.64786228633519813428137373790678D-21
    weight( 51) =   0.82476020965403242936448553126316D-22
    weight( 52) =   0.10005361880214719793491658282977D-22
    weight( 53) =   0.11561395116207304954233181263632D-23
    weight( 54) =   0.12719342731167922655612134264961D-24
    weight( 55) =   0.13316584714165372967340004160814D-25
    weight( 56) =   0.13261218454678944033646108509198D-26
    weight( 57) =   0.12554995447643949807286074138324D-27
    weight( 58) =   0.11294412178579462703240913107219D-28
    weight( 59) =   0.96491020279562119228500608131696D-30
    weight( 60) =   0.78241846768302099396733076955632D-31
    weight( 61) =   0.60181503542219626658249939076636D-32
    weight( 62) =   0.43882482704961741551510518054138D-33
    weight( 63) =   0.30314137647517256304035802501863D-34
    weight( 64) =   0.19826016543944539545224676057020D-35
    weight( 65) =   0.12267623373665926559013654872402D-36
    weight( 66) =   0.71763931692508888943812834967620D-38
    weight( 67) =   0.39659378833836963584113716149270D-39
    weight( 68) =   0.20688970553868040099581951696677D-40
    weight( 69) =   0.10179587017979517245268418427523D-41
    weight( 70) =   0.47200827745986374625714293679649D-43
    weight( 71) =   0.20606828985553374825744353490744D-44
    weight( 72) =   0.84627575907305987245899032156188D-46
    weight( 73) =   0.32661123687088798658026998931647D-47
    weight( 74) =   0.11833939207883162380564134612682D-48
    weight( 75) =   0.40211209123895013807243250164050D-50
    weight( 76) =   0.12799824394111125389430292847476D-51
    weight( 77) =   0.38123877747548846504399051365162D-53
    weight( 78) =   0.10612057542701156767898551949650D-54
    weight( 79) =   0.27571446947200403594113572720812D-56
    weight( 80) =   0.66772544240928492881306904862856D-58
    weight( 81) =   0.15052438383868234954068178600268D-59
    weight( 82) =   0.31538986800113758526689068500772D-61
    weight( 83) =   0.61326614299483180785237418887960D-63
    weight( 84) =   0.11048510030324810567549119229368D-64
    weight( 85) =   0.18410563538091348076979665543900D-66
    weight( 86) =   0.28323926570052832195543883237652D-68
    weight( 87) =   0.40154409843763655508670978777418D-70
    weight( 88) =   0.52351530215683708779772201956106D-72
    weight( 89) =   0.62634476665005100555787696642851D-74
    weight( 90) =   0.68612210535666530365348093803922D-76
    weight( 91) =   0.68651298840956019297134099761855D-78
    weight( 92) =   0.62581388433728084867318704240915D-80
    weight( 93) =   0.51833271237514904046803469968027D-82
    weight( 94) =   0.38893621571918443533108973497673D-84
    weight( 95) =   0.26357711379476932781525533730623D-86
    weight( 96) =   0.16078851293917979699005509638883D-88
    weight( 97) =   0.87978042070968939637972577886624D-91
    weight( 98) =   0.43013405077495109903408697802188D-93
    weight( 99) =   0.18713435881342838527144321803729D-95
    weight(100) =   0.72125744708060471675805761366523D-98
    weight(101) =   0.24508746062177874383231742333023D-100
    weight(102) =   0.73042094619470875777647865078327D-103
    weight(103) =   0.18983290818383463537886818579820D-105
    weight(104) =   0.42757400244246684123093264825902D-108
    weight(105) =   0.82894681420515755691423485228897D-111
    weight(106) =   0.13729432219324400013067050156048D-113
    weight(107) =   0.19265464126404973222043166489406D-116
    weight(108) =   0.22693344503301354826140809941334D-119
    weight(109) =   0.22209290603717355061909071271535D-122
    weight(110) =   0.17851087685544512662856555121755D-125
    weight(111) =   0.11630931990387164467431190485525D-128
    weight(112) =   0.60524443584652392290952805077893D-132
    weight(113) =   0.24729569115063528647628375096400D-135
    weight(114) =   0.77789065006489410364997205809045D-139
    weight(115) =   0.18409738662712607039570678274636D-142
    weight(116) =   0.31900921131079114970179071968597D-146
    weight(117) =   0.39179487139174199737617666077555D-150
    weight(118) =   0.32782158394188697053774429820559D-154
    weight(119) =   0.17793590713138888062819640128739D-158
    weight(120) =   0.58882353408932623157467835381214D-163
    weight(121) =   0.10957236509071169877747203273886D-167
    weight(122) =   0.10281621114867000898285076975760D-172
    weight(123) =   0.41704725557697758145816510853967D-178
    weight(124) =   0.58002877720316101774638319601971D-184
    weight(125) =   0.18873507745825517106171619101120D-190
    weight(126) =   0.69106601826730911682786705950895D-198
    weight(127) =   0.43506813201105855628383313334402D-207

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_LOOKUP_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    write ( *, '(a)' ) '  Legal values are 1, 3, 7, 15, 31, 63 and 127.'
    stop

  end if

  return
end
subroutine laguerre_ss_compute ( order, x, w )

!*****************************************************************************80
!
!! LAGUERRE_SS_COMPUTE computes a Laguerre quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ 0, +oo ).
!
!    The weight function is w(x) = exp ( -x ).
!
!
!    If the integral to approximate is:
!
!      Integral ( 0 <= X < +oo ) exp ( - X ) * F(X) dX
!
!    then the quadrature rule is:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!
!    If the integral to approximate is:
!
!      Integral ( 0 <= X < +oo ) F(X) dX
!
!    then the quadrature rule is:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * exp ( X(I) ) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2009
!
!  Author:
!
!    Original FORTRAN77 by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(order)
  real ( kind = 8 ) x0

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_SS_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if
!
!  Set the recursion coefficients.
!
  do i = 1, order
    b(i) = ( real ( 2 * i - 1, kind = 8 ) )
  end do

  do i = 1, order
    c(i) = real ( i - 1, kind = 8 ) * ( real ( i - 1, kind = 8 ) )
  end do

  cc = product ( c(2:order) )

  do i = 1, order
!
!  Compute an estimate for the root.
!
    if ( i == 1 ) then

      x0 = 3.0D+00 / ( 1.0D+00 + 2.4D+00 * real ( order, kind = 8 ) )

    else if ( i == 2 ) then

      x0 = x0 + 15.0D+00 / ( 1.0D+00 + 2.5D+00 * real ( order, kind = 8 ) )

    else

      r1 = ( 1.0D+00 + 2.55D+00 * real ( i - 2, kind = 8 ) ) &
        / ( 1.9D+00 * real ( i - 2, kind = 8 ) )

      x0 = x0 + r1 * ( x0 - x(i-2) )

    end if
!
!  Use iteration to find the root.
!
    call laguerre_ss_root ( x0, order, dp2, p1, b, c )
!
!  Set the abscissa and weight.
!
    x(i) = x0
!
!  Because of the huge values involved, this calculation breaks down
!  for ORDER = 127.
!
!  It was originally w(i) = cc / dp2 / p1, which breaks down sooner.
!
    w(i) = ( 1.0D+00 / dp2 )
    do j = 2, order
      w(i) = w(i) * real ( j - 1, kind = 8 )
    end do
    w(i) = w(i) / p1
    do j = 2, order
      w(i) = w(i) * real ( j - 1, kind = 8 )
    end do

!   w(i) = ( cc / dp2 ) / p1

  end do

  return
end
subroutine laguerre_ss_recur ( p2, dp2, p1, x, order, b, c )

!*****************************************************************************80
!
!! LAGUERRE_SS_RECUR finds the value and derivative of a Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of L(ORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of L'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Input, real ( kind = 8 ) B(ORDER), C(ORDER), the recursion
!    coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x - 1.0D+00
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine laguerre_ss_root ( x, order, dp2, p1, b, c )

!*****************************************************************************80
!
!! LAGUERRE_SS_ROOT improves an approximate root of a Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial.
!
!    Output, real ( kind = 8 ) DP2, the value of L'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) B(ORDER), C(ORDER), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = r8_epsilon ( )

  do step = 1, step_max

    call laguerre_ss_recur ( p2, dp2, p1, x, order, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine legendre_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE: Legendre quadrature rule by the Elhay-Kautsky method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
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
  zemu = 2.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i * i, kind = 8 ) / real ( 4 * i * i - 1, kind = 8 )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine legendre_compute_points ( order, x )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_POINTS computes abscissas of a Legendre quadrature rule.
!
!  Discussion:
!
!    This is just a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2009
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) x(order)
  real ( kind = 8 ) w(order)

  call legendre_compute ( order, x, w )

  return
end
subroutine legendre_compute_points_np ( order, np, p, x )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_POINTS_NP computes abscissas of a Legendre quadrature rule.
!
!  Discussion:
!
!    This is just a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) x(order)
  real ( kind = 8 ) w(order)

  call legendre_compute ( order, x, w )

  return
end
subroutine legendre_compute_weights ( order, w )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_WEIGHTS computes weights of a Legendre quadrature rule.
!
!  Discussion:
!
!    This is just a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2009.
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) x(order)
  real ( kind = 8 ) w(order)

  call legendre_compute ( order, x, w )

  return
end
subroutine legendre_compute_weights_np ( order, np, p, w )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_WEIGHTS_NP computes weights of a Legendre quadrature rule.
!
!  Discussion:
!
!    This is just a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2009
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) x(order)
  real ( kind = 8 ) w(order)

  call legendre_compute ( order, x, w )

  return
end
subroutine legendre_dr_compute ( order, x, w )

!*****************************************************************************80
!
!! LEGENDRE_DR_COMPUTE computes a Legendre quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2pn
  real ( kind = 8 ) d3pn
  real ( kind = 8 ) d4pn
  real ( kind = 8 ) dp
  real ( kind = 8 ) dpn
  real ( kind = 8 ) e1
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iback
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1mi
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) nmove
  real ( kind = 8 ) p
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pk
  real ( kind = 8 ) pkm1
  real ( kind = 8 ) pkp1
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x(order)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xtemp
  real ( kind = 8 ) w(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_DR_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if

  e1 = real ( order * ( order + 1 ), kind = 8 )

  m = ( order + 1 ) / 2

  do i = 1, m

    mp1mi = m + 1 - i

    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * order + 2, kind = 8 )

    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( order, kind = 8 ) ) &
      / real ( 8 * order * order, kind = 8 ) )

    pkm1 = 1.0D+00
    pk = x0

    do k = 2, order
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do

    d1 = real ( order, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 * x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 * x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) / &
      ( 1.0D+00 - x0 * x0 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp

    xtemp = x0 + h

    x(mp1mi) = xtemp

    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

    w(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )

  end do

  if ( mod ( order, 2 ) == 1 ) then
    x(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = ( order + 1 ) / 2
  ncopy = order - nmove

  do i = 1, nmove
    iback = order + 1 - i
    x(iback) = x(iback-ncopy)
    w(iback) = w(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, order - nmove
    x(i) = - x(order+1-i)
    w(i) = w(order+1-i)
  end do

  return
end
subroutine legendre_integral ( expon, exact )

!*****************************************************************************80
!
!! LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.
!
!  Discussion:
!
!    To test a Legendre quadrature rule, we use it to approximate the
!    integral of a monomial:
!
!      integral ( -1 <= x <= +1 ) x^n dx
!
!    This routine is given the value of the exponent, and returns the
!    exact value of the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Output, real ( kind = 8 ) EXACT, the value of the exact integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
!
!  Get the exact value of the integral.
!
  if ( mod ( expon, 2 ) == 0 ) then

    exact = 2.0D+00 / real ( expon + 1, kind = 8 )

  else

    exact = 0.0D+00

  end if

  return
end
subroutine legendre_zeros ( n, x )

!*****************************************************************************80
!
!! LEGENDRE_ZEROS computes the zeros of the Legendre polynomial of degree N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2011
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    0 < N.
!
!    Output, real ( kind = 8 ) X(N), the locations of the zeros.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2pn
  real ( kind = 8 ) d3pn
  real ( kind = 8 ) d4pn
  real ( kind = 8 ) dp
  real ( kind = 8 ) dpn
  real ( kind = 8 ) e1
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iback
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1mi
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) nmove
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pk
  real ( kind = 8 ) pkm1
  real ( kind = 8 ) pkp1
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xtemp

  e1 = real ( n * ( n + 1 ), kind = 8 )

  m = ( n + 1 ) / 2

  do i = 1, m

    mp1mi = m + 1 - i

    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * n + 2, kind = 8 )

    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( n, kind = 8 ) ) &
      / real ( 8 * n * n, kind = 8 ) )

    pkm1 = 1.0D+00
    pk = x0

    do k = 2, n
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do

    d1 = real ( n, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 ) &
      / ( 1.0D+00 + x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) &
      / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp

    xtemp = x0 + h

    x(mp1mi) = xtemp

    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

  end do

  if ( mod ( n, 2 ) == 1 ) then
    x(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = ( n + 1 ) / 2
  ncopy = n - nmove

  do i = 1, nmove
    iback = n + 1 - i
    x(iback) = x(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, n - nmove
    x(i) = - x(n+1-i)
  end do

  return
end
subroutine level_growth_to_order ( dim_num, level, rule, growth, order )

!*****************************************************************************80
!
!! LEVEL_GROWTH_TO_ORDER: convert Level and Growth to Order.
!
!  Discussion:
!
!    This function is given level, rule, and growth information
!    for each dimension of a quadrature rule, and determines the
!    corresponding order of the rule in each dimension.
!
!    This is a revised version of the LEVEL_GROWTH_TO_ORDER function.
!
!    In particular, it revises the interpretation of the RULE vector as 
!    far as the values 10, 11, and 12 are concerned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL(DIM_NUM), the 1D levels.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
!    11, "UO",  User supplied Open, presumably Non Nested.
!    12, "UC",  User supplied Closed, presumably Non Nested.
!
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the desired growth in
!    each dimension.
!    0, "DF", default growth associated with this quadrature rule
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) o
!
!  The order of the final HGK rule can be one of 35, 37, 41 or 43.
!
  integer ( kind = 4 ) :: o_hgk(0:4) = (/ 1, 3, 9, 19, 43 /)
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) p
!
!  The precision of the final HGK rule must correspondingly be 51, 55, 63 or 67.
!
  integer ( kind = 4 ) :: p_hgk(0:4) = (/ 1, 5, 15, 29, 67 /)
  integer ( kind = 4 ) rule(dim_num)
!
!  Check the input.
!
  do dim = 1, dim_num

    if ( level(dim) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
      write ( *, '(a)' ) '  Negative value of LEVEL(DIM)!'
      write ( *, '(a,i8,a,i8)' ) '  LEVEL(', dim, ') = ', level(dim)
      stop
    end if

    if ( rule(dim) < 1 .or. 12 < rule(dim) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of RULE(DIM)!'
      write ( *, '(a,i8,a,i8)' ) '  RULE(', dim, ') = ', rule(dim)
      stop
    end if

    if ( growth(dim) < 0 .or. 6 < growth(dim) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of GROWTH(DIM)!'
      write ( *, '(a,i8,a,i8)' ) '  GROWTH(', dim, ') = ', growth(dim)
      stop
    end if

  end do
!
!  Compute the order vector.
!
  do dim = 1, dim_num
!
!  CC
!  Default is Moderate Exponential Growth.
!
    if ( rule(dim) == 1 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        if ( level(dim) == 0 ) then
          o = 1
        else
          o = 2
          do while ( o < 2 * level(dim) + 1 )
            o = 2 * ( o - 1 ) + 1
          end do
        end if
      else if ( growth(dim) == 5 .or. growth(dim) == 0 ) then
        if ( level(dim) == 0 ) then
          o = 1
        else
          o = 2
          do while ( o < 4 * level(dim) + 1 )
            o = 2 * ( o - 1 ) + 1
          end do
        end if
      else if ( growth(dim) == 6 ) then
        if ( level(dim) == 0 ) then
          o = 1
        else
          o = 2 ** level(dim) + 1
        end if
      end if
!
!  F2
!  Default is Moderate Exponential Growth.
!
    else if ( rule(dim) == 2 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        o = 1
        do while ( o < 2 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 5 .or. growth(dim) == 0 ) then
        o = 1
        do while ( o < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  GP
!  Default is Moderate Exponential Growth.
!
    else if ( rule(dim) == 3 ) then
      if ( growth(dim) == 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Growth rate 1 for rule 3 not available.'
        stop
      else if ( growth(dim) == 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Growth rate 2 for rule 3 not available.'
        stop
      else if ( growth(dim) == 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Growth rate 3 for rule 3 not available.'
        stop
      else if ( growth(dim) == 4 ) then
        if ( level(dim) == 0 ) then
          o = 1
        else
          p = 5
          o = 3
          do while ( p < 2 * level(dim) + 1 )
            p = 2 * p + 1
            o = 2 * o + 1
          end do
        end if
      else if ( growth(dim) == 5 .or. growth(dim) == 0 ) then
        if ( level(dim) == 0 ) then
          o = 1
        else
          p = 5
          o = 3
          do while ( p < 4 * level(dim) + 1 )
            p = 2 * p + 1
            o = 2 * o + 1
          end do
        end if
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  GL
!  Default is Moderate Linear Growth.
!
    else if ( rule(dim) == 4 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 .or. growth(dim) == 0 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        o = 1
        do while ( 2 * o - 1 < 2 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 5 ) then
        o = 1
        do while ( 2 * o - 1 < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  GH
!  Default is Moderate Linear Growth.
!
    else if ( rule(dim) == 5 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 .or. growth(dim) == 0 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        o = 1
        do while ( 2 * o - 1 < 2 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 5 ) then
        o = 1
        do while ( 2 * o - 1 < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  GGH
!  Default is Moderate Linear Growth.
!
    else if ( rule(dim) == 6 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 .or. growth(dim) == 0 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        o = 1
        do while ( 2 * o - 1 < 2 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 5 ) then
        o = 1
        do while ( 2 * o - 1 < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  LG
!  Default is Moderate Linear Growth.
!
    else if ( rule(dim) == 7 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 .or. growth(dim) == 0 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        o = 1
        do while ( 2 * o - 1 < 2 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 5 ) then
        o = 1
        do while ( 2 * o - 1 < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  GLG
!  Default is Moderate Linear Growth.
!
    else if ( rule(dim) == 8 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 .or. growth(dim) == 0 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        o = 1
        do while ( 2 * o - 1 < 2 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 5 ) then
        o = 1
        do while ( 2 * o - 1 < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  GJ
!  Default is Moderate Linear Growth.
!
    else if ( rule(dim) == 9 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 .or. growth(dim) == 0 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        o = 1
        do while ( 2 * o - 1 < 2 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 5 ) then
        o = 1
        do while ( 2 * o - 1 < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  HGK
!  Default is Moderate Exponential Growth.
!  Exponential growth is interpreted to mean simply take successive rules.
!
    else if ( rule(dim) == 10 ) then

      if ( growth(dim) == 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Growth rate 1 for rule 10 not available.'
        stop
      else if ( growth(dim) == 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Growth rate 2 for rule 10 not available.'
        stop
      else if ( growth(dim) == 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Growth rate 3 for rule 10 not available.'
        stop
      else if ( growth(dim) == 4 ) then
        l = 0
        p = p_hgk(l)
        o = o_hgk(l)
        do while ( p < 2 * level(dim) + 1 )
          l = l + 1
          if ( 5 < l ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
            write ( *, '(a)' ) '  Hermite Genz-Keister maximum level exceeded.'
            stop
          end if
          p = p_hgk(l)
          o = o_hgk(l)
        end do
      else if ( growth(dim) == 5 .or. growth(dim) == 0 ) then
        l = 0
        p = p_hgk(l)
        o = o_hgk(l)
        do while ( p < 4 * level(dim) + 1 )
          l = l + 1
          if ( 5 < l ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
            write ( *, '(a)' ) '  Hermite Genz-Keister maximum level exceeded.'
            stop
          end if
          p = p_hgk(l)
          o = o_hgk(l)
        end do
      else if ( growth(dim) == 6 ) then
        l = level(dim);
        l = max ( l, 0 );
        if ( 5 < l ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LEVEL_GROWTH_TO_ORDER - Fatal error!'
          write ( *, '(a)' ) '  Hermite Genz-Keister maximum level exceeded.'
          stop
        end if
        o = o_hgk(l)
      end if
!
!  UO
!  Default is Moderate Linear Growth.
!  We'll assume the rule is of OPEN type and that it
!  has a precision typical of Gauss rules.
!
    else if ( rule(dim) == 11 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 .or. growth(dim) == 0 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        o = 1
        do while ( 2 * o - 1 < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 5 ) then
        o = 1
        do while ( 2 * o - 1 < 4 * level(dim) + 1 )
          o = 2 * o + 1
        end do
      else if ( growth(dim) == 6 ) then
        o =  2 ** ( level(dim) + 1 ) - 1
      end if
!
!  UC
!  Default is Moderate Linear Growth.
!  We'll assume the rule is of CLOSED type and that it
!  has a precision typical of Clenshaw-Curtis rules.
!
    else if ( rule(dim) == 12 ) then
      if ( growth(dim) == 1 ) then
        o = level(dim) + 1
      else if ( growth(dim) == 2 ) then
        o = 2 * ( ( level(dim) + 1 ) / 2 ) + 1
      else if ( growth(dim) == 3 .or. growth(dim) == 0 ) then
        o = 2 * level(dim) + 1
      else if ( growth(dim) == 4 ) then
        if ( level(dim) == 0 ) then
          o = 1
        else
          o = 2
          do while ( o < 2 * level(dim) + 1 )
            o = 2 * ( o - 1 ) + 1
          end do
        end if
      else if ( growth(dim) == 5 ) then
        if ( level(dim) == 0 ) then
          o = 1
        else
          o = 2
          do while ( o < 4 * level(dim) + 1 )
            o = 2 * ( o - 1 ) + 1
          end do
        end if
      else if ( growth(dim) == 6 ) then
        if ( level(dim) == 0 ) then
          o = 1
        else
          o =  2 ** level(dim) + 1
        end if
      end if

    end if

    order(dim) = o

  end do

  return
end
subroutine level_to_order_default ( dim_num, level, rule, order )

!*****************************************************************************80
!
!! LEVEL_TO_ORDER_DEFAULT: default growth.
!
!  Discussion:
!
!    This function uses:
!
!    * exponential growth rates for fully nested rules,
!      ( "CC", "F2", "GP");
!
!    * linear growth rates for most rules:
!      ( "GL", "GH", "GGH", "LG", "GLG", "GJ", "GW" ).
!
!    * slow exponential growth alternative for fully nested rules:
!      ("CC_SE", "F2_SE", "GP_SE").
!
!    * moderate exponential growth alternative for fully nested rules:
!      ("CC_ME", "F2_ME", "GP_ME").
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL(DIM_NUM), the 1D levels.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
!    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested.
!    12, "F2_SE", Fejer Type 2 Slow Exponential, Closed Fully Nested.
!    13, "GP_SE", Gauss Patterson Slow Exponential, Closed Fully Nested.
!    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested.
!    15, "F2_ME", Fejer Type 2 Moderate Exponential, Closed Fully Nested.
!    16, "GP_ME", Gauss Patterson Moderate Exponential, Closed Fully Nested.
!    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
!
!    Output, integer ( kind = 4 ) ORDER(DIM_NUM), the 1D orders
!    (number of points).
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) o
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) rule(dim_num)

  do dim = 1, dim_num

    if ( level(dim) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEVEL_TO_ORDER_DEFAULT - Fatal error!'
      write ( *, '(a)' ) '  LEVEL(DIM) < 0.'
      write ( *, '(a,i8,a,i8)' ) '  LEVEL(', dim, ') = ', level(dim)
      stop
    else if ( rule(dim) == 1 ) then
      if ( level(dim) == 0 ) then
        order(dim) = 1
      else
        order(dim) = ( 2**level(dim) ) + 1
      end if
    else if ( rule(dim) == 2 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 3 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 4 ) then
      order(dim) = 2 * level(dim) + 1
    else if ( rule(dim) == 5 ) then
      order(dim) = 2 * level(dim) + 1
    else if ( rule(dim) == 6 ) then
      order(dim) = 2 * level(dim) + 1
    else if ( rule(dim) == 7 ) then
      order(dim) = 2 * level(dim) + 1
    else if ( rule(dim) == 8 ) then
      order(dim) = 2 * level(dim) + 1
    else if ( rule(dim) == 9 ) then
      order(dim) = 2 * level(dim) + 1
    else if ( rule(dim) == 10 ) then
      order(dim) = 2 * level(dim) + 1
    else if ( rule(dim) == 11 ) then
      if ( level(dim) == 0 ) then
        o = 1
      else
        o = 2
        do while ( o < 2 * level(dim) + 1 )
          o = 2 * ( o - 1 ) + 1
        end do
      end if
      order(dim) = o
    else if ( rule(dim) == 12 ) then
      o = 1
      do while ( o < 2 * level(dim) + 1 )
        o = 2 * o + 1
      end do
      order(dim) = o
    else if ( rule(dim) == 13 ) then
      if ( level(dim) == 0 ) then
        order(dim) = 1
      else
        p = 5
        o = 3
        do while ( p < 2 * level(dim) + 1 )
          p = 2 * p + 1
          o = 2 * o + 1
        end do
        order(dim) = o
      end if
    else if ( rule(dim) == 14 ) then
      if ( level(dim) == 0 ) then
        o = 1
      else
        o = 2
        do while ( o < 4 * level(dim) + 1 )
          o = 2 * ( o - 1 ) + 1
        end do
      end if
      order(dim) = o
    else if ( rule(dim) == 15 ) then
      o = 1
      do while ( o < 4 * level(dim) + 1 )
        o = 2 * o + 1
      end do
      order(dim) = o
    else if ( rule(dim) == 16 ) then

      if ( level(dim) == 0 ) then
        order(dim) = 1
      else
        p = 5
        o = 3
        do while ( p < 4 * level(dim) + 1 )
          p = 2 * p + 1
          o = 2 * o + 1
        end do
        order(dim) = o
      end if
    else if ( rule(dim) == 17 ) then
      order(dim) = 2 * level(dim) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEVEL_TO_ORDER_DEFAULT - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) &
        '  Unexpected value of RULE(', dim, ') = ', rule(dim)
      stop
    end if

  end do

  return
end
subroutine level_to_order_exponential ( dim_num, level, rule, order )

!*****************************************************************************80
!
!! LEVEL_TO_ORDER_EXPONENTIAL: exponential growth.
!
!  Discussion:
!
!    In exponential growth, from one level to the next, the number of points
!    essentially doubles.
!
!    Closed rules:
!
!      O(0) = 1
!      O(L) = 2^L + 1;
!
!      O = 1, 3, 5, 9, 17, 33, ...
!
!    Open rules:
!
!      O(L) = 2^(L+1) - 1;
!
!      O = 1, 3, 7, 15, 31, 63, ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL(DIM_NUM), the 1D levels.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
!    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested.
!    12, "F2_SE", Fejer Type 2 Slow Exponential, Closed Fully Nested.
!    13, "GP_SE", Gauss Patterson Slow Exponential, Closed Fully Nested.
!    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested.
!    15, "F2_ME", Fejer Type 2 Moderate Exponential, Closed Fully Nested.
!    16, "GP_ME", Gauss Patterson Moderate Exponential, Closed Fully Nested.
!    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
!
!    Output, integer ( kind = 4 ) ORDER(DIM_NUM), the 1D orders
!    (number of points).
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) rule(dim_num)

  do dim = 1, dim_num

    if ( level(dim) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEVEL_TO_ORDER_EXPONENTIAL - Fatal error!'
      write ( *, '(a)' ) '  LEVEL(DIM) < 0.'
      write ( *, '(a,i8,a,i8)' ) '  LEVEL(', dim, ') = ', level(dim)
      stop
    else if ( rule(dim) == 1 ) then
      if ( level(dim) == 0 ) then
        order(dim) = 1
      else
        order(dim) = ( 2**level(dim) ) + 1
      end if
    else if ( rule(dim) == 2 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 3 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 4 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 5 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 6 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 7 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 8 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 9 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 10 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 11 ) then
      if ( level(dim) == 0 ) then
        order(dim) = 1
      else
        order(dim) = ( 2**level(dim) ) + 1
      end if
    else if ( rule(dim) == 12 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 13 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 14 ) then
      if ( level(dim) == 0 ) then
        order(dim) = 1
      else
        order(dim) = ( 2**level(dim) ) + 1
      end if
    else if ( rule(dim) == 15 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 16 ) then
      order(dim) = 2**( level(dim) + 1 ) - 1
    else if ( rule(dim) == 17 ) then
      if ( level(dim) == 0 ) then
        order(dim) = 1
      else
        order(dim) = ( 2**level(dim) ) + 1
      end if
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEVEL_TO_ORDER_EXPONENTIAL - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) &
        '  Unexpected value of RULE(', dim, ') = ', rule(dim)
      stop
    end if

  end do

  return
end
subroutine level_to_order_exponential_slow ( dim_num, level, rule, order )

!*****************************************************************************80
!
!! LEVEL_TO_ORDER_EXPONENTIAL_SLOW: slow exponential growth.
!
!  Discussion:
!
!    We seek a sequence of quadrature rules with two opposing constraints:
!    * a measured rise in polynomial precision with increasing level;
!    * a control on the increase in (new) points per level;
!
!    Essentially, we are trying to keep some of the advantages of nesting,
!    while moderating the cost of the explosive growth in order that occurs
!    due to the repeated order doubling of nesting.
!
!    We wish the number of points at a given level L to be "about" 2 * L + 1,
!    but we also wish the rules to be completely nested.
!
!    One way to do this is to start with a nested family of rules, whose
!    order will tend to grow exponentially (doubling from one to the next),
!    but simply to REPEAT each rule as many times as possible.  We move to
!    the next rule only when the desired precision 2 * L + 1 exceeds the
!    precision of the current rule.
!
!    For both the Clenshaw Curtis and Fejer Type 2 rules, the order and
!    precision are the same if the order is odd.   That is, an 11 point rule
!    will integrate exactly all polynomials up to and including degree 11.
!
!    For Gauss Patterson rules, the relationship between order and precision
!    is somewhat more complicated.  For that rule, we take the philosophy
!    that at each level L, we wish to choose the rule of smallest order
!    so that the precision of 2 * L + 1 is guaranteed.
!
!     L    2*L+1  CC Order    F2 Order    GP Order/Precision
!
!     0        1         1           1        1/1
!     1        3         3           3        3/5
!     2        5         5           7        3/5
!     3        7         9           7        7/11
!     4        9         9          15        7/11
!     5       11        17          15        7/11
!     6       13        17          15       15/23
!     7       15        17          15       15/23
!     8       17        17          31       15/23
!     9       19        33          31       15/23
!    10       21        33          31       15/23
!    11       23        33          31       15/23
!    12       25        33          31       31/47
!    13       27        33          31       31/47
!    14       29        33          31       31/47
!    15       31        33          31       31/47
!    16       33        33          63       31/47
!    17       35        65          63       31/47
!    18       37        65          63       31/47
!    19       39        65          63       31/47
!    20       41        65          63       31/47
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Knut Petras,
!    Smolyak Cubature of Given Polynomial Degree with Few Nodes
!    for Increasing Dimension,
!    Numerische Mathematik,
!    Volume 93, Number 4, February 2003, pages 729-753.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL(DIM_NUM), the 1D levels.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
!    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested.
!    12, "F2_SE", Fejer Type 2 Slow Exponential, Closed Fully Nested.
!    13, "GP_SE", Gauss Patterson Slow Exponential, Closed Fully Nested.
!    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested.
!    15, "F2_ME", Fejer Type 2 Moderate Exponential, Closed Fully Nested.
!    16, "GP_ME", Gauss Patterson Moderate Exponential, Closed Fully Nested.
!    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
!
!    Output, integer ( kind = 4 ) ORDER(DIM_NUM), the 1D orders
!    (number of points).
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) o
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) rule(dim_num)

  if ( any ( level(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVEL_TO_ORDER_EXPONENTIAL_SLOW - Fatal error!'
    write ( *, '(a)' ) '  Some entry of LEVEL is negative.'
    stop
  end if

  do dim = 1, dim_num

    if ( rule(dim) == 1 .or. rule(dim) == 11 .or. rule(dim) == 14 .or. &
         rule(dim) == 17 ) then

      if ( level(dim) == 0 ) then
        o = 1
      else
        o = 2
        do while ( o < 2 * level(dim) + 1 )
          o = 2 * ( o - 1 ) + 1
        end do
      end if

    else if ( rule(dim) == 3 .or. rule(dim) == 13 .or. rule(dim) == 16 ) then

      if ( level(dim) == 0 ) then
        o = 1
      else
        p = 5
        o = 3
        do while ( p < 2 * level(dim) + 1 )
          p = 2 * p + 1
          o = 2 * o + 1
        end do
      end if

    else

      o = 1
      do while ( o < 2 * level(dim) + 1 )
        o = 2 * o + 1
      end do

    end if

    order(dim) = o

  end do

  return
end
subroutine level_to_order_linear ( dim_num, level, rule, order )

!*****************************************************************************80
!
!! LEVEL_TO_ORDER_LINEAR: linear growth.
!
!  Discussion:
!
!      O(L) = 2 * L + 1;
!
!      O = 1, 3, 5, 7, 9, ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL(DIM_NUM), the 1D levels.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
!    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested.
!    12, "F2_SE", Fejer Type 2 Slow Exponential, Closed Fully Nested.
!    13, "GP_SE", Gauss Patterson Slow Exponential, Closed Fully Nested.
!    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested.
!    15, "F2_ME", Fejer Type 2 Moderate Exponential, Closed Fully Nested.
!    16, "GP_ME", Gauss Patterson Moderate Exponential, Closed Fully Nested.
!    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
!
!    Output, integer ( kind = 4 ) ORDER(DIM_NUM), the 1D orders
!    (number of points).
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)
  integer ( kind = 4 ) rule(dim_num)

  if ( any ( level(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVEL_TO_ORDER_LINEAR - Fatal error!'
    write ( *, '(a)' ) '  Some entry of LEVEL is negative.'
    stop
  end if

  order(1:dim_num) = 2 * level(1:dim_num) + 1

  return
end
subroutine nc_compute ( n, x_min, x_max, x, w )

!*****************************************************************************80
!
!! NC_COMPUTE computes a Newton-Cotes quadrature rule.
!
!  Discussion:
!
!    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule
!    estimates
!
!      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
!
!    using N abscissas X and weights W:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
!
!    For the CLOSED rule, the abscissas include the end points.
!    For an OPEN rule, the abscissas do not include the end points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X_MIN, X_MAX, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  real ( kind = 8 ) yvala
  real ( kind = 8 ) yvalb

  do i = 1, n
!
!  Compute the Lagrange basis polynomial which is 1 at X(I),
!  and zero at the other nodes.
!
    d(1:n) = 0.0D+00
    d(i) = 1.0D+00

    do j = 2, n
      do k = j, n
        d(n+j-k) = ( d(n+j-k-1) - d(n+j-k) ) / ( x(n+1-k) - x(n+j-k) )
      end do
    end do

    do j = 1, n - 1
      do k = 1, n - j
        d(n-k) = d(n-k) - x(n-k-j+1) * d(n-k+1)
      end do
    end do
!
!  Evaluate the antiderivative of the polynomial at the endpoints.
!
    yvala = d(n) / real ( n, kind = 8 )
    do j = n - 1, 1, -1
      yvala = yvala * x_min + d(j) / real ( j, kind = 8 )
    end do
    yvala = yvala * x_min

    yvalb = d(n) / real ( n, kind = 8 )
    do j = n - 1, 1, -1
      yvalb = yvalb * x_max + d(j) / real ( j, kind = 8 )
    end do
    yvalb = yvalb * x_max

    w(i) = yvalb - yvala

  end do

  return
end
subroutine ncc_compute_points ( n, x )

!*****************************************************************************80
!
!! NCC_COMPUTE_POINTS: Newton-Cotes Closed points
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  if ( n == 1 ) then

    x(1) = ( x_max + x_min ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = 8 ) * x_min   &
             + real (     i - 1, kind = 8 ) * x_max ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine ncc_compute_weights ( n, w )

!*****************************************************************************80
!
!! NCC_COMPUTE_WEIGHTS: Newton-Cotes Closed weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  if ( n == 1 ) then

    w(1) = x_max - x_min

  else

    call ncc_compute_points ( n, x )

    call nc_compute ( n, x_min, x_max, x, w )

  end if

  return
end
subroutine nco_compute_points ( n, x )

!*****************************************************************************80
!
!! NCO_COMPUTE_POINTS: points for a Newton-Cotes Open quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  do i = 1, n
    x(i) = ( real ( n - i + 1, kind = 8 ) * x_min   &
           + real (     i,     kind = 8 ) * x_max ) &
           / real ( n     + 1, kind = 8 )
  end do

  return
end
subroutine nco_compute_weights ( n, w )

!*****************************************************************************80
!
!! NCO_COMPUTE_WEIGHTS: weights for a Newton-Cotes Open quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  call nco_compute_points ( n, x )

  call nc_compute ( n, x_min, x_max, x, w )

  return
end
subroutine ncoh_compute_points ( n, x )

!*****************************************************************************80
!
!! NCOH_COMPUTE_POINTS: points for a Newton-Cotes Open Half quadrature rule.
!
!  Discussion:
!
!    The input value N is used to define N equal subintervals of [-1,+1].
!    The I-th abscissa is the center of the I-th subinterval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  do i = 1, n
    x(i) = ( real ( 2 * n - 2 * i + 1, kind = 8 ) * x_min   &
           + real (         2 * i - 1, kind = 8 ) * x_max ) &
           / real ( 2 * n,             kind = 8 )
  end do

  return
end
subroutine ncoh_compute_weights ( n, w )

!*****************************************************************************80
!
!! NCOH_COMPUTE_WEIGHTS: weights for a Newton-Cotes Open Half quadrature rule.
!
!  Discussion:
!
!    The input value N is used to define N equal subintervals of [-1,+1].
!    The I-th abscissa is the center of the I-th subinterval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  call ncoh_compute_points ( n, x )

  call nc_compute ( n, x_min, x_max, x, w )

  return
end
subroutine patterson_lookup ( order, points, weights )

!*****************************************************************************80
!
!! PATTERSON_LOOKUP returns the abscissas and weights of a Patterson rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1],
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    Legal values are 1, 3, 7, 15, 31, 63, 127, 255 and 511.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHTS(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) points(order)
  real ( kind = 8 ) weights(order)

  call patterson_lookup_points ( order, points )
  call patterson_lookup_weights ( order, weights );

  return
end
subroutine patterson_lookup_points ( n, points )

!*****************************************************************************80
!
!! PATTERSON_LOOKUP_POINTS returns the abscissas of a Patterson rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1],
!
!    These rules constitute a nested family.  The rules can integrate exactly
!    any polynomial of degree 1, 5, 11, 23, 47, 95, 191, 383 or 767, 
!    respectively.
!
!    The data for N = 511 was supplied by Dirk Laurie, and is derived
!    from a NAG Library function d01arf.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!    NAG Library Documentation,
!    D01ARF,
!    The Numerical Algorithms Group.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 7, 15, 31, 63, 127, 255 and 511.
!
!    Output, real ( kind = 8 ) POINTS(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) points(n)
  real ( kind = 8 ), save, dimension ( 1 ) :: x_001 = (/ &
     0.0D+00 /)
  real ( kind = 8 ), save, dimension ( 3 ) :: x_003 = (/ &
    -0.77459666924148337704D+00, &
     0.0D+00, &
     0.77459666924148337704D+00 /)
  real ( kind = 8 ), save, dimension ( 7 ) :: x_007 = (/ &
    -0.96049126870802028342D+00, &
    -0.77459666924148337704D+00, &
    -0.43424374934680255800D+00, &
     0.0D+00, &
     0.43424374934680255800D+00, &
     0.77459666924148337704D+00, &
     0.96049126870802028342D+00 /)
  real ( kind = 8 ), save, dimension ( 15 ) :: x_015 = (/ &
    -0.99383196321275502221D+00, &
    -0.96049126870802028342D+00, &
    -0.88845923287225699889D+00, &
    -0.77459666924148337704D+00, &
    -0.62110294673722640294D+00, &
    -0.43424374934680255800D+00, &
    -0.22338668642896688163D+00, &
     0.0D+00, &
     0.22338668642896688163D+00, &
     0.43424374934680255800D+00, &
     0.62110294673722640294D+00, &
     0.77459666924148337704D+00, &
     0.88845923287225699889D+00, &
     0.96049126870802028342D+00, &
     0.99383196321275502221D+00 /)
  real ( kind = 8 ), save, dimension ( 31 ) :: x_031 = (/ &
    -0.99909812496766759766D+00, &
    -0.99383196321275502221D+00, &
    -0.98153114955374010687D+00, &
    -0.96049126870802028342D+00, &
    -0.92965485742974005667D+00, &
    -0.88845923287225699889D+00, &
    -0.83672593816886873550D+00, &
    -0.77459666924148337704D+00, &
    -0.70249620649152707861D+00, &
    -0.62110294673722640294D+00, &
    -0.53131974364437562397D+00, &
    -0.43424374934680255800D+00, &
    -0.33113539325797683309D+00, &
    -0.22338668642896688163D+00, &
    -0.11248894313318662575D+00, &
     0.0D+00, &
     0.11248894313318662575D+00, &
     0.22338668642896688163D+00, &
     0.33113539325797683309D+00, &
     0.43424374934680255800D+00, &
     0.53131974364437562397D+00, &
     0.62110294673722640294D+00, &
     0.70249620649152707861D+00, &
     0.77459666924148337704D+00, &
     0.83672593816886873550D+00, &
     0.88845923287225699889D+00, &
     0.92965485742974005667D+00, &
     0.96049126870802028342D+00, &
     0.98153114955374010687D+00, &
     0.99383196321275502221D+00, &
     0.99909812496766759766D+00 /)
  real ( kind = 8 ), save, dimension ( 63 ) :: x_063 = (/ &
    -0.99987288812035761194D+00, &
    -0.99909812496766759766D+00, &
    -0.99720625937222195908D+00, &
    -0.99383196321275502221D+00, &
    -0.98868475754742947994D+00, &
    -0.98153114955374010687D+00, &
    -0.97218287474858179658D+00, &
    -0.96049126870802028342D+00, &
    -0.94634285837340290515D+00, &
    -0.92965485742974005667D+00, &
    -0.91037115695700429250D+00, &
    -0.88845923287225699889D+00, &
    -0.86390793819369047715D+00, &
    -0.83672593816886873550D+00, &
    -0.80694053195021761186D+00, &
    -0.77459666924148337704D+00, &
    -0.73975604435269475868D+00, &
    -0.70249620649152707861D+00, &
    -0.66290966002478059546D+00, &
    -0.62110294673722640294D+00, &
    -0.57719571005204581484D+00, &
    -0.53131974364437562397D+00, &
    -0.48361802694584102756D+00, &
    -0.43424374934680255800D+00, &
    -0.38335932419873034692D+00, &
    -0.33113539325797683309D+00, &
    -0.27774982202182431507D+00, &
    -0.22338668642896688163D+00, &
    -0.16823525155220746498D+00, &
    -0.11248894313318662575D+00, &
    -0.056344313046592789972D+00, &
     0.0D+00, &
     0.056344313046592789972D+00, &
     0.11248894313318662575D+00, &
     0.16823525155220746498D+00, &
     0.22338668642896688163D+00, &
     0.27774982202182431507D+00, &
     0.33113539325797683309D+00, &
     0.38335932419873034692D+00, &
     0.43424374934680255800D+00, &
     0.48361802694584102756D+00, &
     0.53131974364437562397D+00, &
     0.57719571005204581484D+00, &
     0.62110294673722640294D+00, &
     0.66290966002478059546D+00, &
     0.70249620649152707861D+00, &
     0.73975604435269475868D+00, &
     0.77459666924148337704D+00, &
     0.80694053195021761186D+00, &
     0.83672593816886873550D+00, &
     0.86390793819369047715D+00, &
     0.88845923287225699889D+00, &
     0.91037115695700429250D+00, &
     0.92965485742974005667D+00, &
     0.94634285837340290515D+00, &
     0.96049126870802028342D+00, &
     0.97218287474858179658D+00, &
     0.98153114955374010687D+00, &
     0.98868475754742947994D+00, &
     0.99383196321275502221D+00, &
     0.99720625937222195908D+00, &
     0.99909812496766759766D+00, &
     0.99987288812035761194D+00 /)
  real ( kind = 8 ), save, dimension ( 127 ) :: x_127 = (/ &
    -0.99998243035489159858D+00, &
    -0.99987288812035761194D+00, &
    -0.99959879967191068325D+00, &
    -0.99909812496766759766D+00, &
    -0.99831663531840739253D+00, &
    -0.99720625937222195908D+00, &
    -0.99572410469840718851D+00, &
    -0.99383196321275502221D+00, &
    -0.99149572117810613240D+00, &
    -0.98868475754742947994D+00, &
    -0.98537149959852037111D+00, &
    -0.98153114955374010687D+00, &
    -0.97714151463970571416D+00, &
    -0.97218287474858179658D+00, &
    -0.96663785155841656709D+00, &
    -0.96049126870802028342D+00, &
    -0.95373000642576113641D+00, &
    -0.94634285837340290515D+00, &
    -0.93832039777959288365D+00, &
    -0.92965485742974005667D+00, &
    -0.92034002547001242073D+00, &
    -0.91037115695700429250D+00, &
    -0.89974489977694003664D+00, &
    -0.88845923287225699889D+00, &
    -0.87651341448470526974D+00, &
    -0.86390793819369047715D+00, &
    -0.85064449476835027976D+00, &
    -0.83672593816886873550D+00, &
    -0.82215625436498040737D+00, &
    -0.80694053195021761186D+00, &
    -0.79108493379984836143D+00, &
    -0.77459666924148337704D+00, &
    -0.75748396638051363793D+00, &
    -0.73975604435269475868D+00, &
    -0.72142308537009891548D+00, &
    -0.70249620649152707861D+00, &
    -0.68298743109107922809D+00, &
    -0.66290966002478059546D+00, &
    -0.64227664250975951377D+00, &
    -0.62110294673722640294D+00, &
    -0.59940393024224289297D+00, &
    -0.57719571005204581484D+00, &
    -0.55449513263193254887D+00, &
    -0.53131974364437562397D+00, &
    -0.50768775753371660215D+00, &
    -0.48361802694584102756D+00, &
    -0.45913001198983233287D+00, &
    -0.43424374934680255800D+00, &
    -0.40897982122988867241D+00, &
    -0.38335932419873034692D+00, &
    -0.35740383783153215238D+00, &
    -0.33113539325797683309D+00, &
    -0.30457644155671404334D+00, &
    -0.27774982202182431507D+00, &
    -0.25067873030348317661D+00, &
    -0.22338668642896688163D+00, &
    -0.19589750271110015392D+00, &
    -0.16823525155220746498D+00, &
    -0.14042423315256017459D+00, &
    -0.11248894313318662575D+00, &
    -0.084454040083710883710D+00, &
    -0.056344313046592789972D+00, &
    -0.028184648949745694339D+00, &
     0.0D+00, &
     0.028184648949745694339D+00, &
     0.056344313046592789972D+00, &
     0.084454040083710883710D+00, &
     0.11248894313318662575D+00, &
     0.14042423315256017459D+00, &
     0.16823525155220746498D+00, &
     0.19589750271110015392D+00, &
     0.22338668642896688163D+00, &
     0.25067873030348317661D+00, &
     0.27774982202182431507D+00, &
     0.30457644155671404334D+00, &
     0.33113539325797683309D+00, &
     0.35740383783153215238D+00, &
     0.38335932419873034692D+00, &
     0.40897982122988867241D+00, &
     0.43424374934680255800D+00, &
     0.45913001198983233287D+00, &
     0.48361802694584102756D+00, &
     0.50768775753371660215D+00, &
     0.53131974364437562397D+00, &
     0.55449513263193254887D+00, &
     0.57719571005204581484D+00, &
     0.59940393024224289297D+00, &
     0.62110294673722640294D+00, &
     0.64227664250975951377D+00, &
     0.66290966002478059546D+00, &
     0.68298743109107922809D+00, &
     0.70249620649152707861D+00, &
     0.72142308537009891548D+00, &
     0.73975604435269475868D+00, &
     0.75748396638051363793D+00, &
     0.77459666924148337704D+00, &
     0.79108493379984836143D+00, &
     0.80694053195021761186D+00, &
     0.82215625436498040737D+00, &
     0.83672593816886873550D+00, &
     0.85064449476835027976D+00, &
     0.86390793819369047715D+00, &
     0.87651341448470526974D+00, &
     0.88845923287225699889D+00, &
     0.89974489977694003664D+00, &
     0.91037115695700429250D+00, &
     0.92034002547001242073D+00, &
     0.92965485742974005667D+00, &
     0.93832039777959288365D+00, &
     0.94634285837340290515D+00, &
     0.95373000642576113641D+00, &
     0.96049126870802028342D+00, &
     0.96663785155841656709D+00, &
     0.97218287474858179658D+00, &
     0.97714151463970571416D+00, &
     0.98153114955374010687D+00, &
     0.98537149959852037111D+00, &
     0.98868475754742947994D+00, &
     0.99149572117810613240D+00, &
     0.99383196321275502221D+00, &
     0.99572410469840718851D+00, &
     0.99720625937222195908D+00, &
     0.99831663531840739253D+00, &
     0.99909812496766759766D+00, &
     0.99959879967191068325D+00, &
     0.99987288812035761194D+00, &
     0.99998243035489159858D+00 /)

  real ( kind = 8 ), save, dimension ( 255 ) :: x_255 = (/ &
    -0.99999759637974846462D+00, &
    -0.99998243035489159858D+00, &
    -0.99994399620705437576D+00, &
    -0.99987288812035761194D+00, &
    -0.99976049092443204733D+00, &
    -0.99959879967191068325D+00, &
    -0.99938033802502358193D+00, &
    -0.99909812496766759766D+00, &
    -0.99874561446809511470D+00, &
    -0.99831663531840739253D+00, &
    -0.99780535449595727456D+00, &
    -0.99720625937222195908D+00, &
    -0.99651414591489027385D+00, &
    -0.99572410469840718851D+00, &
    -0.99483150280062100052D+00, &
    -0.99383196321275502221D+00, &
    -0.99272134428278861533D+00, &
    -0.99149572117810613240D+00, &
    -0.99015137040077015918D+00, &
    -0.98868475754742947994D+00, &
    -0.98709252795403406719D+00, &
    -0.98537149959852037111D+00, &
    -0.98351865757863272876D+00, &
    -0.98153114955374010687D+00, &
    -0.97940628167086268381D+00, &
    -0.97714151463970571416D+00, &
    -0.97473445975240266776D+00, &
    -0.97218287474858179658D+00, &
    -0.96948465950245923177D+00, &
    -0.96663785155841656709D+00, &
    -0.96364062156981213252D+00, &
    -0.96049126870802028342D+00, &
    -0.95718821610986096274D+00, &
    -0.95373000642576113641D+00, &
    -0.95011529752129487656D+00, &
    -0.94634285837340290515D+00, &
    -0.94241156519108305981D+00, &
    -0.93832039777959288365D+00, &
    -0.93406843615772578800D+00, &
    -0.92965485742974005667D+00, &
    -0.92507893290707565236D+00, &
    -0.92034002547001242073D+00, &
    -0.91543758715576504064D+00, &
    -0.91037115695700429250D+00, &
    -0.90514035881326159519D+00, &
    -0.89974489977694003664D+00, &
    -0.89418456833555902286D+00, &
    -0.88845923287225699889D+00, &
    -0.88256884024734190684D+00, &
    -0.87651341448470526974D+00, &
    -0.87029305554811390585D+00, &
    -0.86390793819369047715D+00, &
    -0.85735831088623215653D+00, &
    -0.85064449476835027976D+00, &
    -0.84376688267270860104D+00, &
    -0.83672593816886873550D+00, &
    -0.82952219463740140018D+00, &
    -0.82215625436498040737D+00, &
    -0.81462878765513741344D+00, &
    -0.80694053195021761186D+00, &
    -0.79909229096084140180D+00, &
    -0.79108493379984836143D+00, &
    -0.78291939411828301639D+00, &
    -0.77459666924148337704D+00, &
    -0.76611781930376009072D+00, &
    -0.75748396638051363793D+00, &
    -0.74869629361693660282D+00, &
    -0.73975604435269475868D+00, &
    -0.73066452124218126133D+00, &
    -0.72142308537009891548D+00, &
    -0.71203315536225203459D+00, &
    -0.70249620649152707861D+00, &
    -0.69281376977911470289D+00, &
    -0.68298743109107922809D+00, &
    -0.67301883023041847920D+00, &
    -0.66290966002478059546D+00, &
    -0.65266166541001749610D+00, &
    -0.64227664250975951377D+00, &
    -0.63175643771119423041D+00, &
    -0.62110294673722640294D+00, &
    -0.61031811371518640016D+00, &
    -0.59940393024224289297D+00, &
    -0.58836243444766254143D+00, &
    -0.57719571005204581484D+00, &
    -0.56590588542365442262D+00, &
    -0.55449513263193254887D+00, &
    -0.54296566649831149049D+00, &
    -0.53131974364437562397D+00, &
    -0.51955966153745702199D+00, &
    -0.50768775753371660215D+00, &
    -0.49570640791876146017D+00, &
    -0.48361802694584102756D+00, &
    -0.47142506587165887693D+00, &
    -0.45913001198983233287D+00, &
    -0.44673538766202847374D+00, &
    -0.43424374934680255800D+00, &
    -0.42165768662616330006D+00, &
    -0.40897982122988867241D+00, &
    -0.39621280605761593918D+00, &
    -0.38335932419873034692D+00, &
    -0.37042208795007823014D+00, &
    -0.35740383783153215238D+00, &
    -0.34430734159943802278D+00, &
    -0.33113539325797683309D+00, &
    -0.31789081206847668318D+00, &
    -0.30457644155671404334D+00, &
    -0.29119514851824668196D+00, &
    -0.27774982202182431507D+00, &
    -0.26424337241092676194D+00, &
    -0.25067873030348317661D+00, &
    -0.23705884558982972721D+00, &
    -0.22338668642896688163D+00, &
    -0.20966523824318119477D+00, &
    -0.19589750271110015392D+00, &
    -0.18208649675925219825D+00, &
    -0.16823525155220746498D+00, &
    -0.15434681148137810869D+00, &
    -0.14042423315256017459D+00, &
    -0.12647058437230196685D+00, &
    -0.11248894313318662575D+00, &
    -0.098482396598119202090D+00, &
    -0.084454040083710883710D+00, &
    -0.070406976042855179063D+00, &
    -0.056344313046592789972D+00, &
    -0.042269164765363603212D+00, &
    -0.028184648949745694339D+00, &
    -0.014093886410782462614D+00, &
    0.0D+00, &
    0.014093886410782462614D+00, &
    0.028184648949745694339D+00, &
    0.042269164765363603212D+00, &
    0.056344313046592789972D+00, &
    0.070406976042855179063D+00, &
    0.084454040083710883710D+00, &
    0.098482396598119202090D+00, &
    0.11248894313318662575D+00, &
    0.12647058437230196685D+00, &
    0.14042423315256017459D+00, &
    0.15434681148137810869D+00, &
    0.16823525155220746498D+00, &
    0.18208649675925219825D+00, &
    0.19589750271110015392D+00, &
    0.20966523824318119477D+00, &
    0.22338668642896688163D+00, &
    0.23705884558982972721D+00, &
    0.25067873030348317661D+00, &
    0.26424337241092676194D+00, &
    0.27774982202182431507D+00, &
    0.29119514851824668196D+00, &
    0.30457644155671404334D+00, &
    0.31789081206847668318D+00, &
    0.33113539325797683309D+00, &
    0.34430734159943802278D+00, &
    0.35740383783153215238D+00, &
    0.37042208795007823014D+00, &
    0.38335932419873034692D+00, &
    0.39621280605761593918D+00, &
    0.40897982122988867241D+00, &
    0.42165768662616330006D+00, &
    0.43424374934680255800D+00, &
    0.44673538766202847374D+00, &
    0.45913001198983233287D+00, &
    0.47142506587165887693D+00, &
    0.48361802694584102756D+00, &
    0.49570640791876146017D+00, &
    0.50768775753371660215D+00, &
    0.51955966153745702199D+00, &
    0.53131974364437562397D+00, &
    0.54296566649831149049D+00, &
    0.55449513263193254887D+00, &
    0.56590588542365442262D+00, &
    0.57719571005204581484D+00, &
    0.58836243444766254143D+00, &
    0.59940393024224289297D+00, &
    0.61031811371518640016D+00, &
    0.62110294673722640294D+00, &
    0.63175643771119423041D+00, &
    0.64227664250975951377D+00, &
    0.65266166541001749610D+00, &
    0.66290966002478059546D+00, &
    0.67301883023041847920D+00, &
    0.68298743109107922809D+00, &
    0.69281376977911470289D+00, &
    0.70249620649152707861D+00, &
    0.71203315536225203459D+00, &
    0.72142308537009891548D+00, &
    0.73066452124218126133D+00, &
    0.73975604435269475868D+00, &
    0.74869629361693660282D+00, &
    0.75748396638051363793D+00, &
    0.76611781930376009072D+00, &
    0.77459666924148337704D+00, &
    0.78291939411828301639D+00, &
    0.79108493379984836143D+00, &
    0.79909229096084140180D+00, &
    0.80694053195021761186D+00, &
    0.81462878765513741344D+00, &
    0.82215625436498040737D+00, &
    0.82952219463740140018D+00, &
    0.83672593816886873550D+00, &
    0.84376688267270860104D+00, &
    0.85064449476835027976D+00, &
    0.85735831088623215653D+00, &
    0.86390793819369047715D+00, &
    0.87029305554811390585D+00, &
    0.87651341448470526974D+00, &
    0.88256884024734190684D+00, &
    0.88845923287225699889D+00, &
    0.89418456833555902286D+00, &
    0.89974489977694003664D+00, &
    0.90514035881326159519D+00, &
    0.91037115695700429250D+00, &
    0.91543758715576504064D+00, &
    0.92034002547001242073D+00, &
    0.92507893290707565236D+00, &
    0.92965485742974005667D+00, &
    0.93406843615772578800D+00, &
    0.93832039777959288365D+00, &
    0.94241156519108305981D+00, &
    0.94634285837340290515D+00, &
    0.95011529752129487656D+00, &
    0.95373000642576113641D+00, &
    0.95718821610986096274D+00, &
    0.96049126870802028342D+00, &
    0.96364062156981213252D+00, &
    0.96663785155841656709D+00, &
    0.96948465950245923177D+00, &
    0.97218287474858179658D+00, &
    0.97473445975240266776D+00, &
    0.97714151463970571416D+00, &
    0.97940628167086268381D+00, &
    0.98153114955374010687D+00, &
    0.98351865757863272876D+00, &
    0.98537149959852037111D+00, &
    0.98709252795403406719D+00, &
    0.98868475754742947994D+00, &
    0.99015137040077015918D+00, &
    0.99149572117810613240D+00, &
    0.99272134428278861533D+00, &
    0.99383196321275502221D+00, &
    0.99483150280062100052D+00, &
    0.99572410469840718851D+00, &
    0.99651414591489027385D+00, &
    0.99720625937222195908D+00, &
    0.99780535449595727456D+00, &
    0.99831663531840739253D+00, &
    0.99874561446809511470D+00, &
    0.99909812496766759766D+00, &
    0.99938033802502358193D+00, &
    0.99959879967191068325D+00, &
    0.99976049092443204733D+00, &
    0.99987288812035761194D+00, &
    0.99994399620705437576D+00, &
    0.99998243035489159858D+00, &
    0.99999759637974846462D+00 /)

  real ( kind = 8 ), save, dimension ( 511 ) :: x_511 = (/ &
    -0.999999672956734384381E+00, &
    -0.999997596379748464620E+00, &
    -0.999992298136257588028E+00, &
    -0.999982430354891598580E+00, &
    -0.999966730098486276883E+00, &
    -0.999943996207054375764E+00, &
    -0.999913081144678282800E+00, &
    -0.999872888120357611938E+00, &
    -0.999822363679787739196E+00, &
    -0.999760490924432047330E+00, &
    -0.999686286448317731776E+00, &
    -0.999598799671910683252E+00, &
    -0.999497112467187190535E+00, &
    -0.999380338025023581928E+00, &
    -0.999247618943342473599E+00, &
    -0.999098124967667597662E+00, &
    -0.998931050830810562236E+00, &
    -0.998745614468095114704E+00, &
    -0.998541055697167906027E+00, &
    -0.998316635318407392531E+00, &
    -0.998071634524930323302E+00, &
    -0.997805354495957274562E+00, &
    -0.997517116063472399965E+00, &
    -0.997206259372221959076E+00, &
    -0.996872143485260161299E+00, &
    -0.996514145914890273849E+00, &
    -0.996131662079315037786E+00, &
    -0.995724104698407188509E+00, &
    -0.995290903148810302261E+00, &
    -0.994831502800621000519E+00, &
    -0.994345364356723405931E+00, &
    -0.993831963212755022209E+00, &
    -0.993290788851684966211E+00, &
    -0.992721344282788615328E+00, &
    -0.992123145530863117683E+00, &
    -0.991495721178106132399E+00, &
    -0.990838611958294243677E+00, &
    -0.990151370400770159181E+00, &
    -0.989433560520240838716E+00, &
    -0.988684757547429479939E+00, &
    -0.987904547695124280467E+00, &
    -0.987092527954034067190E+00, &
    -0.986248305913007552681E+00, &
    -0.985371499598520371114E+00, &
    -0.984461737328814534596E+00, &
    -0.983518657578632728762E+00, &
    -0.982541908851080604251E+00, &
    -0.981531149553740106867E+00, &
    -0.980486047876721339416E+00, &
    -0.979406281670862683806E+00, &
    -0.978291538324758539526E+00, &
    -0.977141514639705714156E+00, &
    -0.975955916702011753129E+00, &
    -0.974734459752402667761E+00, &
    -0.973476868052506926773E+00, &
    -0.972182874748581796578E+00, &
    -0.970852221732792443256E+00, &
    -0.969484659502459231771E+00, &
    -0.968079947017759947964E+00, &
    -0.966637851558416567092E+00, &
    -0.965158148579915665979E+00, &
    -0.963640621569812132521E+00, &
    -0.962085061904651475741E+00, &
    -0.960491268708020283423E+00, &
    -0.958859048710200221356E+00, &
    -0.957188216109860962736E+00, &
    -0.955478592438183697574E+00, &
    -0.953730006425761136415E+00, &
    -0.951942293872573589498E+00, &
    -0.950115297521294876558E+00, &
    -0.948248866934137357063E+00, &
    -0.946342858373402905148E+00, &
    -0.944397134685866648591E+00, &
    -0.942411565191083059813E+00, &
    -0.940386025573669721370E+00, &
    -0.938320397779592883655E+00, &
    -0.936214569916450806625E+00, &
    -0.934068436157725787999E+00, &
    -0.931881896650953639345E+00, &
    -0.929654857429740056670E+00, &
    -0.927387230329536696843E+00, &
    -0.925078932907075652364E+00, &
    -0.922729888363349241523E+00, &
    -0.920340025470012420730E+00, &
    -0.917909278499077501636E+00, &
    -0.915437587155765040644E+00, &
    -0.912924896514370590080E+00, &
    -0.910371156957004292498E+00, &
    -0.907776324115058903624E+00, &
    -0.905140358813261595189E+00, &
    -0.902463227016165675048E+00, &
    -0.899744899776940036639E+00, &
    -0.896985353188316590376E+00, &
    -0.894184568335559022859E+00, &
    -0.891342531251319871666E+00, &
    -0.888459232872256998890E+00, &
    -0.885534668997285008926E+00, &
    -0.882568840247341906842E+00, &
    -0.879561752026556262568E+00, &
    -0.876513414484705269742E+00, &
    -0.873423842480859310192E+00, &
    -0.870293055548113905851E+00, &
    -0.867121077859315215614E+00, &
    -0.863907938193690477146E+00, &
    -0.860653669904299969802E+00, &
    -0.857358310886232156525E+00, &
    -0.854021903545468625813E+00, &
    -0.850644494768350279758E+00, &
    -0.847226135891580884381E+00, &
    -0.843766882672708601038E+00, &
    -0.840266795261030442350E+00, &
    -0.836725938168868735503E+00, &
    -0.833144380243172624728E+00, &
    -0.829522194637401400178E+00, &
    -0.825859458783650001088E+00, &
    -0.822156254364980407373E+00, &
    -0.818412667287925807395E+00, &
    -0.814628787655137413436E+00, &
    -0.810804709738146594361E+00, &
    -0.806940531950217611856E+00, &
    -0.803036356819268687782E+00, &
    -0.799092290960841401800E+00, &
    -0.795108445051100526780E+00, &
    -0.791084933799848361435E+00, &
    -0.787021875923539422170E+00, &
    -0.782919394118283016385E+00, &
    -0.778777615032822744702E+00, &
    -0.774596669241483377036E+00, &
    -0.770376691217076824278E+00, &
    -0.766117819303760090717E+00, &
    -0.761820195689839149173E+00, &
    -0.757483966380513637926E+00, &
    -0.753109281170558142523E+00, &
    -0.748696293616936602823E+00, &
    -0.744245161011347082309E+00, &
    -0.739756044352694758677E+00, &
    -0.735229108319491547663E+00, &
    -0.730664521242181261329E+00, &
    -0.726062455075389632685E+00, &
    -0.721423085370098915485E+00, &
    -0.716746591245747095767E+00, &
    -0.712033155362252034587E+00, &
    -0.707282963891961103412E+00, &
    -0.702496206491527078610E+00, &
    -0.697673076273711232906E+00, &
    -0.692813769779114702895E+00, &
    -0.687918486947839325756E+00, &
    -0.682987431091079228087E+00, &
    -0.678020808862644517838E+00, &
    -0.673018830230418479199E+00, &
    -0.667981708447749702165E+00, &
    -0.662909660024780595461E+00, &
    -0.657802904699713735422E+00, &
    -0.652661665410017496101E+00, &
    -0.647486168263572388782E+00, &
    -0.642276642509759513774E+00, &
    -0.637033320510492495071E+00, &
    -0.631756437711194230414E+00, &
    -0.626446232611719746542E+00, &
    -0.621102946737226402941E+00, &
    -0.615726824608992638014E+00, &
    -0.610318113715186400156E+00, &
    -0.604877064481584353319E+00, &
    -0.599403930242242892974E+00, &
    -0.593898967210121954393E+00, &
    -0.588362434447662541434E+00, &
    -0.582794593837318850840E+00, &
    -0.577195710052045814844E+00, &
    -0.571566050525742833992E+00, &
    -0.565905885423654422623E+00, &
    -0.560215487612728441818E+00, &
    -0.554495132631932548866E+00, &
    -0.548745098662529448608E+00, &
    -0.542965666498311490492E+00, &
    -0.537157119515795115982E+00, &
    -0.531319743644375623972E+00, &
    -0.525453827336442687395E+00, &
    -0.519559661537457021993E+00, &
    -0.513637539655988578507E+00, &
    -0.507687757533716602155E+00, &
    -0.501710613415391878251E+00, &
    -0.495706407918761460170E+00, &
    -0.489675444004456155436E+00, &
    -0.483618026945841027562E+00, &
    -0.477534464298829155284E+00, &
    -0.471425065871658876934E+00, &
    -0.465290143694634735858E+00, &
    -0.459130011989832332874E+00, &
    -0.452944987140767283784E+00, &
    -0.446735387662028473742E+00, &
    -0.440501534168875795783E+00, &
    -0.434243749346802558002E+00, &
    -0.427962357921062742583E+00, &
    -0.421657686626163300056E+00, &
    -0.415330064175321663764E+00, &
    -0.408979821229888672409E+00, &
    -0.402607290368737092671E+00, &
    -0.396212806057615939183E+00, &
    -0.389796704618470795479E+00, &
    -0.383359324198730346916E+00, &
    -0.376901004740559344802E+00, &
    -0.370422087950078230138E+00, &
    -0.363922917266549655269E+00, &
    -0.357403837831532152376E+00, &
    -0.350865196458001209011E+00, &
    -0.344307341599438022777E+00, &
    -0.337730623318886219621E+00, &
    -0.331135393257976833093E+00, &
    -0.324522004605921855207E+00, &
    -0.317890812068476683182E+00, &
    -0.311242171836871800300E+00, &
    -0.304576441556714043335E+00, &
    -0.297893980296857823437E+00, &
    -0.291195148518246681964E+00, &
    -0.284480308042725577496E+00, &
    -0.277749822021824315065E+00, &
    -0.271004054905512543536E+00, &
    -0.264243372410926761945E+00, &
    -0.257468141491069790481E+00, &
    -0.250678730303483176613E+00, &
    -0.243875508178893021593E+00, &
    -0.237058845589829727213E+00, &
    -0.230229114119222177156E+00, &
    -0.223386686428966881628E+00, &
    -0.216531936228472628081E+00, &
    -0.209665238243181194766E+00, &
    -0.202786968183064697557E+00, &
    -0.195897502711100153915E+00, &
    -0.188997219411721861059E+00, &
    -0.182086496759252198246E+00, &
    -0.175165714086311475707E+00, &
    -0.168235251552207464982E+00, &
    -0.161295490111305257361E+00, &
    -0.154346811481378108692E+00, &
    -0.147389598111939940054E+00, &
    -0.140424233152560174594E+00, &
    -0.133451100421161601344E+00, &
    -0.126470584372301966851E+00, &
    -0.119483070065440005133E+00, &
    -0.112488943133186625746E+00, &
    -0.105488589749541988533E+00, &
    -0.984823965981192020903E-01, &
    -0.914707508403553909095E-01, &
    -0.844540400837108837102E-01, &
    -0.774326523498572825675E-01, &
    -0.704069760428551790633E-01, &
    -0.633773999173222898797E-01, &
    -0.563443130465927899720E-01, &
    -0.493081047908686267156E-01, &
    -0.422691647653636032124E-01, &
    -0.352278828084410232603E-01, &
    -0.281846489497456943394E-01, &
    -0.211398533783310883350E-01, &
    -0.140938864107824626142E-01, &
    -0.704713845933674648514E-02, &
    +0.000000000000000000000E+00, &
    +0.704713845933674648514E-02, &
    +0.140938864107824626142E-01, &
    +0.211398533783310883350E-01, &
    +0.281846489497456943394E-01, &
    +0.352278828084410232603E-01, &
    +0.422691647653636032124E-01, &
    +0.493081047908686267156E-01, &
    +0.563443130465927899720E-01, &
    +0.633773999173222898797E-01, &
    +0.704069760428551790633E-01, &
    +0.774326523498572825675E-01, &
    +0.844540400837108837102E-01, &
    +0.914707508403553909095E-01, &
    +0.984823965981192020903E-01, &
    +0.105488589749541988533E+00, &
    +0.112488943133186625746E+00, &
    +0.119483070065440005133E+00, &
    +0.126470584372301966851E+00, &
    +0.133451100421161601344E+00, &
    +0.140424233152560174594E+00, &
    +0.147389598111939940054E+00, &
    +0.154346811481378108692E+00, &
    +0.161295490111305257361E+00, &
    +0.168235251552207464982E+00, &
    +0.175165714086311475707E+00, &
    +0.182086496759252198246E+00, &
    +0.188997219411721861059E+00, &
    +0.195897502711100153915E+00, &
    +0.202786968183064697557E+00, &
    +0.209665238243181194766E+00, &
    +0.216531936228472628081E+00, &
    +0.223386686428966881628E+00, &
    +0.230229114119222177156E+00, &
    +0.237058845589829727213E+00, &
    +0.243875508178893021593E+00, &
    +0.250678730303483176613E+00, &
    +0.257468141491069790481E+00, &
    +0.264243372410926761945E+00, &
    +0.271004054905512543536E+00, &
    +0.277749822021824315065E+00, &
    +0.284480308042725577496E+00, &
    +0.291195148518246681964E+00, &
    +0.297893980296857823437E+00, &
    +0.304576441556714043335E+00, &
    +0.311242171836871800300E+00, &
    +0.317890812068476683182E+00, &
    +0.324522004605921855207E+00, &
    +0.331135393257976833093E+00, &
    +0.337730623318886219621E+00, &
    +0.344307341599438022777E+00, &
    +0.350865196458001209011E+00, &
    +0.357403837831532152376E+00, &
    +0.363922917266549655269E+00, &
    +0.370422087950078230138E+00, &
    +0.376901004740559344802E+00, &
    +0.383359324198730346916E+00, &
    +0.389796704618470795479E+00, &
    +0.396212806057615939183E+00, &
    +0.402607290368737092671E+00, &
    +0.408979821229888672409E+00, &
    +0.415330064175321663764E+00, &
    +0.421657686626163300056E+00, &
    +0.427962357921062742583E+00, &
    +0.434243749346802558002E+00, &
    +0.440501534168875795783E+00, &
    +0.446735387662028473742E+00, &
    +0.452944987140767283784E+00, &
    +0.459130011989832332874E+00, &
    +0.465290143694634735858E+00, &
    +0.471425065871658876934E+00, &
    +0.477534464298829155284E+00, &
    +0.483618026945841027562E+00, &
    +0.489675444004456155436E+00, &
    +0.495706407918761460170E+00, &
    +0.501710613415391878251E+00, &
    +0.507687757533716602155E+00, &
    +0.513637539655988578507E+00, &
    +0.519559661537457021993E+00, &
    +0.525453827336442687395E+00, &
    +0.531319743644375623972E+00, &
    +0.537157119515795115982E+00, &
    +0.542965666498311490492E+00, &
    +0.548745098662529448608E+00, &
    +0.554495132631932548866E+00, &
    +0.560215487612728441818E+00, &
    +0.565905885423654422623E+00, &
    +0.571566050525742833992E+00, &
    +0.577195710052045814844E+00, &
    +0.582794593837318850840E+00, &
    +0.588362434447662541434E+00, &
    +0.593898967210121954393E+00, &
    +0.599403930242242892974E+00, &
    +0.604877064481584353319E+00, &
    +0.610318113715186400156E+00, &
    +0.615726824608992638014E+00, &
    +0.621102946737226402941E+00, &
    +0.626446232611719746542E+00, &
    +0.631756437711194230414E+00, &
    +0.637033320510492495071E+00, &
    +0.642276642509759513774E+00, &
    +0.647486168263572388782E+00, &
    +0.652661665410017496101E+00, &
    +0.657802904699713735422E+00, &
    +0.662909660024780595461E+00, &
    +0.667981708447749702165E+00, &
    +0.673018830230418479199E+00, &
    +0.678020808862644517838E+00, &
    +0.682987431091079228087E+00, &
    +0.687918486947839325756E+00, &
    +0.692813769779114702895E+00, &
    +0.697673076273711232906E+00, &
    +0.702496206491527078610E+00, &
    +0.707282963891961103412E+00, &
    +0.712033155362252034587E+00, &
    +0.716746591245747095767E+00, &
    +0.721423085370098915485E+00, &
    +0.726062455075389632685E+00, &
    +0.730664521242181261329E+00, &
    +0.735229108319491547663E+00, &
    +0.739756044352694758677E+00, &
    +0.744245161011347082309E+00, &
    +0.748696293616936602823E+00, &
    +0.753109281170558142523E+00, &
    +0.757483966380513637926E+00, &
    +0.761820195689839149173E+00, &
    +0.766117819303760090717E+00, &
    +0.770376691217076824278E+00, &
    +0.774596669241483377036E+00, &
    +0.778777615032822744702E+00, &
    +0.782919394118283016385E+00, &
    +0.787021875923539422170E+00, &
    +0.791084933799848361435E+00, &
    +0.795108445051100526780E+00, &
    +0.799092290960841401800E+00, &
    +0.803036356819268687782E+00, &
    +0.806940531950217611856E+00, &
    +0.810804709738146594361E+00, &
    +0.814628787655137413436E+00, &
    +0.818412667287925807395E+00, &
    +0.822156254364980407373E+00, &
    +0.825859458783650001088E+00, &
    +0.829522194637401400178E+00, &
    +0.833144380243172624728E+00, &
    +0.836725938168868735503E+00, &
    +0.840266795261030442350E+00, &
    +0.843766882672708601038E+00, &
    +0.847226135891580884381E+00, &
    +0.850644494768350279758E+00, &
    +0.854021903545468625813E+00, &
    +0.857358310886232156525E+00, &
    +0.860653669904299969802E+00, &
    +0.863907938193690477146E+00, &
    +0.867121077859315215614E+00, &
    +0.870293055548113905851E+00, &
    +0.873423842480859310192E+00, &
    +0.876513414484705269742E+00, &
    +0.879561752026556262568E+00, &
    +0.882568840247341906842E+00, &
    +0.885534668997285008926E+00, &
    +0.888459232872256998890E+00, &
    +0.891342531251319871666E+00, &
    +0.894184568335559022859E+00, &
    +0.896985353188316590376E+00, &
    +0.899744899776940036639E+00, &
    +0.902463227016165675048E+00, &
    +0.905140358813261595189E+00, &
    +0.907776324115058903624E+00, &
    +0.910371156957004292498E+00, &
    +0.912924896514370590080E+00, &
    +0.915437587155765040644E+00, &
    +0.917909278499077501636E+00, &
    +0.920340025470012420730E+00, &
    +0.922729888363349241523E+00, &
    +0.925078932907075652364E+00, &
    +0.927387230329536696843E+00, &
    +0.929654857429740056670E+00, &
    +0.931881896650953639345E+00, &
    +0.934068436157725787999E+00, &
    +0.936214569916450806625E+00, &
    +0.938320397779592883655E+00, &
    +0.940386025573669721370E+00, &
    +0.942411565191083059813E+00, &
    +0.944397134685866648591E+00, &
    +0.946342858373402905148E+00, &
    +0.948248866934137357063E+00, &
    +0.950115297521294876558E+00, &
    +0.951942293872573589498E+00, &
    +0.953730006425761136415E+00, &
    +0.955478592438183697574E+00, &
    +0.957188216109860962736E+00, &
    +0.958859048710200221356E+00, &
    +0.960491268708020283423E+00, &
    +0.962085061904651475741E+00, &
    +0.963640621569812132521E+00, &
    +0.965158148579915665979E+00, &
    +0.966637851558416567092E+00, &
    +0.968079947017759947964E+00, &
    +0.969484659502459231771E+00, &
    +0.970852221732792443256E+00, &
    +0.972182874748581796578E+00, &
    +0.973476868052506926773E+00, &
    +0.974734459752402667761E+00, &
    +0.975955916702011753129E+00, &
    +0.977141514639705714156E+00, &
    +0.978291538324758539526E+00, &
    +0.979406281670862683806E+00, &
    +0.980486047876721339416E+00, &
    +0.981531149553740106867E+00, &
    +0.982541908851080604251E+00, &
    +0.983518657578632728762E+00, &
    +0.984461737328814534596E+00, &
    +0.985371499598520371114E+00, &
    +0.986248305913007552681E+00, &
    +0.987092527954034067190E+00, &
    +0.987904547695124280467E+00, &
    +0.988684757547429479939E+00, &
    +0.989433560520240838716E+00, &
    +0.990151370400770159181E+00, &
    +0.990838611958294243677E+00, &
    +0.991495721178106132399E+00, &
    +0.992123145530863117683E+00, &
    +0.992721344282788615328E+00, &
    +0.993290788851684966211E+00, &
    +0.993831963212755022209E+00, &
    +0.994345364356723405931E+00, &
    +0.994831502800621000519E+00, &
    +0.995290903148810302261E+00, &
    +0.995724104698407188509E+00, &
    +0.996131662079315037786E+00, &
    +0.996514145914890273849E+00, &
    +0.996872143485260161299E+00, &
    +0.997206259372221959076E+00, &
    +0.997517116063472399965E+00, &
    +0.997805354495957274562E+00, &
    +0.998071634524930323302E+00, &
    +0.998316635318407392531E+00, &
    +0.998541055697167906027E+00, &
    +0.998745614468095114704E+00, &
    +0.998931050830810562236E+00, &
    +0.999098124967667597662E+00, &
    +0.999247618943342473599E+00, &
    +0.999380338025023581928E+00, &
    +0.999497112467187190535E+00, &
    +0.999598799671910683252E+00, &
    +0.999686286448317731776E+00, &
    +0.999760490924432047330E+00, &
    +0.999822363679787739196E+00, &
    +0.999872888120357611938E+00, &
    +0.999913081144678282800E+00, &
    +0.999943996207054375764E+00, &
    +0.999966730098486276883E+00, &
    +0.999982430354891598580E+00, &
    +0.999992298136257588028E+00, &
    +0.999997596379748464620E+00, &
    +0.999999672956734384381E+00 /)

  if ( n == 1 ) then

    points(1:n) = x_001(1:n)

  else if ( n == 3 ) then

    points(1:n) = x_003(1:n)

  else if ( n == 7 ) then

    points(1:n) = x_007(1:n)

  else if ( n == 15 ) then

    points(1:n) = x_015(1:n)

  else if ( n == 31 ) then

    points(1:n) = x_031(1:n)

  else if ( n == 63 ) then

    points(1:n) = x_063(1:n)

  else if ( n == 127 ) then

    points(1:n) = x_127(1:n)

  else if ( n == 255 ) then

    points(1:n) = x_255(1:n)

  else if ( n == 511 ) then

    points(1:n) = x_511(1:n)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PATTERSON_LOOKUP_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Unexpected value of N = ', n
    stop

  end if

  return
end
subroutine patterson_lookup_points_np ( order, np, p, points )

!*****************************************************************************80
!
!! PATTERSON_LOOKUP_POINTS_NP returns the abscissas of a Patterson rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1],
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    Legal values are 1, 3, 7, 15, 31, 63, 127, 255 and 511.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) POINTS(ORDER), the abscissas.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) points(order)

  call patterson_lookup_points ( order, points )

  return
end
subroutine patterson_lookup_weights ( n, w )

!*****************************************************************************80
!
!! PATTERSON_LOOKUP_WEIGHTS sets weights for a Patterson rule.
!
!  Discussion:
!
!    The zeroth rule, of order 1, is the standard Legendre rule.
!
!    The first rule, of order 3, is the standard Legendre rule.
!
!    The second rule, of order 7, includes the abscissas of the previous
!    rule.
!
!    Each subsequent rule is nested in a similar way.  Rules are available
!    of orders 1, 3, 7, 15, 31, 63, 127, 255, and 511.
!
!    These rules constitute a nested family.  The rules can integrate exactly
!    any polynomial of degree 1, 5, 11, 23, 47, 95, 191, 383 or 767, 
!    respectively.
!
!    The data for N = 511 was supplied by Dirk Laurie, and is derived
!    from a NAG Library function d01arf.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!    NAG Library Documentation,
!    D01ARF,
!    The Numerical Algorithms Group.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 7, 15, 31, 63, 127, 255 or 511.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)

  if ( n == 1 ) then

    w(1) = 2.0D+00

  else if ( n == 3 ) then

    w(1) = 0.555555555555555555556D+00
    w(2) = 0.888888888888888888889D+00
    w(3) = 0.555555555555555555556D+00

  else if ( n == 7 ) then

    w(1) = 0.104656226026467265194D+00
    w(2) = 0.268488089868333440729D+00
    w(3) = 0.401397414775962222905D+00
    w(4) = 0.450916538658474142345D+00
    w(5) = 0.401397414775962222905D+00
    w(6) = 0.268488089868333440729D+00
    w(7) = 0.104656226026467265194D+00

  else if ( n == 15 ) then

    w( 1) = 0.0170017196299402603390D+00
    w( 2) = 0.0516032829970797396969D+00
    w( 3) = 0.0929271953151245376859D+00
    w( 4) = 0.134415255243784220360D+00
    w( 5) = 0.171511909136391380787D+00
    w( 6) = 0.200628529376989021034D+00
    w( 7) = 0.219156858401587496404D+00
    w( 8) = 0.225510499798206687386D+00
    w( 9) = 0.219156858401587496404D+00
    w(10) = 0.200628529376989021034D+00
    w(11) = 0.171511909136391380787D+00
    w(12) = 0.134415255243784220360D+00
    w(13) = 0.0929271953151245376859D+00
    w(14) = 0.0516032829970797396969D+00
    w(15) = 0.0170017196299402603390D+00

  else if ( n == 31 ) then

    w( 1) = 0.00254478079156187441540D+00
    w( 2) = 0.00843456573932110624631D+00
    w( 3) = 0.0164460498543878109338D+00
    w( 4) = 0.0258075980961766535646D+00
    w( 5) = 0.0359571033071293220968D+00
    w( 6) = 0.0464628932617579865414D+00
    w( 7) = 0.0569795094941233574122D+00
    w( 8) = 0.0672077542959907035404D+00
    w( 9) = 0.0768796204990035310427D+00
    w(10) = 0.0857559200499903511542D+00
    w(11) = 0.0936271099812644736167D+00
    w(12) = 0.100314278611795578771D+00
    w(13) = 0.105669893580234809744D+00
    w(14) = 0.109578421055924638237D+00
    w(15) = 0.111956873020953456880D+00
    w(16) = 0.112755256720768691607D+00
    w(17) = 0.111956873020953456880D+00
    w(18) = 0.109578421055924638237D+00
    w(19) = 0.105669893580234809744D+00
    w(20) = 0.100314278611795578771D+00
    w(21) = 0.0936271099812644736167D+00
    w(22) = 0.0857559200499903511542D+00
    w(23) = 0.0768796204990035310427D+00
    w(24) = 0.0672077542959907035404D+00
    w(25) = 0.0569795094941233574122D+00
    w(26) = 0.0464628932617579865414D+00
    w(27) = 0.0359571033071293220968D+00
    w(28) = 0.0258075980961766535646D+00
    w(29) = 0.0164460498543878109338D+00
    w(30) = 0.00843456573932110624631D+00
    w(31) = 0.00254478079156187441540D+00

  else if ( n == 63 ) then

    w( 1) = 0.000363221481845530659694D+00
    w( 2) = 0.00126515655623006801137D+00
    w( 3) = 0.00257904979468568827243D+00
    w( 4) = 0.00421763044155885483908D+00
    w( 5) = 0.00611550682211724633968D+00
    w( 6) = 0.00822300795723592966926D+00
    w( 7) = 0.0104982469096213218983D+00
    w( 8) = 0.0129038001003512656260D+00
    w( 9) = 0.0154067504665594978021D+00
    w(10) = 0.0179785515681282703329D+00
    w(11) = 0.0205942339159127111492D+00
    w(12) = 0.0232314466399102694433D+00
    w(13) = 0.0258696793272147469108D+00
    w(14) = 0.0284897547458335486125D+00
    w(15) = 0.0310735511116879648799D+00
    w(16) = 0.0336038771482077305417D+00
    w(17) = 0.0360644327807825726401D+00
    w(18) = 0.0384398102494555320386D+00
    w(19) = 0.0407155101169443189339D+00
    w(20) = 0.0428779600250077344929D+00
    w(21) = 0.0449145316536321974143D+00
    w(22) = 0.0468135549906280124026D+00
    w(23) = 0.0485643304066731987159D+00
    w(24) = 0.0501571393058995374137D+00
    w(25) = 0.0515832539520484587768D+00
    w(26) = 0.0528349467901165198621D+00
    w(27) = 0.0539054993352660639269D+00
    w(28) = 0.0547892105279628650322D+00
    w(29) = 0.0554814043565593639878D+00
    w(30) = 0.0559784365104763194076D+00
    w(31) = 0.0562776998312543012726D+00
    w(32) = 0.0563776283603847173877D+00
    w(33) = 0.0562776998312543012726D+00
    w(34) = 0.0559784365104763194076D+00
    w(35) = 0.0554814043565593639878D+00
    w(36) = 0.0547892105279628650322D+00
    w(37) = 0.0539054993352660639269D+00
    w(38) = 0.0528349467901165198621D+00
    w(39) = 0.0515832539520484587768D+00
    w(40) = 0.0501571393058995374137D+00
    w(41) = 0.0485643304066731987159D+00
    w(42) = 0.0468135549906280124026D+00
    w(43) = 0.0449145316536321974143D+00
    w(44) = 0.0428779600250077344929D+00
    w(45) = 0.0407155101169443189339D+00
    w(46) = 0.0384398102494555320386D+00
    w(47) = 0.0360644327807825726401D+00
    w(48) = 0.0336038771482077305417D+00
    w(49) = 0.0310735511116879648799D+00
    w(50) = 0.0284897547458335486125D+00
    w(51) = 0.0258696793272147469108D+00
    w(52) = 0.0232314466399102694433D+00
    w(53) = 0.0205942339159127111492D+00
    w(54) = 0.0179785515681282703329D+00
    w(55) = 0.0154067504665594978021D+00
    w(56) = 0.0129038001003512656260D+00
    w(57) = 0.0104982469096213218983D+00
    w(58) = 0.00822300795723592966926D+00
    w(59) = 0.00611550682211724633968D+00
    w(60) = 0.00421763044155885483908D+00
    w(61) = 0.00257904979468568827243D+00
    w(62) = 0.00126515655623006801137D+00
    w(63) = 0.000363221481845530659694D+00

  else if ( n == 127 ) then

    w(  1) = 0.0000505360952078625176247D+00
    w(  2) = 0.000180739564445388357820D+00
    w(  3) = 0.000377746646326984660274D+00
    w(  4) = 0.000632607319362633544219D+00
    w(  5) = 0.000938369848542381500794D+00
    w(  6) = 0.00128952408261041739210D+00
    w(  7) = 0.00168114286542146990631D+00
    w(  8) = 0.00210881524572663287933D+00
    w(  9) = 0.00256876494379402037313D+00
    w( 10) = 0.00305775341017553113613D+00
    w( 11) = 0.00357289278351729964938D+00
    w( 12) = 0.00411150397865469304717D+00
    w( 13) = 0.00467105037211432174741D+00
    w( 14) = 0.00524912345480885912513D+00
    w( 15) = 0.00584344987583563950756D+00
    w( 16) = 0.00645190005017573692280D+00
    w( 17) = 0.00707248999543355546805D+00
    w( 18) = 0.00770337523327974184817D+00
    w( 19) = 0.00834283875396815770558D+00
    w( 20) = 0.00898927578406413572328D+00
    w( 21) = 0.00964117772970253669530D+00
    w( 22) = 0.0102971169579563555237D+00
    w( 23) = 0.0109557333878379016480D+00
    w( 24) = 0.0116157233199551347270D+00
    w( 25) = 0.0122758305600827700870D+00
    w( 26) = 0.0129348396636073734547D+00
    w( 27) = 0.0135915710097655467896D+00
    w( 28) = 0.0142448773729167743063D+00
    w( 29) = 0.0148936416648151820348D+00
    w( 30) = 0.0155367755558439824399D+00
    w( 31) = 0.0161732187295777199419D+00
    w( 32) = 0.0168019385741038652709D+00
    w( 33) = 0.0174219301594641737472D+00
    w( 34) = 0.0180322163903912863201D+00
    w( 35) = 0.0186318482561387901863D+00
    w( 36) = 0.0192199051247277660193D+00
    w( 37) = 0.0197954950480974994880D+00
    w( 38) = 0.0203577550584721594669D+00
    w( 39) = 0.0209058514458120238522D+00
    w( 40) = 0.0214389800125038672465D+00
    w( 41) = 0.0219563663053178249393D+00
    w( 42) = 0.0224572658268160987071D+00
    w( 43) = 0.0229409642293877487608D+00
    w( 44) = 0.0234067774953140062013D+00
    w( 45) = 0.0238540521060385400804D+00
    w( 46) = 0.0242821652033365993580D+00
    w( 47) = 0.0246905247444876769091D+00
    w( 48) = 0.0250785696529497687068D+00
    w( 49) = 0.0254457699654647658126D+00
    w( 50) = 0.0257916269760242293884D+00
    w( 51) = 0.0261156733767060976805D+00
    w( 52) = 0.0264174733950582599310D+00
    w( 53) = 0.0266966229274503599062D+00
    w( 54) = 0.0269527496676330319634D+00
    w( 55) = 0.0271855132296247918192D+00
    w( 56) = 0.0273946052639814325161D+00
    w( 57) = 0.0275797495664818730349D+00
    w( 58) = 0.0277407021782796819939D+00
    w( 59) = 0.0278772514766137016085D+00
    w( 60) = 0.0279892182552381597038D+00
    w( 61) = 0.0280764557938172466068D+00
    w( 62) = 0.0281388499156271506363D+00
    w( 63) = 0.0281763190330166021307D+00
    w( 64) = 0.0281888141801923586938D+00
    w( 65) = 0.0281763190330166021307D+00
    w( 66) = 0.0281388499156271506363D+00
    w( 67) = 0.0280764557938172466068D+00
    w( 68) = 0.0279892182552381597038D+00
    w( 69) = 0.0278772514766137016085D+00
    w( 70) = 0.0277407021782796819939D+00
    w( 71) = 0.0275797495664818730349D+00
    w( 72) = 0.0273946052639814325161D+00
    w( 73) = 0.0271855132296247918192D+00
    w( 74) = 0.0269527496676330319634D+00
    w( 75) = 0.0266966229274503599062D+00
    w( 76) = 0.0264174733950582599310D+00
    w( 77) = 0.0261156733767060976805D+00
    w( 78) = 0.0257916269760242293884D+00
    w( 79) = 0.0254457699654647658126D+00
    w( 80) = 0.0250785696529497687068D+00
    w( 81) = 0.0246905247444876769091D+00
    w( 82) = 0.0242821652033365993580D+00
    w( 83) = 0.0238540521060385400804D+00
    w( 84) = 0.0234067774953140062013D+00
    w( 85) = 0.0229409642293877487608D+00
    w( 86) = 0.0224572658268160987071D+00
    w( 87) = 0.0219563663053178249393D+00
    w( 88) = 0.0214389800125038672465D+00
    w( 89) = 0.0209058514458120238522D+00
    w( 90) = 0.0203577550584721594669D+00
    w( 91) = 0.0197954950480974994880D+00
    w( 92) = 0.0192199051247277660193D+00
    w( 93) = 0.0186318482561387901863D+00
    w( 94) = 0.0180322163903912863201D+00
    w( 95) = 0.0174219301594641737472D+00
    w( 96) = 0.0168019385741038652709D+00
    w( 97) = 0.0161732187295777199419D+00
    w( 98) = 0.0155367755558439824399D+00
    w( 99) = 0.0148936416648151820348D+00
    w(100) = 0.0142448773729167743063D+00
    w(101) = 0.0135915710097655467896D+00
    w(102) = 0.0129348396636073734547D+00
    w(103) = 0.0122758305600827700870D+00
    w(104) = 0.0116157233199551347270D+00
    w(105) = 0.0109557333878379016480D+00
    w(106) = 0.0102971169579563555237D+00
    w(107) = 0.00964117772970253669530D+00
    w(108) = 0.00898927578406413572328D+00
    w(109) = 0.00834283875396815770558D+00
    w(110) = 0.00770337523327974184817D+00
    w(111) = 0.00707248999543355546805D+00
    w(112) = 0.00645190005017573692280D+00
    w(113) = 0.00584344987583563950756D+00
    w(114) = 0.00524912345480885912513D+00
    w(115) = 0.00467105037211432174741D+00
    w(116) = 0.00411150397865469304717D+00
    w(117) = 0.00357289278351729964938D+00
    w(118) = 0.00305775341017553113613D+00
    w(119) = 0.00256876494379402037313D+00
    w(120) = 0.00210881524572663287933D+00
    w(121) = 0.00168114286542146990631D+00
    w(122) = 0.00128952408261041739210D+00
    w(123) = 0.000938369848542381500794D+00
    w(124) = 0.000632607319362633544219D+00
    w(125) = 0.000377746646326984660274D+00
    w(126) = 0.000180739564445388357820D+00
    w(127) = 0.0000505360952078625176247D+00

  else if ( n == 255 ) then

    w(1) = 0.69379364324108267170D-05
    w(2) = 0.25157870384280661489D-04
    w(3) = 0.53275293669780613125D-04
    w(4) = 0.90372734658751149261D-04
    w(5) = 0.13575491094922871973D-03
    w(6) = 0.18887326450650491366D-03
    w(7) = 0.24921240048299729402D-03
    w(8) = 0.31630366082226447689D-03
    w(9) = 0.38974528447328229322D-03
    w(10) = 0.46918492424785040975D-03
    w(11) = 0.55429531493037471492D-03
    w(12) = 0.64476204130572477933D-03
    w(13) = 0.74028280424450333046D-03
    w(14) = 0.84057143271072246365D-03
    w(15) = 0.94536151685852538246D-03
    w(16) = 0.10544076228633167722D-02
    w(17) = 0.11674841174299594077D-02
    w(18) = 0.12843824718970101768D-02
    w(19) = 0.14049079956551446427D-02
    w(20) = 0.15288767050877655684D-02
    w(21) = 0.16561127281544526052D-02
    w(22) = 0.17864463917586498247D-02
    w(23) = 0.19197129710138724125D-02
    w(24) = 0.20557519893273465236D-02
    w(25) = 0.21944069253638388388D-02
    w(26) = 0.23355251860571608737D-02
    w(27) = 0.24789582266575679307D-02
    w(28) = 0.26245617274044295626D-02
    w(29) = 0.27721957645934509940D-02
    w(30) = 0.29217249379178197538D-02
    w(31) = 0.30730184347025783234D-02
    w(32) = 0.32259500250878684614D-02
    w(33) = 0.33803979910869203823D-02
    w(34) = 0.35362449977167777340D-02
    w(35) = 0.36933779170256508183D-02
    w(36) = 0.38516876166398709241D-02
    w(37) = 0.40110687240750233989D-02
    w(38) = 0.41714193769840788528D-02
    w(39) = 0.43326409680929828545D-02
    w(40) = 0.44946378920320678616D-02
    w(41) = 0.46573172997568547773D-02
    w(42) = 0.48205888648512683476D-02
    w(43) = 0.49843645647655386012D-02
    w(44) = 0.51485584789781777618D-02
    w(45) = 0.53130866051870565663D-02
    w(46) = 0.54778666939189508240D-02
    w(47) = 0.56428181013844441585D-02
    w(48) = 0.58078616599775673635D-02
    w(49) = 0.59729195655081658049D-02
    w(50) = 0.61379152800413850435D-02
    w(51) = 0.63027734490857587172D-02
    w(52) = 0.64674198318036867274D-02
    w(53) = 0.66317812429018878941D-02
    w(54) = 0.67957855048827733948D-02
    w(55) = 0.69593614093904229394D-02
    w(56) = 0.71224386864583871532D-02
    w(57) = 0.72849479805538070639D-02
    w(58) = 0.74468208324075910174D-02
    w(59) = 0.76079896657190565832D-02
    w(60) = 0.77683877779219912200D-02
    w(61) = 0.79279493342948491103D-02
    w(62) = 0.80866093647888599710D-02
    w(63) = 0.82443037630328680306D-02
    w(64) = 0.84009692870519326354D-02
    w(65) = 0.85565435613076896192D-02
    w(66) = 0.87109650797320868736D-02
    w(67) = 0.88641732094824942641D-02
    w(68) = 0.90161081951956431600D-02
    w(69) = 0.91667111635607884067D-02
    w(70) = 0.93159241280693950932D-02
    w(71) = 0.94636899938300652943D-02
    w(72) = 0.96099525623638830097D-02
    w(73) = 0.97546565363174114611D-02
    w(74) = 0.98977475240487497440D-02
    w(75) = 0.10039172044056840798D-01
    w(76) = 0.10178877529236079733D-01
    w(77) = 0.10316812330947621682D-01
    w(78) = 0.10452925722906011926D-01
    w(79) = 0.10587167904885197931D-01
    w(80) = 0.10719490006251933623D-01
    w(81) = 0.10849844089337314099D-01
    w(82) = 0.10978183152658912470D-01
    w(83) = 0.11104461134006926537D-01
    w(84) = 0.11228632913408049354D-01
    w(85) = 0.11350654315980596602D-01
    w(86) = 0.11470482114693874380D-01
    w(87) = 0.11588074033043952568D-01
    w(88) = 0.11703388747657003101D-01
    w(89) = 0.11816385890830235763D-01
    w(90) = 0.11927026053019270040D-01
    w(91) = 0.12035270785279562630D-01
    w(92) = 0.12141082601668299679D-01
    w(93) = 0.12244424981611985899D-01
    w(94) = 0.12345262372243838455D-01
    w(95) = 0.12443560190714035263D-01
    w(96) = 0.12539284826474884353D-01
    w(97) = 0.12632403643542078765D-01
    w(98) = 0.12722884982732382906D-01
    w(99) = 0.12810698163877361967D-01
    w(100) = 0.12895813488012114694D-01
    w(101) = 0.12978202239537399286D-01
    w(102) = 0.13057836688353048840D-01
    w(103) = 0.13134690091960152836D-01
    w(104) = 0.13208736697529129966D-01
    w(105) = 0.13279951743930530650D-01
    w(106) = 0.13348311463725179953D-01
    w(107) = 0.13413793085110098513D-01
    w(108) = 0.13476374833816515982D-01
    w(109) = 0.13536035934956213614D-01
    w(110) = 0.13592756614812395910D-01
    w(111) = 0.13646518102571291428D-01
    w(112) = 0.13697302631990716258D-01
    w(113) = 0.13745093443001896632D-01
    w(114) = 0.13789874783240936517D-01
    w(115) = 0.13831631909506428676D-01
    w(116) = 0.13870351089139840997D-01
    w(117) = 0.13906019601325461264D-01
    w(118) = 0.13938625738306850804D-01
    w(119) = 0.13968158806516938516D-01
    w(120) = 0.13994609127619079852D-01
    w(121) = 0.14017968039456608810D-01
    w(122) = 0.14038227896908623303D-01
    w(123) = 0.14055382072649964277D-01
    w(124) = 0.14069424957813575318D-01
    w(125) = 0.14080351962553661325D-01
    w(126) = 0.14088159516508301065D-01
    w(127) = 0.14092845069160408355D-01
    w(128) = 0.14094407090096179347D-01
    w(129) = 0.14092845069160408355D-01
    w(130) = 0.14088159516508301065D-01
    w(131) = 0.14080351962553661325D-01
    w(132) = 0.14069424957813575318D-01
    w(133) = 0.14055382072649964277D-01
    w(134) = 0.14038227896908623303D-01
    w(135) = 0.14017968039456608810D-01
    w(136) = 0.13994609127619079852D-01
    w(137) = 0.13968158806516938516D-01
    w(138) = 0.13938625738306850804D-01
    w(139) = 0.13906019601325461264D-01
    w(140) = 0.13870351089139840997D-01
    w(141) = 0.13831631909506428676D-01
    w(142) = 0.13789874783240936517D-01
    w(143) = 0.13745093443001896632D-01
    w(144) = 0.13697302631990716258D-01
    w(145) = 0.13646518102571291428D-01
    w(146) = 0.13592756614812395910D-01
    w(147) = 0.13536035934956213614D-01
    w(148) = 0.13476374833816515982D-01
    w(149) = 0.13413793085110098513D-01
    w(150) = 0.13348311463725179953D-01
    w(151) = 0.13279951743930530650D-01
    w(152) = 0.13208736697529129966D-01
    w(153) = 0.13134690091960152836D-01
    w(154) = 0.13057836688353048840D-01
    w(155) = 0.12978202239537399286D-01
    w(156) = 0.12895813488012114694D-01
    w(157) = 0.12810698163877361967D-01
    w(158) = 0.12722884982732382906D-01
    w(159) = 0.12632403643542078765D-01
    w(160) = 0.12539284826474884353D-01
    w(161) = 0.12443560190714035263D-01
    w(162) = 0.12345262372243838455D-01
    w(163) = 0.12244424981611985899D-01
    w(164) = 0.12141082601668299679D-01
    w(165) = 0.12035270785279562630D-01
    w(166) = 0.11927026053019270040D-01
    w(167) = 0.11816385890830235763D-01
    w(168) = 0.11703388747657003101D-01
    w(169) = 0.11588074033043952568D-01
    w(170) = 0.11470482114693874380D-01
    w(171) = 0.11350654315980596602D-01
    w(172) = 0.11228632913408049354D-01
    w(173) = 0.11104461134006926537D-01
    w(174) = 0.10978183152658912470D-01
    w(175) = 0.10849844089337314099D-01
    w(176) = 0.10719490006251933623D-01
    w(177) = 0.10587167904885197931D-01
    w(178) = 0.10452925722906011926D-01
    w(179) = 0.10316812330947621682D-01
    w(180) = 0.10178877529236079733D-01
    w(181) = 0.10039172044056840798D-01
    w(182) = 0.98977475240487497440D-02
    w(183) = 0.97546565363174114611D-02
    w(184) = 0.96099525623638830097D-02
    w(185) = 0.94636899938300652943D-02
    w(186) = 0.93159241280693950932D-02
    w(187) = 0.91667111635607884067D-02
    w(188) = 0.90161081951956431600D-02
    w(189) = 0.88641732094824942641D-02
    w(190) = 0.87109650797320868736D-02
    w(191) = 0.85565435613076896192D-02
    w(192) = 0.84009692870519326354D-02
    w(193) = 0.82443037630328680306D-02
    w(194) = 0.80866093647888599710D-02
    w(195) = 0.79279493342948491103D-02
    w(196) = 0.77683877779219912200D-02
    w(197) = 0.76079896657190565832D-02
    w(198) = 0.74468208324075910174D-02
    w(199) = 0.72849479805538070639D-02
    w(200) = 0.71224386864583871532D-02
    w(201) = 0.69593614093904229394D-02
    w(202) = 0.67957855048827733948D-02
    w(203) = 0.66317812429018878941D-02
    w(204) = 0.64674198318036867274D-02
    w(205) = 0.63027734490857587172D-02
    w(206) = 0.61379152800413850435D-02
    w(207) = 0.59729195655081658049D-02
    w(208) = 0.58078616599775673635D-02
    w(209) = 0.56428181013844441585D-02
    w(210) = 0.54778666939189508240D-02
    w(211) = 0.53130866051870565663D-02
    w(212) = 0.51485584789781777618D-02
    w(213) = 0.49843645647655386012D-02
    w(214) = 0.48205888648512683476D-02
    w(215) = 0.46573172997568547773D-02
    w(216) = 0.44946378920320678616D-02
    w(217) = 0.43326409680929828545D-02
    w(218) = 0.41714193769840788528D-02
    w(219) = 0.40110687240750233989D-02
    w(220) = 0.38516876166398709241D-02
    w(221) = 0.36933779170256508183D-02
    w(222) = 0.35362449977167777340D-02
    w(223) = 0.33803979910869203823D-02
    w(224) = 0.32259500250878684614D-02
    w(225) = 0.30730184347025783234D-02
    w(226) = 0.29217249379178197538D-02
    w(227) = 0.27721957645934509940D-02
    w(228) = 0.26245617274044295626D-02
    w(229) = 0.24789582266575679307D-02
    w(230) = 0.23355251860571608737D-02
    w(231) = 0.21944069253638388388D-02
    w(232) = 0.20557519893273465236D-02
    w(233) = 0.19197129710138724125D-02
    w(234) = 0.17864463917586498247D-02
    w(235) = 0.16561127281544526052D-02
    w(236) = 0.15288767050877655684D-02
    w(237) = 0.14049079956551446427D-02
    w(238) = 0.12843824718970101768D-02
    w(239) = 0.11674841174299594077D-02
    w(240) = 0.10544076228633167722D-02
    w(241) = 0.94536151685852538246D-03
    w(242) = 0.84057143271072246365D-03
    w(243) = 0.74028280424450333046D-03
    w(244) = 0.64476204130572477933D-03
    w(245) = 0.55429531493037471492D-03
    w(246) = 0.46918492424785040975D-03
    w(247) = 0.38974528447328229322D-03
    w(248) = 0.31630366082226447689D-03
    w(249) = 0.24921240048299729402D-03
    w(250) = 0.18887326450650491366D-03
    w(251) = 0.13575491094922871973D-03
    w(252) = 0.90372734658751149261D-04
    w(253) = 0.53275293669780613125D-04
    w(254) = 0.25157870384280661489D-04
    w(255) = 0.69379364324108267170D-05

  else if ( n == 511 ) then

    w(  1) = 0.945715933950007048827D-06
    w(  2) = 0.345456507169149134898D-05
    w(  3) = 0.736624069102321668857D-05
    w(  4) = 0.125792781889592743525D-04
    w(  5) = 0.190213681905875816679D-04
    w(  6) = 0.266376412339000901358D-04
    w(  7) = 0.353751372055189588628D-04
    w(  8) = 0.451863674126296143105D-04
    w(  9) = 0.560319507856164252140D-04
    w( 10) = 0.678774554733972416227D-04
    w( 11) = 0.806899228014035293851D-04
    w( 12) = 0.944366322532705527066D-04
    w( 13) = 0.109085545645741522051D-03
    w( 14) = 0.124606200241498368482D-03
    w( 15) = 0.140970302204104791413D-03
    w( 16) = 0.158151830411132242924D-03
    w( 17) = 0.176126765545083195474D-03
    w( 18) = 0.194872642236641146532D-03
    w( 19) = 0.214368090034216937149D-03
    w( 20) = 0.234592462123925204879D-03
    w( 21) = 0.255525589595236862014D-03
    w( 22) = 0.277147657465187357459D-03
    w( 23) = 0.299439176850911730874D-03
    w( 24) = 0.322381020652862389664D-03
    w( 25) = 0.345954492129903871350D-03
    w( 26) = 0.370141402122251665232D-03
    w( 27) = 0.394924138246873704434D-03
    w( 28) = 0.420285716355361231823D-03
    w( 29) = 0.446209810101403247488D-03
    w( 30) = 0.472680758429262691232D-03
    w( 31) = 0.499683553312800484519D-03
    w( 32) = 0.527203811431658386125D-03
    w( 33) = 0.555227733977307579715D-03
    w( 34) = 0.583742058714979703847D-03
    w( 35) = 0.612734008012225209294D-03
    w( 36) = 0.642191235948505088403D-03
    w( 37) = 0.672101776960108194646D-03
    w( 38) = 0.702453997827572321358D-03
    w( 39) = 0.733236554224767912055D-03
    w( 40) = 0.764438352543882784191D-03
    w( 41) = 0.796048517297550871506D-03
    w( 42) = 0.828056364077226302608D-03
    w( 43) = 0.860451377808527848128D-03
    w( 44) = 0.893223195879324912340D-03
    w( 45) = 0.926361595613111283368D-03
    w( 46) = 0.959856485506936206261D-03
    w( 47) = 0.993697899638760857945D-03
    w( 48) = 0.102787599466367326179D-02
    w( 49) = 0.106238104885340071375D-02
    w( 50) = 0.109720346268191941940D-02
    w( 51) = 0.113233376051597664917D-02
    w( 52) = 0.116776259302858043685D-02
    w( 53) = 0.120348074001265964881D-02
    w( 54) = 0.123947911332878396534D-02
    w( 55) = 0.127574875977346947345D-02
    w( 56) = 0.131228086370221478128D-02
    w( 57) = 0.134906674928353113127D-02
    w( 58) = 0.138609788229672549700D-02
    w( 59) = 0.142336587141720519900D-02
    w( 60) = 0.146086246895890987689D-02
    w( 61) = 0.149857957106456636214D-02
    w( 62) = 0.153650921735128916170D-02
    w( 63) = 0.157464359003212166189D-02
    w( 64) = 0.161297501254393423070D-02
    w( 65) = 0.165149594771914570655D-02
    w( 66) = 0.169019899554346019117D-02
    w( 67) = 0.172907689054461607168D-02
    w( 68) = 0.176812249885838886701D-02
    w( 69) = 0.180732881501808930079D-02
    w( 70) = 0.184668895851282540913D-02
    w( 71) = 0.188619617015808475394D-02
    w( 72) = 0.192584380831993546204D-02
    w( 73) = 0.196562534503150547732D-02
    w( 74) = 0.200553436203751169944D-02
    w( 75) = 0.204556454679958293446D-02
    w( 76) = 0.208570968849203942640D-02
    w( 77) = 0.212596367401472533045D-02
    w( 78) = 0.216632048404649142727D-02
    w( 79) = 0.220677418916003329194D-02
    w( 80) = 0.224731894601603393082D-02
    w( 81) = 0.228794899365195972378D-02
    w( 82) = 0.232865864987842738864D-02
    w( 83) = 0.236944230779380495146D-02
    w( 84) = 0.241029443242563417382D-02
    w( 85) = 0.245120955750556483923D-02
    w( 86) = 0.249218228238276930060D-02
    w( 87) = 0.253320726907925325750D-02
    w( 88) = 0.257427923948908888092D-02
    w( 89) = 0.261539297272236109225D-02
    w( 90) = 0.265654330259352828314D-02
    w( 91) = 0.269772511525294586667D-02
    w( 92) = 0.273893334695947541201D-02
    w( 93) = 0.278016298199139435045D-02
    w( 94) = 0.282140905069222207923D-02
    w( 95) = 0.286266662764757868253D-02
    w( 96) = 0.290393082998878368175D-02
    w( 97) = 0.294519681581857582284D-02
    w( 98) = 0.298645978275408290247D-02
    w( 99) = 0.302771496658198544480D-02
    w(100) = 0.306895764002069252174D-02
    w(101) = 0.311018311158427546158D-02
    w(102) = 0.315138672454287935858D-02
    w(103) = 0.319256385597434736790D-02
    w(104) = 0.323370991590184336368D-02
    w(105) = 0.327482034651233969564D-02
    w(106) = 0.331589062145094394706D-02
    w(107) = 0.335691624518616761342D-02
    w(108) = 0.339789275244138669739D-02
    w(109) = 0.343881570768790591876D-02
    w(110) = 0.347968070469521146972D-02
    w(111) = 0.352048336613417922682D-02
    w(112) = 0.356121934322919357659D-02
    w(113) = 0.360188431545532431869D-02
    w(114) = 0.364247399027690353194D-02
    w(115) = 0.368298410292403911967D-02
    w(116) = 0.372341041620379550870D-02
    w(117) = 0.376374872034296338241D-02
    w(118) = 0.380399483285952829161D-02
    w(119) = 0.384414459846013158917D-02
    w(120) = 0.388419388896099560998D-02
    w(121) = 0.392413860322995774660D-02
    w(122) = 0.396397466714742455513D-02
    w(123) = 0.400369803358421688562D-02
    w(124) = 0.404330468239442998549D-02
    w(125) = 0.408279062042157838350D-02
    w(126) = 0.412215188151643401528D-02
    w(127) = 0.416138452656509745764D-02
    w(128) = 0.420048464352596631772D-02
    w(129) = 0.423944834747438184434D-02
    w(130) = 0.427827178065384480959D-02
    w(131) = 0.431695111253279479928D-02
    w(132) = 0.435548253986604343679D-02
    w(133) = 0.439386228676004195260D-02
    w(134) = 0.443208660474124713206D-02
    w(135) = 0.447015177282692726900D-02
    w(136) = 0.450805409759782158001D-02
    w(137) = 0.454578991327213285488D-02
    w(138) = 0.458335558178039420335D-02
    w(139) = 0.462074749284080687482D-02
    w(140) = 0.465796206403469754658D-02
    w(141) = 0.469499574088179046532D-02
    w(142) = 0.473184499691503264714D-02
    w(143) = 0.476850633375474925263D-02
    w(144) = 0.480497628118194150483D-02
    w(145) = 0.484125139721057135214D-02
    w(146) = 0.487732826815870573054D-02
    w(147) = 0.491320350871841897367D-02
    w(148) = 0.494887376202437487201D-02
    w(149) = 0.498433569972103029914D-02
    w(150) = 0.501958602202842039909D-02
    w(151) = 0.505462145780650125058D-02
    w(152) = 0.508943876461803986674D-02
    w(153) = 0.512403472879005351831D-02
    w(154) = 0.515840616547381084096D-02
    w(155) = 0.519254991870341614863D-02
    w(156) = 0.522646286145300596306D-02
    w(157) = 0.526014189569259311205D-02
    w(158) = 0.529358395244259896547D-02
    w(159) = 0.532678599182711857974D-02
    w(160) = 0.535974500312596681161D-02
    w(161) = 0.539245800482555593606D-02
    w(162) = 0.542492204466865704951D-02
    w(163) = 0.545713419970309863995D-02
    w(164) = 0.548909157632945623482D-02
    w(165) = 0.552079131034778706457D-02
    w(166) = 0.555223056700346326850D-02
    w(167) = 0.558340654103215637610D-02
    w(168) = 0.561431645670402467678D-02
    w(169) = 0.564495756786715368885D-02
    w(170) = 0.567532715799029830087D-02
    w(171) = 0.570542254020497332312D-02
    w(172) = 0.573524105734693719020D-02
    w(173) = 0.576478008199711142954D-02
    w(174) = 0.579403701652197628421D-02
    w(175) = 0.582300929311348057702D-02
    w(176) = 0.585169437382850155033D-02
    w(177) = 0.588008975062788803205D-02
    w(178) = 0.590819294541511788161D-02
    w(179) = 0.593600151007459827614D-02
    w(180) = 0.596351302650963502011D-02
    w(181) = 0.599072510668009471472D-02
    w(182) = 0.601763539263978131522D-02
    w(183) = 0.604424155657354634589D-02
    w(184) = 0.607054130083414983949D-02
    w(185) = 0.609653235797888692923D-02
    w(186) = 0.612221249080599294931D-02
    w(187) = 0.614757949239083790214D-02
    w(188) = 0.617263118612191922727D-02
    w(189) = 0.619736542573665996342D-02
    w(190) = 0.622178009535701763157D-02
    w(191) = 0.624587310952490748541D-02
    w(192) = 0.626964241323744217671D-02
    w(193) = 0.629308598198198836688D-02
    w(194) = 0.631620182177103938227D-02
    w(195) = 0.633898796917690165912D-02
    w(196) = 0.636144249136619145314D-02
    w(197) = 0.638356348613413709795D-02
    w(198) = 0.640534908193868098342D-02
    w(199) = 0.642679743793437438922D-02
    w(200) = 0.644790674400605734710D-02
    w(201) = 0.646867522080231481688D-02
    w(202) = 0.648910111976869964292D-02
    w(203) = 0.650918272318071200827D-02
    w(204) = 0.652891834417652442012D-02
    w(205) = 0.654830632678944064054D-02
    w(206) = 0.656734504598007641819D-02
    w(207) = 0.658603290766824937794D-02
    w(208) = 0.660436834876456498276D-02
    w(209) = 0.662234983720168509457D-02
    w(210) = 0.663997587196526532519D-02
    w(211) = 0.665724498312454708217D-02
    w(212) = 0.667415573186258997654D-02
    w(213) = 0.669070671050613006584D-02
    w(214) = 0.670689654255504925648D-02
    w(215) = 0.672272388271144108036D-02
    w(216) = 0.673818741690825799086D-02
    w(217) = 0.675328586233752529078D-02
    w(218) = 0.676801796747810680683D-02
    w(219) = 0.678238251212300746082D-02
    w(220) = 0.679637830740619795480D-02
    w(221) = 0.681000419582894688374D-02
    w(222) = 0.682325905128564571420D-02
    w(223) = 0.683614177908911221841D-02
    w(224) = 0.684865131599535812903D-02
    w(225) = 0.686078663022780697951D-02
    w(226) = 0.687254672150094831613D-02
    w(227) = 0.688393062104341470995D-02
    w(228) = 0.689493739162046825872D-02
    w(229) = 0.690556612755588354803D-02
    w(230) = 0.691581595475321433825D-02
    w(231) = 0.692568603071643155621D-02
    w(232) = 0.693517554456992049848D-02
    w(233) = 0.694428371707782549438D-02
    w(234) = 0.695300980066273063177D-02
    w(235) = 0.696135307942366551493D-02
    w(236) = 0.696931286915342540213D-02
    w(237) = 0.697688851735519545845D-02
    w(238) = 0.698407940325846925786D-02
    w(239) = 0.699088493783425207545D-02
    w(240) = 0.699730456380953992594D-02
    w(241) = 0.700333775568106572820D-02
    w(242) = 0.700898401972830440494D-02
    w(243) = 0.701424289402572916425D-02
    w(244) = 0.701911394845431165171D-02
    w(245) = 0.702359678471225911031D-02
    w(246) = 0.702769103632498213858D-02
    w(247) = 0.703139636865428709508D-02
    w(248) = 0.703471247890678765907D-02
    w(249) = 0.703763909614153052319D-02
    w(250) = 0.704017598127683066242D-02
    w(251) = 0.704232292709631209597D-02
    w(252) = 0.704407975825415053266D-02
    w(253) = 0.704544633127951476780D-02
    w(254) = 0.704642253458020417748D-02
    w(255) = 0.704700828844548013730D-02
    w(256) = 0.704720354504808967346D-02
    w(257) = 0.704700828844548013730D-02
    w(258) = 0.704642253458020417748D-02
    w(259) = 0.704544633127951476780D-02
    w(260) = 0.704407975825415053266D-02
    w(261) = 0.704232292709631209597D-02
    w(262) = 0.704017598127683066242D-02
    w(263) = 0.703763909614153052319D-02
    w(264) = 0.703471247890678765907D-02
    w(265) = 0.703139636865428709508D-02
    w(266) = 0.702769103632498213858D-02
    w(267) = 0.702359678471225911031D-02
    w(268) = 0.701911394845431165171D-02
    w(269) = 0.701424289402572916425D-02
    w(270) = 0.700898401972830440494D-02
    w(271) = 0.700333775568106572820D-02
    w(272) = 0.699730456380953992594D-02
    w(273) = 0.699088493783425207545D-02
    w(274) = 0.698407940325846925786D-02
    w(275) = 0.697688851735519545845D-02
    w(276) = 0.696931286915342540213D-02
    w(277) = 0.696135307942366551493D-02
    w(278) = 0.695300980066273063177D-02
    w(279) = 0.694428371707782549438D-02
    w(280) = 0.693517554456992049848D-02
    w(281) = 0.692568603071643155621D-02
    w(282) = 0.691581595475321433825D-02
    w(283) = 0.690556612755588354803D-02
    w(284) = 0.689493739162046825872D-02
    w(285) = 0.688393062104341470995D-02
    w(286) = 0.687254672150094831613D-02
    w(287) = 0.686078663022780697951D-02
    w(288) = 0.684865131599535812903D-02
    w(289) = 0.683614177908911221841D-02
    w(290) = 0.682325905128564571420D-02
    w(291) = 0.681000419582894688374D-02
    w(292) = 0.679637830740619795480D-02
    w(293) = 0.678238251212300746082D-02
    w(294) = 0.676801796747810680683D-02
    w(295) = 0.675328586233752529078D-02
    w(296) = 0.673818741690825799086D-02
    w(297) = 0.672272388271144108036D-02
    w(298) = 0.670689654255504925648D-02
    w(299) = 0.669070671050613006584D-02
    w(300) = 0.667415573186258997654D-02
    w(301) = 0.665724498312454708217D-02
    w(302) = 0.663997587196526532519D-02
    w(303) = 0.662234983720168509457D-02
    w(304) = 0.660436834876456498276D-02
    w(305) = 0.658603290766824937794D-02
    w(306) = 0.656734504598007641819D-02
    w(307) = 0.654830632678944064054D-02
    w(308) = 0.652891834417652442012D-02
    w(309) = 0.650918272318071200827D-02
    w(310) = 0.648910111976869964292D-02
    w(311) = 0.646867522080231481688D-02
    w(312) = 0.644790674400605734710D-02
    w(313) = 0.642679743793437438922D-02
    w(314) = 0.640534908193868098342D-02
    w(315) = 0.638356348613413709795D-02
    w(316) = 0.636144249136619145314D-02
    w(317) = 0.633898796917690165912D-02
    w(318) = 0.631620182177103938227D-02
    w(319) = 0.629308598198198836688D-02
    w(320) = 0.626964241323744217671D-02
    w(321) = 0.624587310952490748541D-02
    w(322) = 0.622178009535701763157D-02
    w(323) = 0.619736542573665996342D-02
    w(324) = 0.617263118612191922727D-02
    w(325) = 0.614757949239083790214D-02
    w(326) = 0.612221249080599294931D-02
    w(327) = 0.609653235797888692923D-02
    w(328) = 0.607054130083414983949D-02
    w(329) = 0.604424155657354634589D-02
    w(330) = 0.601763539263978131522D-02
    w(331) = 0.599072510668009471472D-02
    w(332) = 0.596351302650963502011D-02
    w(333) = 0.593600151007459827614D-02
    w(334) = 0.590819294541511788161D-02
    w(335) = 0.588008975062788803205D-02
    w(336) = 0.585169437382850155033D-02
    w(337) = 0.582300929311348057702D-02
    w(338) = 0.579403701652197628421D-02
    w(339) = 0.576478008199711142954D-02
    w(340) = 0.573524105734693719020D-02
    w(341) = 0.570542254020497332312D-02
    w(342) = 0.567532715799029830087D-02
    w(343) = 0.564495756786715368885D-02
    w(344) = 0.561431645670402467678D-02
    w(345) = 0.558340654103215637610D-02
    w(346) = 0.555223056700346326850D-02
    w(347) = 0.552079131034778706457D-02
    w(348) = 0.548909157632945623482D-02
    w(349) = 0.545713419970309863995D-02
    w(350) = 0.542492204466865704951D-02
    w(351) = 0.539245800482555593606D-02
    w(352) = 0.535974500312596681161D-02
    w(353) = 0.532678599182711857974D-02
    w(354) = 0.529358395244259896547D-02
    w(355) = 0.526014189569259311205D-02
    w(356) = 0.522646286145300596306D-02
    w(357) = 0.519254991870341614863D-02
    w(358) = 0.515840616547381084096D-02
    w(359) = 0.512403472879005351831D-02
    w(360) = 0.508943876461803986674D-02
    w(361) = 0.505462145780650125058D-02
    w(362) = 0.501958602202842039909D-02
    w(363) = 0.498433569972103029914D-02
    w(364) = 0.494887376202437487201D-02
    w(365) = 0.491320350871841897367D-02
    w(366) = 0.487732826815870573054D-02
    w(367) = 0.484125139721057135214D-02
    w(368) = 0.480497628118194150483D-02
    w(369) = 0.476850633375474925263D-02
    w(370) = 0.473184499691503264714D-02
    w(371) = 0.469499574088179046532D-02
    w(372) = 0.465796206403469754658D-02
    w(373) = 0.462074749284080687482D-02
    w(374) = 0.458335558178039420335D-02
    w(375) = 0.454578991327213285488D-02
    w(376) = 0.450805409759782158001D-02
    w(377) = 0.447015177282692726900D-02
    w(378) = 0.443208660474124713206D-02
    w(379) = 0.439386228676004195260D-02
    w(380) = 0.435548253986604343679D-02
    w(381) = 0.431695111253279479928D-02
    w(382) = 0.427827178065384480959D-02
    w(383) = 0.423944834747438184434D-02
    w(384) = 0.420048464352596631772D-02
    w(385) = 0.416138452656509745764D-02
    w(386) = 0.412215188151643401528D-02
    w(387) = 0.408279062042157838350D-02
    w(388) = 0.404330468239442998549D-02
    w(389) = 0.400369803358421688562D-02
    w(390) = 0.396397466714742455513D-02
    w(391) = 0.392413860322995774660D-02
    w(392) = 0.388419388896099560998D-02
    w(393) = 0.384414459846013158917D-02
    w(394) = 0.380399483285952829161D-02
    w(395) = 0.376374872034296338241D-02
    w(396) = 0.372341041620379550870D-02
    w(397) = 0.368298410292403911967D-02
    w(398) = 0.364247399027690353194D-02
    w(399) = 0.360188431545532431869D-02
    w(400) = 0.356121934322919357659D-02
    w(401) = 0.352048336613417922682D-02
    w(402) = 0.347968070469521146972D-02
    w(403) = 0.343881570768790591876D-02
    w(404) = 0.339789275244138669739D-02
    w(405) = 0.335691624518616761342D-02
    w(406) = 0.331589062145094394706D-02
    w(407) = 0.327482034651233969564D-02
    w(408) = 0.323370991590184336368D-02
    w(409) = 0.319256385597434736790D-02
    w(410) = 0.315138672454287935858D-02
    w(411) = 0.311018311158427546158D-02
    w(412) = 0.306895764002069252174D-02
    w(413) = 0.302771496658198544480D-02
    w(414) = 0.298645978275408290247D-02
    w(415) = 0.294519681581857582284D-02
    w(416) = 0.290393082998878368175D-02
    w(417) = 0.286266662764757868253D-02
    w(418) = 0.282140905069222207923D-02
    w(419) = 0.278016298199139435045D-02
    w(420) = 0.273893334695947541201D-02
    w(421) = 0.269772511525294586667D-02
    w(422) = 0.265654330259352828314D-02
    w(423) = 0.261539297272236109225D-02
    w(424) = 0.257427923948908888092D-02
    w(425) = 0.253320726907925325750D-02
    w(426) = 0.249218228238276930060D-02
    w(427) = 0.245120955750556483923D-02
    w(428) = 0.241029443242563417382D-02
    w(429) = 0.236944230779380495146D-02
    w(430) = 0.232865864987842738864D-02
    w(431) = 0.228794899365195972378D-02
    w(432) = 0.224731894601603393082D-02
    w(433) = 0.220677418916003329194D-02
    w(434) = 0.216632048404649142727D-02
    w(435) = 0.212596367401472533045D-02
    w(436) = 0.208570968849203942640D-02
    w(437) = 0.204556454679958293446D-02
    w(438) = 0.200553436203751169944D-02
    w(439) = 0.196562534503150547732D-02
    w(440) = 0.192584380831993546204D-02
    w(441) = 0.188619617015808475394D-02
    w(442) = 0.184668895851282540913D-02
    w(443) = 0.180732881501808930079D-02
    w(444) = 0.176812249885838886701D-02
    w(445) = 0.172907689054461607168D-02
    w(446) = 0.169019899554346019117D-02
    w(447) = 0.165149594771914570655D-02
    w(448) = 0.161297501254393423070D-02
    w(449) = 0.157464359003212166189D-02
    w(450) = 0.153650921735128916170D-02
    w(451) = 0.149857957106456636214D-02
    w(452) = 0.146086246895890987689D-02
    w(453) = 0.142336587141720519900D-02
    w(454) = 0.138609788229672549700D-02
    w(455) = 0.134906674928353113127D-02
    w(456) = 0.131228086370221478128D-02
    w(457) = 0.127574875977346947345D-02
    w(458) = 0.123947911332878396534D-02
    w(459) = 0.120348074001265964881D-02
    w(460) = 0.116776259302858043685D-02
    w(461) = 0.113233376051597664917D-02
    w(462) = 0.109720346268191941940D-02
    w(463) = 0.106238104885340071375D-02
    w(464) = 0.102787599466367326179D-02
    w(465) = 0.993697899638760857945D-03
    w(466) = 0.959856485506936206261D-03
    w(467) = 0.926361595613111283368D-03
    w(468) = 0.893223195879324912340D-03
    w(469) = 0.860451377808527848128D-03
    w(470) = 0.828056364077226302608D-03
    w(471) = 0.796048517297550871506D-03
    w(472) = 0.764438352543882784191D-03
    w(473) = 0.733236554224767912055D-03
    w(474) = 0.702453997827572321358D-03
    w(475) = 0.672101776960108194646D-03
    w(476) = 0.642191235948505088403D-03
    w(477) = 0.612734008012225209294D-03
    w(478) = 0.583742058714979703847D-03
    w(479) = 0.555227733977307579715D-03
    w(480) = 0.527203811431658386125D-03
    w(481) = 0.499683553312800484519D-03
    w(482) = 0.472680758429262691232D-03
    w(483) = 0.446209810101403247488D-03
    w(484) = 0.420285716355361231823D-03
    w(485) = 0.394924138246873704434D-03
    w(486) = 0.370141402122251665232D-03
    w(487) = 0.345954492129903871350D-03
    w(488) = 0.322381020652862389664D-03
    w(489) = 0.299439176850911730874D-03
    w(490) = 0.277147657465187357459D-03
    w(491) = 0.255525589595236862014D-03
    w(492) = 0.234592462123925204879D-03
    w(493) = 0.214368090034216937149D-03
    w(494) = 0.194872642236641146532D-03
    w(495) = 0.176126765545083195474D-03
    w(496) = 0.158151830411132242924D-03
    w(497) = 0.140970302204104791413D-03
    w(498) = 0.124606200241498368482D-03
    w(499) = 0.109085545645741522051D-03
    w(500) = 0.944366322532705527066D-04
    w(501) = 0.806899228014035293851D-04
    w(502) = 0.678774554733972416227D-04
    w(503) = 0.560319507856164252140D-04
    w(504) = 0.451863674126296143105D-04
    w(505) = 0.353751372055189588628D-04
    w(506) = 0.266376412339000901358D-04
    w(507) = 0.190213681905875816679D-04
    w(508) = 0.125792781889592743525D-04
    w(509) = 0.736624069102321668857D-05
    w(510) = 0.345456507169149134898D-05
    w(511) = 0.945715933950007048827D-06

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PATTERSON_LOOKUP_WEIGHTS - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 7, 15, 31, 63, 127, 255 or 511.'
    stop

  end if

  return
end
subroutine patterson_lookup_weights_np ( order, np, p, w )

!*****************************************************************************80
!
!! PATTERSON_LOOKUP_WEIGHTS_NP sets weights for a Patterson rule.
!
!  Discussion:
!
!    The zeroth rule, of order 1, is the standard Legendre rule.
!
!    The first rule, of order 3, is the standard Legendre rule.
!
!    The second rule, of order 7, includes the abscissas of the previous
!    rule.
!
!    Each subsequent rule is nested in a similar way.  
!    Rules are available of orders 1, 3, 7, 15, 31, 63, 127, 255 and 511.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    Leval values are 1, 3, 7, 15, 31, 63, 127, 255 and 511.
!
!    Input, integer ( kind = 4 ) NP, the number of parameters.
!
!    Input, real ( kind = 8 ) P(NP), the parameters.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) order

  real ( kind = 8 ) p(*)
  real ( kind = 8 ) w(order)

  call patterson_lookup_weights ( order, w )

  return
end
subroutine point_radial_tol_unique_count ( m, n, a, tol, seed, unique_num )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_COUNT counts the tolerably unique points.
!
!  Discussion:
!
!    The input data is an M x N array A, representing the M-dimensional
!    coordinates of N points.
!
!    The output is the number of tolerably unique points in the list.
!
!    This program performs the same task as POINT_TOL_UNIQUE_COUNT.
!    But that program is guaranteed to use N^2 comparisons.
!
!    It is hoped that this function, on the other hand, will tend
!    to use O(N) comparisons after an O(NLog(N)) sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of tolerably
!    unique points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tol
  logical unique(n)
  integer ( kind = 4 ) unique_num
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) w_sum
  real ( kind = 8 ) z(m)

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if
!
!  Assign a base point Z randomly in the convex hull.
!
  call r8vec_uniform_01 ( n, seed, w )
  w_sum = sum ( w(1:n) )

  w(1:n) = w(1:n) / w_sum

  z = matmul ( a(1:m,1:n), w(1:n) )
!
!  Compute the radial distance R of each point to Z.
!
  do j = 1, n
    r(j) = sqrt ( sum ( ( a(1:m,j) - z(1:m) )**2 ) )
  end do
!
!  Implicitly sort the R array.
!
  call r8vec_sort_heap_index_a ( n, r, indx )
!
!  To determine if a point I is tolerably unique, we only have to check
!  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
!
  unique_num = 0
  unique(1:n) = .true.

  do i = 1, n

    if ( unique(indx(i)) ) then
!
!  Point INDX(I) is unique, in that no earlier point is near it.
!
      unique_num = unique_num + 1
!
!  Look for later points which are close to point INDX(I)
!  in terms of R.
!
      hi = i

      do while ( hi < n )

        if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
          exit
        end if
        hi = hi + 1

      end do
!
!  Points INDX(I+1) through INDX(HI) have an R value close to
!  point INDX(I).  Are they truly close to point INDEX(I)?
!
      do j = i + 1, hi

        if ( unique(indx(j)) ) then
          dist = sqrt ( sum ( ( a(1:m,indx(i)) - a(1:m,indx(j)) )**2 ) )

          if ( dist <= tol ) then
            unique(indx(j)) = .false.
          end if
        end if

      end do

    end if
  end do

  return
end
subroutine point_radial_tol_unique_count_old ( m, n1, a1, n2, a2, tol, seed, &
  unique_num1, unique_num2 )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_COUNT counts the tolerably unique points.
!
!  Discussion:
!
!    The input data includes an M x N1 array A1 and an M x N2 array A2,
!    representing the M-dimensional coordinates of a set of N1
!    "permanent" points and N2 "provisional" points.
!
!    This is an "incremental" version of POINT_RADIAL_TOL_UNIQUE_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N1, the number of permanent points.
!
!    Input, real ( kind = 8 ) A1(M,N1), the permanent points.
!
!    Input, integer ( kind = 4 ) N2, the number of provisional points.
!
!    Input, real ( kind = 8 ) A2(M,N2), the provisional points.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM1, the number of tolerably
!    unique permanent points.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM2, the number of tolerably
!    unique points when the temporary points are included.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a1(m,n1)
  real ( kind = 8 ) a2(m,n2)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n1+n2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(n1+n2)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num1
  integer ( kind = 4 ) unique_num2
  logical unique1(n1)
  logical unique2(n2)
  real ( kind = 8 ) w(n1)
  real ( kind = 8 ) z(m)

  n = n1 + n2

  if ( n <= 0 ) then
    unique_num1 = 0
    unique_num2 = 0
    return
  end if
!
!  Assign a base point Z randomly in the convex hull of the permanent points.
!
  call r8vec_uniform_01 ( n1, seed, w )

  do i = 1, m
    z(i) = dot_product ( a1(i,1:n1), w(1:n1) ) / sum ( w(1:n1) )
  end do
!
!  STEP 1:
!  Compare PERMANENT POINTS to PERMANENT POINTS.
!

!
!  Compute the radial distance R of each permanent point to Z.
!
  do j = 1, n1
    r(j) = sqrt ( sum ( ( a1(1:m,j) - z(1:m) )**2 ) )
  end do
!
!  Implicitly sort the R array.
!
  call r8vec_sort_heap_index_a ( n1, r, indx )
!
!  To determine if a point I is tolerably unique, we only have to check
!  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
!
  unique_num1 = 0
  unique1(1:n1) = .true.

  do i = 1, n1

    if ( indx(i) <= n1 ) then

      if ( unique1(indx(i)) ) then
!
!  Point INDX(I) is unique, in that no earlier point is near it.
!
        unique_num1 = unique_num1 + 1
!
!  Look for later points which are close to point INDX(I) in terms of R.
!
        hi = i

        do while ( hi < n1 )
          if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
            exit
          end if
          hi = hi + 1
        end do
!
!  Points INDX(I+1) through INDX(HI) have an R value close to point INDX(I).
!  Are they permanent points?
!  Did we think they were unique?
!  Are they truly close to point INDEX(I)?
!
        do j = i + 1, hi
          if ( indx(j) <= n1 ) then
            if ( unique1(indx(j)) ) then
              dist = sqrt ( sum ( ( a1(1:m,indx(i)) - a1(1:m,indx(j)) )**2 ) )
              if ( dist <= tol ) then
                unique1(indx(j)) = .false.
              end if
            end if
          end if
        end do

      end if

    end if

  end do
!
!  STEP 2:
!  Compare TEMPORARY POINTS to PERMANENT POINTS.
!
  do j = 1, n2
    r(j+n1) = sqrt ( sum ( ( a2(1:m,j) - z(1:m) )**2 ) )
  end do

  call r8vec_sort_heap_index_a ( n, r, indx )

  unique2(1:n2) = .true.

  do i = 1, n

    if ( indx(i) <= n1 ) then

      if ( unique1(indx(i)) ) then
!
!  Point INDX(I) is a unique permanent point.
!  Look for later points which are close to point INDX(I) in terms of R.
!
        hi = i

        do while ( hi < n )
          if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
            exit
          end if
          hi = hi + 1
        end do
!
!  Points INDX(I+1) through INDX(HI) have an R value close to point INDX(I).
!  Are they temporary points?
!  Did we think they were unique?
!  Are they truly close to permanent point INDEX(I)?
!
        do j = i + 1, hi
          if ( n1 < indx(j) ) then
            if ( unique2(indx(j)) ) then
              dist = sqrt ( sum ( ( a1(1:m,indx(i)) &
                                  - a2(1:m,indx(j)-n1) )**2 ) )
              if ( dist <= tol ) then
                unique2(indx(j)) = .false.
             end if
            end if
          end if
        end do

      end if

    end if

  end do
!
!  STEP 3:
!  Compare TEMPORARY POINTS to TEMPORARY POINTS.
!  (but also look backwards at nearby permanent points).
!
  unique_num2 = unique_num1

  do i = 1, n

    if ( n1 < indx(i) ) then

      if ( unique2(indx(i)) ) then
!
!  Point INDX(I) is unique, in that no earlier point is near it.
!
        unique_num2 = unique_num2 + 1
!
!  Look for later points which are close to point INDX(I) in terms of R.
!
        hi = i

        do while ( hi < n )
          if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
            exit
          end if
          hi = hi + 1
        end do
!
!  First check whether there are any nearby permanent points.
!  If so, cancel this "unique" temporary point.
!  This seemingly needless check occurs because normally we only
!  check "to the right" of a point as we proceed, which would
!  normally be enough to catch a pair of close points.  But if
!  the temporary point occurs to the LEFT of the permanent point,
!  then when we check to the right of all the permanent points, we
!  don't see the temporary point.  When we check to the right of
!  the temporary point, we see the permanent point, but we prefer
!  to mark the permanent point as unique, not the temporary one!
!
        do j = i + 1, hi

          if ( indx(j) <= n1 ) then

            if ( unique1(indx(j)) ) then

              dist = sqrt ( sum ( ( a2(1:m,indx(i)-n1) &
                                  - a1(1:m,indx(j)) )**2 ) )

              if ( dist <= tol ) then
                unique2(indx(i)) = .false.
                unique_num2 = unique_num2 - 1
                hi = i
                exit
              end if

            end if

          end if
        end do
!
!  If INDX(I) is still counted as tolerably unique, we now look for
!  nearby temporary points that can be marked as nonunique.
!
        do j = i + 1, hi
          if ( n1 < indx(j) ) then
            if ( unique2(indx(j)) ) then
              dist = sqrt ( sum ( ( a2(1:m,indx(i)-n1) &
                                  - a2(1:m,indx(j)-n1) )**2 ) )
              if ( dist <= tol ) then
                unique2(indx(j)) = .false.
              end if
            end if
          end if
        end do

      end if

    end if

  end do

  return
end
subroutine point_radial_tol_unique_count_inc1 ( m, n1, a1, tol, seed, z, &
  r1, indx1, unique1, unique_num1 )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_COUNT_INC1 counts the tolerably unique points.
!
!  Discussion:
!
!    The input data includes an M x N1 array A1 of a set of N1
!    "permanent" points and N2 "temporary" points.
!
!    This is an two step version of POINT_RADIAL_TOL_UNIQUE_COUNT_INC.
!
!    This means that we want to identify the tolerably unique points
!    among the permanent points before processing the temporary points.
!
!    If many sets of temporary data are considered, this function will
!    do a lot of unnecessary work resorting the permanent data; it would
!    be possible to avoid repetitions of that work at the expense of saving
!    various work vectors.  This function accepts the overhead of the
!    repeated calculations for the benefit of only having to "remember"
!    the number of unique points discovered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N1, the number of permanent points.
!
!    Input, real ( kind = 8 ) A1(M,N1), the permanent points.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) Z(M), a random base vector used to
!    linearly sort the data.
!
!    Output, real ( kind = 8 ) R1(N1), the scalar values assigned to
!    the data for sorting.
!
!    Output, integer ( kind = 4 ) INDX1(N1), the ascending sort index
!    for A1.
!
!    Output, logical UNIQUE1(N1), is TRUE for each unique permanent point.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM1, the number of tolerably
!    unique permanent points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1

  real ( kind = 8 ) a1(m,n1)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx1(n1)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k1
  real ( kind = 8 ) r1(n1)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num1
  logical unique1(n1)
  real ( kind = 8 ) w(n1)
  real ( kind = 8 ) z(m)
!
!  Assign a base point Z randomly in the convex hull of the permanent points.
!
  call r8vec_uniform_01 ( n1, seed, w )

  do i = 1, m
    z(i) = dot_product ( a1(i,1:n1), w(1:n1) ) / sum ( w(1:n1) )
  end do
!
!  Initialize the permanent point data.
!
  do j1 = 1, n1
    r1(j1) = sqrt ( sum ( ( a1(1:m,j1) - z(1:m) )**2 ) )
  end do

  call r8vec_sort_heap_index_a ( n1, r1, indx1 )

  unique_num1 = 0
  unique1(1:n1) = .true.
!
!  STEP 1:
!  Compare PERMANENT POINTS to PERMANENT POINTS.
!
  do j1 = 1, n1

    if ( unique1(indx1(j1)) ) then

      unique_num1 = unique_num1 + 1

      hi = j1

      do while ( hi < n1 )
        if ( r1(indx1(j1)) + tol < r1(indx1(hi+1)) ) then
          exit
        end if
        hi = hi + 1
      end do

      do k1 = j1 + 1, hi
        if ( unique1(indx1(k1)) ) then
          dist = sqrt ( sum ( ( a1(1:m,indx1(j1)) - a1(1:m,indx1(k1)) )**2 ) )
          if ( dist <= tol ) then
            unique1(indx1(k1)) = .false.
          end if
        end if
      end do

    end if

  end do

  return
end
subroutine point_radial_tol_unique_count_inc2 ( m, n1, a1, n2, a2, tol, &
  z, r1, indx1, unique1, unique_num2 )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_COUNT_INC2 counts the tolerably unique points.
!
!  Discussion:
!
!    The input data includes an M x N1 array A1 and an M x N2 array A2,
!    representing the M-dimensional coordinates of a set of N1
!    "permanent" points and N2 "temporary" points.
!
!    This is a two step version of POINT_RADIAL_TOL_UNIQUE_COUNT_INC.
!
!    This means that we want to identify the tolerably unique points
!    among the permanent points before processing the temporary points.
!
!    If many sets of temporary data are considered, this function will
!    do a lot of unnecessary work resorting the permanent data; it would
!    be possible to avoid repetitions of that work at the expense of saving
!    various work vectors.  This function accepts the overhead of the
!    repeated calculations for the benefit of only having to "remember"
!    the number of unique points discovered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N1, the number of permanent points.
!
!    Input, real ( kind = 8 ) A1(M,N1), the permanent points.
!
!    Input, integer ( kind = 4 ) N2, the number of temporary points.
!
!    Input, real ( kind = 8 ) A2(M,N2), the temporary points.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input, real ( kind = 8 ) Z(M), a random base vector used to
!    linearly sort the data.
!
!    Input, real ( kind = 8 ) R1(N1), the scalar values assigned to
!    the data for sorting.
!
!    Input, integer ( kind = 4 ) INDX1(N1), the ascending sort index
!    for A1.
!
!    Input, logical UNIQUE1(N1), is TRUE for each unique permanent point.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM2, the number of additional
!    tolerably unique points if the temporary points are included.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a1(m,n1)
  real ( kind = 8 ) a2(m,n2)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx1(n1)
  integer ( kind = 4 ) indx2(n2)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2_hi
  integer ( kind = 4 ) j2_lo
  integer ( kind = 4 ) k2
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
  real ( kind = 8 ) r1(n1)
  real ( kind = 8 ) r2(n2)
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num2
  logical unique1(n1)
  logical unique2(n2)
  real ( kind = 8 ) z(m)
!
!  Initialize the temporary point data.
!
  do j2 = 1, n2
    r2(j2) = sqrt ( sum ( ( a2(1:m,j2) - z(1:m) )**2 ) )
  end do

  call r8vec_sort_heap_index_a ( n2, r2, indx2 )

  unique2(1:n2) = .true.

  unique_num2 = 0
!
!  STEP 2:
!  Use PERMANENT points to eliminate TEMPORARY points.
!
  do j1 = 1, n1

    if ( unique1(indx1(j1)) ) then

      r_lo = r1(indx1(j1)) - tol
      r_hi = r1(indx1(j1)) + tol

      call r8vec_index_sorted_range ( n2, r2, indx2, r_lo, r_hi, &
        j2_lo, j2_hi )

      do j2 = j2_lo, j2_hi
        if ( unique2(indx2(j2)) ) then
          dist = sqrt ( sum ( ( a1(1:m,indx1(j1)) &
                              - a2(1:m,indx2(j2)) )**2 ) )
          if ( dist <= tol ) then
            unique2(indx2(j2)) = .false.
          end if
        end if
      end do

    end if

  end do
!
!  STEP 3:
!  Use TEMPORARY points to eliminate TEMPORARY points.
!
  do j2 = 1, n2

    if ( unique2(indx2(j2)) ) then

      unique_num2 = unique_num2 + 1

      hi = j2

      do while ( hi < n2 )
        if ( r2(indx2(j2)) + tol < r2(indx2(hi+1)) ) then
          exit
        end if
        hi = hi + 1
      end do

      do k2 = j2 + 1, hi
        if ( unique2(indx2(k2)) ) then
          dist = sqrt ( sum ( ( a2(1:m,indx2(j2)) &
                              - a2(1:m,indx2(k2)) )**2 ) )
          if ( dist <= tol ) then
            unique2(indx2(k2)) = .false.
          end if
        end if
      end do

    end if

  end do

  return
end
subroutine point_radial_tol_unique_index ( m, n, a, tol, seed, unique_num, &
  undx, xdnu )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_INDEX indexes the tolerably unique points.
!
!  Discussion:
!
!    The input data is an M x N array A, representing the M-dimensional
!    coordinates of N points.
!
!    The output is:
!    * the number of tolerably unique points in the list
!    * the index, in the list of unique items, of the representatives
!      of each point
!    * the index, in A, of the tolerably unique representatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of tolerably
!    unique points.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the index, in A, of the
!    tolerably unique points.
!
!    Output, integer ( kind = 4 ) XDNU(N), the index, in UNDX, of the
!    tolerably unique point that "represents" this point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(unique_num)
  logical unique(n)
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) xdnu(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) w_sum
  real ( kind = 8 ) z(m)

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if
!
!  Assign a base point Z randomly in the convex hull.
!
  call r8vec_uniform_01 ( n, seed, w )
  w_sum = sum ( w(1:n) )
  w(1:n) = w(1:n) / w_sum

  z = matmul ( a(1:m,1:n), w(1:n) )
!
!  Compute the radial distance R of each point to Z.
!
  do j = 1, n
    r(j) = sqrt ( sum ( ( a(1:m,j) - z(1:m) )**2 ) )
  end do
!
!  Implicitly sort the R array.
!
  call r8vec_sort_heap_index_a ( n, r, indx )
!
!  To determine if a point I is tolerably unique, we only have to check
!  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
!
  unique_num = 0
  unique(1:n) = .true.

  do i = 1, n

    if ( unique(indx(i)) ) then
!
!  Point INDX(I) is unique, in that no earlier point is near it.
!
      unique_num = unique_num + 1
      xdnu(indx(i)) = unique_num
      undx(unique_num) = indx(i)
!
!  Look for later points which are close to point INDX(I)
!  in terms of R.
!
      hi = i

      do while ( hi < n )

        if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
          exit
        end if
        hi = hi + 1

      end do
!
!  Points INDX(I+1) through INDX(HI) have an R value close to
!  point INDX(I).  Are they truly close to point INDEX(I)?
!
      do j = i + 1, hi

        if ( unique(indx(j)) ) then
          dist = sqrt ( sum ( ( a(1:m,indx(i)) - a(1:m,indx(j)) )**2 ) )

          if ( dist <= tol ) then
            unique(indx(j)) = .false.
            xdnu(indx(j)) = xdnu(indx(i))
          end if
        end if

      end do

    end if
  end do

  return
end
subroutine point_radial_tol_unique_index_old ( m, n1, a1, n2, a2, tol, seed, &
  unique_num1, unique_num2, undx, xdnu )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_INDEX_OLD indexes the tolerably unique points.
!
!  Discussion:
!
!    The input data includes an M x N1 array A1 and an M x N2 array A2,
!    representing the M-dimensional coordinates of a set of N1
!    "permanent" points and N2 "provisional" points.
!
!    For notation, we use "A" to describe the M x (N1+N2) array that would be
!    formed by starting with A1 and appending A2.
!
!    The output is:
!    * the number of tolerably unique points in the list;
!    * the index, in the list of unique items, of the representatives
!      of each point;
!    * the index, in A, of the tolerably unique representatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N1, the number of permanent points.
!
!    Input, real ( kind = 8 ) A1(M,N1), the permanent points.
!
!    Input, integer ( kind = 4 ) N2, the number of provisional points.
!
!    Input, real ( kind = 8 ) A2(M,N2), the provisional points.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM1, UNIQUE_NUM2, the number
!    of tolerably unique points with just the permanent points, or with
!    the permanent points incremented with the provisional points.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM2), the index, in A, of the
!    tolerably unique points.
!
!    Output, integer ( kind = 4 ) XDNU(N1+N2), the index, in UNDX, of the
!    tolerably unique point that "represents" this point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a1(m,n1)
  real ( kind = 8 ) a2(m,n2)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n1+n2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(n1+n2)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tol
  logical unique(n1+n2)
  integer ( kind = 4 ) undx(n1+n2)
  integer ( kind = 4 ) unique_num1
  integer ( kind = 4 ) unique_num2
  real ( kind = 8 ) w(n1)
  integer ( kind = 4 ) xdnu(n1+n2)
  real ( kind = 8 ) z(m)

  n = n1 + n2

  if ( n <= 0 ) then
    unique_num1 = 0
    unique_num2 = 0
    return
  end if
!
!  Assign a base point Z randomly in the convex hull of the permanent points.
!
  call r8vec_uniform_01 ( n1, seed, w )

  do i = 1, m
    z(i) = dot_product ( a1(i,1:n1), w(1:n1) ) / sum ( w(1:n1) )
  end do
!
!  Compute the radial distance R of each point to Z.
!
  do j = 1, n1
    r(j) = sqrt ( sum ( ( a1(1:m,j) - z(1:m) )**2 ) )
  end do

  do j = 1, n2
    r(j+n1) = sqrt ( sum ( ( a2(1:m,j) - z(1:m) )**2 ) )
  end do
!
!  Implicitly sort the array.
!
  call r8vec_sort_heap_index_a ( n, r, indx )
!
!  To determine if a point I is tolerably unique, we only have to check
!  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
!
  unique_num1 = 0
  unique(1:n) = .true.
!
!  Process the N1 permanent points.
!
  do i = 1, n

    if ( indx(i) <= n1 ) then

      if ( unique(indx(i)) ) then
!
!  Point INDX(I) is unique, in that no earlier point is near it.
!
        unique_num1 = unique_num1 + 1
        xdnu(indx(i)) = unique_num1
        undx(unique_num1) = indx(i)
!
!  Look for later points which are close to point INDX(I) in terms of R.
!
        hi = i

        do while ( hi < n )
          if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
            exit
          end if
          hi = hi + 1
        end do
!
!  Points INDX(I+1) through INDX(HI) have an R value close to point INDX(I).
!  Are they truly close to point INDEX(I)?
!
        do j = i + 1, hi
          if ( unique(indx(j)) ) then
            if ( indx(j) <= n1 ) then
              dist = sqrt ( sum ( ( a1(1:m,indx(i)) - a1(1:m,indx(j)) )**2 ) )
            else
              dist = &
                sqrt ( sum ( ( a1(1:m,indx(i)) - a2(1:m,indx(j)-n1) )**2 ) )
            end if
            if ( dist <= tol ) then
              unique(indx(j)) = .false.
              xdnu(indx(j)) = xdnu(indx(i))
            end if
          end if
        end do

      end if

    end if

  end do
!
!  Process the N2 temporary points.
!  No temporary point is allowed to make a permanent point nonunique.
!
  unique_num2 = unique_num1

  do i = 1, n

    if ( n1 < indx(i) ) then

      if ( unique(indx(i)) ) then
!
!  Point INDX(I) is unique, in that no earlier point is near it.
!
        unique_num2 = unique_num2 + 1
        xdnu(indx(i)) = unique_num2
        undx(unique_num2) = indx(i)
!
!  Look for later points which are close to point INDX(I) in terms of R.
!
        hi = i

        do while ( hi < n )
          if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
            exit
          end if
          hi = hi + 1
        end do
!
!  First check whether there are any nearby permanent points.
!  If so, cancel this "unique" temporary point.
!  This seemingly needless check occurs because normally we only
!  check "to the right" of a point as we proceed, which would
!  normally be enough to catch a pair of close points.  But if
!  the temporary point occurs to the LEFT of the permanent point,
!  then when we check to the right of all the permanent points, we
!  don't see the temporary point.  When we check to the right of
!  the temporary point, we see the permanent point, but we prefer
!  to mark the permanent point as unique, not the temporary one!
!
        do j = i + 1, hi
          if ( unique(indx(j)) ) then
            if ( indx(j) <= n1 ) then
              dist = sqrt ( sum ( ( a2(1:m,indx(i)-n1) &
                                  - a1(1:m,indx(j)) )**2 ) )
            else
              dist = sqrt ( sum ( ( a2(1:m,indx(i)-n1) &
                                  - a2(1:m,indx(j)-n1) )**2 ) )
            end if
            if ( dist <= tol ) then
              if ( indx(j) <= n1 ) then
                unique(indx(i)) = .false.
                xdnu(indx(i)) = xdnu(indx(j))
                undx(unique_num2) = - 1
                unique_num2 = unique_num2 - 1
                hi = i
                exit
              end if
            end if
          end if
        end do
!
!  If INDX(I) is still counted as tolerably unique, we now look for
!  nearby temporary points that can be marked as nonunique.
!
        do j = i + 1, hi
          if ( unique(indx(j)) ) then
            if ( n1 < indx(j) ) then
              dist = sqrt ( sum ( ( a2(1:m,indx(i)-n1) &
                                  - a2(1:m,indx(j)-n1) )**2 ) )
              if ( dist <= tol ) then
                unique(indx(j)) = .false.
                xdnu(indx(j)) = xdnu(indx(i))
              end if
            end if
          end if
        end do

      end if

    end if

  end do

  return
end
subroutine point_radial_tol_unique_index_inc1 ( m, n1, a1, tol, seed, z, &
  r1, indx1, unique1, unique_num1, undx1, xdnu1 )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_INDEX_INC1 indexes the tolerably unique points.
!
!  Discussion:
!
!    The input data includes an M x N1 array A1
!    representing the M-dimensional coordinates of a set of N1
!    "permanent" points.
!
!    The output is:
!    * the number of tolerably unique points in the list;
!    * the index, in the list of unique items, of the representatives
!      of each point;
!    * the index, in A, of the tolerably unique representatives.
!
!    In addition, in order to allow for the temporary inclusion of
!    an augmenting set of data, we must output:
!    * a vector Z used to define a sorting scalar;
!    * a vector R1 containing the sorting scalar for each data item;
!    * a vector INDX1 which ascending sorts the data according to R1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N1, the number of permanent points.
!
!    Input, real ( kind = 8 ) A1(M,N1), the permanent points.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) Z(M), a random base vector used to
!    linearly sort the data.
!
!    Output, real ( kind = 8 ) R1(N1), the scalar values assigned to
!    the data for sorting.
!
!    Output, integer ( kind = 4 ) INDX1(N1), the ascending sort index
!    for A1.
!
!    Output, logical UNIQUE1(N1), is TRUE for each unique permanent point.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM1, the number
!    of tolerably unique points with just the permanent points.
!
!    Output, integer ( kind = 4 ) UNDX1(UNIQUE_NUM1),
!    the index in A1 of the tolerably unique points.
!
!    Output, integer ( kind = 4 ) XDNU1(N1), the index in UNDX1
!    of the tolerably unique point that "represents" this point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1

  real ( kind = 8 ) a1(m,n1)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx1(n1)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2_hi
  integer ( kind = 4 ) j2_lo
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
  real ( kind = 8 ) r1(n1)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tol
  logical unique1(n1)
  integer ( kind = 4 ) undx1(n1)
  integer ( kind = 4 ) unique_num1
  real ( kind = 8 ) w(n1)
  integer ( kind = 4 ) xdnu1(n1)
  real ( kind = 8 ) z(m)
!
!  Assign a base point Z randomly in the convex hull of the permanent points.
!
  call r8vec_uniform_01 ( n1, seed, w )

  do i = 1, m
    z(i) = dot_product ( a1(i,1:n1), w(1:n1) ) / sum ( w(1:n1) )
  end do
!
!  Initialize the permanent point data.
!
  do j1 = 1, n1
    r1(j1) = sqrt ( sum ( ( a1(1:m,j1) - z(1:m) )**2 ) )
  end do

  call r8vec_sort_heap_index_a ( n1, r1, indx1 )

  unique_num1 = 0
  unique1(1:n1) = .true.
!
!  STEP 1:
!  Compare PERMANENT POINTS to PERMANENT POINTS.
!
  do j1 = 1, n1

    if ( unique1(indx1(j1)) ) then

      unique_num1 = unique_num1 + 1
      xdnu1(indx1(j1)) = unique_num1
      undx1(unique_num1) = indx1(j1)
!
!  Look for later points which are close to point INDX(I) in terms of R.
!
      hi = j1

      do while ( hi < n1 )
        if ( r1(indx1(j1)) + tol < r1(indx1(hi+1)) ) then
          exit
        end if
        hi = hi + 1
      end do

      do k1 = j1 + 1, hi
        if ( unique1(indx1(k1)) ) then
          dist = sqrt ( sum ( ( a1(1:m,indx1(j1)) - a1(1:m,indx1(k1)) )**2 ) )
          if ( dist <= tol ) then
            unique1(indx1(k1)) = .false.
            xdnu1(indx1(k1)) = xdnu1(indx1(j1))
          end if
        end if
      end do

    end if

  end do

  return
end
subroutine point_radial_tol_unique_index_inc2 ( m, n1, a1, n2, a2, tol, &
  z, r1, indx1, unique1, unique_num1, undx1, xdnu1, r2, indx2, unique2, &
  unique_num2, undx2, xdnu2 )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_INDEX_INC2 indexes unique temporary points.
!
!  Discussion:
!
!    The input data includes an M x N1 array A1 and an M x N2 array A2,
!    representing the M-dimensional coordinates of a set of N1
!    "permanent" points and N2 "temporary" points.
!
!    For notation, we use "A" to describe the M x (N1+N2) array that would be
!    formed by starting with A1 and appending A2.
!
!    The output includes:
!    * the number of tolerably unique points in the list;
!    * the index, in the list of unique items, of the representatives
!      of each point;
!    * the index, in A, of the tolerably unique representatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N1, the number of permanent points.
!
!    Input, real ( kind = 8 ) A1(M,N1), the permanent points.
!
!    Input, integer ( kind = 4 ) N2, the number of temporary points.
!
!    Input, real ( kind = 8 ) A2(M,N2), the temporary points.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input, real ( kind = 8 ) Z(M), a random base vector used to
!    linearly sort the data.
!
!    Input, real ( kind = 8 ) R1(N1), the scalar values assigned to
!    A1 for sorting.
!
!    Input, integer ( kind = 4 ) INDX1(N1), the ascending sort index
!    for A1.
!
!    Input, logical UNIQUE1(N1), is TRUE for each unique permanent point.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM1, the number
!    of tolerably unique permanent points.
!
!    Input, integer ( kind = 4 ) UNDX1(UNIQUE_NUM1),
!    the index in A1 of the tolerably unique permanent points.
!
!    Input, integer ( kind = 4 ) XDNU1(N1), the index in UNDX1
!    of the tolerably unique permanent point that "represents" this point.
!
!    Output, real ( kind = 8 ) R2(N2), the scalar values assigned to
!    A2 for sorting.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the ascending sort index
!    for A2.
!
!    Output, logical UNIQUE2(N2), is TRUE for unique temporary points.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM2, the number
!    of tolerably unique temporary points.
!
!    Output, integer ( kind = 4 ) UNDX2(UNIQUE_NUM2),
!    the index in A2 of the tolerably unique temporary points, incremented
!    by N1.
!
!    Output, integer ( kind = 4 ) XDNU2(N2), the index, in UNDX1
!    or UNDX2, of the tolerably unique point that "represents" this
!    temporary point.  If the value represents an index in UNDX2, this
!    can be inferred by the fact that its value is greater than UNIQUE_NUM1.
!    To reference UNDX2, the value should then be decremented by
!    UNIQUE_NUM1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a1(m,n1)
  real ( kind = 8 ) a2(m,n2)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx1(n1)
  integer ( kind = 4 ) indx2(n2)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2_hi
  integer ( kind = 4 ) j2_lo
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
  real ( kind = 8 ) r1(n1)
  real ( kind = 8 ) r2(n2)
  real ( kind = 8 ) tol
  logical unique1(n1)
  logical unique2(n2)
  integer ( kind = 4 ) undx1(n1)
  integer ( kind = 4 ) undx2(n2)
  integer ( kind = 4 ) unique_num1
  integer ( kind = 4 ) unique_num2
  real ( kind = 8 ) w(n1)
  integer ( kind = 4 ) xdnu1(n1)
  integer ( kind = 4 ) xdnu2(n2)
  real ( kind = 8 ) z(m)
!
!  Initialize the temporary point data.
!
  do j2 = 1, n2
    r2(j2) = sqrt ( sum ( ( a2(1:m,j2) - z(1:m) )**2 ) )
  end do

  call r8vec_sort_heap_index_a ( n2, r2, indx2 )

  unique2(1:n2) = .true.

  unique_num2 = 0
!
!  STEP 2:
!  Use PERMANENT points to eliminate TEMPORARY points.
!
  do j1 = 1, n1

    if ( unique1(indx1(j1)) ) then

      r_lo = r1(indx1(j1)) - tol
      r_hi = r1(indx1(j1)) + tol

      call r8vec_index_sorted_range ( n2, r2, indx2, r_lo, r_hi, &
        j2_lo, j2_hi )

      do j2 = j2_lo, j2_hi
        if ( unique2(indx2(j2)) ) then
          dist = sqrt ( sum ( ( a1(1:m,indx1(j1)) &
                              - a2(1:m,indx2(j2)) )**2 ) )
          if ( dist <= tol ) then
            unique2(indx2(j2)) = .false.
            xdnu2(indx2(j2)) = xdnu1(indx1(j1))
          end if
        end if
      end do

    end if

  end do
!
!  STEP 3:
!  Use TEMPORARY points to eliminate TEMPORARY points.
!
  do j2 = 1, n2

    if ( unique2(indx2(j2)) ) then

      unique_num2 = unique_num2 + 1
      xdnu2(indx2(j2)) = unique_num1 + unique_num2
      undx2(unique_num2) = indx2(j2) + n1

      hi = j2

      do while ( hi < n2 )
        if ( r2(indx2(j2)) + tol < r2(indx2(hi+1)) ) then
          exit
        end if
        hi = hi + 1
      end do

      do k2 = j2 + 1, hi
        if ( unique2(indx2(k2)) ) then
          dist = sqrt ( sum ( ( a2(1:m,indx2(j2)) &
                              - a2(1:m,indx2(k2)) )**2 ) )
          if ( dist <= tol ) then
            unique2(indx2(k2)) = .false.
            xdnu2(indx2(k2)) = xdnu2(indx2(j2))
          end if
        end if
      end do

    end if

  end do

  return
end
subroutine point_radial_tol_unique_index_inc3 ( m, &
  n1, a1, r1, indx1, unique1, unique_num1, undx1, xdnu1, &
  n2, a2, r2, indx2, unique2, unique_num2, undx2, xdnu2, &
  n3, a3, r3, indx3, unique3, unique_num3, undx3, xdnu3 )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_INDEX_INC3 merges index data.
!
!  Discussion:
!
!    This function may be called after *INDEX_INC1 has created index
!    information for the permanent data, and *INDEX_INC2 has created
!    augmenting information for a set of temporary data which now is
!    to be merged with the permanent data.
!
!    The function merges the data and index information to create a
!    new "permanent" data set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N1, the number of permanent points.
!
!    Input, real ( kind = 8 ) A1(M,N1), the permanent points.
!
!    Input, real ( kind = 8 ) R1(N1), the scalar values assigned to
!    the data for sorting.
!
!    Input, integer ( kind = 4 ) INDX1(N1), the ascending sort index
!    for A1.
!
!    Input, logical UNIQUE1(N1), is TRUE for each unique permanent point.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM1, the number
!    of tolerably unique points with just the permanent points.
!
!    Input, integer ( kind = 4 ) UNDX1(UNIQUE_NUM1),
!    the index in A1 of the tolerably unique points.
!
!    Input, integer ( kind = 4 ) XDNU1(N1), the index in UNDX1
!    of the tolerably unique point that "represents" this point.
!
!    Input, integer ( kind = 4 ) N2, the number of temporary points.
!
!    Input, real ( kind = 8 ) A2(M,N2), the temporary points.
!
!    Input, real ( kind = 8 ) R2(N2), the scalar values assigned to
!    the data for sorting.
!
!    Input, integer ( kind = 4 ) INDX2(N2), the ascending sort index
!    for A2.
!
!    Input, logical UNIQUE2(N2), is TRUE for each unique temporary point.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM2, the number
!    of tolerably unique temporary points.
!
!    Input, integer ( kind = 4 ) UNDX2(UNIQUE_NUM2),
!    the index in A2 of the tolerably unique points, incremented by UNIQUE_NUM1.
!
!    Input, integer ( kind = 4 ) XDNU2(N2), the index in UNDX1 or UNDX2
!    of the tolerably unique point that "represents" this point.
!
!    Output, integer ( kind = 4 ) N3, the number of permanent points.
!
!    Output, real ( kind = 8 ) A3(M,N3), the permanent points.
!
!    Output, real ( kind = 8 ) R3(N3), the scalar values assigned to
!    the data for sorting.
!
!    Output, integer ( kind = 4 ) INDX3(N3), the ascending sort index
!    for A3.
!
!    Output, logical UNIQUE3(N3), is TRUE for each unique permanent point.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM3, the number
!    of tolerably unique points.
!
!    Output, integer ( kind = 4 ) UNDX3(UNIQUE_NUM3),
!    the index in A3 of the tolerably unique points.
!
!    Output, integer ( kind = 4 ) XDNU3(N3), the index in UNDX3
!    of the tolerably unique point that "represents" this point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a1(m,n1)
  real ( kind = 8 ) a2(m,n2)
  real ( kind = 8 ) a3(m,n1+n2)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) indx1(n1)
  integer ( kind = 4 ) indx2(n2)
  integer ( kind = 4 ) indx3(n1+n2)
  integer ( kind = 4 ) n3
  real ( kind = 8 ) r1(n1)
  real ( kind = 8 ) r2(n2)
  real ( kind = 8 ) r3(n1+n2)
  real ( kind = 8 ) r8_huge
  logical unique1(n1)
  logical unique2(n2)
  logical unique3(n1+n2)
  integer ( kind = 4 ) undx1(n1)
  integer ( kind = 4 ) undx2(n2)
  integer ( kind = 4 ) undx3(n1+n2)
  integer ( kind = 4 ) unique_num1
  integer ( kind = 4 ) unique_num2
  integer ( kind = 4 ) unique_num3
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  integer ( kind = 4 ) xdnu1(n1)
  integer ( kind = 4 ) xdnu2(n2)
  integer ( kind = 4 ) xdnu3(n1+n2)

  n3 = n1 + n2

  a3(1:m,1:n1)       = a1(1:m,1:n1)
  a3(1:m,n1+1:n1+n2) = a2(1:m,1:n2)

  r3(1:n1)       = r1(1:n1)
  r3(n1+1:n1+n2) = r2(1:n2)
!
!  Interleave the two INDX arrays so that INDX3 presents the entries
!  of A3 in ascending R3 order.
!
  i1 = 1
  i2 = 1

  do i3 = 1, n3

    if ( i1 <= n1 ) then
      v1 = r1(indx1(i1))
    else
      v1 = r8_huge ( )
    end if

    if ( i2 <= n2 ) then
      v2 = r2(indx2(i2))
    else
      v2 = r8_huge ( )
    end if

    if ( v1 <= v2 ) then
      indx3(i3) = indx1(i1)
      i1 = i1 + 1
    else
      indx3(i3) = indx2(i2) + n1
      i2 = i2 + 1
    end if

  end do

  unique_num3 = unique_num1 + unique_num2

  unique3(1:n1) = unique1(1:n1)
  unique3(n1+1:n1+n2) = unique2(1:n2)
!
!  The entries in UNDX2 were already incremented by N2 if they pointed
!  to an entry of A2, so all entries in UNDX2 correctly index A3.
!
  undx3(1:unique_num1)                         = undx1(1:unique_num1)
  undx3(unique_num1+1:unique_num1+unique_num2) = undx2(1:unique_num2)
!
!  Note that the entries of XDNU2 were already incremented by N2
!  so that they correctly index A3, not A2.
!
  xdnu3(1:n1)       = xdnu1(1:n1)
  xdnu3(n1+1:n1+n2) = xdnu2(1:n2)

  return
end
subroutine point_unique_index ( m, n, a, unique_num, undx, xdnu )

!*****************************************************************************80
!
!! POINT_UNIQUE_INDEX indexes unique points.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points to the unique elements of A, in sorted order,
!    and a vector XDNU, which identifies, for each entry of A, the index of
!    the unique sorted element of A.
!
!    This is all done with index vectors, so that the elements of
!    A are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of A, which creates arrays INDX and XDNI.  (If all the entries
!    of A are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of A in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the vector A could be
!    replaced by a compressed vector XU, containing the unique entries
!    of A in sorted order, using the formula
!
!      XU(*) = A(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector A, or
!    any element of it, by index, as follows:
!
!      A(I) = XU(XDNU(I)).
!
!    We could then replace A by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of A, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector A, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I     A  Indx  Xdni       XU  Undx  Xdnu
!    ----+-----+-----+-----+--------+-----+-----+
!      0 | 11.     0     0 |    11.     0     0
!      1 | 22.     2     4 |    22.     1     1
!      2 | 11.     5     1 |    33.     3     0
!      3 | 33.     8     7 |    55.     4     2
!      4 | 55.     1     8 |                  3
!      5 | 11.     6     2 |                  0
!      6 | 22.     7     5 |                  1
!      7 | 22.     3     6 |                  1
!      8 | 11.     4     3 |                  0
!
!    INDX(2) = 3 means that sorted item(2) is A(3).
!    XDNI(2) = 5 means that A(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
!    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = A(I).
!    XU(I)        = A(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values,
!
!    Input, real ( kind = 8 ) A(M,N), the data values.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM, the number of unique values in A.
!    This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(N), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) unique_num

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) undx(unique_num)
  integer ( kind = 4 ) xdnu(n)
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Walk through the implicitly sorted array.
!
  i = 1
  j = 1
  undx(j) = indx(i)
  xdnu(indx(i)) = j

  do i = 2, n
    diff = 0.0D+00
    do k = 1, m
      diff = max ( diff, abs ( a(k,indx(i)) - a(k,undx(j)) ) )
    end do
    if ( 0.0D+00 < diff ) then
      j = j + 1
      undx(j) = indx(i)
    end if
    xdnu(indx(i)) = j
  end do

  return
end
subroutine product_mixed_weight ( dim_num, order_1d, order_nd, rule, alpha, &
  beta, weight_nd )

!*****************************************************************************80
!
!! PRODUCT_MIXED_WEIGHT computes the weights of a mixed product rule.
!
!  Discussion:
!
!    This routine computes the weights for a quadrature rule which is
!    a product of 1D rules of varying order and kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ORDER_1D(DIM_NUM), the order of the 1D rules.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the order of the product rule.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
!    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested.
!    12, "F2_SE", Fejer Type 2 Slow Exponential, Closed Fully Nested.
!    13, "GP_SE", Gauss Patterson Slow Exponential, Closed Fully Nested.
!    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested.
!    15, "F2_ME", Fejer Type 2 Moderate Exponential, Closed Fully Nested.
!    16, "GP_ME", Gauss Patterson Moderate Exponential, Closed Fully Nested.
!    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
!
!    Input, real ( kind = 8 ) ALPHA(DIM_NUM), BETA(DIM_NUM), parameters used for
!    Generalized Gauss Hermite, Generalized Gauss Laguerre, and Gauss
!    Jacobi rules.
!
!    Output, real ( kind = 8 ) WEIGHT_ND(ORDER_ND), the product rule weights.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  real ( kind = 8 ) alpha(dim_num)
  real ( kind = 8 ) beta(dim_num)
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) rule(dim_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight_1d
  real ( kind = 8 ) weight_nd(order_nd)

  weight_nd(1:order_nd) = 1.0D+00

  do dim = 1, dim_num

    allocate ( weight_1d(1:order_1d(dim) ) )

    if ( rule(dim) == 1 ) then
      call clenshaw_curtis_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 2 ) then
      call fejer2_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 3 ) then
      call patterson_lookup_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 4 ) then
      call legendre_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 5 ) then
      call hermite_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 6 ) then
      call gen_hermite_compute_weights ( order_1d(dim), alpha(dim), weight_1d )
    else if ( rule(dim) == 7 ) then
      call laguerre_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 8 ) then
      call gen_laguerre_compute_weights ( order_1d(dim), alpha(dim), weight_1d )
    else if ( rule(dim) == 9 ) then
      call jacobi_compute_weights ( order_1d(dim), alpha(dim), beta(dim), &
        weight_1d )
    else if ( rule(dim) == 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT - Fatal error!'
      write ( *, '(a,i8)' ) '  Do not know how to set weights for rule 10.'
      stop
    else if ( rule(dim) == 11 ) then
      call clenshaw_curtis_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 12 ) then
      call fejer2_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 13 ) then
      call patterson_lookup_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 14 ) then
      call clenshaw_curtis_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 15 ) then
      call fejer2_compute_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 16 ) then
      call patterson_lookup_weights ( order_1d(dim), weight_1d )
    else if ( rule(dim) == 17 ) then
      call ccn_compute_weights ( order_1d(dim), weight_1d )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT - Fatal error!'
      write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
      stop
    end if

    call r8vec_direct_product2 ( dim, order_1d(dim), weight_1d, &
      dim_num, order_nd, weight_nd )

    deallocate ( weight_1d )

  end do

  return
end
function r8_ceiling ( r )

!*****************************************************************************80
!
!! R8_CEILING rounds an R8 "up" (towards +oo) to the next integer.
!
!  Example:
!
!    R     Value
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the value to be rounded up.
!
!    Output, integer ( kind = 4 ) R8_CEILING, the rounded value.
!
  implicit none

  real ( kind = 8 ) r
  integer ( kind = 4 ) r8_ceiling
  integer ( kind = 4 ) value

  value = int ( r )
  if ( real ( value, kind = 8 ) < r ) then
    value = value + 1
  end if

  r8_ceiling = value

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
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON, the R8 round-off unit.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) d_test
  real ( kind = 8 ) r8_epsilon

  d = 1.0D+00
  d_test = 1.0D+00 + d / 2.0D+00

  do while ( 1.0D+00 < d_test )
    d = d / 2.0D+00
    d_test = 1.0D+00 + d / 2.0D+00
  end do

  r8_epsilon = d

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial function.
!
!  Discussion:
!
!    factorial ( N ) = N! = product ( 1 <= I <= N ) I
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
!    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL2, the value of N!!.
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
function r8_floor ( r )

!*****************************************************************************80
!
!! R8_FLOOR rounds an R8 "down" (towards -infinity) to the next integer.
!
!  Example:
!
!    R     Value
!
!   -1.1  -2
!   -1.0  -1
!   -0.9  -1
!    0.0   0
!    5.0   5
!    5.1   5
!    5.9   5
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the value to be rounded down.
!
!    Output, integer ( kind = 4 ) R8_FLOOR, the rounded value.
!
  implicit none

  real ( kind = 8 ) r
  integer ( kind = 4 ) r8_floor
  integer ( kind = 4 ) value

  value = int ( r )
  if ( r < real ( value, kind = 8 ) ) then
    value = value - 1
  end if

  r8_floor = value

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: one = 1.0D+00
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
  real ( kind = 8 ), parameter :: twelve = 12.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
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
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  parity = .false.
  fact = one
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= zero ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= zero ) then

      if ( y1 /= aint ( y1 * half ) * two ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + one

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
      res = one / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < twelve ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < one ) then

      z = y
      y = y + one
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - one

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = zero
    xden = one
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + one
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
        y = y + one
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
      sum = sum + ( y - half ) * log ( y )
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

  if ( fact /= one ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
subroutine r8_hyper_2f1 ( a_input, b_input, c_input, x_input, hf )

!*****************************************************************************80
!
!! R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
!
!  Discussion:
!
!    A minor bug was corrected.  The HW variable, used in several places as
!    the "old" value of a quantity being iteratively improved, was not
!    being initialized.  JVB, 11 February 2008.
!
!    The original version of this program allowed the input arguments to
!    be modified, although they were restored to their input values before exit.
!    This is unacceptable if the input arguments are allowed to be constants.
!    The code has been modified so that the input arguments are never modified.
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program provided that the copyright
!    is acknowledged.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2008
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_INPUT, B_INPUT, C_INPUT, X_INPUT,
!    the arguments of the function.  The user is allowed to pass these
!    values as constants or variables.
!    C_INPUT must not be equal to a nonpositive integer.
!    X_INPUT < 1.
!
!    Output, real ( kind = 8 ) HF, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_input
  real ( kind = 8 ) a0
  real ( kind = 8 ) aa
  real ( kind = 8 ) b
  real ( kind = 8 ) b_input
  real ( kind = 8 ) bb
  real ( kind = 8 ) c
  real ( kind = 8 ) c_input
  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ), parameter :: el = 0.5772156649015329D+00
  real ( kind = 8 ) eps
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) g0
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) g3
  real ( kind = 8 ) ga
  real ( kind = 8 ) gabc
  real ( kind = 8 ) gam
  real ( kind = 8 ) gb
  real ( kind = 8 ) gbm
  real ( kind = 8 ) gc
  real ( kind = 8 ) gca
  real ( kind = 8 ) gcab
  real ( kind = 8 ) gcb
  real ( kind = 8 ) gm
  real ( kind = 8 ) hf
  real ( kind = 8 ) hw
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nm
  real ( kind = 8 ) pa
  real ( kind = 8 ) pb
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r0
  real ( kind = 8 ) r1
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) rm
  real ( kind = 8 ) rp
  real ( kind = 8 ) sm
  real ( kind = 8 ) sp
  real ( kind = 8 ) sp0
  real ( kind = 8 ) x
  real ( kind = 8 ) x_input
  real ( kind = 8 ) x1
!
!  Immediately copy the input arguments!
!
  a = a_input
  b = b_input
  c = c_input
  x = x_input

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
  l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )

  if ( l0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  Integral C < 0.'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    stop
  end if

  if ( l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    write ( *, '(a)' ) '  1 - X < 0, C - A - B < 0.'
    stop
  end if

  if ( 0.95D+00 < x ) then
    eps = 1.0D-08
  else
    eps = 1.0D-15
  end if

  if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

    hf = 1.0D+00
    return

  else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then

    gc = r8_gamma ( c )
    gcab = r8_gamma ( c - a - b )
    gca = r8_gamma ( c - a )
    gcb = r8_gamma ( c - b )
    hf = gc * gcab / ( gca * gcb )
    return

  else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then

    g0 = sqrt ( pi ) * 2.0D+00**( - a )
    g1 = r8_gamma ( c )
    g2 = r8_gamma ( 1.0D+00 + a / 2.0D+00 - b )
    g3 = r8_gamma ( 0.5D+00 + 0.5D+00 * a )
    hf = g0 * g1 / ( g2 * g3 )
    return

  else if ( l2 .or. l3 ) then

    if ( l2 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l3 ) then
      nm = int ( abs ( b ) )
    end if

    hf = 1.0D+00
    r = 1.0D+00

    do k = 1, nm
      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do

    return

  else if ( l4 .or. l5 ) then

    if ( l4 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l5 ) then
      nm = int ( abs ( c - b ) )
    end if

    hf = 1.0D+00
    r  = 1.0D+00
    do k = 1, nm
      r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do
    hf = ( 1.0D+00 - x )**( c - a - b ) * hf
    return

  end if

  aa = a
  bb = b
  x1 = x

  if ( x < 0.0D+00 ) then
    x = x / ( x - 1.0D+00 )
    if ( a < c .and. b < a .and. 0.0D+00 < b ) then
      a = bb
      b = aa
    end if
    b = c - b
  end if

  if ( 0.75D+00 <= x ) then

    gm = 0.0D+00

    if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then

      m = int ( c - a - b )
      ga = r8_gamma ( a )
      gb = r8_gamma ( b )
      gc = r8_gamma ( c )
      gam = r8_gamma ( a + m )
      gbm = r8_gamma ( b + m )

      pa = r8_psi ( a )
      pb = r8_psi ( b )

      if ( m /= 0 ) then
        gm = 1.0D+00
      end if

      do j = 1, abs ( m ) - 1
        gm = gm * j
      end do

      rm = 1.0D+00
      do j = 1, abs ( m )
        rm = rm * j
      end do

      f0 = 1.0D+00
      r0 = 1.0D+00
      r1 = 1.0D+00
      sp0 = 0.0D+00
      sp = 0.0D+00

      if ( 0 <= m ) then

        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

        do k = 1, m - 1
          r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
            + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = 8 )
        end do

        f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + ( 1.0D+00 - a ) &
              / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
              + 1.0D+00 / ( b + j + k - 1.0D+00 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      else if ( m < 0 ) then

        m = - m
        c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

        do k = 1, m - 1
          r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / real ( k, kind = 8 )
        end do

        f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) &
            / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + 1.0D+00 / real ( j + k, kind = 8 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      end if

    else

      ga = r8_gamma ( a )
      gb = r8_gamma ( b )
      gc = r8_gamma ( c )
      gca = r8_gamma ( c - a )
      gcb = r8_gamma ( c - b )
      gcab = r8_gamma ( c - a - b )
      gabc = r8_gamma ( a + b - c )
      c0 = gc * gcab / ( gca * gcb )
      c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
      hf = 0.0D+00
      hw = hf
      r0 = c0
      r1 = c1

      do k = 1, 250

        r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
          / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

        r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
          / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

        hf = hf + r0 + r1

        if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
          exit
        end if

        hw = hf

      end do

      hf = hf + c0 + c1

    end if

  else

    a0 = 1.0D+00

    if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then

      a0 = ( 1.0D+00 - x )**( c - a - b )
      a = c - a
      b = c - b

    end if

    hf = 1.0D+00
    hw = hf
    r = 1.0D+00

    do k = 1, 250

      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x

      hf = hf + r

      if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
        exit
      end if

      hw = hf

    end do

    hf = a0 * hf

  end if

  if ( x1 < 0.0D+00 ) then
    x = x1
    c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
    hf = c0 * hf
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Warning!'
    write ( *, '(a)' ) '  A large number of iterations were needed.'
    write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8 value.
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
function r8_psi ( xx )

!*****************************************************************************80
!
!! R8_PSI evaluates the function Psi(X).
!
!  Discussion:
!
!    This routine evaluates the logarithmic derivative of the
!    Gamma function,
!
!      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
!             = d/dX LN ( GAMMA(X) )
!
!    for real X, where either
!
!      - XMAX1 < X < - XMIN, and X is not a negative integer,
!
!    or
!
!      XMIN < X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, Number 121, January 1973, pages 123-127.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XX, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_PSI, the value of the function.
!
  implicit none

  real ( kind = 8 ) aug
  real ( kind = 8 ) den
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: fourth = 0.25D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nq
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), dimension ( 9 ) :: p1 = (/ &
   4.5104681245762934160D-03, &
   5.4932855833000385356D+00, &
   3.7646693175929276856D+02, &
   7.9525490849151998065D+03, &
   7.1451595818951933210D+04, &
   3.0655976301987365674D+05, &
   6.3606997788964458797D+05, &
   5.8041312783537569993D+05, &
   1.6585695029761022321D+05 /)
  real ( kind = 8 ), dimension ( 7 ) :: p2 = (/ &
  -2.7103228277757834192D+00, &
  -1.5166271776896121383D+01, &
  -1.9784554148719218667D+01, &
  -8.8100958828312219821D+00, &
  -1.4479614616899842986D+00, &
  -7.3689600332394549911D-02, &
  -6.5135387732718171306D-21 /)
  real ( kind = 8 ), parameter :: piov4 = 0.78539816339744830962D+00
  real ( kind = 8 ), dimension ( 8 ) :: q1 = (/ &
   9.6141654774222358525D+01, &
   2.6287715790581193330D+03, &
   2.9862497022250277920D+04, &
   1.6206566091533671639D+05, &
   4.3487880712768329037D+05, &
   5.4256384537269993733D+05, &
   2.4242185002017985252D+05, &
   6.4155223783576225996D-08 /)
  real ( kind = 8 ), dimension ( 6 ) :: q2 = (/ &
   4.4992760373789365846D+01, &
   2.0240955312679931159D+02, &
   2.4736979003315290057D+02, &
   1.0742543875702278326D+02, &
   1.7463965060678569906D+01, &
   8.8427520398873480342D-01 /)
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) sgn
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ) upper
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: x01 = 187.0D+00
  real ( kind = 8 ), parameter :: x01d = 128.0D+00
  real ( kind = 8 ), parameter :: x02 = 6.9464496836234126266D-04
  real ( kind = 8 ), parameter :: xinf = 1.70D+38
  real ( kind = 8 ), parameter :: xlarge = 2.04D+15
  real ( kind = 8 ), parameter :: xmax1 = 3.60D+16
  real ( kind = 8 ), parameter :: xmin1 = 5.89D-39
  real ( kind = 8 ), parameter :: xsmall = 2.05D-09
  real ( kind = 8 ) xx
  real ( kind = 8 ) z
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  x = xx
  w = abs ( x )
  aug = zero
!
!  Check for valid arguments, then branch to appropriate algorithm.
!
  if ( xmax1 <= - x .or. w < xmin1 ) then

    if ( zero < x ) then
      r8_psi = - xinf
    else
      r8_psi = xinf
    end if

    return
  end if

  if ( x < half ) then
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
    if ( w <= xsmall ) then

      aug = - one / x
!
!  Argument reduction for cotangent.
!
    else

      if ( x < zero ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - real ( int ( w ), kind = 8 )
      nq = int ( w * four )
      w = four * ( w - real ( nq, kind = 8 ) * fourth )
!
!  W is now related to the fractional part of 4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
      n = nq / 2

      if ( n + n /= nq ) then
        w = one - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) /= 0 ) then
        sgn = - sgn
      end if
!
!  Determine the final value for  -pi * cotan(pi*x).
!
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) == 0 ) then
!
!  Check for singularity.
!
        if ( z == zero ) then

          if ( zero < x ) then
            r8_psi = -xinf
          else
            r8_psi = xinf
          end if

          return
        end if

        aug = sgn * ( four / tan ( z ) )

      else

        aug = sgn * ( four * tan ( z ) )

      end if

    end if

    x = one - x

  end if
!
!  0.5 <= X <= 3.0.
!
  if ( x <= three ) then

    den = x
    upper = p1(1) * x
    do i = 1, 7
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do
    den = ( upper + p1(9) ) / ( den + q1(8) )
    x = ( x - x01 / x01d ) - x02
    r8_psi = den * x + aug
    return

  end if
!
!  3.0 < X.
!
  if ( x < xlarge ) then
    w = one / ( x * x )
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do
    aug = ( upper + p2(7) ) / ( den + q2(6) ) - half / x + aug
  end if

  r8_psi = aug + log ( x )

  return
end
subroutine r8col_compare ( m, n, a, i, j, value )

!*****************************************************************************80
!
!! R8COL_COMPARE compares columns in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      VALUE = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) VALUE, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) value
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  end if

  value = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      value = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      value = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8col_sort_heap_a ( m, n, a )

!*****************************************************************************80
!
!! R8COL_SORT_HEAP_A ascending heapsorts an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r8col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call r8col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8col_sort_heap_index_a ( m, n, a, indx )

!*****************************************************************************80
!
!! R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2)
!    is negative.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(*,INDX(*)) is sorted,
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in each column of A.
!
!    Input, integer ( kind = 4 )  N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The I-th element
!    of the sorted array is column INDX(I).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) column(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = ( n / 2 ) + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      column(1:m) = a(1:m,indxt)

    else

      indxt = indx(ir)
      column(1:m) = a(1:m,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then

        call r8vec_compare ( m, a(1:m,indx(j)), a(1:m,indx(j+1)), isgn )

        if ( isgn < 0 ) then
          j = j + 1
        end if

      end if

      call r8vec_compare ( m, column, a(1:m,indx(j)), isgn )

      if ( isgn < 0 ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8col_sorted_tol_undex ( m, n, a, unique_num, tol, undx, xdnu )

!*****************************************************************************80
!
!! R8COL_SORTED_TOL_UNDEX indexes tolerably unique entries in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the tolerably unique elements of A, in sorted order,
!    and a vector XDNU, which identifies, for each entry of A, the index of
!    the unique sorted element of A.
!
!    This is all done with index vectors, so that the elements of
!    A are never moved.
!
!    Assuming A is already sorted, we examine the entries of A in order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the vector A could be
!    replaced by a compressed vector XU, containing the unique entries
!    of A in sorted order, using the formula
!
!      XU(*) = A(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector A, or
!    any element of it, by index, as follows:
!
!      A(I) = XU(XDNU(I)).
!
!    We could then replace A by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of A, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector A, the unique sort and
!    inverse unique sort vectors and the compressed unique sorted vector.
!
!      I      A      XU  Undx  Xdnu
!    ----+------+------+-----+-----+
!      1 | 11.0 |  11.0    1     1
!      2 | 11.0 |  22.0    5     1
!      3 | 11.0 |  33.0    8     1
!      4 | 11.0 |  55.0    9     1
!      5 | 22.0 |                2
!      6 | 22.0 |                2
!      7 | 22.0 |                2
!      8 | 33.0 |                3
!      9 | 55.0 |                4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(M,N), the data values.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM, the number of unique values
!    in A.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(N), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) unique_num

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(unique_num)
  logical unique
  integer ( kind = 4 ) xdnu(n)
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!  Set the representative of item I to the K-th unique item.
!
  i = 1
  k = 1
  undx(k) = i
  xdnu(i) = k
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K, set the
!  K-th unique item to I, and set the representative of I to K.
!
!  If it is not unique, set the representative of item I to a
!  previously determined unique item that is close to it.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      i2 = undx(j)
      diff = maxval ( abs ( a(1:m,i) - a(1:m,i2) ) )
      if ( diff <= tol ) then
        unique = .false.
        xdnu(i) = j
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = i
      xdnu(i) = k
    end if

  end do

  return
end
subroutine r8col_sorted_tol_unique_count ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_TOL_UNIQUE_COUNT: tolerably unique elements in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!    If the tolerance is large enough, then the concept of uniqueness
!    can become ambiguous.  If we have a tolerance of 1.5, then in the
!    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
!    one unique entry?  That would be because 1 may be regarded as unique,
!    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
!    be unique and so on.
!
!    This seems wrongheaded.  So I prefer the idea that an item is not
!    unique under a tolerance only if it is close to something that IS unique.
!    Thus, the unique items are guaranteed to cover the space if we include
!    a disk of radius TOL around each one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(n)
  logical unique
  integer unique_num
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!
  i = 1
  k = 1
  undx(k) = i
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K and set the
!  K-th unique item to I.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      i2 = undx(j)
      diff = maxval ( abs ( a(1:m,i) - a(1:m,i2) ) )
      if ( diff <= tol ) then
        unique = .false.
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = i
    end if

  end do

  unique_num = k

  return
end
subroutine r8col_sorted_unique_count ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!    If the tolerance is large enough, then the concept of uniqueness
!    can become ambiguous.  If we have a tolerance of 1.5, then in the
!    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
!    one unique entry?  That would be because 1 may be regarded as unique,
!    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
!    be unique and so on.
!
!    This seems wrongheaded.  So I prefer the idea that an item is not
!    unique under a tolerance only if it is close to something that IS unique.
!    Thus, the unique items are guaranteed to cover the space if we include
!    a disk of radius TOL around each one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( tol < maxval ( abs ( a(1:m,j1) - a(1:m,j2) ) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine r8col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! R8COL_SWAP swaps columns I and J of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      A = (
!        1.  4.  3.  2.
!        5.  8.  7.  6.
!        9. 12. 11. 10. )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) col(m)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i8)' ) '  J1 =    ', j1
    write ( *, '(a,i8)' ) '  J2 =    ', j2
    write ( *, '(a,i8)' ) '  NCOL = ', n
    stop
  end if

  if ( j1 == j2 ) then
    return
  end if

  col(1:m) = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)

  return
end
subroutine r8col_tol_undex ( m, n, a, unique_num, tol, undx, xdnu )

!*****************************************************************************80
!
!! R8COL_TOL_UNDEX indexes tolerably unique entries of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of A, in sorted order,
!    and a vector XDNU, which identifies, for each entry of A, the index of
!    the unique sorted element of A.
!
!    This is all done with index vectors, so that the elements of
!    A are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of A, which creates arrays INDX and XDNI.  (If all the entries
!    of A are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of A in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the object X could be
!    replaced by a compressed object XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(*) = A(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector A, or
!    any element of it, by index, as follows:
!
!      A(I) = XU(XDNU(I)).
!
!    We could then replace A by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of A, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector A, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I    A   Indx  Xdni      XU   Undx  Xdnu
!    ----+-----+-----+-----+--------+-----+-----+
!      1 | 11.     1     1 |    11.     1     1
!      2 | 22.     3     5 |    22.     2     2
!      3 | 11.     6     2 |    33.     4     1
!      4 | 33.     9     8 |    55.     5     3
!      5 | 55.     2     9 |                  4
!      6 | 11.     7     3 |                  1
!      7 | 22.     8     6 |                  2
!      8 | 22.     4     7 |                  2
!      9 | 11.     5     4 |                  1
!
!    INDX(2) = 3 means that sorted item(2) is A(3).
!    XDNI(2) = 5 means that A(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
!    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = A(I).
!    XU(I)        = A(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(M,N), the data values.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM, the number of unique values
!    in A.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(N), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) unique_num

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(unique_num)
  logical unique
  integer ( kind = 4 ) xdnu(n)
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!  Set the representative of item I to the K-th unique item.
!
  i = 1
  k = 1
  undx(k) = indx(i)
  xdnu(indx(i)) = k
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K, set the
!  K-th unique item to I, and set the representative of I to K.
!
!  If it is not unique, set the representative of item I to a
!  previously determined unique item that is close to it.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      diff = maxval ( abs ( a(1:m,indx(i)) - a(1:m,undx(j)) ) )
      if ( diff <= tol ) then
        unique = .false.
        xdnu(indx(i)) = j
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = indx(i)
      xdnu(indx(i)) = k
    end if

  end do

  return
end
subroutine r8col_tol_unique_count ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8COL_TOL_UNIQUE_COUNT counts tolerably unique entries in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    If the tolerance is large enough, then the concept of uniqueness
!    can become ambiguous.  If we have a tolerance of 1.5, then in the
!    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
!    one unique entry?  That would be because 1 may be regarded as unique,
!    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
!    be unique and so on.
!
!    This seems wrongheaded.  So I prefer the idea that an item is not
!    unique under a tolerance only if it is close to something that IS unique.
!    Thus, the unique items are guaranteed to cover the space if we include
!    a disk of radius TOL around each one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(n)
  logical unique
  integer ( kind = 4 ) unique_num
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!  Set the representative of item I to the K-th unique item.
!
  i = 1
  k = 1
  undx(k) = indx(i)
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K, set the
!  K-th unique item to I, and set the representative of I to K.
!
!  If it is not unique, set the representative of item I to a
!  previously determined unique item that is close to it.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      diff = maxval ( abs ( a(1:m,indx(i)) - a(1:m,undx(j)) ) )
      if ( diff <= tol ) then
        unique = .false.
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = indx(i)
    end if

  end do

  unique_num = k

  return
end
subroutine r8col_undex ( x_dim, x_num, x_val, x_unique_num, tol, undx, xdnu )

!*****************************************************************************80
!
!! R8COL_UNDEX returns unique sorted indexes for an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of X, in sorted order,
!    and a vector XDNU, which identifies, for each entry of X, the index of
!    the unique sorted element of X.
!
!    This is all done with index vectors, so that the elements of
!    X are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of X, which creates arrays INDX and XDNI.  (If all the entries
!    of X are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of X in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the object X could be
!    replaced by a compressed object XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(*) = X(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector X, or
!    any element of it, by index, as follows:
!
!      X(I) = XU(XDNU(I)).
!
!    We could then replace X by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of X, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector X, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I    X   Indx  Xdni      XU   Undx  Xdnu
!    ----+-----+-----+-----+--------+-----+-----+
!      1 | 11.     1     1 |    11.     1     1
!      2 | 22.     3     5 |    22.     2     2
!      3 | 11.     6     2 |    33.     4     1
!      4 | 33.     9     8 |    55.     5     3
!      5 | 55.     2     9 |                  4
!      6 | 11.     7     3 |                  1
!      7 | 22.     8     6 |                  2
!      8 | 22.     4     7 |                  2
!      9 | 11.     5     4 |                  1
!
!    INDX(2) = 3 means that sorted item(2) is X(3).
!    XDNI(2) = 5 means that X(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
!    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = X(I).
!    XU(I)        = X(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_DIM, the dimension of the data values.
!    (the number of rows in the R8COL).
!
!    Input, integer ( kind = 4 ) X_NUM, the number of data values.
!    (the number of columns in the R8COL).
!
!    Input, real ( kind = 8 ) X_VAL(X_DIM,X_NUM), the data values.
!
!    Input, integer ( kind = 4 ) X_UNIQUE_NUM, the number of unique values
!    in X_VAL.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNDX(X_UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(X_NUM), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) x_dim
  integer ( kind = 4 ) x_num
  integer ( kind = 4 ) x_unique_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(x_num)
  integer ( kind = 4 ) j
  real ( kind = 8 ) tol
  integer ( kind = 4 ) undx(x_unique_num)
  real ( kind = 8 ) x_val(x_dim,x_num)
  integer ( kind = 4 ) xdnu(x_num)
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( x_dim, x_num, x_val, indx )
!
!  Walk through the implicitly sorted array X.
!
  i = 1

  j = 1
  undx(j) = indx(i)

  xdnu(indx(i)) = j

  do i = 2, x_num

    if ( tol < &
         maxval ( abs ( x_val(1:x_dim,indx(i)) - x_val(1:x_dim,undx(j)) ) ) &
    ) then
      j = j + 1
      undx(j) = indx(i)
    end if

    xdnu(indx(i)) = j

  end do

  return
end
subroutine r8col_unique_index ( m, n, a, tol, unique_index )

!*****************************************************************************80
!
!! R8COL_UNIQUE_INDEX indexes the unique occurrence of values in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness
!    index of A(1:M,J).  That is, if A_UNIQUE contains the unique elements
!    of A, gathered in order, then
!
!      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!    The length of an "element" of A, and the number of "elements".
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_INDEX(N), the first occurrence index.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_index(n)
  integer ( kind = 4 ) unique_num

  unique_index(1:n) = -1
  unique_num = 0

  do j1 = 1, n

    if ( unique_index(j1) == -1 ) then

      unique_num = unique_num + 1
      unique_index(j1) = unique_num

      do j2 = j1 + 1, n
        if ( maxval ( abs ( a(1:m,j1) - a(1:m,j2) ) ) <= tol ) then
          unique_index(j2) = unique_num
        end if
      end do

    end if

  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
!    14 June 2004
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
!    14 June 2004
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine r8poly_ant_val ( n, c, xv, yv )

!*****************************************************************************80
!
!! R8POLY_ANT_VAL evaluates the antiderivative of a polynomial in standard form.
!
!  Discussion:
!
!    The constant term of the antiderivative is taken to be zero.
!
!    Thus, if 
!      p(x) = c(1)     + c(2) * x   + c(3) * x^2   + ... + c(n) * x^(n-1)
!    then q(x), the antiderivative, is taken to be:
!      q(x) = c(1) * x + c(2) * x/2 + c(3) * x^3/3 + ... + c(n) * x^n/n
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) C(N), the polynomial coefficients.
!    C(1) is the constant term, and C(N) is the coefficient of X^(N-1).
!
!    Input, real ( kind = 8 ) XV, the evaluation point.
!
!    Output, real ( kind = 8 ) YV, the value of the antiderivative of 
!    the polynomial at XV.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) xv
  real ( kind = 8 ) yv

  yv = 0.0D+00

  do i = n, 1, -1
    yv = ( yv + c(i) / real ( i, kind = 8 ) ) * xv
  end do

  return
end
subroutine r8vec_chebyshev ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_CHEBYSHEV creates a vector of Chebyshev spaced values.
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
!    08 June 2011
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
!    Output, real ( kind = 8 ) A(N), a vector of Chebyshev spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n

      theta = real ( n - i, kind = 8 ) * pi &
            / real ( n - 1, kind = 8 )

      c = cos ( theta )

      if ( mod ( n, 2 ) == 1 ) then
        if ( 2 * i - 1 == n ) then
          c = 0.0D+00
        end if
      end if

      a(i) = ( ( 1.0D+00 - c ) * a_first  &
             + ( 1.0D+00 + c ) * a_last ) &
             /   2.0D+00

    end do

  end if

  return
end
subroutine r8vec_compare ( n, a1, a2, isgn )

!*****************************************************************************80
!
!! R8VEC_COMPARE compares two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2.0, 6.0, 2.0 )
!      A2 = ( 2.0, 8.0, 12.0 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A1 > A2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) k

  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = -1
      return
    else if ( a2(k) < a1(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of values
!    to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer SKIP, the distance from the current value of START
!    to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine r8vec_index_sorted_range ( n, r, indx, r_lo, r_hi, i_lo, i_hi )

!*****************************************************************************80
!
!! R8VEC_INDEX_SORTED_RANGE: search index sorted vector for elements in a range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!
!    Input, real ( kind = 8 ) R(N), the index sorted vector.
!
!    Input, integer ( kind = 4 ) INDX(N), the vector used to sort R.
!    The vector R(INDX(*)) is sorted.
!
!    Input, real ( kind = 8 ) R_LO, R_HI, the limits of the range.
!
!    Output, integer ( kind = 4 ) I_LO, I_HI, the range of indices
!    so that I_LO <= I <= I_HI => R_LO <= R(INDX(I)) <= R_HI.  If no
!    values in R lie in the range, then I_HI < I_LO will be returned.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
!
!  Cases we can handle immediately.
!
  if ( r(indx(n)) < r_lo ) then
    i_lo = n + 1
    i_hi = n
    return
  end if

  if ( r_hi < r(indx(1)) ) then
    i_lo = 1
    i_hi = 0
    return
  end if
!
!  Are there are least two intervals?
!
  if ( n == 1 ) then
    if ( r_lo <= r(indx(1)) .and. r(indx(1)) <= r_hi ) then
      i_lo = 1
      i_hi = 1
    else
      i_lo = 0
      i_hi = -1
    end if
    return
  end if
!
!  Bracket R_LO.
!
  if ( r_lo <= r(indx(1)) ) then

    i_lo = 1

  else
!
!  R_LO is in one of the intervals spanned by R(INDX(J1)) to R(INDX(J2)).
!  Examine the intermediate interval [R(INDX(I1)), R(INDX(I1+1))].
!  Does R_LO lie here, or below or above?
!
    j1 = 1
    j2 = n
    i1 = ( j1 + j2 - 1 ) / 2
    i2 = i1 + 1

    do

      if ( r_lo < r(indx(i1)) ) then
        j2 = i1
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else if ( r(indx(i2)) < r_lo ) then
        j1 = i2
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else
        i_lo = i1
        exit
      end if

    end do

  end if
!
!  Bracket R_HI.
!
  if ( r(indx(n)) <= r_hi ) then

    i_hi = n

  else

    j1 = i_lo
    j2 = n
    i1 = ( j1 + j2 - 1 ) / 2
    i2 = i1 + 1

    do

      if ( r_hi < r(indx(i1)) ) then
        j2 = i1
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else if ( r(indx(i2)) < r_hi ) then
        j1 = i2
        i1 = ( j1 + j2 - 1 ) / 2
        i2 = i1 + 1
      else
        i_hi = i2
        exit
      end if

    end do

  end if
!
!  We expect to have computed the largest I_LO and smallest I_HI such that
!    R(INDX(I_LO)) <= R_LO <= R_HI <= R(INDX(I_HI))
!  but what we want is actually
!    R_LO <= R(INDX(I_LO)) <= R(INDX(I_HI)) <= R_HI
!  which we can usually get simply by incrementing I_LO and decrementing I_HI.
!
  if ( r(indx(i_lo)) < r_lo ) then
    i_lo = i_lo + 1
    if ( n < i_lo ) then
      i_hi = i_lo - 1
    end if
  end if

  if ( r_hi < r(indx(i_hi)) ) then
    i_hi = i_hi - 1
    if ( i_hi < 1 ) then
      i_lo = i_hi + 1
    end if
  end if

  return
end
subroutine r8vec_legendre ( n, x_first, x_last, x )

!*****************************************************************************80
!
!! R8VEC_LEGENDRE creates a vector of Legendre-spaced values.
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
!    17 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X_FIRST, X_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of Legendre-spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_first
  real ( kind = 8 ) x_last

  call legendre_zeros ( n, x )

  x(1:n) = ( ( 1.0D+00 - x(1:n) ) * x_first  &
           + ( 1.0D+00 + x(1:n) ) * x_last ) &
           /   2.0D+00

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
function r8vec_min_pos ( n, a )

!*****************************************************************************80
!
!! R8VEC_MIN_POS returns the minimum positive value of an R8VEC.
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
!    08 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) R8VEC_MIN_POS, the smallest positive entry,
!    or R8_HUGE if no entry is positive.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_huge = 1.0D+30
  real ( kind = 8 ) r8vec_min_pos
  real ( kind = 8 ) value

  value = r8_huge

  do i = 1, n
    if ( 0.0D+00 < a(i) ) then
      value = min ( value, a(i) )
    end if
  end do

  r8vec_min_pos = value

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

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, indx, a )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

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
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

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
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
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
subroutine vec_colex_next3 ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_COLEX_NEXT3 generates vectors in colex order.
!
!  Discussion:
!
!    The vectors are produced in colexical order, starting with
!
!    (1,        1,        ...,1),
!    (2,        1,        ...,1),
!     ...
!    (BASE(1),  1,        ...,1)
!
!    (1,        2,        ...,1)
!    (2,        2,        ...,1)
!    ...
!    (BASE(1),  2,        ...,0)
!
!    (1,        3,        ...,1)
!    (2,        3,        ...,1)
!    ...
!    (BASE(1),  BASE(2),  ...,BASE(DIM_NUM)).
!
!  Example:
!
!    DIM_NUM = 2,
!    BASE = ( 3, 3 )
!
!    1   1
!    2   1
!    3   1
!    1   2
!    2   2
!    3   2
!    1   3
!    2   3
!    3   3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases to be used in each
!    dimension.  In dimension I, entries will range from 1 to BASE(I).
!
!    Input/output, integer ( kind = 4 ) A(DIM_NUM).
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output
!    vector and stop calling the routine.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 1
    more = .true.

  else

    do i = 1, dim_num

      a(i) = a(i) + 1

      if ( a(i) <= base(i) ) then
        return
      end if

      a(i) = 1

    end do

    more = .false.

  end if

  return
end
