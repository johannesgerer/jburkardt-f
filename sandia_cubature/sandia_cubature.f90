subroutine c1_geg_monomial_integral ( alpha, expon, value )

!*****************************************************************************80
!
!! C1_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on C1.
!
!  Discussion:
!
!    C1_GEG is the interval [-1,+1] with the Gegenbauer weight function
!
!      w(alpha;x) = (1-x^2)^alpha
!
!    with -1.0 < alpha.
!
!    value = integral ( -1 <= x <= +1 ) x^expon (1-x^2)^alpha dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X^2).
!    - 1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
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

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C1_GEG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( expon < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C1_GEG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  EXPON < 0.'
    stop
  end if

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
subroutine c1_jac_monomial_integral ( alpha, beta, expon, value )

!*****************************************************************************80
!
!! C1_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over C1.
!
!  Discussion:
!
!    C1_JAC is the interval [-1,+1] with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = (1-x)^beta (1+x)^alpha.
!
!    with -1.0 < alpha, -1.0 < beta.
!
!    value = integral ( -1 <= x <= +1 ) x^expon (1-x)^alpha (1+x)^beta dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X) in the weight factor.
!    - 1.0 < ALPHA.
!
!    Input, real ( kind = 8 ) BETA, the exponent of (1+X) in the weight factor.
!    - 1.0 < BETA.
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
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

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C1_JAC_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C1_JAC_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  if ( expon < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C1_JAC_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  EXPON < 0.'
    stop
  end if

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
subroutine c1_leg_monomial_integral ( expon, value )

!*****************************************************************************80
!
!! C1_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on C1.
!
!  Discussion:
!
!    C1_LEG is the interval [-1,+1] with the Legendre weight function
!
!      w(x) = 1.
!
!    value = integral ( -1 <= x <= +1 ) x^expon dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2010
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
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) expon
  real ( kind = 8 ) value

  if ( expon < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C1_LEG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  EXPON < 0.'
    stop
  end if

  if ( mod ( expon, 2 ) == 1 ) then
    value = 0.0D+00
    return
  end if

  value = 2.0D+00 / real ( expon + 1, kind = 8 )

  return
end
subroutine cn_geg_00_1 ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! CN_GEG_00_1 implements the midpoint rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call c1_geg_monomial_integral ( alpha, expon, volume_1d )
  volume = volume_1d ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = volume

  return
end
subroutine cn_geg_00_1_size ( n, alpha, o )

!*****************************************************************************80
!
!! CN_GEG_00_1_SIZE sizes the midpoint rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_00_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine cn_geg_01_1 ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! CN_GEG_01_1 implements a precision 1 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_01_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call c1_geg_monomial_integral ( alpha, expon, volume_1d )
  volume = volume_1d ** n

  expon = 1
  call c1_geg_monomial_integral ( alpha, expon, value2 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = value2 / volume_1d
  w(k) = volume

  return
end
subroutine cn_geg_01_1_size ( n, alpha, o )

!*****************************************************************************80
!
!! CN_GEG_01_1_SIZE sizes a precision 1 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_01_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine cn_geg_02_xiu ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! CN_GEG_02_XIU implements the Xiu precision 2 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_02_XIU - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = 1.0D+00
  delta0 = 0.0D+00
  c1 = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call c1_geg_monomial_integral ( alpha, expon, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_geg_02_xiu_size ( n, alpha, o )

!*****************************************************************************80
!
!! CN_GEG_02_XIU_SIZE sizes the Xiu precision 2 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_02_XIU_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = n + 1

  return
end
subroutine cn_geg_03_xiu ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! CN_GEG_03_XIU implements the Xiu precision 3 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_03_XIU - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call c1_geg_monomial_integral ( alpha, expon, volume_1d )
  volume = volume_1d ** n

  do j = 1, o

    i = 0 
    do r = 1, n / 2
      arg = real ( ( 2 * r - 1 ) * j, kind = 8 ) * pi / real ( n, kind = 8 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) &
        / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) &
        / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * r8_mop ( j ) &
        / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
      if ( n == 1 ) then
        x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_geg_03_xiu_size ( n, alpha, o )

!*****************************************************************************80
!
!! CN_GEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_03_XIU_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 2 * n

  return
end
subroutine cn_geg_monomial_integral ( n, alpha, expon, value )

!*****************************************************************************80
!
!! CN_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on CN.
!
!  Discussion:
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!    value = integral ( CN ) 
!      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i)^2)^alpha dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X).
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  value = 1.0D+00
  do i = 1, n
    call c1_geg_monomial_integral ( alpha, expon(i), value2 )
    value = value * value2
  end do

  return
end
subroutine cn_jac_00_1 ( n, alpha, beta, o, x, w )

!*****************************************************************************80
!
!! CN_JAC_00_1 implements the midpoint rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1.0 < alpha, -1.0 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_00_1 - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  expon = 0
  call c1_jac_monomial_integral ( alpha, beta, expon, volume_1d )
  volume = volume_1d ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = volume

  return
end
subroutine cn_jac_00_1_size ( n, alpha, beta, o )

!*****************************************************************************80
!
!! CN_JAC_00_1_SIZE sizes the midpoint rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1.0 < alpha, -1.0 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_00_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_00_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine cn_jac_01_1 ( n, alpha, beta, o, x, w )

!*****************************************************************************80
!
!! CN_JAC_01_1 implements a precision 1 rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha. 
!
!    with -1.0 < alpha, -1.0 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_01_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_01_1 - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  expon = 0
  call c1_jac_monomial_integral ( alpha, beta, expon, volume_1d )
  volume = volume_1d ** n

  expon = 1
  call c1_jac_monomial_integral ( alpha, beta, expon, value2 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = value2 / volume_1d
  w(k) = volume

  return
end
subroutine cn_jac_01_1_size ( n, alpha, beta, o )

!*****************************************************************************80
!
!! CN_JAC_01_1_SIZE sizes a precision 1 rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha. 
!
!    with -1.0 < alpha, -1.0 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_01_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_01_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine cn_jac_02_xiu ( n, alpha, beta, o, x, w )

!*****************************************************************************80
!
!! CN_JAC_02_XIU implements the Xiu precision 2 rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1.0 < alpha, -1.0 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) beta
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_02_XIU - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_02_XIU - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = ( alpha + beta + 2.0D+00 ) / 2.0D+00
  delta0 = ( alpha - beta ) / 2.0D+00
  c1 = 2.0D+00 * ( alpha + 1.0D+00 ) * ( beta + 1.0D+00 ) &
    / ( alpha + beta + 3.0D+00 ) / ( alpha + beta + 2.0D+00 )

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call c1_jac_monomial_integral ( alpha, beta, expon, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_jac_02_xiu_size ( n, alpha, beta, o )

!*****************************************************************************80
!
!! CN_JAC_02_XIU_SIZE sizes the Xiu precision 2 rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1.0 < alpha, -1.0 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_02_XIU_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_02_XIU_SIZE - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  o = n + 1

  return
end
subroutine cn_jac_monomial_integral ( n, alpha, beta, expon, value )

!*****************************************************************************80
!
!! CN_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over CN.
!
!  Discussion:
!
!    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1.0 < alpha, -1.0 < beta.
!
!    value = integral ( CN ) 
!      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i))^alpha (1+x(i))^beta dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X) in the weight factor.
!    - 1.0 < ALPHA.
!
!    Input, real ( kind = 8 ) BETA, the exponent of (1+X) in the weight factor.
!    - 1.0 < BETA.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  value = 1.0D+00
  do i = 1, n
    call c1_jac_monomial_integral ( alpha, beta, expon(i), value2 )
    value = value * value2
  end do

  return
end
subroutine cn_leg_01_1 ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_01_1 implements the midpoint rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call c1_leg_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  w(k) = volume

  return
end
subroutine cn_leg_01_1_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_01_1_SIZE sizes the midpoint rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine cn_leg_02_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_02_XIU implements the Xiu precision 2 rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = 1.0D+00
  delta0 = 0.0D+00
  c1 = 1.0D+00 / 3.0D+00

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call c1_leg_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_leg_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_02_XIU_SIZE sizes the Xiu precision 2 rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
subroutine cn_leg_03_1 ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_03_1 implements the Stroud rule CN:3-1 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!    The necessary treatment of the final coordinate of points when
!    N is odd seems to vary from what Stroud declares! 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit  none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call c1_leg_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  do j = 1, o

    i = 0

    do r = 1, ( n / 2 )

      arg = real ( ( 2 * r - 1 ) * ( j + 1 ), kind = 8 ) * pi &
        / real ( n, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) / sqrt ( 3.0D+00 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) / sqrt ( 3.0D+00 )

    end do
!
!  The following code does not correspond to what Stroud declares.
!
    if ( i < n ) then

      i = i + 1
      if ( n == 1 ) then
        x(i,j) =                r8_mop ( j + 1 ) / sqrt ( 3.0D+00 )
      else
        x(i,j) = sqrt ( 2.0 ) * r8_mop ( j + 1 ) / sqrt ( 3.0D+00 )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_leg_03_1_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_03_1_SIZE sizes the Stroud rule CN:3-1 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) CN_LEG_03_1_SIZE, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine cn_leg_03_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_03_XIU implements the Xiu precision 3 rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call c1_leg_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  do j = 1, o

    i = 0 
    do r = 1, n / 2
      arg = real ( ( 2 * r - 1 ) * j, kind = 8 ) * pi / real ( n, kind = 8 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) / sqrt ( 3.0D+00 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) / sqrt ( 3.0D+00 )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * r8_mop ( j ) / sqrt ( 3.0D+00 )
      if ( n == 1 ) then
        x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_leg_03_xiu_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine cn_leg_05_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_05_1 implements the Stroud rule CN:5-1 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N must be 4, 5, or 6.
!
!    Input, integer ( kind = 4 ) OPTION, is only used in case N = 4 or 5.
!    In that case, OPTION should be 1 or 2 to select the
!    two available variants of the rule.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) eta
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) k
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) option
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)
  real ( kind = 8 ) xsi

  if ( n < 4 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_LEG_05_1 - Fatal error!'
    write ( *, '(a)' ) '  The value of N must be 4, 5, or 6.'
    stop
  end if

  if ( n == 4 .or. n == 5 ) then
    if ( option < 1 .or. 2 < option ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CN_LEG_05_1 - Fatal error!'
      write ( *, '(a)' ) '  When N = 4 or 5, OPTION must be 1 or 2.'
      stop
    end if
  end if

  expon = 0
  call c1_leg_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  if ( n == 4 .and. option == 1 ) then

    eta    =   0.778984505799815D+00
    lambda =   1.284565137874656D+00
    xsi =    - 0.713647298819253D+00
    mu =     - 0.715669761974162D+00
    gamma =    0.217089151000943D+00
    a =        0.206186096875899D-01 * volume
    b =        0.975705820221664D-02 * volume
    c =        0.733921929172573D-01 * volume

  else if ( n == 4 .and. option == 2 ) then

    eta    =   0.546190755827425D+00
    lambda =   0.745069130115661D+00
    xsi =    - 0.413927294508700D+00
    mu =     - 0.343989637454535D+00
    gamma =    1.134017894600344D+00
    a =        0.853094758323323D-01 * volume
    b =        0.862099000096395D-01 * volume
    c =        0.116418206881849D-01 * volume

  else if ( n == 5 .and. option == 1 ) then

    eta    =   0.522478547481276D+00
    lambda =   0.936135175985774D+00
    xsi =    - 0.246351362101519D+00
    mu =     - 0.496308106093758D+00
    gamma =    0.827180176822930D+00
    a =        0.631976901960153D-01 * volume
    b =        0.511464127430166D-01 * volume
    c =        0.181070246088902D-01 * volume

  else if ( n == 5 .and. option == 2 ) then

    eta    =   0.798317301388741D+00
    lambda =   0.637344273885728D+00
    xsi =    - 0.455245909918377D+00
    mu =     - 1.063446229997311D+00
    gamma =    0.354482076665770D+00
    a =        0.116952384292206D-01 * volume
    b =        0.701731258612708D-01 * volume
    c =        0.137439132264426D-01 * volume

  else if ( n == 6 ) then

    eta    =   0.660225291773525D+00
    lambda =   1.064581294844754D+00
    xsi =      0.000000000000000D+00
    mu =     - 0.660225291773525D+00
    gamma =    0.660225291773525D+00
    a =        0.182742214532872D-01 * volume
    b =        0.346020761245675D-01 * volume
    c =        0.182742214532872D-01 * volume

  end if

  k = 0

  k = k + 1
  do i = 1, n
    x(i,k) = eta
  end do
  w(k) = a

  k = k + 1
  do i = 1, n
    x(i,k) = - eta
  end do
  w(k) = a

  do i1 = 1, n
    k = k + 1
    do i = 1, n
      x(i,k) = xsi
    end do
    x(i1,k) = lambda
    w(k) = b
  end do

  do i1 = 1, n
    k = k + 1
    do i = 1, n
      x(i,k) = - xsi
    end do
    x(i1,k) = - lambda
    w(k) = b
  end do

  do i1 = 1, n - 1
    do i2 = i1 + 1, n
      k = k + 1
      do i = 1, n
        x(i,k) = gamma
      end do
      x(i1,k) = mu
      x(i2,k) = mu
      w(k) = c
    end do
  end do

  do i1 = 1, n - 1
    do i2 = i1 + 1, n
      k = k + 1
      do i = 1, n
        x(i,k) = - gamma
      end do
      x(i1,k) = - mu
      x(i2,k) = - mu
      w(k) = c
    end do
  end do

  return
end
subroutine cn_leg_05_1_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_05_1_SIZE sizes the Stroud rule CN:5-1 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n * n + n + 2

  return o
end
subroutine cn_leg_05_2 ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_05_2 implements the Stroud rule CN:5-2 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 N^2 + 1.
!
!    The rule has precision P = 5.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N must be at least 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_LEG_05_2 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if

  expon = 0
  call c1_leg_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  b0 = real ( 25 * n * n - 115 * n + 162, kind = 8 ) * volume / 162.0D+00
  b1 = real ( 70 - 25 * n, kind = 8 ) * volume / 162.0D+00
  b2 = 25.0D+00 * volume / 324.0D+00

  r = sqrt ( 3.0D+00 / 5.0D+00 )

  k = 0

  k = k + 1
  do i = 1, n
    x(i,k) = 0.0D+00
  end do
  w(k) = b0

  do i1 = 1, n

    k = k + 1
    do i = 1, n
      x(i,k) = 0.0D+00
    end do
    x(i1,k) = + r
    w(k) = b1

    k = k + 1
    do i = 1, n
      x(i,k) = 0.0D+00
    end do
    x(i1,k) = - r
    w(k) = b1

  end do

  do i1 = 1, n - 1
    do i2 = i1 + 1, n

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      x(i1,k) = + r
      x(i2,k) = + r
      w(k) = b2

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      x(i1,k) = + r
      x(i2,k) = - r
      w(k) = b2

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      x(i1,k) = - r
      x(i2,k) = + r
      w(k) = b2

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      x(i1,k) = - r
      x(i2,k) = - r
      w(k) = b2

    end do
  end do

  return
end
subroutine cn_leg_05_2_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_05_2_SIZE sizes the Stroud rule CN:5-2 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 N^2 + 1.
!
!    The rule has precision P = 5.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n * n + 1

  return
end
subroutine cn_leg_monomial_integral ( n, expon, value )

!*****************************************************************************80
!
!! CN_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on CN.
!
!  Discussion:
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!    value = integral ( CN ) product ( 1 <= i <= n ) x(I)^expon(i) dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  value = 1.0D+00
  do i = 1, n
    call c1_leg_monomial_integral ( expon(i), value2 )
    value = value * value2
  end do

  return
end
subroutine en_her_01_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_01_1 implements the Stroud rule 1.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = volume

  return
end
subroutine en_her_01_1_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_01_1_SIZE sizes the Stroud rule 1.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine en_her_02_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_02_XIU implements the Xiu precision 2 rule for region EN_HER.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) delta0
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = 2.0D+00
  delta0 = 0.0D+00
  c1 = 1.0D+00

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine en_her_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_02_XIU_SIZE sizes the Xiu precision 2 rule for region EN_HER.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
subroutine en_her_03_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_03_1 implements the Stroud rule 3.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  a = volume / real ( o, kind = 8 )
  r = sqrt ( real ( n, kind = 8 ) / 2.0D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = a
    k = k + 1
    x(i,k) = + r
    w(k) = a
  end do

  return
end
subroutine en_her_03_1_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_03_1_SIZE sizes the Stroud rule 3.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine en_her_03_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_03_2 implements the Stroud rule 3.2 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^N.
!
!    The rule has precision P = 3.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  logical              more
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  a = volume / real ( o, kind = 8 )
  r = sqrt ( 0.5D+00 );

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - r
  w(k) = a
  more = .true.

  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < 0.0D+00 ) then
        k = k + 1
        x(1:n,k) = x(1:n,k-1)
        x(i,k)     = + r
        x(i+1:n,k) = - r
        w(k) = a
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_her_03_2_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_03_2_SIZE sizes the Stroud rule 3.2 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^N.
!
!    The rule has precision P = 3.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 ** n;

  return
end
subroutine en_her_03_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_03_XIU implements the Xiu precision 3 rule for region EN_HER.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  do j = 1, o

    i = 0 
    do r = 1, n / 2
      arg = real ( ( 2 * r - 1 ) * j, kind = 8 ) * pi / real ( n, kind = 8 )
      i = i + 1
      x(i,j) = cos ( arg )
      i = i + 1
      x(i,j) = sin ( arg )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j )
      if ( n == 1 ) then
        x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine en_her_03_xiu_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_03_XIU_SIZE sizes the Xiu precision 3 rule for region EN_HER.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine en_her_05_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_HER_05_1 implements the Stroud rule 5.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
!    the OPTION variable to 1 or 2.
!
!    Versions of this rule are only available for N = 2 through 7.
!
!    There is a typographical error in the reference.
!    For the second version of the rule for N = 2, the line
!      gamma =    0.313300683022281D+00
!    should read
!      gamma =    0.312200683022281D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    2 <= N <= 7.
!
!    Input, integer ( kind = 4 ) OPTION, selects option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) eta
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)
  real ( kind = 8 ) xsi

  if ( n < 2 .or. 7 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_05_1 - Fatal error!'
    write ( *, '(a)' ) '  2 <= N <= 7 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_05_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 5 .and. n /= 6 ) then 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_HER_05_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION = 2 requires N = 3, 5 or 6.'
      stop
    end if
  end if

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  if ( n == 2 ) then
    eta =      0.446103183094540D+00
    lambda =   0.136602540378444D+01
    xsi =    - 0.366025403784439D+00
    mu =       0.198167882945871D+01
    gamma =    0.000000000000000D+00
    a =        0.328774019778636D+00 * volume
    b =        0.833333333333333D-01 * volume
    c =        0.455931355469736D-02 * volume
  else if ( n == 3 .and. option == 1 ) then
    eta =      0.476731294622796D+00
    lambda =   0.935429018879534D+00
    xsi =    - 0.731237647787132D+00
    mu =       0.433155309477649D+00
    gamma =    0.266922328697744D+01
    a =        0.242000000000000D+00 * volume
    b =        0.810000000000000D-01 * volume
    c =        0.500000000000000D-02 * volume
!
!  The value of gamma that follows corrects an error in the reference.
!
  else if ( n == 3 .and. option == 2 ) then
    eta =      0.476731294622796D+00
    lambda =   0.128679320334269D+01
    xsi =    - 0.379873463323979D+00
    mu =     - 0.192386729447751D+01
    gamma =    0.312200683022281D+00
    a =        0.242000000000000D+00 * volume
    b =        0.810000000000000D-01 * volume
    c =        0.500000000000000D-02 * volume
  else if ( n == 4 ) then
    eta =      0.523945658287507D+00
    lambda =   0.119433782552719D+01
    xsi =    - 0.398112608509063D+00
    mu =     - 0.318569372920112D+00
    gamma =    0.185675837424096D+01
    a =        0.155502116982037D+00 * volume
    b =        0.777510584910183D-01 * volume
    c =        0.558227484231506D-02 * volume
  else if ( n == 5 .and. option == 1 ) then
    eta =      0.214972564378798D+01
    lambda =   0.464252986016289D+01
    xsi =    - 0.623201054093728D+00
    mu =     - 0.447108700673434D+00
    gamma =    0.812171426076311D+00
    a =        0.487749259189752D-03 * volume
    b =        0.487749259189752D-03 * volume
    c =        0.497073504444862D-01 * volume
  else if ( n == 5 .and. option == 2 ) then
    eta =      0.615369528365158D+00
    lambda =   0.132894698387445D+01
    xsi =    - 0.178394363877324D+00
    mu =     - 0.745963266507289D+00
    gamma =    0.135503972310817D+01
    a =        0.726415024414905D-01 * volume
    b =        0.726415024414905D-01 * volume
    c =        0.641509853510569D-02 * volume
  else if ( n == 6 .and. option == 1 ) then
    eta =      0.100000000000000D+01
    lambda =   0.141421356237309D+01
    xsi =      0.000000000000000D+00
    mu =     - 0.100000000000000D+01
    gamma =    0.100000000000000D+01
    a =        0.781250000000000D-02 * volume
    b =        0.625000000000000D-01 * volume
    c =        0.781250000000000D-02 * volume
  else if ( n == 6 .and. option == 2 ) then
    eta =      0.100000000000000D+01
    lambda =   0.942809041582063D+00
    xsi =    - 0.471404520791032D+00
    mu =     - 0.166666666666667D+01
    gamma =    0.333333333333333D+00
    a =        0.781250000000000D-02 * volume
    b =        0.625000000000000D-01 * volume
    c =        0.781250000000000D-02 * volume
  else if ( n == 7 ) then
    eta =      0.000000000000000D+00
    lambda =   0.959724318748357D+00
    xsi =    - 0.772326488820521D+00
    mu =     - 0.141214270131942D+01
    gamma =    0.319908106249452D+00
    a =        0.111111111111111D+00 * volume
    b =        0.138888888888889D-01 * volume
    c =        0.138888888888889D-01 * volume
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2 points.
!
  k = k + 1
  x(1:n,k) = - eta
  w(k) = a
  k = k + 1
  x(1:n,k) = + eta
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(1:n,k) = - xsi
    x(i,k) = - lambda
    w(k) = b
    k = k + 1
    x(1:n,k) = + xsi
    x(i,k) = + lambda
    w(k) = b
  end do
!
!  2 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(1:n,k) = - gamma
      x(i,k) = - mu
      x(j,k) = - mu
      w(k) = c
      k = k + 1
      x(1:n,k) = + gamma
      x(i,k) = + mu
      x(j,k) = + mu
      w(k) = c
    end do
  end do

  return
end
subroutine en_her_05_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_HER_05_1_SIZE sizes the Stroud rule 5.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
!    the OPTION variable to 1 or 2.
!
!    Versions of this rule are only available for N = 2 through 7.
!
!    There is a typographical error in the reference.
!    For the second version of the rule for N = 2, the line
!      gamma =    0.313300683022281D+00
!    should read
!      gamma =    0.312200683022281D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    2 <= N <= 7.
!
!    Input, integer ( kind = 4 ) OPTION, selects option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( n < 2 .or. 7 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_05_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  2 <= N <= 7 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_05_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 5 .and. n /= 6 ) then 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_HER_05_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION = 2 requires N = 3, 5 or 6.'
      stop
    end if
  end if

  o = n * n + n + 2

  return
end
subroutine en_her_05_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_05_2 implements the Stroud rule 5.2 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2 * N^2 + 1.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  a = 2.0D+00 * volume / real ( n + 2, kind = 8 )
  b = real ( 4 - n, kind = 8 ) * volume / 2.0D+00 &
    / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )
  c = volume / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )

  r = sqrt ( real ( n + 2, kind = 8 ) / 2.0D+00 )
  s = sqrt ( real ( n + 2, kind = 8 ) / 4.0D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = b
    k = k + 1
    x(i,k) = + r
    w(k) = b
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - s
      x(j,k) = - s
      w(k) = c
      k = k + 1
      x(i,k) = - s
      x(j,k) = + s
      w(k) = c
      k = k + 1
      x(i,k) = + s
      x(j,k) = - s
      w(k) = c
      k = k + 1
      x(i,k) = + s
      x(j,k) = + s
      w(k) = c
    end do
  end do

  return
end
subroutine en_her_05_2_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_05_2_SIZE sizes the Stroud rule 5.2 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2 * N^2 + 1.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n * n + 1

  return
end
subroutine en_her_05_3 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_05_3 implements the Stroud rule 5.3 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 3 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  logical              more
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_05_3 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  a = 4.0D+00 * volume / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )
  b = real ( ( n - 2 ) * ( n - 2 ), kind = 8 ) * volume &
    / real ( 2**n, kind = 8 ) / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )
  r = sqrt ( real ( n + 2, kind = 8 ) / 4.0D+00 )
  s = sqrt ( real ( n + 2, kind = 8 ) / 2.0D+00 / real ( n - 2, kind = 8 ) )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = a
    k = k + 1
    x(i,k) = + r
    w(k) = a
  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - s
  w(k) = b
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < 0.0D+00 ) then
        k = k + 1
        x(1:n,k) = x(1:n,k-1)
        x(i,k)     = + s
        x(i+1:n,k) = - s
        w(k) = b
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_her_05_3_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_05_3_SIZE sizes the Stroud rule 5.3 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 3 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_05_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  o = 2**n + 2 * n

  return
end
subroutine en_her_05_4 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_05_4 implements the Stroud rule 5.4 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) - 1.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              more
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  s = sqrt ( 0.5D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2^N + 2^(N-1) + 2^(N-2) + ... + 1 = 2^(N+1)-1 points.
!  but do the last point separately.
!
  do i = 1, n

    r = sqrt ( real ( i + 2, kind = 8 ) / 2.0D+00 )
    b = 2.0D+00 ** ( i - n ) * volume / real ( i + 1, kind = 8 ) &
      / real ( i + 2, kind = 8 )

    k = k + 1
    x(i,k) = - r
    x(i+1:n,k) = - s
    w(k) = b
    more = .true.
    do while ( more )
      more = .false.
      do j = n, i, -1
        if ( x(j,k) < 0.0D+00 ) then
          k = k + 1
          x(1:n,k) = x(1:n,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = b
          more = .true.
          exit
        end if
      end do
    end do

  end do
!
!  Last point.
!
  k = k + 1
! x(1:n,k) = 0.0
  w(k) = 2.0D+00 * volume / real ( n + 2, kind = 8 )

  return
end
subroutine en_her_05_4_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_05_4_SIZE sizes the Stroud rule 5.4 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) - 1.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 ** ( n + 1 ) - 1

  return
end
subroutine en_her_05_5 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_05_5 implements the Stroud rule 5.5 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = N * 2^N + 1.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There is a second version of this rule however it results in
!    complex abscissas, and so it has been disabled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              more
  real ( kind = 8 ) n_r8
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  n_r8 = real ( n, kind = 8 )

  a = 2.0D+00 * volume / ( n_r8 + 2.0D+00 )
  b =           volume / ( n_r8 + 2.0D+00 ) / ( 2.0D+00 ** n )

  option = 1

  if ( option == 1 ) then
    r = sqrt ( ( n_r8 + 2.0D+00 &
      + ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
      / 2.0D+00 / n_r8 )
    s = sqrt ( ( n_r8 + 2.0D+00 &
      -                      sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
      / 2.0D+00 / n_r8 )
  else if ( option == 2 ) then
    r = sqrt ( ( n_r8 + 2.0D+00 &
      - ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
      / 2.0D+00 / n_r8 )
    s = sqrt ( ( n_r8 + 2.0D+00 &
      +                      sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
      / 2.0D+00 / n_r8 )
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = a
!
!  N * 2^N points:
!  N choices for location of R, 2^N choices of sign pattern.
!
  do i = 1, n

    k = k + 1
    x(1:n,k) = - s
    x(i,k)   = - r
    w(k) = b

    more = .true.

    do while ( more )
      more = .false.
      do j = n, 1, -1
        if ( x(j,k) < 0.0D+00 ) then
          k = k + 1
          x(1:n,k) = x(1:n,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = b
          more = .true.
          exit
        end if
      end do
    end do

  end do

  return
end
subroutine en_her_05_5_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_05_5_SIZE sizes the Stroud rule 5.5 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = N * 2^N + 1.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There is a second version of this rule however it results in
!    complex abscissas, and so it has been disabled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n * 2 ** n + 1

  return
end
subroutine en_her_05_6 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_05_6 implements the Stroud rule 5.6 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = ( N + 1 ) * 2^N.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 5 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    5 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              more
  real ( kind = 8 ) n_r8
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_05_6 - Fatal error!'
    write ( *, '(a)' ) '  5 <= N is required.'
    stop
  end if

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  n_r8 = real ( n, kind = 8 )

  a = volume / ( 2.0D+00 ** n ) / ( n_r8 + 1.0D+00 )

  r = sqrt ( ( n_r8 - sqrt ( 2.0D+00 ) &
    + ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 1.0D+00 ) ) ) &
    / 2.0D+00 / n_r8 )
  s = sqrt ( ( n_r8 - sqrt ( 2.0D+00 ) &
    -                      sqrt ( 2.0D+00 * ( n_r8 + 1.0D+00 ) ) ) &
    / 2.0D+00 / n_r8 )
  t = sqrt ( ( 1.0D+00 + sqrt ( 2.0D+00 ) ) / 2.0D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  N * 2^N points.
!
  do i = 1, n

    k = k + 1
    x(1:n,k) = - s
    x(i,k)   = - r
    w(k) = a

    more = .true.

    do while ( more )
      more = .false.
      do j = n, 1, -1
        if ( x(j,k) < 0.0D+00 ) then
          k = k + 1
          x(1:n,k) = x(1:n,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = a
          more = .true.
          exit
        end if
      end do
    end do

  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - t
  w(k) = a
  more = .true.
  do while ( more )
    more = .false.
    do j = n, 1, -1
      if ( x(j,k) < 0.0D+00 ) then
        k = k + 1
        x(1:n,k) = x(1:n,k-1)
        x(j,k)     =   abs ( x(j,k) )
        x(j+1:n,k) = - abs ( x(j+1:n,k) )
        w(k) = a
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_her_05_6_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_05_6_SIZE sizes the Stroud rule 5.6 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = ( N + 1 ) * 2^N.
!
!    The rule has precision P = 5.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 5 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    5 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( n < 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_05_6_SIZE - Fatal error!'
    write ( *, '(a)' ) '  5 <= N is required.'
    stop
  end if

  o = ( 2 ** n ) * ( n + 1 )

  return
end
subroutine en_her_07_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_HER_07_1 implements the Stroud rule 7.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N^2 + 1.
!
!    The rule has precision P = 7.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of the rule, chosen by setting the
!    OPTION variable to 1 or 2.  
!
!    Option 1 is only valid for N = 3, 4, 6 or 7.
!    Option 2 is only valid for N = 3 or 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N = 3, 4, 6 or 7.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              more
  real ( kind = 8 ) n_r8
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_07_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 1 ) then
    if ( n /= 3 .and. n /= 4 .and. n /= 6 .and. n /= 7 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_HER_07_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION 1 requires N =  3, 4, 6 or 7.'
      stop
    end if
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_HER_07_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION 2 requires N =  3 or 4.'
      stop
    end if
  end if

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  n_r8 = real ( n, kind = 8 )

  if ( option == 1 ) then
    r = sqrt ( ( 3.0D+00 * ( 8.0D+00 - n_r8 ) - ( n_r8 - 2.0D+00 ) &
      * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 / ( 5.0D+00 - n_r8 ) )
    s = sqrt ( ( 3.0D+00 *             n_r8   -          2.0D+00   &
      * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 &
      / ( 3.0D+00 * n_r8 - 8.0D+00 ) )
    t = sqrt ( ( 6.0D+00 + sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 )
  else if ( option == 2 ) then
    r = sqrt ( ( 3.0D+00 * ( 8.0D+00 - n_r8 ) + ( n_r8 - 2.0D+00 ) &
      * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 / ( 5.0D+00 - n_r8 ) )
    s = sqrt ( ( 3.0D+00 *             n_r8   +          2.0D+00   &
      * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 &
      / ( 3.0D+00 * n_r8 - 8.0D+00 ) )
    t = sqrt ( ( 6.0D+00 - sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 )
  end if

  b = ( 8.0D+00 - n_r8 ) * volume / 8.0D+00 / r ** 6
  c = volume / 2.0D+00 ** ( n + 3 ) / s ** 6
  d = volume / 16.0D+00 / t ** 6
  a = volume - 2.0D+00 * n_r8 * b - 2.0D+00 ** n * c - 2.0D+00 * n_r8 &
    * ( n_r8 - 1.0D+00 ) * d

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = b
    k = k + 1
    x(i,k) = + r
    w(k) = b
  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - s
  w(k) = c
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < 0.0D+00 ) then
        k = k + 1
        x(1:n,k) = x(1:n,k-1)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = c
        more = .true.
        exit
      end if
    end do
  end do
!
!  2 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - t
      x(j,k) = - t
      w(k) = d
      k = k + 1
      x(i,k) = - t
      x(j,k) = + t
      w(k) = d
      k = k + 1
      x(i,k) = + t
      x(j,k) = - t
      w(k) = d
      k = k + 1
      x(i,k) = + t
      x(j,k) = + t
      w(k) = d
    end do
  end do

  return
end
subroutine en_her_07_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_HER_07_1_SIZE sizes the Stroud rule 7.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N^2 + 1.
!
!    The rule has precision P = 7.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of the rule, chosen by setting the
!    OPTION variable to 1 or 2.  
!
!    Option 1 is only valid for N = 3, 4, 6 or 7.
!    Option 2 is only valid for N = 3 or 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N = 3, 4, 6 or 7.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_07_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 1 ) then
    if ( n /= 3 .and. n /= 4 .and. n /= 6 .and. n /= 7 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_HER_07_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION 1 requires N =  3, 4, 6 or 7.'
      stop
    end if
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_HER_07_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION 2 requires N =  3 or 4.'
      stop
    end if
  end if

  o = 2 ** n + 2 * n ** 2 + 1

  return
end
subroutine en_her_07_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_HER_07_2 implements the Stroud rule 7.2 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) + 4 * N^2.
!
!    The rule has precision P = 7.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 3 <= N.
!
!    The reference has a typographical error in the description of this rule.
!    The formula:
!
!      (t,t,t,...,t)FS
!
!    should read
!
!      (t,t,0,...,0)FS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              more
  real ( kind = 8 ) n_r8
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) rho1
  real ( kind = 8 ) rho2
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_07_2 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  n_r8 = real ( n, kind = 8 )

  rho1 = sqrt ( ( n_r8 + 2.0D+00 - sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
    / 2.0D+00 )
  rho2 = sqrt ( ( n_r8 + 2.0D+00 + sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
    / 2.0D+00 )
  a1 = ( n_r8 + 2.0D+00 + sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) / 2.0D+00 &
    / ( n_r8 + 2.0D+00 )
  a2 = ( n_r8 + 2.0D+00 - sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) / 2.0D+00 &
    / ( n_r8 + 2.0D+00 )

  r = 1.0D+00
  s = sqrt ( 1.0D+00 / n_r8 )
  t = sqrt ( 0.5D+00 )
  b = ( 8.0D+00 - n_r8 ) * volume / n_r8 / ( n_r8 + 2.0D+00 ) &
    / ( n_r8 + 4.0D+00 )
  c = n_r8 ** 3 * volume / 2.0D+00 ** n / n_r8 / ( n_r8 + 2.0D+00 ) &
    / ( n_r8 + 4.0D+00 )
  d = 4.0D+00 * volume / n_r8 / ( n_r8 + 2.0D+00 ) / ( n_r8 + 4.0D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2 * 2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - rho1 * r
    w(k) = a1 * b
    k = k + 1
    x(i,k) = - rho2 * r
    w(k) = a2 * b
    k = k + 1
    x(i,k) = + rho1 * r
    w(k) = a1 * b
    k = k + 1
    x(i,k) = + rho2 * r
    w(k) = a2 * b
  end do
!
!  2 * 2^N points.
!
  k = k + 1
  x(1:n,k) = - rho1 * s
  w(k) = a1 * c
  k = k + 1
  x(1:n,k) = - rho2 * s
  w(k) = a2 * c
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < 0.0D+00 ) then
        k = k + 1
        x(1:n,k) =     x(1:n,k-2)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = a1 * c
        k = k + 1
        x(1:n,k) =     x(1:n,k-2)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = a2 * c
        more = .true.
        exit
      end if
    end do
  end do
!
!  2 * 4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - rho1 * t
      x(j,k) = - rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = - rho1 * t
      x(j,k) = + rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = + rho1 * t
      x(j,k) = - rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = + rho1 * t
      x(j,k) = + rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = - rho2 * t
      x(j,k) = - rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = - rho2 * t
      x(j,k) = + rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = + rho2 * t
      x(j,k) = - rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = + rho2 * t
      x(j,k) = + rho2 * t
      w(k) = a2 * d
    end do
  end do

  return
end
subroutine en_her_07_2_size ( n, o )

!*****************************************************************************80
!
!! EN_HER_07_2_SIZE sizes the Stroud rule 7.2 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) + 4 * N^2.
!
!    The rule has precision P = 7.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 3 <= N.
!
!    The reference has a typographical error in the description of this rule.
!    The formula:
!
!      (t,t,t,...,t)FS
!
!    should read
!
!      (t,t,0,...,0)FS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_07_2_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  o = 2 ** ( n + 1 ) + 4 * n * n

  return
end
subroutine en_her_07_3 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_HER_07_3 implements the Stroud rule 7.3 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
!
!    The rule has precision P = 7.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   45
!     4   97
!     5  181
!     6  305
!
!    The reference has a typographical error for N = 5, OPTION 1, B4:
!
!      -(1)0.736330882774831
!
!    should read
!
!      (-1)0.736330882774831
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) b4
  real ( kind = 8 ) b5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  logical              more
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_07_3 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_07_3 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  if ( n == 3 .and. option == 1 ) then
    u =    0.524647623275290D+00
    v =    0.165068012388578D+01
    b0 = - 0.166705761599566D+02
    b1 =   0.100296981655678D+02
    b2 =   0.161699246687754D+00
    b3 = - 0.604719151221535D+01
    b4 =   0.234381399489666D-01
    b5 =   0.417194501880647D+01
  else if ( n == 3 .and. option == 2 ) then
    u =    0.165068012388578D+01
    v =    0.524647623275290D+00
    b0 =   0.166705761599566D+02
    b1 =   0.178903161957074D+00
    b2 = - 0.665808190965810D+01
    b3 =   0.148361823143070D-01
    b4 =   0.229669852539758D+01
    b5 =   0.430097881732984D-02
  else if ( n == 4 .and. option == 1 ) then
    u  =   0.524647623275290D+00
    v  =   0.165068012388578D+01
    b0 = - 0.167539329651562D+03
    b1 =   0.687922329603575D+02
    b2 =   0.203518409659014D+00
    b3 = - 0.255075279116885D+02
    b4 =   0.415430214106084D-01
    b5 =   0.739458001434961D+01
  else if ( n == 4 .and. option == 2 ) then
    u =    0.165068012388578D+01
    v =    0.524647623275290D+00
    b0 =   0.688432856406677D+02
    b1 =   0.294997847268286D+00
    b2 = - 0.199427272118378D+02
    b3 =   0.110498755408511D-01
    b4 =   0.407079214570997D+01
    b5 =   0.762328646743931D-02
  else if ( n == 5 .and. option == 1 ) then
    u  =   0.524647623275290D+00
    v  =   0.165068012388578D+01
    b0 = - 0.826940846964452D+03
    b1 =   0.264779097660331D+03
    b2 =   0.213460812375320D+00
    b3 = - 0.714240197186780D+02
    b4 =   0.736330882774831D-01
    b5 =   0.131065518222629D+02
  else if ( n == 5 .and. option == 2 ) then
    u =    0.165068012388578D+01
    v =    0.524647623275290D+00
    b0 =   0.220502344940121D+03
    b1 =   0.537746975313769D+00
    b2 = - 0.497781460739792D+02
    b3 = - 0.743845245712926D-02
    b4 =   0.721529121489956D+01
    b5 =   0.135119234557687D-01
  else if ( n == 6 .and. option == 1 ) then
    u  =   0.524647623275290D+00
    v  =   0.165068012388578D+01
    b0 = - 0.309679578630802E+04
    b1 =   0.815423321880237D+03
    b2 =   0.117326937169073D+00
    b3 = - 0.173057295296448D+03
    b4 =   0.130511250871491D+00
    b5 =   0.232307582494626D+02
  else if ( n == 6 .and. option == 2 ) then
    u =    0.165068012388578D+01
    v =    0.524647623275290D+00
    b0 =   0.616293651884027D+03
    b1 =   0.107529736766179D+01
    b2 = - 0.113807008098269D+03
    b3 = - 0.610828352270520D-01
    b4 =   0.127887706992535D+02
    b5 =   0.239492607623178D-01
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b3
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b4
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b5
      end do
    end do
  end do

  return
end
subroutine en_her_07_3_size ( n, option, o )

!*****************************************************************************80
!
!! EN_HER_07_3_SIZE sizes the Stroud rule 7.3 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
!
!    The rule has precision P = 7.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   45
!     4   97
!     5  181
!     6  305
!
!    The reference has a typographical error for N = 5, OPTION 1, B4:
!
!      -(1)0.736330882774831
!
!    should read
!
!      (-1)0.736330882774831
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_07_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_07_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  o = ( 4 * n ** 3 + 8 * n + 3 ) / 3

  return
end
subroutine en_her_09_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_HER_09_1 implements the Stroud rule 9.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
!
!    The rule has precision P = 9.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the 
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   77
!     4  193
!     5  421
!     6  825
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) b4
  real ( kind = 8 ) b5
  real ( kind = 8 ) b6
  real ( kind = 8 ) b7
  real ( kind = 8 ) b8
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical              more
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_09_1 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_09_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  if ( n == 3 ) then
    u =    0.202018287045609D+01
    v =    0.958572464613819D+00
    b0 =   0.676448734429924D+00
    b1 =   0.511989106291551D-02
    b2 =   0.448595723493744D+00
    b3 =   0.235223454595606D-03
    b4 =   0.915390713080005D-01
    b5 =   0.139208199920793D-01
    b6 =   0.235223454595606D-03
    b7 =   0.915390713080008D-01
    b8 =   0.000000000000000D+00
  else if ( n == 4 .and. option == 1 ) then
    u =    0.202018287045609D+01
    v =    0.958572464613819D+00
    b0 = - 0.860452945007048D+00
    b1 = - 0.405511998533795D-01
    b2 =   0.107026475449715D+01
    b3 =   0.138974239307092D-03
    b4 = - 0.162248779448181D+00
    b5 =   0.246740110027234D-01
    b6 =   0.138974239307094D-03
    b7 =   0.162248779448181D+00
    b8 =   0.138974239307094D-03
  else if ( n == 4 .and. option == 2 ) then
    u =    0.958572464613819D+00
    v =    0.202018287045609D+01
    b0 =   0.265029088766810D-02
    b1 =   0.637601342635332D+00
    b2 = - 0.394394059389228D-01
    b3 =   0.540829264827264D-01
    b4 = - 0.416922717921281D-03
    b5 =   0.246740110027234D-01
    b6 =   0.540829264827270D-01
    b7 =   0.416922717921281D-03
    b8 =   0.540829264827269D-01
  else if ( n == 5 .and. option == 1 ) then
    u =    0.202018287045609D+01
    v =    0.958572464613819D+00
    b0 = - 0.827347006200826D+01
    b1 = - 0.160820174530905D+00
    b2 =   0.353499863758467D+01
    b3 =   0.738976276909564D-03
    b4 = - 0.862735421812943D+00
    b5 =   0.437335458190621D-01
    b6 = - 0.246325425636523D-03
    b7 =   0.287578473937648D+00
    b8 =   0.246325425636523D-03
  else if ( n == 5 .and. option == 2 ) then
    u =    0.958572464613819D+00
    v =    0.202018287045609D+01
    b0 = - 0.624416791055272D+00
    b1 =   0.467494915583104D+00
    b2 = - 0.152937760910536D+00
    b3 =   0.287578473937646D+00
    b4 = - 0.221692883072871D-02
    b5 =   0.437335458190621D-01
    b6 = - 0.958594913125490D-01
    b7 =   0.738976276909568D-03
    b8 =   0.958594913125492D-01
  else if ( n == 6 .and. option == 1 ) then
    u =    0.202018287045609D+01
    v =    0.958572464613819D+00
    b0 = - 0.361840434143098D+02
    b1 = - 0.447936529138517D+00
    b2 =   0.112077863004144D+02
    b3 =   0.392940404320855D-02
    b4 = - 0.254859786784158D+01
    b5 =   0.775156917007496D-01
    b6 = - 0.130980134773619D-02
    b7 =   0.509719573568315D+00
    b8 =   0.436600449245395D-03
  else if ( n == 6 .and. option == 2 ) then
    u =    0.958572464613819D+00
    v =    0.202018287045609D+01
    b0 =   0.448873836333650D+01
    b1 = - 0.238473566140736D+01
    b2 = - 0.413008493198885D+00
    b3 =   0.152915872070494D+01
    b4 = - 0.654900673868093D-02
    b5 =   0.775156917007496D-01
    b6 = - 0.509719573568314D+00
    b7 =   0.130980134773618D-02
    b8 =   0.169906524522772D+00
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b3
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b4
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = - u
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = + u
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = + u
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = - u
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = + u
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = - u
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = + u
      w(k) = b5
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b6
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b7
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
        end do
      end do
    end do
  end do

  return
end
subroutine en_her_09_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_HER_09_1_SIZE sizes the Stroud rule 9.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
!
!    The rule has precision P = 9.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the 
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   77
!     4  193
!     5  421
!     6  825
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_09_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_09_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  o = ( 2 * n ** 4 - 4 * n ** 3 + 22 * n ** 2 - 8 * n + 3 ) / 3

  return
end
subroutine en_her_11_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_HER_11_1 implements the Stroud rule 11.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order 
!
!      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
!
!    The rule has precision P = 11.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 5.
!
!     N    O
!    __  ___
!     3  151
!     4  417
!     5  983
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 5.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) b4
  real ( kind = 8 ) b5
  real ( kind = 8 ) b6
  real ( kind = 8 ) b7
  real ( kind = 8 ) b8
  real ( kind = 8 ) b9
  real ( kind = 8 ) b10
  real ( kind = 8 ) b11
  real ( kind = 8 ) b12
  real ( kind = 8 ) b13
  real ( kind = 8 ) b14
  real ( kind = 8 ) b15
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i5
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical              more
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w2
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 .or. 5 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_11_1 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 5 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_11_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  if ( n == 3 .and. option == 1 ) then
    u =     0.235060497367449D+01
    v =     0.436077411927617D+00
    w2 =    0.133584907401370D+01
    b0 =  - 0.881591029957858D+01
    b1 =  - 0.751996143360650D-01
    b2 =    0.621743189471515D+01
    b3 =    0.241426451456494D+00
    b4 =  - 0.120709739276065D-02
    b5 =  - 0.427751221210138D+01
    b6 =    0.550169924840163D-01
    b7 =    0.237084999634707D-01
    b8 =  - 0.169791992887741D-02
    b9 =  - 0.252266276123350D-04
    b10 =   0.326777873717691D+01
    b11 =   0.968469949206802D-02
    b12 =   0.789754514877422D-03
    b13 =   0.000000000000000D+00
    b14 =   0.000000000000000D+00
    b15 =   0.000000000000000D+00
  else if ( n == 3 .and. option == 2 ) then
    u =     0.235060497367449D+01
    v =     0.133584907401370D+01
    w2 =    0.436077411927617D+00
    b0 =  - 0.141214037032900D+02
    b1 =  - 0.803730274707282D-01
    b2 =    0.235546545595906D+00
    b3 =    0.888123191556611D+01
    b4 =    0.142467131155533D-03
    b5 =    0.582993124006494D-01
    b6 =  - 0.561099173155661D+01
    b7 =  - 0.204028691521686D-02
    b8 =    0.252880089932256D-01
    b9 =  - 0.814378678627283D-04
    b10 =   0.804353953375146D-02
    b11 =   0.393451849690453D+01
    b12 =   0.171183493169724D-03
    b13 =   0.000000000000000D+00
    b14 =   0.000000000000000D+00
    b15 =   0.000000000000000D+00
  else if ( n == 4 .and. option == 1 ) then
    u =     0.235060497367449D+01
    v =     0.436077411927617D+00
    w2 =    0.133584907401370D+01
    b0 =    0.241502736147339D+03
    b1 =  - 0.196095938531478D+00
    b2 =  - 0.128675737999280D+03
    b3 =    0.307568784278696D+00
    b4 =  - 0.480908422319460D-02
    b5 =    0.698087019367085D+02
    b6 =    0.631837143743771D-01
    b7 =    0.392226151971179D-01
    b8 =  - 0.300948471646799D-02
    b9 =  - 0.650235306755170D-04
    b10 = - 0.386951974646715D+02
    b11 =   0.171656829095787D-01
    b12 =   0.139980343116450D-02
    b13 =   0.101552487093372D-04
    b14 =   0.222435922356439D+02
    b15 =   0.000000000000000D+00
  else if ( n == 4 .and. option == 2 ) then
    u =     0.235060497367449D+01
    v =     0.133584907401370D+01
    w2 =    0.436077411927617D+00
    b0 =  - 0.151944464736584D+03
    b1 =  - 0.223498438689039D+00
    b2 =    0.243574919068010D+00
    b3 =    0.634373877008693D+02
    b4 =  - 0.782065187814018D-04
    b5 =    0.911833754536616D-01
    b6 =  - 0.238927288245914D+02
    b7 =  - 0.422314408318853D-02
    b8 =    0.448218289217760D-01
    b9 =  - 0.138053374667391D-03
    b10 =   0.607473265800655D-02
    b11 =   0.697375246129742D+01
    b12 =   0.303414841680135D-03
    b13 = - 0.314574391771792D-05
    b14 =   0.409103498175100D-02
    b15 =   0.000000000000000D+00
  else if ( n == 5 .and. option == 1 ) then
    u =     0.235060497367449D+01
    v =     0.436077411927617D+00
    w2 =    0.133584907401370D+01
    b0 =    0.255885269311763E+04
    b1 =  - 0.439598677491526D+00
    b2 =  - 0.106541406144610E+04
    b3 =    0.453540909054264D+00
    b4 =  - 0.132100905623778D-01
    b5 =    0.418606568954203D+03
    b6 =    0.511394563043680D-01
    b7 =    0.645581013845604D-01
    b8 =  - 0.533417277494500D-02
    b9 =  - 0.137981626254496D-03
    b10 = - 0.147436933189884D+03
    b11 =   0.304253807765057D-01
    b12 =   0.248108698207828D-02
    b13 =   0.113652094546015D-04
    b14 =   0.394257407160391D+02
    b15 =   0.331725011358320D-05
  else if ( n == 5 .and. option == 2 ) then
    u =     0.235060497367449D+01
    v =     0.133584907401370D+01
    w2 =    0.436077411927617D+00
    b0 =  - 0.761305347548192D+03
    b1 =  - 0.536360805019297D+00
    b2 =    0.110669832078736D+00
    b3 =    0.246421088923968D+03
    b4 =  - 0.773649327968607D-03
    b5 =    0.169088641205970D+00
    b6 =  - 0.670700680243651D+02
    b7 =  - 0.856090560229205D-02
    b8 =    0.794446232770302D-01
    b9 =  - 0.220272863263544D-03
    b10 = - 0.373515812228225D-02
    b11 =   0.123606544052884D+02
    b12 =   0.537788804557843D-03
    b13 = - 0.122101861480881D-04
    b14 =   0.725117070759373D-02
    b15 =   0.331725011358320D-05
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - w2
    w(k) = b3
    k = k + 1
    x(i,k) = + w2
    w(k) = b3
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b4
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b4
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b4
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b4
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b5
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - w2
      x(j,k) = - w2
      w(k) = b6
      k = k + 1
      x(i,k) = - w2
      x(j,k) = + w2
      w(k) = b6
      k = k + 1
      x(i,k) = + w2
      x(j,k) = - w2
      w(k) = b6
      k = k + 1
      x(i,k) = + w2
      x(j,k) = + w2
      w(k) = b6
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - v
      w(k) = b7
      k = k + 1
      x(i,k) = - u
      x(j,k) = + v
      w(k) = b7
      k = k + 1
      x(i,k) = + u
      x(j,k) = - v
      w(k) = b7
      k = k + 1
      x(i,k) = + u
      x(j,k) = + v
      w(k) = b7
      k = k + 1
      x(i,k) = - v
      x(j,k) = - u
      w(k) = b7
      k = k + 1
      x(i,k) = - v
      x(j,k) = + u
      w(k) = b7
      k = k + 1
      x(i,k) = + v
      x(j,k) = - u
      w(k) = b7
      k = k + 1
      x(i,k) = + v
      x(j,k) = + u
      w(k) = b7
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - w2
      w(k) = b8
      k = k + 1
      x(i,k) = - u
      x(j,k) = + w2
      w(k) = b8
      k = k + 1
      x(i,k) = + u
      x(j,k) = - w2
      w(k) = b8
      k = k + 1
      x(i,k) = + u
      x(j,k) = + w2
      w(k) = b8
      k = k + 1
      x(i,k) = - w2
      x(j,k) = - u
      w(k) = b8
      k = k + 1
      x(i,k) = - w2
      x(j,k) = + u
      w(k) = b8
      k = k + 1
      x(i,k) = + w2
      x(j,k) = - u
      w(k) = b8
      k = k + 1
      x(i,k) = + w2
      x(j,k) = + u
      w(k) = b8
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b9
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b10
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - w2
        x(j,k) = - w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = - w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = + w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = + w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = - w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = - w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = + w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = + w2
        x(l,k) = + w2
        w(k) = b11
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 2 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b12
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
        end do
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
        end do
      end do
    end do
  end do
!
!  All quintuples UUUUU with 32 sign combinations.
!
  do i1 = 1, n - 4
    do i2 = i1 + 1, n - 3
      do i3 = i2 + 1, n - 2
        do i4 = i3 + 1, n - 1
          do i5 = i4 + 1, n
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
          end do
        end do
      end do
    end do
  end do

  return
end
subroutine en_her_11_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_HER_11_1_SIZE sizes the Stroud rule 11.1 for region EN_HER.
!
!  Discussion:
!
!    The rule has order 
!
!      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
!
!    The rule has precision P = 11.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 5.
!
!     N    O
!    __  ___
!     3  151
!     4  417
!     5  983
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 5.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( n < 3 .or. 5 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_11_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 5 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_11_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  o = ( 4 * n ** 5 - 20 * n ** 4 + 140 * n ** 3 - 130 * n ** 2 &
    + 96 * n + 15 ) / 15


  return
end
subroutine en_her_monomial_integral ( n, expon, value )

!*****************************************************************************80
!
!! EN_HER_MONOMIAL_INTEGRAL evaluates monomial integrals in EN_HER.
!
!  Discussion:
!
!    EXPON is the set of polynomial exponents.
!
!    EN_HER is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The integral to be evaluated is
!
!      value = integral ( EN ) x(1)^expon(1) * x(2)^expon(2) * ... 
!        * x(n)^expon(n) * w(x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(N), the polynomial exponents.
!    0 <= ALPHA(*).
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) value

  if ( any ( expon(1:n) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_HER_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  Some EXPON(I) < 0.'
    stop
  else if ( any ( mod ( expon(1:n), 2 ) == 1 ) ) then
    value = 0.0D+00
  else
    value = 1.0D+00
    do i = 1, n
      arg = ( real ( expon(i) + 1, kind = 8 ) ) / 2.0D+00
      value = value * r8_gamma ( arg )
    end do
  end if

  return
end
subroutine ep1_glg_monomial_integral ( expon, alpha, exact )

!*****************************************************************************80
!
!! EP1_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EP1.
!
!  Discussion:
!
!    EP1_GLG is the interval [0,+oo) with generalized Laguerre weight function:
!
!      w(alpha;x) = x^alpha exp ( - x )
!
!    value = integral ( 0 <= x < +oo ) x^expon x^alpha exp ( - x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2010
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
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EP1_GLG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( expon < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EP1_GLG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  EXPON < 0.'
    stop
  end if

  arg = alpha + real ( expon + 1, kind = 8 )

  exact = r8_gamma ( arg )

  return
end
subroutine ep1_lag_monomial_integral ( expon, value )

!*****************************************************************************80
!
!! EP1_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EP1.
!
!  Discussion:
!
!    EP1_LAG is the interval [0,+oo) with exponential or Laguerre 
!    weight function:
!
!      w(x) = exp ( - x )
!
!    value = integral ( 0 <= x < +oo ) x^expon exp ( - x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2010
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
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) value

  if ( expon < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EP1_LAG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  EXPON < 0.'
    stop
  end if

  value = r8_factorial ( expon )

  return
end
subroutine epn_glg_00_1 ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! EPN_GLG_00_1 implements the "midpoint rule" for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call ep1_glg_monomial_integral ( expon, alpha, volume_1d )
  volume = volume_1d ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = 1.0D+00
  w(k) = volume

  return
end
subroutine epn_glg_00_1_size ( n, alpha, o )

!*****************************************************************************80
!
!! EPN_GLG_00_1_SIZE sizes the midpoint rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_00_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine epn_glg_01_1 ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! EPN_GLG_01_1 implements a precision 1 rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_01_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call ep1_glg_monomial_integral ( expon, alpha, volume_1d )
  volume = volume_1d ** n

  expon = 1
  call ep1_glg_monomial_integral ( expon, alpha, value2 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = value2 / volume_1d
  w(k) = volume

  return
end
subroutine epn_glg_01_1_size ( n, alpha, o )

!*****************************************************************************80
!
!! EPN_GLG_01_1_SIZE sizes a precision 1 rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_01_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine epn_glg_02_xiu ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! EPN_GLG_02_XIU implements the Xiu precision 2 rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_02_XIU - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = - 1.0D+00
  delta0 = alpha + 1.0D+00
  c1 = - alpha - 1.0D+00

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call ep1_glg_monomial_integral ( expon, alpha, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine epn_glg_02_xiu_size ( n, alpha, o )

!*****************************************************************************80
!
!! EPN_GLG_02_XIU_SIZE sizes the Xiu precision 2 rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_02_XIUI_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = n + 1

  return
end
subroutine epn_glg_monomial_integral ( n, expon, alpha, value )

!*****************************************************************************80
!
!! EPN_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EPN.
!
!  Discussion:
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!    value = integral ( EPN ) 
!      product ( 1 <= i <= n ) x(I)^expon(i) x(i)^alpha exp ( - x(i) ) dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  value = 1.0D+00
  do i = 1, n
    call ep1_glg_monomial_integral ( expon(i), alpha, value2 )
    value = value * value2
  end do

  return
end
subroutine epn_lag_00_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EPN_LAG_00_1 implements the "midpoint rule" for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call ep1_lag_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = 1.0D+00
  w(k) = volume

  return
end
subroutine epn_lag_00_1_size ( n, o )

!*****************************************************************************80
!
!! EPN_LAG_00_1_SIZE sizes the midpoint rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine epn_lag_01_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EPN_LAG_01_1 implements a precision 1 rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call ep1_lag_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  expon = 1
  call ep1_lag_monomial_integral ( expon, value2 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = value2 / volume_1d
  w(k) = volume

  return
end
subroutine epn_lag_01_1_size ( n, o )

!*****************************************************************************80
!
!! EPN_LAG_01_1_SIZE sizes a precision 1 rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine epn_lag_02_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! EPN_LAG_02_XIU implements the Xiu precision 2 rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = - 1.0D+00
  delta0 = 1.0D+00
  c1 = - 1.0D+00

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call ep1_lag_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine epn_lag_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! EPN_LAG_02_XIU_SIZE sizes the Xiu precision 2 rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
subroutine epn_lag_monomial_integral ( n, expon, value )

!*****************************************************************************80
!
!! EPN_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EPN.
!
!  Discussion:
!
!    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!    value = integral ( EPN ) 
!      product ( 1 <= i <= n ) x(I)^expon(i) exp ( -x(i) ) dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  value = 1.0D+00
  do i = 1, n
    call ep1_lag_monomial_integral ( expon(i), value2 )
    value = value * value2
  end do

  return
end
subroutine gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )

!*****************************************************************************80
!
!! GW_02_XIU implements the Golub-Welsch version of the Xiu rule.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    It is assumed that the integral is over an N-dimensional region,
!    and has the form
!
!      Integral f(x) w(x) dx
!
!    where w(x) is separable into identical and independent components:
!
!      w(x) = v(x1) * v(x2) * ... * v(xn)
!
!    Associated with the weight function v(x), we assume there is a
!    family of orthogonal polynomials satisfying a three-term recurrence
!    of the form:
!
!      x P(n,x) = An * P(n+1,x) + Bn * P(n,x) + Cn * P(n-1,x)
!
!    with P(0,x) = 1, and P(-1,x) = 0.
!
!    This routine can construct the desired quadrature rule by knowing
!    the values of C1, used in the definition of P2, the values
!    GAMMA0 = 1/A0 and DELTA0 = - B0/A0, for which it is the case that
!    P(1,X) = GAMMA0 * X + DELTA0, and the value of VOLUME_1D, that is,
!    the 1D integral of v(x) over the region.
!
!    Note the values for the following standard polynomial families:
!
!    Chebyshev Type 1
!      V(X) =      1 / sqrt ( 1 - X^2 )
!      Interval =  [-1,+1]
!      GAMMA0 =    1.0
!      DELTA0 =    0.0
!      C1 =        1/2
!      VOLUME_1D = PI
!
!    Chebyshev Type 2
!      V(X) =      sqrt ( 1 - X^2 )
!      Interval =  [-1,+1]
!      GAMMA0 =    2.0
!      DELTA0 =    0.0
!      C1 =        1/2
!      VOLUME_1D = PI / 2
!
!    Gegenbauer
!      V(X) =      ( 1 - X^2 )^A
!      Interval =  [-1,+1]
!      GAMMA0 =    2 * A + 1
!      DELTA0 =    0.0
!      C1 =        ( 2 * A + 1 ) / ( 2 A + 3 )
!      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
!
!    Gegenbauer* (Removes singularity at ALPHA = -0.5):
!      V(X) =      ( 1 - X^2 )^A
!      Interval =  [-1,+1]
!      GAMMA0 =    1
!      DELTA0 =    0.0
!      C1 =        1 / ( 2 A + 3 )
!      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
!
!    Generalized Hermite
!      V(X) = |x|^A exp ( - x^2 )
!      Interval = (-oo,+oo)
!      GAMMA0 =    2
!      DELTA0 =    0
!      C1 =        2+2A
!      VOLUME_1D = Gamma((A+1)/2)
!
!    Generalized Laguerre
!      V(X) =       x^A exp ( - x )
!      Interval =  [0,+oo)
!      GAMMA0 =    -1.0
!      DELTA0 =     A+1.0
!      C1 =        -A-1.0
!      VOLUME_1D =  Gamma(A+1)
!
!    Hermite (physicist)
!      V(X) =      exp ( - x^2 )
!      Interval =  (-oo,+oo)
!      GAMMA0 =    2.0
!      DELTA0 =    0.0
!      C1 =        1.0
!      VOLUME_1D = sqrt ( PI )
!
!    Hermite (probabilist)
!      V(X) =      exp ( - x^2 / 2 )
!      Interval =  (-oo,+oo)
!      GAMMA0 =    1.0
!      DELTA0 =    0.0
!      C1 =        1.0
!      VOLUME_1D = sqrt ( 2 PI )
!
!    Jacobi
!      V(X) =      (1-x)^A (1+x)^B
!      Interval =  [-1,+1]
!      GAMMA0 =    (A+B+2)/2  
!      DELTA0 =    (A-B)/2
!      C1 =        2(A+1)(B+1)/(A+B+3)/(A+B+2)
!      VOLUME_1D = 2^(A+B+1) * Gamma(A+1) * Gamma(B+1) / ( A+B+1) / Gamma(A+B+1)
!
!    Laguerre
!      V(X) =       exp ( - x )
!      Interval =  [0,+oo)
!      GAMMA0 =    -1.0
!      DELTA0 =     1.0
!      C1 =        -1.0
!      VOLUME_1D =  1.0
!
!    Legendre
!      V(X) =      1.0
!      Interval =  [-1,+1]
!      GAMMA0 =    1.0
!      DELTA0 =    0.0
!      C1 =        1/3
!      VOLUME_1D = 2.0
!                                  
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Input, real ( kind = 8 ) GAMMA0, the ratio 1 / A0.
!
!    Input, real ( kind = 8 ) DELTA0, the ratio B0 / A0.
!
!    Input, real ( kind = 8 ) C1, the coefficient of P(0,X) in 
!    the definition of P(2,X).
!
!    Input, real ( kind = 8 ) VOLUME_1D, the 1D integral of V(X).
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) delta0
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  do j = 1, o

    i = 0
    do r = 1, ( n / 2 )
      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do
!
!  Adjust for the GW rule.
!
  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0
!
!  The weights are equal.
!
  w(1:o) = volume_1d ** n  / real ( o, kind = 8 )

  return
end
subroutine gw_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! GW_02_XIU_SIZE sizes the Golub Welsch version of the Xiu rule.
!
!  Discussion:
!
!    The rule has order O = N + 1;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
subroutine monomial_value ( dim_num, point_num, x, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the point coordinates.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the value of the monomial.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) expon(dim_num)
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  value(1:point_num) = 1.0D+00

  do dim = 1, dim_num
    if ( 0 /= expon(dim) ) then
      value(1:point_num) = value(1:point_num) * x(dim,1:point_num)**expon(dim)
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
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program provided that the copyright
!    is acknowledged.
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
  logical              l0
  logical              l1
  logical              l2
  logical              l3
  logical              l4
  logical              l5
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
    write ( *, '(a,g14.6)' ) '  C = ', c
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    stop
  end if

  if ( l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    write ( *, '(a)' ) '  1 - X < 0, C - A - B < 0.'
    write ( *, '(a,g14.6)' ) '  A = ', a
    write ( *, '(a,g14.6)' ) '  B = ', b
    write ( *, '(a,g14.6)' ) '  C = ', c
    write ( *, '(a,g14.6)' ) '  X = ', x
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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nq
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

  x = xx
  w = abs ( x )
  aug = 0.0D+00
!
!  Check for valid arguments, then branch to appropriate algorithm.
!
  if ( xmax1 <= - x .or. w < xmin1 ) then

    if ( 0.0D+00 < x ) then
      r8_psi = - xinf
    else
      r8_psi = xinf
    end if

    return
  end if

  if ( x < 0.5D+00 ) then
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
    if ( w <= xsmall ) then

      aug = - 1.0D+00 / x
!
!  Argument reduction for cotangent.
!
    else

      if ( x < 0.0D+00 ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - real ( int ( w ), kind = 8 )
      nq = int ( w * 4.0D+00 )
      w = 4.0D+00 * ( w - real ( nq, kind = 8 ) * 0.25D+00 )
!
!  W is now related to the fractional part of 4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
      n = nq / 2

      if ( n + n /= nq ) then
        w = 1.0D+00 - w
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
        if ( z == 0.0D+00 ) then

          if ( 0.0D+00 < x ) then
            r8_psi = - xinf
          else
            r8_psi = xinf
          end if

          return
        end if

        aug = sgn * ( 4.0D+00 / tan ( z ) )

      else

        aug = sgn * ( 4.0D+00 * tan ( z ) )

      end if

    end if

    x = 1.0D+00 - x

  end if
!
!  0.5 <= X <= 3.0.
!
  if ( x <= 3.0D+00 ) then

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
    w = 1.0D+00 / ( x * x )
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do
    aug = ( upper + p2(7) ) / ( den + q2(6) ) - 0.5D+00 / x + aug
  end if

  r8_psi = aug + log ( x )

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
