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
subroutine pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u )

!*****************************************************************************80
!
!! PCE_ODE_HERMITE applies the polynomial chaos expansion to a scalar ODE.
!
!  Discussion:
!
!    The deterministic equation is
!
!      du/dt = - alpha * u,
!      u(0) = u0
!
!    In the stochastic version, it is assumed that the decay coefficient
!    ALPHA is a Gaussian random variable with mean value ALPHA_MU and variance
!    ALPHA_SIGMA^2.
!
!    The exact expected value of the stochastic equation will be
!
!      u(t) = u0 * exp ( t^2/2)
!
!    This should be matched by the first component of the polynomial chaos
!    expansion.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TI, TF, the initial and final times.
!
!    Input, integer ( kind = 4 ) NT, the number of output points.
!
!    Input, real ( kind = 8 ) UI, the initial condition.
!
!    Input, integer ( kind = 4 ) NP, the degree of the expansion.  Polynomials 
!    of degree 0 through NP will be used.
!
!    Input, real ( kind = 8 ) ALPHA_MU, ALPHA_SIGMA, the mean and standard 
!    deviation of the decay coefficient.
!
!    Output, real ( kind = 8 ) T(0:NT), U(0:NT,0:NP), the times and the PCE 
!    coefficients at the successive time steps.
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) nt

  real ( kind = 8 ) alpha_mu
  real ( kind = 8 ) alpha_sigma
  real ( kind = 8 ) dp
  real ( kind = 8 ) dt
  real ( kind = 8 ) he_double_product_integral
  real ( kind = 8 ) he_triple_product_integral
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) t(0:nt)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) term
  real ( kind = 8 ) tf
  real ( kind = 8 ) ti
  real ( kind = 8 ) tp
  real ( kind = 8 ) u(0:nt,0:np)
  real ( kind = 8 ) u1(0:np)
  real ( kind = 8 ) u2(0:np)
  real ( kind = 8 ) ui

  dt = ( tf - ti ) / real ( nt, kind = 8 )
!
!  Set the PCE coefficients for the initial time.
!
  t1 = ti
  u1(0) = ui
  u1(1:np) = 0.0D+00
!
!  Copy into the output arrays.
!
  t(0) = t1
  u(0,0:np) = u1(0:np)
!
!  Time integration.
!
  do it = 1, nt

    t2 = ( real ( nt - it, kind = 8 ) * ti   &
         + real (      it, kind = 8 ) * tf ) &
         / real ( nt,      kind = 8 )

    do k = 0, np

      dp = he_double_product_integral ( k, k )

      term = - alpha_mu * u1(k)

      i = 1
      do j = 0, np
        tp = he_triple_product_integral ( i, j, k )
        term = term - alpha_sigma * u1(j) * tp / dp
      end do

      u2(k) = u1(k) + dt * term

    end do
!
!  Prepare for next step.
!
    t1 = t2
    u1(0:np) = u2(0:np)
!
!  Copy into the output arrays.
!
    t(it) = t1
    u(it,0:np) = u1(0:np)

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
