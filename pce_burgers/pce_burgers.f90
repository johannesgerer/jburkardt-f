program main

!*****************************************************************************80
!
!! MAIN is the main program for PCE_BURGERS.
!
!  Discussion:
!
!    The time-dependent viscous Burgers equation to be solved is:
!
!      du/dt = - d ( u*(1/2-u)) /dx + nu d2u/dx2
!
!    with boundary conditions
!
!      u(-3.0) = 0.0, u(+3.0) = 1.0.
!
!    The viscosity nu is assumed to be an uncertain quantity with
!    normal distribution of known mean and variance.
!
!    A polynomial chaos expansion is to be used, with Hermite polynomial
!    basis functions h(i,x), 0 <= i <= n.
!
!    Because the first two Hermite polynomials are simply 1 and x, 
!    we have that 
!
!      nu = nu_mean * h(0,x) + nu_variance * h(1,x).
!
!    We replace the time derivative by an explicit Euler approximation,
!    so that the equation now describes the value of U(x,t+dt) in terms
!    of known data at time t.
!
!    Now assume that the solution U(x,t) can be approximated
!    by the truncated expansion:
!
!      U(x,t) = sum ( 0 <= i <= n ) c(i,t) * h(i,x)
!
!    In the equation, we replace U by its expansion, and then multiply
!    successively by each of the basis functions h(*,x) to get a set of
!    n+1 equations that can be used to determine the values of c(i,t+dt).
!
!    This process is repeated until the desired final time is reached.
!
!    At any time, the coefficients c(0,t) contain information definining
!    the expected value of u(x,t) at that time, while the higher order 
!    coefficients can be used to deterimine higher moments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2012
!
!  Author:
!
!    Original FORTRAN90 version by Gianluca Iaccarino.
!    This FORTRAN90 version is by John Burkardt.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) DT, the timestep.
!
!    Local, real ( kind = 8 ) DX, the spacing between grid points.
!
!    Local, integer ( kind = 4 ) N, the number of intervals in the spatial domain.
!
!    Local, real ( kind = 8 ) NUMEAN, the mean of viscosity.
!
!    Local, real ( kind = 8 ) NUVARIANCE, the variance of viscosity.
!
!    Local, integer ( kind = 4 ) P, the order of the PC expansion.
!
!    Local, real ( kind = 8 ) T, the current time.
!
!    Local, real ( kind = 8 ) TF, the final integration time.
!
!    Local, real ( kind = 8 ) U1(N+1,P+1), the PCE representation at the 
!    current time.
!
!    Local, real ( kind = 8 ) U2(N+1,P+1), the PCE representation for the next time.
!
!    Local, real ( kind = 8 ) X(N+1,1), the grid points.
!
  implicit none

  real ( kind = 8 ) conv
  real ( kind = 8 ) dp
  real ( kind = 8 ) dt
  real ( kind = 8 ) dx
  real ( kind = 8 ) he_double_product_integral
  real ( kind = 8 ) he_triple_product_integral
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nt
  real ( kind = 8 ) numean
  real ( kind = 8 ) nuvariance
  character ( len = 80 ) output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) p
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) term1
  real ( kind = 8 ) term2
  real ( kind = 8 ) tf
  real ( kind = 8 ) ti
  real ( kind = 8 ) tp
  real ( kind = 8 ), allocatable :: u1(:,:)
  real ( kind = 8 ), allocatable :: u2(:,:)
  real ( kind = 8 ), allocatable :: umean(:)
  real ( kind = 8 ), allocatable :: uvariance(:)
  real ( kind = 8 ) :: visc(2)
  real ( kind = 8 ), allocatable :: x(:)

  p = 5 
  n = 32
  nt = 2000
  ti = 0.0D+00
  tf = 2.0D+00
  dt = ( tf - ti ) / real ( nt, kind = 8 )
  numean = 0.25D+00
  nuvariance = 0.08D+00

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PCE_BURGERS:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Polynomial Chaos Solution'
  write ( *, '(a)' ) '  1D Burgers equation'
  write ( *, '(a)' ) '  Original version by Gianluca Iaccarino'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  PCE order       = ', p
  write ( *, '(a,i8)' ) '  Number of cells = ', n
  write ( *, '(a,g14.6)' ) '  Time step       = ', dt
  write ( *, '(a,g14.6)' ) '  Initial time    = ', ti
  write ( *, '(a,g14.6)' ) '  Final time      = ', tf
  write ( *, '(a,g14.6)' ) '  Viscosity Mean  = ', numean
  write ( *, '(a,g14.6)' ) '  Viscosity Var   = ', nuvariance
  write ( *, '(a)' ) ' '

  allocate ( u1(n+1,p+1) )
  allocate ( u2(n+1,p+1) )
  allocate ( x(n+1) )
!
!  Define some numerical parameters
!
  dx = 6.0D+00 / real ( n, kind = 8 )
  conv = dt / ( 2.0D+00 * dx )
!
!  The expansion for viscosity stops at the linear term.
!
  visc(1) = numean * dt / ( dx * dx )
  visc(2) = nuvariance * dt / ( dx * dx )
!
!  Define a uniform grid
!
  do i = 1, n + 1
    x(i) = ( real ( n - i + 1, kind = 8 ) * ( -3.0D+00 )   &
           + real (     i - 1, kind = 8 ) * ( +3.0D+00 ) ) &
           / real ( n,         kind = 8 )
  end do
!
!  Set the initial conditions
!
  u1(1:n+1,1:p+1) = 0.0D+00
  u1(1:n+1,1) = 0.5D+00 + x(1:n+1) / 6.0D+00
!
!  Write the current solution.
!
  output_filename = 'burgers.history.txt'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace' )

  write ( output_unit, '(a)' ) '----------'
  write ( output_unit, '(a,g14.6)' ) 'T = ', t1
  do i = 1, n + 1
    write ( output_unit, '(10(2x,g14.6))' ) x(i), u1(i,1:p+1)
  end do
!
!  Time integration
!
  t1 = ti

  do it = 1, nt

    t2 = ( real ( nt - it, kind = 8 ) * ti   &
         + real (      it, kind = 8 ) * tf ) &
         / real ( nt,      kind = 8 )
!
!  Boundary conditions.
!
    u2(1,1:p+1) = 0.0D+00
    u2(n+1,1) = 1.0D+00
    u2(n+1,2:p+1) = 0.0D+00

    do k = 1, p + 1

      dp = he_double_product_integral ( k - 1, k - 1 )

      do ix = 2, n
!
!  Viscous term.
!
        term1 = visc(1) * ( u1(ix+1,k) - 2.0D+00 * u1(ix,k) + u1(ix-1,k) )
        i = 2
        do j = 1, p + 1
          tp = he_triple_product_integral ( i - 1, j - 1, k - 1 )
          term1 = term1 + visc(i) * ( u1(ix+1,j) - 2.0D+00 * u1(ix,j) + u1(ix-1,j) ) * tp / dp
        end do
!
!  Convective term.
!
        term2 = - conv * 0.5D+00 * ( u1(ix+1,k) - u1(ix-1,k) )
        do j = 1, p + 1
          do i = 1, p + 1
            tp = he_triple_product_integral ( i - 1, j - 1, k - 1 )
            term2 = term2 + ( conv * u1(ix,i) * ( u1(ix+1,j) - u1(ix-1,j) ) * tp ) / dp
          end do
        end do

        u2(ix,k) = u1(ix,k) + term1 + term2

      end do
    end do

    t1 = t2
    u1(1:n+1,1:p+1) = u2(1:n+1,1:p+1)
!
!  Print solution every 100 time steps.
!
    if ( mod ( it, 100 ) == 0 ) then
      write ( output_unit, '(a)' ) '----------'
      write ( output_unit, '(a,g14.6)' ) 'T = ', t1
      do i = 1, n + 1
        write ( output_unit, '(20(2x,g14.6))' ) x(i), u1(i,1:p+1)
      end do
    end if

  end do

  close ( unit = output_unit )
  write ( *, '(a)' ) '  Time history in "' // trim ( output_filename ) // '".'
!
!  Compute the mean and variance.
!
  allocate ( umean(n+1) )
  allocate ( uvariance(n+1) )

  umean(1:n+1) = u1(1:n+1,1)

  do i = 1, n + 1
    uvariance(i) = 0.0D+00
    do j = 2, p + 1
      dp = he_double_product_integral ( j - 1, j - 1 )
      uvariance(i) = uvariance(i) + u1(i,j)**2 * dp
    end do
  end do
!
!  Write data about the solution at the final time.
!
  output_filename = 'burgers.moments.txt'
  open ( unit = output_unit, file = output_filename, status = 'replace' )
  write ( output_unit, '(a)' ) 'X E[U] Var[U]'
  do i = 1, n + 1
    write ( output_unit, '(20(g18.8,2x))' ) x(i), umean(i), uvariance(i)
  end do
  close ( unit = output_unit )
  write ( *, '(a)' ) '  Moments written to "' // trim ( output_filename ) // '".'

  output_filename = 'burgers.modes.txt'
  open ( unit = output_unit, file = output_filename, status = 'replace' )
  write ( output_unit, '(a)' ) 'X U_0 ... U_P'
  do i = 1, n + 1
    write ( output_unit, '(20(2x,g14.6))' ) x(i), u1(i,1:p+1)
  end do
  close ( unit = output_unit )
  write ( *, '(a)' ) '  Final modes written to "' // trim ( output_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PCE_BURGERS:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
