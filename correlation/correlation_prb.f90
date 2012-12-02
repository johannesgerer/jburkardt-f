program main

!*****************************************************************************80
!
!! CORRELATION_PRB tests the CORRELATION library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CORRELATION_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the CORRELATION library.'

  call correlation_test01 ( )
  call correlation_test02 ( )
  call correlation_test03 ( )
  call correlation_test04 ( )
  call correlation_test05 ( )
  call correlation_test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CORRELATION_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine correlation_test01 ( )

!*****************************************************************************80
!
!! CORRELATION_TEST01 plots the correlation functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: c(:)
  real ( kind = 8 ) e
  integer ( kind = 4 ) n
  real ( kind = 8 ) nu
  real ( kind = 8 ), allocatable :: rho(:)
  real ( kind = 8 ) rho0

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CORRELATION_TEST01'
  write ( *, '(a)' ) '  Make plots of correlation functions.'
  write ( *, '(a)' ) ' '

  n = 101
  allocate ( rho(1:n) )
  allocate ( c(1:n) )
!
!  besselj
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -8.0D+00, 8.0D+00, rho )
  call correlation_besselj ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'besselj', 'Bessel J correlation' )
!
!  besselk
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -4.0D+00, 4.0D+00, rho )
  call correlation_besselk ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'besselk', 'Bessel K correlation' )
!
!  circular
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_circular ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'circular', 'Circular correlation' )
!
!  constant
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_constant ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'constant', 'Constant correlation' )
!
!  cubic
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_cubic ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'cubic', 'Cubic correlation' )
!
!  damped_cosine
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -6.0D+00, 6.0D+00, rho )
  call correlation_damped_cosine ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'damped_cosine', 'Damped cosine correlation' )
!
!  damped_sine
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -12.0D+00, 12.0D+00, rho )
  call correlation_damped_sine ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'damped_sine', 'Damped sine correlation' )
!
!  exponential
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_exponential ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'exponential', 'Exponential correlation' )
!
!  gaussian
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_gaussian ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'gaussian', 'Gaussian correlation' )
!
!  hole
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -6.0D+00, 6.0D+00, rho )
  call correlation_hole ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'hole', 'Hole correlation' )
!
!  linear
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_linear ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'linear', 'Linear correlation' )
!
!  matern, nu = 2.5
!
  rho0 = 1.0D+00
  nu = 2.5D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
! call correlation_matern ( n, rho, rho0, nu, c )
  call correlation_matern ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'matern', 'Matern correlation (NU = 2.5)' )
!
!  pentaspherical
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_pentaspherical ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'pentaspherical', 'Pentaspherical correlation' )
!
!  power, e = 2.0
!
  rho0 = 1.0D+00
  e = 2.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
! call correlation_power ( n, rho, rho0, e, c )
  call correlation_power ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'power', 'Power correlation' )
!
!  rational_quadratic
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -4.0D+00, 4.0D+00, rho )
  call correlation_rational_quadratic ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'rational_quadratic', 'Rational quadratic correlation' )
!
!  spherical
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_spherical ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'spherical', 'Spherical correlation' )
!
!  white_noise
!
  rho0 = 1.0D+00
  call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
  call correlation_white_noise ( n, rho, rho0, c )
  call correlation_plot ( n, rho, c, 'white_noise', 'White noise correlation' )

  deallocate ( c )
  deallocate ( rho )

  return
end
subroutine correlation_test02 ( )

!*****************************************************************************80
!
!! CORRELATION_TEST02 plots sample paths with SAMPLE_PATHS_CHOLESKY/EIGEN.
!
!  Discussion:
!
!    Most paths will be blue, but make the LAST one red so that there will
!    always be one distinguished path that is easy to follow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  external correlation_besselj
  external correlation_besselk
  external correlation_circular
  external correlation_constant
  external correlation_cubic
  external correlation_damped_cosine
  external correlation_damped_sine
  external correlation_exponential
  external correlation_gaussian
  external correlation_hole
  external correlation_linear
  external correlation_matern
  external correlation_pentaspherical
  external correlation_power
  external correlation_rational_quadratic
  external correlation_spherical
  external correlation_white_noise
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 8 ), allocatable :: rho(:)
  real ( kind = 8 ), parameter :: rho0 = 1.0D+00
  real ( kind = 8 ), parameter :: rhomax = 10.0D+00
  real ( kind = 8 ), parameter :: rhomin = 0.0D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CORRELATION_TEST02'
  write ( *, '(a)' ) '  SAMPLE_PATHS_CHOLESKY generates sample paths from the'
  write ( *, '(a)' ) '  correlation matrix, factored using the Cholesky factor.'
  write ( *, '(a)' ) '  It requires that the correlation matrix is nonnegative definite.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SAMPLE_PATHS_EIGEN generates sample paths from the'
  write ( *, '(a)' ) '  correlation matrix, factored using the eigen factorization.'
  write ( *, '(a)' ) '  If the correlation matrix is not nonnegative definite,'
  write ( *, '(a)' ) '  we simply suppress negative eigenvalues.'
  write ( *, '(a)' ) ''

  n = 101
  n2 = 3
  allocate ( rho(1:n) )
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( x(n,n2) )
!
!  besselj
!  Use EIGEN, because CHOLESKY fails.
!
  seed = 123456789
  call sample_paths_eigen ( n, n2, rhomax, rho0, correlation_besselj, seed, x )
  call paths_plot ( n, n2, rho, x, 'besselj', 'Bessel J correlation' )
!
!  besselk
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_besselk, seed, x )
  call paths_plot ( n, n2, rho, x, 'besselk', 'Bessel K correlation' )
!
!  circular
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_circular, seed, x )
  call paths_plot ( n, n2, rho, x, 'circular', 'Circular correlation' )
!
!  constant
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_constant, seed, x )
  call paths_plot ( n, n2, rho, x, 'constant', 'Constant correlation' )
!
!  cubic
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_cubic, seed, x )
  call paths_plot ( n, n2, rho, x, 'cubic', 'Cubic correlation' )
!
!  damped_cosine
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_damped_cosine, seed, x )
  call paths_plot ( n, n2, rho, x, 'damped_cosine', 'Damped cosine correlation' )
!
!  damped_sine
!  Use EIGEN, because CHOLESKY fails.
!
  seed = 123456789
  call sample_paths_eigen ( n, n2, rhomax, rho0, correlation_damped_sine, seed, x )
  call paths_plot ( n, n2, rho, x, 'damped_sine', 'Damped sine correlation' )
!
!  exponential
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_exponential, seed, x )
  call paths_plot ( n, n2, rho, x, 'exponential', 'Exponential correlation' )
!
!  gaussian
!  Use EIGEN, because CHOLESKY fails.
!
  seed = 123456789
  call sample_paths_eigen ( n, n2, rhomax, rho0, correlation_gaussian, seed, x )
  call paths_plot ( n, n2, rho, x, 'gaussian', 'Gaussian correlation' )
!
!  hole
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_hole, seed, x )
  call paths_plot ( n, n2, rho, x, 'hole', 'Hole correlation' )
!
!  linear
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_linear, seed, x )
  call paths_plot ( n, n2, rho, x, 'linear', 'Linear correlation' )
!
!  matern ( nu = 2.5 )
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_matern, seed, x )
  call paths_plot ( n, n2, rho, x, 'matern', 'Matern correlation (nu=2.5)' )
!
!  pentaspherical
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_pentaspherical, seed, x )
  call paths_plot ( n, n2, rho, x, 'pentaspherical', 'Pentaspherical correlation' )
!
!  power ( e = 2.0 )
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_power, seed, x )
  call paths_plot ( n, n2, rho, x, 'power', 'Power correlation (e=2.0)' )
!
!  rational_quadratic
!  Use EIGEN, because CHOLESKY fails.
!
  seed = 123456789
  call sample_paths_eigen ( n, n2, rhomax, rho0, correlation_rational_quadratic, seed, x )
  call paths_plot ( n, n2, rho, x, 'rational_quadratic', 'Rational quadratic correlation' )
!
!  spherical
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_spherical, seed, x )
  call paths_plot ( n, n2, rho, x, 'spherical', 'Spherical correlation' )
!
!  white_noise
!
  seed = 123456789
  call sample_paths_cholesky ( n, n2, rhomax, rho0, correlation_white_noise, seed, x )
  call paths_plot ( n, n2, rho, x, 'white_noise', 'White noise correlation' )

  return
end
subroutine correlation_test03 ( )

!*****************************************************************************80
!
!! CORRELATION_TEST03 plots a correlation function for several values of RH00.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: c(:,:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 8 ), allocatable :: rho(:)
  real ( kind = 8 ), allocatable :: rho0(:)
  real ( kind = 8 ) rhomax
  real ( kind = 8 ) rhomin

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CORRELATION_TEST03'
  write ( *, '(a)' ) '  Make plots of correlation functions for'
  write ( *, '(a)' ) '  a range of correlation lengths.'
  write ( *, '(a)' ) ''
!
!  besselj
!
  n = 101
  n2 = 5
  allocate ( rho0(1:n2) )
  rho0 = (/ 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -8.0D+00
  rhomax = +8.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_besselj ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'besselj', 'Bessel J correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  besselk
!
  n = 101
  n2 = 5
  allocate ( rho0(1:n2) )
  rho0 = (/ 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -4.0D+00
  rhomax = +4.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_besselk ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'besselk', 'Bessel K correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  circular
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -4.0D+00
  rhomax = +4.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_circular ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'circular', 'Circular correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  constant
!  1 plot is enough
!
  n = 101
  n2 = 1
  allocate ( rho0(1:n2) )
  rho0 = (/ 1.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_constant ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'constant', 'Constant correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  cubic
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -8.0D+00
  rhomax = +8.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_cubic ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'cubic', 'Cubic correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  damped_cosine
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -6.0D+00
  rhomax = +6.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_damped_cosine ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'damped_cosine', 'Damped cosine correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  damped_sine
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -8.0D+00
  rhomax = +8.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_damped_sine ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'damped_sine', 'Damped sine correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  exponential
!
  n = 101
  n2 = 7
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.25, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_exponential ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'exponential', 'Exponential correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  gaussian
!
  n = 101
  n2 = 7
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.25, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_gaussian ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'gaussian', 'Gaussian correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  hole
!
  n = 101
  n2 = 7
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.25, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_hole ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'hole', 'Hole correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  linear
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_linear ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'linear', 'Linear correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  matern, nu = 2.5
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_matern ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'matern', 'Matern correlation (NU = 2.5)' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  pentaspherical
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_pentaspherical ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'pentaspherical', 'Pentaspherical correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  power, e = 2.0
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_power ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'power', 'Power correlation (E = 2.0)' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  rational_quadratic
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -4.0D+00
  rhomax = +4.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_rational_quadratic ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'rational_quadratic', 'Rational quadratic correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  spherical
!
  n = 101
  n2 = 6
  allocate ( rho0(1:n2) )
  rho0 = (/ 0.5, 1.0, 1.5, 2.0, 4.0, 8.0 /)
  allocate ( rho(1:n) )
  rhomin = -8.0D+00
  rhomax = +8.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_spherical ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'spherical', 'Spherical correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )
!
!  white_noise
!  1 plot is enough
!
  n = 101
  n2 = 1
  allocate ( rho0(1:n2) )
  rho0 = (/ 1.0 /)
  allocate ( rho(1:n) )
  rhomin = -2.0D+00
  rhomax = +2.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  allocate ( c(n,n2) )
  do j = 1, n2
    call correlation_white_noise ( n, rho, rho0(j), c(1:n,j) )
  end do
  call correlation_plots ( n, n2, rho, rho0, c, 'white_noise', 'White noise correlation' )
  deallocate ( c )
  deallocate ( rho )
  deallocate ( rho0 )

  return
end
subroutine correlation_test04 ( )

!*****************************************************************************80
!
!! CORRELATION_TEST04 converts between covariance and correlation matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: c(:,:)
  real ( kind = 8 ), allocatable :: k(:,:)
  real ( kind = 8 ), allocatable :: k2(:,:)
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: sigma(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CORRELATION_TEST04'
  write ( *, '(a)' ) '  Convert between a correlation and a covariance matrix.'

  n = 5
  allocate ( k(1:n,1:n) )
  call minij ( n, n, k )

  call r8mat_print ( n, n, k, '  Covariance matrix K:' )

  allocate ( c(1:n,1:n) )
  allocate ( sigma(1:n) )

  call covariance_to_correlation ( n, k, c, sigma )

  call r8mat_print ( n, n, c, '  Correlation matrix C:' )

  call r8vec_print ( n, sigma, '  Variances:' )

  allocate ( k2(1:n,1:n) )

  call correlation_to_covariance ( n, c, sigma, k2 )

  call r8mat_print ( n, n, k2, '  Recovered covariance matrix K2:' )

  deallocate ( c )
  deallocate ( k )
  deallocate ( k2 )
  deallocate ( sigma )

  return
end
subroutine correlation_test05 ( )

!*****************************************************************************80
!
!! CORRELATION_TEST05 calls CORRELATION_BROWNIAN_DISPLAY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CORRELATION_TEST05'
  write ( *, '(a)' ) '  CORRELATION_BROWNIAN_DISPLAY displays 4 slices of'
  write ( *, '(a)' ) '  the Brownian correlation function.'

  call correlation_brownian_display ( )

  return
end
subroutine correlation_test06 ( )

!*****************************************************************************80
!
!! CORRELATION_TEST06 plots sample paths with SAMPLE_PATHS2_CHOLESKY/EIGEN/FFT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  external correlation_brownian
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 8 ), allocatable :: rho(:)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhomax
  real ( kind = 8 ) rhomin
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CORRELATION_TEST06'
  write ( *, '(a)' ) '  For non-stationary correlation functions:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  SAMPLE_PATHS2_CHOLESKY generates sample paths from the'
  write ( *, '(a)' ) '  correlation matrix, factored using the Cholesky factor.'
  write ( *, '(a)' ) '  It requires that the correlation matrix is nonnegative definite.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  SAMPLE_PATHS2_EIGEN generates sample paths from the'
  write ( *, '(a)' ) '  correlation matrix, factored using the eigen factorization.';
  write ( *, '(a)' ) '  If the correlation matrix is not nonnegative definite,'
  write ( *, '(a)' ) '  we simply suppress negative eigenvalues.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  SAMPLE_PATHS2_FFT generates sample paths from the'
  write ( *, '(a)' ) '  correlation matrix, factored using the FFT factorization';
  write ( *, '(a)' ) '  of the correlation matrix after embedding in a circulant.'
  write ( *, '(a)' ) ''
!
!  brownian
!
  n = 101
  n2 = 3
  rhomin = 0.0D+00
  rhomax = 10.0D+00
  rho0 = 1.0D+00
  seed = 123456789
  allocate ( x(1:n,1:n2) )
  call sample_paths2_cholesky ( n, n2, rhomin, rhomax, rho0, correlation_brownian, seed, x )
  allocate ( rho(1:n) )
  call r8vec_linspace ( n, rhomin, rhomax, rho )
  call paths_plot ( n, n2, rho, x, 'brownian', 'Brownian correlation' )

  deallocate ( rho )
  deallocate ( x )

  return
end

