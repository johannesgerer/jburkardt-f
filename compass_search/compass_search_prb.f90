program main

!*****************************************************************************80
!
!! COMPASS_SEARCH_TEST tests COMPASS_SEARCH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMPASS_SEARCH_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the COMPASS_SEARCH library.'

  call beale_test ( )
  call bohach1_test ( )
  call bohach2_test ( )
  call broyden_test ( )
  call extended_rosenbrock_test ( )
  call goldstein_price_test ( )
  call himmelblau_test ( )
  call local_test ( )
  call mckinnon_test ( )
  call powell_test ( )
  call rosenbrock_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMPASS_SEARCH_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine beale_test ( )

!*****************************************************************************80
!
!! BEALE_TEST works with the Beale function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), external :: beale
  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BEALE_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Beale function.'
  delta_tol = 0.00001D+00
  delta = 0.1D+00
  k_max = 20000

  x0 = (/ 1.0D+00, 1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', beale ( m, x0 )

  call compass_search ( beale, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i10)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Repeat with more difficult start.
!
  x0 = (/ 1.0D+00, 4.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', beale ( m, x0 )

  call compass_search ( beale, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate correct minimizer.
!
  x = (/ 3.0D+00, 0.5D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', beale ( m, x )

  return
end
subroutine bohach1_test ( )

!*****************************************************************************80
!
!! BOHACH1_TEST works with the Bohachevsky function #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), external :: bohach1
  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BOHACH1_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Bohachevsky function #1.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ 0.5D+00, 1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', bohach1 ( m, x0 )

  call compass_search ( bohach1, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate correct minimizer.
!
  x = (/ 0.0D+00, 0.0D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', bohach1 ( m, x )

  return
end
subroutine bohach2_test ( )

!*****************************************************************************80
!
!! BOHACH2_TEST works with the Bohachevsky function #2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), external :: bohach2
  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BOHACH2_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Bohachevsky function #2.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ 0.6D+00, 1.3D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', bohach2 ( m, x0 )

  call compass_search ( bohach2, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate correct minimizer.
!
  x = (/ 0.0D+00, 0.0D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', bohach2 ( m, x )

  return
end
subroutine broyden_test ( )

!*****************************************************************************80
!
!! BROYDEN_TEST works with the Broyden function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), external :: broyden
  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BROYDEN_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Broyden function.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ -0.9D+00, -1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', broyden ( m, x0 )

  call compass_search ( broyden, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate correct minimizer.
!
  x = (/ -0.511547141775014D+00, -0.398160951772036D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', broyden ( m, x )

  return
end
subroutine extended_rosenbrock_test ( )

!*****************************************************************************80
!
!! EXTENDED_ROSENBROCK_TEST works with the extended Rosenbrock function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4

  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ), external :: extended_rosenbrock
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXTENDED_ROSENBROCK_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the extended Rosenbrock function.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ - 1.2D+00, 1.0D+00,  -1.5D+00, 1.2D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', extended_rosenbrock ( m, x0 )

  call compass_search ( extended_rosenbrock, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate correct minimizer.
!
  x = (/ 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', extended_rosenbrock ( m, x )

  return
end
subroutine goldstein_price_test ( )

!*****************************************************************************80
!
!! GOLDSTEIN_PRICE_TEST works with the Goldstein-Price function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  real ( kind = 8 ), external :: goldstein_price
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GOLDSTEIN_PRICE_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Goldstein-Price function.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ -0.5D+00, 0.25D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', goldstein_price ( m, x0 )

  call compass_search ( goldstein_price, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Repeat with more difficult start.
!
  x0 = (/ -4.0D+00, 5.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', goldstein_price ( m, x0 )

  call compass_search ( goldstein_price, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate correct minimizer.
!
  x = (/ 0.0D+00, -1.0D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', goldstein_price ( m, x )

  return
end
subroutine himmelblau_test ( )

!*****************************************************************************80
!
!! HIMMELBLAU_TEST works with the Himmelblau function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  real ( kind = 8 ), external :: himmelblau
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HIMMELBLAU_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Himmelblau function.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ 1.0D+00, 1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', himmelblau ( m, x0 )

  call compass_search ( himmelblau, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k

  x0 = (/ -1.0D+00, 1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', himmelblau ( m, x0 )

  call compass_search ( himmelblau, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k

  x0 = (/ -1.0D+00, -1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', himmelblau ( m, x0 )

  call compass_search ( himmelblau, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k

  x0 = (/ 1.0D+00, -1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', himmelblau ( m, x0 )

  call compass_search ( himmelblau, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate Himmelblau minimizers.
!
  x = (/ 3.0D+00, 2.0D+00 /)
  call r8vec_print ( m, x, '  Correct Himmelblau minimizer X1*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', himmelblau ( m, x )

  x = (/ 3.58439D+00, -1.84813D+00 /)
  call r8vec_print ( m, x, '  Correct Himmelblau minimizer X2*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', himmelblau ( m, x )

  x = (/ -3.77934D+00, -3.28317D+00 /)
  call r8vec_print ( m, x, '  Correct Himmelblau minimizer X3*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', himmelblau ( m, x )

  x = (/ -2.80512D+00,  3.13134D+00 /)
  call r8vec_print ( m, x, '  Correct Himmelblau minimizer X4*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', himmelblau ( m, x )

  return
end
subroutine local_test ( )

!*****************************************************************************80
!
!! LOCAL_TEST works with the Local function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ), external :: local
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOCAL_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Local function.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ 1.0D+00, 1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', local ( m, x0 )

  call compass_search ( local, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate local minimizer.
!
  x = (/ 0.2858054412D+00, 0.2793263206D+00 /)
  call r8vec_print ( m, x, '  Correct local minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', local ( m, x )
!
!  Try for global minimizer.
!
  x0 = (/ -15.0D+00, -35.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', local ( m, x0 )

  call compass_search ( local, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate global minimizer.
!
  x = (/ -21.02667179D+00, -36.75997872D+00 /)
  call r8vec_print ( m, x, '  Correct global minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', local ( m, x )

  return
end
subroutine mckinnon_test ( )

!*****************************************************************************80
!
!! MCKINNON_TEST works with the McKinnon function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ), external :: mckinnon
  real ( kind = 8 ) phi
  real ( kind = 8 ) tau
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MCKINNON_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the McKinnon function.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000
!
!  Test 1
!
  a = ( 1.0D+00 + sqrt ( 33.0D+00 ) ) / 8.0D+00
  b = ( 1.0D+00 - sqrt ( 33.0D+00 ) ) / 8.0D+00

  phi = 10.0D+00
  tau = 1.0D+00
  theta = 15.0D+00

  call mckinnon_parameters ( 'set', phi, tau, theta )

  x0 = (/ a, b /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) '  PHI = ', phi, ' TAU = ', tau, ' THETA = ', theta
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', mckinnon ( m, x0 )

  call compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k

  x = (/ 0.0D+00, -0.5D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', mckinnon ( m, x )
!
!  Test 2
!
  a = ( 1.0D+00 + sqrt ( 33.0D+00 ) ) / 8.0D+00
  b = ( 1.0D+00 - sqrt ( 33.0D+00 ) ) / 8.0D+00

  phi = 60.0D+00
  tau = 2.0D+00
  theta = 6.0D+00

  call mckinnon_parameters ( 'set', phi, tau, theta )

  x0 = (/ a, b /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) '  PHI = ', phi, ' TAU = ', tau, ' THETA = ', theta
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', mckinnon ( m, x0 )

  call compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k

  x = (/ 0.0D+00, -0.5D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', mckinnon ( m, x )
!
!  Test 3
!
  a = ( 1.0D+00 + sqrt ( 33.0D+00 ) ) / 8.0D+00
  b = ( 1.0D+00 - sqrt ( 33.0D+00 ) ) / 8.0D+00

  phi = 4000.0D+00
  tau = 3.0D+00
  theta = 6.0D+00

  call mckinnon_parameters ( 'set', phi, tau, theta )

  x0 = (/ a, b /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) '  PHI = ', phi, ' TAU = ', tau, ' THETA = ', theta
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', mckinnon ( m, x0 )

  call compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k

  x = (/ 0.0D+00, -0.5D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', mckinnon ( m, x )

  return
end
subroutine powell_test ( )

!*****************************************************************************80
!
!! POWELL_TEST works with the Powell function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4

  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ), external :: powell
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POWELL_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Powell function.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ 3.0D+00, -1.0D+00, 0.0D+00, 1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', powell ( m, x0 )

  call compass_search ( powell, m, x0, delta_tol, delta, k_max, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate correct minimizer.
!
  x = (/ 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', powell ( m, x )

  return
end
subroutine rosenbrock_test ( )

!*****************************************************************************80
!
!! ROSENBROCK_TEST works with the Rosenbrock function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ), external :: rosenbrock
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ROSENBROCK_TEST:'
  write ( *, '(a)' ) '  Test COMPASS_SEARCH with the Rosenbrock function.'
  delta_tol = 0.00001D+00
  delta = 0.3D+00
  k_max = 20000

  x0 = (/ - 1.2D+00, 1.0D+00 /)
  call r8vec_print ( m, x0, '  Initial point X0:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X0) = ', rosenbrock ( m, x0 )

  call compass_search ( rosenbrock, m, x0, delta_tol, delta, k_max, x, fx, k, x, fx, k )
  call r8vec_print ( m, x, '  Estimated minimizer X1:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,i12)' ) '  F(X1) = ', fx, ' number of steps = ', k
!
!  Demonstrate correct minimizer.
!
  x = (/ 1.0D+00, 1.0D+00 /)
  call r8vec_print ( m, x, '  Correct minimizer X*:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', rosenbrock ( m, x )

  return
end
function beale ( m, x )

!*****************************************************************************80
!
!! BEALE computes the Beale function.
!
!  Discussion:
!
!    This function has a global minimizer:
!
!      X = ( 3.0, 0.5 )
!
!    for which
!
!      F(X) = 0.
!
!    For a relatively easy computation, start the search at
!
!      X = ( 1.0, 1.0 )
!
!    A harder computation starts at
!
!      X = ( 1.0, 4.0 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Evelyn Beale,
!    On an Iterative Method for Finding a Local Minimum of a Function
!    of More than One Variable,
!    Technical Report 25, 
!    Statistical Techniques Research Group,
!    Princeton University, 1958.
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) BEALE, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) beale
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) x(m)

  f1 = 1.5D+00   - x(1) * ( 1.0D+00 - x(2)    )
  f2 = 2.25D+00  - x(1) * ( 1.0D+00 - x(2)**2 )
  f3 = 2.625D+00 - x(1) * ( 1.0D+00 - x(2)**3 )

  f = f1**2 + f2**2 + f3**2
 
  beale = f

  return
end
function bohach1 ( m, x )

!*****************************************************************************80
!
!! BOHACH1 evaluates the Bohachevsky function #1.
!
!  Discussion:
!
!    The minimizer is
!
!      X* = [ 0.0, 0.0 ]
!      F(X*) = 0.0
!
!    Suggested starting point:
!
!      X = [ 0.5, 1.0 ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) BOHACH1, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) bohach1
  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m)

  f =           x(1) * x(1) - 0.3D+00 * cos ( 3.0D+00 * pi * x(1) ) &
    + 2.0D+00 * x(2) * x(2) - 0.4D+00 * cos ( 4.0D+00 * pi * x(2) ) &
    + 0.7D+00

  bohach1 = f

  return
end
function bohach2 ( m, x )

!*****************************************************************************80
!
!! BOHACH2 evaluates the Bohachevsky function #2.
!
!  Discussion:
!
!    The minimizer is
!
!      X* = [ 0.0, 0.0 ]
!      F(X*) = 0.0
!
!    Suggested starting point:
!
!      X = [ 0.6, 1.3 ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) BOHACH2, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) bohach2
  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m)

  f =           x(1) * x(1) & 
    + 2.0D+00 * x(2) * x(2) &
    - 0.3D+00 * cos ( 3.0D+00 * pi * x(1) ) * cos ( 4.0D+00 * pi * x(2) ) &
    + 0.3D+00

  bohach2 = f

  return
end
function broyden ( m, x )

!*****************************************************************************80
!
!! BROYDEN computes the two dimensional modified Broyden function.
!
!  Discussion:
!
!    This function has a global minimizer:
!
!      X = ( -0.511547141775014, -0.398160951772036 )
!
!    for which
!
!      F(X) = 1.44E-04
!
!    Start the search at
!
!      X = ( -0.9, -1.0 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Broyden,
!    A class of methods for solving nonlinear simultaneous equations,
!    Mathematics of Computation,
!    Volume 19, 1965, pages 577-593.
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    Testing unconstrained optimization software,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 1, March 1981, pages 17-41. 
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) BROYDEN, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) broyden
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) p
  real ( kind = 8 ) x(m)

  f1 = abs ( ( 3.0D+00 -           x(1) ) * x(1) - 2.0D+00 * x(2) + 1.0D+00 )
  f2 = abs ( ( 3.0D+00 - 2.0D+00 * x(2) ) * x(2) -           x(1) + 1.0D+00 )

  p = 3.0D+00 / 7.0D+00

  f = f1**p + f2**p
 
  broyden = f

  return
end
function extended_rosenbrock ( m, x )

!*****************************************************************************80
!
!! EXTENDED_ROSENBROCK computes the extended Rosenbrock function.
!
!  Discussion:
!
!    The number of dimensions is arbitrary, except that it must be even.
!
!    There is a global minimum at X* = (1,1,&), F(X*) = 0.
!
!    The contours are sharply twisted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Rosenbrock,
!    An Automatic Method for Finding the Greatest or Least Value of a Function,
!    Computer Journal,
!    Volume 3, 1960, pages 175-184.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) EXTENDED_ROSENBROCK, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) extended_rosenbrock
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(m)

  f = 0.0D+00

  do i = 1, m - 1, 2
    f1 = 1.0D+00 - x(i)
    f2 = 10.0D+00 * ( x(i+1) - x(i)**2 )
    f = f + f1**2 + f2**2
  end do

  extended_rosenbrock = f

  return
end
function goldstein_price ( m, x )

!*****************************************************************************80
!
!! GOLDSTEIN_PRICE evaluates the Goldstein-Price polynomial.
!
!  Discussion:
!
!    The minimizer is
!
!      X* = [ 0.0, -1.0 ]
!      F(X*) = 3.0
!
!    Suggested starting point:
!
!      X = [ -0.5, 0.25 ] (easy convergence)
!      X = [ -4.0, 5.00 ] (harder convergence)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Zbigniew Michalewicz,
!    Genetic Algorithms + Data Structures = Evolution Programs,
!    Third Edition,
!    Springer Verlag, 1996,
!    ISBN: 3-540-60676-9,
!    LC: QA76.618.M53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) GOLDSTEIN_PRICE, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) goldstein_price
  real ( kind = 8 ) x(m)

  a = x(1) + x(2) + 1.0D+00

  b = 19.0D+00 - 14.0D+00 * x(1) + 3.0D+00 * x(1) * x(1) - 14.0D+00 * x(2) &
    + 6.0D+00 * x(1) * x(2) + 3.0D+00 * x(2) * x(2)

  c = 2.0D+00 * x(1) - 3.0D+00 * x(2)

  d = 18.0D+00 - 32.0D+00 * x(1) + 12.0D+00 * x(1) * x(1) + 48.0D+00 * x(2) &
    - 36.0D+00 * x(1) * x(2) + 27.0D+00 * x(2) * x(2)

  f = ( 1.0D+00 + a * a * b ) * ( 30.0D+00 + c * c * d )

  goldstein_price = f

  return
end
function himmelblau ( m, x )

!*****************************************************************************80
!
!! HIMMELBLAU computes the Himmelblau function.
!
!  Discussion:
!
!    This function has 4 global minima:
!
!      X* = (  3,        2       ), F(X*) = 0.
!      X* = (  3.58439, -1.84813 ), F(X*) = 0.
!      X* = ( -3.77934, -3.28317 ), F(X*) = 0.
!      X* = ( -2.80512,  3.13134 ), F(X*) = 0.
!
!    Suggested starting points are
!
!      (+1,+1),
!      (-1,+1),
!      (-1,-1),
!      (+1,-1),
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Himmelblau,
!    Applied Nonlinear Programming,
!    McGraw Hill, 1972,
!    ISBN13: 978-0070289215,
!    LC: T57.8.H55.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) HIMMELBLAU, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) f
  real ( kind = 8 ) himmelblau
  real ( kind = 8 ) x(m)

  f = ( x(1)**2 + x(2) - 11.0D+00 )**2 + ( x(1) + x(2)**2 - 7.0D+00 )**2

  himmelblau = f
 
  return
end
function local ( m, x )

!*****************************************************************************80
!
!! LOCAL computes the local function.
!
!  Discussion:
!
!    This function has a local minimizer:
!
!      X* = ( 0.2858054412&, 0.2793263206&), F(X*) = 5.9225&
!
!    and a global minimizer:
!
!      X* = ( -21.02667179&, -36.75997872&), F(X*) = 0.
!
!    Suggested starting point for local minimizer:
!
!      X = ( 1, 1 ), F(X) = 3.33 * 10^6.
!
!    Suggested starting point for global minimizer:
!
!      X = ( -15, -35), F(X) = 1.49 * 10^8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Himmelblau,
!    Applied Nonlinear Programming,
!    McGraw Hill, 1972,
!    ISBN13: 978-0070289215,
!    LC: T57.8.H55.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) LOCAL, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) f
  real ( kind = 8 ) local
  real ( kind = 8 ) x(m)

  f = ( x(1)**2 + 12.0D+00 * x(2) - 1.0D+00 )**2 &
    + ( 49.0D+00 * x(1)**2 + 49.0D+00 * x(2)**2 + 84.0D+00 * x(1) &
    + 2324.0D+00 * x(2) - 681.0D+00 )**2
 
  local = f

  return
end
function mckinnon ( m, x )

!*****************************************************************************80
!
!! MCKINNON computes the McKinnon function.
!
!  Discussion:
!
!    This function has a global minimizer:
!
!      X* = ( 0.0, -0.5 ), F(X*) = -0.25.
!
!    There are three parameters, TAU, THETA and PHI.
!
!    1 < TAU, then F is strictly convex.
!             and F has continuous first derivatives.
!    2 < TAU, then F has continuous second derivatives.
!    3 < TAU, then F has continuous third derivatives.
!
!    However, this function can cause the Nelder-Mead optimization
!    algorithm to "converge" to a point which is not the minimizer
!    of the function F.
!
!    Sample parameter values which cause problems for Nelder-Mead 
!    include:
!
!      PHI = 10,  TAU = 1, THETA = 15
!      PHI = 60,  TAU = 2, THETA =  6
!      PHI = 400, TAU = 3, THETA =  6
!
!    To get the bad behavior, we also assume the initial simplex has the form
!
!      X1 = (0,0),
!      X2 = (1,1),
!      X3 = (A,B), 
!
!    where 
!
!      A = (1+sqrt(33))/8 =  0.84307&
!      B = (1-sqrt(33))/8 = -0.59307&
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ken McKinnon,
!    Convergence of the Nelder-Mead simplex method to a nonstationary point,
!    SIAM Journal on Optimization,
!    Volume 9, Number 1, 1998, pages 148-158.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) MCKINNON, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) f
  real ( kind = 8 ) mckinnon
  real ( kind = 8 ) phi
  real ( kind = 8 ) tau
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(m)

  call mckinnon_parameters ( 'get', phi, tau, theta )

  if ( x(1) <= 0.0D+00 ) then
    f = theta * phi * abs ( x(1) )**tau + x(2) * ( 1.0D+00 + x(2) )
  else
    f = theta       *       x(1)**tau   + x(2) * ( 1.0D+00 + x(2) )
  end if

  mckinnon = f

  return
end
subroutine mckinnon_parameters ( action, phi, tau, theta )

!*****************************************************************************80
!
!! MCKINNON_PARAMETERS sets or gets McKinnon function parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.  
!    'S', set the parameters.
!    'G', get the parameters.
!
!    Input/output, real ( kind = 8 ) PHI, TAU, THETA, the parameter values.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ) phi
  real ( kind = 8 ), save :: phi_save = 60.0D+00
  real ( kind = 8 ) tau
  real ( kind = 8 ), save :: tau_save = 2.0D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ), save :: theta_save = 6.0D+00

  if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then
    phi_save = phi
    tau_save = tau
    theta_save = theta
  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then
    phi = phi_save
    tau = tau_save
    theta = theta_save
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MCKINNON_PARAMETERS - Fatal error!'
    write ( *, '(a)' ) '  Unexpected value of ACTION.'
    stop
  end if

  return
end
function powell ( m, x )

!*****************************************************************************80
!
!! POWELL computes the Powell singular quartic function.
!
!  Discussion:
!
!    This function has a global minimizer:
!
!      X* = ( 0.0, 0.0, 0.0, 0.0 ), F(X*) = 0.
!
!    Start the search at
!
!      X = ( 3.0, -1.0, 0.0, 1.0 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Powell,
!    An Iterative Method for Finding Stationary Values of a Function
!    of Several Variables,
!    Computer Journal,
!    Volume 5, 1962, pages 147-151.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) POWELL, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ) powell
  real ( kind = 8 ) x(m)

  f1 = x(1) + 10.0D+00 * x(2)
  f2 =                                x(3) - x(4)
  f3 =               x(2) - 2.0D+00 * x(3)
  f4 = x(1)                                - x(4)

  f = f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4
 
  powell = f

  return
end
function rosenbrock ( m, x )

!*****************************************************************************80
!
!! ROSENBROCK computes the Rosenbrock function.
!
!  Discussion:
!
!    There is a global minimum at X* = (1,1), F(X*) = 0.
!
!    The starting point X = [ -1.2, 1.0 ] is suggested.
!
!    The contours are sharply twisted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Rosenbrock,
!    An Automatic Method for Finding the Greatest or Least Value of a Function,
!    Computer Journal,
!    Volume 3, 1960, pages 175-184.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X(M), the argument of the function.
!
!    Output, real ( kind = 8 ) ROSENBROCK, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) f
  real ( kind = 8 ) rosenbrock
  real ( kind = 8 ) x(m)

  f = ( 1.0D+00 - x(1) )**2 + 100.0D+00 * ( x(2) - x(1) * x(1) )**2

  rosenbrock = f

  return
end
