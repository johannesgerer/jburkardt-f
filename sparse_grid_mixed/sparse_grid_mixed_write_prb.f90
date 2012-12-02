program main

!****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_MIXED_WRITE_PRB.
!
!  Discussion:
!
!    SPARSE_GRID_MIXED_WRITE_PRB tests the SPARSE_GRID_MIXED_WRITE routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
  implicit none

  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) tol

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_WRITE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'

  tol = r8_epsilon ( )
  call sparse_grid_mixed_write_tests ( tol )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_WRITE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine sparse_grid_mixed_write_tests ( tol )

!****************************************************************************80
!
!! SPARSE_GRID_MIXED_WRITE_TESTS calls SPARSE_GRID_MIXED_WRITE_TEST with various arguments.
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
!    Input, real ( kind = 8 ) TOL, a tolerance for point equality.
!    A value of sqrt ( eps ) is reasonable, and will allow the code to
!    consolidate points which are equal, or very nearly so.  A value of
!    -1.0, on the other hand, will force the code to use every point, regardless
!    of duplication.
!
  implicit none

  real ( kind = 8  ), allocatable, dimension ( : ) :: alpha
  real ( kind = 8  ), allocatable, dimension ( : ) :: beta
  integer   ( kind = 4  ) dim_num
  character ( len = 255 ) file_name
  integer   ( kind = 4  ) level_max
  integer   ( kind = 4  ) level_max_max
  integer   ( kind = 4  ) level_max_min
  integer   ( kind = 4  ), allocatable, dimension ( : ) :: order_1d
  integer   ( kind = 4  ) order_nd
  integer   ( kind = 4  ), allocatable, dimension ( : ) :: rule
  real ( kind = 8  ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_WRITE_TESTS'
  write ( *, '(a)' ) '  Call SPARSE_GRID_MIXED_WRITE_TEST with various arguments.'
  write ( *, '(a,g14.6)' ) '  All tests will use a point equality tolerance of ', tol

  dim_num = 2
  level_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l2_ccxcc'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max = 3
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 3 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l3_ccxgp'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 4 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l2_ccxgl'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 7 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l2_ccxlg'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 8 /)
  alpha = (/ 0.0D+00, 1.5D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l2_ccxglg'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 2, 9 /)
  alpha = (/ 0.0D+00, 0.5D+00 /)
  beta = (/ 0.0D+00, 1.5D+00 /)
  file_name = 'sparse_grid_mixed_d2_l2_f2xgj'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 6, 4 /)
  alpha = (/ 2.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l2_gghxgl'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 3
  level_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 2, 5 /)
  alpha = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d3_l2_ccxf2xgh'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )
!
!  Dimension 2, Level 4, Rule 3
!
  dim_num = 2
  level_max = 4
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 3, 3 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l4_gpxgp'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )
!
!  Dimension 2, Level 4, Rule 13
!
  dim_num = 2
  level_max = 4
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 13, 13 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l4_gpsexgpse'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )
!
!  Dimension 2, Level 4, Rule 16
!
  dim_num = 2
  level_max = 4
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 16, 16 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l4_gpmexgpme'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )
!
!  Dimension 2, Level 4, Rule 17
!
  dim_num = 2
  level_max = 4
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 17, 17 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  file_name = 'sparse_grid_mixed_d2_l4_ccnxccn'
  call sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha, &
    beta, tol, file_name )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  return
end
subroutine sparse_grid_mixed_write_test ( dim_num, level_max, rule, &
  alpha, beta, tol, file_name )

!****************************************************************************80
!
!! SPARSE_GRID_MIXED_WRITE_TEST tests SPARSE_GRID_MIXED_WRITE.
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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level that defines the grid.
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
!    Generalized Gauss Hermite, Generalized Gauss Laguerre, and Gauss Jacobi rules.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for point equality.
!
!    Input, character ( len = * ) FILE_NAME, the main name of the
!    output files.
!
  implicit none

  integer   ( kind = 4 ) dim_num

  real ( kind = 8 ) alpha(dim_num)
  real ( kind = 8 ) beta(dim_num)
  character ( len = *  ) file_name
  integer   ( kind = 4 ) level_max
  integer   ( kind = 4 ) point_num
  integer   ( kind = 4 ) point_total_num
  integer   ( kind = 4 ) rule(dim_num)
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_index
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_point
  integer   ( kind = 4 ), allocatable, dimension(:) :: sparse_unique_index
  real ( kind = 8 ), allocatable, dimension ( : ) :: sparse_weight
  real ( kind = 8 ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_WRITE_TEST'
  write ( *, '(a)' ) '  SPARSE_GRID_MIXED_WRITE writes a sparse grid rule'
  write ( *, '(a)' ) '  to A, B, R, W, and X files.'
!
!  Compute necessary data.
!
  call sparse_grid_mixed_size_total ( dim_num, level_max, rule, &
    point_total_num )

  call sparse_grid_mixed_size ( dim_num, level_max, rule, alpha, beta, &
    tol, point_num )

  allocate ( sparse_unique_index(1:point_total_num) )

  call sparse_grid_mixed_unique_index ( dim_num, level_max, rule, alpha, &
    beta, tol, point_num, point_total_num, sparse_unique_index )

  allocate ( sparse_order(1:dim_num,1:point_num) )
  allocate ( sparse_index(1:dim_num,1:point_num) )

  call sparse_grid_mixed_index ( dim_num, level_max, rule, point_num, &
    point_total_num, sparse_unique_index, sparse_order, sparse_index )
!
!  Compute points and weights.
!
  allocate ( sparse_point(1:dim_num,1:point_num) )

  call sparse_grid_mixed_point ( dim_num, level_max, rule, alpha, beta, &
    point_num, sparse_order, sparse_index, sparse_point )

  allocate ( sparse_weight(1:point_num) )

  call sparse_grid_mixed_weight ( dim_num, level_max, rule, alpha, beta, &
    point_num, point_total_num, sparse_unique_index, sparse_weight )
!
!  Write points and weights to files.
!
  call sparse_grid_mixed_write ( dim_num, rule, alpha, beta, point_num, &
    sparse_weight, sparse_point, file_name )

  deallocate ( sparse_index )
  deallocate ( sparse_order )
  deallocate ( sparse_point )
  deallocate ( sparse_unique_index )
  deallocate ( sparse_weight )

  return
end
