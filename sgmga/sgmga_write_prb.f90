program main

!*****************************************************************************80
!
!! MAIN is the main program for SGMGA_WRITE_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_WRITE_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SGMA_WRITE function.'

  call sgmga_write_tests  ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_WRITE_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sgmga_write_tests ( )

!****************************************************************************80
!
!! SGMGA_WRITE_TESTS calls SGMGA_WRITE_TEST with various arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) TOL, a tolerance for point equality.
!    A value of sqrt ( eps ) is reasonable, and will allow the code to
!    consolidate points which are equal, or very nearly so.  A value of
!    -1.0, on the other hand, will force the code to use every point, regardless
!    of duplication.
!
  implicit none

  integer ( kind = 4 )  dim
  integer ( kind = 4 )  dim_num
  character ( len = 255 ) file_name
  integer ( kind = 4 ), allocatable :: growth(:)
  real ( kind = 8 ), allocatable :: importance(:)
  integer ( kind = 4 )  level_max
  integer ( kind = 4 )  level_max_max
  integer ( kind = 4 )  level_max_min
  real ( kind = 8 ), allocatable :: level_weight(:)
  integer ( kind = 4 ), allocatable :: np(:)
  integer ( kind = 4 )  np_sum
  integer ( kind = 4 ), allocatable :: order_1d(:)
  integer ( kind = 4 )  order_nd
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 )  r8_epsilon
  integer ( kind = 4 ), allocatable :: rule(:)
  real ( kind = 8 )  tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_WRITE_TESTS'
  write ( *, '(a)' ) '  Call SGMGA_WRITE_TEST with various arguments.'
!
!  Set the point equality tolerance.
!
  tol = sqrt ( r8_epsilon ( ) )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  All tests will use a point equality tolerance of ', tol

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = 1.0D+00
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l2_ccxcc_iso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l2_ccxcc_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 3
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = 1.0D+00
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d3_l2_ccxccxcc_iso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 3
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d3_l2_ccxccxcc_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 3
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 3 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l3_ccxgp_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 4 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l2_ccxgl_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 7 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l2_ccxlg_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 8 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 1 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 1.5D+00
  file_name = 'sgmga_d2_l2_ccxglg_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 2, 9 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 2 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 0.5D+00
  p(2) = 1.5D+00
  file_name = 'sgmga_d2_l2_f2xgj_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 6, 10 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 3, 4 /)
  allocate ( np(1:dim_num) )
  np = (/ 1, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 1.0D+00
  file_name = 'sgmga_d2_l2_gghxhgk_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  LEVEL_MAX 1
!
  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 1
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l1_ccxcc_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  LEVEL_MAX 2 already done
!
!  LEVEL_MAX 3
!
  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 3
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l3_ccxcc_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  LEVEL_MAX 4
!
  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 4
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l4_ccxcc_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  LEVEL_MAX 5
!
  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 5
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l5_ccxcc_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  Dimension 3
!
  dim_num = 3
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 4, 5 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3, 3 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d3_l2_ccxglxgh_aniso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  Rule 3, LEVEL_MAX 4
!
  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = 1.0D+00
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 4
  allocate ( rule(1:dim_num) )
  rule = (/ 3, 3 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l4_gpxgp_iso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  Rule 3, Growth 4, LEVEL_MAX 4
!
  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = 1.0D+00
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 4
  allocate ( rule(1:dim_num) )
  rule = (/ 3, 3 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 4, 4 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l4_gpsexgpse_iso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  Rule 3, Growth 5, LEVEL_MAX 4
!
  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = 1.0D+00
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 4
  allocate ( rule(1:dim_num) )
  rule = (/ 3, 3 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 5, 5 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  file_name = 'sgmga_d2_l4_gpmexgpme_iso'
  call sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
    np, p, tol, file_name )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  return
end
subroutine sgmga_write_test ( dim_num, level_weight, level_max, rule, growth, &
  np, p, tol, file_name )

!****************************************************************************80
!
!! SGMGA_WRITE_TEST tests SGMGA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weight
!    for each dimension.
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
!    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
!    11, "UO",  User supplied Open, presumably Non Nested.
!    12, "UC",  User supplied Closed, presumably Non Nested.
!
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind = 4 ) NP(DIM_NUM), the number of parameters used by
!    each rule.
!
!    Input, real ( kind = 8 ) P(*), the parameters needed by each rule.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for point equality.
!
!    Input, character ( len = * ) FILE_NAME, the main name of the
!    output files.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * )  file_name
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) level_max
  real ( kind = 8 ) level_weight(dim_num)
  integer ( kind = 4 ) np(dim_num)
  real ( kind = 8 ) p(*)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num
  integer ( kind = 4 ) rule(dim_num)
  integer ( kind = 4 ), allocatable :: sparse_index(:,:)
  integer ( kind = 4 ), allocatable :: sparse_order(:,:)
  real ( kind = 8 ), allocatable :: sparse_point(:,:)
  integer ( kind = 4 ), allocatable :: sparse_unique_index(:)
  real ( kind = 8 ), allocatable :: sparse_weight(:)
  real ( kind = 8 ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_WRITE_TEST'
  write ( *, '(a)' ) '  SGMGA_WRITE writes a sparse grid rule to files.'
!
!  Compute necessary data.
!
  call sgmga_size_total ( dim_num, level_weight, level_max, rule, growth, &
    point_total_num )

  call sgmga_size ( dim_num, level_weight, level_max, rule, growth, np, &
    p, tol, point_num )

  allocate ( sparse_unique_index(1:point_total_num) )

  call sgmga_unique_index ( dim_num, level_weight, level_max, &
    rule, growth, np, p, tol, point_num, point_total_num, sparse_unique_index )

  allocate ( sparse_order(1:dim_num,1:point_num) )
  allocate ( sparse_index(1:dim_num,1:point_num) )

  call sgmga_index ( dim_num, level_weight, level_max, rule, growth, point_num, &
    point_total_num, sparse_unique_index, sparse_order, sparse_index )
!
!  Compute points and weights.
!
  allocate ( sparse_point(1:dim_num,1:point_num) )

  call sgmga_point ( dim_num, level_weight, level_max, rule, growth, &
    np, p, point_num, sparse_order, sparse_index, sparse_point )

  allocate ( sparse_weight(1:point_num) )

  call sgmga_weight ( dim_num, level_weight, level_max, rule, growth, &
    np, p, point_num, point_total_num, sparse_unique_index, &
    sparse_weight )
!
!  Write points and weights to files.
!
  call sgmga_write ( dim_num, level_weight, rule, growth, np, p, &
    point_num, sparse_weight, sparse_point, file_name )

  deallocate ( sparse_index )
  deallocate ( sparse_order )
  deallocate ( sparse_point )
  deallocate ( sparse_unique_index )
  deallocate ( sparse_weight )

  return
end
