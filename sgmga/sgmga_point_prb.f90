program main

!*****************************************************************************80
!
!! MAIN is the main program for SGMGA_POINT_PRB.
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
  write ( *, '(a)' ) 'SGMGA_POINT_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SGMGA_POINT function.'

  call sgmga_point_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_POINT_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sgmga_point_tests ( )

!****************************************************************************80
!
!! SGMGA_POINT_TESTS calls SGMGA_POINT_TEST with various arguments.
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

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable :: growth(:)
  real ( kind = 8 ), allocatable :: importance(:)
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  real ( kind = 8 ), allocatable :: level_weight(:)
  integer ( kind = 4 ), allocatable :: np(:)
  integer ( kind = 4 ) np_sum
  integer ( kind = 4 ), allocatable :: order_1d(:)
  integer ( kind = 4 ) order_nd
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ), allocatable :: rule(:)
  real ( kind = 8 ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_POINT_TESTS'
  write ( *, '(a)' ) '  Call SGMGA_POINT_TEST with various arguments.'
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6, 6 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6, 6 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 3 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 4 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 7 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 8 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 1 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 1.5D+00
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 2, 9 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 2 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 0.5D+00
  p(2) = 1.5D+00
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 6, 10 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 3, 4 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 1, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 2.0D+00
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
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
  level_max_min = 0
  level_max_max = 2
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 4, 5 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 3, 3 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  Look at a case of interest to Mike.
!
  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim_num + 1 - dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max_min = 0
  level_max_max = 5
  allocate ( rule(1:dim_num) )
  rule = (/ 5, 5 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 3, 3 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )
!
!  Look at a case that includes a "0" importance dimension.
!
  dim_num = 3
  allocate ( importance(1:dim_num) )
  importance(1:dim_num) = (/ 1.0D+00, 0.0D+00, 1.0D+00 /)
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max_min = 0
  level_max_max = 3
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1, 1 /)
  allocate ( growth(1:dim_num) )
  growth = (/ 6, 6, 6 /)
  allocate ( np(1:dim_num) )
  np(1:dim_num) = (/ 0, 0, 0/)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_point_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max, rule, growth, np, p, tol )
  deallocate ( growth )
  deallocate ( importance )
  deallocate ( level_weight )
  deallocate ( np )
  deallocate ( p )
  deallocate ( rule )

  return
end
subroutine sgmga_point_test ( dim_num, importance, level_weight, &
  level_max_min, level_max_max, rule, growth, np, p, tol )

!****************************************************************************80
!
!! SGMGA_POINT_TEST tests SGMGA_POINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2010
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
!    Input, integer ( kind = 4 ) LEVEL_MAX_MIN, LEVEL_MAX_MAX, the minimum and
!    maximum values of LEVEL_MAX.
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
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the desired growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind = 4 ) NP(DIM_NUM), the number of parameters used by each rule.
!
!    Input, real ( kind = 8 ) P(*), the parameters needed by each rule.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for point equality.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) growth(dim_num)
  real ( kind = 8 ) importance(dim_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  real ( kind = 8 ) level_weight(dim_num)
  integer ( kind = 4 ) np(dim_num)
  integer ( kind = 4 ) np_sum
  real ( kind = 8 ) p(*)
  integer ( kind = 4 ) p_index
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num
  integer ( kind = 4 ) rule(dim_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_index
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_point
  integer ( kind = 4 ), allocatable, dimension ( : ) :: sparse_unique_index
  real ( kind = 8 ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_POINT_TEST'
  write ( *, '(a)' ) '  SGMGA_POINT returns an array of the points'
  write ( *, '(a)' ) '  forming a multidimensional sparse grid with mixed factors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each sparse grid is of spatial dimension DIM_NUM,'
  write ( *, '(a)' ) '  and is made up of product grids of levels up to LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IMPORTANCE:'
  write ( *, '(5g14.6)' ) importance(1:dim_num)
  write ( *, '(a)' ) '  LEVEL_WEIGHT:'
  write ( *, '(5g14.6)' ) level_weight(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Dimension      Rule       Growth     Parameters'
  write ( *, '(a)' ) ' '

  p_index = 1

  do dim = 1, dim_num

    if ( rule(dim) == 1 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    else if ( rule(dim) ==  2 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    else if ( rule(dim) ==  3 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    else if ( rule(dim) ==  4 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    else if ( rule(dim) ==  5 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    else if ( rule(dim) ==  6 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim), p(p_index)
    else if ( rule(dim) ==  7 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    else if ( rule(dim) ==  8 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim), p(p_index)
    else if ( rule(dim) ==  9 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim), p(p_index), p(p_index+1)
    else if ( rule(dim) == 10 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    else if ( rule(dim) == 11 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    else if ( rule(dim) == 12 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), growth(dim)
    end if

    p_index = p_index + np(dim)

  end do

  do level_max = level_max_min, level_max_max

    call sgmga_size_total ( dim_num, level_weight, level_max, rule, growth, &
      point_total_num )

    call sgmga_size ( dim_num, level_weight, level_max, rule, growth, np, p, &
      tol, point_num )

    allocate ( sparse_unique_index(1:point_total_num) )

    call sgmga_unique_index ( dim_num, level_weight, level_max, rule, growth, np, &
      p, tol, point_num, point_total_num, sparse_unique_index )

    allocate ( sparse_index(1:dim_num,1:point_num) )
    allocate ( sparse_order(1:dim_num,1:point_num) )

    call sgmga_index ( dim_num, level_weight, level_max, rule, growth, point_num, &
      point_total_num, sparse_unique_index, sparse_order, sparse_index )

    allocate ( sparse_point(1:dim_num,1:point_num) )

    call sgmga_point ( dim_num, level_weight, level_max, rule, growth, np, p, &
      point_num, sparse_order, sparse_index, sparse_point )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  For LEVEL_MAX = ', level_max
    write ( *, '(a)' ) ' '
    do point = 1, point_num
      write ( *, '(2x,i4,2x,5g12.6)' ) &
        point, sparse_point(1:dim_num,point)
    end do

    deallocate ( sparse_index )
    deallocate ( sparse_order )
    deallocate ( sparse_point )
    deallocate ( sparse_unique_index )

  end do

  return
end
