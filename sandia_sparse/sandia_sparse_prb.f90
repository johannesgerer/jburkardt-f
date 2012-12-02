program main

!****************************************************************************80
!
!! MAIN is the main program for SANDIA_SPARSE_PRB.
!
!  Discussion:
!
!    SANDIA_SPARSE_PRB tests the SANDIA_SPARSE routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), parameter :: rule_max = 7

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_SPARSE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SANDIA_SPARSE library.'
!
!  Test LEVELS_INDEX_SIZE for one example each of CFN, OFN, OWN and ONN rules.
!
  if ( .true. ) then

    rule = 1

    dim_min = 1
    dim_max = 1
    level_max_min = 0
    level_max_max = 10

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 1
    dim_max = 6
    level_max_min = 0
    level_max_max = 6

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 6
    dim_max = 10
    level_max_min = 0
    level_max_max = 5

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 100
    dim_max = 100
    level_max_min = 0
    level_max_max = 2

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    rule = 2

    dim_min = 1
    dim_max = 1
    level_max_min = 0
    level_max_max = 10

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 1
    dim_max = 6
    level_max_min = 0
    level_max_max = 6

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 6
    dim_max = 10
    level_max_min = 0
    level_max_max = 5

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 100
    dim_max = 100
    level_max_min = 0
    level_max_max = 2

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    rule = 5

    dim_min = 1
    dim_max = 1
    level_max_min = 0
    level_max_max = 10

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 1
    dim_max = 6
    level_max_min = 0
    level_max_max = 6

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 6
    dim_max = 10
    level_max_min = 0
    level_max_max = 5

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 100
    dim_max = 100
    level_max_min = 0
    level_max_max = 2

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    rule = 7

    dim_min = 1
    dim_max = 1
    level_max_min = 0
    level_max_max = 10

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 1
    dim_max = 6
    level_max_min = 0
    level_max_max = 6

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 6
    dim_max = 10
    level_max_min = 0
    level_max_max = 5

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

    dim_min = 100
    dim_max = 100
    level_max_min = 0
    level_max_max = 2

    call levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
      level_max_max )

  end if
!
!  Test LEVELS_INDEX for one example each of CFN, OFN, OWN and ONN rules.
!
  if ( .true. ) then

    rule = 1

    dim_num = 2
    level_max = 1
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 2
    level_max = 3
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 3
    level_max = 0
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 3
    level_max = 2
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 6
    level_max = 2
    call levels_index_test ( rule, dim_num, level_max )

    rule = 2

    dim_num = 2
    level_max = 1
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 2
    level_max = 3
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 3
    level_max = 0
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 3
    level_max = 2
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 6
    level_max = 2
    call levels_index_test ( rule, dim_num, level_max )

    rule = 5

    dim_num = 2
    level_max = 1
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 2
    level_max = 3
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 3
    level_max = 0
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 3
    level_max = 2
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 6
    level_max = 2
    call levels_index_test ( rule, dim_num, level_max )

    rule = 7

    dim_num = 2
    level_max = 1
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 2
    level_max = 3
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 3
    level_max = 0
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 3
    level_max = 2
    call levels_index_test ( rule, dim_num, level_max )

    dim_num = 6
    level_max = 2
    call levels_index_test ( rule, dim_num, level_max )

  end if
!
!  Test SPARSE_GRID by having it compute a few sparse grids based on each rule.
!
  if ( .true. ) then

    do rule = 1, rule_max

      dim_num = 2
      level_max = 1
      call sparse_grid_compute_test ( rule, dim_num, level_max )

      dim_num = 2
      level_max = 2
      call sparse_grid_compute_test ( rule, dim_num, level_max )

      dim_num = 3
      level_max = 1
      call sparse_grid_compute_test ( rule, dim_num, level_max )

    end do

  end if
!
!  Test SPARSE_GRID by having it compute a few sparse grids based on each rule,
!  and comparing the sum of the quadrature weights to the expected sum.
!
  if ( .true. ) then

    do rule = 1, rule_max

      dim_num = 2
      level_max = 4
      call sparse_grid_weight_test ( rule, dim_num, level_max )

      dim_num = 3
      level_max = 0
      call sparse_grid_weight_test ( rule, dim_num, level_max )

      dim_num = 3
      level_max = 1
      call sparse_grid_weight_test ( rule, dim_num, level_max )

      dim_num = 3
      level_max = 6
      call sparse_grid_weight_test ( rule, dim_num, level_max )

      dim_num = 10
      level_max = 3
      call sparse_grid_weight_test ( rule, dim_num, level_max )

    end do

  end if
!
!  Test SPARSE_GRID by having it compute a few sparse grids based on each rule,
!  and comparing estimated and exact monomial integrals.
!
  if ( .true. ) then

    do rule = 1, rule_max

      dim_num = 2
      level_max = 0
      degree_max = 3
      call sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max )

      dim_num = 2
      level_max = 1
      degree_max = 5
      call sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max )

      dim_num = 2
      level_max = 2
      degree_max = 7
      call sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max )

      dim_num = 2
      level_max = 3
      degree_max = 9
      call sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max )

      dim_num = 3
      level_max = 0
      degree_max = 2
      call sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max )

      dim_num = 3
      level_max = 1
      degree_max = 4
      call sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max )

      dim_num = 3
      level_max = 2
      degree_max = 6
      call sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max )

    end do

  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_SPARSE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine levels_index_size_test ( rule, dim_min, dim_max, level_max_min, &
  level_max_max )

!****************************************************************************80
!
!! LEVELS_INDEX_SIZE_TEST tests LEVELS_INDEX_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) DIM_MIN, the minimum spatial dimension.
!
!    Input, integer ( kind = 4 ) DIM_MAX, the maximum spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MIN, the minimum value of LEVEL_MAX.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MAX, the maximum value of LEVEL_MAX.
!
  implicit none

  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)
  integer ( kind = 4 ) rule

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEVELS_INDEX_SIZE_TEST'
  write ( *, '(a)' ) '  LEVELS_INDEX_SIZE returns the number of distinct'
  write ( *, '(a)' ) '  points in a sparse grid derived from a 1D rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  We are looking at rules like rule ', rule
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each sparse grid is of spatial dimension DIM,'
  write ( *, '(a)' ) '  and is made up of product grids such that'
  write ( *, '(a)' ) '  LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i8))' ) '   DIM: ', point_num(dim_min:dim_max)
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) '   ---------'

  do level_max = level_max_min, level_max_max

    do dim_num = dim_min, dim_max
      call levels_index_size ( dim_num, level_max, rule, point_num(dim_num) )
    end do

    write ( *, '(a4,i4,6(2x,i8))' ) '    ', level_max, point_num(dim_min:dim_max)

  end do

  return
end
subroutine levels_index_test ( rule, dim_num, level_max )

!****************************************************************************80
!
!! LEVELS_INDEX_TEST tests LEVELS_INDEX.
!
!  Discussion:
!
!    The routine computes the indices of the unique points used in a sparse
!    multidimensional grid whose size is controlled by a parameter LEVEL_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_base
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) rule

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEVELS_INDEX_TEST'
  write ( *, '(a)' ) '  LEVELS_INDEX returns all grid indexes'
  write ( *, '(a)' ) '  whose level value satisfies'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) '  Here, LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  and the order of the rule is 2^LEVEL + 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  We are looking at rules like rule ', rule
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max

  call levels_index_size ( dim_num, level_max, rule, point_num )

  write ( *, '(a,i8)' ) '  Unique points in the grid = ', point_num
!
!  Compute the orders and points.
!
  allocate ( grid_base(1:dim_num,1:point_num) )
  allocate ( grid_index(1:dim_num,1:point_num) )

  call levels_index ( dim_num, level_max, rule, point_num, grid_index, &
    grid_base )
!
!  Now we're done.  Print the merged grid data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Point     Grid indices:'
  write ( *, '(a)' ) '           Grid bases:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,10i6)' ) point, grid_index(1:dim_num,point)
    write ( *, '(2x,4x,2x,10i6)' )        grid_base(1:dim_num,point)
  end do

  deallocate ( grid_base )
  deallocate ( grid_index )

  return
end
subroutine sparse_grid_compute_test ( rule, dim_num, level_max )

!****************************************************************************80
!
!! SPARSE_GRID_COMPUTE_TEST computes and prints a sparse grid rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) rule

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_COMPUTE_TEST:'
  write ( *, '(a)' ) '  SPARSE_GRID can make a sparse grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max
  write ( *, '(a,i8)' ) '  1D quadrature index RULE =  ', rule
!
!  Determine the number of points.
!
  call levels_index_size ( dim_num, level_max, rule, point_num )

  write ( *, '(a,i8)' ) '  Unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( grid_weight(1:point_num) )
  allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid ( dim_num, level_max, rule, point_num, grid_weight, &
    grid_point )
!
!  Print them out.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid weights:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,g14.6)' ) point, grid_weight(point)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,5g14.6)' ) point, grid_point(1:dim_num,point)
  end do

  deallocate ( grid_point )
  deallocate ( grid_weight )

  return
end
subroutine sparse_grid_weight_test ( rule, dim_num, level_max )

!****************************************************************************80
!
!! SPARSE_GRID_WEIGHT_TEST checks the sum of the quadrature weights.
!
!  Discussion:
!
!    This routine gets the sparse grid indices and determines the
!    corresponding sparse grid abscissas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level.
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) rule
  real ( kind = 8 ) weight_sum
  real ( kind = 8 ) weight_sum_error
  real ( kind = 8 ) weight_sum_exact

  level_min = max ( 0, level_max + 1 - dim_num )

  if ( 1 <= rule .and. rule <= 5 ) then
    weight_sum_exact = 2.0D+00**dim_num
  else if ( 6 == rule ) then
    weight_sum_exact = sqrt ( pi**dim_num )
  else if ( 7 == rule ) then
    weight_sum_exact = 1.0D+00
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_WEIGHT_TEST:'
  write ( *, '(a)' ) '  Compute the weights of a sparse grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  As a simple test, sum these weights.'
  write ( *, '(a,g14.6)' ) '  They should sum to exactly ', weight_sum_exact
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max
  write ( *, '(a,i8)' ) '  1D quadrature index RULE =  ', rule
!
!  Determine the number of points.
!
  call levels_index_size ( dim_num, level_max, rule, point_num )

  write ( *, '(a,i8)' ) '  Unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( grid_weight(1:point_num) )
  allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid ( dim_num, level_max, rule, point_num, grid_weight, &
    grid_point )
!
!  Sum the weights.
!
  weight_sum = sum ( grid_weight(1:point_num) )

  weight_sum_error = abs ( weight_sum - weight_sum_exact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Weight sum  Expected sum    Difference'
  write ( *, '(a)' ) ' '
  write ( *, '(a2,g14.6,a2,g14.6,a2,g14.6)' ) &
    '  ', weight_sum, '  ', weight_sum_exact, '  ', weight_sum_error

  deallocate ( grid_point )
  deallocate ( grid_weight )

  return
end
subroutine sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max )

!****************************************************************************80
!
!! SPARSE_GRID_MONOMIAL_TEST tests monomial exactness of the sparse grid rules.
!
!  Discussion:
!
!    This test is going to check EVERY monomial of total degree DEGREE_MAX
!    or less.  Even for a moderately high dimension of DIM_NUM = 10, you
!    do NOT want to use a large value of DEGREE_MAX, since there are
!
!      1         monomials of total degree 0,
!      DIM_NUM   monomials of total degree 1,
!      DIM_NUM^2 monomials of total degree 2,
!      DIM_NUM^3 monomials of total degree 3, and so on.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level.
!
!    Input, integer ( kind = 4 ) DEGREE_MAX, the maximum monomial total
!    degree to check.
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: expon
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight
  integer ( kind = 4 ) h
  integer ( kind = 4 ) last
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) quad_error
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) t
  real ( kind = 8 ) volume

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MONOMIAL_TEST'
  write ( *, '(a)' ) '  Check the exactness of a sparse grid quadrature rule,'
  write ( *, '(a)' ) '  applied to all monomials of orders 0 to DEGREE_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For cases where the dimension is greater than 1,'
  write ( *, '(a)' ) '  many sparse grid of this level have accuracy through'
  write ( *, '(a,i8)' ) '  monomials of total degree   ', 2 * level_max + 1
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max
  write ( *, '(a,i8)' ) '  1D quadrature index RULE =  ', rule
  write ( *, '(a,i8)' ) '  Check up to DEGREE_MAX =    ', degree_max
!
!  Determine the number of points in the rule.
!
  call levels_index_size ( dim_num, level_max, rule, point_num )

  write ( *, '(a,i8)' ) '  Unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( grid_weight(1:point_num) )
  allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid ( dim_num, level_max, rule, point_num, grid_weight, &
    grid_point )
!
!  Compare exact and estimated values of the integrals of various monomials.
!
  allocate ( expon(1:dim_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Error      Total   Monomial'
  write ( *, '(a)' ) '                 Degree  Exponents'
  write ( *, '(a)' ) ' '

  do degree = 0, degree_max

    more = .false.

    do

      call comp_next ( degree, dim_num, expon, more, h, t )

      call monomial_quadrature ( dim_num, expon, point_num, &
        grid_weight, grid_point, rule, quad_error )

      write ( *, '(a,g14.6,a,i2,a,10i2)' ) &
        '  ', quad_error, '  ', degree, '    ', expon(1:dim_num)

      if ( .not. more ) then
        exit
      end if

    end do

    if ( 1 < dim_num ) then
      write ( *, '(a)' ) ' '
    end if

  end do

  deallocate ( expon )
  deallocate ( grid_point )
  deallocate ( grid_weight )

  return
end
