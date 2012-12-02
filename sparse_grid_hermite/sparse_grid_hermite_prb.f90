program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_HERMITE_PRB.
!
!  Discussion:
!
!    SPARSE_GRID_HERMITE_PRB tests the SPARSE_GRID_HERMITE routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2012
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

  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_HERMITE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPARSE_GRID_HERMITE library.'
!
!  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
!
  if ( .true. ) then

  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10

  call test01 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5

  call test01 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 100
  dim_max = 100
  level_max_min = 0
  level_max_max = 2

  call test01 ( dim_min, dim_max, level_max_min, level_max_max )

  end if
!
!  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
!
  if ( .true. ) then

  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10

  call test02 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5

  call test02 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 100
  dim_max = 100
  level_max_min = 0
  level_max_max = 2

  call test02 ( dim_min, dim_max, level_max_min, level_max_max )

  end if
!
!  Compute abstract grid indices of sparse grid points as selected from product grid
!  for DIMENSION, LEVEL_MAX.
!
  if ( .true. ) then

  dim_num = 2
  level_max = 3
  call test03 ( dim_num, level_max )

  dim_num = 2
  level_max = 4
  call test03 ( dim_num, level_max )

  dim_num = 3
  level_max = 0
  call test03 ( dim_num, level_max )

  dim_num = 3
  level_max = 2
  call test03 ( dim_num, level_max )

  dim_num = 6
  level_max = 2
  call test03 ( dim_num, level_max )

  end if
!
!  Compute sparse Gauss-Hermite rule for DIMENSION, LEVEL_MAX.
!
  if ( .true. ) then

  dim_num = 2
  level_max = 0
  call test04 ( dim_num, level_max )

  dim_num = 2
  level_max = 3
  call test04 ( dim_num, level_max )

  dim_num = 2
  level_max = 4
  call test04 ( dim_num, level_max )

  dim_num = 3
  level_max = 0
  call test04 ( dim_num, level_max )

  dim_num = 3
  level_max = 2
  call test04 ( dim_num, level_max )

  end if
!
!  Test sum of weights for DIMENSION, LEVEL_MAX.
!
  if ( .true. ) then

  dim_num = 2
  level_max = 4
  call test05 ( dim_num, level_max )

  dim_num = 3
  level_max = 0
  call test05 ( dim_num, level_max )

  dim_num = 3
  level_max = 1
  call test05 ( dim_num, level_max )

  dim_num = 3
  level_max = 6
  call test05 ( dim_num, level_max )

  dim_num = 10
  level_max = 3
  call test05 ( dim_num, level_max )

  end if
!
!  Test monomial exactness for DIMENSION, LEVEL_MAX, DEGREE_MAX.
!
  if ( .true. ) then

  dim_num = 2
  level_max = 0
  degree_max = 3
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 1
  degree_max = 5
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 2
  degree_max = 8
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 3
  degree_max = 10
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 4
  degree_max = 13
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 5
  degree_max = 13
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 3
  level_max = 0
  degree_max = 2
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 3
  level_max = 1
  degree_max = 9
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 3
  level_max = 2
  degree_max = 9
  call test06 ( dim_num, level_max, degree_max )

  dim_num = 3
  level_max = 3
  degree_max = 9
  call test06 ( dim_num, level_max, degree_max )

  end if
!
!  Show how to write a rule to a file.
!
  if ( .true. ) then

  dim_num = 2
  level_max = 3

  call test07 ( dim_num, level_max )

  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_HERMITE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST01 tests SPARSE_GRID_HERMITE_SIZE_TOTAL.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_MIN, the minimum spatial dimension
!    to consider.
!
!    Input, integer ( kind = 4 ) DIM_MAX, the maximum spatial dimension
!    to consider.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MIN, the minimum value of
!    LEVEL_MAX to consider.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MAX, the maximum value of
!    LEVEL_MAX to consider.
!
  implicit none

  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SPARSE_GRID_HERMITE_SIZE_TOTAL returns the number of'
  write ( *, '(a)' ) '  points in a Gauss-Hermite sparse grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each sparse grid is of spatial dimension DIM,'
  write ( *, '(a)' ) '  and is made up of all product grids of levels up to LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a10,8(2x,i10))' ) '   DIM:   ', ( dim_num, dim_num = dim_min, dim_max )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max

    do dim_num = dim_min, dim_max
      call sparse_grid_hermite_size ( dim_num, level_max, point_num, &
        point_total_num(dim_num) )
    end do

    write ( *, '(8(2x,i10))' ) level_max, point_total_num(dim_min:dim_max)

  end do

  return
end
subroutine test02 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST02 tests SPARSE_GRID_HERMITE_SIZE.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_MIN, the minimum spatial dimension
!    to consider.
!
!    Input, integer ( kind = 4 ) DIM_MAX, the maximum spatial dimension
!    to consider.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MIN, the minimum value of
!    LEVEL_MAX to consider.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MAX, the maximum value of
!    LEVEL_MAX to consider.
!
  implicit none

  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)
  integer ( kind = 4 ) point_total_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  SPARSE_GRID_HERMITE_SIZE returns the number of distinct'
  write ( *, '(a)' ) '  points in a Gauss-Hermite sparse grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each sparse grid is of spatial dimension DIM,'
  write ( *, '(a)' ) '  and is made up of all product grids of levels up to LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a10,8(2x,i10))' ) '   DIM:   ', ( dim_num, dim_num = dim_min, dim_max )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max

    do dim_num = dim_min, dim_max
      call sparse_grid_hermite_size ( dim_num, level_max, &
        point_num(dim_num), point_total_num )
    end do

    write ( *, '(8(2x,i10))' ) level_max, point_num(dim_min:dim_max)

  end do

  return
end
subroutine test03 ( dim_num, level_max )

!****************************************************************************80
!
!! TEST03 tests SPARSE_GRID_HERMITE_INDEX.
!
!  Discussion:
!
!    The routine computes abstract indices that describe the sparse grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_index
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_order
  integer ( kind = 4 ), allocatable, dimension ( : ) :: sparse_unique_index
  real ( kind = 8 ) tol

  tol = 1.0D-07
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  SPARSE_GRID_HERMITE_INDEX returns abstract indices for the'
  write ( *, '(a)' ) '  points that make up a Hermite sparse grid.'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num

  call sparse_grid_hermite_size ( dim_num, level_max, point_num, point_total_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points in the grid        = ', point_total_num
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Unique indices.
!
  allocate ( sparse_unique_index(1:point_total_num) )

  call sparse_grid_hermite_unique_index ( dim_num, level_max, tol, point_num, &
    point_total_num, sparse_unique_index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sparse grid unique index:'
  write ( *, '(a)' ) ' '
  do point = 1, point_total_num
    write ( *, '(2x,i4,2x,10i6)' ) point, sparse_unique_index(point)
  end do
!
!  Compute the indices.
!
  allocate ( sparse_index(1:dim_num,1:point_num) )
  allocate ( sparse_order(1:dim_num,1:point_num) )

  call sparse_grid_hermite_index ( dim_num, level_max, point_num, point_total_num, &
    sparse_unique_index, sparse_order, sparse_index )
!
!  Print the indices.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sparse grid index/order:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,10i6)' ) point, sparse_index(1:dim_num,point)
    write ( *, '(2x,4x,2x,10i6)' )        sparse_order(1:dim_num,point)
  end do

  deallocate ( sparse_index )
  deallocate ( sparse_order )
  deallocate ( sparse_unique_index )

  return
end
subroutine test04 ( dim_num, level_max )

!****************************************************************************80
!
!! TEST04 computes the points and weights of a Hermite sparse grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: sparse_weight

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  SPARSE_GRID_HERMITE computes Hermite sparse grid points and weights.'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
!
!  Determine the number of points.
!
  call sparse_grid_hermite_size ( dim_num, level_max, point_num, point_total_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points in the grid        = ', point_total_num
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( sparse_weight(1:point_num) )
  allocate ( sparse_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid_hermite ( dim_num, level_max, point_total_num, point_num, &
    sparse_weight, sparse_point )
!
!  Print them out.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid weights:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,f10.6)' ) point, sparse_weight(point)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,5f10.6)' ) point, sparse_point(1:dim_num,point)
  end do

  deallocate ( sparse_point )
  deallocate ( sparse_weight )

  return
end
subroutine test05 ( dim_num, level_max )

!****************************************************************************80
!
!! TEST05 sums the weights and compares them to SQRT(PI**DIM_NUM).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: sparse_weight
  real ( kind = 8 ) weight_sum
  real ( kind = 8 ) weight_sum_error
  real ( kind = 8 ) weight_sum_exact

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  Compute the weights of a Gauss-Hermite sparse grid .'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  As a simple test, sum these weights.'
  write ( *, '(a)' ) '  They should sum to SQRT(PI^DIM_NUM).'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
!
!  Determine the number of points.
!
  call sparse_grid_hermite_size ( dim_num, level_max, point_num, point_total_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points in the grid        = ', point_total_num
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( sparse_weight(1:point_num) )
  allocate ( sparse_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid_hermite ( dim_num, level_max, point_total_num, point_num, &
    sparse_weight, sparse_point )
!
!  Sum the weights.
!
  weight_sum = sum ( sparse_weight(1:point_num) )

  weight_sum_exact = sqrt ( pi**dim_num )

  weight_sum_error = abs ( weight_sum - weight_sum_exact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Weight sum     Exact sum    Difference'
  write ( *, '(a)' ) ' '
  write ( *, '(a2,g14.6,a2,g14.6,a2,g14.6)' ) &
    '  ', weight_sum, '  ', weight_sum_exact, '  ', weight_sum_error

  deallocate ( sparse_point )
  deallocate ( sparse_weight )

  return
end
subroutine test06 ( dim_num, level_max, degree_max )

!****************************************************************************80
!
!! TEST06 tests a Gauss-Hermite sparse grid rule for monomial exactness.
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
!    15 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  real ( kind = 8 ) exact
  integer ( kind = 4 ), allocatable, dimension ( : ) :: expon
  integer ( kind = 4 ) h
  integer ( kind = 4 ) last
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num
  real ( kind = 8 ) quad
  real ( kind = 8 ) quad_error
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: sparse_weight
  integer ( kind = 4 ) t
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Check the exactness of a Gauss-Hermite sparse '
  write ( *, '(a)' ) '  grid quadrature rule, applied to all monomials '
  write ( *, '(a)' ) '  of orders 0 to DEGREE_MAX.'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  The maximum total degree to be checked is DEGREE_MAX = ', degree_max
!
!  Determine the number of points.
!
  call sparse_grid_hermite_size ( dim_num, level_max, point_num, point_total_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points in the grid        = ', point_total_num
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( sparse_weight(1:point_num) )
  allocate ( sparse_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid_hermite ( dim_num, level_max, point_total_num, point_num, &
    sparse_weight, sparse_point )
!
!  Explore the monomials.
!
  allocate ( expon(1:dim_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Exact        Estimate         Error       Total   Monomial'
  write ( *, '(a)' ) '                                                Degree  Exponents'
  write ( *, '(a)' ) '  --------------  --------------  --------------  --    --------------------'
  write ( *, '(a)' ) ' '

  do degree = 0, degree_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( degree, dim_num, expon, more, h, t )

      call monomial_quadrature_hermite ( dim_num, expon, point_num, &
        sparse_weight, sparse_point, exact, quad, quad_error )

      write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,i2,2x,10i3)' ) &
        exact, quad, quad_error, degree, expon(1:dim_num)

      if ( .not. more ) then
        exit
      end if

    end do

    write ( *, '(a)' ) ' '

  end do

  deallocate ( expon )
  deallocate ( sparse_point )
  deallocate ( sparse_weight )

  return
end
subroutine test07 ( dim_num, level_max )

!****************************************************************************80
!
!! TEST07 creates a sparse Gauss-Hermite grid and writes it to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  character ( len = 255 ) filename
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: sparse_weight

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  Call SPARSE_GRID_HERMITE_WRITE to write'
  write ( *, '(a)' ) '  sparse grid data to a set of files.'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
!
!  Determine the number of points.
!
  call sparse_grid_hermite_size ( dim_num, level_max, point_num, point_total_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points in the grid        = ', point_total_num
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( sparse_weight(1:point_num) )
  allocate ( sparse_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid_hermite ( dim_num, level_max, point_total_num, point_num, &
    sparse_weight, sparse_point )
!
!  Write the data out.
!
  filename = 'test07'

  call sparse_grid_hermite_write ( dim_num, point_num, sparse_weight, &
    sparse_point, filename )

  deallocate ( sparse_weight )
  deallocate ( sparse_point )

  return
end
