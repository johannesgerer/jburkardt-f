program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_OPEN_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2010
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

  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_OPEN_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPARSE_GRID_OPEN library.'
!
!  Point growth table for general Open Fully Nested "OFN" rules.
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test01 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test01 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  Point growth table for general Open Non Nested "ONN" rules.
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test011 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test011 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  Point growth table for general Open Weaky Nested "OWN" rules.
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test012 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test012 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  Point growth table for Fejer Type 2 Slow rules.
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test013 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test013 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  Point growth table for Gauss-Patterson-Slow rules.
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test015 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test015 ( dim_min, dim_max, level_max_min, level_max_max )

  call test02 ( 2, 2 )
  call test02 ( 2, 3 )
  call test02 ( 2, 4 )
  call test02 ( 3, 2 )
  call test02 ( 6, 2 )

  call test04 ( 2, 3 )

  call test05 ( 2, 3 )

  call test06 ( 2, 4 )

  call test08 ( 2, 1 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_OPEN_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST01 tests SPARSE_GRID_OFN_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SPARSE_GRID_OFN_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open fully nested '
  write ( *, '(a)' ) '  quadrature rules.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sparse grid is the sum of all product grids'
  write ( *, '(a)' ) '  of order LEVEL, with'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  the order of the 1D rule is 2^(LEVEL+1) - 1,'
  write ( *, '(a)' ) '  the region is [-1,1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this kind of rule, there is complete nesting,'
  write ( *, '(a)' ) '  that is, a sparse grid of a given level includes'
  write ( *, '(a)' ) '  ALL the points on grids of lower levels.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_ofn_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test011 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST011 tests SPARSE_GRID_ONN_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  SPARSE_GRID_ONN_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open non-nested '
  write ( *, '(a)' ) '  quadrature rules.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sparse grid is the sum of all product grids'
  write ( *, '(a)' ) '  of order LEVEL, with'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  the order of the 1D rule is 2^(LEVEL+1) - 1,'
  write ( *, '(a)' ) '  the region is [-1,1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this kind of rule, there is no nesting.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_onn_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test012 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST012 tests SPARSE_GRID_OWN_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  SPARSE_GRID_OWN_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open weakly nested '
  write ( *, '(a)' ) '  quadrature rules.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sparse grid is the sum of all product grids'
  write ( *, '(a)' ) '  of order LEVEL, with'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  the order of the 1D rule is 2^(LEVEL+1) - 1,'
  write ( *, '(a)' ) '  the region is [-1,1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this kind of rule, there is weak nesting,'
  write ( *, '(a)' ) '  that is, a sparse grid of a given level includes'
  write ( *, '(a)' ) '  only the abscissa 0 from the previous level.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_own_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test013 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST013 tests SPARSE_GRID_F2S_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  SPARSE_GRID_F2S_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid made up of product'
  write ( *, '(a)' ) '  grids formed from 1D Fejer Type 2 Slow rules.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sparse grid is the sum of all product grids'
  write ( *, '(a)' ) '  of order LEVEL, with'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  the order of the 1D rule is 2^(LEVEL+1) - 1,'
  write ( *, '(a)' ) '  the region is [-1,1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this kind of rule, there is complete nesting,'
  write ( *, '(a)' ) '  that is, a sparse grid of a given level includes'
  write ( *, '(a)' ) '  ALL the points on grids of lower levels.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_f2s_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test015 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST015 tests SPARSE_GRID_GPS_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  SPARSE_GRID_GPS_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid made up of product'
  write ( *, '(a)' ) '  grids formed from 1D Gauss-Patterson-Slow rules.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sparse grid is the sum of all product grids'
  write ( *, '(a)' ) '  of order LEVEL, with'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  the order of the 1D rule is 2^(LEVEL+1) - 1,'
  write ( *, '(a)' ) '  the region is [-1,1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this kind of rule, there is complete nesting,'
  write ( *, '(a)' ) '  that is, a sparse grid of a given level includes'
  write ( *, '(a)' ) '  ALL the points on grids of lower levels.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_gps_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test02 ( dim_num, level_max )

!*****************************************************************************80
!
!! TEST02 tests LEVELS_OPEN_INDEX.
!
!  Discussion:
!
!    The routine under study computes the indices of the unique points
!    used in a sparse multidimensional grid whose size is controlled
!    by a parameter LEVEL.
!
!    Once these indices are returned, they can be interpreted as the
!    indices of points in a product grid based on Fejer Type 2,
!    Newton Cotes Open, or Gauss Patterson rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  integer ( kind = 4 ) grid_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point_num

  level_min = level_max + 1 - dim_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  LEVELS_OPEN_INDEX returns all grid indexes'
  write ( *, '(a)' ) '  whose level value satisfies'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) '  Here, LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  and the order of the rule is 2^(LEVEL+1) - 1.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max

  call sparse_grid_ofn_size ( dim_num, level_max, point_num )

  write ( *, '(a,i8)' ) '  Unique points in the grid = ', point_num
!
!  Allocate the space.
!
  allocate ( grid_index(1:dim_num,1:point_num) )
!
!  Compute the grid index values.
!
  call levels_open_index ( dim_num, level_max, point_num, grid_index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid index:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i4,2x,6i6)' ) j, grid_index(1:dim_num,j)
  end do

  deallocate ( grid_index)

  return
end
subroutine test04 ( dim_num, level_max )

!*****************************************************************************80
!
!! TEST04 tests LEVELS_OPEN_INDEX to make a Fejer type 2 grid.
!
!  Discussion:
!
!    This routine gets the sparse grid indices and determines the
!    corresponding sparse grid abscissas for a Fejer type 2 scheme.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) f2_abscissa
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  integer ( kind = 4 ) grid_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_num

  level_min = level_max + 1 - dim_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  Make a Fejer type 2 sparse grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEVELS_OPEN_INDEX returns all grid indexes'
  write ( *, '(a)' ) '  whose level value satisfies'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) '  Here, LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  and the order of the rule is 2^(LEVEL+1) - 1..'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now we demonstrate how to convert grid indices'
  write ( *, '(a)' ) '  into physical grid points.  In this case, we'
  write ( *, '(a)' ) '  want points on [-1,+1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max

  call sparse_grid_ofn_size ( dim_num, level_max, point_num )

  write ( *, '(a,i8)' ) '  Unique points in the grid = ', point_num
!
!  Allocate the space.
!
  allocate ( grid_index(1:dim_num,1:point_num) )
!
!  Compute the grid index values.
!
  call levels_open_index ( dim_num, level_max, point_num, grid_index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid index:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i4,2x,6i6)' ) j, grid_index(1:dim_num,j)
  end do
!
!  Convert index information to physical information.
!
  order_max = 2**( level_max + 1 ) - 1

  allocate ( grid_point(1:dim_num,1:point_num) )

  do j = 1, point_num
    do dim = 1, dim_num
      grid_point(dim,j) = f2_abscissa ( order_max, grid_index(dim,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i8,2x,6f10.6)' ) j, grid_point(1:dim_num,j)
  end do

  deallocate ( grid_index )
  deallocate ( grid_point )

  return
end
subroutine test05 ( dim_num, level_max )

!*****************************************************************************80
!
!! TEST05 tests LEVELS_OPEN_INDEX to make a Gauss Patterson grid.
!
!  Discussion:
!
!    This routine gets the sparse grid indices and determines the
!    corresponding sparse grid abscissas for a Gauss-Patterson scheme.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) gp_abscissa
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  integer ( kind = 4 ) grid_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_num

  level_min = level_max + 1 - dim_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  Make a Gauss-Patterson sparse grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEVELS_OPEN_INDEX returns all grid indexes'
  write ( *, '(a)' ) '  whose level value satisfies'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) '  Here, LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  and the order of the rule is 2^(LEVEL+1) - 1..'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now we demonstrate how to convert grid indices'
  write ( *, '(a)' ) '  into physical grid points.  In this case, we'
  write ( *, '(a)' ) '  want points on [-1,+1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max

  call sparse_grid_ofn_size ( dim_num, level_max, point_num )

  write ( *, '(a,i8)' ) '  Unique points in the grid = ', point_num
!
!  Allocate the space.
!
  allocate ( grid_index(1:dim_num,1:point_num) )
!
!  Compute the grid index values.
!
  write ( * , * ) 'DEBUG: Call LEVELS_OPEN_INDEX'

  call levels_open_index ( dim_num, level_max, point_num, grid_index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid index:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i4,2x,6i6)' ) j, grid_index(1:dim_num,j)
  end do
!
!  Convert index information to physical information.
!  Note that GP_ABSCISSA expects the LEVEL value, not the ORDER.
!
  allocate ( grid_point(1:dim_num,1:point_num) )

  order_max = 2**(level_max+1) - 1

  do j = 1, point_num
    do dim = 1, dim_num
      grid_point(dim,j) = gp_abscissa ( order_max, grid_index(dim,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i8,2x,6f10.6)' ) j, grid_point(1:dim_num,j)
  end do

  deallocate ( grid_index )
  deallocate ( grid_point )

  return
end
subroutine test06 ( dim_num, level_max )

!*****************************************************************************80
!
!! TEST06 tests LEVELS_OPEN_INDEX to make a Newton Cotes Open grid.
!
!  Discussion:
!
!    This routine gets the sparse grid indices and determines the
!    corresponding sparse grid abscissas for a Newton Cotes Open scheme.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  integer ( kind = 4 ) grid_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  real ( kind = 8 ) nco_abscissa
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_num

  level_min = level_max + 1 - dim_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  Make a Newton Cotes Open sparse grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEVELS_OPEN_INDEX returns all grid indexes'
  write ( *, '(a)' ) '  whose level value satisfies'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) '  Here, LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  and the order of the rule is 2^(LEVEL+1) - 1..'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now we demonstrate how to convert grid indices'
  write ( *, '(a)' ) '  into physical grid points.  In this case, we'
  write ( *, '(a)' ) '  want points on [0,+1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max

  call sparse_grid_ofn_size ( dim_num, level_max, point_num )

  write ( *, '(a,i8)' ) '  Unique points in the grid = ', point_num
!
!  Allocate the space.
!
  allocate ( grid_index(1:dim_num,1:point_num) )
!
!  Compute the grid index values.
!
  call levels_open_index ( dim_num, level_max, point_num, grid_index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid index:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i4,2x,6i6)' ) j, grid_index(1:dim_num,j)
  end do
!
!  Convert index information to physical information.
!
  order_max = 2**( level_max + 1 ) - 1

  allocate ( grid_point(1:dim_num,1:point_num) )

  do j = 1, point_num
    do dim = 1, dim_num
      grid_point(dim,j) = nco_abscissa ( order_max, grid_index(dim,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i8,2x,6f10.6)' ) j, grid_point(1:dim_num,j)
  end do

  deallocate ( grid_index )
  deallocate ( grid_point )

  return
end
subroutine test08 ( dim_num, level_max )

!*****************************************************************************80
!
!! TEST08 creates and writes sparse grid files of all types.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ) dim
  integer   ( kind = 4 ) dim_num
  real      ( kind = 8 ) f2_abscissa
  real      ( kind = 8 ) gp_abscissa
  character ( len = 80 ) file_name
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  integer   ( kind = 4 ) grid_num
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) level
  integer   ( kind = 4 ) level_max
  integer   ( kind = 4 ) level_min
  real      ( kind = 8 ) nco_abscissa
  integer   ( kind = 4 ) order_max
  integer   ( kind = 4 ) point_num

  level_min = level_max + 1 - dim_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  Make sparse grids and write to files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  LEVEL_MIN =                 ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX =                 ', level_max

  call sparse_grid_ofn_size ( dim_num, level_max, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Allocate the space.
!
  allocate ( grid_index(1:dim_num,1:point_num) )
!
!  Compute the orders and points.
!
  call levels_open_index ( dim_num, level_max, point_num, grid_index )
!
!  Now we're done.  Print the merged grid data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid index:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i4,2x,6i6)' ) j, grid_index(1:dim_num,j)
  end do
!
!  Convert index information to physical information.
!
  order_max = 2**( level_max + 1 ) - 1

  allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Create F2 data and write to file.
!
  do j = 1, point_num
    do dim = 1, dim_num
      grid_point(dim,j) = f2_abscissa ( order_max, grid_index(dim,j)  )
    end do
  end do

  write ( file_name, '(a,i2,a,i2,a)' ) &
    'f2_d', dim_num, '_level', level_max, '.txt'

  call s_blank_delete ( file_name )

  call r8mat_write ( file_name, dim_num, point_num, grid_point )

  write ( *, '(a)' ) '  Wrote file "' // trim ( file_name ) // '".'
!
!  Create GP data and write to file.
!  Note that GP_ABSCISSA wants the value of LEVEL_MAX, not ORDER_MAX!
!
  do j = 1, point_num
    do dim = 1, dim_num
      grid_point(dim,j) = gp_abscissa ( order_max, grid_index(dim,j) )
    end do
  end do

  write ( file_name, '(a,i2,a,i2,a)' ) &
    'gp_d', dim_num, '_level', level_max, '.txt'

  call s_blank_delete ( file_name )

  call r8mat_write ( file_name, dim_num, point_num, grid_point )

  write ( *, '(a)' ) '  Wrote file "' // trim ( file_name ) // '".'
!
!  Create NCO data and write to file.
!
  do j = 1, point_num
    do dim = 1, dim_num
      grid_point(dim,j) = nco_abscissa ( order_max, grid_index(dim,j) )
    end do
  end do

  write ( file_name, '(a,i2,a,i2,a)' ) &
    'nco_d', dim_num, '_level', level_max, '.txt'

  call s_blank_delete ( file_name )

  call r8mat_write ( file_name, dim_num, point_num, grid_point )

  write ( *, '(a)' ) '  Wrote file "' // trim ( file_name ) // '".'

  deallocate ( grid_index )
  deallocate ( grid_point )

  return
end
