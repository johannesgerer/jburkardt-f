program main

!*****************************************************************************80
!
!!   MAIN is the main program for SPARSE_GRID_GL_PRB.
!
!  Discussion:
!
!    SPARSE_GRID_GL_PRB tests the SPARSE_GRID_GL routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2010
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
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_GL_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPARSE_GRID_GL library.'
!
!  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
!
  dim_min = 1
  dim_max = 1
  level_max_min = 0
  level_max_max = 10

  call test01 ( dim_min, dim_max, level_max_min, level_max_max )

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

  dim_min = 11
  dim_max = 15
  level_max_min = 0
  level_max_max = 5

  call test01 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 100
  dim_max = 100
  level_max_min = 0
  level_max_max = 2

  call test01 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  Compute abstract grid indices of sparse grid points as selected from product grid
!  for DIMENSION, LEVEL_MAX.
!
  dim_num = 2
  level_max = 3
  call test02 ( dim_num, level_max )

  dim_num = 2
  level_max = 4
  call test02 ( dim_num, level_max )

  dim_num = 3
  level_max = 0
  call test02 ( dim_num, level_max )

  dim_num = 3
  level_max = 2
  call test02 ( dim_num, level_max )

  dim_num = 6
  level_max = 2
  call test02 ( dim_num, level_max )
!
!  Compute sparse Gauss-Legendre rule for DIMENSION, LEVEL_MAX.
!
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
!
!  Test sum of weights for DIMENSION, LEVEL_MAX.
!
  dim_num = 2
  level_max = 4
  call test04 ( dim_num, level_max )

  dim_num = 3
  level_max = 0
  call test04 ( dim_num, level_max )

  dim_num = 3
  level_max = 1
  call test04 ( dim_num, level_max )

  dim_num = 3
  level_max = 6
  call test04 ( dim_num, level_max )

  dim_num = 10
  level_max = 3
  call test04 ( dim_num, level_max )
!
!  Test monomial exactness for DIMENSION, LEVEL_MAX, DEGREE_MAX.
!
  dim_num = 2
  level_max = 0
  degree_max = 3
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 1
  degree_max = 5
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 2
  degree_max = 7
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 3
  degree_max = 9
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 4
  degree_max = 11
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 2
  level_max = 5
  degree_max = 13
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 3
  level_max = 0
  degree_max = 2
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 3
  level_max = 1
  degree_max = 4
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 3
  level_max = 2
  degree_max = 6
  call test05 ( dim_num, level_max, degree_max )

  dim_num = 3
  level_max = 3
  degree_max = 8
  call test05 ( dim_num, level_max, degree_max )
!
!  Show how to write a rule to a file.
!
  dim_num = 2
  level_max = 3

  call test06 ( dim_num, level_max )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_GL_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST01 tests SPARSE_GRID_GL_SIZE.
!
!  Discussion:
!
!    Note that, unlike most sparse grids, a sparse grid
!    based on Gauss-Legendre points is almost entirely
!    NOT nested.
!
!    Hence the point counts should be much higher than for a grid of
!    the same level, but using rules such as Fejer1 or Fejer2 or
!    Gauss-Patterson or Newton-Cotes-Open or Newton-Cotes-Open-Half.
!
!    Each sparse grid is of spatial dimension DIM,
!    and is made up of all product grids of levels up to LEVEL_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2010
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SPARSE_GRID_GL_SIZE returns the number of distinct'
  write ( *, '(a)' ) '  points in a Gauss-Legendre sparse grid.'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,a8,8(2x,i8))' ) '  DIM: ', ( dim_num, dim_num = dim_min, dim_max )
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max

    do dim_num = dim_min, dim_max
      call sparse_grid_gl_size ( dim_num, level_max, point_num(dim_num) )
    end do

    write ( *, '(8(2x,i8))' ) level_max, point_num(dim_min:dim_max)

  end do

  return
end
subroutine test02 ( dim_num, level_max )

!****************************************************************************80
!
!! TEST02 tests SPARSE_GRID_GL_INDEX.
!
!  Discussion:
!
!    The routine computes abstract indices that describe the sparse grid
!    of GL points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2007
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
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_base
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  SPARSE_GRID_GL_INDEX returns abstract indices for the'
  write ( *, '(a)' ) '  points that make up a Gauss-Legendre sparse grid.'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num

  call sparse_grid_gl_size ( dim_num, level_max, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Compute the indices.
!
  allocate ( grid_index(1:dim_num,1:point_num) )
  allocate ( grid_base(1:dim_num,1:point_num) )

  call sparse_grid_gl_index ( dim_num, level_max, point_num, grid_index, &
    grid_base )
!
!  Print the indices.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid index/base:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,10i6)' ) point, grid_index(1:dim_num,point)
    write ( *, '(2x,4x,2x,10i6)' )        grid_base(1:dim_num,point)
  end do

  deallocate ( grid_index )
  deallocate ( grid_base )

  return
end
subroutine test03 ( dim_num, level_max )

!****************************************************************************80
!
!! TEST03 call SPARSE_GRID_GL to create a sparse Gauss-Legendre grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2007
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
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  SPARSE_GRID_GL makes a sparse Gauss-Legendre grid.'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
!
!  Determine the number of points.
!
  call sparse_grid_gl_size ( dim_num, level_max, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( grid_weight(1:point_num) )
  allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid_gl ( dim_num, level_max, point_num, grid_weight, grid_point )
!
!  Print them out.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid weights:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,f10.6)' ) point, grid_weight(point)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do point = 1, point_num
    write ( *, '(2x,i4,2x,5f10.6)' ) point, grid_point(1:dim_num,point)
  end do

  deallocate ( grid_point )
  deallocate ( grid_weight )

  return
end
subroutine test04 ( dim_num, level_max )

!****************************************************************************80
!
!! TEST04 sums the weights and compares them to 2^DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2007
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
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) weight_sum
  real ( kind = 8 ) weight_sum_error
  real ( kind = 8 ) weight_sum_exact

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  Compute the weights of a Gauss-Legendre sparse grid .'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  As a simple test, sum these weights.'
  write ( *, '(a)' ) '  They should sum to exactly 2^DIM_NUM.'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
!
!  Determine the number of points.
!
  call sparse_grid_gl_size ( dim_num, level_max, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( grid_weight(1:point_num) )
  allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid_gl ( dim_num, level_max, point_num, grid_weight, grid_point )
!
!  Sum the weights.
!
  weight_sum = sum ( grid_weight(1:point_num) )

  weight_sum_exact = 2.0D+00**dim_num

  weight_sum_error = abs ( weight_sum - weight_sum_exact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Weight sum     Exact sum    Difference'
  write ( *, '(a)' ) ' '
  write ( *, '(a2,g14.6,a2,g14.6,a2,g14.6)' ) &
    '  ', weight_sum, '  ', weight_sum_exact, '  ', weight_sum_error

  deallocate ( grid_point )
  deallocate ( grid_weight )

  return
end
subroutine test05 ( dim_num, level_max, degree_max )

!****************************************************************************80
!
!! TEST05 tests a Gauss-Legendre sparse grid rule for monomial exactness.
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
!    04 July 2008
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
  integer ( kind = 4 ) t
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Check the exactness of a Gauss-Legendre sparse '
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
!  Determine the number of points in the rule.
!
  call sparse_grid_gl_size ( dim_num, level_max, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique points in the grid = ', point_num
!
!  Allocate space for the weights and points.
!
  allocate ( grid_weight(1:point_num) )
  allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  call sparse_grid_gl ( dim_num, level_max, point_num, grid_weight, &
    grid_point )
!
!  Rescale the weights, and translate the abscissas.
!
  volume = 2.0D+00**dim_num

  grid_weight(1:point_num) = grid_weight(1:point_num) / volume

  grid_point(1:dim_num,1:point_num) = &
    ( grid_point(1:dim_num,1:point_num) + 1.0D+00 ) / 2.0D+00
!
!  Explore the monomials.
!
  allocate ( expon(1:dim_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Error      Total   Monomial'
  write ( *, '(a)' ) '                 Degree  Exponents'
  write ( *, '(a)' ) ' '

  do degree = 0, degree_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( degree, dim_num, expon, more, h, t )

      call monomial_quadrature ( dim_num, expon, point_num, grid_weight, &
        grid_point, quad_error )

      write ( *, '(a,g14.6,a,i2,a,10i3)' ) &
        '  ', quad_error, '  ', degree, '    ', expon(1:dim_num)

      if ( .not. more ) then
        exit
      end if

    end do

    write ( *, '(a)' ) ' '

  end do

  deallocate ( expon )
  deallocate ( grid_point )
  deallocate ( grid_weight )

  return
end
subroutine test06 ( dim_num, level_max )

!****************************************************************************80
!
!! TEST06 creates a sparse Gauss-Legendre grid and writes it to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2007
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
  logical header
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_2 = 2
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  character ( len = 80 ) r_filename
  real ( kind = 8 ), allocatable, dimension ( : )    :: w
  character ( len = 80 ) w_filename
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  character ( len = 80 ) x_filename

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  Call SPARSE_GRID_GL to make a sparse Gauss-Legendre grid.'
  write ( *, '(a)' ) '  Write the data to a set of quadrature files.'

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
!
!  Determine the number of points.
!
  call sparse_grid_gl_size ( dim_num, level_max, point_num )
!
!  Allocate space for the weights and points.
!
  allocate ( r(1:dim_num,1:2) )
  allocate ( w(1:point_num) )
  allocate ( x(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  r(1:dim_num,1) = -1.0D+00
  r(1:dim_num,2) = +1.0D+00

  call sparse_grid_gl ( dim_num, level_max, point_num, w, x )
!
!  Write the data out.
!
  header = .false.

  write ( r_filename, '(a,i2,a,i3,a)' ) 'gl_d', dim_num, '_level', level_max, '_r.txt'
  write ( w_filename, '(a,i2,a,i3,a)' ) 'gl_d', dim_num, '_level', level_max, '_w.txt'
  write ( x_filename, '(a,i2,a,i3,a)' ) 'gl_d', dim_num, '_level', level_max, '_x.txt'

  call s_blank_delete ( r_filename )
  call s_blank_delete ( w_filename )
  call s_blank_delete ( x_filename )

  call r8mat_write ( r_filename, dim_num, 2,         r )
  call r8mat_write ( w_filename, 1,       point_num, w )
  call r8mat_write ( x_filename, dim_num, point_num, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R data written to "' // trim ( r_filename ) // '".'
  write ( *, '(a)' ) '  W data written to "' // trim ( w_filename ) // '".'
  write ( *, '(a)' ) '  X data written to "' // trim ( x_filename ) // '".'

  deallocate ( r )
  deallocate ( w )
  deallocate ( x )

  return
end
