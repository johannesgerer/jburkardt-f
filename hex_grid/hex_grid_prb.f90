program main

!*****************************************************************************80
!
!! MAIN is the main program for HEX_GRID_PRB.
!
!  Discussion:
!
!    HEX_GRID_PRB calls a set of problems for HEX_GRID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEX_GRID_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HEX_GRID library.'
!
!  Tests for unit square.
!
  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Tests for arbitrary size coordinate box.
!
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
  call test11 ( )
  call test12 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEX_GRID_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests HEX_GRID_01_LAYERS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 17

  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: nodes_per_layer_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101, 1001, 10001 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a hexagonal grid of points in the unit square,'
  write ( *, '(a)' ) '  given NODES_PER_LAYER, the number of grid points'
  write ( *, '(a)' ) '  along the first layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_01_LAYERS computes LAYERS, the number of'
  write ( *, '(a)' ) '  layers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   NODES  LAYERS'
  write ( *, '(a)' ) '     PER'
  write ( *, '(a)' ) '   LAYER'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nodes_per_layer = nodes_per_layer_test ( test )
    call hex_grid_01_layers ( nodes_per_layer, layers )
    write ( *, '(2x,i6,2x,i6)' ) nodes_per_layer, layers
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HEX_GRID_01_H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 17

  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: nodes_per_layer_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101, 1001, 10001 /)
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For a hexagonal grid of points in the unit square,'
  write ( *, '(a)' ) '  given NODES_PER_LAYER, the number of grid points'
  write ( *, '(a)' ) '  along the first layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_01_H computes HX and HY, the spacings'
  write ( *, '(a)' ) '  in the row and column directions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    NODES    LAYERS   HX          HY'
  write ( *, '(a)' ) '      PER'
  write ( *, '(a)' ) '    LAYER'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nodes_per_layer = nodes_per_layer_test ( test )
    call hex_grid_01_layers ( nodes_per_layer, layers )
    call hex_grid_01_h ( nodes_per_layer, hx, hy )
    write ( *, '(2x,i6,2x,i6,2x,f10.6,2x,f10.6)' ) &
      nodes_per_layer, layers, hx, hy
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LAYERS is chosen so that LAYERS-1 layers just fit'
  write ( *, '(a)' ) '  inside the unit square, but LAYERS layers do not.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LAYERS      HY     (LAYERS-1)*HY    LAYERS*HY'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nodes_per_layer = nodes_per_layer_test ( test )
    call hex_grid_01_layers ( nodes_per_layer, layers )
    call hex_grid_01_h ( nodes_per_layer, hx, hy )

    temp1 = real ( layers - 1, kind = 8 ) * hy
    temp2 = real ( layers, kind = 8 ) * hy

    write ( *, '(2x,i6,2x,f10.6,2x,f10.6,2x,f10.6)' ) &
      layers, hy, temp1, temp2
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests HEX_GRID_01_N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 17

  integer ( kind = 4 ) n
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: nodes_per_layer_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101, 1001, 10001 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a hexagonal grid of points in the unit square,'
  write ( *, '(a)' ) '  given NODES_PER_LAYER, the number of grid points'
  write ( *, '(a)' ) '  along the first layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_01_N computes N, the total number of grid'
  write ( *, '(a)' ) '  points in the unit square.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    NODES   LAYERS           N'
  write ( *, '(a)' ) '      PER'
  write ( *, '(a)' ) '    LAYER'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    nodes_per_layer = nodes_per_layer_test ( test )

    call hex_grid_01_layers ( nodes_per_layer, layers )

    call hex_grid_01_n ( nodes_per_layer, n )

    write ( *, '(2x,i6,2x,i6,2x,i12)' ) nodes_per_layer, layers, n

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests HEX_GRID_01_POINTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 15

  real ( kind = 8 ) box(dim_num,2)
  character ( len = 80 ) file_name
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) n
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: nodes_per_layer_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 /)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a hexagonal grid of points in the unit square,'
  write ( *, '(a)' ) '  given NODES_PER_LAYER, the number of grid points'
  write ( *, '(a)' ) '  along the first layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_01_POINTS computes P, the coordinates'
  write ( *, '(a)' ) '  of the points of the grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_WRITE writes the data to a file.'

  box(1:dim_num,1:2) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00 /), (/ dim_num, 2 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   NODES  LAYERS    N    Filename'
  write ( *, '(a)' ) '     PER'
  write ( *, '(a)' ) '   LAYER'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    nodes_per_layer = nodes_per_layer_test ( test )
    call hex_grid_01_layers ( nodes_per_layer, layers )
    call hex_grid_01_h ( nodes_per_layer, hx, hy )
    call hex_grid_01_n ( nodes_per_layer, n )

    allocate ( p(1:dim_num,1:n) )

    call hex_grid_01_points ( nodes_per_layer, layers, n, p )

    write ( file_name, '(a,i3,a,i3,a,i6,a)' ) &
      'hex_grid_01_', nodes_per_layer, '_', layers, '_', n, '.txt'

    call s_blank_delete ( file_name )

    write ( *, '(2x,i6,2x,i3,2x,i6,4x,a)' ) &
      nodes_per_layer, layers, n, trim ( file_name )

    call hex_grid_write ( n, nodes_per_layer, layers, hx, hy, box, &
      p, file_name )

    deallocate ( p )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests HEX_GRID_01_APPROXIMATE_N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_goal
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: n_goal_test = (/ &
    100, 200, 500, 10000 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For a hexagonal grid of points in the unit box,'
  write ( *, '(a)' ) '  HEX_GRID_01_APPROXIMATE_N seeks the value of'
  write ( *, '(a)' ) '  NODES_PER_LAYER that produces a mesh of about N points.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N_GOAL  NODES_PER_LAYER       N'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    n_goal = n_goal_test ( test )

    call hex_grid_01_approximate_n ( n_goal, nodes_per_layer, n )

    write ( *, '(2x,i6,2x,9x,i6,2x,i6)' ) n_goal, nodes_per_layer, n

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests HEX_GRID_01_APPROXIMATE_H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) h
  real ( kind = 8 ) h_goal
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nodes_per_layer
  real ( kind = 8 ), dimension (test_num) :: h_goal_test = (/ &
    0.10D+00, 0.01D+00, 0.0001D+00 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For a hexagonal grid of points in the unit box,'
  write ( *, '(a)' ) '  HEX_GRID_01_APPROXIMATE_H seeks the value of'
  write ( *, '(a)' ) '  NODES_PER_LAYER that produces a mesh with spacing'
  write ( *, '(a)' ) '  that is H or less.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      H_GOAL      NODES_PER_LAYER      H                      N'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    h_goal = h_goal_test ( test )

    call hex_grid_01_approximate_h ( h_goal, nodes_per_layer, h )

    call hex_grid_01_n ( nodes_per_layer, n )

    write ( *, '(2x,g14.6,2x,9x,i6,2x,g14.6,2x,i12)' ) &
      h_goal, nodes_per_layer, h, n

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests HEX_GRID_LAYERS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 15

  real ( kind = 8 ) box(dim_num,2)
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: nodes_per_layer_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For a hexagonal grid of points in a coordinate box,'
  write ( *, '(a)' ) '  given NODES_PER_LAYER, the number of grid points'
  write ( *, '(a)' ) '  along the first layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_LAYERS computes LAYERS, the number of'
  write ( *, '(a)' ) '  layers.'

  box(1:dim_num,1:2) = reshape ( (/ &
    1.0D+00, 2.0D+00, &
    4.0D+00, 7.0D+00 /), (/ dim_num, 2 /) )

  call box_print_2d ( box )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   NODES  LAYERS'
  write ( *, '(a)' ) '     PER'
  write ( *, '(a)' ) '   LAYER'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nodes_per_layer = nodes_per_layer_test ( test )
    call hex_grid_layers ( nodes_per_layer, box, layers )
    write ( *, '(2x,i6,2x,i6)' ) nodes_per_layer, layers
  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests HEX_GRID_H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 15

  real ( kind = 8 ) box(dim_num,2)
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: nodes_per_layer_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 /)
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  For a hexagonal grid of points in a coordinate box,'
  write ( *, '(a)' ) '  given NODES_PER_LAYER, the number of grid points'
  write ( *, '(a)' ) '  along the first layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_H computes HX and HY, the spacings'
  write ( *, '(a)' ) '  in the row and column directions.'

  box(1:dim_num,1:2) = reshape ( (/ &
    1.0D+00, 2.0D+00, &
    4.0D+00, 7.0D+00 /), (/ dim_num, 2 /) )

  call box_print_2d ( box )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    NODES    LAYERS   HX          HY'
  write ( *, '(a)' ) '      PER'
  write ( *, '(a)' ) '    LAYER'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nodes_per_layer = nodes_per_layer_test ( test )
    call hex_grid_layers ( nodes_per_layer, box, layers )
    call hex_grid_h ( nodes_per_layer, box, hx, hy )
    write ( *, '(2x,i6,2x,i6,2x,f10.6,2x,f10.6)' ) &
      nodes_per_layer, layers, hx, hy
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LAYERS is chosen so that LAYERS-1 layers just fit'
  write ( *, '(a)' ) '  inside the unit square, but LAYERS layers do not.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LAYERS      HY     (LAYERS-1)*HY    LAYERS*HY'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nodes_per_layer = nodes_per_layer_test ( test )
    call hex_grid_layers ( nodes_per_layer, box, layers )
    call hex_grid_h ( nodes_per_layer, box, hx, hy )

    temp1 = real ( layers - 1, kind = 8 ) * hy
    temp2 = real ( layers, kind = 8 ) * hy

    write ( *, '(2x,i6,2x,f10.6,2x,f10.6,2x,f10.6)' ) &
      layers, hy, temp1, temp2
  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests HEX_GRID_N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 15

  real ( kind = 8 ) box(dim_num,2)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: nodes_per_layer_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For a hexagonal grid of points in a coordinate box,'
  write ( *, '(a)' ) '  given NODES_PER_LAYER, the number of grid points'
  write ( *, '(a)' ) '  along the first layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_N computes N, the total number of grid'
  write ( *, '(a)' ) '  points in the coordinate box.'

  box(1:dim_num,1:2) = reshape ( (/ &
    1.0D+00, 2.0D+00, &
    4.0D+00, 7.0D+00 /), (/ dim_num, 2 /) )

  call box_print_2d ( box )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    NODES   LAYERS           N'
  write ( *, '(a)' ) '      PER'
  write ( *, '(a)' ) '    LAYER'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    nodes_per_layer = nodes_per_layer_test ( test )

    call hex_grid_layers ( nodes_per_layer, box, layers )

    call hex_grid_n ( nodes_per_layer, box, n )

    write ( *, '(2x,i6,2x,i6,2x,i12)' ) nodes_per_layer, layers, n

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests HEX_GRID_POINTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 12

  real ( kind = 8 ) box(dim_num,2)
  character ( len = 80 ) file_name
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) n
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: nodes_per_layer_test = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21 /)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For a hexagonal grid of points in a coordinate box,'
  write ( *, '(a)' ) '  given NODES_PER_LAYER, the number of grid points'
  write ( *, '(a)' ) '  along the first layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_POINTS computes P, the coordinates'
  write ( *, '(a)' ) '  of the points of the grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HEX_GRID_WRITE writes the data to a file.'

  box(1:dim_num,1:2) = reshape ( (/ &
    1.0D+00, 2.0D+00, &
    4.0D+00, 7.0D+00 /), (/ dim_num, 2 /) )

  call box_print_2d ( box )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   NODES  LAYERS    N    Filename'
  write ( *, '(a)' ) '     PER'
  write ( *, '(a)' ) '   LAYER'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    nodes_per_layer = nodes_per_layer_test ( test )
    call hex_grid_layers ( nodes_per_layer, box, layers )
    call hex_grid_h ( nodes_per_layer, box, hx, hy )
    call hex_grid_n ( nodes_per_layer, box, n )

    allocate ( p(1:dim_num,1:n) )

    call hex_grid_points ( nodes_per_layer, layers, n, box, p )

    write ( file_name, '(a,i3,a,i3,a,i6,a)' ) &
      'hex_grid_', nodes_per_layer, '_', layers, '_', n, '.txt'

    call s_blank_delete ( file_name )

    write ( *, '(2x,i6,2x,i3,2x,i6,4x,a)' ) &
      nodes_per_layer, layers, n, trim ( file_name )

    call hex_grid_write ( n, nodes_per_layer, layers, hx, hy, box, &
      p, file_name )

    deallocate ( p )

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests HEX_GRID_APPROXIMATE_N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) box(dim_num,2)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_goal
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ), dimension (test_num) :: n_goal_test = (/ &
    100, 200, 500, 10000 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  For a hexagonal grid of points in a coordinate box,'
  write ( *, '(a)' ) '  HEX_GRID_APPROXIMATE_N seeks the value of'
  write ( *, '(a)' ) '  NODES_PER_LAYER that produces a mesh of about N points.'

  box(1:dim_num,1:2) = reshape ( (/ &
    1.0D+00, 2.0D+00, &
    4.0D+00, 7.0D+00 /), (/ dim_num, 2 /) )

  call box_print_2d ( box )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N_GOAL  NODES_PER_LAYER       N'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    n_goal = n_goal_test ( test )

    call hex_grid_approximate_n ( box, n_goal, nodes_per_layer, n )

    write ( *, '(2x,i6,2x,9x,i6,2x,i6)' ) n_goal, nodes_per_layer, n

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests HEX_GRID_APPROXIMATE_H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) box(dim_num,2)
  real ( kind = 8 ) h
  real ( kind = 8 ) h_goal
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nodes_per_layer
  real ( kind = 8 ), dimension (test_num) :: h_goal_test = (/ &
    0.10D+00, 0.01D+00, 0.0001D+00 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For a hexagonal grid of points in the unit box,'
  write ( *, '(a)' ) '  HEX_GRID_APPROXIMATE_H seeks the value of'
  write ( *, '(a)' ) '  NODES_PER_LAYER that produces a mesh with spacing'
  write ( *, '(a)' ) '  that is H or less.'
  write ( *, '(a)' ) ' '

  box(1:dim_num,1:2) = reshape ( (/ &
    1.0D+00, 2.0D+00, &
    4.0D+00, 7.0D+00 /), (/ dim_num, 2 /) )

  call box_print_2d ( box )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      H_GOAL      NODES_PER_LAYER      H                      N'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    h_goal = h_goal_test ( test )

    call hex_grid_approximate_h ( box, h_goal, nodes_per_layer, h )

    call hex_grid_n ( nodes_per_layer, box, n )

    write ( *, '(2x,g14.6,2x,9x,i6,2x,g14.6,2x,i12)' ) &
      h_goal, nodes_per_layer, h, n

  end do

  return
end
