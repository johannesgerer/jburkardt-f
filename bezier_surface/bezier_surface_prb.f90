program main

!*****************************************************************************80
!
!! MAIN is the main program for BEZIER_SURFACE_PRB.
!
!  Discussion:
!
!    BEZIER_SURFACE_PRB tests routines from the BEZIER_SURFACE library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BEZIER_SURFACE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BEZIER_SURFACE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BEZIER_SURFACE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests routines to read a Bezier surface definition.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 100 ) :: node_file_name = 'teapot_nodes.txt'
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz
  character ( len = 100 ) :: rectangle_file_name = 'teapot_rectangles.txt'
  integer ( kind = 4 ) rectangle_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: rectangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  BEZIER_SURFACE_NODE_SIZE determines the number of'
  write ( *, '(a)' ) '    nodes in a Bezier surface node file.'
  write ( *, '(a)' ) '  BEZIER_SURFACE_NODE_READ reads the'
  write ( *, '(a)' ) '    nodes in a Bezier surface node file.'
  write ( *, '(a)' ) '  BEZIER_SURFACE_RECTANGLES_SIZE determines the number of'
  write ( *, '(a)' ) '    rectangles in a Bezier surface rectangle file.'
  write ( *, '(a)' ) '  BEZIER_SURFACE_RECTANGLES_READ reads the'
  write ( *, '(a)' ) '    rectangles in a Bezier surface rectangle file.'
!
!  Get the number of nodes, allocate space for them, and read them in.
!
  call bezier_surface_node_size ( node_file_name, node_num )
  allocate ( node_xyz(1:3,1:node_num) )
  call bezier_surface_node_read ( node_file_name, node_num, node_xyz )
  call bezier_surface_node_print ( node_num, node_xyz )
!
!  Get the number of rectangles, allocate space for them, and read them in.
!
  call bezier_surface_rectangle_size ( rectangle_file_name, rectangle_num )
  allocate ( rectangle_node(1:16,1:rectangle_num) )
  call bezier_surface_rectangle_read ( rectangle_file_name, rectangle_num, &
    rectangle_node )
  call bezier_surface_rectangle_print ( rectangle_num, rectangle_node )
!
!  Free the memory.
!
  deallocate ( node_xyz )
  deallocate ( rectangle_node )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BEZIER_SURFACE_NEIGHBORS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 100 ) :: rectangle_file_name = 'teapot_rectangles.txt'
  integer ( kind = 4 ) rectangle_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: rectangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: rectangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  BEZIER_SURFACE_NEIGHBORS determines patch neighbors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Note that, for this example, the teapot, there are'
  write ( *, '(a)' ) '  cases where more than two patches meet at a'
  write ( *, '(a)' ) '  (degenerate) side.  This routine will not handle'
  write ( *, '(a)' ) '  such cases completely.'
!
!  Get the number of rectangles, allocate space for them, and read them in.
!
  call bezier_surface_rectangle_size ( rectangle_file_name, rectangle_num )
  allocate ( rectangle_node(1:16,1:rectangle_num) )
  call bezier_surface_rectangle_read ( rectangle_file_name, rectangle_num, &
    rectangle_node )
! call bezier_surface_rectangle_print ( rectangle_num, rectangle_node )
!
!  Compute and print the neighbor array.
!
  allocate ( rectangle_neighbor(1:4,1:rectangle_num) )
  call bezier_surface_neighbors ( rectangle_num, rectangle_node, &
    rectangle_neighbor )
  call i4mat_transpose_print ( 4, rectangle_num, rectangle_neighbor, &
    '  Bezier patch neighbors:' )
!
!  Free the memory.
!
  deallocate ( rectangle_neighbor )
  deallocate ( rectangle_node )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests BEZIER_PATCH_EVALUATE.
!
!  Discussion:
!
!    For simplicity, we set up a Bezier surface of a single patch.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 16
  integer ( kind = 4 ), parameter :: point_num = 16
  integer ( kind = 4 ), parameter :: rectangle_num = 1

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz(3,node_num)
  integer ( kind = 4 ) patch
  integer ( kind = 4 ) point
  real ( kind = 8 ) point_uv(2,point_num)
  real ( kind = 8 ) point_xyz(3,point_num)
  integer ( kind = 4 ), dimension(16,rectangle_num) :: rectangle_node = reshape ( (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 /), (/ 16, 1 /) )
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  BEZIER_PATCH_EVALUATE evaluates points in one'
  write ( *, '(a)' ) '  patch of a Bezier surface.'

  node = 0
  do j = 1, 4
    y = real ( j - 1, kind = 8 ) / 3.0D+00
    do i = 1, 4
      x = real ( i - 1, kind = 8 ) / 3.0D+00
      node = node + 1
      node_xyz(1,node) = x
      node_xyz(2,node) = y
      node_xyz(3,node) = x * ( 1.0D+00 - x ) * y * y
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nodal coordinates:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i2,2x,5g12.4)' ) node, node_xyz(1:3,node)
  end do

  patch = 1

  point = 0
  do j = 1, 4
    do i = 1, 4
      point = point + 1
      point_uv(1,point) = real ( i - 1, kind = 8 ) / 3.0D+00
      point_uv(2,point) = real ( j - 1, kind = 8 ) / 3.0D+00
    end do
  end do

  call bezier_patch_evaluate ( node_num, node_xyz, rectangle_num, &
    rectangle_node, patch, point_num, point_uv, point_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (U,V) --> (X,Y,Z) coordinates:'
  write ( *, '(a)' ) ' '

  do point = 1, point_num
    write ( *, '(2x,i2,2x,5g12.4)' ) point, point_uv(1:2,point), &
      point_xyz(1:3,point)
  end do

  return
end
