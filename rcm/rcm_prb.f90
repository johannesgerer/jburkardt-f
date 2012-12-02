program main

!*****************************************************************************80
!
!! MAIN is the main program for RCM_PRB.
!
!  Discussion:
!
!    RCM_PRB runs the RCM tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RCM_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the RCM library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
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
  write ( *, '(a)' ) 'RCM_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ADJ_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 10
  integer ( kind = 4 ), parameter :: adj_max = node_num * ( node_num - 1 )

  integer ( kind = 4 ) adj(adj_max)
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n_calls
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  ADJ_SET sets up an adjacency matrix incrementally.'

  n_calls = i4_uniform ( 1, adj_max, seed )

  call adj_set ( node_num, adj_max, adj_num, adj_row, adj, -1, -1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Creating and recording adjacency information:'
  write ( *, '(a)' ) ' '

  do k = 1, n_calls

    i = i4_uniform ( 1, node_num, seed )
    j = i4_uniform ( 1, node_num, seed )

    write ( *, '(2x,i8,2x,i8)' ) i, j

    call adj_set ( node_num, adj_max, adj_num, adj_row, adj, i, j )

  end do

  call adj_print ( node_num, adj_num, adj_row, adj, &
    '  Random adjacency matrix:' )

  call adj_show ( node_num, adj_num, adj_row, adj )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests GENRCM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_bandwidth
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_perm_bandwidth
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) i
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm_inv

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  GENRCM reorders the nodes in a graph using'
  write ( *, '(a)' ) '  the Reverse Cuthill McKee algorithm.'

  call graph_01_size ( node_num, adj_num )

  allocate ( adj_row(node_num+1) )
  allocate ( adj(adj_num) )
  allocate ( perm(node_num) )
  allocate ( perm_inv(node_num) )

  call graph_01_adj ( node_num, adj_num, adj_row, adj )

  call adj_print ( node_num, adj_num, adj_row, adj, '  Adjacency matrix:' )

  call adj_show ( node_num, adj_num, adj_row, adj )

  bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ADJ bandwidth = ', bandwidth

  call genrcm ( node_num, adj_num, adj_row, adj, perm )

  call perm_inverse3 ( node_num, perm, perm_inv )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The RCM permutation and inverse:'
  write ( *, '(a)' ) ' '

  do i = 1, node_num
    write ( *, '(2x,3i8)' ) i, perm(i), perm_inv(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Permuted adjacency matrix:'
  write ( *, '(a)' ) ' '

  call adj_perm_show ( node_num, adj_num, adj_row, adj, perm, perm_inv )

  bandwidth = adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, &
    perm, perm_inv )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ADJ (permuted) bandwidth = ', bandwidth

  deallocate ( adj )
  deallocate ( adj_row )
  deallocate ( perm )
  deallocate ( perm_inv )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests GENRCM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_bandwidth
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_perm_bandwidth
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm_inv
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  GENRCM generates the Reverse Cuthill McKee ordering.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do the test twice.  On the second test, randomly'
  write ( *, '(a)' ) '  permute the initial nodes.'

  call triangulation_order3_example2_size ( node_num, triangle_num, hole_num )

  do test = 1, 2

    allocate ( node_xy(1:2,1:node_num) )
    allocate ( triangle_node(1:3,1:triangle_num) )
    allocate ( triangle_neighbor(1:3,1:triangle_num) )

    call triangulation_order3_example2 ( node_num, triangle_num, node_xy, &
      triangle_node, triangle_neighbor )
!
!  Randomly permute the nodes.
!
    if ( test == 2 ) then

      seed = 123456789

      allocate ( perm(1:node_num) )

      call perm_uniform ( node_num, seed, perm )

      call i4vec_print ( node_num, perm, '  The random permutation:' )

      do i = 1, 3
        do j = 1, triangle_num
          node = triangle_node(i,j)
          triangle_node(i,j) = perm ( node )
        end do
      end do

      deallocate ( perm )

    end if

    call i4mat_transpose_print ( 3, triangle_num, triangle_node, &
      '  TRIANGLE_NODE:' )

    allocate ( adj_row(1:node_num+1) )

    call triangulation_order3_adj_count ( node_num, triangle_num, &
      triangle_node, triangle_neighbor, adj_num, adj_row )

    allocate ( adj(1:adj_num) )

    call triangulation_order3_adj_set ( node_num, triangle_num, triangle_node, &
      triangle_neighbor, adj_num, adj_row, adj )

    call adj_print ( node_num, adj_num, adj_row, adj, '  ADJ array:' )

    bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ADJ bandwidth = ', bandwidth

    allocate ( perm(1:node_num) )

    call genrcm ( node_num, adj_num, adj_row, adj, perm )

    call i4vec_print ( node_num, perm, '  The RCM permutation:' )

    allocate ( perm_inv(1:node_num) )

    call perm_inverse3 ( node_num, perm, perm_inv )

    bandwidth = adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, &
      perm, perm_inv )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Permuted ADJ bandwidth = ', bandwidth

    deallocate ( adj )
    deallocate ( adj_row )
    deallocate ( node_xy )
    deallocate ( perm )
    deallocate ( perm_inv )
    deallocate ( triangle_neighbor )
    deallocate ( triangle_node )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests GENRCM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_bandwidth
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_perm_bandwidth
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm_inv
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node
  integer ( kind = 4 ), parameter :: triangle_order = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  GENRCM generates the Reverse Cuthill McKee ordering.'

  call triangulation_order6_example2_size ( node_num, triangle_num, hole_num )

  allocate ( node_xy(1:2,1:node_num) )
  allocate ( triangle_node(1:triangle_order,1:triangle_num) )
  allocate ( triangle_neighbor(1:3,1:triangle_num) )

  call triangulation_order6_example2 ( node_num, triangle_num, node_xy, &
    triangle_node, triangle_neighbor )
!
!  Randomly permute the nodes.
!
  seed = 123456789

  allocate ( perm(1:node_num) )

  call perm_uniform ( node_num, seed, perm )

  call i4vec_print ( node_num, perm, '  The random permutation:' )

  do i = 1, triangle_order
    do j = 1, triangle_num
      node = triangle_node(i,j)
      triangle_node(i,j) = perm ( node )
    end do
  end do

  call i4mat_transpose_print ( triangle_order, triangle_num, triangle_node, &
    '  Permuted TRIANGLE_NODE' )

  deallocate ( perm )

  allocate ( adj_row(1:node_num+1) )

  call triangulation_order6_adj_count ( node_num, triangle_num, triangle_node, &
    triangle_neighbor, adj_num, adj_row )

  allocate ( adj(1:adj_num) )

  call triangulation_order6_adj_set ( node_num, triangle_num, triangle_node, &
    triangle_neighbor, adj_num, adj_row, adj )

  call adj_print ( node_num, adj_num, adj_row, adj, '  ADJ array:' )

  bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ADJ bandwidth = ', bandwidth

  allocate ( perm(1:node_num) )

  call genrcm ( node_num, adj_num, adj_row, adj, perm )

  call i4vec_print ( node_num, perm, '  The RCM permutation:' )

  allocate ( perm_inv(1:node_num) )

  call perm_inverse3 ( node_num, perm, perm_inv )

  bandwidth = adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, &
    perm, perm_inv )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Permuted ADJ bandwidth = ', bandwidth

  deallocate ( adj )
  deallocate ( adj_row )
  deallocate ( node_xy )
  deallocate ( perm )
  deallocate ( perm_inv )
  deallocate ( triangle_neighbor )
  deallocate ( triangle_node )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests GRAPH_01_ADJ and GRAPH_01_SIZE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  GRAPH_01_SIZE returns the sizes for graph 1.'
  write ( *, '(a)' ) '  GRAPH_01_ADJ returns the adjacency for graph 1.'
  write ( *, '(a)' ) '  ADJ_PRINT prints the adjacency information.'

  call graph_01_size ( node_num, adj_num )

  allocate ( adj_row(node_num+1) )
  allocate ( adj(adj_num) )

  call graph_01_adj ( node_num, adj_num, adj_row, adj )

  call adj_print ( node_num, adj_num, adj_row, adj, &
    '  Adjacency for GRAPH_01:' )

  call adj_show ( node_num, adj_num, adj_row, adj )

  deallocate ( adj )
  deallocate ( adj_row )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests LEVEL_SET;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ), allocatable, dimension ( : ) :: level
  integer ( kind = 4 ) level_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: level_row
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mask
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) root
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  LEVEL_SET computes the level sets of a graph,'
  write ( *, '(a)' ) '  given a root node (which defines level 1).'

  call graph_01_size ( node_num, adj_num )

  allocate ( adj_row(node_num+1) )
  allocate ( adj(adj_num) )

  call graph_01_adj ( node_num, adj_num, adj_row, adj )

  call adj_print ( node_num, adj_num, adj_row, adj, '  Adjacency matrix:' )

  call adj_show ( node_num, adj_num, adj_row, adj )
!
!  Choose different roots.
!
  allocate ( level(node_num) )
  allocate ( level_row(node_num+1) )
  allocate ( mask(node_num) )

  do i = 1, 3

    root = i4_uniform ( 1, node_num, seed )

    mask(1:node_num) = 1

    call level_set ( root, adj_num, adj_row, adj, mask, level_num, &
      level_row, level, node_num )

    call level_set_print ( node_num, level_num, level_row, level )

  end do

  deallocate ( adj )
  deallocate ( adj_row )
  deallocate ( level )
  deallocate ( level_row )
  deallocate ( mask )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests ROOT_FIND;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable, dimension ( : ) :: level
  integer ( kind = 4 ) level_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: level_row
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mask
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) root
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  ROOT_FIND is given a node in the graph,'
  write ( *, '(a)' ) '  and returns a better node to use as a starting'
  write ( *, '(a)' ) '  point for reordering.'

  call graph_01_size ( node_num, adj_num )

  allocate ( adj_row(node_num+1) )
  allocate ( adj(adj_num) )

  call graph_01_adj ( node_num, adj_num, adj_row, adj )

  call adj_print ( node_num, adj_num, adj_row, adj, '  Adjacency matrix:' )

  call adj_show ( node_num, adj_num, adj_row, adj )

  allocate ( level(node_num) )
  allocate ( level_row(node_num+1) )
  allocate ( mask(node_num) )

  do i = 1, node_num

    root = i

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Starting root =    ', root

    mask(1:node_num) = 1

    call root_find ( root, adj_num, adj_row, adj, mask, level_num, &
      level_row, level, node_num )

    write ( *, '(a,i8)' ) '  Suggested root =   ', root
    write ( *, '(a,i8)' ) '  Number of levels = ', level_num

  end do

  deallocate ( adj )
  deallocate ( adj_row )
  deallocate ( level )
  deallocate ( level_row )
  deallocate ( mask )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests TRIANGULATION_ORDER3_ADJ_COUNT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_ADJ_COUNT counts the (lower)'
  write ( *, '(a)' ) '  adjacencies defined by a triangulation.'

  call triangulation_order3_example2_size ( node_num, triangle_num, hole_num )

  allocate ( node_xy(1:2,1:node_num) )
  allocate ( triangle_node(1:3,1:triangle_num) )
  allocate ( triangle_neighbor(1:3,1:triangle_num) )

  call triangulation_order3_example2 ( node_num, triangle_num, node_xy, &
    triangle_node, triangle_neighbor )

  call i4mat_transpose_print ( 3, triangle_num, triangle_node, &
    '  TRIANGLE_NODE' )

  allocate ( adj_row(1:node_num+1) )

  call triangulation_order3_adj_count ( node_num, triangle_num, triangle_node, &
    triangle_neighbor, adj_num, adj_row )

  call i4vec_print ( node_num+1, adj_row, '  ADJ_ROW' )

  deallocate ( adj_row )
  deallocate ( node_xy )
  deallocate ( triangle_neighbor )
  deallocate ( triangle_node )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests TRIANGULATION_ORDER3_ADJ_SET
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_bandwidth
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_ADJ_SET sets the (lower)'
  write ( *, '(a)' ) '  adjacencies defined by a triangulation.'

  call triangulation_order3_example2_size ( node_num, triangle_num, hole_num )

  allocate ( node_xy(1:2,1:node_num) )
  allocate ( triangle_node(1:3,1:triangle_num) )
  allocate ( triangle_neighbor(1:3,1:triangle_num) )

  call triangulation_order3_example2 ( node_num, triangle_num, node_xy, &
    triangle_node, triangle_neighbor )

  call i4mat_transpose_print ( 3, triangle_num, triangle_node, &
    '  TRIANGLE_NODE' )

  allocate ( adj_row(1:node_num+1) )

  call triangulation_order3_adj_count ( node_num, triangle_num, triangle_node, &
    triangle_neighbor, adj_num, adj_row )

  allocate ( adj(1:adj_num) )

  call triangulation_order3_adj_set ( node_num, triangle_num, triangle_node, &
    triangle_neighbor, adj_num, adj_row, adj )

  call adj_print ( node_num, adj_num, adj_row, adj, '  ADJ array:' )

  bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ADJ bandwidth = ', bandwidth

  deallocate ( adj )
  deallocate ( adj_row )
  deallocate ( node_xy )
  deallocate ( triangle_neighbor )
  deallocate ( triangle_node )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests TRIANGULATION_NEIGHBOR_TRIANGLES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: triangle_num = 16
  integer ( kind = 4 ), parameter :: triangle_order = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension (triangle_order,triangle_num) :: triangle_node = reshape ( (/ &
     3,   4,   1, &
     3,   1,   2, &
     3,   2,   8, &
     2,   1,   5, &
     8,   2,  13, &
     8,  13,   9, &
     3,   8,   9, &
    13,   2,   5, &
     9,  13,   7, &
     7,  13,   5, &
     6,   7,   5, &
     9,   7,   6, &
    10,   9,   6, &
     6,   5,  12, &
    11,   6,  12, &
    10,   6,  11 /), (/ triangle_order, triangle_num /) )
  integer ( kind = 4 ), dimension (3,triangle_num) :: triangle_neighbor

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For a triangulation of a set of nodes,'
  write ( *, '(a)' ) '  TRIANGULATION_NEIGHBOR_TRIANGLES determines the'
  write ( *, '(a)' ) '    adjacency relationships between triangles.'

  call i4mat_transpose_print ( triangle_order, triangle_num, triangle_node, &
    '  Triangles:' )

  call triangulation_neighbor_triangles ( triangle_order, triangle_num, &
    triangle_node, triangle_neighbor )

  call i4mat_transpose_print ( 3, triangle_num, triangle_neighbor, &
    '  Triangle neighbors:' )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests TRIANGULATION_ORDER6_ADJ_COUNT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_ADJ_COUNT counts the (lower)'
  write ( *, '(a)' ) '  adjacencies defined by a triangulation.'

  call triangulation_order6_example2_size ( node_num, triangle_num, hole_num )

  allocate ( node_xy(1:2,1:node_num) )
  allocate ( triangle_node(1:6,1:triangle_num) )
  allocate ( triangle_neighbor(1:3,1:triangle_num) )

  call triangulation_order6_example2 ( node_num, triangle_num, node_xy, &
    triangle_node, triangle_neighbor )

  allocate ( adj_row(1:node_num+1) )

  call triangulation_order6_adj_count ( node_num, triangle_num, triangle_node, &
    triangle_neighbor, adj_num, adj_row )

  call i4vec_print ( node_num+1, adj_row, '  ADJ_ROW' )

  deallocate ( adj_row )
  deallocate ( node_xy )
  deallocate ( triangle_neighbor )
  deallocate ( triangle_node )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests TRIANGULATION_ORDER6_ADJ_SET
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_bandwidth
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_ADJ_SET sets the (lower)'
  write ( *, '(a)' ) '  adjacencies defined by a triangulation.'

  call triangulation_order6_example2_size ( node_num, triangle_num, hole_num )

  allocate ( node_xy(1:2,1:node_num) )
  allocate ( triangle_node(1:6,1:triangle_num) )
  allocate ( triangle_neighbor(1:3,1:triangle_num) )

  call triangulation_order6_example2 ( node_num, triangle_num, node_xy, &
    triangle_node, triangle_neighbor )

  call i4mat_transpose_print ( 6, triangle_num, triangle_node, &
    '  TRIANGLE_NODE' )

  allocate ( adj_row(1:node_num+1) )

  call triangulation_order6_adj_count ( node_num, triangle_num, triangle_node, &
    triangle_neighbor, adj_num, adj_row )

  allocate ( adj(1:adj_num) )

  call triangulation_order6_adj_set ( node_num, triangle_num, triangle_node, &
    triangle_neighbor, adj_num, adj_row, adj )

  call adj_print ( node_num, adj_num, adj_row, adj, '  ADJ array:' )

  bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ADJ bandwidth = ', bandwidth

  deallocate ( adj )
  deallocate ( adj_row )
  deallocate ( node_xy )
  deallocate ( triangle_neighbor )
  deallocate ( triangle_node )

  return
end
