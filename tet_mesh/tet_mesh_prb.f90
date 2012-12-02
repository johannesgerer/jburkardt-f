program main

!*****************************************************************************80
!
!! MAIN is the main program for TET_MESH_PRB.
!
!  Discussion:
!
!    TET_MESH_PRB tests the routines in TET_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TET_MESH_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TET_MESH library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test004 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TET_MESH_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests R8MAT_SOLVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: rhs_num = 2
!
!  Each row of this definition is a COLUMN of the matrix.
!
  real ( kind = 8 ), dimension (n,n+rhs_num) :: a = reshape ( &
    (/ 1.0D+00,  4.0D+00,  7.0D+00, &
       2.0D+00,  5.0D+00,  8.0D+00, &
       3.0D+00,  6.0D+00,  0.0D+00, &
      14.0D+00, 32.0D+00, 23.0D+00, &
       7.0D+00, 16.0D+00,  7.0D+00 /), &
    (/ n, n+rhs_num /) )
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  R8MAT_SOLVE solves linear systems.'
!
!  Print out the matrix to be inverted.
!
  call r8mat_print ( n, n+rhs_num, a, '  The linear system:' )
!
!  Solve the systems.
!
  call r8mat_solve ( n, rhs_num, a, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input matrix was singular.'
    write ( *, '(a)' ) '  The solutions could not be computed.'
    write ( *, '(a)' ) ' '
    return
  end if

  call r8mat_print ( n, rhs_num, a(1:n,n+1:n+rhs_num), &
    '  The computed solutions' )

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) phy(3,n)
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) ref(3,n)
  real ( kind = 8 ) ref2(3,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(3,4) :: tet = reshape ( (/ &
    5.0D+00, 0.0D+00, 0.0D+00, &
    8.0D+00, 0.0D+00, 0.0D+00, &
    5.0D+00, 2.0D+00, 0.0D+00, &
    6.0D+00, 1.0D+00, 2.0D+00 /), (/ 3, 4 /) )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  For an order 4 tetrahedron,'
  write ( *, '(a)' ) '  TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE '
  write ( *, '(a)' ) '    maps a physical point to a reference point.'
  write ( *, '(a)' ) '  TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL '
  write ( *, '(a)' ) '    maps a reference point to a physical point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     ( R, S, T )          ==>  ( X, Y, Z )           ==> ( R2, S2, T2 )'
  write ( *, '(a)' ) ' '

  call tetrahedron_reference_sample ( n, seed, ref )

  call tetrahedron_order4_reference_to_physical ( tet, n, ref, phy )
  call tetrahedron_order4_physical_to_reference ( tet, n, phy, ref2 )

  do j = 1, n

    write ( *, '(2x,3f8.4,2x,3f8.4,2x,3f8.4)' ) &
      ref(1:3,j), phy(1:3,j), ref2(1:3,j)

  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests TETRAHEDRON_ORDER10_TO_ORDER4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tet_node1
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tet_node2
  integer ( kind = 4 ) tet_num1
  integer ( kind = 4 ) tet_num2
  integer ( kind = 4 ), parameter :: tet_order1 = 10
  integer ( kind = 4 ), parameter :: tet_order2 = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  For an order 10 tet mesh,'
  write ( *, '(a)' ) '  TETRAHEDRON_ORDER10_TO_ORDER4 '
  write ( *, '(a)' ) '    makes a linear (order 4) tet mesh by using'
  write ( *, '(a)' ) '    the existing nodes, and replacing each'
  write ( *, '(a)' ) '    quadratic tetrahedron by 8 linear tetrahedrons.'

  call tet_mesh_order10_example_size ( node_num1, tet_num1 )

  allocate ( node_xyz(1:3,1:node_num1) )
  allocate ( tet_node1(1:tet_order1,1:tet_num1) )

  call tet_mesh_order10_example_set ( node_num1, tet_num1, &
    node_xyz, tet_node1 )

  call i4mat_transpose_print_some ( tet_order1, tet_num1, tet_node1, &
    1, 1, tet_order1, 5, '  First 5 quadratic tetrahedrons:' )

  call tet_mesh_order10_to_order4_size ( node_num1, tet_num1, &
    node_num2, tet_num2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Quadratic mesh size is       ', tet_num1
  write ( *, '(a,i8)' ) '  Linearized mesh size will be ', tet_num2

  allocate ( tet_node2(1:tet_order2,1:tet_num2) )

  call tet_mesh_order10_to_order4_compute ( tet_num1, tet_node1, &
    tet_num2, tet_node2 )

  call i4mat_transpose_print_some ( tet_order2, tet_num2, tet_node2, &
    1, 1, tet_order2, 5, '  First 5 linear tetrahedrons:' )

  deallocate ( node_xyz )
  deallocate ( tet_node1 )
  deallocate ( tet_node2 )

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests TET_MESH_NODE_ORDER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: node_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tet_node
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ), parameter :: tet_order = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  TET_MESH_NODE_ORDER determines the order of '
  write ( *, '(a)' ) '  each node in a tet mesh.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The order of a node is the number of tetrahedrons'
  write ( *, '(a)' ) '  that use the node as part of their definition.'

  call tet_mesh_order10_example_size ( node_num, tet_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  This mesh has tetrahedron order ', tet_order
  write ( *, '(a,i8)' ) '  The number of tetrahedrons is   ', tet_num

  allocate ( node_xyz(1:3,1:node_num) )
  allocate ( tet_node(1:tet_order,1:tet_num) )

  call tet_mesh_order10_example_set ( node_num, tet_num, &
    node_xyz, tet_node )

  call i4mat_transpose_print ( tet_order, tet_num, tet_node, &
    '  The tet mesh:' )

  allocate ( node_order(1:node_num) )

  call tet_mesh_node_order ( tet_order, tet_num, tet_node, &
    node_num, node_order )

  call i4vec_print ( node_num, node_order, '  Node orders:' );

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Check that the following are equal:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  Number of tetrahedrons * order = ', tet_num * tet_order
  write ( *, '(a,i8)' ) &
    '  Sum of node orders             = ', sum ( node_order(1:node_num) )

  deallocate ( node_xyz )
  deallocate ( node_order )
  deallocate ( tet_node )

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests TETRAHEDRON_BARYCENTRIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c1(4)
  real ( kind = 8 ) c1_sum
  real ( kind = 8 ) c2(4)
  real ( kind = 8 ) p(3)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test1
  integer ( kind = 4 ), parameter :: test1_num = 3
  integer ( kind = 4 ) test2
  integer ( kind = 4 ), parameter :: test2_num = 5
  real ( kind = 8 ) tet_xyz(3,4)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  TETRAHEDRON_BARYCENTRIC computes the barycentric'
  write ( *, '(a)' ) '  coordinates of a point.'
!
!  Choose a random tetrahedron.
!
  do test1 = 1, test1_num

    call r8mat_uniform_01 ( 3, 4, seed, tet_xyz )

    call r8mat_transpose_print ( 3, 4, tet_xyz, '  Random tetrahedron:' )
!
!  Choose barycentric coordinates C1 at random.
!
!  Define a point P.
!
!  Have TETRAHEDRON_BARYCENTRIC compute C2, the barycentric coordinates of P.
!
    do test2 = 1, test2_num

      call r8vec_uniform_01 ( 4, seed, c1 )
      c1_sum = sum ( c1(1:4) )
      c1(1:4) = c1(1:4) / c1_sum

      p(1:3) = matmul ( tet_xyz(1:3,1:4), c1(1:4) )

      call tetrahedron_barycentric ( tet_xyz, p, c2 )

      write ( *, '(a)' ) ' '
      write ( *, '(2x,a,4(2x,g14.6))' ) 'C1 = ', c1(1:4)
      write ( *, '(2x,a,4(2x,g14.6))' ) 'C2 = ', c2(1:4)

    end do
  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests TET_MESH_TET_NEIGHBORS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tet_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tet_node
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ), parameter :: tet_order = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  TET_MESH_TET_NEIGHBORS computes the 4 neighboring'
  write ( *, '(a)' ) '  tetrahedrons of each tetrahedron in a tet mesh.'
  write ( *, '(a)' ) '  containing a point.'
!
!  Set up the example tetrahedron mesh.
!
  call tet_mesh_order4_example_size ( node_num, tet_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  This mesh has tetrahedron order ', tet_order
  write ( *, '(a,i8)' ) '  The number of tetrahedrons is   ', tet_num

  allocate ( node_xyz(1:3,1:node_num) )
  allocate ( tet_node(1:tet_order,1:tet_num) )

  call tet_mesh_order4_example_set ( node_num, tet_num, node_xyz, tet_node )
!
!  Print the tets.
!
  call i4mat_transpose_print_some ( tet_order, tet_num, tet_node, &
  1, 1, tet_order, 10, '  First 10 Tets:' )
!
!  The TET_NEIGHBOR array is needed by TET_MESH_DELAUNAY_SEARCH.
!
  allocate ( tet_neighbor(4,tet_num) )

  call tet_mesh_neighbor_tets ( tet_order, tet_num, tet_node, &
    tet_neighbor )

  call i4mat_transpose_print_some ( 4, tet_num, tet_neighbor, &
    1, 1, 4, 10, '  First 10 Tet Neighbors:' )

  deallocate ( node_xyz )
  deallocate ( tet_neighbor )
  deallocate ( tet_node )

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests TET_MESH_SEARCH_NAIVE and TET_MESH_SEARCH_DELAUNAY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha(4)
  integer ( kind = 4 ) face
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz
  real ( kind = 8 ) p(3)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  integer ( kind = 4 ) tet1
  integer ( kind = 4 ) tet2
  integer ( kind = 4 ) tet3
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tet_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tet_node
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ), parameter :: tet_order = 4

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  TET_MESH_SEARCH_NAIVE uses a naive algorithm'
  write ( *, '(a)' ) '  to search a tetrahedral mesh for the tetrahedron'
  write ( *, '(a)' ) '  containing a point.'
  write ( *, '(a)' ) '  TET_MESH_SEARCH_DELAUNAY uses the Delaunay search algorithm'
  write ( *, '(a)' ) '  to search a Delaunay tetrahedral mesh for the tetrahedron'
  write ( *, '(a)' ) '  containing a point.'
!
!  Set up the example tetrahedron mesh.
!
  call tet_mesh_order4_example_size ( node_num, tet_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  This mesh has tetrahedron order ', tet_order
  write ( *, '(a,i8)' ) '  The number of tetrahedrons is   ', tet_num

  allocate ( node_xyz(1:3,1:node_num) )
  allocate ( tet_node(1:tet_order,1:tet_num) )

  call tet_mesh_order4_example_set ( node_num, tet_num, node_xyz, tet_node )
!
!  Initializing TET3 to -1 helps TET_MESH_DELAUNAY_SEARCH get started.
!
  tet3 = - 1
!
!  The TET_NEIGHBOR array is needed by TET_MESH_DELAUNAY_SEARCH.
!
  allocate ( tet_neighbor(4,tet_num) )

  call tet_mesh_neighbor_tets ( tet_order, tet_num, tet_node, &
    tet_neighbor )

  do test = 1, test_num
!
!  Choose a tetrahedron at random.
!
    tet1 = i4_uniform ( 1, tet_num, seed )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Point was chosen from tetrahedron    ', tet1
!
!  Choose a point in the tetrahedron at random.
!
    call tetrahedron_sample ( node_xyz(1:3,tet_node(1:4,tet1)), 1, seed, p )
!
!  Naive search.
!
    call tet_mesh_search_naive ( node_num, node_xyz, tet_order, tet_num, &
      tet_node, p, tet2, step_num )

    write ( *, '(a,i8,a,i8)' ) '  Naive search ended in tetrahedron    ', tet2, &
      ', number of steps = ', step_num
!
!  Delaunay search.
! 
    call tet_mesh_search_delaunay ( node_num, node_xyz, tet_order, &
      tet_num, tet_node, tet_neighbor, p, tet3, face, step_num )

    write ( *, '(a,i8,a,i8)' ) '  Delaunay search ended in tetrahedron ', tet3, &
      ', number of steps = ', step_num

  end do

  deallocate ( node_xyz )
  deallocate ( tet_neighbor )
  deallocate ( tet_node )

  return
end
