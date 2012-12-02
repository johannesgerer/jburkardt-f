program main

!*****************************************************************************80
!
!! MAIN is the main program for GRAFPACK.
!
!  Discussion:
!
!    GRAFPACK_PRB calls the GRAFPACK test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAFPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GRAFPACK library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test004 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
  call test008 ( )
  call test009 ( )
  call test0095 ( )
  call test010 ( )
  call test0105 ( )

  call test011 ( )
  call test012 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
  call test0155 ( )
  call test016 ( )
  call test017 ( )
  call test018 ( )
  call test019 ( )
  call test020 ( )

  call test021 ( )
  call test022 ( )
  call test023 ( )
  call test024 ( )
  call test025 ( )
  call test026 ( )
  call test027 ( )
  call test028 ( )
  call test029 ( )
  call test030 ( )

  call test031 ( )
  call test032 ( )
  call test033 ( )
  call test034 ( )
  call test035 ( )
  call test0335 ( )
  call test036 ( )
  call test0365 ( )
  call test0366 ( )
  call test037 ( )
  call test0375 ( )
  call test038 ( )
  call test039 ( )
  call test040 ( )

  call test041 ( )
  call test042 ( )
  call test043 ( )
  call test044 ( )
  call test045 ( )
  call test046 ( )
  call test047 ( )
  call test048 ( )
  call test049 ( )
  call test050 ( )

  call test051 ( )
  call test052 ( )
  call test053 ( )
  call test054 ( )
  call test055 ( )
  call test056 ( )
  call test057 ( )
  call test058 ( )
  call test059 ( )
  call test060 ( )

  call test061 ( )
  call test062 ( )
  call test063 ( )
  call test064 ( )
  call test065 ( )
  call test066 ( )
  call test0665 ( )
  call test067 ( )
  call test068 ( )
  call test069 ( )
  call test0695 ( )
  call test0696 ( )
  call test0697 ( )
  call test070 ( )

  call test071 ( )
  call test072 ( )
  call test073 ( )
  call test074 ( )
  call test075 ( )
  call test076 ( )
  call test077 ( )
  call test078 ( )
  call test079 ( )
  call test080 ( )

  call test081 ( )
  call test082 ( )
  call test083 ( )
  call test084 ( )
  call test085 ( )
  call test086 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAFPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests COLOR_DIGRAPH_ADJ_RANDOM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) mcolor
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) seed

  seed = 123456789
  ncolor = 3
  nedge = 15

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  COLOR_DIGRAPH_ADJ_RANDOM returns a random '
  write ( *, '(a)' ) '  color digraph.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random object is to have:'
  write ( *, '(a,i8)' ) '  Number of colors = ', ncolor
  write ( *, '(a,i8)' ) '  Number of nodes = ', nnode
  write ( *, '(a,i8)' ) '  Number of edges = ', nedge

  call color_digraph_adj_random ( nnode, ncolor, nedge, seed, adj )

  call color_digraph_adj_print ( adj, nnode, nnode, '  The color digraph:' )
!
!  Count the edges.
!
  call color_digraph_adj_edge_count ( adj, nnode, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge
!
!  Count the colors.
!
  call color_graph_adj_color_count ( adj, nnode, nnode, mcolor, ncolor )

  write ( *, '(a,i8)' ) '  Number of colors is    ', ncolor
  write ( *, '(a,i8)' ) '  Maximum color index is ', mcolor

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests COLOR_GRAPH_ADJ_CONNECT_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) mcolor
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seed

  ncolor = 3
  nedge = 8
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  COLOR_GRAPH_ADJ_CONNECT_RANDOM returns a random ' // &
    'connected color graph;'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random object is to have:'
  write ( *, '(a,i8)' ) '  Number of colors = ', ncolor
  write ( *, '(a,i8)' ) '  Number of nodes = ', nnode
  write ( *, '(a,i8)' ) '  Number of edges = ', nedge

  call color_graph_adj_connect_random ( lda, nnode, nedge, ncolor, seed, adj )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
!
!  Check connectedness.
!
  call graph_adj_is_edge_connected ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT edgewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS edgewise connected.'
  end if

  call graph_adj_is_node_connected ( adj, lda, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT nodewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS nodewise connected.'
  end if
!
!  Count the edges.
!
  call color_graph_adj_edge_count ( adj, lda, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge
!
!  Count the colors.
!
  call color_graph_adj_color_count ( adj, lda, nnode, mcolor, ncolor )
 
  write ( *, '(a,i8)' ) '  Number of colors is    ', ncolor
  write ( *, '(a,i8)' ) '  Maximum color index is ', mcolor

  return
end 
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests COLOR_GRAPH_ADJ_RANDOM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) mcolor
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) seed

  ncolor = 3
  nedge = 7
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  COLOR_GRAPH_ADJ_RANDOM returns a random color digraph.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random object is to have:'
  write ( *, '(a,i8)' ) '  Number of colors = ', ncolor
  write ( *, '(a,i8)' ) '  Number of nodes = ', nnode
  write ( *, '(a,i8)' ) '  Number of edges = ', nedge

  call color_graph_adj_random ( lda, nnode, ncolor, nedge, seed, adj )

  call color_graph_adj_print ( adj, lda, nnode, '  The color graph:' )
!
!  Count the edges.
!
  call color_graph_adj_edge_count ( adj, lda, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge
!
!  Count the colors.
!
  call color_graph_adj_color_count ( adj, lda, nnode, mcolor, ncolor )
 
  write ( *, '(a,i8)' ) '  Number of colors is    ', ncolor
  write ( *, '(a,i8)' ) '  Maximum color index is ', mcolor

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests DEGREE_SEQ_IS_GRAPHIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: ntest = 5

  integer ( kind = 4 ) degree_seq(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  DEGREE_SEQ_IS_GRAPHIC reports whether'
  write ( *, '(a)' ) '  a given sequence can represent the degree'
  write ( *, '(a)' ) '  sequence of a graph.'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    call i4vec_uniform ( nnode, 1, nnode-1, seed, degree_seq )

    call i4vec_sort_heap_d ( nnode, degree_seq )

    call i4vec_print ( nnode, degree_seq, '  The degree sequence:' )

    call degree_seq_is_graphic ( nnode, degree_seq, result )

    write ( *, '(a)' ) ' '
    if ( result == 0 ) then
      write ( *, '(a)' ) '  The sequence is NOT graphic.'
    else if ( result == 1 ) then
      write ( *, '(a)' ) '  The sequence IS graphic.'
    end if

  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests DEGREE_SEQ_TO_GRAPH_ADJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), dimension ( nnode ) :: seq = (/ 5, 5, 4, 3, 3, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  DEGREE_SEQ_TO_GRAPH_ADJ is given a degree'
  write ( *, '(a)' ) '    sequence, and constructs the adjancency'
  write ( *, '(a)' ) '    matrix of a corresponding graph.'

  call i4vec_print ( nnode, seq, '  The degree sequence:' )

  call degree_seq_to_graph_adj ( nnode, seq, lda, adj, ierror )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests DIGRAPH_ADJ_CLOSURE and DIGRAPH_ADJ_REDUCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 13
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,6) = 1
  adj(1,7) = 1

  adj(3,1) = 1

  adj(4,6) = 1

  adj(5,4) = 1

  adj(6,5) = 1

  adj(7,3) = 1
  adj(7,5) = 1
  adj(7,10) = 1

  adj(8,7) = 1
  adj(8,9) = 1

  adj(9,8) = 1

  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1

  adj(12,7) = 1
  adj(12,13) = 1

  adj(13,12) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_CLOSURE finds the transitive '
  write ( *, '(a)' ) '    closure of a digraph;'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_REDUCE finds the transitive '
  write ( *, '(a)' ) '    reduction of a digraph.'
  write ( *, '(a)' ) ' '

  call digraph_adj_print ( adj, lda, nnode, '  Adjacency matrix for G:' )

  call digraph_adj_closure ( adj, lda, nnode )

  call digraph_adj_print ( adj, lda, nnode, &
    '  Adjacency matrix for H, the transitive closure of G:' )

  call digraph_adj_reduce ( adj, nnode )

  call digraph_adj_print ( adj, lda, nnode, &
    '  Adjacency matrix for G2, the transitive reduction of H:' )

  call digraph_adj_closure ( adj, lda, nnode )

  call digraph_adj_print ( adj, lda, nnode, &
    '  Adjacency matrix for H2, the transitive closure of G2:' )

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests DIGRAPH_ADJ_COMPONENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 13
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) comp(nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_COMPONENTS finds strongly connected'
  write ( *, '(a)' ) '    components of a directed graph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,11) = 1

  adj(2,3) = 1
  adj(2,6) = 1

  adj(3,4) = 1
  adj(3,5) = 1

  adj(4,3) = 1

  adj(5,4) = 1

  adj(6,7) = 1
  adj(6,8) = 1

  adj(7,6) = 1

  adj(8,9) = 1
  adj(8,10) = 1

  adj(9,7) = 1

  adj(10,9) = 1

  adj(11,12) = 1
  adj(11,13) = 1

  adj(12,1) = 1

  adj(13,1) = 1
  adj(13,12) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph' )
 
  call digraph_adj_components ( adj, lda, nnode, ncomp, comp, dad, order )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of components = ', ncomp
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  Node, Dad, Component, Order'
  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(5i8)' ) i, dad(i), comp(i), order(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The correct components are:'
  write ( *, '(a)' ) ' (1,11,12,13), (2), (3,4,5), (6,7,8,9,10).'
!
!  Compute a reordering of the nodes.
!
  do i = 1, nnode
    order(i) = i
  end do

  do i = 2, nnode
    do j = 1, i - 1
      if ( comp(j) > comp(i) .or. &
         ( comp(j) == comp(i) .and. order(j) > order(i) ) ) then
        call i4_swap ( comp(j), comp(i) )
        call i4_swap ( order(j), order(i) )
      end if
    end do
  end do

  call i4vec2_print ( nnode, comp, order, '  I, Component(I), Node(I)' )

  call perm_inv ( nnode, order )

  call i4mat_perm ( lda, nnode, adj, order )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  return
end 
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests DIGRAPH_ADJ_CYCLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 9

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) adj2(lda,lda)
  integer ( kind = 4 ) dad(lda)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) order(lda)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_CYCLE searches for cycles in a digraph.'

  call digraph_adj_example_cycler ( adj, lda, nnode )

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )
!
!  Count the edges.
!
  call digraph_adj_edge_count ( adj, lda, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of edges is ', nedge

  call digraph_adj_cycle ( adj, lda, nnode, adj2, dad, order )

  call i4vec2_print ( nnode, dad, order, '  Node, Dad, Order' )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  Adjacency matrix with cycles marked.'
  write ( *, '(a)' ) ' '
  
  do i = 1, nnode
    write ( *, '(10i3)' ) adj2(i,1:nnode)
  end do

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests DIGRAPH_ADJ_DEGREE, DIGRAPH_ADJ_DEGREE_MAX, DIGRAPH_ADJ_DEGREE_SEQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 10

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) indegree(lda)
  integer ( kind = 4 ) indegree_max
  integer ( kind = 4 ) indegree_seq(lda)
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) outdegree(lda)
  integer ( kind = 4 ) outdegree_max
  integer ( kind = 4 ) outdegree_seq(lda)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  For a directed graph:'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_DEGREE computes the degree of the nodes;'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_DEGREE_MAX computes the maximum'
  write ( *, '(a)' ) '    degree of the nodes;'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_DEGREE_SEQ computes the degree'
  write ( *, '(a)' ) '    sequence;'

  call digraph_adj_example_cycler ( adj, lda, nnode )

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  call digraph_adj_degree ( adj, lda, nnode, indegree, outdegree )

  call i4vec2_print ( nnode, indegree, outdegree, '  Node, In/Outdegree' )

  call digraph_adj_degree_seq ( adj, lda, nnode, indegree_seq, outdegree_seq )

  call i4vec2_print ( nnode, indegree_seq, outdegree_seq, &
    '  Node, In/Outdegree sequence' )

  call digraph_adj_degree_max ( adj, lda, nnode, indegree_max, outdegree_max, &
    degree_max )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a,i8)' ) '  Maximum  indegree is ',  indegree_max
  write ( *, '(a,i8)' ) '  Maximum outdegree is ', outdegree_max
  write ( *, '(a,i8)' ) '  Maximum    degree is ',    degree_max
  write ( *, '(a)' ) ' '

  return
end
subroutine test0095 ( )

!*****************************************************************************80
!
!! TEST0095 tests DIGRAPH_ADJ_EIGEN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 9

  integer ( kind = 4 ) adj(lda,lda)
  real ( kind = 8 ) eigeni(lda)
  real ( kind = 8 ) eigenr(lda)
  integer ( kind = 4 ) neigen
  integer ( kind = 4 ) nnode

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0095'
  write ( *, '(a)' ) '  For a digraph:'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_EIGEN computes the eigenvalues.'

  call digraph_adj_example_cycler ( adj, lda, nnode )

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  call digraph_adj_eigen ( adj, lda, nnode, neigen, eigenr, eigeni )

  call r8vec2_print ( neigen, eigenr, eigeni, &
    '  Real and imaginary parts of eigenvalues:' )

  if ( neigen < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Warning!  Not all eigenvalues were computed.'
  end if

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests DIGRAPH_ADJ_HAM_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 20
  integer ( kind = 4 ), parameter :: lda = nnode
  integer ( kind = 4 ), parameter :: maxstack = 100

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical more
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_HAM_NEXT produces Hamilton circuits;'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0

  adj(1,8) = 1
  adj(1,2) = 1
  adj(1,20) = 1
  adj(2,3) = 1
  adj(2,15) = 1
  adj(3,7) = 1
  adj(3,4) = 1
  adj(4,5) = 1
  adj(4,14) = 1
  adj(5,6) = 1
  adj(5,12) = 1
  adj(6,10) = 1
  adj(6,7) = 1
  adj(7,8) = 1
  adj(8,9) = 1
  adj(9,10) = 1
  adj(9,19) = 1
  adj(10,11) = 1
  adj(11,12) = 1
  adj(11,18) = 1
  adj(12,13) = 1
  adj(13,14) = 1
  adj(13,17) = 1
  adj(14,15) = 1
  adj(15,16) = 1
  adj(16,17) = 1
  adj(16,20) = 1
  adj(17,18) = 1
  adj(18,19) = 1
  adj(19,20) = 1
 
  do i = 1, nnode-1
    do j = i+1, nnode
      if ( adj(i,j) == 1 ) then
        adj(j,i) = 1
      end if
    end do
  end do

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0

  more = .false.

  do

    call digraph_adj_ham_next ( adj, lda, nnode, circuit, stack, maxstack, &
      ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine test0105

!*****************************************************************************80
!
!! TEST0105 tests DIGRAPH_ADJ_HAM_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 9
  integer ( kind = 4 ), parameter :: lda = nnode
  integer ( kind = 4 ), parameter :: maxstack = 100

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0105'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_HAM_NEXT produces Hamilton circuits;'
  write ( *, '(a)' ) ' '

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,6) = 1
 
  adj(2,3) = 1
  adj(2,5) = 1
 
  adj(3,4) = 1
 
  adj(4,1) = 1
  adj(4,5) = 1
  adj(4,8) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
  adj(5,3) = 1
  adj(5,4) = 1
  adj(5,7) = 1
  adj(5,8) = 1
  adj(5,9) = 1
 
  adj(6,3) = 1
  adj(6,5) = 1
  adj(6,8) = 1
 
  adj(7,2) = 1
  adj(7,4) = 1
  adj(7,5) = 1
 
  adj(8,4) = 1
  adj(8,5) = 1
  adj(8,6) = 1
  adj(8,9) = 1
 
  adj(9,5) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0

  more = .false.

  do

    call digraph_adj_ham_next ( adj, lda, nnode, circuit, stack, maxstack, &
      ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine test011

!*****************************************************************************80
!
!! TEST011 tests DIGRAPH_ADJ_HAM_NEXT_BRUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 9
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iset

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_HAM_NEXT_BRUTE seeks circuits'
  write ( *, '(a)' ) '    in a directed graph which visit every node.'
  write ( *, '(a)' ) '  A brute force algorithm is used.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,6) = 1
 
  adj(2,3) = 1
  adj(2,5) = 1
 
  adj(3,4) = 1
 
  adj(4,1) = 1
  adj(4,5) = 1
  adj(4,8) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
  adj(5,3) = 1
  adj(5,4) = 1
  adj(5,7) = 1
  adj(5,8) = 1
  adj(5,9) = 1
 
  adj(6,3) = 1
  adj(6,5) = 1
  adj(6,8) = 1
 
  adj(7,2) = 1
  adj(7,4) = 1
  adj(7,5) = 1
 
  adj(8,4) = 1
  adj(8,5) = 1
  adj(8,6) = 1
  adj(8,9) = 1
 
  adj(9,5) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )
 
  iset = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0

  do
 
    call digraph_adj_ham_next_brute ( adj, lda, nnode, circuit, iset )
 
    if ( iset == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more circuits were found.'
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine test012

!*****************************************************************************80
!
!! TEST012 tests DIGRAPH_ADJ_HAM_PATH_NEXT_BRUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 4
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) iset
  integer ( kind = 4 ) j
  integer ( kind = 4 ) path(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_HAM_PATH_NEXT_BRUTE seeks paths in a'
  write ( *, '(a)' ) '    digraph which visit every node once.'
  write ( *, '(a)' ) '  A brute force algorithm is used.'
!
!  Initialize the adjacency matrix to the identity.
!
  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj(i,j) = 1
      else
        adj(i,j) = 0
      end if
    end do
  end do
!
!  Add entries for specific edges.  This is a directed graph.
!  ADJ(I, j) = 1 means there's a edge from I to J.
!
  adj(1,2) = 1
  adj(1,4) = 1
 
  adj(2,4) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
 
  adj(4,2) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )
 
  iset = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Paths:'
  write ( *, '(a)' ) ' '
  i = 0

  do
 
    call digraph_adj_ham_path_next_brute ( adj, lda, nnode, path, iset )
 
    if ( iset == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more paths were found.'
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, path(1:nnode)

  end do
 
  return
end
subroutine test013

!*****************************************************************************80
!
!! TEST013 tests DIGRAPH_ADJ_IS_EDGE_CONNECTED;
!
!  Discussion:
!
!    Here is a picture of the digraph.
!
!    1-->--2
!    |     |
!    A     A
!    |     |
!    4--<--3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 4
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_IS_EDGE_CONNECTED reports if a'
  write ( *, '(a)' ) '    digraph is edgewise connected;'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(3,2) = 1
  adj(3,4) = 1
  adj(4,1) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  call digraph_adj_is_edge_connected ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT edgewise connected.'
  else
    write ( *, '(a)' ) '  The digraph IS edgewise connected.'
  end if

  return
end 
subroutine test014

!*****************************************************************************80
!
!! TEST014 tests DIGRAPH_ADJ_IS_EULERIAN;
!
!  Discussion:
!
!    Here is a picture of the digraph:
!
!    1->---2-->---3
!        A V      V
!      6<--5--<---4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_IS_EULERIAN reports if a digraph is Eulerian;'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(3,4) = 1
  adj(4,5) = 1
  adj(5,6) = 1
  adj(6,2) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  call digraph_adj_is_eulerian ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT Eulerian.'
  else if ( result == 1 ) then
    write ( *, '(a)' ) '  The digraph IS path Eulerian.'
  else if ( result == 2 ) then
    write ( *, '(a)' ) '  The digraph IS circuit Eulerian.'
  end if

  return
end
subroutine test015

!*****************************************************************************80
!
!! TEST015 tests DIGRAPH_ADJ_IS_STRONG_CONNECTED;
!
!  Discussion:
!
!    Here are pictures of the digraphs:
!
!  1)
!
!    1-->--2
!    |     |
!    A     A
!    |     |
!    4--<--3
!
!  2)
!
!    1-->--2-->--3-->--4
!    |     |     |     |
!    A     V     A     V
!    |     |     |     |
!    5--<--6     7--<--8
!
!  3)
!
!    1-->--2-->--3-->--4
!    |     |     |     |
!    A     V     A     V
!    |     |     |     |
!    5--<--6--<--7--<--8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 8

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_IS_STRONG_CONNECTED reports if a'
  write ( *, '(a)' ) '    digraph is strongly connected;'

  nnode = 4

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(3,2) = 1
  adj(3,4) = 1
  adj(4,1) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  call digraph_adj_is_strong_connected ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT strongly connected.'
  else
    write ( *, '(a)' ) '  The digraph IS strongly connected.'
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(2,6) = 1
  adj(6,5) = 1
  adj(5,1) = 1
  adj(3,4) = 1
  adj(4,8) = 1
  adj(8,7) = 1
  adj(7,3) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  call digraph_adj_is_strong_connected ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT strongly connected.'
  else
    write ( *, '(a)' ) '  The digraph IS strongly connected.'
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(2,6) = 1
  adj(6,5) = 1
  adj(5,1) = 1
  adj(3,4) = 1
  adj(4,8) = 1
  adj(8,7) = 1
  adj(7,3) = 1
  adj(7,6) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  call digraph_adj_is_strong_connected ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT strongly connected.'
  else
    write ( *, '(a)' ) '  The digraph IS strongly connected.'
  end if

  return
end 
subroutine test0155

!*****************************************************************************80
!
!! TEST0155 tests DIGRAPH_ADJ_TOURNAMENT_RANDOM, DIGRAPH_ADJ_IS_TOURNAMENT;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seed

  seed = 123456789
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0155'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_TOURNAMENT_RANDOM returns a random'
  write ( *, '(a)' ) '    tourname digraph.'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_IS_TOURNAMENT reports if a'
  write ( *, '(a)' ) '    digraph is a tournament.'

  call digraph_adj_tournament_random ( lda, nnode, seed, adj )

  call digraph_adj_print ( adj, lda, nnode, '  A random tournament digraph:' )

  call digraph_adj_is_tournament ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT a tournament.'
  else
    write ( *, '(a)' ) '  The digraph IS a tournament.'
  end if

  return
end 
subroutine test016

!*****************************************************************************80
!
!! TEST016 tests DIGRAPH_ADJ_IS_TRANSITIVE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 12

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_IS_TRANSITIVE reports if a'
  write ( *, '(a)' ) '    digraph is transitive;'

  call digraph_adj_example_sixty ( adj, lda, nnode )

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  call digraph_adj_is_transitive ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT transitive.'
  else
    write ( *, '(a)' ) '  The digraph IS transitive.'
  end if

  return
end 
subroutine test017

!*****************************************************************************80
!
!! TEST017 tests DIGRAPH_ADJ_RANDOM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) seed

  seed = 123456789
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_RANDOM returns a random digraph.'
  write ( *, '(a)' ) ' '

  nedge = 10
  write ( *, '(a,i8)' ) '  Number of edges requested = ', nedge

  call digraph_adj_random ( lda, nnode, nedge, seed, adj )

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )
!
!  Count the edges.
!
  call digraph_adj_edge_count ( adj, lda, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge

  return
end
subroutine test018

!*****************************************************************************80
!
!! TEST018 tests DIGRAPH_ADJ_TO_DIGRAPH_ARC;
!
!    1->---2-->---3
!    A     V      V
!    6--<--5--<---4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode
  integer ( kind = 4 ), parameter :: maxarc = 10

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) inode(maxarc)
  integer ( kind = 4 ) jnode(maxarc)
  integer ( kind = 4 ) narc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_TO_DIGRAPH_ARC converts a digraph in'
  write ( *, '(a)' ) '  adjacency form to arc list form;'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(3,4) = 1
  adj(4,5) = 1
  adj(5,6) = 1
  adj(6,2) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph in adjacency form:' )

  call digraph_adj_to_digraph_arc ( adj, lda, nnode, maxarc, narc, &
    inode, jnode )

  call digraph_arc_print ( narc, inode, jnode, &
    '  The digraph in arc list form:' )

  return
end
subroutine test019

!*****************************************************************************80
!
!! TEST019 tests DIGRAPH_ADJ_TO_DIGRAPH_INC;
!
!    1->---2-->---3
!        A V      V
!      6<--5--<---4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode
  integer ( kind = 4 ), parameter :: maxarc = 10

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) inc(lda,maxarc)
  integer ( kind = 4 ) narc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_TO_DIGRAPH_INC converts a digraph in'
  write ( *, '(a)' ) '  adjacency form to incidence matrix form;'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(3,4) = 1
  adj(4,5) = 1
  adj(5,6) = 1
  adj(6,2) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph in adjacency form:' )

  call digraph_adj_to_digraph_inc ( adj, lda, nnode, maxarc, narc, inc )

  call digraph_inc_print ( lda, nnode, narc, inc, &
    '  The digraph in incidence form:' )

  return
end
subroutine test020

!*****************************************************************************80
!
!! TEST020 tests DIGRAPH_ADJ_TOP_SORT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 13

  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) node_list(nnode)
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) visit(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  DIGRAPH_ADJ_TOP_SORT does a topological sort'
  write ( *, '(a)' ) '    of an acyclic digraph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,6) = 1

  adj(5,4) = 1

  adj(6,4) = 1
  adj(6,5) = 1

  adj(7,3) = 1
  adj(7,5) = 1
  adj(7,8) = 1

  adj(8,9) = 1

  adj(10,7) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1

  adj(12,7) = 1
  adj(12,13) = 1

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )
 
  call digraph_adj_top_sort ( adj, lda, nnode, dad, visit, node_list )

  call i4vec_print ( nnode, dad, '  Nodes and "Dads":' )

  call i4vec_print ( nnode, visit, '  Nodes and order of visit:' )

  call i4vec_print ( nnode, node_list, '  Nodes and reverse topological order:' )
!
!  Invert the listing to get a permutation.
!
  order(1:nnode) = node_list(1:nnode)

  call perm_inv ( nnode, order )
!
!  Apply reordering and print adjacency matrix.
!
  call i4mat_perm ( lda, nnode, adj, order )

  call digraph_adj_print ( adj, lda, nnode, '  The reordered digraph:' )

  return
end 
subroutine test021

!*****************************************************************************80
!
!! TEST021 tests DIGRAPH_ARC_DEGREE.
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |  
!         4--7  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 11
  integer ( kind = 4 ), parameter :: nnode = 10

  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) outdegree(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  For a digraph described by an arc list:'
  write ( *, '(a)' ) '  DIGRAPH_ARC_DEGREE computes the degree of the nodes;'

  inode = (/ 1, 1,  1, 2,  2, 3, 3, 4, 4, 6,  8 /)
  jnode = (/ 3, 7, 10, 5, 10, 6, 9, 7, 8, 9, 10 /)

  call digraph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call digraph_arc_degree ( nnode, nedge, inode, jnode, indegree, outdegree )

  call i4vec2_print ( nnode, indegree, outdegree, '  Node, Indegree, Outdegree' )

  return
end 
subroutine test022

!*****************************************************************************80
!
!! TEST022 tests DIGRAPH_ARC_EULER_CIRC_NEXT, DIGRAPH_ARC_IS_EULERIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxstack = 130
  integer ( kind = 4 ), parameter :: nedge = 10
  integer ( kind = 4 ), parameter :: nnode = 5

  integer ( kind = 4 ) circuit(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 1, 3, 1, 5, 2, 4, 2, 4, 3, 5 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2, 1, 4, 1, 3, 2, 5, 3, 5, 4 /)
  logical more
  integer ( kind = 4 ) ncan(nedge)
  integer ( kind = 4 ) outdegree(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  For a digraph described by an arc list:'
  write ( *, '(a)' ) '  DIGRAPH_ARC_IS_EULERIAN checks if a graph'
  write ( *, '(a)' ) '    has an Euler circuit.'
  write ( *, '(a)' ) '  DIGRAPH_ARC_EULER_CIRC_NEXT finds the next'
  write ( *, '(a)' ) '    Euler circuit of a graph.'
  write ( *, '(a)' ) ' '

  call digraph_arc_print ( nedge, inode, jnode, '  The digraph:' )

  call digraph_arc_is_eulerian ( nnode, nedge, inode, jnode, indegree, &
    outdegree, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph is NOT eulerian.'
    return
  else if ( result == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph has an eulerian path,'
    write ( *, '(a)' ) '  but not an eulerian circuit.'
  else if ( result == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph has an eulerian circuit.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0

  more = .false.

  do

    call digraph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
      maxstack, ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nedge)

  end do

  return
end
subroutine test023

!*****************************************************************************80
!
!! TEST023 tests DIGRAPH_ARC_TO_DIGRAPH_ADJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxedge = 20
  integer ( kind = 4 ), parameter :: maxnode = 20
  integer ( kind = 4 ), parameter :: lda = maxnode

  integer ( kind = 4 ) adj(lda,maxnode)
  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  DIGRAPH_ARC_TO_DIGRAPH_ADJ converts an arclist'
  write ( *, '(a)' ) '    digraph to an adjacency digraph.'
  write ( *, '(a)' ) ' '

  call digraph_arc_example_cycler ( maxedge, nedge, inode, jnode )

  call digraph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call digraph_arc_to_digraph_adj ( nedge, inode, jnode, adj, lda, nnode )

  call digraph_adj_print ( adj, lda, nnode, '  The digraph:' )

  return
end
subroutine test024

!*****************************************************************************80
!
!! TEST024 tests FACE_CHECK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_edge = 30
  integer ( kind = 4 ), parameter :: max_order = 4
  integer ( kind = 4 ), parameter :: max_face = 10

  integer ( kind = 4 ) edge(4,max_edge)
  integer ( kind = 4 ) face(max_order,max_face)
  integer ( kind = 4 ) face_object(max_face)
  integer ( kind = 4 ) face_order(max_face)
  integer ( kind = 4 ) face_rank(max_face)
  integer ( kind = 4 ) face_tier(max_face)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) num_edge
  integer ( kind = 4 ) num_face
  integer ( kind = 4 ) num_object

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024'
  write ( *, '(a)' ) '  FACE_CHECK checks faces;'
!
!  Get the problem data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  max_face =  ', max_face
  write ( *, '(a,i8)' ) '  max_order = ', max_order
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get a test example'

  call face_example_pieces ( face, face_order, max_face, max_order, num_face )
!
!  List the problem data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Order, Nodes'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(10i3)' ) i, face_order(i), ( face(j,i), j = 1, face_order(i) )
  end do
!
!  Check the problem data.
!
  call face_check ( edge, face, face_object, face_order, face_rank, &
    face_tier, max_edge, max_order, num_edge, num_face, num_object )

  return
end
subroutine test025

!*****************************************************************************80
!
!! TEST025 tests GRAPH_ADJ_BFS.
!
!  This example is from page 22 of
!
!  Alan Gibbons,
!  Algorithmic Graph Theory,
!  Cambridge University Press, 1985
!  ISBN 0-521-28881-9
!
!  The correct result is
!
!  Node  Idad   Ideep Iorder
!
!   1      0       1    1
!   2      1       2    2
!   3      1       2    3
!   4      1       2    4
!   5      1       2    5
!   6      1       2    6
!   7      1       2    7
!   8      1       2    8
!   9      0       3    9
!  10      9       4   10
!  11     10       5   12
!  12     10       5   13
!  13      9       4   11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 13
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) deep(nnode)
  integer ( kind = 4 ) order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  GRAPH_ADJ_BFS sets up a breadth-first'
  write ( *, '(a)' ) '    traversal of a graph.'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,7) = 1
  adj(1,8) = 1
 
  adj(2,1) = 1
  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,8) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
  adj(3,7) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  adj(7,1) = 1
  adj(7,3) = 1
 
  adj(8,1) = 1
  adj(8,2) = 1
 
  adj(9,10) = 1
  adj(9,13) = 1
 
  adj(10,9) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1
 
  adj(11,10) = 1
  adj(11,12) = 1
 
  adj(12,10) = 1
  adj(12,11) = 1
 
  adj(13,9) = 1
  adj(13,10) = 1

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
 
  call graph_adj_bfs ( adj, lda, nnode, dad, deep, order )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, dad(i), deep(i), order(i)'
  write ( *, '(a)' ) ' '
 
  do i = 1, nnode
    write ( *, '(4i8)' )  i, dad(i), deep(i), order(i)
  end do
 
  return
end
subroutine test026

!*****************************************************************************80
!
!! TEST026 tests GRAPH_ADJ_BIPARTITE_RANDOM, GRAPH_ADJ_IS_BIPARTITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode1 = 4
  integer ( kind = 4 ), parameter :: nnode2 = 6
  integer ( kind = 4 ), parameter :: nnode = nnode1 + nnode2
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  GRAPH_ADJ_BIPARTITE_RANDOM returns a random ' // &
    'bipartite graph;'
  write ( *, '(a)' ) '  GRAPH_ADJ_IS_BIPARTITE reports if a graph is bipartite.'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of nodes in set 1 is ', nnode1
  write ( *, '(a,i8)' ) '  Number of nodes in set 2 is ', nnode2

  call graph_adj_bipartite_random ( lda, nnode1, nnode2, seed, nedge, adj )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  call graph_adj_is_bipartite ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT bipartite.'
  else
    write ( *, '(a)' ) '  The graph IS bipartite.'
  end if
!
!  Count the edges.
!
  call graph_adj_edge_count ( adj, lda, nnode, nedge2 )

  write ( *, '(a,i8)' ) '  Total number of edges is   ', nedge
  write ( *, '(a,i8)' ) '  Counted number of edges is ', nedge2

  return
end 
subroutine test027

!*****************************************************************************80
!
!! TEST027 tests GRAPH_ADJ_BLOCK.
!
!  The correct result is
!
!  3 blocks
!
!  Node  Idad  Iorder
!
!     1     0      -1
!     2     1       2
!     3     4       5
!     4     1      -4
!     5     4       6
!     6     2       3
!
!  Revised adjacency matrix:
!
!    0 1 0 3 3 1
!    1 0 0 0 0 1
!    0 0 0 2 0 0
!    3 0 2 0 3 0
!    3 0 0 3 0 0
!    1 1 0 0 0 0
!
!  The three blocks are defined by the edges:
!
!  1: (6,1), (2,6), (1,2)
!
!  2: (4,3)
!
!  3: (1,4), (4,5), (5,1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)

  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) nblock
  integer ( kind = 4 ) stack(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  GRAPH_ADJ_BLOCK finds the blocks in a graph.'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0
 
  adj(1,2) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
 
  adj(2,1) = 1
  adj(2,6) = 1
 
  adj(3,4) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
  adj(4,5) = 1
 
  adj(5,1) = 1
  adj(5,4) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  call graph_adj_block ( adj, lda, nnode, dad, order, stack, nblock )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of blocks = ', nblock

  call i4vec2_print ( nnode, dad, order, '  I, DAD(I), ORDER(I)' )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
 
  return
end
subroutine test028

!*****************************************************************************80
!
!! TEST028 tests GRAPH_ADJ_CLOSURE, GRAPH_ADJ_REDUCE.
!
!    1--5      2
!    | /|
!    |/ |      8--3--7
!    4  6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 8
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj(i,j) = 1
      else
        adj(i,j) = 0
      end if
    end do
  end do

  adj(1,4) = 1
  adj(1,5) = 1

  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,1) = 1
  adj(4,5) = 1

  adj(5,1) = 1
  adj(5,4) = 1
  adj(5,6) = 1

  adj(6,5) = 1

  adj(7,3) = 1

  adj(8,3) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  GRAPH_ADJ_CLOSURE finds the transitive closure '
  write ( *, '(a)' ) '    of a graph;'
  write ( *, '(a)' ) '  GRAPH_ADJ_REDUCE finds the transitive reduction'
  write ( *, '(a)' ) '    of a graph.'

  call graph_adj_print ( adj, lda, nnode, '  The adjacency matrix for G:' )

  call graph_adj_closure ( adj, lda, nnode )

  call graph_adj_print ( adj, lda, nnode, &
    '  Adjacency matrix for H, the transitive closure of G:' )

  call graph_adj_reduce ( adj, nnode )

  call graph_adj_print ( adj, lda, nnode, &
    '  Adjacency matrix for G2, the transitive reduction of H:' )

  call graph_adj_closure ( adj, lda, nnode )

  call graph_adj_print ( adj, lda, nnode, &
    '  Adjacency matrix for H2, the transitive closure of G2:' )

  return
end
subroutine test029

!*****************************************************************************80
!
!! TEST029 tests GRAPH_ADJ_COLOR_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 4
  integer ( kind = 4 ), parameter :: lda = nnode
  integer ( kind = 4 ), parameter :: maxstack = 20

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) color(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical more
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) :: ncolor = 3
  integer ( kind = 4 ) stack(maxstack)

  data ( ( adj(i,j), j = 1, nnode ), i = 1, nnode) / &
    0, 1, 0, 1, &
    1, 0, 1, 0, &
    0, 1, 0, 1, &
    1, 0, 1, 0 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029'
  write ( *, '(a)' ) '  GRAPH_ADJ_COLOR_NEXT produces colorings of a graph'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of colors available is ', ncolor

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Possible node colorings:'
  write ( *, '(a)' ) ' '

  more = .false.

  do

    call graph_adj_color_next ( adj, lda, nnode, ncolor, color, stack, &
      maxstack, ncan, more )

    if ( .not. more ) then
      exit
    end if

    write ( *, '(19i4)' ) color(1:nnode)

  end do

  return
end
subroutine test030

!*****************************************************************************80
!
!! TEST030 tests GRAPH_ADJ_CONNECT_RANDOM, GRAPH_ADJ_IS_EDGE_CONNECTED, GRAPH_ADJ_IS_NODE_CONNECTED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 8
  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seed

  seed = 123456789
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST030'
  write ( *, '(a)' ) '  GRAPH_ADJ_CONNECT_RANDOM returns a random connected graph;'
  write ( *, '(a)' ) '  GRAPH_ADJ_IS_EDGE_CONNECTED reports if a'
  write ( *, '(a)' ) '    graph is edgewise connected;'
  write ( *, '(a)' ) '  GRAPH_ADJ_IS_NODE_CONNECTED reports if a'
  write ( *, '(a)' ) '    graph is node connected;'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge

  call graph_adj_connect_random ( lda, nnode, nedge, seed, adj )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
!
!  Check connectedness.
!
  call graph_adj_is_edge_connected ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT edgewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS edgewise connected.'
  end if

  call graph_adj_is_node_connected ( adj, lda, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT nodewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS nodewise connected.'
  end if

  return
end 
subroutine test031

!*****************************************************************************80
!
!! TEST031 tests GRAPH_ADJ_CONNECT_RANDOM, GRAPH_ADJ_IS_EDGE_CONNECTED, and
!                GRAPH_ADJ_IS_NODE_CONNECTED, GRAPH_ADJ_IS_TREE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: nedge = nnode - 1
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seed

  seed = 123456789
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031'
  write ( *, '(a)' ) '  GRAPH_ADJ_CONNECT_RANDOM returns a random connected graph;'
  write ( *, '(a)' ) '  GRAPH_ADJ_IS_EDGE_CONNECTED reports if a'
  write ( *, '(a)' ) '    graph is edgewise connected;'
  write ( *, '(a)' ) '  GRAPH_ADJ_IS_NODE_CONNECTED reports if a'
  write ( *, '(a)' ) '    graph is node connected;'
  write ( *, '(a)' ) '  GRAPH_ADJ_IS_TREE reports if a graph is a tree.'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge

  call graph_adj_connect_random ( lda, nnode, nedge, seed, adj )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
!
!  Check connectedness.
!
  call graph_adj_is_edge_connected ( adj, lda, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT edgewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS edgewise connected.'
  end if

  call graph_adj_is_node_connected ( adj, lda, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT nodewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS nodewise connected.'
  end if
!
!  Check arboricity.
!
  call graph_adj_is_tree ( adj, lda, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT a tree.'
  else
    write ( *, '(a)' ) '  The graph IS a tree.'
  end if

  return
end 
subroutine test032

!*****************************************************************************80
!
!! TEST032 tests GRAPH_ADJ_CYCLE.
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |  
!         4--7  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxstack = 100
  integer ( kind = 4 ), parameter :: nnode = 10

  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032'
  write ( *, '(a)' ) '  GRAPH_ADJ_CYCLE searches for cycles in a graph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,7) = 1
  adj(1,10) = 1

  adj(2,5) = 1
  adj(2,10) = 1

  adj(3,1) = 1
  adj(3,6) = 1
  adj(3,9) = 1

  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,2) = 1

  adj(6,3) = 1
  adj(6,9) = 1

  adj(7,1) = 1
  adj(7,4) = 1
 
  adj(8,4) = 1
  adj(8,10) = 1

  adj(9,3) = 1
  adj(9,6) = 1

  adj(10,1) = 1
  adj(10,2) = 1
  adj(10,8) = 1
   
  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  call graph_adj_cycle ( adj, lda, nnode, dad, order, maxstack, stack )

  call i4vec2_print ( nnode, dad, order, '  Node, Dad, Order' )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  Adjacency matrix with cycles marked.'
  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(10i3)') adj(i,1:nnode)
  end do

  return
end 
subroutine test033

!*****************************************************************************80
!
!! TEST033 tests GRAPH_ADJ_DEGREE, GRAPH_ADJ_DEGREE_MAX, GRAPH_ADJ_DEGREE_SEQ.
!
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |  
!         4--7  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 10
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) degree_seq(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033'
  write ( *, '(a)' ) '  For a graph:'
  write ( *, '(a)' ) '  GRAPH_ADJ_DEGREE computes the degree of the nodes;'
  write ( *, '(a)' ) '  GRAPH_ADJ_DEGREE_MAX computes the maximum'
  write ( *, '(a)' ) '    degree of the nodes;'
  write ( *, '(a)' ) '  GRAPH_ADJ_DEGREE_SEQ computes the degree sequence;'

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,7) = 1
  adj(1,10) = 1

  adj(2,5) = 1
  adj(2,10) = 1

  adj(3,1) = 1
  adj(3,6) = 1
  adj(3,9) = 1

  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,2) = 1

  adj(6,3) = 1
  adj(6,9) = 1

  adj(7,1) = 1
  adj(7,4) = 1
 
  adj(8,4) = 1
  adj(8,10) = 1

  adj(9,3) = 1
  adj(9,6) = 1

  adj(10,1) = 1
  adj(10,2) = 1
  adj(10,8) = 1
   
  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  call graph_adj_degree ( adj, lda, nnode, degree )

  call i4vec_print ( nnode, degree, '  Node degrees:' )

  call graph_adj_degree_seq ( adj, lda, nnode, degree_seq )

  call i4vec_print ( nnode, degree_seq, '  Degree sequence:' )

  call graph_adj_degree_max ( adj, lda, nnode, degree_max )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a,i8)' ) '  Maximum node degree is ', degree_max
  write ( *, '(a)' ) ' '

  return
end 
subroutine test034

!*****************************************************************************80
!
!! TEST034 tests GRAPH_ADJ_DFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 13
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034'
  write ( *, '(a)' ) '  GRAPH_ADJ_DFS does depth first search of graph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,6) = 1
  adj(1,7) = 1

  adj(5,4) = 1
  adj(5,7) = 1

  adj(6,5) = 1

  adj(8,9) = 1

  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1

  adj(12,13) = 1

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  call graph_adj_dfs ( adj, lda, nnode, dad, order )

  call i4vec2_print ( nnode, dad, order, '  Node, Dad, Order' )

  return
end 
subroutine test0335

!*****************************************************************************80
!
!! TEST0335 tests GRAPH_ADJ_EIGEN.
!
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |
!         4--7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 10
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  real ( kind = 8 ) eigen(nnode)
  integer ( kind = 4 ) neigen

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0335'
  write ( *, '(a)' ) '  For a graph:'
  write ( *, '(a)' ) '  GRAPH_ADJ_EIGEN computes the eigenvalues.'

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,7) = 1
  adj(1,10) = 1

  adj(2,5) = 1
  adj(2,10) = 1

  adj(3,1) = 1
  adj(3,6) = 1
  adj(3,9) = 1

  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,2) = 1

  adj(6,3) = 1
  adj(6,9) = 1

  adj(7,1) = 1
  adj(7,4) = 1

  adj(8,4) = 1
  adj(8,10) = 1

  adj(9,3) = 1
  adj(9,6) = 1

  adj(10,1) = 1
  adj(10,2) = 1
  adj(10,8) = 1

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  call graph_adj_eigen ( adj, lda, nnode, neigen, eigen )

  call r8vec_print ( neigen, eigen, '  The eigenvalues:' )

  if ( neigen < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Warning!  Not all eigenvalues were computed.'
  end if

  return
end
subroutine test035

!*****************************************************************************80
!
!! TEST035 tests GRAPH_ADJ_DFS_2.
!
!  Discussion:
!
!    This example is from page 22 of
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985
!    ISBN 0-521-28881-9
!
!    The correct result is
!
!    Node  Idad  Iorder
!
!     1      0      1
!     2      1      2
!     3      1      6
!     4      3      7
!     5      2      3
!     6      2      4
!     7      3      8
!     8      2      5
!     9      0      9
!    10      9     10
!    11     10     11
!    12     10     12
!    13     10     13
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 13
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  GRAPH_ADJ_DFS_2 sets up depth-first traversal'
  write ( *, '(a)' ) '    of a graph described by an adjacency matrix.'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,7) = 1
  adj(1,8) = 1
 
  adj(2,1) = 1
  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,8) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
  adj(3,7) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  adj(7,1) = 1
  adj(7,3) = 1
 
  adj(8,1) = 1
  adj(8,2) = 1
 
  adj(9,10) = 1
  adj(9,13) = 1
 
  adj(10,9) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1
 
  adj(11,10) = 1
  adj(11,12) = 1
 
  adj(12,10) = 1
  adj(12,11) = 1
 
  adj(13,9) = 1
  adj(13,10) = 1

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
 
  call graph_adj_dfs_2 ( adj, lda, nnode, dad, order )
 
  call i4vec2_print ( nnode, dad, order, '  I, DAD(I), ORDER(I)' )
 
  return
end
subroutine test036

!*****************************************************************************80
!
!! TEST036 tests GRAPH_ADJ_HAM_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 20
  integer ( kind = 4 ), parameter :: lda = nnode
  integer ( kind = 4 ), parameter :: maxstack = 100

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical more
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036'
  write ( *, '(a)' ) '  GRAPH_ADJ_HAM_NEXT produces Hamilton circuits;'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,8) = 1
  adj(1,20) = 1

  adj(2,1) = 1
  adj(2,3) = 1
  adj(2,15) = 1

  adj(3,2) = 1
  adj(3,7) = 1
  adj(3,4) = 1

  adj(4,3) = 1
  adj(4,5) = 1
  adj(4,14) = 1

  adj(5,4) = 1
  adj(5,6) = 1
  adj(5,12) = 1

  adj(6,10) = 1
  adj(6,7) = 1

  adj(7,3) = 1
  adj(7,6) = 1
  adj(7,8) = 1

  adj(8,1) = 1
  adj(8,7) = 1
  adj(8,9) = 1

  adj(9,8) = 1
  adj(9,10) = 1
  adj(9,19) = 1

  adj(10,6) = 1
  adj(10,9) = 1
  adj(10,11) = 1

  adj(11,10) = 1
  adj(11,12) = 1
  adj(11,18) = 1

  adj(12,5) = 1
  adj(12,11) = 1
  adj(12,13) = 1

  adj(13,12) = 1
  adj(13,14) = 1
  adj(13,17) = 1

  adj(14,4) = 1
  adj(14,13) = 1
  adj(14,15) = 1

  adj(15,2) = 1
  adj(15,14) = 1
  adj(15,16) = 1

  adj(16,15) = 1
  adj(16,17) = 1
  adj(16,20) = 1

  adj(17,13) = 1
  adj(17,16) = 1
  adj(17,18) = 1

  adj(18,11) = 1
  adj(18,17) = 1
  adj(18,19) = 1

  adj(19,9) = 1
  adj(19,18) = 1
  adj(19,20) = 1

  adj(20,1) = 1
  adj(20,16) = 1
  adj(20,19) = 1
 
  do i = 1, nnode-1
    do j = i+1, nnode
      if ( adj(i,j) == 1 ) then
        adj(j,i) = 1
      end if
    end do
  end do

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '

  i = 0

  more = .false.

  do

    call graph_adj_ham_next ( adj, lda, nnode, circuit, stack, maxstack, &
      ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine test0365

!*****************************************************************************80
!
!! TEST0365 tests GRAPH_ADJ_HAM_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 9
  integer ( kind = 4 ), parameter :: lda = nnode
  integer ( kind = 4 ), parameter :: maxstack = 100

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0365'
  write ( *, '(a)' ) '  GRAPH_ADJ_HAM_NEXT produces Hamilton circuits;'
  write ( *, '(a)' ) ' '

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,4) = 1
  adj(1,6) = 1
 
  adj(2,1) = 1
  adj(2,3) = 1
  adj(2,7) = 1
 
  adj(3,2) = 1
  adj(3,4) = 1
  adj(3,6) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
  adj(4,7) = 1
 
  adj(5,6) = 1
  adj(5,7) = 1
  adj(5,9) = 1
 
  adj(6,1) = 1
  adj(6,3) = 1
  adj(6,5) = 1
  adj(6,8) = 1
 
  adj(7,2) = 1
  adj(7,4) = 1
  adj(7,5) = 1
 
  adj(8,6) = 1
  adj(8,9) = 1
 
  adj(9,5) = 1
  adj(9,8) = 1

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '

  i = 0

  more = .false.

  do

    call graph_adj_ham_next ( adj, lda, nnode, circuit, stack, maxstack, &
      ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(2x,i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine test0366

!*****************************************************************************80
!
!! TEST0366 tests GRAPH_ADJ_HAM_NEXT_BRUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 9
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iset

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0366'
  write ( *, '(a)' ) '  GRAPH_ADJ_HAM_NEXT_BRUTE seeks circuits'
  write ( *, '(a)' ) '    in a graph which visit every node.'
  write ( *, '(a)' ) '  A brute force algorithm is used.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,4) = 1
  adj(1,6) = 1
 
  adj(2,1) = 1
  adj(2,3) = 1
  adj(2,7) = 1
 
  adj(3,2) = 1
  adj(3,4) = 1
  adj(3,6) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
  adj(4,7) = 1
 
  adj(5,6) = 1
  adj(5,7) = 1
  adj(5,9) = 1
 
  adj(6,1) = 1
  adj(6,3) = 1
  adj(6,5) = 1
  adj(6,8) = 1
 
  adj(7,2) = 1
  adj(7,4) = 1
  adj(7,5) = 1
 
  adj(8,6) = 1
  adj(8,9) = 1
 
  adj(9,5) = 1
  adj(9,8) = 1

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
 
  iset = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '

  i = 0
  
  do
 
    call graph_adj_ham_next_brute ( adj, lda, nnode, circuit, iset )
 
    if ( iset == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more circuits were found.'
      exit
    end if

    i = i + 1
    write ( *, '(2x,i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine test037

!*****************************************************************************80
!
!! TEST037 tests GRAPH_ADJ_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 10
  integer ( kind = 4 ), parameter :: nnode = 6

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037'
  write ( *, '(a)' ) '  GRAPH_ADJ_RANDOM returns a random graph;'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of edges requested = ', nedge

  call graph_adj_random ( nnode, nedge, seed, adj )

  call graph_adj_print ( adj, nnode, nnode, '  The graph:' )

  return
end
subroutine test0375

!*****************************************************************************80
!
!! TEST0375 tests GRAPH_ADJ_RANDOM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 20
  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) adj(nnode,nnode)
  real ( kind = 8 ) eigen(nnode)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) neigen
  real ( kind = 8 ) prob
  real ( kind = 8 ), dimension ( test_num ) :: prob_test = (/ &
    0.25D+00, 0.40D+00,  0.65D+00 /)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0375'
  write ( *, '(a)' ) '  GRAPH_ADJ_RANDOM2 returns a random graph, for which'
  write ( *, '(a)' ) '  edges are generated with a given probability.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we show the effect of increasing connectivity'
  write ( *, '(a)' ) '  on the singularity of the adjacency matrix.'

  do test = 1, test_num

    prob = prob_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Probability of edge generation = ', prob

    call graph_adj_random2 ( nnode, prob, seed, nedge, adj )

    write ( *, '(a,i8)' ) '  Number of edges generated = ', nedge
    write ( *, '(a,g14.6)' ) '  Ratio = ', &
      real ( nedge, kind = 8 ) / real ( ( nnode * ( nnode - 1 ) ) / 2, kind = 8 )

    call graph_adj_print ( adj, nnode, nnode, '  The graph:' )

    call graph_adj_eigen ( adj, nnode, nnode, neigen, eigen )

    call r8vec_print ( neigen, eigen, '  The eigenvalues:' )

  end do

  return
end 
subroutine test038

!*****************************************************************************80
!
!! TEST038 tests GRAPH_ADJ_SPAN_TREE, GRAPH_ADJ_SPAN_TREE_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 13
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) tree_num

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,7) = 1
  adj(1,8) = 1
 
  adj(2,1) = 1
  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,8) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
  adj(3,7) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  adj(7,1) = 1
  adj(7,3) = 1
 
  adj(8,1) = 1
  adj(8,2) = 1
  adj(8,9) = 1

  adj(9,8) = 1 
  adj(9,10) = 1
  adj(9,13) = 1
 
  adj(10,9) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1
 
  adj(11,10) = 1
  adj(11,12) = 1
 
  adj(12,10) = 1
  adj(12,11) = 1
 
  adj(13,9) = 1
  adj(13,10) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038'
  write ( *, '(a)' ) '  GRAPH_ADJ_SPAN_TREE constructs a spanning tree of a graph.'
  write ( *, '(a)' ) '  GRAPH_ADJ_SPAN_TREE_ENUM enumerates the spanning trees'
  write ( *, '(a)' ) '     of a graph.'

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  call graph_adj_span_tree_enum ( adj, lda, nnode, tree_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of spanning trees is ', tree_num

  call graph_adj_span_tree ( adj, lda, nnode, inode, jnode )

  call graph_arc_print ( nnode-1, inode, jnode, '  The spanning tree:' )

  return
end
subroutine test039

!*****************************************************************************80
!
!! TEST039 tests GRAPH_ARC_EDGE_CON2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 17
  integer ( kind = 4 ), parameter :: nnode = 9

  integer ( kind = 4 ) edge_con
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 6,2,3,6,7,1,4,7,3,4,9,6,5,4,2,9,4 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 8,5,1,3,2,8,3,5,8,1,2,1,9,8,6,7,2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST039'
  write ( *, '(a)' ) '  GRAPH_ARC_EDGE_CON2 finds graph edge connectivity.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  call graph_arc_edge_con2 ( nnode, nedge, inode, jnode, edge_con )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The computed edge connectivity is ', edge_con

  return
end
subroutine test040

!*****************************************************************************80
!
!  TEST040 tests GRAPH_ARC_MATCH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: nedge = 14
  integer ( kind = 4 ), parameter :: nnode = 12

  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 6, 9, 3,  4, 11, 6, 4, 5,  6, 10, 3, 4, 1, 3 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 2, 7, 7, 10,  5, 8, 6, 7, 12,  2, 1, 2, 5, 5 /)
  integer ( kind = 4 ), dimension ( nnode ) :: match
  integer ( kind = 4 ), dimension ( nnode ) :: type = (/ &
    1, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040'
  write ( *, '(a)' ) '  GRAPH_ARC_MATCH finds a maximal matching in a graph.'

  call graph_arc_print ( nedge, inode, jnode, '  The edge list of the graph:' )

  call i4vec_print ( nnode, type, '  Nodes and their types:' )

  call graph_arc_match ( nnode, nedge, inode, jnode, type, match )

  call i4vec_print ( nnode, match, '  Node and matching node:' )

  return
end
subroutine test041

!*****************************************************************************80
!
!! TEST041 tests GRAPH_ARC_MIN_PATH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: lda = nnode
  integer ( kind = 4 ), parameter :: nedge = 6

  real ( kind = 8 ), save, dimension ( nedge ) :: cost = (/ &
    1.0D+00, 1.0D+00, 3.0D+00, 2.0D+00, 2.0D+00, 5.0D+00 /)
  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ), save, dimension ( nedge ) :: inode = (/ 1, 1, 2, 2, 3, 3 /)
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) istop
  integer ( kind = 4 ), save, dimension ( nedge ) :: jnode = (/ 2, 3, 3, 5, 4, 5 /)
  integer ( kind = 4 ) num_path
  integer ( kind = 4 ) path(nnode)
  real ( kind = 8 ) path_length

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041'
  write ( *, '(a)' ) '  GRAPH_ARC_MIN_PATH computes the shortest path from one'
  write ( *, '(a)' ) '  node to another.'
  write ( *, '(a)' ) ' '

  call graph_arc_weight_print ( nedge, inode, jnode, cost, &
    '  The weighted graph:' )

  dist(1:nnode,1:nnode) = 0.0D+00

  do istart = 1, nnode
    do istop = istart+1, nnode
      call graph_arc_min_path ( nnode, nedge, inode, jnode, cost, istart, &
        istop, num_path, path, path_length )
      dist(istart,istop) = path_length
      dist(istop,istart) = path_length
    end do
  end do

  call graph_dist_print ( dist, lda, nnode, &
    '  The distance matrix constructed by GRAPH_ARC_MIN_PATH:' )
 
  istart = 4
  istop = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The routine actually also computes the path.'
  write ( *, '(a,i8)' ) '  For instance, starting at node ', istart
  write ( *, '(a,i8)' ) '  we compute the shortest path to node ', istop

  call graph_arc_min_path ( nnode, nedge, inode, jnode, cost, istart, &
    istop, num_path, path, path_length )

  call i4vec_print ( num_path, path, '  The path:' )

  return
end
subroutine test042

!*****************************************************************************80
!
!! TEST042 tests GRAPH_ARC_MIN_SPAN_TREE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 10
  integer ( kind = 4 ), parameter :: nnode = 5

  real ( kind = 8 ), dimension ( nedge ) :: cost = &
    (/ 100.0, 125.0, 120.0, 110.0, 40.0, 65.0, 60.0, 45.0, 55.0, 50.0 /)
  real ( kind = 8 ), dimension ( nnode-1) :: ctree
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)
  integer ( kind = 4 ) jtree(nnode-1)
  real ( kind = 8 ) tree_cost

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST042'
  write ( *, '(a)' ) '  GRAPH_ARC_MIN_SPAN_TREE finds a minimum length'
  write ( *, '(a)' ) '  spanning tree.'
  write ( *, '(a)' ) ' '

  call graph_arc_weight_print ( nedge, inode, jnode, cost, &
    '  The weighted graph:' )

  call graph_arc_min_span_tree ( nnode, nedge, inode, jnode, cost, &
    itree, jtree, tree_cost )

  do i = 1, nnode-1
    ctree(i) = 0.0D+00
    do j = 1, nedge
      if ( ( inode(j) == itree(i) .and. jnode(j) == jtree(i) ) .or. &
           ( inode(j) == jtree(i) .and. jnode(j) == itree(i) ) ) then
        ctree(i) = cost(j)
        exit
      end if
    end do
  end do

  call graph_arc_weight_print ( nnode-1, itree, jtree, ctree, &
    '  The minimal spanning tree:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( ctree )

  return
end
subroutine test043

!*****************************************************************************80
!
!! TEST043 tests GRAPH_ARC_SPAN_FOREST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 14
  integer ( kind = 4 ), parameter :: nedge = 10

  integer ( kind = 4 ) component(nnode)
  integer ( kind = 4 ), save, dimension ( nedge ) :: inode = &
    (/  2,  4,  1,  7,  5,  2,  6,  2,  3,  4 /)
  integer ( kind = 4 ), save, dimension ( nedge ) :: jnode = &
    (/  3,  7,  9, 11,  8,  5, 10,  8,  8, 11 /)
  integer ( kind = 4 ) ncomp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST043'
  write ( *, '(a)' ) '  GRAPH_ARC_SPAN_FOREST'
  write ( *, '(a)' ) '  computes a spanning forest for a graph'
 
  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_span_forest ( nnode, nedge, inode, jnode, ncomp, component )

  call graph_arc_print ( nedge, inode, jnode, &
    '  The reordered endpoint array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of connected components = ', ncomp

  call i4vec_print ( nnode, component, '  Node component membership:' )

  return
end
subroutine test044

!*****************************************************************************80
!
!! TEST044 tests GRAPH_ARC_TO_DIGRAPH_ARC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 8
  integer ( kind = 4 ), parameter :: maxarc = 2 * nedge

  integer ( kind = 4 ) iarc(maxarc)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 1, 1, 1, 2, 3, 4, 2, 4 /)
  integer ( kind = 4 ) jarc(maxarc)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2, 1, 4, 1, 2, 1, 3, 2 /)
  integer ( kind = 4 ) narc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST044'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_TO_DIGRAPH_ARC makes a directed graph'
  write ( *, '(a)' ) '    from an undirected one.'

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_to_digraph_arc ( iarc, jarc, inode, jnode, maxarc, narc, &
    nedge )

  call digraph_arc_print ( narc, iarc, jarc, '  The digraph:' )

  return
end
subroutine test045

!*****************************************************************************80
!
!! TEST045 tests GRAPH_ARC_TO_GRAPH_ADJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 8
  integer ( kind = 4 ), parameter :: maxnode = 5
  integer ( kind = 4 ), parameter :: lda = maxnode

  integer ( kind = 4 ) adj(lda,maxnode)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 1, 1, 1, 2, 3, 4, 2, 4 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2, 1, 4, 1, 2, 1, 3, 2 /)
  integer ( kind = 4 ) nnode

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_TO_GRAPH_ADJ converts an arclist'
  write ( *, '(a)' ) '    graph to an adjacency graph.'
  write ( *, '(a)' ) ' '

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_to_graph_adj ( nedge, inode, jnode, adj, lda, nnode )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )

  return
end
subroutine test046

!*****************************************************************************80
!
!! TEST046 tests GRAPH_ARC_COMPLEMENT, GRAPH_ARC_EDGE_SORT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxedge = 90
  integer ( kind = 4 ), parameter :: maxnode = 10

  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) inode2(maxedge)
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) jnode2(maxedge)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) nnode
  real ( kind = 8 ) x(maxnode)
  real ( kind = 8 ) y(maxnode)
  real ( kind = 8 ) z(maxnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST046'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_COMPLEMENT computes the complement'
  write ( *, '(a)' ) '    of a graph described by its edge array;'
  write ( *, '(a)' ) '  GRAPH_ARC_EDGE_SORT sorts the edge array.'

  call graph_arc_example_diamond ( inode, jnode, maxedge, nedge, nnode, x, y, z )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges in original graph is ', nedge
  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
 
  call graph_arc_edge_sort ( nedge, inode, jnode )
 
  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )
 
  call graph_arc_complement ( inode, jnode, inode2, jnode2, maxedge, nedge, &
    nedge2, nnode )
 
  write ( *, '(a,i8)' ) 'Number of edges in complement is ', nedge2

  call graph_arc_edge_sort ( nedge2, inode2, jnode2 )
 
  call graph_arc_print ( nedge, inode, jnode, '  The complement graph:' )
 
  return
end
subroutine test047

!*****************************************************************************80
!
!! TEST047 tests GRAPH_ARC_DEGREE.
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |  
!         4--7  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 11
  integer ( kind = 4 ), parameter :: nnode = 10

  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST047'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_DEGREE computes the degree of the nodes;'

  inode = (/ 1, 1,  1, 2,  2, 3, 3, 4, 4, 6,  8 /)
  jnode = (/ 3, 7, 10, 5, 10, 6, 9, 7, 8, 9, 10 /)

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )

  call i4vec_print ( nnode, degree, '  The node degrees:' )

  return
end 
subroutine test048

!*****************************************************************************80
!
!! TEST048 tests GRAPH_ARC_DEGREE.
!
!
!  5--2--100-1--3--0
!         |  |  | /
!        88  |  9
!         |  |  
!      (-4)--7  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 11

  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 1, 1,   1, 2,   2, 3, 3, -4, -4, 0,  88 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 3, 7, 100, 5, 100, 0, 9,  7, 88, 9, 100 /)
  integer ( kind = 4 ) mnode
  integer ( kind = 4 ) nnode

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST048'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_NODE_COUNT counts the nodes and'
  write ( *, '(a)' ) '  finds the highest label.'

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_node_count ( nedge, inode, jnode, mnode, nnode )

  write ( *, '(a,i8)' ) '  Number of nodes is    ', nnode
  write ( *, '(a,i8)' ) '  Maximum node label is ', mnode

  return
end 
subroutine test049

!*****************************************************************************80
!
!! TEST049 tests GRAPH_ARC_EULER_CIRC_NEXT, GRAPH_ARC_IS_EULERIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxstack = 130
  integer ( kind = 4 ), parameter :: nedge = 10
  integer ( kind = 4 ), parameter :: nnode = 5

  integer ( kind = 4 ) circuit(nedge)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)
  logical more
  integer ( kind = 4 ) ncan(nedge)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST049'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_IS_EULERIAN checks if a graph has an'
  write ( *, '(a)' ) '    Euler circuit.'
  write ( *, '(a)' ) '  GRAPH_ARC_EULER_CIRC_NEXT finds the next'
  write ( *, '(a)' ) '    Euler circuit of a graph.'
  write ( *, '(a)' ) ' '

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_is_eulerian ( nnode, nedge, inode, jnode, degree, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The graph is NOT eulerian.'
    return
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The graph is eulerian.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0
  more = .false.

  do

    call graph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
      maxstack, ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nedge)

  end do

  return
end
subroutine test050

!*****************************************************************************80
!
!! TEST050 tests GRAPH_ARC_EULER_CIRC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 10
  integer ( kind = 4 ), parameter :: nnode = 5

  integer ( kind = 4 ) circuit(nedge)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)
  integer ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST050'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_IS_EULERIAN determines if a graph'
  write ( *, '(a)' ) '    is Eulerian;'
  write ( *, '(a)' ) '  GRAPH_ARC_EULER_CIRC returns an Euler circuit'
  write ( *, '(a)' ) '    of a graph.'

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_is_eulerian ( nnode, nedge, inode, jnode, degree, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The graph is NOT eulerian.'
    return
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The graph is eulerian.'
  end if

  call graph_arc_euler_circ ( nnode, nedge, inode, jnode, circuit )

  call i4vec_print ( nedge, circuit, '  The nodes in the Euler circuit:' )

  return
end
subroutine test051

!*****************************************************************************80
!
!! TEST051 tests GRAPH_ARC_SPAN_TREE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 18
  integer ( kind = 4 ), parameter :: nnode = 13

  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)

  inode(1) = 1
  jnode(1) = 2
  inode(2) = 1
  jnode(2) = 3
  inode(3) = 1
  jnode(3) = 4
  inode(4) = 1
  jnode(4) = 5
  inode(5) = 1
  jnode(5) = 6
  inode(6) = 1
  jnode(6) = 7
  inode(7) = 1
  jnode(7) = 8
 
  inode(8) = 2
  jnode(8) = 5
  inode(9) = 2
  jnode(9) = 6
  inode(10) = 2
  jnode(10) = 8
 
  inode(11) = 3
  jnode(11) = 4
  inode(12) = 3
  jnode(12) = 7
 
  inode(13) = 9
  jnode(13) = 10
  inode(14) = 9
  jnode(14) = 13
 
  inode(15) = 10
  jnode(15) = 11
  inode(16) = 10
  jnode(16) = 12
  inode(17) = 10
  jnode(17) = 13
 
  inode(18) = 11
  jnode(18) = 12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST051'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_SPAN_TREE constructs a spanning tree.'
  write ( *, '(a)' ) ' '

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_span_tree ( nedge, inode, jnode, nnode, dad )

  call i4vec_print ( nnode, dad, '  Nodes and Parent Nodes:' )

  return
end
subroutine test052

!*****************************************************************************80
!
!! TEST052 tests GRAPH_CHRO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6 
  integer ( kind = 4 ), parameter :: nedge = 12
  integer ( kind = 4 ), parameter :: maxstack = nnode * nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(nnode)
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jarray(nnode)
  integer ( kind = 4 ) karray(nnode)
  integer ( kind = 4 ) stack(2,maxstack)

  data ( ( iendpt(i,j), i = 1, 2 ), j = 1, nedge ) / &
    1,2, 1,3, 1,4, 1,5, 2,3, 2,4, 2,6, 3,5, 3,6, 4,5, 4,6, 5,6 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052'
  write ( *, '(a)' ) '  GRAPH_CHRO finds the chromatic polynomial of a graph.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The end point array:'
  write ( *, '(a)' ) ' '
  write ( *, '(19i4)' ) ( iendpt(1,i), i = 1, nedge )
  write ( *, '(19i4)' ) ( iendpt(2,i), i = 1, nedge )
 
  call graph_chro ( nnode, nedge, iendpt, iarray, jarray, karray, &
    stack, maxstack )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The chromatic polynomial:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Power sum form:'
  write ( *, '(19i4)' ) iarray(1:nnode)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tutte or tree form:'
  write ( *, '(19i4)' ) jarray(1:nnode)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Stirling form:'
  write ( *, '(19i4)' ) karray(1:nnode)
 
  return
end
subroutine test053

!*****************************************************************************80
!
!! TEST053 tests GRAPH_DIST_ALL.
!
!  The graph is:
!
!  N3 --3-- N2 --4-- N4 --5-- N5
!
!     \      |      /
!       6    2     1
!        \   |    /
!
!            N1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dinfin
  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) i
  real ( kind = 8 ) path_dist(lda,nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST053'
  write ( *, '(a)' ) '  GRAPH_DIST_ALL computes the distance between'
  write ( *, '(a)' ) '    all pairs of nodes.'
  write ( *, '(a)' ) ' '

  dinfin = 1000.0D+00
 
  dist(1:nnode,1:nnode) = dinfin

  do i = 1, nnode
    dist(i,i) = 0.0D+00
  end do
 
  dist(1,2) = 2.0D+00
  dist(1,3) = 6.0D+00
  dist(1,4) = 1.0D+00

  dist(2,1) = 2.0D+00
  dist(2,3) = 3.0D+00
  dist(2,4) = 4.0D+00

  dist(3,1) = 6.0D+00
  dist(3,2) = 3.0D+00

  dist(4,1) = 1.0D+00
  dist(4,2) = 4.0D+00
  dist(4,5) = 5.0D+00

  dist(5,4) = 5.0D+00
 
  call graph_dist_print ( dist, lda, nnode, &
    '  Immediate node distance matrix:' )

  call graph_dist_all ( dist, dinfin, lda, nnode, path_dist )
 
  call graph_dist_print ( path_dist, lda, nnode, &
    '  Total node distance matrix:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Note that "infinity" is represented by ', dinfin
 
  return
end
subroutine test054

!*****************************************************************************80
!
!! TEST054 tests GRAPH_DIST_CHECK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 15
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) a(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j

  data ( ( a(i,j), j = 1, nnode ), i = 1, nnode ) / &
     0., 29., 82., 46., 68., 52., 72., 42., 51., 55., 29., 74., 23., 72., 46., &
    29.,  0., 55., 46., 42., 43., 43., 23., 23., 31., 41., 51., 11., 52., 21., &
    82., 55.,  0., 68., 46., 55., 23., 43., 41., 29., 79., 21., 64., 31., 51., &
    46., 46., 68.,  0., 82., 15., 72., 31., 62., 42., 21., 51., 51., 43., 64., &
    68., 42., 46., 82.,  0., 74., 23., 52., 21., 46., 82., 58., 46., 65., 23., &
    52., 43., 55., 15., 74.,  0., 61., 23., 55., 31., 33., 37., 51., 29., 59., &
    72., 43., 23., 72., 23., 61.,  0., 42., 23., 31., 77., 37., 51., 46., 33., &
    42., 23., 43., 31., 52., 23., 42.,  0., 33., 15., 37., 33., 33., 31., 37., &
    51., 23., 41., 62., 21., 55., 23., 33.,  0., 29., 62., 46., 29., 51., 11., &
    55., 31., 29., 42., 46., 31., 31., 15., 29.,  0., 51., 21., 41., 23., 37., &
    29., 41., 79., 21., 82., 33., 77., 37., 62., 51.,  0., 65., 42., 59., 61., &
    74., 51., 21., 51., 58., 37., 37., 33., 46., 21., 65.,  0., 61., 11., 55., &
    23., 11., 64., 51., 46., 51., 51., 33., 29., 41., 42., 61.,  0., 62., 23., &
    72., 52., 31., 43., 65., 29., 46., 31., 51., 23., 59., 11., 62.,  0., 59., &
    46., 21., 51., 64., 23., 59., 33., 37., 11., 37., 61., 55., 23., 59.,  0. /
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054'
  write ( *, '(a)' ) '  GRAPH_DIST_CHECK checks a distance matrix.'

  call graph_dist_check ( a, lda, nnode, ierror )

  if ( ierror == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The distance matrix passed all tests.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'The distance matrix failed test ', ierror
  end if

  return
end
subroutine test055

!*****************************************************************************80
!
!! TEST055 tests GRAPH_DIST_MIN_SPAN_TREE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtree(nnode-1)
  real ( kind = 8 ) wtree(nnode-1)

  data ( ( dist(i,j), i = 1, nnode ), j = 1, nnode ) / &
       0.0D+00, 100.0D+00, 125.0D+00, 120.0D+00, 110.0D+00, &
     100.0D+00,   0.0D+00,  40.0D+00,  65.0D+00,  60.0D+00, &
     125.0D+00,  40.0D+00,   0.0D+00,  45.0D+00,  55.0D+00, &
     120.0D+00,  65.0D+00,  45.0D+00,   0.0D+00,  50.0D+00, &
     110.0D+00,  60.0D+00,  55.0D+00,  50.0D+00,   0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055'
  write ( *, '(a)' ) '  For a graph defined by a distance matrix,'
  write ( *, '(a)' ) '  GRAPH_DIST_MIN_SPAN_TREE finds a minimum spanning tree.'
  write ( *, '(a)' ) ' '
 
  call graph_dist_print ( dist, lda, nnode, '  The graph:' )

  call graph_dist_min_span_tree ( lda, nnode, dist, itree, jtree )
 
  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do
 
  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The minimal spanning tree:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )
 
  return
end
subroutine test056

!*****************************************************************************80
!
!! TEST056 tests GRAPH_DIST_MIN_SPAN_TREE2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: lda = nnode

  integer ( kind = 4 ) class(nnode)
  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtree(nnode-1)
  real ( kind = 8 ) wtree(nnode-1)

  data ( ( dist(i,j), i = 1, nnode ), j = 1, nnode ) / &
       0.0, 100.0, 125.0, 120.0, 110.0, &
     100.0,   0.0,  40.0,  65.0,  60.0, &
     125.0,  40.0,   0.0,  45.0,  55.0, &
     120.0,  65.0,  45.0,   0.0,  50.0, &
     110.0,  60.0,  55.0,  50.0,   0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST056'
  write ( *, '(a)' ) '  For a graph defined by a distance matrix,'
  write ( *, '(a)' ) '  GRAPH_DIST_MIN_SPAN_TREE2 finds a minimum spanning tree.'
  write ( *, '(a)' ) ' '
 
  call graph_dist_print ( dist, lda, nnode, '  The graph:' )

  call graph_dist_min_span_tree2 ( lda, nnode, dist, class, itree, jtree )
 
  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do

  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The minimal spanning tree:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )
 
  return
end
subroutine test057

!*****************************************************************************80
!
!! TEST057 tests GRAPH_DIST_MIN_SPAN_TREE3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) jtree(nnode-1)
  integer ( kind = 4 ) j
  real ( kind = 8 ) wtree(nnode-1)

  data ( ( dist(i,j), i = 1, nnode ), j = 1, nnode ) / &
       0.0, 100.0, 125.0, 120.0, 110.0, &
     100.0,   0.0,  40.0,  65.0,  60.0, &
     125.0,  40.0,   0.0,  45.0,  55.0, &
     120.0,  65.0,  45.0,   0.0,  50.0, &
     110.0,  60.0,  55.0,  50.0,   0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST057'
  write ( *, '(a)' ) '  For a graph defined by a distance matrix,'
  write ( *, '(a)' ) '  GRAPH_DIST_MIN_SPAN_TREE3 finds a minimum spanning tree.'
  write ( *, '(a)' ) ' '
 
  call graph_dist_print ( dist, lda, nnode, '  The graph:' )

  call graph_dist_min_span_tree3 ( lda, nnode, dist, itree, jtree )

  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do

  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The minimal spanning tree:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )
 
  return
end
subroutine test058

!*****************************************************************************80
!
!! TEST058 tests GRAPH_DIST_MIN_SPAN_TREE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 57
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dist(lda,nnode)
  character ( len = 80 ) :: file_name = '57_city_distances.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jtree(nnode-1)
  real ( kind = 8 ) wtree(nnode-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST058'
  write ( *, '(a)' ) '  GRAPH_DIST_MIN_SPAN_TREE finds a minimum '
  write ( *, '(a)' ) '    spanning tree.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read distance data for 57 cities from file.'
!
!  Read the data.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Problems opening the file: ' // trim ( file_name )
    write ( *, '(a)' ) '  The test was abandoned.'
    return
  end if

  do i = 1, nnode

    read ( iunit, *, iostat = ios ) dist(i,1:nnode)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problems reading the data.'
      write ( *, '(a)' ) '  The test was abandoned.'
      return
    end if

  end do

  close ( unit = iunit )
!
!  Compute the tree.
!
  call graph_dist_min_span_tree ( lda, nnode, dist, itree, jtree )
 
  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do
 
  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The weighted tree:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )

  return
end
subroutine test059

!*****************************************************************************80
!
!! TEST059 tests GRAPH_DIST_ONE.
!
!  Discussion:
!
!    This example appears on page 15 of the reference book by Gibbons.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dinfin
  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idad(nnode)
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) path(nnode)
  integer ( kind = 4 ) itemp(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) length
  real ( kind = 8 ) path_dist(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059'
  write ( *, '(a)' ) '  GRAPH_DIST_ONE computes the distance from one'
  write ( *, '(a)' ) '    node to all others in a graph.'
  write ( *, '(a)' ) ' '

  dinfin = 1000.0D+00
 
  do i = 1, nnode
    do j = 1, nnode
      dist(i,j) = dinfin
    end do
    dist(i,i) = 0.0D+00
  end do
 
  dist(1,2) = 1.0D+00
  dist(1,3) = 3.0D+00
 
  dist(2,1) = 2.0D+00
  dist(2,3) = 1.0D+00
  dist(2,5) = 2.0D+00
 
  dist(3,4) = 2.0D+00
  dist(3,5) = 3.0D+00
 
  dist(4,3) = 1.0D+00
 
  dist(5,1) = 1.0D+00
  dist(5,2) = 3.0D+00
  dist(5,4) = 6.0D+00

  call graph_dist_print ( dist, lda, nnode, '  Edge Distance Matrix:' )

  inode = 5
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'The starting node is ', inode
  write ( *, '(a)' ) ' '

  call graph_dist_one ( dist, dinfin, path_dist, idad, inode, path, &
    lda, nnode )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node    Distance   Path Idad'
  write ( *, '(a)' ) ' '
 
  do i = 1, nnode
    write ( *, '(i5,g14.6,2i5)' ) i, path_dist(i), path(i), idad(i)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Note that "infinity" is represented by ', dinfin
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here are the paths for each node:'
  write ( *, '(a)' ) ' '
 
  do i = 1, nnode

    length = 1
    itemp(length) = i
 
    do while ( itemp(length) /= inode )
      length = length+1
      itemp(length) = idad(itemp(length-1))
    end do
 
    write ( *, '(5i5)' ) itemp(1:length)
 
  end do
 
  return
end
subroutine test060

!*****************************************************************************80
!
!! TEST060 tests VLA_TO_GRAPH_ARC, GRAPH_ARC_FACE, FACE_TO_IV;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxedge = 1000
  integer ( kind = 4 ), parameter :: maxface = 2000
  integer ( kind = 4 ), parameter :: maxnode = 1000
  integer ( kind = 4 ), parameter :: maxorder = 20

  integer ( kind = 4 ) face(maxorder,maxface)
  integer ( kind = 4 ) face_count(maxedge)
  integer ( kind = 4 ) face_order(maxface)
  character ( len = 80 ) :: file_in = 'fish_lines.vla'
  character ( len = 80 ) :: file_out = 'fish_faces.iv'
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface(maxedge)
  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) jface(maxedge)
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nnode
  real ( kind = 8 ) x(maxnode)
  real ( kind = 8 ) y(maxnode)
  real ( kind = 8 ) z(maxnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060'
  write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC converts VLA edge data to a'
  write ( *, '(a)' ) '    graph defined by arcs;'
  write ( *, '(a)' ) '  GRAPH_ARC_FACE constructs the faces of an orientable graph.'
  write ( *, '(a)' ) '  FACE_TO_IV writes face data to an IV file.'
!
!  Get the edge array for the graph.
!
  call vla_to_graph_arc ( file_in, maxedge, maxnode, nedge, nnode, inode, &
    jnode, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) 'TEST060 - Error!'
    write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC returned an error.'
    return
  end if
!
!  Sort the edge array.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sort the edges:'

  call graph_arc_edge_sort ( nedge, inode, jnode )
!
!  Determine the faces.
!
  write ( *, '(a)' ) '  Determine the faces:'

  call graph_arc_face ( face, face_count, face_order, iface, jface, &
    inode, jnode, maxface, maxorder, nedge, nface, nnode )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces found was ', nface
  write ( *, '(a,i8)' ) '  Euler predicted ', nedge + 2 - nnode
!
!  Write the faces to an IV file.
!
  call face_to_iv ( file_out, face, face_order, inode, jnode, &
    nedge, maxnode, maxface, maxorder, nnode, nface, x, y, z )

  return
end
subroutine test061

!*****************************************************************************80
!
!! TEST061 tests GRF_READ, GRAPH_ARC_TO_PS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxedge = 500
  integer ( kind = 4 ), parameter :: maxnode = 100

  integer ( kind = 4 ), parameter :: lda = maxnode

  integer ( kind = 4 ) adj(lda,maxnode)
  character ( len = 80 ) :: file_grf = 'knightstour.grf'
  character ( len = 80 ) :: file_ps = 'knightstour.eps'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  real ( kind = 8 ) x(maxnode)
  real ( kind = 8 ) y(maxnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST061'
  write ( *, '(a)' ) '  GRF_READ reads a GRF file,'
  write ( *, '(a)' ) '  GRAPH_ARC_TO_PS writes a PostScript version of it.'

  call grf_read ( file_grf, inode, jnode, maxedge, maxnode, nedge, nnode, x, y )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Node, X, Y'
  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(i8,2g14.6)' ) i, x(i), y(i)
  end do

  call graph_arc_to_graph_adj ( nedge, inode, jnode, adj, lda, nnode )

  call graph_adj_print ( adj, lda, nnode, '  The graph:' )
!
!  Now write out a PostScript version.
!
  call graph_arc_to_ps ( file_ps, inode, jnode, nedge, nnode, x, y )

  return
end
subroutine test062

!*****************************************************************************80
!
!! TEST062 tests GREEDY.
!
!  Discussion:
!
!    Random data is used in setting up the problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 15

  real ( kind = 8 ) dist
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) nodeb(nnode)
  integer ( kind = 4 ) nodeb1
  integer ( kind = 4 ) noder(nnode)
  integer ( kind = 4 ) noder1
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tol
  real ( kind = 8 ) total
  real ( kind = 8 ) xb(nnode)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xr(nnode)
  real ( kind = 8 ) yb(nnode)
  real ( kind = 8 ) yhi
  real ( kind = 8 ) ylo
  real ( kind = 8 ) yr(nnode)

  seed = 123456789
!
!  IDO just tells us if this is the first or later trials.
!
  ido = 1
!
!  Set the maximum number of iterations.
!
  maxit = 10
!
!  Set the range of the X and Y coordinates.
!
  xhi = 10.0D+00
  xlo = 0.0D+00
  yhi = 5.0D+00
  ylo = 3.0D+00
!
!  Set the relative tolerance for the stepwise distance decrease.
!
  tol = 0.05D+00
!
!  Say hello.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062'

  write ( *, '(a)' ) '  GREEDY tries to minimize the total distance'
  write ( *, '(a)' ) '    in a pairing of black and red nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Try to find a pairing of two sets of nodes'
  write ( *, '(a)' ) '  with a low discrepancy.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Relative tolerance for step decrease = ', tol
  write ( *, '(a,i8)' ) '  Maximum number of steps = ', maxit
  write ( *, '(a,g14.6,a,g14.6)' ) '  X range is ', xlo,' to ', xhi
  write ( *, '(a,g14.6,a,g14.6)' ) '  Y range is ', ylo,' to ', yhi
!
!  Make an arbitrary pairing of the nodes.
!
  do indx = 1, nnode
    nodeb(indx) = indx
    noder(indx) = indx
  end do
!
!  Make up a random set of X, Y coordinates for the nodes.
!
  call r8vec_uniform ( nnode, xlo, xhi, seed, xb )
  call r8vec_uniform ( nnode, xlo, xhi, seed, xr )
  call r8vec_uniform ( nnode, ylo, yhi, seed, yb )
  call r8vec_uniform ( nnode, ylo, yhi, seed, yr )
!
!  We will jump back here if we restart with a permuted NODER.
!
  do ido = 1, 2
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Initial black node coordinates:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I   Black   X             Y'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      write ( *, '(2i8,2g14.6)' ) indx, nodeb(indx), xb(indx), yb(indx)
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Initial red node coordinates:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I    Red    X             Y'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      write ( *, '(2i8,2g14.6)' ) indx, noder(indx), xr(indx), yr(indx)
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Initial pairing of nodes:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I   Black  Red    Distance'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      nodeb1 = nodeb(indx)
      noder1 = noder(indx)
      dist = sqrt ( ( xb(nodeb1) - xr(noder1) )**2 + &
                    ( yb(nodeb1) - yr(noder1) )**2 )

      write ( *, '(3i8,g14.6)' ) indx, nodeb1, noder1, dist
    end do
 
    total = 0.0D+00
    do indx = 1, nnode
      nodeb1 = nodeb(indx)
      noder1 = noder(indx)
      total = total + sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                           + ( yb(nodeb1) - yr(noder1) )**2 )
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) 'Total discrepancy of initial pairing = ', total
!
!  Call GREEDY to seek a better pairing.
!
    call greedy ( maxit, nodeb, noder, nnode, tol, xb, xr, yb, yr )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Final black node coordinates:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I   Black   X             Y'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      write ( *, '(2i8,2g14.6)' ) indx, nodeb(indx), xb(indx), yb(indx)
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Final red node coordinates:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I    Red    X             Y'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      write ( *, '(2i8,2g14.6)' ) indx, noder(indx), xr(indx), yr(indx)
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Final pairing of nodes:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I   Black  Red    Distance'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode

      nodeb1 = nodeb(indx)
      noder1 = noder(indx)

      dist = sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                  + ( yb(nodeb1) - yr(noder1) )**2 )

      write ( *, '(3i8,g14.6)') indx, nodeb1, noder1, dist

    end do
 
    total = 0.0D+00
    do indx = 1, nnode
      nodeb1 = nodeb(indx)
      noder1 = noder(indx)
      dist = sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                  + ( yb(nodeb1) - yr(noder1) )**2 )

      total = total + dist
 
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Total discrepancy of final pairing = ', total
!
!  On the second try, reverse the ordering of the red nodes.
!  Any random permutation would be worth trying.
!
    if ( ido == 1 ) then
 
      do indx = 1, nnode / 2
        call i4_swap ( noder(indx), noder(nnode+1-indx) )
      end do
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Reversing NODER!'
 
    end if

  end do
 
  return
end
subroutine test063

!*****************************************************************************80
!
!! TEST063 tests MAZE_DIAM, MAZE_PATH, MAZE_PRINT, MAZE_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) bar(m,n+1)
  integer ( kind = 4 ) dad(m,n)
  integer ( kind = 4 ) degree(m,n)
  integer ( kind = 4 ) diam
  integer ( kind = 4 ) flat(m+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstart
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) path(m,n)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST063'
  write ( *, '(a)' ) '  MAZE_RANDOM: generate a random maze;'
  write ( *, '(a)' ) '  MAZE_DIAM: find two far apart cells;'
  write ( *, '(a)' ) '  MAZE_PATH: generate a path.'
  write ( *, '(a)' ) '  MAZE_PRINT: print a maze.'
!
!  Print out the cell numbers for the maze.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cell numbers for the maze:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(20i3)' ) ( (j-1)*m+i, j = 1, n )
  end do
!
!  Get a random maze and print it.
!
  call maze_random ( m, n, seed, bar, dad, flat )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A random maze:'
  write ( *, '(a,i8)' ) '    Number of rows =   ', m
  write ( *, '(a,i8)' ) '    Number of columns = ', n

  istart = 0
  jstart = 0

  istop = 0
  jstop = 0

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rooted tree representation:'
  write ( *, '(a)' ) '  (0 is the root.  All other cells print the'
  write ( *, '(a)' ) '  cell number of their parent on the tree.)'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(20i3)' ) dad(i,1:n)
  end do
!
!  Get start and end points that are far apart and print the maze.
!
  call maze_diam ( bar, degree, diam, flat, m, n, path, istart, jstart, &
    istop, jstop )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random maze with far apart ends:'
  write ( *, '(a,i8)' ) '    Diameter = ', diam
  write ( *, '(a,2i8)' ) '    Starting cell = ', istart, jstart
  write ( *, '(a,2i8)' ) '    Stopping cell = ', istop, jstop

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )
!
!  Find a path and print it.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random maze with path from start to stop:'

  call maze_path ( bar, flat, m, n, istart, jstart, istop, jstop )

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze' )

  return
end
subroutine test064

!*****************************************************************************80
!
!! TEST064 tests MAZE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 3
 
  integer ( kind = 4 ), parameter :: INDEF = -1
  integer ( kind = 4 ), parameter :: WALL = 0
  integer ( kind = 4 ), parameter :: OPEN = 1

  integer ( kind = 4 ) bar(m,n+1)
  integer ( kind = 4 ) flat(m+1,n)
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) jstart
  integer ( kind = 4 ) jstop

  bar(1:m,1:n+1) = WALL
  flat(1:m+1,1:n) = WALL

  bar(1,2) = OPEN
  bar(1,4) = INDEF
  bar(2,3) = OPEN

  flat(1,3) = INDEF
  flat(2,1) = OPEN
  flat(2,2) = OPEN
  flat(2,3) = OPEN
  flat(3,1) = OPEN

  istart = 2
  jstart = 1

  istop = 1
  jstop = 3
!
!  Now mark the path.
!
  flat(2,1) = 2
  bar(1,2) = 2
  flat(2,2) = 2
  bar(2,3) = 2
  flat(2,3) = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST064'
  write ( *, '(a)' ) '  MAZE_PRINT prints a maze with path marked.'
  write ( *, '(a)' ) ' '

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )

  return
end
subroutine test065

!*****************************************************************************80
!
!! TEST065 tests NETWORK_FLOW_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: nedge = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icut(nnode)
  integer ( kind = 4 ) icpflo(2,nedge)
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) :: isink = 6
  integer ( kind = 4 ) :: isorce = 1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node_flow(nnode)

  data ( ( iendpt(i,j), j = 1, nedge ), i = 1, 2 ) / &
    1,2, 1,3, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5, 4,6, 5,6, &
    2,1, 3,1, 3,2, 4,2, 5,2, 4,3, 5,3, 5,4, 6,4, 6,5 /
 
  data ( ( icpflo(i,j), j = 1, nedge ), i = 1, 2 ) / &
    3,0,7,0,2,0,5,0,4,0,1,0,4,0,2,0,8,0,3,0, &
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  NETWORK_FLOW_MAX finds the maximum flow on a network.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The source is node ', isorce
  write ( *, '(a,i8)' ) '  The sink is node   ', isink
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Endpoint array:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) ( iendpt(1,i), i = 1, nedge )
  write ( *, '(20i3)' ) ( iendpt(2,i), i = 1, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input edge capacity array:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) ( icpflo(1,i), i = 1, nedge)
 
  call network_flow_max ( nnode, nedge, iendpt, icpflo, isorce, &
    isink, icut, node_flow )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reordered endpoint array:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) ( iendpt(1,i), i = 1, nedge )
  write ( *, '(20i3)' ) ( iendpt(2,i), i = 1, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output edge capacity/flow array:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) ( icpflo(1,i), i = 1, nedge )
  write ( *, '(20i3)' ) ( icpflo(2,i), i = 1, nedge )

  call i4vec_print ( nnode, icut, '  Minimal node cut vector:' )

  call i4vec_print ( nnode, node_flow, '  Nodal flow vector:' )

  return
end
subroutine test066

!*****************************************************************************80
!
!! TEST066 tests NODE_RELAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_cor3 = 100
  integer ( kind = 4 ), parameter :: max_face = 100
  integer ( kind = 4 ), parameter :: max_order = 5

  real ( kind = 8 ) cor3(3,max_cor3)
  real ( kind = 8 ) cor3_new(3,max_cor3)
  integer ( kind = 4 ) cor3_num(max_cor3)
  integer ( kind = 4 ) face(max_order,max_face)
  integer ( kind = 4 ) face_order(max_face)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) num_cor3
  integer ( kind = 4 ) num_face

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST066'
  write ( *, '(a)' ) '  NODE_RELAX smooths a surface.'

  num_cor3 = 8

  cor3(1,1) = 0.0D+00
  cor3(2,1) = 0.0D+00
  cor3(3,1) = 0.0D+00

  cor3(1,2) = 1.0D+00
  cor3(2,2) = 0.0D+00
  cor3(3,2) = 0.0D+00

  cor3(1,3) = 1.0D+00
  cor3(2,3) = 1.0D+00
  cor3(3,3) = 0.0D+00

  cor3(1,4) = 0.0D+00
  cor3(2,4) = 1.0D+00
  cor3(3,4) = 0.0D+00

  cor3(1,5) = 0.0D+00
  cor3(2,5) = 0.0D+00
  cor3(3,5) = 1.0D+00

  cor3(1,6) = 1.0D+00
  cor3(2,6) = 0.0D+00
  cor3(3,6) = 1.0D+00

  cor3(1,7) = 1.0D+00
  cor3(2,7) = 1.0D+00
  cor3(3,7) = 1.0D+00

  cor3(1,8) = 0.0D+00
  cor3(2,8) = 1.0D+00
  cor3(3,8) = 1.0D+00

  num_face = 6

  face(1,1) = 1
  face(2,1) = 4
  face(3,1) = 3
  face(4,1) = 2

  face(1,2) = 2
  face(2,2) = 6
  face(3,2) = 7
  face(4,2) = 3

  face(1,3) = 3
  face(2,3) = 7
  face(3,3) = 8
  face(4,3) = 4

  face(1,4) = 4
  face(2,4) = 8
  face(3,4) = 5
  face(4,4) = 1

  face(1,5) = 1
  face(2,5) = 5
  face(3,5) = 6
  face(4,5) = 2

  face(1,6) = 5
  face(2,6) = 8
  face(3,6) = 7
  face(4,6) = 6

  face_order(1:num_face) = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Old coordinates'
  write ( *, '(a)' ) ' '
  do j = 1, num_cor3
    write ( *, '(i4, 3g14.6)' ) j, cor3(1:3,j)
  end do

  call node_relax ( cor3, cor3_new, cor3_num, face, face_order, max_cor3, &
    max_face, max_order, num_cor3, num_face )

  write ( *, '(a)' ) ' '
  write ( *, '(a)') '  After 1 step'
  write ( *, '(a)' ) ' '

  do j = 1, num_cor3
    write ( *, '(i4, 3g14.6)' ) j, cor3_new(1:3,j)
  end do

  cor3(1:3,1:num_cor3) = cor3_new(1:3,1:num_cor3)

  call node_relax ( cor3, cor3_new, cor3_num, face, face_order, max_cor3, &
    max_face, max_order, num_cor3, num_face )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After 2 steps'
  write ( *, '(a)' ) ' '

  do j = 1, num_cor3
    write ( *, '(i4, 3g14.6)' ) j, cor3_new(1:3,j)
  end do

  cor3(1:3,1:num_cor3) = cor3_new(1:3,1:num_cor3)

  call node_relax ( cor3, cor3_new, cor3_num, face, face_order, max_cor3, &
    max_face, max_order, num_cor3, num_face )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After 3 steps'
  write ( *, '(a)' ) ' '

  do j = 1, num_cor3
    write ( *, '(i4, 3g14.6)' ) j, cor3_new(1:3,j)
  end do

  return
end
subroutine test0665

!*****************************************************************************80
!
!! TEST0665 tests PERM_INC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) perm(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0665'
  write ( *, '(a)' ) '  PERM_INC increments a permutation.'
  write ( *, '(a)' ) ' '

  i = 0
  ipos = 0

  do

    call perm_inc ( perm, ipos, n )

    if ( ipos == 0 ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,4i2)' ) i, perm(1:n)

  end do

  return
end
subroutine test067

!*****************************************************************************80
!
!! TEST067 tests POLY_TO_TRI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_face = 20
  integer ( kind = 4 ), parameter :: max_vert = 5

  integer ( kind = 4 ) face(max_vert,max_face)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) num_face
  integer ( kind = 4 ) num_vert(max_face)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST067'
  write ( *, '(a)' ) '  POLY_TO_TRI replaces a polygonal mesh with a'
  write ( *, '(a)' ) '    triangular one.'

  num_face = 4

  num_vert(1) = 4
  face(1,1) = 1
  face(2,1) = 3
  face(3,1) = 5
  face(4,1) = 7

  num_vert(2) = 3
  face(1,2) = 2
  face(2,2) = 3
  face(3,2) = 9

  num_vert(3) = 5
  face(1,3) = 3
  face(2,3) = 7
  face(3,3) = 8
  face(4,3) = 23
  face(5,3) = 2

  num_vert(4) = 4
  face(1,4) = 4
  face(2,4) = 7
  face(3,4) = 8
  face(4,4) = 23

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces = ', num_face

  call i4vec_print ( num_face, num_vert, '  Faces and number of vertices:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face   Vertices'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(6i8)' ) i, ( face(j,i), j = 1, num_vert(i) )
  end do

  call poly_to_tri ( face, ierror, max_face, max_vert, num_face, num_vert )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The algorithm failed.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of faces = ', num_face

    call i4vec_print ( num_face, num_vert, '  Faces and number of vertices:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Face   Vertices'
    write ( *, '(a)' ) ' '
    do i = 1, num_face
      write ( *, '(6i8)' ) i, ( face(j,i), j = 1, num_vert(i) )
    end do

  end if

  return
end
subroutine test068

!*****************************************************************************80
!
!! TEST068 tests PRUEFER_TO_TREE_ARC.
!
!  The tree is
!
!          5
!          |
!    2-3-6-8-1-9
!      |   |
!      7   4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 9

  integer ( kind = 4 ), save, dimension ( nnode-2 ) :: code = (/ 1, 3, 8, 8, 3, 6, 8 /)
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) jnode(nnode-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST068'
  write ( *, '(a)' ) '  PRUEFER_TO_TREE_ARC computes a tree from its Pruefer code.'
 
  call i4vec_print ( nnode-2, code, '  The Pruefer code:' )

  call pruefer_to_tree_arc ( nnode, code, inode, jnode )
 
  call graph_arc_print ( nnode-1, inode, jnode, '  The graph:' )

  return
end
subroutine test069

!*****************************************************************************80
!
!! TEST069 tests PRUEFER_TO_TREE_2.
!
!  The tree is
!
!          5
!          |
!    2-3-6-8-1-9
!      |   |
!      7   4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 9

  integer ( kind = 4 ), save, dimension ( nnode ) :: code = (/ 1, 3, 8, 8, 3, 6, 8, 0, 0 /)
  integer ( kind = 4 ) itree(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST069'
  write ( *, '(a)' ) '  PRUEFER_TO_TREE_2 produces a tree from its Pruefer code'

  call i4vec_print ( nnode-2, code, '  The Pruefer code:' )

  call pruefer_to_tree_2 ( nnode, code, itree )
 
  call i4vec_print ( nnode-1, itree, '  The edge list of the tree:' )
 
  return
end
subroutine test0695

!*****************************************************************************80
!
!! TEST0695 tests VLA_TO_GRAPH_ARC, SHAPE_3D_NODES_TO_PS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: max_edge = 1000
  integer   ( kind = 4 ), parameter :: max_node = 1000

  character ( len = 80 ) :: file_in = 'fish_lines.vla'
  character ( len = 80 ) :: file_out = 'fish_nodes.ps'
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) inode(max_edge)
  integer   ( kind = 4 ) jnode(max_edge)
  integer   ( kind = 4 ) num_edge
  integer   ( kind = 4 ) num_node
  real      ( kind = 8 ) x(max_node)
  real      ( kind = 8 ) y(max_node)
  real      ( kind = 8 ) z(max_node)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0695'
  write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC reads a VLA file and converts it'
  write ( *, '(a)' ) '    to a graph defined by an arc list.'
  write ( *, '(a)' ) '  SHAPE_3D_NODES_TO_PS writes the nodes to a PostScript file.'
!
!  Get the edge array for the graph.
!
  call vla_to_graph_arc ( file_in, max_edge, max_node, num_edge, &
    num_node, inode, jnode, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)') '  VLA_TO_GRAPH_ARC returned an error.'
    return
  end if

  call shape_3d_nodes_to_ps ( file_out, num_node, x, y, z )

  return
end
subroutine test0696

!*****************************************************************************80
!
!! TEST0696 tests VLA_TO_GRAPH_ARC, SHAPE_3D_EDGES_TO_PS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: max_edge = 1000
  integer   ( kind = 4 ), parameter :: max_face = 2000
  integer   ( kind = 4 ), parameter :: max_node = 1000
  integer   ( kind = 4 ), parameter :: max_order = 20

  integer   ( kind = 4 ) face(max_order,max_face)
  integer   ( kind = 4 ) face_count(max_edge)
  integer   ( kind = 4 ) face_order(max_face)
  character ( len = 80 ) :: file_in = 'fish_lines.vla'
  character ( len = 80 ) :: file_out = 'fish_edges.ps'
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iface(max_edge)
  integer   ( kind = 4 ) inode(max_edge)
  integer   ( kind = 4 ) jface(max_edge)
  integer   ( kind = 4 ) jnode(max_edge)
  integer   ( kind = 4 ) num_edge
  integer   ( kind = 4 ) num_face
  integer   ( kind = 4 ) num_node
  real      ( kind = 8 ) x(max_node)
  real      ( kind = 8 ) y(max_node)
  real      ( kind = 8 ) z(max_node)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0696'
  write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC reads a VLA file and converts it'
  write ( *, '(a)' ) '    to a graph defined by an arc list.'
  write ( *, '(a)' ) '  SHAPE_3D_EDGES_TO_PS writes the edges to a PostScript file.'
!
!  Get the edge array for the graph.
!
  call vla_to_graph_arc ( file_in, max_edge, max_node, num_edge, &
    num_node, inode, jnode, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC returned an error.'
    return
  end if
!
!  Sort the edge array.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sort the edges:'

  call graph_arc_edge_sort ( num_edge, inode, jnode )
!
!  Determine the faces.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Determine the faces:'

  call graph_arc_face ( face, face_count, face_order, iface, jface, inode, &
    jnode, max_face, max_order, num_edge, num_face, num_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The faces were determined.'
  write ( *, '(a,i8)' ) '    Number of faces found was ', num_face
  write ( *, '(a,i8)' ) '    Euler predicted ', num_edge + 2 - num_node

  call shape_3d_edges_to_ps ( file_out, max_order, num_face, num_node, &
    face, face_order, x, y, z )

  return
end
subroutine test0697

!*****************************************************************************80
!
!! TEST0697 tests VLA_TO_GRAPH_ARC, SHAPE_3D_FACES_TO_PS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_edge = 1000
  integer ( kind = 4 ), parameter :: max_face = 2000
  integer ( kind = 4 ), parameter :: max_node = 500
  integer ( kind = 4 ), parameter :: max_order = 20

  integer ( kind = 4 ) face(max_order,max_face)
  integer ( kind = 4 ) face_count(max_edge)
  integer ( kind = 4 ) face_order(max_face)
  character ( len = 80 ) :: file_in = 'fish_lines.vla'
  character ( len = 80 ) :: file_out = 'fish_faces.ps'
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface(max_edge)
  integer ( kind = 4 ) inode(max_edge)
  integer ( kind = 4 ) jface(max_edge)
  integer ( kind = 4 ) jnode(max_edge)
  integer ( kind = 4 ) num_edge
  integer ( kind = 4 ) num_face
  integer ( kind = 4 ) num_node
  real ( kind = 8 ) x(max_node)
  real ( kind = 8 ) y(max_node)
  real ( kind = 8 ) z(max_node)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0697'
  write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC reads a VLA file and converts it'
  write ( *, '(a)' ) '    to a graph defined by an arc list.'
  write ( *, '(a)' ) '  SHAPE_3D_FACES_TO_PS writes the faces to a PostScript file.'
!
!  Get the edge array for the graph.
!
  call vla_to_graph_arc ( file_in, max_edge, max_node, num_edge, &
    num_node, inode, jnode, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) 'TEST0697 - Error!'
    write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC returned an error.'
    return
  end if
!
!  Sort the edge array.

  call graph_arc_edge_sort ( num_edge, inode, jnode )
!
!  Determine the faces.
!
  call graph_arc_face ( face, face_count, face_order, iface, jface, inode, &
    jnode, max_face, max_order, num_edge, num_face, num_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces found was ', num_face
  write ( *, '(a,i8)' ) '  Euler predicted ', num_edge + 2 - num_node

  call shape_3d_faces_to_ps ( file_out, max_order, num_face, num_node, &
    face, face_order, x, y, z )

  return
end
subroutine test070

!*****************************************************************************80
!
!! TEST070 tests SORT_HEAP_EXTERNAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST070'
  write ( *, '(a)' ) '  SORT_HEAP_EXTERNAL sorts objects externally.'
  write ( *, '(a)' ) ' '

  indx = 0
  i = 0
  j = 0
  isgn = 0

  call i4vec_uniform ( n, 1, n, seed, a )

  call i4vec_print ( n, a, '  Before sorting:' ) 
 
  do

    call sort_heap_external ( n, indx, i, j, isgn )
 
    if ( indx < 0 ) then
      isgn = 1
      if ( a(i) <= a(j) ) then
        isgn = -1
      end if
    else if ( indx > 0 ) then
      call i4_swap ( a(i), a(j) )
    else
      exit
    end if

  end do

  call i4vec_print ( n, a, '  After sorting:' )
 
  return
end
subroutine test071

!*****************************************************************************80
!
!! TEST071 tests SPAN_FOREST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 14
  integer ( kind = 4 ), parameter :: nedge = 10

  integer ( kind = 4 ) component(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  data ( ( iendpt(i,j), i = 1, 2 ), j = 1, nedge ) / &
    2,3, 4,7, 1,9, 7,11, 5,8, 2,5, 6,10, 2,8, 3,8, 4,11 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST071'
  write ( *, '(a)' ) '  SPAN_FOREST: a spanning forest for a graph'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial end point array:'
  write ( *, '(a)' ) ' '
  write ( *, '(19i4)' ) ( iendpt(1,j), j = 1, nedge )
  write ( *, '(19i4)' ) ( iendpt(2,j), j = 1, nedge )
 
  call span_forest ( nnode, nedge, iendpt, k, component )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reordered endpoint array:'
  write ( *, '(a)' ) ' '
  write ( *, '(19i4)' ) ( iendpt(1,j), j = 1, nedge )
  write ( *, '(19i4)' ) ( iendpt(2,j), j = 1, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of connected components = ', k
 
  call i4vec_print ( nnode, component, '  Node, Component' )

  return
end
subroutine test072

!*****************************************************************************80
!
!! TEST072 tests SPAN_TREE_NEXT;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: nedge = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(nnode-1)
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncan(nnode-1)
  integer ( kind = 4 ) nspan
  integer ( kind = 4 ) signal

  data ( ( iendpt(i,j), i = 1, 2 ), j = 1, nedge ) / &
    1,2, 1,3, 1,4, 1,5, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST072'
  write ( *, '(a)' ) '  SPAN_TREE_NEXT constructs spanning trees'
  write ( *, '(a)' ) '    of a graph using a backtrack search.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Node1   Node2'
  write ( *, '(a)' ) ' '
  do i = 1, nedge
    write ( *, '(3i8)' ) iendpt(1,i), iendpt(2,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Edges in spanning tree:'
  write ( *, '(a)' ) ' '

  nspan = 0
  signal = 0

  do

    call span_tree_next ( signal, nnode, nedge, iendpt, iarray, ncan )

    if ( signal == 0 ) then 
      exit
    end if

    nspan = nspan + 1
    write ( *, '(i4,4x,5i4)' ) nspan, iarray(1:nnode-1)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of spanning trees found was ', nspan

  return
end
subroutine test073

!*****************************************************************************80
!
!! TEST073 tests TREE_ARC_TO_PRUEFER.
!
!  The tree is
!
!          5
!          |
!    2-3-6-8-1-9
!      |   |
!      7   4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 9

  integer ( kind = 4 ) iarray(nnode-2)
  integer ( kind = 4 ), dimension ( nnode - 1 ) :: inode = (/ 2, 3, 3, 6, 8, 8, 8, 1 /)
  integer ( kind = 4 ), dimension ( nnode - 1 ) :: jnode = (/ 3, 7, 6, 8, 4, 5, 1, 9 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST073'
  write ( *, '(a)' ) '  TREE_ARC_TO_PRUEFER: Tree => Pruefer code'

  call graph_arc_print ( nnode-1, inode, jnode, '  The graph:' )
 
  call tree_arc_to_pruefer ( nnode, inode, jnode, iarray )

  call i4vec_print ( nnode-2, iarray, '  The Pruefer code:' )
 
  return
end
subroutine test074

!*****************************************************************************
!
!! TEST074 tests TREE_ARC_CENTER.
!
!  2---3---6---8---1---9
!     /       / \
!    7       5   4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 9

  integer ( kind = 4 ) center(2)
  integer ( kind = 4 ) eccent
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nnode - 1 ) :: inode = (/ 2, 3, 3, 6, 8, 8, 8, 1 /)
  integer ( kind = 4 ), dimension ( nnode - 1 ) :: jnode = (/ 3, 7, 6, 8, 4, 5, 1, 9 /)
  integer ( kind = 4 ) parity

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST074'
  write ( *, '(a)' ) '  TREE_ARC_CENTER computes the center of a tree.'

  call graph_arc_print ( nnode-1, inode, jnode, '  The edge list of the tree:' )

  call tree_arc_center ( nnode, inode, jnode, center, eccent, parity )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Parity = ', parity
  write ( *, '(a,i8)' ) '  Eccentricity is ', eccent

  if ( parity == 0 ) then
    write ( *, '(a)' ) '  No center node (degenerate case).'
  else if ( parity == 1 ) then
    write ( *, '(a,i8)' ) '  Center node: ', center(1)
  else
    write ( *, '(a,2i8)' ) '  Center nodes: ', center(1), center(2)
  end if

  return
end
subroutine test075

!*****************************************************************************
!
!! TEST075 tests TREE_ARC_DIAM.
!
!  2---3---6---8---1---9
!     /       / \
!    7       5   4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: nnode = 9

  integer ( kind = 4 ) diam
  integer ( kind = 4 ), dimension ( nnode-1 ) :: inode = (/ 2, 3, 3, 6, 8, 8, 8, 1 /)
  integer ( kind = 4 ), dimension ( nnode-1 ) :: jnode = (/ 3, 7, 6, 8, 4, 5, 1, 9 /)
  integer ( kind = 4 ) label(nnode)
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST075'
  write ( *, '(a)' ) '  TREE_ARC_DIAM computes the diameter of a tree.'

  call graph_arc_print ( nnode-1, inode, jnode, '  The edge list of the tree:' )

  call tree_arc_diam ( nnode, inode, jnode, diam, label, nnode1, nnode2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  This tree has a diameter of ', diam
  write ( *, '(a,i8,a,i8)' ) '  between nodes ', nnode1, ' and ', nnode2

  call i4vec_print ( nnode, label, '  Nodes and labels:' )

  return
end
subroutine test076

!*****************************************************************************80
!
!! TEST076 tests TREE_ARC_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icode(nnode-2)
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076'
  write ( *, '(a)' ) '  TREE_ARC_RANDOM produces a random labeled'
  write ( *, '(a)' ) '    tree and its Pruefer code.'
  write ( *, '(a)' ) ' '
 
  do i = 1, 5

    call tree_arc_random ( nnode, seed, icode, inode, jnode )

    call graph_arc_print ( nnode-1, inode, jnode, '  The random tree:' )

    call i4vec_print ( nnode-2, icode, '  The Pruefer code:' )

  end do
 
  return
end
subroutine test077

!*****************************************************************************80
!
!! TEST077 tests TREE_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST077'
  write ( *, '(a)' ) '  TREE_ENUM enumerates the labeled trees on a given'
  write ( *, '(a)' ) '  number of nodes.'
  write ( *, '(a)' ) ' '

  do nnode = 0, 10

    call tree_enum ( nnode, num )

    write ( *, '(i8,i10)' ) nnode, num

  end do
 
  return
end
subroutine test078

!*****************************************************************************80
!
!! TEST078 tests TREE_PARENT_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 4

  integer ( kind = 4 ) iarray(nnode)
  integer ( kind = 4 ) icode(nnode)
  integer ( kind = 4 ) itree(nnode)
  logical more

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078'
  write ( *, '(a)' ) '  TREE_PARENT_NEXT finds all labeled trees of a given '
  write ( *, '(a)' ) '    order, and their Pruefer codes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Pruefer code     Tree'
  write ( *, '(a)' ) ' '
 
  more = .false.
 
  do
 
    call tree_parent_next ( nnode, iarray, icode, itree, more )
 
    write ( *, '(2i2,14x,3i2)' ) icode(1:nnode-2), itree(1:nnode-1)

    if ( .not. more ) then
      exit
    end if

  end do
 
  return
end
subroutine test079

!*****************************************************************************80
!
!! TEST079 tests TREE_RB_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST079'
  write ( *, '(a)' ) '  TREE_RB_ENUM enumerates the rooted binary trees on a '
  write ( *, '(a)' ) '    given number of nodes.'
  write ( *, '(a)' ) ' '

  do nnode = 0, 11

    call tree_rb_enum ( nnode, num )

    write ( *, '(2x,i8,2x,i8)' ) nnode, num

  end do
 
  return
end
subroutine test080

!*****************************************************************************80
!
!! TEST080 tests TREE_RB_LEX_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  logical more

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST080'
  write ( *, '(a)' ) '  TREE_RB_LEX_NEXT produces all rooted binary trees with'
  write ( *, '(a)' ) '  a given number of nodes, in lexicographic order, using'
  write ( *, '(a)' ) '  the preorder traversal representation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes N = ', n
  write ( *, '(a)' ) ' '

  more = .false.
  i = 0

  do

    call tree_rb_lex_next ( n, a, more )

    if ( .not. more ) then 
      exit
    end if

    i = i + 1
    write ( *, '(i2,2x,11i1)' ) i, a(1:11)

  end do

  return
end
subroutine test081

!*****************************************************************************80
!
!! TEST081 tests TREE_RB_LEX_NEXT, TREE_RB_TO_PARENT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) parent(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST081'
  write ( *, '(a)' ) '  TREE_RB_LEX_NEXT produces all rooted binary trees with'
  write ( *, '(a)' ) '    a given number of nodes, in lexicographic order,'
  write ( *, '(a)' ) '    using the preorder traversal representation.'
  write ( *, '(a)' ) '  TREE_RB_TO_PARENT converts the preorder traversal form'
  write ( *, '(a)' ) '    to the more comprehensible parent node representation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes N = ', n
  write ( *, '(a)' ) ' '

  more = .false.
  i = 0

  do

    call tree_rb_lex_next ( n, a, more )

    if ( .not. more ) then 
      exit
    end if

    call tree_rb_to_parent ( n, a, parent )

    i = i + 1
    write ( *, '(i2,2x,11i3)' ) i, parent(1:n)

  end do

  return
end
subroutine test082

!*****************************************************************************80
!
!! TEST082 tests TREE_RB_YULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  integer ( kind = 4 ) a(n_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST082'
  write ( *, '(a)' ) '  TREE_RB_YULE carries out one step of the Yule model'
  write ( *, '(a)' ) '    on a rooted binary tree stored in preorder traversal form.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each call adds two children to an arbitary leaf.'

  do i = 1, 5

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Simulation ', i
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Nodes  Preorder code'
    write ( *, '(a)' ) ' '

    n = 0

    do

      call tree_rb_yule ( n, seed, a )

      write ( *, '(i2,2x,11i1)' ) n, a(1:n)

      if ( n + 2 > n_max ) then
        exit
      end if

    end do

  end do

  return
end
subroutine test083

!*****************************************************************************80
!
!! TEST083 tests TREE_ROOTED_CODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 12

  integer ( kind = 4 ) code(2*nnode)
  integer ( kind = 4 ), dimension ( nnode ) :: parent = &
    (/ 0, 1, 1, 2, 2, 2, 3, 3, 5, 5, 6, 10 /) 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST083'
  write ( *, '(a)' ) '  TREE_ROOTED_CODE: code of a rooted tree.'
  write ( *, '(a)' ) ' '

  call i4vec_print ( nnode, parent, '  Parent vector for tree:' )

  call tree_rooted_code ( nnode, parent, code )

  call i4vec_print ( 2*nnode, code, '  The tree code:' )

  return
end
subroutine test084

!*****************************************************************************80
!
!! TEST084 tests TREE_ROOTED_DEPTH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 12

  integer ( kind = 4 ) depth
  integer ( kind = 4 ) depth_node(nnode)
  integer ( kind = 4 ), dimension ( nnode ) :: parent = &
    (/ 0, 1, 1, 2, 2, 2, 3, 3, 5, 5, 6, 10 /) 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST084'
  write ( *, '(a)' ) '  TREE_ROOTED_DEPTH: depth of a rooted tree.'
  write ( *, '(a)' ) ' '

  call i4vec_print ( nnode, parent, '  Parent vector for tree:' )

  call tree_rooted_depth ( nnode, parent, depth, depth_node )

  call i4vec_print ( nnode, depth_node, '  Individual node depths:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Overall rooted tree depth: ', depth

  return
end
subroutine test085

!*****************************************************************************80
!
!! TEST085 tests TREE_ROOTED_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) itree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntree(nnode)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  TREE_ROOTED_RANDOM: random unlabeled rooted trees.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random trees, rooted at 1'
 
  do i = 1, 5

    call tree_rooted_random ( nnode, seed, ntree, itree )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Endpoint array for tree:'
    write ( *, '(19i4)' ) ( j, j = 2, nnode )
    write ( *, '(19i4)' ) itree(2:nnode)

  end do
 
  call i4vec_print ( nnode, ntree, &
    '  Number of trees with given number of nodes:' )
 
  return
end
subroutine test086 ( )

!*****************************************************************************80
!
!! TEST086 tests TREE_ROOTED_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 10

  integer ( kind = 4 ) ntree(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST086'
  write ( *, '(a)' ) '  TREE_ROOTED_ENUM counts unlabeled rooted trees.'

  call tree_rooted_enum ( nnode, ntree )

  call i4vec_print ( nnode, ntree, &
    '  Number of trees with given number of nodes:' )

  return
end
