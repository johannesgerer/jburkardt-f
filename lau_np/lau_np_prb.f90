program main

!*****************************************************************************80
!
!! MAIN is the main program for LAU_NP_PRB.
!
!  Discussion:
!
!    LAU_NP_PRB tests the LAU_NP package.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAU_NP_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LAU_NP library.'

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
  call test13 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAU_NP_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests GRAPH_ARC_EULER_CIRC.
!
!  Modified:
!
!    09 February 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 10
  integer ( kind = 4 ), parameter :: nnode = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)
  integer ( kind = 4 ), dimension ( nedge ) :: loop

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a graph described by an arc list:'
  write ( *, '(a)' ) '  GRAPH_ARC_EULER_CIRC finds an'
  write ( *, '(a)' ) '  Euler circuit of a graph.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The graph has order NNODE = ', nnode

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_euler_circ ( nnode, nedge, inode, jnode, loop )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The nodes in the Euler circuit:'
  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(2x,i8,2x,i8)' ) i, loop(i)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests GRAPH_ARC_MIN_PATH.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: nedge = 6

  real ( kind = 8 ), dimension ( nedge ) :: cost = (/ &
    1.0D+00, 1.0D+00, 3.0D+00, 2.0D+00, 2.0D+00, 5.0D+00 /)
  real ( kind = 8 ) dist(nnode,nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ &
    1, 1, 2, 2, 3, 3 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istop
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ &
    2, 3, 3, 5, 4, 5 /)
  integer ( kind = 4 ) num_path
  integer ( kind = 4 ) path(nnode)
  real ( kind = 8 ) path_length
  integer ( kind = 4 ) start

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  GRAPH_ARC_MIN_PATH computes the shortest path from one'
  write ( *, '(a)' ) '  node to another.'
  write ( *, '(a)' ) ' '

  call graph_arc_weight_print ( nedge, inode, jnode, cost, &
    '  The weighted graph:' )

  dist(1:nnode,1:nnode) = 0.0D+00

  do start = 1, nnode
    do istop = start+1, nnode
      call graph_arc_min_path ( nnode, nedge, inode, jnode, cost, start, &
        istop, num_path, path, path_length )
      dist(start,istop) = path_length
      dist(istop,start) = path_length
    end do
  end do

  call graph_dist_print ( nnode, dist, &
    '  The distance matrix constructed by GRAPH_ARC_MIN_PATH:' )

  start = 4
  istop = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)') '  The routine actually also computes the path.'
  write ( *, '(a,i8)' ) '  For instance, starting at node ', start
  write ( *, '(a,i8)' ) '  we compute the shortest path to node ', istop

  call graph_arc_min_path ( nnode, nedge, inode, jnode, cost, start, &
    istop, num_path, path, path_length )

  write ( *, '(a)' ) ' '

  do i = 1, num_path
    write ( *, '(i8,2x,i8)' ) i, path(i)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests GRAPH_ARC_MIN_SPAN_TREE.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 10
  integer ( kind = 4 ), parameter :: nnode = 5

  real ( kind = 8 ), dimension ( nedge ) :: cost = &
    (/ 100.0D+00, 125.0D+00, 120.0D+00, 110.0D+00, 40.0D+00, &
        65.0D+00, 60.0D+00, 45.0D+00, 55.0D+00, 50.0D+00 /)
  real ( kind = 8 ), dimension ( nnode-1) :: ctree
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)
  integer ( kind = 4 ) jtree(nnode-1)
  real ( kind = 8 ) tree_cost

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  GRAPH_ARC_MIN_SPAN_TREE finds a minimum '
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
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests GRAPH_DIST_MIN_SPAN_TREE3.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5

  real ( kind = 8 ), dimension ( nnode,nnode) :: dist = reshape ( (/ &
       0.0D+00, 100.0D+00, 125.0D+00, 120.0D+00, 110.0D+00, &
     100.0D+00,   0.0D+00,  40.0D+00,  65.0D+00,  60.0D+00, &
     125.0D+00,  40.0D+00,   0.0D+00,  45.0D+00,  55.0D+00, &
     120.0D+00,  65.0D+00,  45.0D+00,   0.0D+00,  50.0D+00, &
     110.0D+00,  60.0D+00,  55.0D+00,  50.0D+00,   0.0D+00 /), &
    (/ nnode, nnode /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) jtree(nnode-1)
  integer ( kind = 4 ) j
  real ( kind = 8 ) wtree(nnode-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a graph defined by a distance matrix,'
  write ( *, '(a)' ) '  GRAPH_DIST_MIN_SPAN_TREE3 finds a minimum '
  write ( *, '(a)' ) '  spanning tree.'
  write ( *, '(a)' ) ' '

  call graph_dist_print ( nnode, dist, '  The graph:' )

  call graph_dist_min_span_tree3 ( nnode, dist, itree, jtree )

  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do

  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The minimal spanning tree:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests INT_LP
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(m+1,n+m+1)
  real ( kind = 8 ), dimension ( m ) :: b = (/ 1.0D+00, -1.0D+00, 4.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: c = (/ &
    -2.0D+00, 1.0D+00, 4.0D+00, -1.0D+00, -3.0D+00 /)
  integer ( kind = 4 ) infs
  real ( kind = 8 ) x(n)

  a(1,1) = - 3.0D+00
  a(1,2) = - 1.0D+00
  a(1,3) =   2.0D+00
  a(1,4) =   3.0D+00
  a(1,5) = - 3.0D+00

  a(2,1) =   0.0D+00
  a(2,2) =   1.0D+00
  a(2,3) = - 1.0D+00
  a(2,4) = - 4.0D+00
  a(2,5) = - 2.0D+00

  a(3,1) =   1.0D+00
  a(3,2) =   0.0D+00
  a(3,3) =   4.0D+00
  a(3,4) =   3.0D+00
  a(3,5) =   0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  INT_LP is a heuristic algorithm for the'
  write ( *, '(a)' ) '  integer ( kind = 4 ) linear programming problem.'

  call int_lp ( m, n, a, b, c, x, infs )

  if ( infs == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The problem is infeasible.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  COMPUTED SOLUTION:'
    write ( *, '(a)' ) ' '
    write ( *, '(5g12.4)' ) x(1:n)
  end if

  x(1:5) = (/ 0.0D+00, 2.0D+00, 1.0D+00, 0.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CORRECT SOLUTION:'
  write ( *, '(a)' ) ' '
  write ( *, '(5g12.4)' ) x(1:n)

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests K_CENTER.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 10
  integer ( kind = 4 ), parameter :: m = ( nnode * ( nnode - 1 ) ) / 2

  real ( kind = 8 ) cost(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) knum
  integer ( kind = 4 ) kset(nnode)

  data ( ( cost(i,j), j = 1, nnode ), i = 1, nnode ) / &
     0.0D+00, 15.0D+00, 72.0D+00, 51.0D+00, 50.0D+00, &
    59.0D+00, 53.0D+00, 68.0D+00, 11.0D+00, 33.0D+00, &
    15.0D+00,  0.0D+00, 66.0D+00, 44.0D+00, 43.0D+00, &
    45.0D+00, 56.0D+00, 65.0D+00,  9.0D+00, 35.0D+00, &
    72.0D+00, 66.0D+00,  0.0D+00,104.0D+00, 23.0D+00, &
    77.0D+00, 38.0D+00, 11.0D+00, 62.0D+00, 44.0D+00, &
    51.0D+00, 44.0D+00,104.0D+00,  0.0D+00, 82.0D+00, &
    41.0D+00, 99.0D+00,106.0D+00, 52.0D+00, 79.0D+00, &
    50.0D+00, 43.0D+00, 23.0D+00, 82.0D+00,  0.0D+00, &
    59.0D+00, 26.0D+00, 25.0D+00, 39.0D+00, 28.0D+00, &
    59.0D+00, 45.0D+00, 77.0D+00, 41.0D+00, 59.0D+00, &
     0.0D+00, 82.0D+00, 83.0D+00, 52.0D+00, 70.0D+00, &
    53.0D+00, 56.0D+00, 38.0D+00, 99.0D+00, 26.0D+00, &
    82.0D+00,  0.0D+00, 28.0D+00, 47.0D+00, 21.0D+00, &
    68.0D+00, 65.0D+00, 11.0D+00,106.0D+00, 25.0D+00, &
    83.0D+00, 28.0D+00,  0.0D+00, 59.0D+00, 37.0D+00, &
    11.0D+00,  9.0D+00, 62.0D+00, 52.0D+00, 39.0D+00, &
    52.0D+00, 47.0D+00, 59.0D+00,  0.0D+00, 27.0D+00, &
    33.0D+00, 35.0D+00, 44.0D+00, 79.0D+00, 28.0D+00, &
    70.0D+00, 21.0D+00, 37.0D+00, 27.0D+00,  0.0D+00 /

  kmax = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  K_CENTER is a heuristic algorithm'
  write ( *, '(a)' ) '  for the K center problem.'

  call k_center ( nnode, cost, kmax, knum, kset )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COMPUTED RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The centers: '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) kset(1:knum)

  knum = 4
  kset(1:4) = (/ 9, 5, 4, 6 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CORRECT RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The centers: '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) kset(1:knum)

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests K_MEDIAN.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 15

  real ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isol(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  data ( ( c(i,j), j = 1, n ), i = 1, m ) / &
    2.0D+00, 3.0D+00, 0.0D+00, 6.0D+00, 5.0D+00, &
    2.0D+00, 3.0D+00, 2.0D+00, 4.0D+00, 9.0D+00, &
    2.0D+00, 0.0D+00, 8.0D+00, 7.0D+00, 3.0D+00, &
    4.0D+00, 8.0D+00, 6.0D+00, 0.0D+00, 8.0D+00, &
    1.0D+00, 3.0D+00, 3.0D+00, 4.0D+00, 1.0D+00, &
    1.0D+00, 9.0D+00, 3.0D+00, 8.0D+00, 3.0D+00, &
    0.0D+00, 8.0D+00, 4.0D+00, 6.0D+00, 2.0D+00, &
    3.0D+00, 2.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
    5.0D+00, 7.0D+00, 9.0D+00, 9.0D+00, 0.0D+00, &
    9.0D+00, 7.0D+00, 2.0D+00, 6.0D+00, 4.0D+00, &
    2.0D+00, 1.0D+00, 5.0D+00, 9.0D+00, 0.0D+00, &
    7.0D+00, 1.0D+00, 1.0D+00, 4.0D+00, 2.0D+00, &
    7.0D+00, 4.0D+00, 5.0D+00, 3.0D+00, 3.0D+00, &
    4.0D+00, 5.0D+00, 0.0D+00, 3.0D+00, 4.0D+00, &
    3.0D+00, 3.0D+00, 8.0D+00, 4.0D+00, 9.0D+00, &
    6.0D+00, 1.0D+00, 1.0D+00, 7.0D+00, 7.0D+00, &
    9.0D+00, 7.0D+00, 4.0D+00, 6.0D+00, 7.0D+00, &
    2.0D+00, 2.0D+00, 1.0D+00, 5.0D+00, 0.0D+00, &
    7.0D+00, 4.0D+00, 0.0D+00, 7.0D+00, 7.0D+00, &
    4.0D+00, 3.0D+00, 2.0D+00, 4.0D+00, 3.0D+00, &
    9.0D+00, 5.0D+00, 1.0D+00, 8.0D+00, 5.0D+00, &
    5.0D+00, 1.0D+00, 5.0D+00, 7.0D+00, 0.0D+00, &
    8.0D+00, 4.0D+00, 6.0D+00, 5.0D+00, 6.0D+00, &
    4.0D+00, 3.0D+00, 5.0D+00, 2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, 5.0D+00, 2.0D+00, 4.0D+00, &
    7.0D+00, 4.0D+00, 7.0D+00, 0.0D+00, 9.0D+00, &
    7.0D+00, 5.0D+00, 2.0D+00, 1.0D+00, 7.0D+00, &
    7.0D+00, 9.0D+00, 0.0D+00, 0.0D+00, 6.0D+00, &
    3.0D+00, 0.0D+00, 8.0D+00, 3.0D+00, 9.0D+00, &
    1.0D+00, 7.0D+00, 1.0D+00, 6.0D+00, 5.0D+00 /

  k = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  K_MEDIAN is a heuristic algorithm'
  write ( *, '(a)' ) '  for the K median problem.'

  call k_median ( m, n, c, k, isol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COMPUTED RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The K rows found: '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) isol(1:k)

  isol(1:4) = (/ 5, 2, 6, 4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CORRECT RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The K rows found: '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) isol(1:k)

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests KNAPSACK.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 15
  integer ( kind = 4 ), parameter :: n = 7

  real ( kind = 8 ), dimension ( n ) :: a = (/ &
    41.0D+00, 50.0D+00, 49.0D+00, 59.0D+00, 55.0D+00, 57.0D+00, 68.0D+00 /)
  real ( kind = 8 ) :: b = 170.0D+00
  real ( kind = 8 ), dimension ( n ) :: c = &
    (/ 442.0D+00, 525.0D+00, 511.0D+00, 593.0D+00, 546.0D+00, &
    564.0D+00, 617.0D+00 /)
  real ( kind = 8 ) :: eps = 0.8D+00
  integer ( kind = 4 ) isol(n)
  integer ( kind = 4 ) numsol
  real ( kind = 8 ) objval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  KNAPSACK is a heuristic algorithm'
  write ( *, '(a)' ) '  for the 0/1 knapsack problem.'

  call knapsack ( n, a, b, c, eps, m, objval, numsol, isol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COMPUTED RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Objective function = '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,g14.6)' ) objval
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nonzero variables = '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) isol(1:numsol)

  objval = 1652.0D+00
  numsol = 3
  isol(1:3) =  (/ 7, 4, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CORRECT RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Objective function = '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,g14.6)' ) objval
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nonzero variables = '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) isol(1:numsol)

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests MULTI_KNAP.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 7

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ), dimension ( m ) :: b = (/ 19.0D+00, 14.0D+00, 17.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: c = (/ &
    31.0D+00, 26.0D+00, 27.0D+00, 29.0D+00, 32.0D+00, 30.0D+00, 28.0D+00 /)
  integer ( kind = 4 ) isol(n)
  integer ( kind = 4 ) numsol
  real ( kind = 8 ) objval

  a(1,1) =   4.0D+00
  a(1,2) =   5.0D+00
  a(1,3) =   3.0D+00
  a(1,4) =   3.0D+00
  a(1,5) =   7.0D+00
  a(1,6) =   8.0D+00
  a(1,7) =   8.0D+00

  a(2,1) =   3.0D+00
  a(2,2) =   7.0D+00
  a(2,3) =   4.0D+00
  a(2,4) =   9.0D+00
  a(2,5) =   8.0D+00
  a(2,6) =   5.0D+00
  a(2,7) =   6.0D+00

  a(3,1) =   3.0D+00
  a(3,2) =   1.0D+00
  a(3,3) =   2.0D+00
  a(3,4) =   5.0D+00
  a(3,5) =   4.0D+00
  a(3,6) =   4.0D+00
  a(3,7) =   6.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  MULTI_KNAP is a heuristic algoritm for the'
  write ( *, '(a)' ) '  multidimensional 0/1 knapsack problem.'

  call multi_knap ( m, n, a, b, c, objval, numsol, isol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COMPUTED ANSWERS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The objective function value:'
  write ( *, '(a)' ) ' '
  write ( *, '(g14.6)' ) objval
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The nonzero variables are:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) isol(1:numsol)

  objval = 88.0D+00
  numsol = 3
  isol(1:3) = (/ 1, 3, 6 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CORRECT ANSWERS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The objective function value:'
  write ( *, '(a)' ) ' '
  write ( *, '(g14.6)' ) objval
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The nonzero variables are:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) isol(1:numsol)

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests PARTITION.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) cost(2*n,2*n)
  integer ( kind = 4 ) i
  logical init
  integer ( kind = 4 ) ip(n)
  integer ( kind = 4 ) iq(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kp(n)
  integer ( kind = 4 ) kq(n)
  real ( kind = 8 ) tcost

  data ( ( cost(i,j), j = 1, 2*n ), i = 1, 2*n ) / &
    0.0D+00, 2.0D+00, 4.0D+00, 7.0D+00, 4.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00, 5.0D+00, 1.0D+00, &
    2.0D+00, 0.0D+00, 3.0D+00, 6.0D+00, 3.0D+00, &
    1.0D+00, 1.0D+00, 0.0D+00, 1.0D+00, 5.0D+00, &
    4.0D+00, 3.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, &
    1.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, &
    7.0D+00, 6.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 1.0D+00, &
    4.0D+00, 3.0D+00, 2.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 2.0D+00, 0.0D+00, 4.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 3.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00, &
    5.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 5.0D+00, 0.0D+00, 1.0D+00, 4.0D+00, &
    1.0D+00, 3.0D+00, 1.0D+00, 1.0D+00, 0.0D+00 /

  init = .true.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  PARTITION is a heuristic algorithm'
  write ( *, '(a)' ) '  for the graph partition problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The order of the sets, and of the items in the set,'
  write ( *, '(a)' ) '  does not matter.'

  call partition ( n, cost, init, ip, iq, kp, kq, tcost )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COMPUTED RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First set:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) kp(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Second set:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) kq(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Total cost = ', tcost

  tcost = 25.0D+00

  kp(1:5) = (/ 9,  3, 1, 4, 2 /)
  kq(1:5) = (/ 5, 10, 7, 6, 8 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CORRECT RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First set:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) kp(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Second set:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) kq(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Total cost = ', tcost

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests PMATCH.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ), dimension ((n*(n-1))/2) :: cost = (/ &
    3.0D+00, 3.0D+00, 5.0D+00, 2.0D+00, 2.0D+00, &
             4.0D+00, 7.0D+00, 5.0D+00, 3.0D+00, &
                      4.0D+00, 7.0D+00, 5.0D+00, &
                               5.0D+00, 8.0D+00, &
                                        3.0D+00  /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pair(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  PMATCH is a heuristic algorithm'
  write ( *, '(a)' ) '  for the graph matching problem.'

  call pmatch ( n, cost, pair )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COMPUTED RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node   Paired Node'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i8)' ) i, pair(i)
  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests STEINER.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 20
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ns = 5

  real ( kind = 8 ), dimension ( m ) :: cost = (/ &
    6.0D+00, 2.0D+00, 6.0D+00, 9.0D+00, 3.0D+00, &
    4.0D+00, 5.0D+00, 9.0D+00, 3.0D+00, 4.0D+00, &
    3.0D+00, 2.0D+00, 3.0D+00, 5.0D+00, 8.0D+00, &
    6.0D+00, 7.0D+00, 4.0D+00, 4.0D+00, 8.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( m ) :: inode = (/ &
    2,  3,  8, 10,  4,  3,  7,  8,  9, 10, &
    4,  7,  9,  5, 10,  5,  7,  5,  7,  9 /)
  integer ( kind = 4 ) istree(n)
  integer ( kind = 4 ), dimension ( m ) :: jnode = (/ &
    1,  2,  3,  6,  2,  4,  3,  7,  8,  9, &
    1,  4,  7,  1,  5,  4,  5,  6,  6,  6 /)
  integer ( kind = 4 ) jstree(n)
  integer ( kind = 4 ) nsp
  logical spoint(n)
  real ( kind = 8 ) xlen

  spoint(1:n) = .false.

  spoint(3) = .true.
  spoint(1) = .true.
  spoint(10) = .true.
  spoint(6) = .true.
  spoint(9) = .true.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  STEINER is a heuristic algorithm'
  write ( *, '(a)' ) '  for the Steiner tree problem.'

  call steiner ( n, m, inode, jnode, cost, ns, spoint, nsp, &
    istree, jstree, xlen )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COMPUTED RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Steiner tree edges:'
  write ( *, '(a)' ) ' '
  do i = 1, nsp
    write ( *, '(4x,2i3)' ) istree(i), jstree(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Total length = ', xlen

  nsp = 6
  istree(1:6) = (/ 3, 4, 4, 7, 7,  9 /)
  jstree(1:6) = (/ 4, 7, 1, 9, 6, 10 /)
  xlen = 20.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CORRECT RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Steiner tree edges:'
  write ( *, '(a)' ) ' '
  do i = 1, nsp
    write ( *, '(4x,2i3)' ) istree(i), jstree(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Total length = ', xlen

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests TSP.
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 15

  real ( kind = 8 ) dist(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) isol(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) length

  data ( ( dist(i,j), j = 1, n ), i = 1, n ) / &
    0.0D+00, 29.0D+00, 82.0D+00, 46.0D+00, 68.0D+00, &
   52.0D+00, 72.0D+00, 42.0D+00, 51.0D+00, 55.0D+00, &
   29.0D+00, 74.0D+00, 23.0D+00, 72.0D+00, 46.0D+00, &
   29.0D+00,  0.0D+00, 55.0D+00, 46.0D+00, 42.0D+00, &
   43.0D+00, 43.0D+00, 23.0D+00, 23.0D+00, 31.0D+00, &
   41.0D+00, 51.0D+00, 11.0D+00, 52.0D+00, 21.0D+00, &
   82.0D+00, 55.0D+00,  0.0D+00, 68.0D+00, 46.0D+00, &
   55.0D+00, 23.0D+00, 43.0D+00, 41.0D+00, 29.0D+00, &
   79.0D+00, 21.0D+00, 64.0D+00, 31.0D+00, 51.0D+00, &
   46.0D+00, 46.0D+00, 68.0D+00,  0.0D+00, 82.0D+00, &
   15.0D+00, 72.0D+00, 31.0D+00, 62.0D+00, 42.0D+00, &
   21.0D+00, 51.0D+00, 51.0D+00, 43.0D+00, 64.0D+00, &
   68.0D+00, 42.0D+00, 46.0D+00, 82.0D+00,  0.0D+00, &
   74.0D+00, 23.0D+00, 52.0D+00, 21.0D+00, 46.0D+00, &
   82.0D+00, 58.0D+00, 46.0D+00, 65.0D+00, 23.0D+00, &
   52.0D+00, 43.0D+00, 55.0D+00, 15.0D+00, 74.0D+00, &
    0.0D+00, 61.0D+00, 23.0D+00, 55.0D+00, 31.0D+00, &
   33.0D+00, 37.0D+00, 51.0D+00, 29.0D+00, 59.0D+00, &
   72.0D+00, 43.0D+00, 23.0D+00, 72.0D+00, 23.0D+00, &
   61.0D+00,  0.0D+00, 42.0D+00, 23.0D+00, 31.0D+00, &
   77.0D+00, 37.0D+00, 51.0D+00, 46.0D+00, 33.0D+00, &
   42.0D+00, 23.0D+00, 43.0D+00, 31.0D+00, 52.0D+00, &
   23.0D+00, 42.0D+00,  0.0D+00, 33.0D+00, 15.0D+00, &
   37.0D+00, 33.0D+00, 33.0D+00, 31.0D+00, 37.0D+00, &
   51.0D+00, 23.0D+00, 41.0D+00, 62.0D+00, 21.0D+00, &
   55.0D+00, 23.0D+00, 33.0D+00,  0.0D+00, 29.0D+00, &
   62.0D+00, 46.0D+00, 29.0D+00, 51.0D+00, 11.0D+00, &
   55.0D+00, 31.0D+00, 29.0D+00, 42.0D+00, 46.0D+00, &
   31.0D+00, 31.0D+00, 15.0D+00, 29.0D+00,  0.0D+00, &
   51.0D+00, 21.0D+00, 41.0D+00, 23.0D+00, 37.0D+00, &
   29.0D+00, 41.0D+00, 79.0D+00, 21.0D+00, 82.0D+00, &
   33.0D+00, 77.0D+00, 37.0D+00, 62.0D+00, 51.0D+00, &
    0.0D+00, 65.0D+00, 42.0D+00, 59.0D+00, 61.0D+00, &
   74.0D+00, 51.0D+00, 21.0D+00, 51.0D+00, 58.0D+00, &
   37.0D+00, 37.0D+00, 33.0D+00, 46.0D+00, 21.0D+00, &
   65.0D+00,  0.0D+00, 61.0D+00, 11.0D+00, 55.0D+00, &
   23.0D+00, 11.0D+00, 64.0D+00, 51.0D+00, 46.0D+00, &
   51.0D+00, 51.0D+00, 33.0D+00, 29.0D+00, 41.0D+00, &
   42.0D+00, 61.0D+00,  0.0D+00, 62.0D+00, 23.0D+00, &
   72.0D+00, 52.0D+00, 31.0D+00, 43.0D+00, 65.0D+00, &
   29.0D+00, 46.0D+00, 31.0D+00, 51.0D+00, 23.0D+00, &
   59.0D+00, 11.0D+00, 62.0D+00,  0.0D+00, 59.0D+00, &
   46.0D+00, 21.0D+00, 51.0D+00, 64.0D+00, 23.0D+00, &
   59.0D+00, 33.0D+00, 37.0D+00, 11.0D+00, 37.0D+00, &
   61.0D+00, 55.0D+00, 23.0D+00, 59.0D+00,  0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  TSP is a heuristic algorithm'
  write ( *, '(a)' ) '  for the traveling salesman problem.'

  call tsp ( n, dist, isol )

  length = 0.0D+00
  do i = 1, n
    i1 = isol(i)
    if ( i < n ) then
      i2 = isol(i+1)
    else
      i2 = isol(1)
    end if
    length = length + dist(i1,i2)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  COMPUTED RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) isol(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Path length = ', length

  isol(1:15) =  (/ &
    1, 13,  2, 15,  9, &
    5,  7,  3, 12, 14, &
   10,  8,  6,  4, 11 /)

  length = 0.0D+00
  do i = 1, n
    i1 = isol(i)
    if ( i < n ) then
      i2 = isol(i+1)
    else
      i2 = isol(1)
    end if
    length = length + dist(i1,i2)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CORRECT RESULTS:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) isol(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Path length = ', length

  return
end
