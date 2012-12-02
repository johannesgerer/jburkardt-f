program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIANGULATION_PRB.
!
!  Discussion:
!
!    TRIANGULATION_PRB tests routines from the TRIANGULATION library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TRIANGULATION library.'

  call test01 ( )
  call test02 ( )
  call test025 ( )
  call test026 ( )
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
  call test125 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
  call test213 ( )
  call test215 ( )
  call test217 ( )
  call test219 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
  call test25 ( )
  call test26 ( )
  call test265 ( )
  call test27 ( )

  call test31 ( )
  call test32 ( )
  call test33 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ALPHA_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha_ave
  real ( kind = 8 ) alpha_area
  real ( kind = 8 ) alpha_min
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), parameter :: element_order = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  ALPHA_MEASURE returns the ALPHA measure of'
  write ( *, '(a)' ) '  quality of a triangulation.'
!
!  Get the sizes.
!
  call triangulation_order3_example1_size ( node_num, element_num, hole_num )
!
!  Allocate space.
!
  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
!
!  Get the triangulation data.
!
  call triangulation_order3_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )
!
!  Compute the triangulation quality.
!
  call alpha_measure ( node_num, node_xy, element_order, element_num, &
    element_node, alpha_min, alpha_ave, alpha_area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  ALPHA_MIN  = ', alpha_min
  write ( *, '(a,f12.6)' ) '  ALPHA_AVE  = ', alpha_ave
  write ( *, '(a,f12.6)' ) '  ALPHA_AREA = ', alpha_area
!
!  Free the memory.
!
  deallocate ( node_xy )
  deallocate ( element_node )
  deallocate ( element_neighbor )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests AREA_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area_ave
  real ( kind = 8 ) area_max
  real ( kind = 8 ) area_min
  real ( kind = 8 ) area_ratio
  real ( kind = 8 ) area_std
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  real ( kind = 8 ) quality
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), parameter :: element_order = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  AREA_MEASURE returns the AREA measure of'
  write ( *, '(a)' ) '  quality of a triangulation.'
!
!  Get the sizes.
!
  call triangulation_order3_example1_size ( node_num, element_num, hole_num )
!
!  Allocate space.
!
  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
!
!  Get the triangulation data.
!
  call triangulation_order3_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )
!
!  Compute the triangulation quality.
!
  call area_measure ( node_num, node_xy, element_order, element_num, &
    element_node, area_min, area_max, area_ratio, area_ave, area_std )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  AREA_MIN   = ', area_min
  write ( *, '(a,f12.6)' ) '  AREA_MAX   = ', area_max
  write ( *, '(a,f12.6)' ) '  AREA_RATIO = ', area_ratio
  write ( *, '(a,f12.6)' ) '  AREA_AVE   = ', area_ave
  write ( *, '(a,f12.6)' ) '  AREA_STD   = ', area_std
!
!  Free the memory.
!
  deallocate ( node_xy )
  deallocate ( element_node )
  deallocate ( element_neighbor )

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests DELAUNAY_SWAP_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 4
  integer ( kind = 4 ), parameter :: element_num = 2
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ) alpha_area
  real ( kind = 8 ) alpha_ave
  real ( kind = 8 ) alpha_min_swapped
  real ( kind = 8 ) alpha_min_unswapped
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) :: seed = 123456789
  logical              swap
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ) element_node(element_order,element_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  DELAUNAY_SWAP_TEST determines whether two triangles'
  write ( *, '(a)' ) '  with a common edge need to "swap" diagonals.'
  write ( *, '(a)' ) '  If swapping is indicated, then ALPHA_MIN should increase.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Swap   ALPHA_MIN   ALPHA_MIN'
  write ( *, '(a)' ) '         Unswapped   Swapped'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
!
!  Generate a random quadrilateral (1,2,3,4).
!
    call quad_convex_random ( seed, node_xy )
!
!  Does it need swapping?
!
    call delaunay_swap_test ( node_xy, swap )
!
!  Compute ALPHA_MIN unswapped.
!
    element_node(1:3,1) = (/ 1, 2, 3 /)
    element_node(1:3,2) = (/ 1, 3, 4 /)

    call alpha_measure ( node_num, node_xy, element_order, element_num, &
      element_node, alpha_min_unswapped, alpha_ave, alpha_area )
!
!  Compute ALPHA_MIN swapped.
!
    element_node(1:3,1) = (/ 1, 2, 4 /)
    element_node(1:3,2) = (/ 2, 3, 4 /)

    call alpha_measure ( node_num, node_xy, element_order, element_num, &
      element_node, alpha_min_swapped, alpha_ave, alpha_area )

    if ( .false. ) then
      call r8mat_transpose_print ( 2, node_num, node_xy, '  Quadrilateral' )
    end if

    write ( *, '(2x,3x,l1,2x,f10.6,2x,f10.6)' ) &
      swap, alpha_min_unswapped, alpha_min_swapped

  end do

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests DIAEDG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 4
  integer ( kind = 4 ), parameter :: element_num = 2
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ) alpha_area
  real ( kind = 8 ) alpha_ave
  real ( kind = 8 ) alpha_min_swapped
  real ( kind = 8 ) alpha_min_unswapped
  integer ( kind = 4 ) diaedg
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) :: seed = 123456789
  logical              swap
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  DIAEDG determines whether two triangles'
  write ( *, '(a)' ) '  with a common edge need to "swap" diagonals.'
  write ( *, '(a)' ) '  If swapping is indicated, then ALPHA_MIN should increase.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Swap   ALPHA_MIN   ALPHA_MIN'
  write ( *, '(a)' ) '         Unswapped   Swapped'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
!
!  Generate a random quadrilateral (1,2,3,4).
!
    call quad_convex_random ( seed, node_xy )
!
!  Does it need swapping?
!
    value = diaedg ( &
      node_xy(1,1), node_xy(2,1), &
      node_xy(1,2), node_xy(2,2), &
      node_xy(1,3), node_xy(2,3), &
      node_xy(1,4), node_xy(2,4) )

    if ( value == 1 ) then
      swap = .false.
    else
      swap = .true.
    end if
!
!  Compute ALPHA_MIN unswapped.
!
    element_node(1:3,1) = (/ 1, 2, 3 /)
    element_node(1:3,2) = (/ 1, 3, 4 /)

    call alpha_measure ( node_num, node_xy, element_order, element_num, &
      element_node, alpha_min_unswapped, alpha_ave, alpha_area )
!
!  Compute ALPHA_MIN swapped.
!
    element_node(1:3,1) = (/ 1, 2, 4 /)
    element_node(1:3,2) = (/ 2, 3, 4 /)

    call alpha_measure ( node_num, node_xy, element_order, element_num, &
      element_node, alpha_min_swapped, alpha_ave, alpha_area )

    if ( .false. ) then
      call r8mat_transpose_print ( 2, node_num, node_xy, '  Quadrilateral' )
    end if

    write ( *, '(2x,3x,l1,2x,f10.6,2x,f10.6)' ) &
      swap, alpha_min_unswapped, alpha_min_swapped

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests NODE_MERGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num = 15
  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_rep(node_num)
  real ( kind = 8 ), dimension (2,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       1.0D+00, 0.0D+00, &
       3.0D+00, 0.0D+00, &
       4.0D+00, 0.0D+00, &
       1.0D+00, 1.0D+00, &
       4.0D+00, 1.0D+00, &
       2.0D+00, 2.0D+00, &
       3.0D+00, 3.0D+00, &
       2.0D+00, 3.5D+00, &
       0.5D+00, 4.0D+00, &
       1.0D+00, 4.0D+00, &
       1.5D+00, 4.0D+00, &
       4.0D+00, 4.0D+00, &
       1.0D+00, 4.5D+00, &
       1.0D+00, 4.5D+00 /), (/ 2, node_num /) )
  integer ( kind = 4 ) rep
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) test
  real ( kind = 8 ) tolerance
  real ( kind = 8 ), dimension ( test_num ) :: tolerance_test = (/ &
    0.01D+00, 0.75D+00, 1.2D+00, 1.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  NODE_MERGE identifies groups of nodes'
  write ( *, '(a)' ) '  that can be merged, with a given tolerance.'

  call r8mat_transpose_print ( 2, node_num, node_xy, '  Node coordinates:' )

  do test = 1, test_num

    tolerance = tolerance_test(test)

    call node_merge ( dim_num, node_num, node_xy, tolerance, node_rep )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  TOLERANCE = ', tolerance
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      Node  Representatives:'
    write ( *, '(a)' ) ' '

    do node = 1, node_num
      write ( *, '(2x,i8,2x,i8)' ) node, node_rep(node)
    end do
!
!  Make a list of the node representatives.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      Rep   Coordinates:'
    write ( *, '(a)' ) ' '

    call i4vec_sort_heap_a ( node_num, node_rep )

    rep_num = 0

    do node = 1, node_num

      if ( 2 <= node ) then
        if ( node_rep(node-1) == node_rep(node) ) then
          cycle
        end if
      end if

      rep = node_rep(node)
      rep_num = rep_num + 1

      write ( *, '(2x,i8,2x,2g14.6)' ) rep_num, node_xy(1:2,rep)

    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests NS_ADJ_COL_SET, NS_ADJ_COUNT and NS_ADJ_ROW_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 15
  integer ( kind = 4 ), parameter :: element_num = 4
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ), parameter :: variable_num = 36

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), dimension ( variable_num+1 ) :: adj_col
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_row
  character ( len = 80 ) :: file_name = 'ns_triangulation.eps'
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_show
  integer ( kind = 4 ), dimension ( node_num ) :: node_p_variable = (/ &
    3, -1,  8, -1, 13, &
   -1, -1, -1, -1, &
   24, -1, 29, &
   -1, -1, &
   36 /)
  integer ( kind = 4 ), dimension ( node_num ) :: node_u_variable = (/ &
    1,  4,  6,  9, 11, &
   14, 16, 18, 20, &
   22, 25, 27, &
   30, 32, &
   34 /)
  integer ( kind = 4 ), dimension ( node_num ) :: node_v_variable = (/ &
    2,  5,  7, 10, 12, &
   15, 17, 19, 21, &
   23, 26, 28, &
   31, 33, &
   35 /)
  real ( kind = 8 ), dimension ( 2, node_num ) :: node_xy = &
    reshape ( (/ &
   0.0D+00, 0.0D+00, &
   0.0D+00, 1.0D+00, &
   0.0D+00, 2.0D+00, &
   0.0D+00, 3.0D+00, &
   0.0D+00, 4.0D+00, &
   1.0D+00, 0.0D+00, &
   1.0D+00, 1.0D+00, &
   1.0D+00, 2.0D+00, &
   1.0D+00, 3.0D+00, &
   2.0D+00, 0.0D+00, &
   2.0D+00, 1.0D+00, &
   2.0D+00, 2.0D+00, &
   3.0D+00, 0.0D+00, &
   3.0D+00, 1.0D+00, &
   4.0D+00, 0.0D+00 /), &
    (/ 2, node_num /) )
  integer ( kind = 4 ) num
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rhi
  integer ( kind = 4 ) rlo
  integer ( kind = 4 ), dimension ( 3, element_num ) :: element_neighbor = &
    reshape ( (/ &
    -1,  2, -1, &
     3,  1,  4, &
     2, -1, -1, &
    -1, -1,  2 /), (/ 3, element_num /) )
  integer ( kind = 4 ), dimension ( element_order, element_num ) :: element_node = &
    reshape ( (/ &
     1, 10,  3,  6,  7,  2, &
    12,  3, 10,  8,  7, 11, &
     3, 12,  5,  8,  9,  4, &
    10, 15, 12, 13, 14, 11 /), (/ element_order, element_num /) )
  integer ( kind = 4 ) element_show
  integer ( kind = 4 ) variable

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For an order 3/order 6 Taylor Hood triangulation'
  write ( *, '(a)' ) '  for Navier Stokes velocity and pressure,'
  write ( *, '(a)' ) '  NS_ADJ_COUNT counts variable adjacencies'
  write ( *, '(a)' ) '  and sets up the sparse compressed column'
  write ( *, '(a)' ) '  column pointer array.'
  write ( *, '(a)' ) '  NS_ADJ_COL_SET sets up the sparse compressed column'
  write ( *, '(a)' ) '  COL vector.'
  write ( *, '(a)' ) '  NS_ADJ_ROW_SET sets up the sparse compressed column'
  write ( *, '(a)' ) '  ROW vector.'
!
!  Plot the example.
!
  node_show = 2
  element_show = 2

  call triangulation_order6_plot ( file_name, node_num, node_xy, &
    element_num, element_node, node_show, element_show )
!
!  Get the count of the variable adjacencies.
!  We don't really need to make this call, since the next
!  call does the calculation as part of getting ADJ_COL.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of variables is ', variable_num

  call ns_adj_count ( node_num, element_num, variable_num, element_node, &
    element_neighbor, node_u_variable, node_v_variable, node_p_variable, &
    adj_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)'    ) '  As computed by NS_ADJ_COUNT,'
  write ( *, '(a,i8)' ) '  Number of variable adjacency entries is ', adj_num
!
!  Get the count of the variable adjacencies and the COL vector.
!
  call ns_adj_col_set ( node_num, element_num, variable_num, element_node, &
    element_neighbor, node_u_variable, node_v_variable, node_p_variable, &
    adj_num, adj_col )

  write ( *, '(a)' ) ' '
  write ( *, '(a)'    ) '  As computed by NS_ADJ_COL_SET,'
  write ( *, '(a,i8)' ) '  Number of variable adjacency entries is ', adj_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Variable adjacency column pointers:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Variable     First      Last    Number'
  write ( *, '(a)' ) ' '

  do variable = 1, variable_num

    num = adj_col(variable+1) - adj_col(variable)

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) &
      variable, adj_col(variable), adj_col(variable+1)-1, num

  end do
!
!  Get the ROW vector.
!
  allocate ( adj_row(1:adj_num) )

  call ns_adj_row_set ( node_num, element_num, variable_num, element_node, &
    element_neighbor, node_u_variable, node_v_variable, node_p_variable, &
    adj_num, adj_col, adj_row )
!
!  This is a huge array.  We only print out the beginning and end.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Variable adjacency row entries:'
  write ( *, '(a)' ) '  (Partial printout only)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Entry     Row       Col'
  write ( *, '(a)' ) ' '

  do variable = 1, variable_num

    rlo = adj_col(variable)
    rhi = adj_col(variable+1)-1

    if ( variable <= 3 .or. variable_num - 3 <= variable ) then

      write ( *, '(a)' ) ' '

      do r = rlo, rhi
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
          r, adj_row(r), variable
      end do

    end if

    if ( variable == 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  (SKIPPING MANY MANY ENTRIES...)'
      write ( *, '(a)' ) ' '
    end if

  end do

  deallocate ( adj_row )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests POINTS_DELAUNAY_NAIVE_2D.
!
!  Diagram:
!
!    !....3&11....
!    !............
!    !............
!    X..9.........
!    !.....5......
!    !...........6
!    !.4.2...10...
!    !.....8...12.
!    V............
!    !..7.........
!    !......1.....
!    !............
!    !............
!    !----V----X--
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtri = 20
  integer ( kind = 4 ), parameter :: node_num = 12
  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) tri(element_order,maxtri)
  integer ( kind = 4 ) element_num
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xy = reshape ( (/ &
     7.0D+00,  3.0D+00, &
     4.0D+00,  7.0D+00, &
     5.0D+00, 13.0D+00, &
     2.0D+00,  7.0D+00, &
     6.0D+00,  9.0D+00, &
    12.0D+00,  8.0D+00, &
     3.0D+00,  4.0D+00, &
     6.0D+00,  6.0D+00, &
     3.0D+00, 10.0D+00, &
     8.0D+00,  7.0D+00, &
     5.0D+00, 13.0D+00, &
    10.0D+00,  6.0D+00 /), (/ dim_num, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  POINTS_DELAUNAY_NAIVE_2D computes the Delaunay'
  write ( *, '(a)' ) '    triangulation of a set of nodes.'

  call r8mat_transpose_print ( dim_num, node_num, node_xy, '  The nodes:' )

  call points_delaunay_naive_2d ( node_num, node_xy, maxtri, element_num, tri )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of triangles is ', element_num

  call i4mat_transpose_print ( element_order, element_num, tri, &
    '  The triangles:' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests POINTS_HULL_2D.
!
!  Diagram:
!
!    !....3.......
!    !............
!    !..9.........
!    !.....5......
!    !...........6
!    !.4.2...10...
!    !.....8......
!    !.........12.
!    !..7.........
!    !......1.....
!    !............
!    !............
!    !-----------
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num = 12

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival(node_num)
  integer ( kind = 4 ) nval
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xy = reshape ( (/ &
       7.0D+00,  3.0D+00, &
       4.0D+00,  7.0D+00, &
       5.0D+00, 13.0D+00, &
       2.0D+00,  7.0D+00, &
       6.0D+00,  9.0D+00, &
      12.0D+00,  8.0D+00, &
       3.0D+00,  4.0D+00, &
       6.0D+00,  6.0D+00, &
       3.0D+00, 10.0D+00, &
       8.0D+00,  7.0D+00, &
       5.0D+00, 13.0D+00, &
      10.0D+00,  6.0D+00 /), (/ dim_num, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  POINTS_HULL_2D computes the convex hull'
  write ( *, '(a)' ) '    of a set of nodes.'

  call r8mat_transpose_print ( dim_num, node_num, node_xy, '  The nodes:' )

  call points_hull_2d ( node_num, node_xy, nval, ival )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The convex hull is formed by connecting:'
  write ( *, '(a)' ) ' '
  do i = 1, nval
    write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) &
      i, ival(i), node_xy(1:dim_num,ival(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The correct sequence of nodes is:'
  write ( *, '(a)' ) '  4, 9, 3, 6, 12, 1, 7, (4).'
  write ( *, '(a)' ) ' '

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests Q_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  real ( kind = 8 ) q_area
  real ( kind = 8 ) q_ave
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  real ( kind = 8 ) quality
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), parameter :: element_order = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Q_MEASURE returns the Q measure of'
  write ( *, '(a)' ) '  quality of a triangulation.'
!
!  Get the sizes.
!
  call triangulation_order3_example1_size ( node_num, element_num, hole_num )
!
!  Allocate space.
!
  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
!
!  Get the triangulation data.
!
  call triangulation_order3_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )
!
!  Compute the triangulation quality.
!
  call q_measure ( node_num, node_xy, element_order, element_num, &
    element_node, q_min, q_max, q_ave, q_area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Q_MIN  = ', q_min
  write ( *, '(a,f12.6)' ) '  Q_MAX  = ', q_max
  write ( *, '(a,f12.6)' ) '  Q_AVE  = ', q_ave
  write ( *, '(a,f12.6)' ) '  Q_AREA = ', q_area
!
!  Free the memory.
!
  deallocate ( node_xy )
  deallocate ( element_node )
  deallocate ( element_neighbor )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests R8TRIS2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num = 9
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ), dimension (dim_num,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       0.2D+00, 0.5D+00, &
       0.3D+00, 0.6D+00, &
       0.4D+00, 0.5D+00, &
       0.6D+00, 0.4D+00, &
       0.6D+00, 0.5D+00, &
       1.0D+00, 0.0D+00, &
       1.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) triangle(element_order,2*node_num)
  integer ( kind = 4 ) element_neighbor(3,2*node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  R8TRIS2 computes the Delaunay triangulation of'
  write ( *, '(a)' ) '    a set of nodes in 2D.'
!
!  Set up the Delaunay triangulation.
!
  call r8tris2 ( node_num, node_xy, element_num, triangle, element_neighbor )

  call triangulation_order3_print ( node_num, element_num, node_xy, &
    triangle, element_neighbor )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) ref2(2,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(2,3) :: t = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, &
    2.0D+00, 5.0D+00 /), (/ 2, 3 /) )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For an order 3 triangle,'
  write ( *, '(a)' ) '  TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE '
  write ( *, '(a)' ) '    maps a physical point to a reference point.'
  write ( *, '(a)' ) '  TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL '
  write ( *, '(a)' ) '    maps a reference point to a physical point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )'
  write ( *, '(a)' ) ' '

  call triangle_reference_sample ( n, seed, ref )

  call element_order3_reference_to_physical ( t, n, ref, phy )

  call element_order3_physical_to_reference ( t, n, phy, ref2 )

  do j = 1, n

    write ( *, '(2x,2f8.4,2x,2f8.4,2x,2f8.4)' ) &
      ref(1:2,j), phy(1:2,j), ref2(1:2,j)

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) ref2(2,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(2,6) :: t = reshape ( (/ &
    7.0D+00, 2.0D+00, &
    9.0D+00, 2.0D+00, &
    7.0D+00, 3.0D+00, &
    8.0D+00, 2.0D+00, &
    8.0D+00, 2.5D+00, &
    7.0D+00, 2.5D+00 /), (/ 2, 6 /) )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For an order 6 triangle,'
  write ( *, '(a)' ) '  TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE '
  write ( *, '(a)' ) '    maps a physical point to a reference point.'
  write ( *, '(a)' ) '  TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL '
  write ( *, '(a)' ) '    maps a reference point to a physical point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )'
  write ( *, '(a)' ) ' '

  call triangle_reference_sample ( n, seed, ref )

  call element_order6_reference_to_physical ( t, n, ref, phy )

  call element_order6_physical_to_reference ( t, n, phy, ref2 )

  do j = 1, n

    write ( *, '(2x,2f8.4,2x,2f8.4,2x,2f8.4)' ) &
      ref(1:2,j), phy(1:2,j), ref2(1:2,j)

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests TRIANGULATION_NODE_ORDER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 36
  integer ( kind = 4 ), parameter :: element_num = 41
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ), dimension ( node_num ) :: node_order
  integer ( kind = 4 ), dimension ( element_order, element_num ) :: element_node &
    = reshape ( (/ &
     1,  8,  7, &
     1,  2,  8, &
     2,  9,  8, &
     2,  3,  9, &
     3, 10,  9, &
     3,  4, 10, &
     4, 11, 10, &
     4,  5, 11, &
     5, 12, 11, &
     5,  6, 12, &
     7, 14, 13, &
     7,  8, 14, &
     8, 15, 14, &
     8,  9, 15, &
    11, 18, 17, &
    11, 12, 18, &
    13, 20, 19, &
    13, 14, 20, &
    14, 21, 20, &
    14, 15, 21, &
    15, 22, 21, &
    15, 16, 22, &
    16, 23, 22, &
    16, 17, 23, &
    17, 24, 23, &
    17, 18, 24, &
    19, 26, 25, &
    19, 20, 26, &
    21, 28, 27, &
    21, 22, 28, &
    25, 30, 29, &
    25, 26, 30, &
    26, 31, 30, &
    27, 32, 31, &
    27, 28, 32, &
    29, 34, 33, &
    29, 30, 34, &
    30, 35, 34, &
    30, 31, 35, &
    31, 36, 35, &
    31, 32, 36 /), (/ element_order, element_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  TRIANGULATION_NODE_ORDER computes the order'
  write ( *, '(a)' ) '  of the nodes in a triangulation.'

  call triangulation_node_order ( element_order, element_num, &
    element_node, node_num, node_order )

  call i4vec_print ( node_num, node_order, '  NODE ORDER:' )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests TRIANGULATION_ORDER3_ADJ_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_col
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) k
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ) element_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_ADJ_COUNT counts node adjacencies'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_ADJ_SET sets node adjacencies.'
!
!  Get the sizes of the example.
!
  call triangulation_order3_example1_size ( node_num, element_num, hole_num )

  allocate ( adj_col(1:node_num+1) )
  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
  allocate ( element_node(1:element_order,1:element_num) )
!
!  Get the data of the example.
!
  call triangulation_order3_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )
!
!  Get the count of the node adjacencies.
!
  call triangulation_order3_adj_count ( node_num, element_num, element_node, &
    element_neighbor, adj_num, adj_col )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of adjacency entries is ', adj_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Adjacency pointers:'
  write ( *, '(a)' ) ' '
  do node = 1, node_num
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) node, adj_col(node), adj_col(node+1)-1
  end do
!
!  Get the node adjacencies.
!
  allocate ( adj(1:adj_num) )

  call triangulation_order3_adj_set ( node_num, element_num, element_node, &
    element_neighbor, adj_num, adj_col, adj )
!
!  Print the node adjacencies.
!
  do node = 1, node_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Nodes adjacent to node ', node
    write ( *, '(a)' ) ' '

    do k = adj_col(node), adj_col(node+1)-1
      write ( *, '(2x,i8)' ) adj(k)
    end do

  end do

  deallocate ( adj )
  deallocate ( adj_col )
  deallocate ( node_xy )
  deallocate ( element_neighbor )
  deallocate ( element_node )

  return
end
subroutine test125 ( )

!*****************************************************************************80
!
!! TEST125 tests TRIANGULATION_ORDER3_ADJ_SET2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_col
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ia
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ja
  integer ( kind = 4 ) k
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ) element_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST125'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_ADJ_COUNT counts node adjacencies'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_ADJ_SET2 sets node adjacencies'
  write ( *, '(a)' ) '  as a pair of vectors IA(*), JA(*).'
!
!  Get the sizes of the example.
!
  call triangulation_order3_example1_size ( node_num, element_num, hole_num )

  allocate ( adj_col(1:node_num+1) )
  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
  allocate ( element_node(1:element_order,1:element_num) )
!
!  Get the data of the example.
!
  call triangulation_order3_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )
!
!  Get the count of the node adjacencies.
!
  call triangulation_order3_adj_count ( node_num, element_num, element_node, &
    element_neighbor, adj_num, adj_col )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of adjacency entries is ', adj_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Adjacency pointers:'
  write ( *, '(a)' ) ' '
  do node = 1, node_num
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) node, adj_col(node), adj_col(node+1)-1
  end do
!
!  Get the node adjacencies.
!
  allocate ( ia(1:adj_num) )
  allocate ( ja(1:adj_num) )

  call triangulation_order3_adj_set2 ( node_num, element_num, element_node, &
    element_neighbor, adj_num, adj_col, ia, ja )
!
!  Print the node adjacencies.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node adjacencies stored in IA(*) and JA(*):'
  write ( *, '(a)' ) ' '

  do adj = 1, adj_num

    write ( *, '(2x,i8,2x,a,i2,a,i2,a)' ) adj, '(', ia(adj), ',', ja(adj), ')'

  end do

  deallocate ( adj_col )
  deallocate ( ia )
  deallocate ( ja )
  deallocate ( node_xy )
  deallocate ( element_neighbor )
  deallocate ( element_node )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num = 36
  integer ( kind = 4 ), parameter :: element_num = 41
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) boundary_edge_num
  character ( len = 80 ) :: file_name = 'triangulation_order3_plot2.eps'
  integer ( kind = 4 ) :: node_show = 2
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xy = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    3.0D+00, 0.0D+00, &
    4.0D+00, 0.0D+00, &
    5.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, &
    2.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, &
    4.0D+00, 1.0D+00, &
    5.0D+00, 1.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, &
    2.0D+00, 2.0D+00, &
    3.0D+00, 2.0D+00, &
    4.0D+00, 2.0D+00, &
    5.0D+00, 2.0D+00, &
    0.0D+00, 3.0D+00, &
    1.0D+00, 3.0D+00, &
    2.0D+00, 3.0D+00, &
    3.0D+00, 3.0D+00, &
    4.0D+00, 3.0D+00, &
    5.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00, &
    1.0D+00, 4.0D+00, &
    2.0D+00, 4.0D+00, &
    3.0D+00, 4.0D+00, &
    0.0D+00, 5.0D+00, &
    1.0D+00, 5.0D+00, &
    2.0D+00, 5.0D+00, &
    3.0D+00, 5.0D+00, &
    0.0D+00, 6.0D+00, &
    1.0D+00, 6.0D+00, &
    2.0D+00, 6.0D+00, &
    3.0D+00, 6.0D+00 /), (/ dim_num, node_num /) )
  integer ( kind = 4 ), dimension ( element_order, element_num ) :: &
    element_node = reshape ( (/ &
     1,  8,  7, &
     1,  2,  8, &
     2,  9,  8, &
     2,  3,  9, &
     3, 10,  9, &
     3,  4, 10, &
     4, 11, 10, &
     4,  5, 11, &
     5, 12, 11, &
     5,  6, 12, &
     7, 14, 13, &
     7,  8, 14, &
     8, 15, 14, &
     8,  9, 15, &
    11, 18, 17, &
    11, 12, 18, &
    13, 20, 19, &
    13, 14, 20, &
    14, 21, 20, &
    14, 15, 21, &
    15, 22, 21, &
    15, 16, 22, &
    16, 23, 22, &
    16, 17, 23, &
    17, 24, 23, &
    17, 18, 24, &
    19, 26, 25, &
    19, 20, 26, &
    21, 28, 27, &
    21, 22, 28, &
    25, 30, 29, &
    25, 26, 30, &
    26, 31, 30, &
    27, 32, 31, &
    27, 28, 32, &
    29, 34, 33, &
    29, 30, 34, &
    30, 35, 34, &
    30, 31, 35, &
    31, 36, 35, &
    31, 32, 36 /), (/ element_order, element_num /) )
  integer ( kind = 4 ) :: element_show = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the'
  write ( *, '(a)' ) '    boundary edges.'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PLOT plots the triangulation.'

  call triangulation_order3_plot ( file_name, node_num, node_xy, &
    element_num, element_node, node_show, element_show )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An Encapsulated PostScript image of this'
  write ( *, '(a)' ) '  triangulation is in "' // trim ( file_name ) // '".'

  call triangulation_order3_boundary_edge_count ( element_num, element_node, &
    boundary_edge_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary edges = ', boundary_edge_num
  write ( *, '(a,i8)' ) '  Correct number =           ', 33

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ), parameter :: hole_num = 2
  integer ( kind = 4 ), parameter :: node_num = 36
  integer ( kind = 4 ), parameter :: element_num = 41

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER'
  write ( *, '(a)' ) '  determines the number of edges that lie on the'
  write ( *, '(a)' ) '  boundary of a region that has been triangulated.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =          ', node_num
  write ( *, '(a,i8)' ) '  Number of triangles =      ', element_num
  write ( *, '(a,i8)' ) '  Number of holes =          ', hole_num

  call triangulation_order3_boundary_edge_count_euler ( node_num, &
    element_num, hole_num, boundary_num )

  write ( *, '(a,i8)' ) '  Number of boundary edges = ', boundary_num

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests TRIANGULATION_ORDER3_BOUNDARY_NODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num = 36
  integer ( kind = 4 ), parameter :: element_num = 41
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) i
  logical node_boundary(node_num)
  integer ( kind = 4 ), dimension ( element_order, element_num ) :: &
    element_node = reshape ( (/ &
     1,  8,  7, &
     1,  2,  8, &
     2,  9,  8, &
     2,  3,  9, &
     3, 10,  9, &
     3,  4, 10, &
     4, 11, 10, &
     4,  5, 11, &
     5, 12, 11, &
     5,  6, 12, &
     7, 14, 13, &
     7,  8, 14, &
     8, 15, 14, &
     8,  9, 15, &
    11, 18, 17, &
    11, 12, 18, &
    13, 20, 19, &
    13, 14, 20, &
    14, 21, 20, &
    14, 15, 21, &
    15, 22, 21, &
    15, 16, 22, &
    16, 23, 22, &
    16, 17, 23, &
    17, 24, 23, &
    17, 18, 24, &
    19, 26, 25, &
    19, 20, 26, &
    21, 28, 27, &
    21, 22, 28, &
    25, 30, 29, &
    25, 26, 30, &
    26, 31, 30, &
    27, 32, 31, &
    27, 28, 32, &
    29, 34, 33, &
    29, 30, 34, &
    30, 35, 34, &
    30, 31, 35, &
    31, 36, 35, &
    31, 32, 36 /), (/ element_order, element_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_BOUNDARY_NODE determines which'
  write ( *, '(a)' ) '  nodes lie on the boundary.'

  call triangulation_order3_boundary_node ( node_num, element_num, &
    element_node, node_boundary )

  call lvec_print ( node_num, node_boundary, '    Node  BN?' )

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests TRIANGULATION_ORDER3_CHECK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: element_num = 16
  integer ( kind = 4 ), parameter :: node_num = 13
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isave
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) element_num2
  integer ( kind = 4 ), dimension (element_order,element_num ) :: &
    element_node = reshape ( (/ &
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
    10,   6,  11 /), (/ element_order, element_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_CHECK checks the triangulation.'

  call i4mat_transpose_print ( element_order, element_num, element_node, &
    '  Triangles:' );
!
!  Pass all tests.
!
  call triangulation_order3_check ( node_num, element_num, element_node, &
    ierror )

  write ( *, '(a,i8)' ) '  Error code = ', ierror
!
!  Fail test 1.
!
  node_num2 = 2

  call triangulation_order3_check ( node_num2, element_num, &
    element_node, ierror )

  write ( *, '(a,i8)' ) '  Error code = ', ierror
!
!  Fail test 2.
!
  element_num2 = 0

  call triangulation_order3_check ( node_num, element_num2, &
    element_node, ierror )

  write ( *, '(a,i8)' ) '  Error code = ', ierror
!
!  Fail test 3.
!
  isave = element_node(2,5)
  element_node(2,5) = 0

  call triangulation_order3_check ( node_num, element_num, element_node, &
    ierror )

  write ( *, '(a,i8)' ) '  Error code = ', ierror
  element_node(2,5) = isave
!
!  Fail test 4.
!
  isave = element_node(3,10)
  element_node(3,10) = 2 * node_num + 1

  call triangulation_order3_check ( node_num, element_num, element_node, &
    ierror )

  write ( *, '(a,i8)' ) '  Error code = ', ierror
  element_node(3,10) = isave
!
!  Fail test 5.
!
  element_node(3,4) = 3
  element_node(3,8) = 3
  element_node(3,10) = 3
  element_node(3,11) = 3
  element_node(2,14) = 3

  call triangulation_order3_check ( node_num, element_num, element_node, &
    ierror )
  write ( *, '(a,i8)' ) '  Error code = ', ierror

  element_node(3,4) = 5
  element_node(3,8) = 5
  element_node(3,10) = 5
  element_node(3,11) = 5
  element_node(2,14) = 5
!
!  Fail test 6.
!
  element_node(1,9) = 7
  call triangulation_order3_check ( node_num, element_num, element_node, &
    ierror )
  write ( *, '(a,i8)' ) '  Error code = ', ierror
  element_node(1,9) = 9
!
!  Fail test 7.
!
  element_node(3,7) = 2
  call triangulation_order3_check ( node_num, element_num, element_node, &
    ierror )
  write ( *, '(a,i8)' ) '  Error code = ', ierror
  element_node(3,7) = 9

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests TRIANGULATION_ORDER3_EXAMPLE1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_EXAMPLE1_SIZE gives the sizes'
  write ( *, '(a)' ) '    for an example triangulation;'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_EXAMPLE1 returns the information'
  write ( *, '(a)' ) '    for an example triangulation;'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PRINT prints a triangulation.'
!
!  Get the sizes.
!
  call triangulation_order3_example1_size ( node_num, element_num, hole_num )

  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
  allocate ( element_node(1:element_order,1:element_num) )

  call triangulation_order3_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )

  call triangulation_order3_print ( node_num, element_num, node_xy, &
    element_node, element_neighbor )

  deallocate ( node_xy )
  deallocate ( element_neighbor )
  deallocate ( element_node )

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests TRIANGULATION_ORDER3_NEIGHBOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: element_num = 16
  integer ( kind = 4 ), parameter :: node_num = 13
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ), dimension (element_order,element_num ) :: &
    element_node = reshape ( (/ &
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
    10,   6,  11 /), (/ element_order, element_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_NEIGHBOR determines the'
  write ( *, '(a)' ) '  triangle neighbors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T1  S1  T2  S2'
  write ( *, '(a)' ) ' '
  do t1 = 1, element_num
    do s1 = 1, 3
      call triangulation_order3_neighbor ( element_num, element_node, &
        t1, s1, t2, s2 )
      write ( *, '(2x,4i4)' ) t1, s1, t2, s2
    end do
  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests TRIANGULATION_NEIGHBOR_ELEMENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: element_num = 16
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension (element_order,element_num) :: &
    element_node = reshape ( (/ &
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
    10,   6,  11 /), (/ element_order, element_num /) )
  integer ( kind = 4 ), dimension (3,element_num) :: element_neighbor

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  TRIANGULATION_NEIGHBOR_ELEMENTS determines the'
  write ( *, '(a)' ) '  adjacency relationships between elements.'

  call i4mat_transpose_print ( element_order, element_num, element_node, &
    '  Elements:' )

  call triangulation_neighbor_elements ( element_order, element_num, &
    element_node, element_neighbor )

  call i4mat_transpose_print ( 3, element_num, element_neighbor, &
    '  Element neighbors:' )

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests TRIANGULATION_ORDER3_PLOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: file_name = 'triangulation_order3_plot.eps'
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) :: node_show = 2
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ) :: element_show = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PLOT can plot a triangulation.'
!
!  Get the sizes.
!
  call triangulation_order3_example1_size ( node_num, element_num, hole_num )
!
!  Allocate space.
!
  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
!
!  Get the example data.
!
  call triangulation_order3_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )
!
!  Make the plot.
!
  call triangulation_order3_plot ( file_name, node_num, node_xy, element_num, &
    element_node, node_show, element_show )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PLOT has created an'
  write ( *, '(a)' ) '  Encapsulated PostScript file (EPS) containing'
  write ( *, '(a)' ) '  an image of the triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This file is called "' // trim ( file_name ) //'".'

  deallocate ( node_xy )
  deallocate ( element_node )
  deallocate ( element_neighbor )

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests TRIANGULATION_ORDER3_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 9
  integer ( kind = 4 ), parameter :: element_num = 12
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ), dimension (2,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       0.2D+00, 0.5D+00, &
       0.3D+00, 0.6D+00, &
       0.4D+00, 0.5D+00, &
       0.6D+00, 0.4D+00, &
       0.6D+00, 0.5D+00, &
       1.0D+00, 0.0D+00, &
       1.0D+00, 1.0D+00 /), (/ 2, node_num /) )
  integer ( kind = 4 ), dimension (element_order,element_num) :: &
    element_node = reshape ( (/ &
       2, 1, 3, &
       3, 1, 6, &
       2, 3, 4, &
       4, 3, 5, &
       7, 4, 5, &
       5, 3, 6, &
       7, 5, 6, &
       9, 4, 7, &
       6, 1, 8, &
       7, 6, 8, &
       7, 8, 9, &
       2, 4, 9 /), (/ element_order, element_num /) )
  integer ( kind = 4 ), dimension (3,element_num) :: element_neighbor = reshape ( (/ &
       -28,   2,  3, &
         1,   9,  6, &
         1,   4, 12, &
         3,   6,  5, &
         8,   4,  7, &
         4,   2,  7, &
         5,   6, 10, &
        12,   5, 11, &
         2, -34, 10, &
         7,   9, 11, &
        10, -38,  8, &
         3,   8, -3 /), (/ 3, element_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PRINT prints out the data.'

  call triangulation_order3_print ( node_num, element_num, node_xy, &
    element_node, element_neighbor )

  return
end
subroutine test213 ( )

!*****************************************************************************80
!
!! TEST213 tests TRIANGULATION_ORDER3_QUAD.
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

  integer ( kind = 4 ), parameter :: quad_num = 6

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n11
  integer ( kind = 4 ) n12
  integer ( kind = 4 ) n21
  integer ( kind = 4 ) n22
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  external quad_fun
  real ( kind = 8 ) quad_value
  real ( kind = 8 ), dimension ( quad_num ) :: quad_w = (/ &
    0.1666666666666666D+00, &
    0.1666666666666666D+00, &
    0.1666666666666666D+00, &
    0.1666666666666666D+00, &
    0.1666666666666666D+00, &
    0.1666666666666666D+00 /)
  real ( kind = 8 ), dimension(2,quad_num) :: quad_xy = reshape ( (/ &
    0.659027622374092D+00,  0.231933368553031D+00, &
    0.659027622374092D+00,  0.109039009072877D+00, &
    0.231933368553031D+00,  0.659027622374092D+00, &
    0.231933368553031D+00,  0.109039009072877D+00, &
    0.109039009072877D+00,  0.659027622374092D+00, &
    0.109039009072877D+00,  0.231933368553031D+00  &
  /), (/ 2, quad_num /) )
  real ( kind = 8 ) region_area
  integer ( kind = 4 ) test
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), parameter :: test_num = 4
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ) element_num
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST213'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_QUAD can apply a quadrature rule'
  write ( *, '(a)' ) '  to every triangle in a triangulated region,'
  write ( *, '(a)' ) '  and estimate the integral of a function over'
  write ( *, '(a)' ) '  that region.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NODE_NUM   TRI_NUM  Integral estim  Area of Region'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
!
!  Set up the grid.
!
    n = 2**( test - 1 )
    node_num = ( n + 1 ) * ( n + 1 )

    allocate ( node_xy(1:2,1:node_num) )

    k = 0
    do j = 1, n + 1
      y = real ( j - 1, kind = 8 ) / real ( n + 1 - 1, kind = 8 )
      do i = 1, n + 1
        x = real ( i - 1, kind = 8 ) / real ( n + 1 - 1, kind = 8 )
        k = k + 1
        node_xy(1:2,k) = (/ x, y /)
      end do
    end do
!
!  Set up the triangulation.
!
    element_num = 2 * n * n

    allocate ( element_node(1:element_order,1:element_num) )

    k = 0
    do j = 1, n
      do i = 1, n

        n11 = i     + ( j     - 1 ) * ( n + 1 )
        n12 = i     + ( j + 1 - 1 ) * ( n + 1 )
        n21 = i + 1 + ( j     - 1 ) * ( n + 1 )
        n22 = i + 1 + ( j + 1 - 1 ) * ( n + 1 )

        k = k + 1
        element_node(1:element_order,k) = (/ n11, n21, n12 /)
        k = k + 1
        element_node(1:element_order,k) = (/ n22, n12, n21 /)

      end do
    end do
!
!  Estimate the integral.
!
    call triangulation_order3_quad ( node_num, node_xy, element_order, &
      element_num, element_node, quad_fun, quad_num, quad_xy, quad_w, &
      quad_value, region_area )

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) &
      node_num, element_num, quad_value, region_area
!
!  Delete allocatables.
!
    deallocate ( node_xy )
    deallocate ( element_node )

  end do

  return
end
subroutine quad_fun ( n, xy_vec, f_vec )

!*****************************************************************************80
!
!! QUAD_FUN is a sample integrand function for TRIANGULATION_QUAD.
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
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XY_VEC(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F_VEC(N), the value of the integrand
!    function at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_vec(n)
  real ( kind = 8 ) xy_vec(2,n)

  f_vec(1:n) = exp ( xy_vec(1,1:n)**2 + xy_vec(2,1:n)**2 )

  return
end
subroutine test215 ( )

!*****************************************************************************80
!
!! TEST215 tests TRIANGULATION_ORDER3_REFINE_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num1 = 5
  integer ( kind = 4 ), parameter :: element_num1 = 3
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_data
  integer ( kind = 4 ) node_num2
  real ( kind = 8 ), dimension (dim_num,node_num1) :: node_xy1 = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       1.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       1.0D+00, 1.0D+00, &
       0.5D+00, 1.5D+00  /), (/ dim_num, node_num1 /) )
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy2
  integer ( kind = 4 ), dimension (element_order,element_num1) :: &
    element_node1 = reshape ( (/ &
       1, 2, 3, &
       4, 3, 2, &
       3, 4, 5 /), (/ element_order, element_num1 /) )
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node2
  integer ( kind = 4 ) element_num2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST215'
  write ( *, '(a)' ) '  For an order3 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_REFINE_SIZE determines the'
  write ( *, '(a)' ) '  size of a refined triangulation.'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_REFINE_COMPUTES computes the'
  write ( *, '(a)' ) '  refined triangulation.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes is ', node_num1
  write ( *, '(a,i8)' ) '  The number of triangles is ', element_num1

  call r8mat_transpose_print ( dim_num, node_num1, node_xy1, &
    '  The nodes' )

  call i4mat_transpose_print ( element_order, element_num1, element_node1, &
    '  The triangles:' )

  allocate ( edge_data(5,3*element_num1) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sizing the refined mesh:'

  call triangulation_order3_refine_size ( node_num1, element_num1, &
    element_node1, node_num2, element_num2, edge_data )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Information about the refined mesh:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes is ', node_num2
  write ( *, '(a,i8)' ) '  The number of triangles is ', element_num2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing the refined mesh:'

  allocate ( node_xy2(dim_num,node_num2) )
  allocate ( element_node2(element_order,element_num2) )

  call triangulation_order3_refine_compute ( node_num1, element_num1, &
    node_xy1, element_node1, node_num2, element_num2, edge_data, node_xy2, &
    element_node2 )

  call r8mat_transpose_print ( dim_num, node_num2, node_xy2, &
    '  The refined nodes' )

  call i4mat_transpose_print ( element_order, element_num2, element_node2, &
    '  The refined triangles:' )

  deallocate ( edge_data )
  deallocate ( node_xy2 )
  deallocate ( element_node2 )

  return
end
subroutine test217 ( )

!*****************************************************************************80
!
!! TEST217 tests TRIANGULATION_SEARCH_DELAUNAY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num = 13
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) dist
  real ( kind = 8 ) dnear
  integer ( kind = 4 ) edge
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) nnear
  real ( kind = 8 ), dimension (dim_num,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       2.0D+00, 2.0D+00, &
      -1.0D+00, 3.0D+00, &
      -2.0D+00, 2.0D+00, &
       8.0D+00, 2.0D+00, &
       9.0D+00, 5.0D+00, &
       7.0D+00, 4.0D+00, &
       5.0D+00, 6.0D+00, &
       6.0D+00, 7.0D+00, &
       8.0D+00, 8.0D+00, &
      11.0D+00, 7.0D+00, &
      10.0D+00, 4.0D+00, &
       6.0D+00, 4.0D+00 /), (/ dim_num, node_num /) )
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) td(test_num)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_node(element_order,2*node_num)
  integer ( kind = 4 ) element_index
  integer ( kind = 4 ) element_neighbor(3,2*node_num)
  real ( kind = 8 ) xd(dim_num,test_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST217'
  write ( *, '(a)' ) '  Given a set of nodes NODE_XY, and a single point XD,'
  write ( *, '(a)' ) '  find the nearest node in NODE_XY to XD.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINTS_POINT_NEAR_NAIVE_ND uses a naive method.'
  write ( *, '(a)' ) '  TRIANGULATION_SEARCH_DELAUNAY finds a triangle'
  write ( *, '(a)' ) '    containing the point.  Often, one of these vertices'
  write ( *, '(a)' ) '    is the closest point.'
!
!  Set up the Delaunay triangulation.
!
  call r8tris2 ( node_num, node_xy, element_num, element_node, &
    element_neighbor )
!
!  Get the test points.
!
  seed = 123456789

  call triangulation_order3_sample ( node_num, node_xy, element_num, &
    element_node, test_num, seed, xd, td )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X        Y     Distance  Index     Steps'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = xd(1:dim_num,test)

    call points_point_near_naive_nd ( dim_num, node_num, node_xy, &
      p, nnear, dnear )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2f8.4   )' ) '  XD       ', p(1:dim_num)
    write ( *, '(a,3f8.4,i8)' ) '  Naive    ', node_xy(1:dim_num,nnear), &
      dnear, nnear

    call triangulation_search_delaunay ( node_num, node_xy, element_order, &
      element_num, element_node, element_neighbor, p, element_index, alpha, &
      beta, gamma, edge, step_num )

    if ( element_index < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Error: the search failed.'
      cycle
    end if

    i1 = element_node(1,element_index)
    d1 = sqrt ( sum ( ( p(1:dim_num) - node_xy(1:dim_num,i1) )**2 ) )

    dist = d1
    nnear = i1

    i2 = element_node(2,element_index)
    d2 = sqrt ( sum ( ( p(1:dim_num) - node_xy(1:dim_num,i2) )**2 ) )

    if ( d2 < dist ) then
      dnear = d2
      nnear = i2
    end if

    i3 = element_node(3,element_index)
    d3 = sqrt ( sum ( ( p(1:dim_num) - node_xy(1:dim_num,i3) )**2 ) )

    if ( d3 < dist ) then
      dnear = d3
      nnear = i3
    end if

    write ( *, '(a,3f8.4,i8,2x,i8)' ) &
      '  Delaunay ', node_xy(1:2,nnear), dnear, nnear, step_num

  end do

  return
end
subroutine test219 ( )

!*****************************************************************************80
!
!! TEST219 tests TRIANGULATION_SEARCH_DELAUNAY, TRIANGULATION_SEARCH_NAIVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num = 13
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) edge
  real ( kind = 8 ) gamma
  real ( kind = 8 ), dimension (dim_num,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       2.0D+00, 2.0D+00, &
      -1.0D+00, 3.0D+00, &
      -2.0D+00, 2.0D+00, &
       8.0D+00, 2.0D+00, &
       9.0D+00, 5.0D+00, &
       7.0D+00, 4.0D+00, &
       5.0D+00, 6.0D+00, &
       6.0D+00, 7.0D+00, &
       8.0D+00, 8.0D+00, &
      11.0D+00, 7.0D+00, &
      10.0D+00, 4.0D+00, &
       6.0D+00, 4.0D+00 /), (/ dim_num, node_num /) )
  real ( kind = 8 ) p_test(dim_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) t_test(test_num)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_node(element_order,2*node_num)
  integer ( kind = 4 ) element_index1
  integer ( kind = 4 ) element_index2
  integer ( kind = 4 ) element_neighbor(3,2*node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST219'
  write ( *, '(a)' ) '  Given a triangulation, and a point P,'
  write ( *, '(a)' ) '  find the triangle T containing to P.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TRIANGULATION_SEARCH_NAIVE uses a naive method.'
  write ( *, '(a)' ) '  TRIANGULATION_SEARCH_DELAUNAY uses a method that will'
  write ( *, '(a)' ) '    work fast if the triangulation is Delaunay.'
!
!  Set up the Delaunay triangulation.
!
  call r8tris2 ( node_num, node_xy, element_num, element_node, &
    element_neighbor )
!
!  Get the test points.
!
  call triangulation_order3_sample ( node_num, node_xy, element_num, &
    element_node, test_num, seed, p_test, t_test )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X         Y     Naive  Delaunay    ' // &
    '   Alpha        Beta       Gamma     Steps'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call triangulation_search_naive ( node_num, node_xy, element_order, &
      element_num, element_node, p_test(1:dim_num,test), element_index1 )

    call triangulation_search_delaunay ( node_num, node_xy, element_order, &
      element_num, element_node, element_neighbor, p_test(1:dim_num,test), &
      element_index2, alpha, beta, gamma, edge, step_num )

    write ( *, &
      '(2x,f8.4,2x,f8.4,2x,i8,2x,i8,2x,f10.4,2x,f10.4,2x,f10.4,2x,i8)' ) &
      p_test(1:dim_num,test), element_index1, element_index2, &
      alpha, beta, gamma, step_num

  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests TRIANGULATION_ORDER6_ADJ_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_col
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) k
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) element_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  For an order6 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_ADJ_COUNT counts node adjacencies'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_ADJ_SET sets node adjacencies.'
!
!  Get the sizes of the example.
!
  call triangulation_order6_example1_size ( node_num, element_num, hole_num )

  allocate ( adj_col(1:node_num+1) )
  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
  allocate ( element_node(1:element_order,1:element_num) )
!
!  Get the data of the example.
!
  call triangulation_order6_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )
!
!  Get the count of the node adjacencies.
!
  call triangulation_order6_adj_count ( node_num, element_num, element_node, &
    element_neighbor, adj_num, adj_col )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of adjacency entries is ', adj_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Adjacency pointers:'
  write ( *, '(a)' ) ' '
  do node = 1, node_num
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) node, adj_col(node), adj_col(node+1)-1
  end do
!
!  Get the node adjacencies.
!
  allocate ( adj(1:adj_num) )

  call triangulation_order6_adj_set ( node_num, element_num, element_node, &
    element_neighbor, adj_num, adj_col, adj )
!
!  Print the node adjacencies.
!
  do node = 1, node_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Nodes adjacent to node ', node
    write ( *, '(a)' ) ' '

    do k = adj_col(node), adj_col(node+1)-1
      write ( *, '(2x,i8)' ) adj(k)
    end do

  end do

  deallocate ( adj )
  deallocate ( adj_col )
  deallocate ( node_xy )
  deallocate ( element_neighbor )
  deallocate ( element_node )

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) boundary_edge_num
  integer ( kind = 4 ) :: dim_num = 2
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  For an order6 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the'
  write ( *, '(a)' ) '    boundary edges.'

  call triangulation_order6_example1_size ( node_num, element_num, hole_num )

  allocate ( node_xy(1:dim_num,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( element_neighbor(1:3,1:element_num) )

  call triangulation_order6_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )

  call triangulation_order6_boundary_edge_count ( element_num, &
    element_node, boundary_edge_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary edges = ', boundary_edge_num
  write ( *, '(a,i8)' ) '  Correct number =           ', 16

  deallocate ( node_xy )
  deallocate ( element_node )
  deallocate ( element_neighbor )

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  For an order6 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER'
  write ( *, '(a)' ) '  determines the number of edges that lie on the'
  write ( *, '(a)' ) '  boundary of a region that has been triangulated.'

  call triangulation_order6_example1_size ( node_num, element_num, hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =          ', node_num
  write ( *, '(a,i8)' ) '  Number of triangles =      ', element_num
  write ( *, '(a,i8)' ) '  Number of holes =          ', hole_num

  call triangulation_order6_boundary_edge_count_euler ( node_num, &
    element_num, hole_num, boundary_num )

  write ( *, '(a,i8)' ) '  Number of boundary edges = ', boundary_num
  write ( *, '(a,i8)' ) '  Correct number =           ', 16

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests TRIANGULATION_ORDER6_BOUNDARY_NODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: file_name = 'triangulation_order6_plot.eps'
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) hole_num
  logical, allocatable, dimension ( : ) :: node_boundary
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) :: node_show = 2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) :: element_show = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  For an order6 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_BOUNDARY_COUNT counts the boundary'
  write ( *, '(a)' ) '    edges.'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_PLOT plots the triangulation.'

  call triangulation_order6_example1_size ( node_num, element_num, hole_num )

  allocate ( node_boundary(1:node_num) )
  allocate ( node_xy(1:dim_num,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( element_neighbor(1:3,1:element_num) )

  call triangulation_order6_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )
!
!  Make the plot.
!
  call triangulation_order6_plot ( file_name, node_num, node_xy, element_num, &
    element_node, node_show, element_show )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An Encapsulated PostScript image of this'
  write ( *, '(a)' ) '  triangulation is in "' // trim ( file_name ) // '".'

  call triangulation_order6_boundary_node ( node_num, element_num, &
    element_node, node_boundary )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Node  BN?'
  write ( *, '(a)' ) ' '

  do i = 1, node_num
    write ( *, '(2x,i8,2x,l1)' ) i, node_boundary(i)
  end do

  deallocate ( node_boundary )
  deallocate ( node_xy )
  deallocate ( element_node )
  deallocate ( element_neighbor )

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests TRIANGULATION_ORDER6_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  For an order6 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_PRINT prints the data.'

  call triangulation_order6_example1_size ( node_num, element_num, hole_num )

  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
  allocate ( element_node(1:element_order,1:element_num) )

  call triangulation_order6_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )

  call triangulation_order6_print ( node_num, element_num, node_xy, &
    element_node, element_neighbor )

  deallocate ( node_xy )
  deallocate ( element_neighbor )
  deallocate ( element_node )

  return
end
subroutine test265 ( )

!*****************************************************************************80
!
!! TEST265 tests TRIANGULATION_ORDER6_REFINE_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num1 = 12
  integer ( kind = 4 ), parameter :: element_num1 = 3
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_data
  integer ( kind = 4 ) node_num2
  real ( kind = 8 ), dimension (dim_num,node_num1) :: node_xy1 = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       2.0D+00, 0.0D+00, &
       0.0D+00, 2.0D+00, &
       2.0D+00, 2.0D+00, &
       1.0D+00, 3.0D+00, &
       1.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       1.0D+00, 1.0D+00, &
       2.0D+00, 1.0D+00, &
       1.0D+00, 2.0D+00, &
       0.5D+00, 2.5D+00, &
       1.5D+00, 2.5D+00  /), (/ dim_num, node_num1 /) )
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy2
  integer ( kind = 4 ), dimension (element_order,element_num1) :: &
    element_node1 = reshape ( (/ &
       1,  2,  3,  6,  8,  7, &
       4,  3,  2,  9, 10,  8, &
       3,  4,  5, 10, 12, 11 /), (/ element_order, element_num1 /) )
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node2
  integer ( kind = 4 ) element_num2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST265'
  write ( *, '(a)' ) '  For an order6 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_REFINE_SIZE determines the'
  write ( *, '(a)' ) '  size of a refined triangulation.'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_REFINE_COMPUTES computes the'
  write ( *, '(a)' ) '  refined triangulation.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes is ', node_num1
  write ( *, '(a,i8)' ) '  The number of triangles is ', element_num1

  call r8mat_transpose_print ( dim_num, node_num1, node_xy1, &
    '  The nodes' )

  call i4mat_transpose_print ( element_order, element_num1, element_node1, &
    '  The triangles:' )

  allocate ( edge_data(5,3*element_num1) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sizing the refined mesh:'

  call triangulation_order6_refine_size ( node_num1, element_num1, &
    element_node1, node_num2, element_num2, edge_data )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Information about the refined mesh:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes is ', node_num2
  write ( *, '(a,i8)' ) '  The number of triangles is ', element_num2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing the refined mesh:'

  allocate ( node_xy2(dim_num,node_num2) )
  allocate ( element_node2(element_order,element_num2) )

  call triangulation_order6_refine_compute ( node_num1, element_num1, &
    node_xy1, element_node1, node_num2, element_num2, edge_data, node_xy2, &
    element_node2 )

  call r8mat_transpose_print ( dim_num, node_num2, node_xy2, &
    '  The refined nodes' )

  call i4mat_transpose_print ( element_order, element_num2, element_node2, &
    '  The refined triangles:' )

  deallocate ( edge_data )
  deallocate ( node_xy2 )
  deallocate ( element_node2 )

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests TRIANGULATION_ORDER6_VERTEX_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) midside_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) vertex_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  For an order6 triangulation:'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER6_VERTEX_COUNT counts the '
  write ( *, '(a)' ) '  vertex nodes and midside nodes.'

  call triangulation_order6_example1_size ( node_num, element_num, hole_num )

  allocate ( node_xy(1:2,1:node_num) )
  allocate ( element_neighbor(1:3,1:element_num) )
  allocate ( element_node(1:element_order,1:element_num) )

  call triangulation_order6_example1 ( node_num, element_num, node_xy, &
    element_node, element_neighbor )

  call triangulation_order6_vertex_count ( node_num, element_num, &
    element_node, vertex_num, midside_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =         ', node_num
  write ( *, '(a,i8)' ) '  Number of vertex nodes =  ', vertex_num
  write ( *, '(a,i8)' ) '  Number of midside nodes = ', midside_num

  deallocate ( node_xy )
  deallocate ( element_neighbor )
  deallocate ( element_node )

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests VORONOI_POLYGON_AREA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: neighbor_num = 4
  integer ( kind = 4 ), parameter :: node_num = 5

  real ( kind = 8 ) area
  real ( kind = 8 ), parameter :: area_correct = 0.5D+00
  integer ( kind = 4 ), parameter :: center = 5
  integer ( kind = 4 ), dimension ( neighbor_num ) :: neighbor_index = (/ 1, 2, 3, 4 /)
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xy = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, &
    0.5D+00, 0.5D+00 /), (/ dim_num, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  VORONOI_POLYGON_AREA computes the area of'
  write ( *, '(a)' ) '  a finite Voronoi polygon.'

  call voronoi_polygon_area ( center, neighbor_num, neighbor_index, &
    node_num, node_xy, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The computed area is ', area
  write ( *, '(a,g14.6)' ) '  The correct area is  ', area_correct

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests VORONOI_POLYGON_CENTROID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: neighbor_num = 4
  integer ( kind = 4 ), parameter :: node_num = 5

  integer ( kind = 4 ), parameter :: center = 5
  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ), dimension (dim_num) :: centroid_exact = (/ &
    0.5D+00, 0.5D+00 /)
  integer ( kind = 4 ), dimension ( neighbor_num ) :: neighbor_index = (/ 1, 2, 3, 4 /)
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xy = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, &
    0.5D+00, 0.5D+00 /), (/ dim_num, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  VORONOI_POLYGON_CENTROID computes the centroid of'
  write ( *, '(a)' ) '  a finite Voronoi polygon.'

  call voronoi_polygon_centroid ( center, neighbor_num, neighbor_index, &
    node_num, node_xy, centroid )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Computed centroid ', centroid(1:dim_num)
  write ( *, '(a,2g14.6)' ) '  Correct centroid  ', centroid_exact(1:dim_num)

  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests VORONOI_POLYGON_VERTICES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: neighbor_num = 4
  integer ( kind = 4 ), parameter :: node_num = 5

  integer ( kind = 4 ), parameter :: center = 5
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( neighbor_num ) :: neighbor_index = (/ 1, 2, 3, 4 /)
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xy = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, &
    0.5D+00, 0.5D+00 /), (/ dim_num, node_num /) )
  real ( kind = 8 ) v(2,neighbor_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  VORONOI_POLYGON_VERTICES computes the vertices of'
  write ( *, '(a)' ) '  a finite Voronoi polygon.'

  call voronoi_polygon_vertices ( center, neighbor_num, neighbor_index, &
    node_num, node_xy, v )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Voronoi Polygon Vertex coordinates:'
  write ( *, '(a)' ) ' '
  do i = 1, neighbor_num
    write ( *, '(2x,i8,4x,2g14.6)' ) i, v(1:2,i)
  end do

  return
end
