program main

!*****************************************************************************80
!
!! MAIN is the main program for CODEPACK_PRB.
!
!  Discussion:
!
!    CODEPACK_PRB calls the CODEPACK test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CODEPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CODEPACK library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test004 ( )
  call test005 ( )
  call test009 ( )
  call test010 ( )

  call test011 ( )
  call test012 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
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
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CODEPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests CDG_CODE_BACK, CDG_CODE_BRUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  For a color digraph code:'
  write ( *, '(a)' ) '  CDG_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) '  CDG_CODE_BRUTE uses brute force;'
!
!  Choose the example.
!
  call cdg_example_cube ( adj, nnode )

  call cdg_print ( adj, nnode, '  The color digraph adjacency matrix:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode
    order1(i) = i
  end do

  call cdg_order_code ( adj, nnode, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Initial node ordering:' )

  call cdg_code_print ( nnode, code1, '  The order-dependent code:' )
!
!  Brute force.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BRUTE FORCE calculation:'

  call cdg_code_brute ( adj, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Induced node ordering:' )

  call cdg_code_print ( nnode, code1, '  The maximal code:' )
!
!  Backtrack.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BACKTRACK calculation:'

  call cdg_code_back ( adj, nnode, code2, order2 )

  call node_order_print ( nnode, order2, '  Induced node ordering:' )

  call cdg_code_print ( nnode, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call cdg_code_compare ( code1, code2, nnode, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SUCCESS: The codes are equal.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FAILURE: The codes are unequal.'
  end if

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests CDG_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 65

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  character results(test_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  CDG_COMPARE compares color digraphs.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare all pairs of test graphs.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    example = i
    call cdg_example_octo ( example, adj1, nnode1, seed )

    do j = i, test_num

      example = j
      call cdg_example_octo ( example, adj2, nnode2, seed )

      call cdg_compare ( adj1, nnode1, adj2, nnode2, order1, &
        order2, result )

      if ( ( i == j .and. result /= 0 ) .or. &
           ( i /= j .and. result == 0 ) ) then

        write ( *, '(a,2i6)' ) '  FAILURE on graphs ', i, j
        call cdg_print ( adj1, nnode1, '  CDG #1:' )
        call cdg_print ( adj2, nnode2, '  CDG#2:' )

      end if

      if ( result < 0 ) then
        results(j,i) = '.'
        write ( results(i,j), '(i1)' ) abs ( result )
      else if ( 0 < result ) then
        results(i,j) = '.'
        write ( results(j,i), '(i1)' ) result
      else if ( result == 0 ) then
        results(i,j) = '='
        results(j,i) = '='
      end if

    end do

  end do

  do i = 1, test_num
    write ( *, '(2x,i2,2x,65a)' ) i, ( results(i,j), j = 1, test_num )
  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests CG_CODE_BACK, CG_CODE_BRUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  For a color graph code:'
  write ( *, '(a)' ) '  CG_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) '  CG_CODE_BRUTE uses brute force;'
!
!  Choose the example.
!
  call cg_example_cube ( adj, nnode )

  call cg_print ( adj, nnode, '  The color graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode
    order1(i) = i
  end do

  call cg_order_code ( adj, nnode, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Initial node ordering:' )

  call cg_code_print ( nnode, code1, '  The order-dependent code:' )
!
!  Brute force.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BRUTE FORCE calculation:'

  call cg_code_brute ( adj, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Induced node ordering:' )

  call cg_code_print ( nnode, code1, '  The maximal code:' )
!
!  Backtrack.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BACKTRACK calculation:'

  call cg_code_back ( adj, nnode, code2, order2 )

  call node_order_print ( nnode, order2, '  Induced node ordering:' )

  call cg_code_print ( nnode, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call cg_code_compare ( code1, code2, nnode, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SUCCESS: The codes are equal.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST003 - FAILURE'
    write ( *, '(a)' ) '  The codes are unequal.'
  end if

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests CG_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 40

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  character results(test_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  CG_COMPARE compares two color graphs.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare all pairs of test graphs.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    example = i
    call cg_example_octo ( example, adj1, nnode1, seed )

    do j = i, test_num

      example = j
      call cg_example_octo ( example, adj2, nnode2, seed )

      call cg_compare ( adj1, nnode1, adj2, nnode2, order1, &
        order2, result )

      if ( i == j .and. result /= 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  FAILURE on graph ', i
        call cg_print ( adj1, nnode1, '  Version #1 of the color graph:' )
        call cg_print ( adj2, nnode2, '  Version #2 of the color graph:' )

      else if ( i /= j .and. result == 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a,2i6)' ) '  FAILURE on graphs ', i, j

      end if

      if ( result < 0 ) then
        results(j,i) = '.'
        write ( results(i,j), '(i1)' ) abs ( result )
      else if ( 0 < result ) then
        results(i,j) = '.'
        write ( results(j,i), '(i1)' ) result
      else if ( result == 0 ) then
        results(i,j) = '='
        results(j,i) = '='
      end if

    end do

  end do

  do i = 1, test_num
    write ( *, '(2x,i2,2x,65a)' ) i, results(i,1:test_num)
  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests CG_CODE_COMPARE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  CG_CODE_COMPARE'
  write ( *, '(a)' ) '    compares two color graph codes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare the color graph codes of the cube and'
  write ( *, '(a)' ) '  the permuted cube.'
!
!  Set the graph to the color cube.
!
  call cg_example_cube ( adj1, nnode1 )

  call cg_print ( adj1, nnode1, '  The color cube:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode1
    order1(i) = i
  end do

  call cg_order_code ( adj1, nnode1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Initial node ordering:' )

  call cg_code_print ( nnode1, code1, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call cg_code_back ( adj1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Induced node ordering:' )

  call cg_code_print ( nnode1, code1, '  The maximal code:' )
!
!  Now permute the nodes of graph 1 to get graph 2, get its code and print it.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now permute the graph:'

  nnode2 = nnode1
!
!  Initialize the node ordering.
!
  do i = 1, nnode2
    order2(i) = i
  end do
!
!  Permute the node ordering.
!
  call i4vec_perm_random ( nnode2, order2, seed )
!
!  Update the adjacency matrix.
!
  do i = 1, nnode2
    ip = order2(i)
    do j = 1, nnode2
      jp = order2(j)
      adj2(i,j) = adj1(ip,jp)
    end do
  end do

  call cg_print ( adj2, nnode2, '  The color graph:' )
!
!  Compute the order dependent code.
!
  call cg_order_code ( adj2, nnode2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Initial node ordering:' )

  call cg_code_print ( nnode2, code2, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call cg_code_back ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call cg_code_print ( nnode2, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call cg_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  SUCCESS: CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)') '  FAILURE: CODE2 < CODE1'
  end if

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests CG_CODE_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) color_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) mcolor
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  CG_CODE_COMPARE'
  write ( *, '(a)' ) '    compares two color graph codes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare the color graph codes of the cube and'
  write ( *, '(a)' ) '  the cube with permuted colors.'
!
!  Set the graph to the color cube.
!
  call cg_example_cube ( adj1, nnode1 )
!
!  Count the colors.
!
  call cg_color_count ( adj1, nnode1, mcolor, ncolor )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of colors =    ', ncolor
  write ( *, '(a,i6)' ) '  Maximum color index = ', mcolor

  call cg_print ( adj1, nnode1, '  The color graph:' )
!
!  Compute the maximal code by backtracking.
!
  call cg_code_back ( adj1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Induced node ordering:' )

  call cg_code_print ( nnode1, code1, '  The maximal code:' )
!
!  Second example should have a higher graph code.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Graph 2 is made by permuting graph 1'
  write ( *, '(a)' ) '  and increasing the color of one node.'

  nnode2 = nnode1

  adj2(1:nnode2,1:nnode2) = adj1(1:nnode2,1:nnode2)
!
!  Initialize the node ordering.
!
  do i = 1, nnode2
    order2(i) = i
  end do
!
!  Get a random permutation for the colors.
!
  call i4vec_perm_random ( nnode2, order2, seed )
!
!  Permute the adjacency matrix.
!
  call i4mat_perm ( nnode2, adj2, order2 )
!
!  Alter one color.
!
  color_min = adj2(1,1)
  i_min = 1
  do i = 2, nnode2
    if ( adj2(i,i) < color_min ) then
      color_min = adj2(i,i)
      i_min = i
    end if
  end do

  adj2(i_min,i_min) = adj2(i_min,i_min) + 1
!
!  Print the matrix.
!
  call cg_print ( adj2, nnode2, '  The color graph:' )
!
!  Compute the maximal code by backtracking.
!
  call cg_code_back ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call cg_code_print ( nnode2, code2, &
    '  The maximal code, using backtracking:' )
!
!  Compare the codes.
!
  call cg_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  SUCCESS: CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  FAILURE: CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE2 < CODE1'
  end if
!
!  Compute the maximal code by brute force.
!
  call cg_code_brute ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call cg_code_print ( nnode2, code2, &
    '  The maximal code, by brute force:' )
!
!  Compare the codes.
!
  call cg_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  SUCCESS: CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  FAILURE: CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE2 < CODE1'
  end if

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests CG_CODE_COMPARE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj1(10,10)
  integer ( kind = 4 ) adj2(10,10)
  integer ( kind = 4 ) code1(10,10)
  integer ( kind = 4 ) code2(10,10)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) order1(10)
  integer ( kind = 4 ) order2(10)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  CG_CODE_COMPARE compares two color'
  write ( *, '(a)' ) '    graph codes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare color graph codes of the cube and a '
  write ( *, '(a)' ) '  graph with same number of nodes, links, and '
  write ( *, '(a)' ) '  colors.'
!
!  Set the graph to the color cube.
!
  call cg_example_cube ( adj1, nnode1 )

  call cg_print ( adj1, nnode1, '  The color graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode1
    order1(i) = i
  end do

  call cg_order_code ( adj1, nnode1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Initial node ordering:' )

  call cg_code_print ( nnode1, code1, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call cg_code_back ( adj1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Induced node ordering:' )

  call cg_code_print ( nnode1, code1, '  The maximal code:' )
!
!  Now choose a random graph on N1 nodes with the same number of links
!  and colors.
!
  nnode2 = nnode1

  call cg_edge_count ( adj1, nnode1, nedge )

  ncolor = 0
  do i = 1, nnode1
    ncolor = max ( ncolor, adj1(i,i) )
  end do

  call cg_random ( adj2, nnode2, ncolor, nedge, seed )

  call cg_print ( adj2, nnode2, '  The color graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode2
    order2(i) = i
  end do

  call cg_order_code ( adj2, nnode2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Initial node ordering:' )

  call cg_code_print ( nnode2, code2, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call cg_code_back ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call cg_code_print ( nnode2, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call cg_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)' ) '  CODE2 < CODE1'
  end if

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests DG_CODE_BACK;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  For a digraph code:'
  write ( *, '(a)' ) '  DG_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compare the digraph codes of'
  write ( *, '(a)' ) '  the cube digraph, and a node-reordered copy of'
  write ( *, '(a)' ) '  the cube digraph.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The codes should be the same.'
!
!  Set up the cube graph.
!
  call g_example_cube ( adj1, nnode1 )

  call g_print ( adj1, nnode1, '  The graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode1
    order1(i) = i
  end do

  call dg_order_code ( adj1, nnode1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Initial node ordering:' )

  call dg_code_print ( nnode1, code1, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call dg_code_back ( adj1,  nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Induced node ordering:' )

  call dg_code_print ( nnode1, code1, '  The maximal code:' )
!
!  Now permute the nodes of digraph 1 to get digraph 2,
!  get its code and print it.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now permute the original digraph:'

  nnode2 = nnode1
!
!  Get a random permutation.
!
  do i = 1, nnode2
    order2(i) = i
  end do

  call i4vec_perm_random ( nnode2, order2, seed )
!
!  Reorder the nodes of the digraph.
!
  do i = 1, nnode2
    ip = order2(i)
    do j = 1, nnode2
      jp = order2(j)
      adj2(i,j) = adj1(ip,jp)
    end do
  end do

  call g_print ( adj2, nnode2, '  The graph:' )
!
!  Compute the order dependent code.
!
  call dg_order_code ( adj2, nnode2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Permuted node ordering:' )

  call dg_code_print ( nnode2, code2, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call dg_code_back ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call dg_code_print ( nnode2, code2, '  The maximal code:' )
!
!  Compare the two codes.
!
  call dg_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  SUCCESS: CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE2 < CODE1'
  end if

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests DG_CODE_BRUTE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) code(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) order(8)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  For a digraph code:'
  write ( *, '(a)' ) '  DG_CODE_BRUTE uses brute force;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute the digraph code of'
  write ( *, '(a)' ) '  the cube digraph by brute force.'
!
!  Set up the cube graph.
!
  call g_example_cube ( adj, nnode )

  call g_print ( adj, nnode, '  The graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode
    order(i) = i
  end do

  call dg_order_code ( adj, nnode, nnode, code, order )

  call node_order_print ( nnode, order, '  Initial node ordering:')

  call dg_code_print ( nnode, code, '  The order-dependent code:' )
!
!  Compute the maximal code by brute force.
!
  call dg_code_brute ( adj, nnode, code, order )

  call node_order_print ( nnode, order, '  Induced node ordering:' )

  call dg_code_print ( nnode, code, '  The maximal code:' )

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests DG_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 13

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  character results(test_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  DG_COMPARE compares two digraphs.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare all pairs of test graphs.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    example = i
    call dg_example_octo ( example, adj1, nnode1, seed )

    do j = i, test_num

      example = j
      call dg_example_octo ( example, adj2, nnode2, seed )

      call dg_compare ( adj1, nnode1, adj2, nnode2, order1, &
        order2, result )

      if ( i == j .and. result /= 0 ) then
        write ( *, '(a,i6)' ) '  FAILURE on graph ', i
      else if ( i /= j .and. result == 0 ) then
        write ( *, '(a,2i6)' ) '  FAILURE on graphs ',i, j
      end if

      if ( result < 0 ) then
        results(j,i) = '.'
        write ( results(i,j), '(i1)' ) abs ( result )
      else if ( 0 < result ) then
        results(i,j) = '.'
        write ( results(j,i), '(i1)' ) result
      else if ( result == 0 ) then
        results(i,j) = '='
        results(j,i) = '='
      end if

    end do

  end do

  do i = 1, test_num
    write ( *, '(2x,i2,2x,60a)' ) i, results(i,1:test_num)
  end do

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests G_CODE_BACK;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) code(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(8)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  For a graph code:'
  write ( *, '(a)' ) '  G_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute the graph code of'
  write ( *, '(a)' ) '  the cube graph by backtracking.'
!
!  Set up the cube graph.
!
  call g_example_cube ( adj, nnode )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of nodes is ', nnode

  call g_print ( adj, nnode, '  The graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode
    order(i) = i
  end do

  call g_order_code ( adj, nnode, npart, code, order )

  call node_order_print ( nnode, order, '  Initial node ordering:' )

  call g_code_print ( nnode, code, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call g_code_back ( adj, nnode, code, order )

  call node_order_print ( nnode, order, '  Induced node ordering:' )

  call g_code_print ( nnode, code, '  The maximal code:' )

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests G_CODE_BACK;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  For a graph code:'
  write ( *, '(a)' ) '  G_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compare the graph codes of'
  write ( *, '(a)' ) '  the cube graph, and a node-reordered copy of'
  write ( *, '(a)' ) '  the cube graph.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The codes should be the same.'
!
!  Set up the cube graph.
!
  call g_example_cube ( adj1, nnode1 )

  call g_print ( adj1, nnode1, '  The graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode1
    order1(i) = i
  end do

  call g_order_code ( adj1, nnode1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Initial node ordering:' )

  call g_code_print ( nnode1, code1, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call g_code_back ( adj1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Induced node ordering:' )

  call g_code_print ( nnode1, code1, '  The maximal code:' )
!
!  Now permute the nodes of graph 1 to get graph 2, get its code and print it.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now permute the original graph:'

  nnode2 = nnode1
!
!  Get a random permutation.
!
  do i = 1, nnode2
    order2(i) = i
  end do

  call i4vec_perm_random ( nnode2, order2, seed )
!
!  Reorder the nodes of the graph.
!
  do i = 1, nnode2
    ip = order2(i)
    do j = 1, nnode2
      jp = order2(j)
      adj2(i,j) = adj1(ip,jp)
    end do
  end do

  call g_print ( adj2, nnode2, '  The graph:' )
!
!  Compute the order dependent code.
!
  call g_order_code ( adj2, nnode2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Permuted node ordering:' )

  call g_code_print ( nnode2, code2, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call g_code_back ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call g_code_print ( nnode2, code2, '  The maximal code:' )
!
!  Compare the two codes.
!
  call g_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  SUCCESS: CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE2 < CODE1'
  end if

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests G_CODE_BACK;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  For a graph code:'
  write ( *, '(a)' ) '  G_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compare the graph codes'
  write ( *, '(a)' ) '  of the cube graph and a random graph.'
!
!  Set up the cube graph.
!
  call g_example_cube ( adj1, nnode1 )

  call g_print ( adj1, nnode1, '  The graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode1
    order1(i) = i
  end do

  call g_order_code ( adj1, nnode1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Initial node ordering:' )

  call g_code_print ( nnode1, code1, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call g_code_back ( adj1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Induced node ordering:' )

  call g_code_print ( nnode1, code1, '  The maximal code:' )
!
!  Now choose a random graph on N1 nodes with the same number of links.
!
  nnode2 = nnode1

  call g_edge_count ( adj1, nnode1, nedge )
!
!  Get the random graph.
!
  call g_random ( adj2, nnode2, nedge, seed )

  call g_print ( adj2, nnode2, '  The graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode2
    order2(i) = i
  end do

  call g_order_code ( adj2, nnode2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Initial node ordering:' )

  call g_code_print ( nnode2, code2, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call g_code_back ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call g_code_print ( nnode2, code2, '  The maximal code:' )
!
!  Compare the graph codes.
!
  call g_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)' ) '  CODE2 < CODE1'
  end if

  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests G_CODE_BRUTE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) code(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) order(8)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  For a graph code:'
  write ( *, '(a)' ) '  G_CODE_BRUTE uses brute force;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute the graph code of'
  write ( *, '(a)' ) '  the cube graph by brute force.'
!
!  Set up the cube graph.
!
  call g_example_cube ( adj, nnode )

  call g_print ( adj, nnode, '  The graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode
    order(i) = i
  end do

  call g_order_code ( adj, nnode, nnode, code, order )

  call node_order_print ( nnode, order, '  Initial node ordering:' )

  call g_code_print ( nnode, code, '  The order-dependent code:' )
!
!  Compute the maximal code by brute force.
!
  call g_code_brute ( adj, nnode, code, order )

  call node_order_print ( nnode, order, '  Induced node ordering:' )

  call g_code_print ( nnode, code, '  The maximal code:' )

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests G_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  character results(test_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018'
  write ( *, '(a)' ) '  G_COMPARE compares two graphs.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare all pairs of test graphs.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    example = i
    call g_example_octo ( example, adj1, nnode1, seed )

    do j = i, test_num

      example = j
      call g_example_octo ( example, adj2, nnode2, seed )

      call g_compare ( adj1, nnode1, adj2, nnode2, order1, order2, &
        result )

      if ( i == j .and. result /= 0 ) then
        write ( *, '(a,i6)' ) '  FAILURE on graph ', i
      else if ( i /= j .and. result == 0 ) then
        write ( *, '(a,2i6)' ) '  FAILURE on graphs ',i, j
      end if

      if ( result < 0 ) then
        results(j,i) = '.'
        write ( results(i,j), '(i1)' ) abs ( result )
      else if ( 0 < result ) then
        results(i,j) = '.'
        write ( results(j,i), '(i1)' ) result
      else if ( result == 0 ) then
        results(i,j) = '='
        results(j,i) = '='
      end if

    end do

  end do

  do i = 1, test_num
    write ( *, '(2x,i2,2x,60a)' ) i, results(i,1:test_num)
  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests G_CODE_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  G_CODE_COMPARE compares two graph codes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare the codes of the cube and'
  write ( *, '(a)' ) '  the cube with permuted nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The codes should be the same.'
!
!  Set the graph to the color cube.
!
  call g_example_cube ( adj1, nnode1 )

  call g_print ( adj1, nnode1, '  The graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode1
    order1(i) = i
  end do

  call g_order_code ( adj1, nnode1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Initial node ordering:' )

  call g_code_print ( nnode1, code1, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call g_code_back ( adj1, nnode1, code1, order1 )

  call node_order_print ( nnode1, order1, '  Induced node ordering:' )

  call g_code_print ( nnode1, code1, '  The maximal code:' )
!
!  Now permute the colors of graph 1 to get graph 2, get its code and print it.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Graph 2 is made by permuting graph 1.'

  nnode2 = nnode1

  call perm_random ( nnode2, order2, seed )

  adj2(1:nnode2,1:nnode2) = adj1(1:nnode2,1:nnode2)

  call i4mat_perm ( nnode2, adj2, order2 )
!
!  Compute the order dependent code.
!
  do i = 1, nnode2
    order2(i) = i
  end do

  call g_order_code ( adj2, nnode2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Initial node ordering:' )

  call g_code_print ( nnode2, code2, '  The order-dependent code:' )
!
!  Compute the maximal code by backtracking.
!
  call g_code_back ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call g_code_print ( nnode2, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call g_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  SUCCESS: CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE2 < CODE1'
  end if
!
!  Compute the maximal code by brute force.
!
  do i = 1, nnode2
    order2(i) = i
  end do

  call g_code_brute ( adj2, nnode2, code2, order2 )

  call node_order_print ( nnode2, order2, '  Induced node ordering:' )

  call g_code_print ( nnode2, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call g_code_compare ( code1, code2, nnode1, nnode2, result )

  write ( *, '(a)' ) ' '
  if ( result == -1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE1 < CODE2'
  else if ( result == 0 ) then
    write ( *, '(a)' ) '  SUCCESS: CODE1 = CODE2'
  else if ( result == +1 ) then
    write ( *, '(a)' ) '  FAILURE: CODE2 < CODE1'
  end if

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests MG_CODE_BACK, MG_CODE_BRUTE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  For a multigraph code:'
  write ( *, '(a)' ) '  MG_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) '  MG_CODE_BRUTE uses brute force;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The results should be the same.'
!
!  Set up a random multigraph.
!
  nedge = 25
  nnode = 8

  call mg_random ( adj, nnode, nedge, seed )

  call mg_print ( adj, nnode, '  The multigraph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode
    order1(i) = i
  end do

  call mg_order_code ( adj, nnode, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Initial node ordering:' )

  call mg_code_print ( nnode, code1, '  The order-dependent code:' )
!
!  Brute force.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BRUTE FORCE calculation:'

  call mg_code_brute ( adj, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Induced node ordering:' )

  call mg_code_print ( nnode, code1, '  The maximal code:' )
!
!  Backtrack.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BACKTRACK calculation:'

  call mg_code_back ( adj, nnode, code2, order2 )

  call node_order_print ( nnode, order2, '  Induced node ordering:' )

  call mg_code_print ( nnode, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call mg_code_compare ( code1, code2, nnode, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SUCCESS: The codes are equal.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FAILURE: The codes are unequal.'
  end if

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests MG_ADJ_MAX_MAX, MG_ADJ_MAX_SEQ, MG_ADJ_SEQ;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) adj_max_max
  integer ( kind = 4 ) adj_max_seq(8)
  integer ( kind = 4 ) adj_seq(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  For a multigraph,'
  write ( *, '(a)' ) '  MG_ADJ_MAX_MAX computes the adjacency'
  write ( *, '(a)' ) '    maximum maximum;'
  write ( *, '(a)' ) '  MG_ADJ_MAX_SEQ computes the adjacency'
  write ( *, '(a)' ) '    maximum sequence;'
  write ( *, '(a)' ) '  MG_ADJ_SEQ computes the adjacency'
  write ( *, '(a)' ) '    sequence;'
!
!  Set up a random multigraph.
!
  nedge = 25
  nnode = 8

  call mg_random ( adj, nnode, nedge, seed )

  call mg_print ( adj, nnode, '  The multigraph:' )

  call mg_adj_max_max ( adj, nnode, adj_max_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The adjacency maximum maximum = ', adj_max_max

  call mg_adj_max_seq ( adj, nnode, adj_max_seq )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The adjacency maximum sequence:'
  write ( *, '(a)' ) ' '
  do i = 1, nnode
    write ( *, '(2i4)' ) i, adj_max_seq(i)
  end do

  call mg_adj_seq ( adj, nnode, adj_seq )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The adjacency sequence:'
  write ( *, '(a)' ) ' '
  do i = 1, nnode
    write ( *, '(10i4)' ) adj_seq(i,1:nnode)
  end do

  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests MG_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  character results(test_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  MG_COMPARE compares two multigraphs.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare pairs of test graphs.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    example = i
    call mg_example_octo ( example, adj1, nnode1, seed )

    do j = i, test_num

      example = j
      call mg_example_octo ( example, adj2, nnode2, seed )

      call mg_compare ( adj1, nnode1, adj2, nnode2, order1, &
        order2, result )

      if ( ( i == j .and. result /= 0 ) .or. &
           ( i /= j .and. result == 0 ) ) then

        write ( *, '(a)') ' '
        write ( *, '(a,2i6)' ) '  FAILURE:'
        write ( *, '(a,2i6)' ) '    Graph #1 = example  ', i
        write ( *, '(a,2i6)' ) '    Graph #2 = example  ', j
        write ( *, '(a,2i6)' ) '    Comparison result = ', result
        call mg_print ( adj1, nnode1, '  Multigraph #1:' )
        call mg_print ( adj2, nnode2, '  Multigraph #2:' )

      end if

      if ( result < 0 ) then
        results(j,i) = '.'
        write ( results(i,j), '(i1)' ) abs ( result )
      else if ( 0 < result ) then
        results(i,j) = '.'
        write ( results(j,i), '(i1)' ) result
      else if ( result == 0 ) then
        results(i,j) = '='
        results(j,i) = '='
      end if

    end do

  end do

  write ( *, '(a)' ) ' '

  do i = 1, test_num
    write ( *, '(2x,i2,2x,60a)' ) i, results(i,1:test_num)
  end do

  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests DMG_CODE_BACK, DMG_CODE_BRUTE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  For a dimultigraph code:'
  write ( *, '(a)' ) '  DMG_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) '  DMG_CODE_BRUTE uses brute force;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The results should be the same.'
!
!  Set up a random dimultigraph.
!
  nedge = 25
  nnode = 8

  call dmg_random ( adj, nnode, nedge, seed )

  call dmg_print ( adj, nnode, '  DM-graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode
    order1(i) = i
  end do

  call dmg_order_code ( adj, nnode, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Initial node ordering:' )

  call dmg_code_print ( nnode, code1, '  The order-dependent code:' )
!
!  Brute force.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BRUTE FORCE calculation:'

  call dmg_code_brute ( adj, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Induced node ordering:' )

  call dmg_code_print ( nnode, code1, '  The maximal code:' )
!
!  Backtrack.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BACKTRACK calculation:'

  call dmg_code_back ( adj, nnode, code2, order2 )

  call node_order_print ( nnode, order2, '  Induced node ordering:' )

  call dmg_code_print ( nnode, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call dmg_code_compare ( code1, code2, nnode, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SUCCESS: The codes are equal.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FAILURE: The codes are unequal.'
  end if

  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests DMG_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character, dimension ( 12 ) :: list = (/ &
    '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C' /)
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  character results(test_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024'
  write ( *, '(a)' ) '  DMG_COMPARE compares two dimultigraphs.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare pairs of test graphs.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    example = i
    call dmg_example_octo ( example, adj1, nnode1, seed )

    do j = i, test_num

      example = j
      call dmg_example_octo ( example, adj2, nnode2, seed )

      call dmg_compare ( adj1, nnode1, adj2, nnode2, order1, &
        order2, result )

      if ( ( i == j .and. result /= 0 ) .or. &
           ( i /= j .and. result == 0 ) ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a,2i6)' ) '  FAILURE:'
        write ( *, '(a,2i6)' ) '    Graph #1 = example ', i
        write ( *, '(a,2i6)' ) '    Graph #2 = example ', j
        write ( *, '(a,i6)'  ) '    Comparison =       ', result

        call dmg_print ( adj1, nnode1, '  DM-graph #1:' )
        call dmg_print ( adj2, nnode2, '  DM-graph #2:' )

      end if

      if ( result < 0 ) then
        results(j,i) = '.'
        results(i,j) = list ( abs ( result ) )
      else if ( 0 < result ) then
        results(i,j) = '.'
        results(j,i) = list ( result )
      else if ( result == 0 ) then
        results(i,j) = '='
        results(j,i) = '='
      end if

    end do

  end do

  write ( *, '(a)' ) ' '
  do i = 1, test_num
    write ( *, '(2x,i2,2x,60a)' ) i, ( results(i,j), j = 1, test_num )
  end do

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests CDMG_CODE_BACK, CDMG_CODE_BRUTE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) code1(8,8)
  integer ( kind = 4 ) code2(8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  For a color dimultigraph code:'
  write ( *, '(a)' ) '  CDMG_CODE_BACK uses backtracking;'
  write ( *, '(a)' ) '  CDMG_CODE_BRUTE uses brute force;'
  write ( *, '(a)' ) ' '
!
!  Set up a random color dimultigraph.
!
  ncolor = 5
  nedge = 25
  nnode = 8

  call cdmg_random ( adj, nnode, ncolor, nedge, seed )

  call cdmg_print ( adj, nnode, '  The random CDM-graph:' )
!
!  Compute the order dependent code.
!
  do i = 1, nnode
    order1(i) = i
  end do

  call cdmg_order_code ( adj, nnode, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Initial node ordering:' )

  call cdmg_code_print ( nnode, code1, '  The order-dependent code:' )
!
!  Brute force.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BRUTE FORCE calculation:'

  call cdmg_code_brute ( adj, nnode, code1, order1 )

  call node_order_print ( nnode, order1, '  Induced node ordering:' )

  call cdmg_code_print ( nnode, code1, '  The maximal code:' )
!
!  Backtrack.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BACKTRACK calculation:'

  call cdmg_code_back ( adj, nnode, code2, order2 )

  call node_order_print ( nnode, order2, '  Induced node ordering:' )

  call cdmg_code_print ( nnode, code2, '  The maximal code:' )
!
!  Compare the codes.
!
  call cdmg_code_compare ( code1, code2, nnode, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SUCCESS: The codes are equal.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FAILURE: he codes are unequal.'
  end if

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests CDMG_COMPARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 12

  integer ( kind = 4 ) adj1(8,8)
  integer ( kind = 4 ) adj2(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character, dimension ( 15 ) :: list = (/ &
    '1', '2', '3', '4', '5', '6', '7', '8', '9', &
    'A', 'B', 'C', 'D', 'E', 'F' /)
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) order1(8)
  integer ( kind = 4 ) order2(8)
  integer ( kind = 4 ) result
  character results(test_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  CDMG_COMPARE compares two color dimultigraphs.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare pairs of test graphs.'

  do i = 1, test_num

    example = i
    call cdmg_example_octo ( example, adj1, nnode1, seed )

    do j = i, test_num

      example = j
      call cdmg_example_octo ( example, adj2, nnode2, seed )

      call cdmg_compare ( adj1, nnode1, adj2, nnode2, order1, &
        order2, result )

      if ( ( i == j .and. result /= 0 ) .or. &
           ( i /= j .and. result == 0 ) ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a,2i6)' ) '  FAILURE:'
        write ( *, '(a,2i6)' ) '    Graph #1 = example ', i
        write ( *, '(a,2i6)' ) '    Graph #2 = example ', j
        write ( *, '(a,i6)'  ) '    Comparison =       ', result

        call cdmg_print ( adj1, nnode1, '  CDM-graph #1:' )
        call cdmg_print ( adj2, nnode2, '  CDM-graph #2:' )

      end if

      if ( result < 0 ) then
        results(j,i) = '.'
        results(i,j) = list ( abs ( result ) )
      else if ( 0 < result ) then
        results(i,j) = '.'
        results(j,i) = list ( result )
      else if ( result == 0 ) then
        results(i,j) = '='
        results(j,i) = '='
      end if

    end do

  end do

  write ( *, '(a)' ) ' '

  do i = 1, test_num
    write ( *, '(2x,i2,2x,60a)' ) i, results(i,1:test_num)
  end do

  return
end
