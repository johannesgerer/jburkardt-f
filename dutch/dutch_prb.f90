program main

!*****************************************************************************80
!
!! MAIN is the main program for DUTCH.
!
!  Discussion:
!
!    DUTCH_PRB tests the DUTCH routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DUTCH_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DUTCH library.'

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
  write ( *, '(a)' ) 'DUTCH_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests POINTS_CONVEX_HULL_CUBIC_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 7

  integer ( kind = 4 ) hull_num
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) hull(node_num)
  real ( kind = 8 ) hull_xy(2,node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  POINTS_CONVEX_HULL_CUBIC_2D computes the convex hull'
  write ( *, '(a)' ) '  of a set of N 2D points, using an algorithm '
  write ( *, '(a)' ) '  that is cubic in N.'

  node_xy(1:2,1:node_num) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, &
    2.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 3.0D+00, &
    2.0D+00, 2.0D+00 /), (/ 2, node_num /) )

  call r8mat_transpose_print ( 2, node_num, node_xy, &
    '  Coordinates of the points:' )

  call points_convex_hull_cubic_2d ( node_num, node_xy, hull_num, hull )

  hull_xy(1:2,1:hull_num) = node_xy(1:2,hull(1:hull_num))

  call r8mat_transpose_print ( 2, hull_num, hull_xy, &
    '  Coordinates of the convex hull:' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests POINTS_CONVEX_HULL_NLOGN_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 7

  integer ( kind = 4 ) hull_num
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) hull(node_num)
  real ( kind = 8 ) hull_xy(2,node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  POINTS_CONVEX_HULL_NLOGN_2D computes the convex hull '
  write ( *, '(a)' ) '  of a set of N 2D points using an algorithm'
  write ( *, '(a)' ) '  that is NlogN in N.'

  node_xy(1:2,1:node_num) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, &
    2.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 3.0D+00, &
    2.0D+00, 2.0D+00 /), (/ 2, node_num /) )

  call r8mat_transpose_print ( 2, node_num, node_xy, &
    '  Coordinates of the points:' )

  call points_convex_hull_nlogn_2d ( node_num, node_xy, hull_num, hull )

  hull_xy(1:2,1:hull_num) = node_xy(1:2,hull(1:hull_num))

  call r8mat_transpose_print ( 2, hull_num, hull_xy, &
    '  Coordinates of the convex hull:' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests POINTS_CONVEX_HULL_NLOGH_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 7

  integer ( kind = 4 ) hull_num
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) hull(node_num)
  real ( kind = 8 ) hull_xy(2,node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  POINTS_CONVEX_HULL_NLOGH_2D computes the convex hull'
  write ( *, '(a)' ) '  of a set of N 2D points using an algorithm'
  write ( *, '(a)' ) '  that is order NlogH.'
  write ( *, '(a)' ) '  (H is the number of points on the convex hull.)'

  node_xy(1:2,1:node_num) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, &
    2.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 3.0D+00, &
    2.0D+00, 2.0D+00 /), (/ 2, node_num /) )

  call r8mat_transpose_print ( 2, node_num, node_xy, &
    '  Coordinates of the points:' )

  call points_convex_hull_nlogh_2d ( node_num, node_xy, hull_num, hull )

  hull_xy(1:2,1:hull_num) = node_xy(1:2,hull(1:hull_num))

  call r8mat_transpose_print ( 2, hull_num, hull_xy, &
    '  Coordinates of the convex hull:' )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests POLYCON_MINKOWSKI_SUM_LINEAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nu = 4
  integer ( kind = 4 ), parameter :: nv = 3

  integer ( kind = 4 ) nw
  real ( kind = 8 ) ux(nu)
  real ( kind = 8 ) uy(nu)
  real ( kind = 8 ) vx(nv)
  real ( kind = 8 ) vy(nv)
  real ( kind = 8 ) wx(nu+nv+2)
  real ( kind = 8 ) wy(nu+nv+2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  POLYCON_MINKOWSKI_SUM_LINEAR computes the Minkowski sum'
  write ( *, '(a)' ) '  of two convex polygons using a linear algorithm.'

  ux(1) = 0.0D+00
  uy(1) = 0.0D+00
  ux(2) = 2.0D+00
  uy(2) = 2.0D+00
  ux(3) = -1.0D+00
  uy(3) = 3.0D+00
  ux(4) = -2.0D+00
  uy(4) = 2.0D+00

  call r8vec2_print ( nu, ux, uy, '  Coordinates of polygon U:' )

  vx(1) = 8.0D+00
  vy(1) = 2.0D+00
  vx(2) = 9.0D+00
  vy(2) = 5.0D+00
  vx(3) = 7.0D+00
  vy(3) = 4.0D+00

  call r8vec2_print ( nv, vx, vy, '  Coordinates of polygon V:' )

  call polycon_minkowski_sum_linear ( nu, ux, uy, nv, vx, vy, nw, wx, wy )

  call r8vec2_print ( nw, wx, wy, '  Coordinates of Minkowski sum polygon W:' )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests POLYCON_MINKOWSKI_SUM_LINEAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nu = 4
  integer ( kind = 4 ), parameter :: nv = 3

  integer ( kind = 4 ) nw
  real ( kind = 8 ) ux(nu)
  real ( kind = 8 ) uy(nu)
  real ( kind = 8 ) vx(nv)
  real ( kind = 8 ) vy(nv)
  real ( kind = 8 ) w_xy(2,nu+nv+2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  POLYCON_MINKOWSKI_SUM_N2LOGN2 computes the Minkowski'
  write ( *, '(a)' ) '  sum of two convex polygons using a linear algorithm.'

  ux(1) = 0.0D+00
  uy(1) = 0.0D+00
  ux(2) = 2.0D+00
  uy(2) = 2.0D+00
  ux(3) = -1.0D+00
  uy(3) = 3.0D+00
  ux(4) = -2.0D+00
  uy(4) = 2.0D+00

  call r8vec2_print ( nu, ux, uy, '  Coordinates of polygon U:' )

  vx(1) = 8.0D+00
  vy(1) = 2.0D+00
  vx(2) = 9.0D+00
  vy(2) = 5.0D+00
  vx(3) = 7.0D+00
  vy(3) = 4.0D+00

  call r8vec2_print ( nv, vx, vy, '  Coordinates of polygon V:' )

  call polycon_minkowski_sum_n2logn2 ( nu, ux, uy, nv, vx, vy, nw, w_xy )

  call r8mat_transpose_print ( 2, nw, w_xy, &
    '  Coordinates of Minkowski sum polygon W:' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests PERM_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  PERM_RANDOM produces a random permutation;'
  write ( *, '(a,i6)' ) '  For this test, N = ', n
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call perm_random ( n, seed, p )
    call perm_print ( n, p, ' ' )
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests POINTS_MINIDISC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  real ( kind = 8 ) px(n)
  real ( kind = 8 ) py(n)
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  POINTS_MINIDISC computes the smallest circle that'
  write ( *, '(a)' ) '  contains a set of N 2D points.'

  px(1) = 0.0D+00
  py(1) = 0.0D+00
  px(2) = 1.0D+00
  py(2) = 2.0D+00
  px(3) = 2.0D+00
  py(3) = 0.0D+00
  px(4) = 1.0D+00
  py(4) = 1.0D+00
  px(5) = 0.0D+00
  py(5) = 2.0D+00
  px(6) = 1.0D+00
  py(6) = 3.0D+00
  px(7) = 2.0D+00
  py(7) = 2.0D+00

  call r8vec2_print ( n, px, py, '  Coordinates of the points:' )

  call points_minidisc_2d ( n, px, py, r, cx, cy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The enclosing circle has:'
  write ( *, '(a,g14.6)' ) '  Radius R =      ', r
  write ( *, '(a,2g14.6)' ) '  Center (CX,CY)= ', cx, cy

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests LINE_SEG_VEC_INT_2D.
!
!  Discussion:
!
!    We are expecting:
!
!    2, 4, 2.0, 1.0
!    2, 8, 2.0, 1.0
!    4, 6, 3.0, 2.0
!    4, 8, 2.0, 1.0
!    6, 8, 2.0, 3.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) xint
  real ( kind = 8 ), dimension ( n ) :: x1 = &
    (/ 4.0D+00, 0.0D+00, 4.0D+00, 1.0D+00, 0.0D+00, 1.0D+00, 3.0D+00, 2.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: x2 = &
    (/ 5.0D+00, 3.0D+00, 5.0D+00, 5.0D+00, 1.0D+00, 4.0D+00, 3.0D+00, 2.0D+00 /)
  real ( kind = 8 ) yint
  real ( kind = 8 ), dimension ( n ) :: y1 = &
    (/ 0.0D+00, 1.0D+00, 2.0D+00, 0.0D+00, 3.0D+00, 4.0D+00, 4.0D+00, 4.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: y2 = &
    (/ 0.0D+00, 1.0D+00, 3.0D+00, 4.0D+00, 2.0D+00, 1.0D+00, 3.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  LINE_SEG_VEC_INT_2D finds the intersections of a set'
  write ( *, '(a)' ) '  of line segments.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, X1(I), Y1(I), X2(I), Y2(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i2,4f7.1)' ) i, x1(i), y1(i), x2(i), y2(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Seg#1, Seg#2, Intersection'
  write ( *, '(a)' ) ' '

  i = 0
  j = 0

  do

    call line_seg_vec_int_2d ( n, x1, y1, x2, y2, i, j, flag, xint, yint )

    if ( flag == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more intersections.'
      exit
    end if

    write ( *, '(2i6,2g14.6)' ) i, j, xint, yint

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests RECT_INT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntest = 8

  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ), dimension ( ntest ) :: x3 = &
    (/ 1.0D+00, 3.0D+00, -2.0D+00, 3.0D+00, 2.0D+00, &
       5.0D+00, 5.0D+00, 6.0D+00 /)
  real ( kind = 8 ), dimension ( ntest ) :: x4 = &
    (/ 2.0D+00, 5.0D+00,  1.0D+00, 9.0D+00, 8.0D+00, &
       8.0D+00, 8.0D+00, 8.0D+00 /)
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ), dimension ( ntest ) :: y3 = (/ &
    2.0D+00, 2.0D+00, 2.0D+00, -2.0D+00, -3.0D+00, &
    1.0D+00, 5.0D+00, 7.0D+00 /)
  real ( kind = 8 ), dimension ( ntest ) :: y4 = (/ &
    4.0D+00, 5.0D+00, 4.0D+00,  1.0D+00,  8.0D+00, &
    3.0D+00, 8.0D+00, 9.0D+00 /)
  real ( kind = 8 ) y5
  real ( kind = 8 ) y6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  RECT_INT_2D finds the intersections of two rectangles.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For all tests, rectangle #1 is (0,0), (5,5).'
  write ( *, '(a)' ) ' '

  x1 = 0.0D+00
  y1 = 0.0D+00

  x2 = 5.0D+00
  y2 = 5.0D+00

  do i = 1, ntest

    call rect_int_2d ( x1, y1, x2, y2, x3(i), y3(i), x4(i), y4(i), &
      flag, x5, y5, x6, y6 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Second rectangle:'
    write ( *, '(a)' ) ' '
    write ( *, '(2f8.4)' ) x3(i), y3(i)
    write ( *, '(2f8.4)' ) x4(i), y4(i)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Intersection rectangle:'
    write ( *, '(a)' ) ' '

    if ( flag == 0 ) then
      write ( *, '(a)' ) '  EMPTY'
    else
      write ( *, '(2f8.4)' ) x5, y5
      write ( *, '(2f8.4)' ) x6, y6
    end if

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests POLY_REORDER_NODES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: npoly = 7
  integer ( kind = 4 ), parameter :: nxy = 8

  integer ( kind = 4 ), dimension ( npoly ) :: poly = (/ 7, 8, 1, 6, 2, 5, 4 /)
  real ( kind = 8 ), dimension ( nxy ) :: x = &
    (/ 4.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, 2.0D+00, 3.0D+00, 5.0D+00 /)
  real ( kind = 8 ), dimension ( nxy ) :: y = &
    (/ 3.0D+00, 3.0D+00, 0.0D+00, 0.0D+00, 2.0D+00, 5.0D+00, 1.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  POLY_REORDER_NODES reorders the nodes of a polygon'
  write ( *, '(a)' ) '  so that node 1 is leftmost (and lowest, in ties).'

  call r8vec2_print ( nxy, x, y, '  The nodes:' )

  call i4vec_print ( npoly, poly, '  The sequence of nodes:' )

  call poly_reorder_nodes ( nxy, x, y, npoly, poly )

  call i4vec_print ( npoly, poly, '  The reordered sequence of nodes:' )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests POLY_TRIANGULATE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 13

  integer ( kind = 4 ) triang(3,n-2)
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    7.0D+00, 7.0D+00, 9.0D+00,  4.0D+00, 4.0D+00, &
    6.0D+00, 6.0D+00, 4.0D+00,  2.0D+00, 2.0D+00, &
    0.0D+00, 2.0D+00, 4.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: y = (/ &
    0.0D+00, 4.0D+00, 5.0D+00, 10.0D+00, 4.0D+00, &
    6.0D+00, 2.0D+00, 2.0D+00,  4.0D+00, 2.0D+00, &
    2.0D+00, 0.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  POLY_TRIANGULATE_2D triangulates a polygon.'

  call r8vec2_print ( n, x, y, '  The nodes of the polygon:' )

  call poly_triangulate_2d ( n, x, y, triang )

  call i4mat_transpose_print ( 3, n-2, triang, '  The triangulation:' )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests TRIANGULATE_TRICOLOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 13

  integer ( kind = 4 ) color(n)
  integer ( kind = 4 ), dimension ( 3, n-2 ) :: triang =  reshape ( (/ &
    11,    12,    10, &
    12,    13,    10, &
    10,    13,     9, &
     9,    13,     8, &
    13,     1,     8, &
     8,     1,     7, &
     5,     6,     4, &
     4,     6,     3, &
     7,     1,     6, &
     6,     2,     3, &
     6,     1,     2 /), (/ 3, n-2 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  TRIANGULATE_TRICOLOR tricolors a triangulation.'

  call i4mat_transpose_print ( 3, n-2, triang, '  The triangulation:' )

  call triangulate_tricolor ( n, triang, color )

  call i4vec_print ( n, color, '  The node coloring' )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests TRIANGULATION_BOUNDARY_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) bound_num
  integer ( kind = 4 ), parameter :: point_num = 13
  integer ( kind = 4 ), parameter :: tri_num = 16

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  TRIANGULATION_BOUNDARY_COUNT determines the number of'
  write ( *, '(a)' ) '  edges that lie on the convex hull of a region that has'
  write ( *, '(a)' ) '  been triangulated.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of points =         ', point_num
  write ( *, '(a,i6)' ) '  Number of triangles =      ', tri_num

  call triangulation_boundary_count ( point_num, tri_num, bound_num )

  write ( *, '(a,i6)' ) '  Number of boundary edges = ', bound_num

  return
end
