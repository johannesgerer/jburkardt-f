program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM2D_PACK_PRB.
!
!  Discussion:
!
!    FEM2D_PACK_PRB calls the various FEM2D_PACK tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_PACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FEM2D_PACK library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test065 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test105 ( )
  call test107 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_PACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BANDWIDTH_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension(:,:) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  BANDWIDTH_MESH computes the geometric bandwidth:'
  write ( *, '(a)' ) '  of a finite element mesh.'

  nelemx = 2
  nelemy = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NELEMX = ', nelemx
  write ( *, '(a,i8)' ) '  NELEMY = ', nelemy

  element_order = 6
  call grid_element_num ( 'T6', nelemx, nelemy, element_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ELEMENT_ORDER = ', element_order
  write ( *, '(a,i8)' ) '  ELEMENT_NUM   = ', element_num

  allocate ( element_node(1:element_order,1:element_num) )

  call grid_t6_element ( nelemx, nelemy, element_node )

  call grid_print ( element_order, element_num, element_node )

  call bandwidth_mesh ( element_order, element_num, element_node, &
    ml, mu, m )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
  write ( *, '(a,i8)' ) '  Total bandwidth M  = ', m

  deallocate ( element_node )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BANDWIDTH_VAR, NS_T6_VAR_COUNT, NS_T6_VAR_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension(:,:) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: var
  integer ( kind = 4 ), allocatable, dimension ( : ) :: var_node
  integer ( kind = 4 ) var_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For the Navier Stokes variables associated with'
  write ( *, '(a)' ) '  a mesh of T6 elements,'
  write ( *, '(a)' ) '  NS_T6_VAR_COUNT counts variables,'
  write ( *, '(a)' ) '  NS_T6_VAR_SET sets them,'
  write ( *, '(a)' ) '  BANDWIDTH_VAR computes the variable bandwidth.'

  nelemx = 2
  nelemy = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NELEMX = ', nelemx
  write ( *, '(a,i8)' ) '  NELEMY = ', nelemy

  element_order = 6
  call grid_element_num ( 'T6', nelemx, nelemy, element_num )
  call grid_node_num ( 'T6', nelemx, nelemy, node_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ELEMENT_ORDER = ', element_order
  write ( *, '(a,i8)' ) '  ELEMENT_NUM   = ', element_num
  write ( *, '(a,i8)' ) '  NODE_NUM      = ', node_num

  allocate ( element_node(1:element_order,1:element_num) )

  call grid_t6_element ( nelemx, nelemy, element_node )

  call grid_print ( element_order, element_num, element_node )

  allocate ( var_node(1:node_num+1) )

  call ns_t6_var_count ( element_num, element_node, node_num, var_node, &
    var_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of variables VAR_NUM = ', var_num

  call i4vec_print ( node_num+1, var_node, '  VAR_NODE pointer vector:' )

  allocate ( var(1:var_num) )

  call ns_t6_var_set ( element_num, element_node, node_num, var_node, &
    var_num, var )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node    U_Var    V_Var    P_Var'
  write ( *, '(a)' ) ' '

  do node = 1, node_num

    ilo = var_node(node)
    ihi = var_node(node+1)-1

    write ( *, '(2x,i8,4x,i8,2x,i8,2x,i8)' ) node, var(ilo:ihi)

  end do

  call bandwidth_var ( element_order, element_num, element_node, &
    node_num, var_node, var_num, var, ml, mu, m )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
  write ( *, '(a,i8)' ) '  Total bandwidth M  = ', m

  deallocate ( element_node )
  deallocate ( var )
  deallocate ( var_node )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests BASIS_11_**_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Test the computation of ONE basis function'
  write ( *, '(a)' ) '  at ONE point in a given element:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BASIS_11_Q4_TEST : Q4 element.'
  write ( *, '(a)' ) '  BASIS_11_T3_TEST : T3 element.'
  write ( *, '(a)' ) '  BASIS_11_T4_TEST : T4 element.'
  write ( *, '(a)' ) '  BASIS_11_T6_TEST : T6 element.'

  call basis_11_q4_test ( )

  call basis_11_t3_test ( )

  call basis_11_t4_test ( )

  call basis_11_t6_test ( )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests BASIS_MN_**_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test the computation of all M basis functions'
  write ( *, '(a)' ) '  at N points in a given element:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BASIS_MN_Q4_TEST : Q4 element.'
  write ( *, '(a)' ) '  BASIS_MN_T3_TEST : T3 element.'
  write ( *, '(a)' ) '  BASIS_MN_T4_TEST : T4 element.'
  write ( *, '(a)' ) '  BASIS_MN_T6_TEST : T6 element.'

  call basis_mn_q4_test ( )

  call basis_mn_t3_test ( )

  call basis_mn_t4_test ( )

  call basis_mn_t6_test ( )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 demonstrates DERIVATIVE_AVERAGE_T3.
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

  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), parameter :: nelemx = 7
  integer ( kind = 4 ), parameter :: nelemy = 5

  real    ( kind = 8 ) angle
  real    ( kind = 8 ), allocatable, dimension ( : ) :: c
  integer ( kind = 4 ) col
  real    ( kind = 8 ), allocatable, dimension ( : ) :: dcdx
  real    ( kind = 8 ) dcdx_exact
  real    ( kind = 8 ), allocatable, dimension ( : ) :: dcdy
  real    ( kind = 8 ) dcdy_exact
  real    ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) element
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  real    ( kind = 8 ) r
  integer ( kind = 4 ) row
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  DERIVATIVE_AVERAGE_T3 averages the spatial derivatives'
  write ( *, '(a)' ) '  of a finite element function at the nodes.'
!
!  How many elements are there?
!
  call grid_t3_element_num ( nelemx, nelemy, element_num )

  allocate ( element_node(element_order,element_num) )
!
!  How many nodes are there?
!
  call grid_t3_node_num ( nelemx, nelemy, node_num )

  allocate ( c(node_num) )
  allocate ( dcdx(node_num) )
  allocate ( dcdy(node_num) )
  allocate ( node_xy(2,node_num) )
!
!  Get the nodes that make up each element.
!
  call grid_t3_element ( nelemx, nelemy, element_node )
!
!  Generate the coordinates of the nodes.
!
  node = 0

  do row = 0, nelemy

    r = ( real ( nelemy - row, kind = 8 ) *  1.0D+00   &
        + real (        + row, kind = 8 ) *  3.0D+00 ) &
        / real ( nelemy,       kind = 8 )

    do col = 0, nelemx

      node = node + 1
      angle = ( real ( nelemx - col, kind = 8 ) * 135.0D+00   &
              + real (        + col, kind = 8 ) *  45.0D+00 ) &
              / real ( nelemx,       kind = 8 )

      angle = degrees_to_radians ( angle )

      node_xy(1,node) = r * cos ( angle )
      node_xy(2,node) = r * sin ( angle )

    end do

  end do
!
!  Set the finite element function.
!
  do node = 1, node_num
    x = node_xy(1,node)
    y = node_xy(2,node)
    c(node) = sin ( x ) * ( 1.0D+00 + y * y )
  end do

  call derivative_average_t3 ( node_num, node_xy, element_num, &
    element_node, c, dcdx, dcdy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C         X               Y'
  write ( *, '(a)' ) '         dCdX(computed)  dCdY(computed)'
  write ( *, '(a)' ) '         dCdX(exact)     dCdY(exact)'
  write ( *, '(a)' ) ' '

  do node = 1, node_num

    x = node_xy(1,node)
    y = node_xy(2,node)

    dcdx_exact = cos ( x ) * ( 1.0D+00 * y * y )
    dcdy_exact = sin ( x ) * 2.0D+00 * y

    write ( *, '(a)' ) ' '
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) c(node), node_xy(1:2,node)
    write ( *, '(2x,  14x,2x,g14.6,2x,g14.6)' ) dcdx(node), dcdy(node)
    write ( *, '(2x,  14x,2x,g14.6,2x,g14.6)' ) dcdx_exact, dcdy_exact

  end do

  deallocate ( c )
  deallocate ( dcdx )
  deallocate ( dcdy )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests DIV_Q4.
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

  integer ( kind = 4 ), parameter :: m = 21
  integer ( kind = 4 ), parameter :: n = 13

  real    ( kind = 8 ) div(m-1,n-1)
  real    ( kind = 8 ) dudx
  real    ( kind = 8 ) dudy
  real    ( kind = 8 ) dvdx
  real    ( kind = 8 ) dvdy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real    ( kind = 8 ), dimension (2,4) :: q = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 4 /) )
  real    ( kind = 8 ) u(m,n)
  real    ( kind = 8 ) v(m,n)
  real    ( kind = 8 ) vort(m-1,n-1)
  real    ( kind = 8 ) x(m,n)
  real    ( kind = 8 ) xm(m-1,n-1)
  real    ( kind = 8 ) y(m,n)
  real    ( kind = 8 ) ym(m-1,n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  DIV_Q4 estimates divergence and vorticity'
  write ( *, '(a)' ) '  using 4 node quadrilateral elements.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Original U, V data forms ', m, ' rows and '
  write ( *, '(a,i8,a)' ) '  ', n, ' columns.'

  ilo = 5
  jlo = 5
  ihi = 8
  jhi = 8

  do i = 1, m
    do j = 1, n

      y(i,j) = ( real ( m - i,     kind = 8 ) * q(2,1)   &
               + real (     i - 1, kind = 8 ) * q(2,3) ) &
               / real ( m     - 1, kind = 8 )

      x(i,j) = ( real ( n - j,     kind = 8 ) * q(1,1)   &
               + real (     j - 1, kind = 8 ) * q(1,3) ) &
               / real ( n     - 1, kind = 8 )
    end do
  end do

  call r8mat_print_some ( m, n, x, ilo, jlo, ihi, jhi, '  Some of X:' )
  call r8mat_print_some ( m, n, y, ilo, jlo, ihi, jhi, '  Some of Y:' )
!
!  Put dummy data into U and V at the data nodes.
!
  do i = 1, m
    do j = 1, n
      u(i,j) = x(i,j) * y(i,j)
      v(i,j) = sin ( x(i,j) * x(i,j) + y(i,j) * y(i,j) )
    end do
  end do

  call r8mat_print_some ( m, n, u, ilo, jlo, ihi, jhi, '  Some of U:' )
  call r8mat_print_some ( m, n, v, ilo, jlo, ihi, jhi, '  Some of V:' )
!
!  Get DIV and VORT.
!
  call div_q4 ( m, n, u, v, q, div, vort )
!
!  Compare computed and known values at the centers of the elements.
!
  do i = 1, m-1

    do j = 1, n-1

      ym(i,j) = ( real ( 2 * m - 2 * i - 1, kind = 8 ) * q(2,1)   &
                + real (         2 * i - 1, kind = 8 ) * q(2,3) ) &
                / real ( 2 * m         - 2, kind = 8 )

      xm(i,j) = ( real ( 2 * n - 2 * j - 1, kind = 8 ) * q(1,1)   &
                + real (         2 * j - 1, kind = 8 ) * q(1,3) ) &
                / real ( 2 * n         - 2, kind = 8 )
    end do
  end do

  call r8mat_print_some ( m-1, n-1, xm, ilo, jlo, ihi, jhi, '  Some of XM:' )
  call r8mat_print_some ( m-1, n-1, ym, ilo, jlo, ihi, jhi, '  Some of YM:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J   XM     YM      DIV(I,J)   Exact DIV        Diff'
  write ( *, '(a)' ) ' '

  do i = ilo, ihi

    do j = jlo, jhi

      dudx = ym(i,j)
      dudy = xm(i,j)
      dvdx = 2.0D+00 * xm(i,j) * cos ( xm(i,j) * xm(i,j) + ym(i,j) * ym(i,j) )
      dvdy = 2.0D+00 * ym(i,j) * cos ( xm(i,j) * xm(i,j) + ym(i,j) * ym(i,j) )

      write ( *, '(2i4,2f7.4,3g14.6)' ) &
        i, j, xm(i,j), ym(i,j), div(i,j), dudx + dvdy, &
        div(i,j) - ( dudx + dvdy )

    end do
    write ( *, '(a)' ) ' '
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   I   J   XM     YM      VORT(I,J)  Exact Vort       Diff'
  write ( *, '(a)' ) ' '

  do i = ilo, ihi

    do j = jlo, jhi

      dudx = ym(i,j)
      dudy = xm(i,j)
      dvdx = 2.0D+00 * xm(i,j) * cos ( xm(i,j) * xm(i,j) + ym(i,j) * ym(i,j) )
      dvdy = 2.0D+00 * ym(i,j) * cos ( xm(i,j) * xm(i,j) + ym(i,j) * ym(i,j) )

      write ( *, '(2i4,2f7.4,3g14.6)' ) &
        i, j, xm(i,j), ym(i,j), vort(i,j), dvdx - dudy, &
        vort(i,j) - ( dvdx - dudy )

    end do
    write ( *, '(a)' ) ' '
  end do

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests DIV_T3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 21
  integer ( kind = 4 ), parameter :: n = 13

  real    ( kind = 8 ) div(2,m-1,n-1)
  real    ( kind = 8 ) dudx
  real    ( kind = 8 ) dudy
  real    ( kind = 8 ) dvdx
  real    ( kind = 8 ) dvdy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real    ( kind = 8 ), dimension (2,4) :: q = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 4 /) )
  real    ( kind = 8 ) u(m,n)
  real    ( kind = 8 ) v(m,n)
  real    ( kind = 8 ) vort(2,m-1,n-1)
  real    ( kind = 8 ) x(m,n)
  real    ( kind = 8 ) xc
  real    ( kind = 8 ) y(m,n)
  real    ( kind = 8 ) yc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  DIV_T3 estimates divergence and vorticity'
  write ( *, '(a)' ) '  using 3 node triangular elements.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Original U, V data forms ', m, ' rows and '
  write ( *, '(a,i8,a)' ) '  ', n, ' columns.'

  ilo = 5
  jlo = 5
  ihi = 8
  jhi = 8
!
!  This mapping does NOT require the region defined by Q to be a rectangle.
!
  do i = 1, m
    do j = 1, n

      x(i,j) = ( &
               real ( n - j,     kind = 8 )  * ( real ( m - i,     kind = 8 ) * q(1,1)   &
                                               + real (     i - 1, kind = 8 ) * q(1,2) ) &
                                               / real ( m     - 1, kind = 8 )            &
             + real (     j - 1, kind = 8 )  * ( real ( m - i,     kind = 8 ) * q(1,4)   &
                                               + real (     i - 1, kind = 8 ) * q(1,3) ) &
                                               / real ( m     - 1, kind = 8 )            &
           ) / real ( n     - 1, kind = 8 )

      y(i,j) = ( &
               real ( n - j,     kind = 8 )  * ( real ( m - i,     kind = 8 ) * q(2,1)   &
                                               + real (     i - 1, kind = 8 ) * q(2,2) ) &
                                               / real ( m     - 1, kind = 8 )            &
             + real (     j - 1, kind = 8 )  * ( real ( m - i,     kind = 8 ) * q(2,4)   &
                                               + real (     i - 1, kind = 8 ) * q(2,3) ) &
                                               / real ( m     - 1, kind = 8 )            &
           ) / real ( n     - 1, kind = 8 )

    end do
  end do

  call r8mat_print_some ( m, n, x, ilo, jlo, ihi, jhi, '  Some of X:' )
  call r8mat_print_some ( m, n, y, ilo, jlo, ihi, jhi, '  Some of Y:' )
!
!  Put data into U and V at the data nodes.
!
  do i = 1, m
    do j = 1, n
      u(i,j) = x(i,j) * y(i,j)
      v(i,j) = sin ( x(i,j) * x(i,j) + y(i,j) * y(i,j) )
    end do
  end do

  call r8mat_print_some ( m, n, u, ilo, jlo, ihi, jhi, '  Some of U:' )
  call r8mat_print_some ( m, n, v, ilo, jlo, ihi, jhi, '  Some of V:' )
!
!  Get DIV and VORT estimated by DIV_T3.
!
  call div_t3 ( m, n, u, v, q, div, vort )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J    X      Y      DIV(I,J)   Exact DIV        Diff'
  write ( *, '(a)' ) ' '

  do i = ilo, ihi

    do j = jlo, jhi

      xc = ( x(i,j) + x(i,j+1) + x(i+1,j+1) ) / 3.0D+00
      yc = ( y(i,j) + y(i,j+1) + y(i+1,j+1) ) / 3.0D+00

      dudx = yc
      dudy = xc
      dvdx = 2.0D+00 * xc * cos ( xc * xc + yc * yc )
      dvdy = 2.0D+00 * yc * cos ( xc * xc + yc * yc )

      write ( *, '(2i4,2f7.4,3g14.6)' ) &
        i, j, xc, yc, div(1,i,j), dudx + dvdy, div(1,i,j) - ( dudx + dvdy )

      xc = ( x(i+1,j+1) + x(i+1,j) + x(i,j) ) / 3.0D+00
      yc = ( y(i+1,j+1) + y(i+1,j) + y(i,j) ) / 3.0D+00

      dudx = yc
      dudy = xc
      dvdx = 2.0D+00 * xc * cos ( xc * xc + yc * yc )
      dvdy = 2.0D+00 * yc * cos ( xc * xc + yc * yc )

      write ( *, '(2i4,2f7.4,3g14.6)' ) &
        i, j, xc, yc, div(2,i,j), dudx + dvdy, div(2,i,j) - ( dudx + dvdy )

    end do
    write ( *, '(a)' ) ' '
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   I   J    X      Y      VORT(I,J)  Exact Vort       Diff'
  write ( *, '(a)' ) ' '

  do i = ilo, ihi

    do j = jlo, jhi

      xc = ( x(i,j) + x(i,j+1) + x(i+1,j+1) ) / 3.0D+00
      yc = ( y(i,j) + y(i,j+1) + y(i+1,j+1) ) / 3.0D+00

      dudx = yc
      dudy = xc
      dvdx = 2.0D+00 * xc * cos ( xc * xc + yc * yc )
      dvdy = 2.0D+00 * yc * cos ( xc * xc + yc * yc )

      write ( *, '(2i4,2f7.4,3g14.6)' ) &
        i, j, xc, yc, vort(1,i,j), dvdx - dudy, vort(1,i,j) - ( dvdx - dudy )

      xc = ( x(i+1,j+1) + x(i+1,j) + x(i,j) ) / 3.0D+00
      yc = ( y(i+1,j+1) + y(i+1,j) + y(i,j) ) / 3.0D+00

      dudx = yc
      dudy = xc
      dvdx = 2.0D+00 * xc * cos ( xc * xc + yc * yc )
      dvdy = 2.0D+00 * yc * cos ( xc * xc + yc * yc )

      write ( *, '(2i4,2f7.4,3g14.6)' ) &
        i, j, xc, yc, vort(2,i,j), dvdx - dudy, vort(2,i,j) - ( dvdx - dudy )

    end do
    write ( *, '(a)' ) ' '
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 demonstrates ELEMENTS_EPS using Q4 elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: element_order = 4
  integer   ( kind = 4 ), parameter :: nelemx = 7
  integer   ( kind = 4 ), parameter :: nelemy = 5

  real      ( kind = 8 ) angle
  character ( len = 4 ) :: code = 'Q4'
  integer   ( kind = 4 ) col
  real      ( kind = 8 ) degrees_to_radians
  integer   ( kind = 4 ) element
  logical, allocatable, dimension ( : ) :: element_mask
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer   ( kind = 4 ) element_num
  character ( len = 80 ) :: file_name = 'fem2d_pack_prb_q4.eps'
  integer   ( kind = 4 ) node
  integer   ( kind = 4 ) node_num
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  real      ( kind = 8 ) r
  integer   ( kind = 4 ) row
  character ( len = 80 ) :: title = 'Grid of Q4 Elements'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  ELEMENTS_EPS creates an Encapsulated PostScript'
  write ( *, '(a)' ) '  file containing an image of a Q4 mesh.'
!
!  How many elements are there?
!
  call grid_q4_element_num ( nelemx, nelemy, element_num )

  allocate ( element_mask(element_num) )
  allocate ( element_node(element_order,element_num) )
!
!  How many nodes are there?
!
  call grid_q4_node_num ( nelemx, nelemy, node_num )

  allocate ( node_xy(2,node_num) )
!
!  Get the nodes that make up each element.
!
  call grid_q4_element ( nelemx, nelemy, element_node )
!
!  Generate the coordinates of the nodes.
!
  node = 0

  do row = 0, nelemy

    r = ( real ( nelemy - row, kind = 8 ) *  1.0D+00   &
        + real (        + row, kind = 8 ) *  3.0D+00 ) &
        / real ( nelemy,       kind = 8 )

    do col = 0, nelemx

      node = node + 1
      angle = ( real ( nelemx - col, kind = 8 ) * 135.0D+00   &
              + real (        + col, kind = 8 ) *  45.0D+00 ) &
              / real ( nelemx,       kind = 8 )

      angle = degrees_to_radians ( angle )

      node_xy(1,node) = r * cos ( angle )
      node_xy(2,node) = r * sin ( angle )

    end do

  end do

  do element = 1, element_num
    element_mask(element) = .true.
  end do

  call elements_eps ( file_name, node_num, node_xy, code, &
    element_order, element_num, element_mask, element_node, title )

  deallocate ( element_mask )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 demonstrates ELEMENTS_EPS using T3 elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: element_order = 3
  integer   ( kind = 4 ), parameter :: nelemx = 7
  integer   ( kind = 4 ), parameter :: nelemy = 5

  real      ( kind = 8 ) angle
  character ( len = 4 ) :: code = 'T3'
  integer   ( kind = 4 ) col
  real      ( kind = 8 ) degrees_to_radians
  integer   ( kind = 4 ) element
  logical, allocatable, dimension ( : ) :: element_mask
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer   ( kind = 4 ) element_num
  character ( len = 80 ) :: file_name = 'fem2d_pack_prb_t3.eps'
  integer   ( kind = 4 ) node
  integer   ( kind = 4 ) node_num
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  real      ( kind = 8 ) r
  integer   ( kind = 4 ) row
  character ( len = 80 ) :: title = 'Grid of T3 Elements'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  ELEMENTS_EPS creates an Encapsulated PostScript'
  write ( *, '(a)' ) '  file containing an image of a T3 mesh.'
!
!  How many elements are there?
!
  call grid_t3_element_num ( nelemx, nelemy, element_num )

  allocate ( element_mask(element_num) )
  allocate ( element_node(element_order,element_num) )
!
!  How many nodes are there?
!
  call grid_t3_node_num ( nelemx, nelemy, node_num )
  allocate ( node_xy(2,node_num) )
!
!  Get the nodes that make up each element.
!
  call grid_t3_element ( nelemx, nelemy, element_node )
!
!  Generate the coordinates of the nodes.
!
  node = 0

  do row = 0, nelemy

    r = ( real ( nelemy - row, kind = 8 ) *  1.0D+00   &
        + real (        + row, kind = 8 ) *  3.0D+00 ) &
        / real ( nelemy,       kind = 8 )

    do col = 0, nelemx

      node = node + 1
      angle = ( real ( nelemx - col, kind = 8 ) * 135.0D+00   &
              + real (        + col, kind = 8 ) *  45.0D+00 ) &
              / real ( nelemx,       kind = 8 )

      angle = degrees_to_radians ( angle )

      node_xy(1,node) = r * cos ( angle )
      node_xy(2,node) = r * sin ( angle )

    end do

  end do

  element_mask(1:element_num) = .true.

  call elements_eps ( file_name, node_num, node_xy, code, &
    element_order, element_num, element_mask, element_node, title )

  deallocate ( element_mask )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 demonstrates ELEMENTS_EPS using T4 elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: element_order = 4
  integer ( kind = 4 ), parameter :: nelemx = 7
  integer ( kind = 4 ), parameter :: nelemy = 5

  real    ( kind = 8 ) angle
  character ( len = 4 ) :: code = 'T4'
  integer ( kind = 4 ) col
  real    ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) element
  logical, allocatable, dimension ( : ) :: element_mask
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  character ( len = 80 ) :: file_name = 'fem2d_pack_prb_t4.eps'
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  real    ( kind = 8 ) r
  integer ( kind = 4 ) row
  character ( len = 80 ) :: title = 'Grid of T4 Elements'
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  ELEMENTS_EPS creates an Encapsulated PostScript'
  write ( *, '(a)' ) '  file containing an image of a T4 mesh.'
!
!  How many elements are there?
!
  call grid_t4_element_num ( nelemx, nelemy, element_num )

  allocate ( element_mask(element_num) )
  allocate ( element_node(element_order,element_num) )
!
!  How many nodes are there?
!
  call grid_t4_node_num ( nelemx, nelemy, node_num )

  allocate ( node_xy(2,node_num) )
!
!  Get the nodes that make up each element.
!
  call grid_t4_element ( nelemx, nelemy, element_node )
!
!  Generate the coordinates of the nodes.
!
  node = 0

  do row = 0, nelemy

    y = ( real ( 3*nelemy - row, kind = 8 ) *  0.0D+00   &
        + real (          + row, kind = 8 ) *  6.0D+00 ) &
        / real ( 3*nelemy,       kind = 8 )

    do col = 0, nelemx

      node = node + 1

      x = ( real ( 2*nelemx - col, kind = 8 ) * 0.0D+00   &
          + real (          + col, kind = 8 ) * 6.0D+00 ) &
          / real ( 2*nelemx,       kind = 8 )

      node_xy(1,node) = x
      node_xy(2,node) = y

    end do
!
!  Skip over the two rows of interior nodes.
!
    node = node + nelemx
    node = node + nelemx

  end do
!
!  The coordinates of interior nodes are the average of the vertices.
!
  do element = 1, element_num
    node = element_node(4,element)
    node_xy(1,node) = sum ( node_xy(1,element_node(1:3,element)) ) / 3.0D+00
    node_xy(2,node) = sum ( node_xy(2,element_node(1:3,element)) ) / 3.0D+00
  end do

  element_mask(1:element_num) = .true.

  call elements_eps ( file_name, node_num, node_xy, code, &
    element_order, element_num, element_mask, element_node, title )

  deallocate ( element_mask )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 demonstrates ELEMENTS_EPS using T6 elements.
!
!  Discussion:
!
!    We generate a big grid of T6 elements, but we only want to
!    look at the six elements shared by node 85.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: element_order = 6
  integer   ( kind = 4 ), parameter :: nelemx = 6
  integer   ( kind = 4 ), parameter :: nelemy = 6

  character ( len = 4 ) :: code = 'T6'
  integer   ( kind = 4 ) col
  integer   ( kind = 4 ) element
  logical, allocatable, dimension ( : ) :: element_mask
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer   ( kind = 4 ) element_num
  character ( len = 80 ) :: file_name = 'fem2d_pack_prb_t6.eps'
  integer   ( kind = 4 ) node
  integer   ( kind = 4 ) node_num
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer   ( kind = 4 ) row
  character ( len = 80 ) :: title = 'Grid of T6 Elements'
  real      ( kind = 8 ) x
  real      ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  ELEMENTS_EPS creates an Encapsulated PostScript'
  write ( *, '(a)' ) '  file containing an image of a T6 mesh.'
!
!  How many elements are there?
!
  call grid_t6_element_num ( nelemx, nelemy, element_num )

  allocate ( element_mask(element_num) )
  allocate ( element_node(element_order,element_num) )
!
!  How many nodes are there?
!
  call grid_t6_node_num ( nelemx, nelemy, node_num )

  allocate ( node_xy(2,node_num) )
!
!  Get the nodes that make up each element.
!
  call grid_t6_element ( nelemx, nelemy, element_node )
!
!  Generate the coordinates of the nodes.
!
  node = 0

  do row = 0, 2*nelemy

    y = ( real ( 2*nelemy - row, kind = 8 ) *  0.0D+00   &
        + real (          + row, kind = 8 ) *  6.0D+00 ) &
        / real ( 2*nelemy,       kind = 8 )

    do col = 0, 2*nelemx

      node = node + 1
      x = ( real ( 2*nelemx - col, kind = 8 ) * 0.0D+00   &
          + real (          + col, kind = 8 ) * 6.0D+00 ) &
          / real ( 2*nelemx,       kind = 8 )

      node_xy(1,node) = x
      node_xy(2,node) = y

    end do

  end do

  element_mask(1:element_num) = .false.

  element_mask(30) = .true.
  element_mask(31) = .true.
  element_mask(32) = .true.
  element_mask(41) = .true.
  element_mask(42) = .true.
  element_mask(43) = .true.

  call elements_eps ( file_name, node_num, node_xy, code, &
    element_order, element_num, element_mask, element_node, title )

  deallocate ( element_mask )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test105 ( )

!*****************************************************************************80
!
!! TEST105 tests GRID_NODES_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: num_x = 5
  integer ( kind = 4 ), parameter :: num_y = 3

  integer ( kind = 4 ), parameter :: node_num = num_x * num_y

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  real    ( kind = 8 ) node_xy(2,node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GRID_NODES_01 creates a regular grid in the unit square.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NUM_X =    ', num_x
  write ( *, '(a,i8)' ) '  NUM_Y =    ', num_y
  write ( *, '(a,i8)' ) '  NODE_NUM = ', node_num
  write ( *, '(a)' ) ' '

  call grid_nodes_01 ( num_x, num_y, node_xy )

  do node = 1, node_num
    write ( *, '(2x,i8,2x,f14.6,2x,f14.6)' ) node, node_xy(1:2,node)
  end do

  return
end
subroutine test107 ( )

!*****************************************************************************80
!
!! TEST107 demonstrates the GRID_T3 routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: element_order = 3
  integer   ( kind = 4 ), parameter :: nelemx = 8
  integer   ( kind = 4 ), parameter :: nelemy = 8

  character ( len = 4 ) :: code = 'T3'
  integer   ( kind = 4 ) col
  integer   ( kind = 4 ) element
  character ( len = 80 ) :: element_filename = 'rectangle_t3_triangles.txt'
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) node
  character ( len = 80 ) :: node_filename = 'rectangle_t3_nodes.txt'
  integer   ( kind = 4 ) node_num
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer   ( kind = 4 ) row
  real      ( kind = 8 ) x
  real      ( kind = 8 ) :: x_max = 8.0D+00
  real      ( kind = 8 ) :: x_min = 0.0D+00
  real      ( kind = 8 ) y
  real      ( kind = 8 ) :: y_max = 8.0D+00
  real      ( kind = 8 ) :: y_min = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST107'
  write ( *, '(a)' ) '  For a T3 mesh of a rectangle,'
  write ( *, '(a)' ) '  GRID_T3_ELEMENT_NUM counte elements;'
  write ( *, '(a)' ) '  GRID_T3_NODE_NUM counte nodes;'
  write ( *, '(a)' ) '  GRID_T3_ELEMENT creates elements;'
!
!  How many elements are there?
!
  call grid_t3_element_num ( nelemx, nelemy, element_num )

  allocate ( element_node(element_order,element_num) )
!
!  How many nodes are there?
!
  call grid_t3_node_num ( nelemx, nelemy, node_num )
  allocate ( node_xy(2,node_num) )
!
!  Get the nodes that make up each element.
!
  call grid_t3_element ( nelemx, nelemy, element_node )
!
!  Generate the coordinates of the nodes.
!
  node = 0

  do row = 0, nelemy

    y = ( real ( nelemy - row, kind = 8 ) *  y_min   &
        + real (        + row, kind = 8 ) *  y_max ) &
        / real ( nelemy,       kind = 8 )

    do col = 0, nelemx

      node = node + 1
      x = ( real ( nelemx - col, kind = 8 ) * x_min   &
          + real (        + col, kind = 8 ) *  x_max ) &
          / real ( nelemx,       kind = 8 )

      node_xy(1,node) = x
      node_xy(2,node) = y

    end do

  end do

  call r8mat_write ( node_filename, 2, node_num, node_xy )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote node file "' // trim ( node_filename ) // '".'

  call i4mat_write ( element_filename, 3, element_num, element_node )
  write ( *, '(a)' ) '  Wrote element file "' // trim ( element_filename ) // '".'

  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests GRID_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  GRID_TEST tests the grid routines.'

  call grid_test ( 'Q4' )

  call grid_test ( 'Q8' )

  call grid_test ( 'Q9' )

  call grid_test ( 'Q12' )

  call grid_test ( 'Q16' )

  call grid_test ( 'QL' )

  call grid_test ( 'T3' )

  call grid_test ( 'T4' )

  call grid_test ( 'T6' )

  call grid_test ( 'T10' )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests INTERP_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  INTERP_TEST tests the interpolating power'
  write ( *, '(a)' ) '  of the element.'

  call interp_test ( 'Q4' )

  call interp_test ( 'Q8' )

  call interp_test ( 'Q9' )

  call interp_test ( 'Q12' )

  call interp_test ( 'Q16' )

  call interp_test ( 'QL' )

  call interp_test ( 'T3' )

  call interp_test ( 'T4' )

  call interp_test ( 'T6' )

  call interp_test ( 'T10' )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests MAP_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  MAP_TEST tests the map routines.'

  call map_test ( 'Q4' )

  call map_test ( 'Q8' )

  call map_test ( 'Q9' )

  call map_test ( 'Q12' )

  call map_test ( 'Q16' )

  call map_test ( 'QL' )

  call map_test ( 'T3' )

  call map_test ( 'T4' )

  call map_test ( 'T6' )

  call map_test ( 'T10' )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 demonstrates MASS_MATRIX_T6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: element_num = 2
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ), parameter :: node_num = 9

  real    ( kind = 8 ) a(node_num,node_num)
  integer ( kind = 4 ), dimension ( element_order, element_num ) :: element_node = &
    reshape ( (/ &
    1, 3, 7, 2, 5, 4, &
    9, 7, 3, 8, 5, 6 /), (/ element_order, element_num /) )
  real    ( kind = 8 ), dimension ( 2, node_num ) :: node_xy = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.5D+00, &
    0.0D+00, 1.0D+00, &
    0.5D+00, 0.0D+00, &
    0.5D+00, 0.5D+00, &
    0.5D+00, 1.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 0.5D+00, &
    1.0D+00, 1.0D+00 /), (/ 2, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  MASS_MATRIX_T6 computes the mass matrix for'
  write ( *, '(a)' ) '  a finite element system using T6 elements'
  write ( *, '(a)' ) '  (quadratic triangles).'

  call mass_matrix_t6 ( node_num, element_num, element_node, node_xy, a )

  call r8mat_print ( node_num, node_num, a, '  The T6 mass matrix:' )

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests PHYSICAL_TO_REFERENCE_T3 and REFERENCE_TO_PHYSICAL_T3.
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

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phy(2,n)
  real    ( kind = 8 ) r8_uniform_01
  real    ( kind = 8 ) ref(2,n)
  real    ( kind = 8 ) ref2(2,n)
  integer ( kind = 4 ) seed
  real    ( kind = 8 ), dimension(2,3) :: t = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, &
    2.0D+00, 5.0D+00 /), (/ 2, 3 /) )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For an order 3 triangle,'
  write ( *, '(a)' ) '  PHYSICAL_TO_REFERENCE_T3 maps a physical point to'
  write ( *, '(a)' ) '    a reference point.'
  write ( *, '(a)' ) '  REFERENCE_TO_PHYSICAL_T3 maps a reference point to'
  write ( *, '(a)' ) '    a physical point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (   XSI    ETA ) ==> ( X     Y  )  ==> ( XSI2   ETA2 )'
  write ( *, '(a)' ) ' '

  do j = 1, n

    do i = 1, 2
      ref(i,j) = r8_uniform_01 ( seed )
    end do

    if ( 1.0D+00 < sum ( ref(1:2,j) ) ) then
      ref(1:2,j) = 1.0D+00 - ref(1:2,j)
    end if

  end do

  call reference_to_physical_t3 ( t, n, ref, phy )
  call physical_to_reference_t3 ( t, n, phy, ref2 )

  do j = 1, n

    write ( *, '(2x,2f8.4,2x,2f8.4,2x,2f8.4)' ) &
      ref(1:2,j), phy(1:2,j), ref2(1:2,j)
  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests REFERENCE_TO_PHYSICAL_T6.
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

  integer ( kind = 4 ), parameter :: n = 16

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phy(2,n)
  real    ( kind = 8 ), dimension (2,n) :: ref = reshape ( (/ &
    0.00D+00, 0.00D+00, &
    1.00D+00, 0.00D+00, &
    0.00D+00, 1.00D+00, &
    0.50D+00, 0.00D+00, &
    0.50D+00, 0.50D+00, &
    0.00D+00, 0.50D+00, &
    0.25D+00, 0.75D+00, &
    0.75D+00, 0.25D+00, &
    0.40D+00, 0.10D+00, &
    0.30D+00, 0.20D+00, &
    0.20D+00, 0.30D+00, &
    0.10D+00, 0.40D+00, &
    0.10D+00, 0.10D+00, &
    0.20D+00, 0.20D+00, &
    0.30D+00, 0.30D+00, &
    0.40D+00, 0.40D+00 /), (/ 2, n /) )
  real    ( kind = 8 ), dimension(2,6) :: t = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    0.0D+00, 4.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 2.0D+00 /), (/ 2, 6 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  For an order 6 triangle,'
  write ( *, '(a)' ) '  REFERENCE_TO_PHYSICAL_T6 maps a reference point to'
  write ( *, '(a)' ) '    a physical point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      XSI     ETA  ==>  X       Y'
  write ( *, '(a)' ) ' '

  call reference_to_physical_t6 ( t, n, ref, phy )

  do j = 1, n
    write ( *, '(2x,2f8.4,2x,2f8.4)' ) ref(1:2,j), phy(1:2,j)
  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 demonstrates S_L2NORM.
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

  integer ( kind = 4 ), parameter :: element_num = 40
  integer ( kind = 4 ), parameter :: psi_num = element_num + 1
  integer ( kind = 4 ), parameter :: quad_num = 5

  integer ( kind = 4 ) element
  real    ( kind = 8 ) element_area(element_num)
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) i
  real    ( kind = 8 ) l2norm
  real    ( kind = 8 ) node_x(psi_num)
  integer ( kind = 4 ) psi
  real    ( kind = 8 ) psi_quad(psi_num,element_num,quad_num)
  integer ( kind = 4 ) quad
  real    ( kind = 8 ) quad_weight(quad_num)
  real    ( kind = 8 ) quad_x(quad_num)
  real    ( kind = 8 ) s_coef(psi_num)
  real    ( kind = 8 ) x
  real    ( kind = 8 ), parameter :: x_max = 1.0D+00
  real    ( kind = 8 ), parameter :: x_min = 0.0D+00
  real    ( kind = 8 ) xl
  real    ( kind = 8 ) xr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  S_L2NORM computes the L2 norm of a scalar function'
  write ( *, '(a)' ) '  S(X) over a region (of any dimension), assuming:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * that there is a set of finite element basis'
  write ( *, '(a)' ) '    functions PSI;'
  write ( *, '(a)' ) '  * that the region is decomposed into a number of'
  write ( *, '(a)' ) '    elements whose areas are known;'
  write ( *, '(a)' ) '  * that the integral is to be computed by a quadrature'
  write ( *, '(a)' ) '    rule applied in the same way to each element;'
  write ( *, '(a)' ) '  * that the value of the basis functions is given'
  write ( *, '(a)' ) '    at every quadrature node in every element;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our example will have one spatial dimension.'
  write ( *, '(a,g14.6)' ) '  XMIN = ', x_min
  write ( *, '(a,g14.6)' ) '  XMAX = ', x_max
  write ( *, '(a,i8)' ) '  The number of intervals will be    ', element_num
  write ( *, '(a,i8)' ) '  The number of basis functions is   ', psi_num
  write ( *, '(a,i8)' ) '  The number of quadrature points is ', quad_num
!
!  Set the nodes.
!
  do i = 1, psi_num
    node_x(i) = ( real ( psi_num - i,     kind = 8 ) * x_min   &
                + real (           i - 1, kind = 8 ) * x_max ) &
                / real ( psi_num     - 1, kind = 8 )
  end do
!
!  Set the element "areas".
!
  do i = 1, element_num
    element_area(i) = node_x(i+1) - node_x(i)
  end do
!
!  Set the quadrature weights and abscissas.
!
  quad_weight(1) = 0.236926885056189087514264040720D+00
  quad_weight(2) = 0.478628670499366468041291514836D+00
  quad_weight(3) = 0.568888888888888888888888888889D+00
  quad_weight(4) = 0.478628670499366468041291514836D+00
  quad_weight(5) = 0.236926885056189087514264040720D+00

  quad_weight(1:5) = quad_weight(1:5) / 2.0D+00

  quad_x(1) = - 0.906179845938663992797626878299D+00
  quad_x(2) = - 0.538469310105683091036314420700D+00
  quad_x(3) =   0.0D+00
  quad_x(4) =   0.538469310105683091036314420700D+00
  quad_x(5) =   0.906179845938663992797626878299D+00
!
!  Set the finite element coefficients of S.  In our formulation,
!  these finite element coefficients are simply the function values
!  at the nodes.
!
  do i = 1, psi_num
    s_coef(i) = sin ( node_x(i) )
  end do
!
!  For each basis function I,
!    for each element J,
!      for each quadrature point K,
!
!  ...evaluate the basis function in the element at the quadrature point.
!
  psi_quad(1:psi_num,1:element_num,1:quad_num) = 0.0D+00

  do psi = 1, psi_num

    do element = 1, element_num

      if ( element < psi - 1 ) then

        cycle

      else if ( element == psi - 1 ) then

        xl = node_x(element)
        xr = node_x(element+1)

        do quad = 1, quad_num

          x = ( ( 1.0D+00 - quad_x(quad) ) * xl   &
              + ( 1.0D+00 + quad_x(quad) ) * xr ) &
              /   2.0D+00

          psi_quad(psi,element,quad) = ( x - xl ) / ( xr - xl )

        end do

      else if ( element == psi ) then

        xl = node_x(element)
        xr = node_x(element+1)

        do quad = 1, quad_num

          x = ( ( 1.0D+00 - quad_x(quad) ) * xl &
              + ( 1.0D+00 + quad_x(quad) ) * xr ) &
              /   2.0D+00

          psi_quad(psi,element,quad) = ( xr - x ) / ( xr - xl )

        end do

      else if ( psi < element ) then
        cycle
      end if

    end do

  end do
!
!  Ask S_L2NORM to compute the integral of L2 norm of S.
!
  call s_l2norm ( psi_num, element_num, quad_num, element_area, &
    quad_weight, psi_quad, s_coef, l2norm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The computed L2 norm is ', l2norm
!
!  The integral of (sin(x))**2 is x/2 - sin(2x)/4
!
  exact = sqrt ( &
      ( 0.5D+00 * x_max - 0.25D+00 * sin ( 2.0D+00 * x_max ) ) &
    - ( 0.5D+00 * x_min - 0.25D+00 * sin ( 2.0D+00 * x_min ) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The exact value is      ', exact

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests SHAPE_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  SHAPE_TEST tests the shape routines.'

  call shape_test ( 'Q4' )

  call shape_test ( 'Q8' )

  call shape_test ( 'Q9' )

  call shape_test ( 'Q12' )

  call shape_test ( 'Q16' )

  call shape_test ( 'QL' )

  call shape_test ( 'T3' )

  call shape_test ( 'T4' )

  call shape_test ( 'T6' )

  call shape_test ( 'T10' )

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests SPHERE_GRID_Q4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) element
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 4
  integer ( kind = 4 ), parameter :: nelemx = 8
  integer ( kind = 4 ), parameter :: nelemy = 8
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  SPHERE_GRID_Q4_ELEMENT sets up a grid of '
  write ( *, '(a)' ) '    Q4 quadrilaterals on a sphere.'
  write ( *, '(a)' ) '  SPHERE_GRID_Q4_ELEMENT_NUM returns the number'
  write ( *, '(a)' ) '    of elements in the grid'
  write ( *, '(a)' ) '  SPHERE_GRID_Q4_NODE_NUM returns the number'
  write ( *, '(a)' ) '    of nodes in the grid.'
  write ( *, '(a)' ) '  SPHERE_GRID_Q4_NODE_XYZ returns the coordinates'
  write ( *, '(a)' ) '    of nodes in the grid.'

  call sphere_grid_q4_element_num ( nelemx, nelemy, element_num )
  call sphere_grid_q4_node_num ( nelemx, nelemy, node_num )

  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( node_xyz(1:3,1:node_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Expected number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Expected number of elements = ', element_num

  call sphere_grid_q4_element ( nelemx, nelemy, element_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The elements and their nodes:'
  write ( *, '(a)' ) ' '

  do element = 1, element_num
    write ( *, '(i4,2x,4i4)' ) element, element_node(1:element_order,element)
  end do

  call sphere_grid_q4_node_xyz ( nelemx, nelemy, node_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The node coordinates:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i4,2x,3g14.6)' ) node, node_xyz(1:3,node)
  end do
!
!  Write the elements and nodes to files.
!
  call r8mat_write ( 'sphere_q4_nodes.txt', 3, node_num, node_xyz )

  call i4mat_write ( 'sphere_q4_elements.txt', element_order, element_num, &
    element_node )

  deallocate ( element_node )
  deallocate ( node_xyz )

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests SPHERE_GRID_Q9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) element
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 9
  integer ( kind = 4 ), parameter :: nelemx = 3
  integer ( kind = 4 ), parameter :: nelemy = 4
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  SPHERE_GRID_Q9_ELEMENT sets up a grid of '
  write ( *, '(a)' ) '    Q9 quadrilaterals on a sphere.'
  write ( *, '(a)' ) '  SPHERE_GRID_Q9_ELEMENT_NUM returns the number'
  write ( *, '(a)' ) '    of elements in the grid'
  write ( *, '(a)' ) '  SPHERE_GRID_Q9_NODE_NUM returns the number'
  write ( *, '(a)' ) '    of nodes in the grid.'
  write ( *, '(a)' ) '  SPHERE_GRID_Q9_NODE_XYZ returns the coordinates'
  write ( *, '(a)' ) '    of nodes in the grid.'

  call sphere_grid_q9_element_num ( nelemx, nelemy, element_num )
  call sphere_grid_q9_node_num ( nelemx, nelemy, node_num )

  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( node_xyz(1:3,1:node_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Expected number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Expected number of elements = ', element_num

  call sphere_grid_q9_element ( nelemx, nelemy, element_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The elements and their nodes:'
  write ( *, '(a)' ) ' '

  do element = 1, element_num
    write ( *, '(i4,2x,9i4)' ) element, element_node(1:element_order,element)
  end do

  call sphere_grid_q9_node_xyz ( nelemx, nelemy, node_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The node coordinates:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i4,2x,3g14.6)' ) node, node_xyz(1:3,node)
  end do
!
!  Write the elements and nodes to files.
!
  call r8mat_write ( 'sphere_q9_nodes.txt', 3, node_num, node_xyz )

  call i4mat_write ( 'sphere_q9_elements.txt', element_order, element_num, &
    element_node )

  deallocate ( element_node )
  deallocate ( node_xyz )

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests SPHERE_GRID_Q16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) element
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 16
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: nelemx = 2
  integer ( kind = 4 ), parameter :: nelemy = 2
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  SPHERE_GRID_Q16_ELEMENT sets up a grid of '
  write ( *, '(a)' ) '    Q16 quadrilaterals on a sphere.'
  write ( *, '(a)' ) '  SPHERE_GRID_Q16_ELEMENT_NUM returns the number'
  write ( *, '(a)' ) '    of elements in the grid'
  write ( *, '(a)' ) '  SPHERE_GRID_Q16_NODE_NUM returns the number'
  write ( *, '(a)' ) '    of nodes in the grid.'
  write ( *, '(a)' ) '  SPHERE_GRID_Q16_NODE_XYZ returns the coordinates'
  write ( *, '(a)' ) '    of nodes in the grid.'

  call sphere_grid_q16_element_num ( nelemx, nelemy, element_num )
  call sphere_grid_q16_node_num ( nelemx, nelemy, node_num )

  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( node_xyz(1:3,1:node_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Expected number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Expected number of elements = ', element_num

  call sphere_grid_q16_element ( nelemx, nelemy, element_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The elements and their nodes, listed in a way'
  write ( *, '(a)' ) '  that suggests their geometry:'
  write ( *, '(a)' ) ' '

  element = element_num

  do j = 1, nelemy
    do i = 1, nelemx
      write ( *, '(a)' ) ' '
      write ( *, '(i4,2x,4i4)' ) element, element_node(13:16,element)
      write ( *, '(4x,2x,4i4)' )          element_node(9:12,element)
      write ( *, '(4x,2x,4i4)' )          element_node(5:8,element)
      write ( *, '(4x,2x,4i4)' )          element_node(1:4,element)
      element = element - 1
    end do
  end do

  call sphere_grid_q16_node_xyz ( nelemx, nelemy, node_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The node coordinates:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i4,2x,3g14.6)' ) node, node_xyz(1:3,node)
  end do
!
!  Write the elements and nodes to files.
!
  call r8mat_write ( 'sphere_q16_nodes.txt', 3, node_num, node_xyz )

  call i4mat_write ( 'sphere_q16_elements.txt', element_order, element_num, &
    element_node )

  deallocate ( element_node )
  deallocate ( node_xyz )

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests SPHERE_GRID_T3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) element
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), parameter :: nelemx = 8
  integer ( kind = 4 ), parameter :: nelemy = 8
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  SPHERE_GRID_T3_ELEMENT sets up a grid of T3 triangles'
  write ( *, '(a)' ) '    on a sphere.'
  write ( *, '(a)' ) '  SPHERE_GRID_T3_ELEMENT_NUM returns the number'
  write ( *, '(a)' ) '    of elements in the grid'
  write ( *, '(a)' ) '  SPHERE_GRID_T3_NODE_NUM returns the number'
  write ( *, '(a)' ) '    of nodes in the grid.'
  write ( *, '(a)' ) '  SPHERE_GRID_T3_NODE_XYZ returns the coordinates'
  write ( *, '(a)' ) '    of nodes in the grid.'

  call sphere_grid_t3_element_num ( nelemx, nelemy, element_num )
  call sphere_grid_t3_node_num ( nelemx, nelemy, node_num )

  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( node_xyz(1:3,1:node_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Expected number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Expected number of elements = ', element_num

  call sphere_grid_t3_element ( nelemx, nelemy, element_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The elements and their nodes:'
  write ( *, '(a)' ) ' '

  do element = 1, element_num
    write ( *, '(i4,2x,3i4)' ) element, element_node(1:element_order,element)
  end do

  call sphere_grid_t3_node_xyz ( nelemx, nelemy, node_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The node coordinates:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i4,2x,3g14.6)' ) node, node_xyz(1:3,node)
  end do
!
!  Write the elements and nodes to files.
!
  call r8mat_write ( 'sphere_t3_nodes.txt', 3, node_num, node_xyz )

  call i4mat_write ( 'sphere_t3_elements.txt', element_order, element_num, &
    element_node )

  deallocate ( element_node )
  deallocate ( node_xyz )

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests SPHERE_GRID_T6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) element
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ), parameter :: nelemx = 3
  integer ( kind = 4 ), parameter :: nelemy = 4
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  SPHERE_GRID_T6_ELEMENT sets up a grid of T6 triangles'
  write ( *, '(a)' ) '    on a sphere.'
  write ( *, '(a)' ) '  SPHERE_GRID_T6_ELEMENT_NUM returns the number'
  write ( *, '(a)' ) '    of elements in the grid'
  write ( *, '(a)' ) '  SPHERE_GRID_T6_NODE_NUM returns the number'
  write ( *, '(a)' ) '    of nodes in the grid.'
  write ( *, '(a)' ) '  SPHERE_GRID_T6_NODE_XYZ returns the coordinates'
  write ( *, '(a)' ) '    of nodes in the grid.'

  call sphere_grid_t6_element_num ( nelemx, nelemy, element_num )
  call sphere_grid_t6_node_num ( nelemx, nelemy, node_num )

  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( node_xyz(1:3,1:node_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Expected number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Expected number of elements = ', element_num

  call sphere_grid_t6_element ( nelemx, nelemy, element_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The elements and their nodes:'
  write ( *, '(a)' ) ' '

  do element = 1, element_num
    write ( *, '(i4,2x,6i4)' ) element, element_node(1:element_order,element)
  end do

  call sphere_grid_t6_node_xyz ( nelemx, nelemy, node_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The node coordinates:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i4,2x,3g14.6)' ) node, node_xyz(1:3,node)
  end do
!
!  Write the elements and nodes to files.
!
  call r8mat_write ( 'sphere_t6_nodes.txt', 3, node_num, node_xyz )

  call i4mat_write ( 'sphere_t6_elements.txt', element_order, element_num, &
    element_node )

  deallocate ( element_node )
  deallocate ( node_xyz )

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests TRIANGLE_UNIT_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 64

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real    ( kind = 8 ) coef
  real    ( kind = 8 ) err
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real    ( kind = 8 ) quad
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), parameter :: rule_max = 20
  integer ( kind = 4 ) triangle_unit_size
  real    ( kind = 8 ) value
  real    ( kind = 8 ) weight(order_max)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xtab(order_max)
  real    ( kind = 8 ) y
  real    ( kind = 8 ) ytab(order_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  TRIANGLE_UNIT_SET sets up a quadrature '
  write ( *, '(a)' ) '    in the unit triangle,'
  write ( *, '(a)' ) ' '

  do a = 0, 10

    do b = 0, 10 - a

      coef = real ( a + b + 2, kind = 8 ) * real ( a + b + 1, kind = 8 )
      do i = 1, b
        coef = coef * real ( a + i, kind = 8 ) / real ( i, kind = 8 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' ) '  A = ', a, '  B = ', b
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rule       QUAD           ERROR'
      write ( *, '(a)' ) ' '

      do rule = 1, rule_max

        order = triangle_unit_size ( rule )

        call triangle_unit_set ( rule, xtab, ytab, weight )

        quad = 0.0D+00

        do i = 1, order

          x = xtab(i)
          y = ytab(i)

          if ( a == 0 .and. b == 0 ) then
            value = coef
          else if ( a == 0 .and. b /= 0 ) then
            value = coef              * ytab(i)**b
          else if ( a /= 0 .and. b == 0 ) then
            value = coef * xtab(i)**a
          else if ( a /= 0 .and. b /= 0 ) then
            value = coef * xtab(i)**a * ytab(i)**b
          end if

          quad = quad + 0.5D+00 * weight(i) * value

        end do

        exact = 1.0D+00
        err = abs ( exact - quad )

        write ( *, '(2x,i4,2x,g14.6,2x,f14.8)' ) rule, quad, err

      end do

    end do

  end do

  return
end
