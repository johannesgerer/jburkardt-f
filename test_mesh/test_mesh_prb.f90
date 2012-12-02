program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_MESH_PRB.
!
!  Discussion:
!
!    TEST_MESH_PRB runs the TEST_MESH tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MESH_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_MESH library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MESH_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests MESH00_ELEMENT_NUM, MESH00_NODE_NUM, MESH00_NUM, MESH00_NAME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mesh_num
  integer ( kind = 4 ) node_num
  character ( len = 80 ) name

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  MESH00_NUM reports the number of meshes.'
  write ( *, '(a)' ) '  MESH00_NAME gives the name of any mesh.'
  write ( *, '(a)' ) '  MESH00_ELEMENT_NUM gives the number of elements.'
  write ( *, '(a)' ) '  MESH00_NODE_NUM gives the number of nodes.'

  call mesh00_num ( mesh_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of meshes available = ', mesh_num

  do i = 1, mesh_num

    call mesh00_name ( i, name )
    call mesh00_element_num ( i, element_num )
    call mesh00_node_num ( i, node_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Mesh number:        ', i
    write ( *, '(a, a)' ) '  Name:              "', trim ( name ) // '"'
    write ( *, '(a,i6)' ) '  Number of elements: ', element_num
    write ( *, '(a,i6)' ) '  Number of nodes:    ', node_num

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests MESH00_ELEMENT_NODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mesh_num
  character ( len = 80 ) name
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) number
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ), allocatable, dimension ( : ) :: y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  MESH00_NAME gets the name of a mesh.'
  write ( *, '(a)' ) '  MESH00_NODE_NUM gets the number of nodes.'
  write ( *, '(a)' ) '  MESH00_NODE_XY gets the node coordinates.'
  write ( *, '(a)' ) '  MESH00_ELEMENT_NUM gets the number of elements.'
  write ( *, '(a)' ) '  MESH00_ELEMENT_NODE gets the element->node data.'

  call mesh00_num ( mesh_num )

  do number = 1, mesh_num

    call mesh00_name ( number, name )

    call mesh00_node_num ( number, node_num )

    call mesh00_element_num ( number, element_num )

    allocate ( x(1:node_num) )
    allocate ( y(1:node_num) )

    call mesh00_node_xy ( number, node_num, x, y )

    allocate ( element_node(1:3,1:element_num) )

    call mesh00_element_node ( number, element_num, element_node )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Mesh name: ' // trim ( name )
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of nodes = ', node_num

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Node coordinates:'
    write ( *, '(a)' ) ' '
    do i = 1, node_num
      write ( *, '(i4,2f14.6)' ) i, x(i), y(i)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of elements = ', element_num

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Element nodes:'
    write ( *, '(a)' ) ' '
    do j = 1, element_num
      write ( *, '(i4,3i6)' ) j, element_node(1:3,j)
    end do

    deallocate ( element_node )
    deallocate ( x )
    deallocate ( y )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests MESH00_NODE_EPS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) file_name
  integer ( kind = 4 ) mesh_num
  integer ( kind = 4 ) number

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  MESH00_NODE_EPS makes an EPS file containing'
  write ( *, '(a)' ) '  an image of the nodes of a mesh.'

  call mesh00_num ( mesh_num )

  file_name = 'mesh00_nodes.eps'

  write ( *, '(a)' ) ' '

  do number = 1, mesh_num

    call file_name_inc ( file_name )

    call mesh00_node_eps ( number, file_name )

    write ( *, '(a)' ) '  Created file ' // trim ( file_name )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests MESH00_ELEMENT3_EPS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) file_name
  integer ( kind = 4 ) mesh_num
  integer ( kind = 4 ) number

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  MESH00_ELEMENT3_EPS makes an EPS file containing'
  write ( *, '(a)' ) '  an image of the elements of a mesh.'

  call mesh00_num ( mesh_num )

  file_name = 'mesh00_elements.eps'

  write ( *, '(a)' ) ' '

  do number = 1, mesh_num

    call file_name_inc ( file_name )

    call mesh00_element3_eps ( number, file_name )

    write ( *, '(a)' ) '  Created file ' // trim ( file_name )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests MESH00_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) file_name
  integer ( kind = 4 ) mesh_num
  integer ( kind = 4 ) number

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  MESH00_POLY makes a POLY file containing the mesh'
  write ( *, '(a)' ) '  data, for input to TRIANGLE.'

  call mesh00_num ( mesh_num )

  file_name = 'mesh00.poly'

  write ( *, '(a)' ) ' '

  do number = 1, mesh_num

    call file_name_inc ( file_name )

    call mesh00_poly ( number, file_name )

    write ( *, '(a)' ) '  Created file ' // trim ( file_name )

  end do

  return
end
