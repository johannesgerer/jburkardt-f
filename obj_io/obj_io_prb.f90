program main

!*****************************************************************************80
!
!! MAIN is the main program for OBJ_IO_PRB.
!
!  Discussion:
!
!    OBJ_IO_PRB runs the tests of the OBJ_IO routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OBJ_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the OBJ_IO library.'

  call test01 ( 'cube.obj' )
  call test02 ( 'cube.obj' )
  call test03  ( 'cube_normals.obj' )
  call test04 ( 'cube_no_normals.obj' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OBJ_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( input_file_name )

!*****************************************************************************80
!
!! TEST01 tests OBJ_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) face_num
  character ( len = * ) input_file_name
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) normal_num
  integer ( kind = 4 ) order_max

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  OBJ_SIZE determines the size of various objects'
  write ( *, '(a)' ) '  in an OBJ file.'

  call obj_size ( input_file_name, node_num, face_num, normal_num, &
    order_max )

  call obj_size_print ( input_file_name, node_num, face_num, &
    normal_num, order_max )

  return
end
subroutine test02 ( input_file_name )

!*****************************************************************************80
!
!! TEST02 tests OBJ_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension(:,:) :: face_node
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension(:) :: face_order
  integer ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension(:,:) :: node_xyz
  real ( kind = 8 ), allocatable, dimension(:,:) :: normal_vector
  integer ( kind = 4 ) normal_num
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ), allocatable, dimension(:,:) :: vertex_normal

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  OBJ_READ reads an object in an OBJ file.'

  call obj_size ( input_file_name, node_num, face_num, normal_num, &
    order_max )

  allocate ( face_node(order_max,face_num) )
  allocate ( face_order(face_num) )
  allocate ( node_xyz(3,node_num) )
  allocate ( normal_vector(3,normal_num) )
  allocate ( vertex_normal(order_max,face_num) )

  call obj_read ( input_file_name, node_num, face_num, normal_num, &
    order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal )

  call obj_face_node_print ( face_num, order_max, face_order, face_node )
  call obj_normal_vector_print ( normal_num, normal_vector )
  call obj_node_xyz_print ( node_num, node_xyz )

  deallocate ( face_node )
  deallocate ( face_order )
  deallocate ( node_xyz )
  deallocate ( normal_vector )
  deallocate ( vertex_normal )

  return
end
subroutine test03 ( output_file_name )

!*****************************************************************************80
!
!! TEST03 tests OBJ_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: face_num = 12
  integer ( kind = 4 ), parameter :: node_num = 8
  integer ( kind = 4 ), parameter :: normal_num = 6
  integer ( kind = 4 ), parameter :: order_max = 3

  integer ( kind = 4 ), dimension ( order_max, face_num ) :: face_node = reshape ( (/ &
    1, 3, 2, &
    2, 3, 4, &
    1, 6, 5, &
    1, 2, 6, &
    3, 7, 4, &
    4, 7, 8, &
    5, 6, 8, &
    5, 8, 7, &
    1, 5, 7, &
    1, 7, 3, &
    2, 4, 6, &
    6, 4, 8 /), (/ order_max, face_num /) )
  integer ( kind = 4 ), dimension ( face_num ) :: face_order = (/ &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
  real ( kind = 8 ), dimension ( 3, normal_num ) :: normal_vector = &
    reshape ( (/ &
    0.0D+00,  0.0D+00,  1.0D+00, &
    0.0D+00,  0.0D+00, -1.0D+00, &
    0.0D+00,  1.0D+00,  0.0D+00, &
    0.0D+00, -1.0D+00,  0.0D+00, &
    1.0D+00,  0.0D+00,  0.0D+00, &
   -1.0D+00,  0.0D+00,  0.0D+00 /), (/ 3, normal_num /) )
  character ( len = * ) :: output_file_name
  real ( kind = 8 ), dimension ( 3, node_num ) :: node_xyz = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00 /), (/ 3, node_num /) )
  integer ( kind = 4 ), dimension ( order_max, face_num ) :: vertex_normal = &
    reshape ( (/ &
    2, 2, 2, &
    2, 2, 2, &
    4, 4, 4, &
    4, 4, 4, &
    3, 3, 3, &
    3, 3, 3, &
    1, 1, 1, &
    1, 1, 1, &
    6, 6, 6, &
    6, 6, 6, &
    5, 5, 5, &
    5, 5, 5 /), (/ order_max, face_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  OBJ_WRITE writes an ASCII OBJ file.'

  call obj_write ( output_file_name, node_num, face_num, normal_num, &
    order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Graphics data was written to the OBJ file "' // &
    trim ( output_file_name ) // '":'

  return
end
subroutine test04 ( output_file_name )

!*****************************************************************************80
!
!! TEST04 tests OBJ_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: face_num = 12
  integer ( kind = 4 ), parameter :: node_num = 8
  integer ( kind = 4 ), parameter :: normal_num = 0
  integer ( kind = 4 ), parameter :: order_max = 3

  integer ( kind = 4 ), dimension ( order_max, face_num ) :: face_node = reshape ( (/ &
    1, 3, 2, &
    2, 3, 4, &
    1, 6, 5, &
    1, 2, 6, &
    3, 7, 4, &
    4, 7, 8, &
    5, 6, 8, &
    5, 8, 7, &
    1, 5, 7, &
    1, 7, 3, &
    2, 4, 6, &
    6, 4, 8 /), (/ order_max, face_num /) )
  integer ( kind = 4 ), dimension ( face_num ) :: face_order = (/ &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
  character ( len = * ) :: output_file_name
  real ( kind = 8 ), dimension ( 3, node_num ) :: node_xyz = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00 /), (/ 3, node_num /) )
  real ( kind = 8 ) normal_vector(1,1)
  integer ( kind = 4 ) vertex_normal(1,1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  OBJ_WRITE writes an ASCII OBJ file.'
  write ( *, '(a)' ) '  Here, we do NOT supply any normal vectors.'

  call obj_write ( output_file_name, node_num, face_num, normal_num, &
    order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Graphics data was written to the OBJ file "' // &
    trim ( output_file_name ) // '":'

  return
end
