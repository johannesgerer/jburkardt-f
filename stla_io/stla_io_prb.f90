program main

!*****************************************************************************80
!
!! MAIN is the main program for STLA_IO_PRB.
!
!  Discussion:
!
!    STLA_IO_PRB runs the tests of the STLA_IO routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STLA_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the STLA_IO library.'

  call test01 ( 'cube.stla' )
  call test02 ( 'cube.stla' )
  call test03 ( 'cube.stla' )
  call test04 ( 'cube_new.stla' )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STLA_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( input_file_name )

!*****************************************************************************80
!
!! TEST01 tests STLA_CHECK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = * ) input_file_name
  logical stla_check

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  STLA_CHECK makes some simple checks on a file.'

  write ( *, '(a)' ) ' '
  if ( stla_check ( input_file_name ) ) then
    write ( *, '(a)' ) '  The file "' // trim ( input_file_name ) // &
      '" seems to be a legal ASCII STL file.'
  else
    write ( *, '(a)' ) '  The file "' // trim ( input_file_name ) // &
      '" does NOT seem to be a legal ASCII STL file.'
  end if

  return
end
subroutine test02 ( input_file_name )

!*****************************************************************************80
!
!! TEST02 tests STLA_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) face_num
  character ( len = * ) input_file_name
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) solid_num
  integer ( kind = 4 ) text_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  STLA_SIZE determines the size of various objects'
  write ( *, '(a)' ) '  in an ASCII STL file.'

  call stla_size ( input_file_name, solid_num, node_num, face_num, text_num )

  call stla_size_print ( input_file_name, solid_num, node_num, face_num, &
    text_num )

  return
end
subroutine test03 ( input_file_name )

!*****************************************************************************80
!
!! TEST03 tests STLA_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension(:,:) :: face_node
  real ( kind = 8 ), allocatable, dimension(:,:) :: face_normal
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension(:,:) :: node_xyz
  integer ( kind = 4 ) solid_num
  integer ( kind = 4 ) text_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  STLA_READ reads an object in an ASCII STL file.'

  call stla_size ( input_file_name, solid_num, node_num, face_num, text_num )

  allocate ( face_node(3,face_num) )
  allocate ( face_normal(3,face_num) )
  allocate ( node_xyz(3,node_num) )

  call stla_read ( input_file_name, node_num, face_num, node_xyz, &
    face_node, face_normal, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  STLA_READ returned IERROR = ', ierror
    return
  end if

  call stla_size_print ( input_file_name, solid_num, node_num, face_num, &
    text_num )

  call stla_face_node_print ( face_num, face_node )
  call stla_face_normal_print ( face_num, face_normal )
  call stla_node_xyz_print ( node_num, node_xyz )

  deallocate ( face_node )
  deallocate ( face_normal )
  deallocate ( node_xyz )

  return
end
subroutine test04 ( output_file_name )

!*****************************************************************************80
!
!! TEST04 tests STLA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: face_num = 12
  integer ( kind = 4 ), parameter :: node_num = 8

  integer ( kind = 4 ), dimension ( 3, face_num ) :: face_node = reshape ( (/ &
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
    6, 4, 8 /), (/ 3, face_num /) )
  real ( kind = 8 ), dimension ( 3, face_num ) :: face_normal = reshape ( (/ &
    0.0D+00,  0.0D+00, -1.0D+00, &
    0.0D+00,  0.0D+00, -1.0D+00, &
    0.0D+00, -1.0D+00,  0.0D+00, &
    0.0D+00, -1.0D+00,  0.0D+00, &
    0.0D+00, +1.0D+00,  0.0D+00, &
    0.0D+00, +1.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00, +1.0D+00, &
    0.0D+00,  0.0D+00, +1.0D+00, &
   -1.0D+00,  0.0D+00,  0.0D+00, &
   -1.0D+00,  0.0D+00,  0.0D+00, &
   +1.0D+00,  0.0D+00,  0.0D+00, &
   +1.0D+00,  0.0D+00,  0.0D+00 /), (/ 3, face_num /) )
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  STLA_WRITE writes an ASCII STL file.'

  call stla_write ( output_file_name, node_num, face_num, node_xyz, &
    face_node, face_normal )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Graphics data was written to the STLA file "' // &
    trim ( output_file_name ) // '":'

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests STLA_FACE_NORMAL_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: face_num = 12
  integer ( kind = 4 ), parameter :: node_num = 8

  real ( kind = 8 ) dot_max
  integer ( kind = 4 ) face
  integer ( kind = 4 ), dimension ( 3, face_num ) :: face_node = reshape ( (/ &
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
    6, 4, 8 /), (/ 3, face_num /) )
  real ( kind = 8 ), dimension ( 3, face_num ) :: face_normal = reshape ( (/ &
    0.0D+00,  0.0D+00, -1.0D+00, &
    0.0D+00,  0.0D+00, -1.0D+00, &
    0.0D+00, -1.0D+00,  0.0D+00, &
    0.0D+00, -1.0D+00,  0.0D+00, &
    0.0D+00, +1.0D+00,  0.0D+00, &
    0.0D+00, +1.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00, +1.0D+00, &
    0.0D+00,  0.0D+00, +1.0D+00, &
   -1.0D+00,  0.0D+00,  0.0D+00, &
   -1.0D+00,  0.0D+00,  0.0D+00, &
   +1.0D+00,  0.0D+00,  0.0D+00, &
   +1.0D+00,  0.0D+00,  0.0D+00 /), (/ 3, face_num /) )
  real ( kind = 8 ), dimension ( 3, face_num ) :: face_normal2
  real ( kind = 8 ), dimension ( 3, node_num ) :: node_xyz = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00 /), (/ 3, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  STLA_FACE_NORMAL_COMPUTE computes the face normal'
  write ( *, '(a)' ) '  vectors for an STLA file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We have an STLA solid, and its exact normals.'
  write ( *, '(a)' ) '  We now call STLA_FACE_NORMAL_COMPUTE to '
  write ( *, '(a)' ) '  recompute the normals.'

  call stla_face_normal_compute ( node_num, face_num, node_xyz, &
    face_node, face_normal2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We print out the maximum error, defined as'
  write ( *, '(a)' ) '    |1 - dot ( n1, n2 )|'
  write ( *, '(a)' ) '  where n1 and n2 are the exact and computed normals.'

  dot_max = 0.0D+00

  do face = 1, face_num

    dot_max = max ( dot_max, &
      abs ( 1.0D+00 - &
      dot_product ( face_normal(1:3,face), face_normal2(1:3,face) ) ) )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum error = ', dot_max

  return
end
