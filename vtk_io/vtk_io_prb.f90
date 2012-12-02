program main

!*****************************************************************************80
!
!! MAIN is the main program for VTK_IO_PRB.
!
!  Discussion:
!
!    VTK_IO_PRB tests the routines in the VTK_IO library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VTK_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the VTK_IO library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VTK_IO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests VTK_PUVW_WRITE.
!
!  Discussion:
!
!    VTK_PUVW_WRITE writes pressure (P) and velocity (UVW) for a 3D
!    fluid flow calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: element_num = 6
  integer ( kind = 4 ), parameter :: element_order = 8
  integer ( kind = 4 ), parameter :: node_num = 24

  integer ( kind = 4 ) :: element_node(element_order,element_num) = reshape ( &
    (/ &
      1,  2,  5,  6, 13, 14, 17, 18, &
      2,  3,  6,  7, 14, 15, 18, 19, &
      3,  4,  7,  8, 15, 16, 19, 20, &
      5,  6,  9, 10, 17, 18, 21, 22, &
      6,  7, 10, 11, 18, 19, 22, 23, &
      7,  8, 11, 12, 19, 20, 23, 24  &
    /), (/ element_order, element_num /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) node
  character ( len = 80 ) output_filename
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) p(node_num)
  character ( len = 80 ) title
  real ( kind = 8 ) uvw(3,node_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) xyz(3,node_num)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  VTK_PUVW_WRITE writes 3d fluid data, pressure and '
  write ( *, '(a)' ) '  velocity, to a VTK file.'

  output_filename = 'puvw_data.txt'
  title = 'Sample data for VTK_PUVW_WRITE.'

  node = 0
  do k = 1, 2
    z = real ( k - 1, kind = 8 )
    do j = 1, 3
      y = real ( j - 1, kind = 8 )
      do i = 1, 4
        x = real ( i - 1, kind = 8 )
        node = node + 1
        xyz(1:3,node) = (/ x, y, z /)
        p(node_num) = 10.0D+00 * x
        uvw(1,node_num) = 2.0D+00 * x
        uvw(2,node_num) = 3.0D+00 * y
        uvw(3,node_num) = 4.0D+00 * z
      end do
    end do
  end do

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace' )

  call vtk_puvw_write ( output_unit, title, node_num, element_num, &
    element_order, xyz, element_node, p, uvw )

  close (  unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VTK_PUVW_WRITE created the file.'

  return
end
