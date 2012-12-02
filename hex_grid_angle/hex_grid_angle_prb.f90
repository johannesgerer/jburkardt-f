program main

!*****************************************************************************80
!
!! MAIN is the main program for HEX_GRID_ANGLE.
!
!  Discussion:
!
!    HEX_GRID_ANGLE_PRB calls a set of problems for HEX_GRID_ANGLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEX_GRID_ANGLE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HEX_GRID_ANGLE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEX_GRID_ANGLE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests HEX_GRID_ANGLE_01_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_uniform_01
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  HEX_GRID_ANGLE_01_SIZE computes the number of'
  write ( *, '(a)' ) '  points in a hexagonal grid in the unit square,'
  write ( *, '(a)' ) '  with grid lines at a given angle ANGLE,'
  write ( *, '(a)' ) '  with a given spacing H between points on a grid line,'
  write ( *, '(a)' ) '  with the coordinates of the center at CENTER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        CENTER           ANGLE       H         N'
  write ( *, '(a)' ) ' '

  do i = 1, 3

    call r8vec_uniform_01 ( dim_num, seed, center )

    angle = 45.0D+00
    h = 0.25D+00

    call hex_grid_angle_01_size ( center, angle, h, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,i6)' ) &
      center(1:dim_num), angle, h, n

  end do

  write ( *, '(a)' ) ' '

  do i = 1, 3

    center(1:dim_num) = 0.5D+00
    angle = 45.0D+00
    h = 0.25D+00 / real ( i, kind = 8 )

    call hex_grid_angle_01_size ( center, angle, h, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,i6)' ) &
      center(1:dim_num), angle, h, n

  end do

  write ( *, '(a)' ) ' '

  do i = 1, 4

    center(1:dim_num) = 0.5D+00
    angle = 180.D+00 * r8_uniform_01 ( seed )
    h = 0.25D+00

    call hex_grid_angle_01_size ( center, angle, h, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,i6)' ) &
      center(1:dim_num), angle, h, n

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HEX_GRID_ANGLE_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ) center(dim_num)
  integer ( kind = 4 ) n
  real ( kind = 8 ) h
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HEX_GRID_ANGLE_01 computes the'
  write ( *, '(a)' ) '  points in a hexagonal grid in the unit square,'
  write ( *, '(a)' ) '  with grid lines at a given angle ANGLE,'
  write ( *, '(a)' ) '  with a given spacing H between points on a grid line,'
  write ( *, '(a)' ) '  with the coordinates of the center at CENTER.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        CENTER           ANGLE       H         N'
  write ( *, '(a)' ) ' '

  center(1) = 0.5D+00
  center(2) = 0.0D+00
  angle = 45.0D+00
  h = 0.25D+00

  call hex_grid_angle_01_size ( center, angle, h, n )

  write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,i6)' ) &
    center(1:dim_num), angle, h, n

  allocate ( r(dim_num,n) )

  call hex_grid_angle_01 ( center, angle, h, n, r )

  call r8mat_transpose_print ( dim_num, n, r, '  Grid points:' )

  deallocate ( r )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests HEX_GRID_ANGLE_01_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ) center(dim_num)
  character ( len = 100 ) :: file_name = 'hex_grid_angle_01_dataset.txt'
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HEX_GRID_ANGLE_01_WRITE writes the points of a'
  write ( *, '(a)' ) '  angled hexagonal grid to a file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        CENTER           ANGLE       H         N'
  write ( *, '(a)' ) ' '

  do i = 1, 3

    call r8vec_uniform_01 ( dim_num, seed, center )

    angle = 180.0D+00 * r8_uniform_01 ( seed )
    h = 0.107457D+00

    call hex_grid_angle_01_size ( center, angle, h, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f10.6,2x,i6)' ) &
      center(1:dim_num), angle, h, n

  end do

  allocate ( r(dim_num,n) )

  call hex_grid_angle_01 ( center, angle, h, n, r )

  call hex_grid_angle_01_write ( center, angle, h, n, r, file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to "' &
    // trim ( file_name ) // '".'

  deallocate ( r )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests HEX_GRID_ANGLE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( dim_num, 2 ) :: box = reshape ( (/ &
    10.0D+00, 2.0D+00, &
    12.0D+00, 2.5D+00 /), (/ dim_num, 2 /) )
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  HEX_GRID_ANGLE_SIZE computes the number of'
  write ( *, '(a)' ) '  points in a hexagonal grid in the unit square,'
  write ( *, '(a)' ) '  with grid lines at a given angle ANGLE,'
  write ( *, '(a)' ) '  with a given spacing H between points on a grid line,'
  write ( *, '(a)' ) '  with the coordinates of the center at CENTER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        CENTER           ANGLE       H         N'
  write ( *, '(a)' ) ' '

  do i = 1, 3

    center(1) = box(1,1) + r8_uniform_01 ( seed ) * ( box(1,2) - box(1,1) )
    center(2) = box(2,1) + r8_uniform_01 ( seed ) * ( box(2,2) - box(2,1) )

    angle = 45.0D+00
    h = 0.25D+00

    call hex_grid_angle_size ( box, center, angle, h, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,i6)' ) &
      center(1:dim_num), angle, h, n

  end do

  write ( *, '(a)' ) ' '

  do i = 1, 3

    center(1:dim_num) = 0.5D+00 * ( box(1:dim_num,1) +  box(1:dim_num,2) )
    angle = 45.0D+00
    h = 0.25D+00 / real ( i, kind = 8 )

    call hex_grid_angle_size ( box, center, angle, h, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,i6)' ) &
      center(1:dim_num), angle, h, n

  end do

  write ( *, '(a)' ) ' '

  do i = 1, 4

    center(1:dim_num) = 0.5D+00 * ( box(1:dim_num,1) +  box(1:dim_num,2) )
    angle = 180.0D+00 * r8_uniform_01 ( seed )
    h = 0.25D+00

    call hex_grid_angle_size ( box, center, angle, h, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,i6)' ) &
      center(1:dim_num), angle, h, n

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests HEX_GRID_ANGLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( dim_num, 2 ) :: box = reshape ( (/ &
    10.0D+00, 2.0D+00, &
    12.0D+00, 2.5D+00 /), (/ dim_num, 2 /) )
  real ( kind = 8 ) center(dim_num)
  integer ( kind = 4 ) n
  real ( kind = 8 ) h
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  HEX_GRID_ANGLE computes the'
  write ( *, '(a)' ) '  points in a hexagonal grid in the unit square,'
  write ( *, '(a)' ) '  with grid lines at a given angle ANGLE,'
  write ( *, '(a)' ) '  with a given spacing H between points on a grid line,'
  write ( *, '(a)' ) '  with the coordinates of the center at CENTER.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        CENTER           ANGLE       H         N'
  write ( *, '(a)' ) ' '

  center(1:dim_num) = 0.5D+00 * ( box(1:dim_num,1) +  box(1:dim_num,2) )
  angle = 45.0D+00
  h = 0.25D+00

  call hex_grid_angle_size ( box, center, angle, h, n )

  write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,i6)' ) &
    center(1:dim_num), angle, h, n

  allocate ( r(dim_num,n) )

  call hex_grid_angle ( box, center, angle, h, n, r )

  call r8mat_transpose_print ( dim_num, n, r, '  Grid points:' )

  deallocate ( r )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests HEX_GRID_ANGLE_01_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( dim_num, 2 ) :: box = reshape ( (/ &
    10.0D+00, 2.0D+00, &
    12.0D+00, 2.5D+00 /), (/ dim_num, 2 /) )
  real ( kind = 8 ) center(dim_num)
  character ( len = 100 ) :: file_name = 'hex_grid_angle_dataset.txt'
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  HEX_GRID_ANGLE_WRITE writes the points of a'
  write ( *, '(a)' ) '  angled hexagonal grid to a file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        CENTER           ANGLE       H         N'
  write ( *, '(a)' ) ' '

  do i = 1, 3

    center(1:dim_num) = 0.5D+00 * ( box(1:dim_num,1) +  box(1:dim_num,2) )
    angle = 180.0D+00 * r8_uniform_01 ( seed )
    h = 0.107457D+00

    call hex_grid_angle_size ( box, center, angle, h, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f10.6,2x,i6)' ) &
      center(1:dim_num), angle, h, n

  end do

  allocate ( r(dim_num,n) )

  call hex_grid_angle ( box, center, angle, h, n, r )

  call hex_grid_angle_write ( box, center, angle, h, n, r, file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to "' &
    // trim ( file_name ) // '".'

  deallocate ( r )

  return
end
