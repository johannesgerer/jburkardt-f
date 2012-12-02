program main

!*****************************************************************************80
!
!! MAIN is the main program for GNUFOR_PRB.
!
!  Discussion:
!
!    GNUFOR_PRB tests the GNUFOR routines.
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

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GNUFOR_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GNUFOR library.'

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
  write ( *, '(a)' ) 'GNUFOR_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the plotting of Y(X) data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 101

  real ( kind = 8 ) angle
  real ( kind = 8 ) area
  character ( len = 255 ) :: command_filename = 'test01_commands.txt'
  character ( len = 255 ) :: data_filename = 'test01_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: xmin = 0.0D+00
  real ( kind = 8 ), parameter :: xmax = 20.0D+00
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  To plot a simple set of (X,Y) data,'
  write ( *, '(a)' ) '  WRITE_XY_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XY_PLOT writes the plot command file.'

  do i = 1, n

    x(i) = ( real ( n - i,     kind = 8 ) * xmin   &
           + real (     i - 1, kind = 8 ) * xmax ) &
           / real ( n     - 1, kind = 8 )

    y(i) = sin ( x(i) ) * sin ( 4.0D+00 * x(i) )

  end do

  call write_xy_data ( data_filename, n, x, y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01'
    write ( *, '(a,i6)' ) '  WRITE_XY_DATA returned IERROR = ', ierror
  end if

  call write_xy_plot ( command_filename, data_filename, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01'
    write ( *, '(a,i6)' ) '  WRITE_XY_PLOT returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_filename )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the plotting of a table of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nrow = 101
  integer ( kind = 4 ), parameter :: ncol = 4

  integer ( kind = 4 ), parameter :: lda = nrow

  real ( kind = 8 ) angle
  real ( kind = 8 ) area
  character ( len = 255 ) :: command_filename = 'test02_commands.txt'
  character ( len = 255 ) :: data_filename = 'test02_data.txt'
  real ( kind = 8 ) height
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  real ( kind = 8 ), parameter :: r = 50.0D+00
  real ( kind = 8 ) width
  real ( kind = 8 ) x(lda,ncol)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  To plot X versus multiple sets of Y data,'
  write ( *, '(a)' ) '  WRITE_XYY_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XYY_PLOT writes the plot command file.'

  do i = 1, nrow

    height = 2.0d+00 * r * real ( i - 1, kind = 8 ) / real ( nrow - 1, kind = 8 )
    width = 2.0D+00 * sqrt ( r**2 - ( r - height )**2 )
    angle = acos ( ( r - height ) / r )
    area = 0.5D+00 * r**2 * 2.0d+00 * acos ( ( r - height ) / r ) &
      - ( r - height ) * sqrt ( height * ( 2.0D+00 * r - height ) )

    x(i,1) = height
    x(i,2) = width
    x(i,3) = angle
    x(i,4) = area

  end do

  call write_xyy_data ( data_filename, lda, nrow, ncol, x, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02'
    write ( *, '(a,i6)' ) '  WRITE_XYY_DATA returned IERROR = ', ierror
  end if

  call write_xyy_plots ( command_filename, data_filename, ncol, &
    ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02'
    write ( *, '(a,i6)' ) '  WRITE_XYY_PLOTS returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_filename )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 plots parameter (X,Y,Z) data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 101

  character ( len = 255 ) :: command_filename = 'test03_commands.txt'
  character ( len = 255 ) :: data_filename = 'test03_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), parameter :: nturn = 5
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r = 5.0D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)
  real ( kind = 8 ), parameter :: zmax = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  To plot a (parametric) set of (X,Y,Z) data,'
  write ( *, '(a)' ) '  WRITE_XYZ_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XYZ_PLOT writes the plot command file.'

  do i = 1, n

    z(i) = zmax * real ( i - 1, kind = 8 ) / real ( n - 1, kind = 8 )

    theta = ( 2.0D+00 * pi ) * z(i) * real ( nturn, kind = 8 ) / zmax

    x(i) = r * cos ( theta )
    y(i) = r * sin ( theta )

  end do

  call write_xyz_data ( data_filename, n, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03'
    write ( *, '(a,i6)' ) '  WRITE_XYZ_DATA returned IERROR = ', ierror
  end if

  call write_xyz_plot ( command_filename, data_filename, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03'
    write ( *, '(a,i6)' ) '  WRITE_XYZ_PLOT returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_filename )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 plots vector data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nx = 21
  integer ( kind = 4 ), parameter :: ny = 21
  integer ( kind = 4 ), parameter :: n = nx * ny

  character ( len = 255 ) :: command_filename = 'test04_commands.txt'
  character ( len = 255 ) :: data_filename = 'test04_data.txt'
  real ( kind = 8 ) dx(n)
  real ( kind = 8 ) dy(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: xmax = 1.0D+00
  real ( kind = 8 ), parameter :: xmin = -1.0D+00
  real ( kind = 8 ) xx
  real ( kind = 8 ) y(n)
  real ( kind = 8 ), parameter :: ymax = 1.0D+00
  real ( kind = 8 ), parameter :: ymin = -1.0D+00
  real ( kind = 8 ) yy

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  To plot a vector field,'
  write ( *, '(a)' ) '  WRITE_VECTOR_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_VECTOR_PLOT writes the plot command file.'

  k = 0

  do i = 1, nx

    do j = 1, ny

      k = k + 1

      xx = ( real ( nx - i,     kind = 8     ) * xmin   &
           + real (      i - 1, kind = 8 ) * xmax ) &
           / real ( nx     - 1, kind = 8 )

      yy = ( real ( ny - j,     kind = 8     ) * ymin   &
           + real (      j - 1, kind = 8 ) * ymax ) &
           / real ( ny     - 1, kind = 8 )

      dx(k) = - 0.10D+00 * yy
      dy(k) =   0.10D+00 * xx

      x(k) = xx - 0.5D+00 * dx(k)
      y(k) = yy - 0.5D+00 * dy(k)

    end do

  end do

  call write_vector_data ( data_filename, n, x, y, dx, dy, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04'
    write ( *, '(a,i6)' ) '  WRITE_VECTOR_DATA returned IERROR = ', ierror
  end if

  call write_vector_plot ( command_filename, data_filename, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04'
    write ( *, '(a,i6)' ) '  WRITE_VECTOR_PLOT returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_filename )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 plots Z(X,Y) grid data as a surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nx = 21
  integer ( kind = 4 ), parameter :: ny = 21

  integer ( kind = 4 ), parameter :: nrow = nx * ny

  character ( len = 255 ) :: command_filename = 'test05_commands.txt'
  character ( len = 255 ) :: data_filename = 'test05_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xmax = 1.0D+00
  real ( kind = 8 ), parameter :: xmin = 0.0D+00
  real ( kind = 8 ) xyz(3,nx,ny)
  real ( kind = 8 ) y
  real ( kind = 8 ), parameter :: ymax = 1.0D+00
  real ( kind = 8 ), parameter :: ymin = 0.0D+00
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  To plot a gridded set of Z(X,Y) data as a surface,'
  write ( *, '(a)' ) '  WRITE_XYZGRID_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XYZGRID_SURFACE writes the plot command file.'

  do i = 1, nx

    x = ( real ( nx - i,     kind = 8     ) * xmin   &
        + real (      i - 1, kind = 8 ) * xmax ) &
        / real ( nx     - 1, kind = 8 )

    do j = 1, ny

      y = ( real ( ny - j,     kind = 8     ) * ymin   &
          + real (      j - 1, kind = 8 ) * ymax ) &
          / real ( ny     - 1, kind = 8 )

      z = sin ( 64.0D+00 * ( x - 0.5D+00 )**2 * ( y - 0.5D+00 )**2 )

      xyz(1,i,j) = x
      xyz(2,i,j) = y
      xyz(3,i,j) = z

    end do

  end do

  call write_xyzgrid_data ( data_filename, nx, ny, xyz, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05'
    write ( *, '(a,i6)' ) '  WRITE_XYZGRID_DATA returned IERROR = ', ierror
  end if

  call write_xyzgrid_surface ( command_filename, data_filename, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05'
    write ( *, '(a,i6)' ) '  WRITE_XYZGRID_SURFACE returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_filename )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 plots Z(X,Y) grid data as contours.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nx = 41
  integer ( kind = 4 ), parameter :: ny = 41

  character ( len = 255 ) :: command_filename = 'test06_commands.txt'
  character ( len = 255 ) :: data_filename = 'test06_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xmax = 1.0D+00
  real ( kind = 8 ), parameter :: xmin = 0.0D+00
  real ( kind = 8 ) xyz(3,nx,ny)
  real ( kind = 8 ) y
  real ( kind = 8 ), parameter :: ymax = 1.0D+00
  real ( kind = 8 ), parameter :: ymin = 0.0D+00
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  To plot gridded Z(X,Y) data as contours,'
  write ( *, '(a)' ) '  WRITE_XYZGRID_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XYZGRID_CONTOUR writes the plot command file.'

  do i = 1, nx

    x = ( real ( nx - i,     kind = 8     ) * xmin   &
        + real (      i - 1, kind = 8 ) * xmax ) &
        / real ( nx     - 1, kind = 8 )

    do j = 1, ny

      y = ( real ( ny - j,     kind = 8     ) * ymin   &
          + real (      j - 1, kind = 8 ) * ymax ) &
          / real ( ny     - 1, kind = 8 )

      z = sin ( 64.0D+00 * ( x - 0.5D+00 )**2 * ( y - 0.5D+00 )**2 )

      xyz(1,i,j) = x
      xyz(2,i,j) = y
      xyz(3,i,j) = z

    end do

  end do

  call write_xyzgrid_data ( data_filename, nx, ny, xyz, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06'
    write ( *, '(a,i6)' ) '  WRITE_XYZGRID_DATA returned IERROR = ', ierror
  end if

  call write_xyzgrid_contour ( command_filename, data_filename, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06'
    write ( *, '(a,i6)' ) '  WRITE_XYZGRID_CONTOUR returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_filename )

  return
end
