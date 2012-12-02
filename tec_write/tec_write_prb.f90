program main

!*****************************************************************************80
!
!! MAIN is the main program for TEC_WRITE_PRB.
!
!  Discussion:
!
!    TEC_WRITE_PRB tests the routines in TEC_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_WRITE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEC_WRITE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_WRITE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 reads data from a file, and writes a 3D cylinder version for Tecplot.
!
!  Discussion:
!
!    The data represents flow in a cylinder.  The data is assumed to be
!    cylindrically symmetric, and data (VR,VZ,VT) was computed only for
!    a single plane of THETA = 0.
!
!    The plane of data was computed for a number of times.
!
!    To make a "nice" 3D version of the data, the plane of data is written
!    out at NT equally spaced values of THETA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nr = 17
  integer, parameter :: nz = 17

  integer ierror
  character ( len = 255 ) :: input_file = 'mydata.txt'
  integer ios
  integer iunit1
  integer iunit2
  integer, parameter :: nt = 17
  character ( len = 255 ) :: output_file = 'cylinder_3d.dat'
  real ( kind = 8 ) r(nr)
  real ( kind = 8 ), parameter :: rlo = 0.0D+00
  real ( kind = 8 ), parameter :: rhi = 1.0D+00
  integer time_step
  real ( kind = 8 ) vr(nr,nz)
  real ( kind = 8 ) vt(nr,nz)
  real ( kind = 8 ) vz(nr,nz)
  real ( kind = 8 ) z(nz)
  real ( kind = 8 ), parameter :: zlo = 0.0D+00
  real ( kind = 8 ), parameter :: zhi = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Read my data file.'
  write ( *, '(a)' ) '  Write a Tecplot file with 3D cylinder data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of copies of plane data = ', nt
!
!  The R and Z data is simply evenly spaced between given limits.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Create evenly spaced R and Z data.'

  call r8vec_even ( rlo, rhi, nr, r )
  call r8vec_even ( zlo, zhi, nz, z )

  call get_unit ( iunit1 )

  open ( unit = iunit1, file = input_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input data file.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the Tecplot file : ' // trim ( output_file )

  call tec_write_open ( output_file, iunit2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the Tecplot header.'

  call tec_write_header ( iunit2, &
    'Cylindrical Flow Velocity Data', &
    '"X","Y","Z","VX","VY", "VZ"' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a set of data records for a sequence of times.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each data record is a (VR,VZ,VT) velocity field'
  write ( *, '(a)' ) '  in the plane (R,Z,0).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The routine TEC_WRITE_CYL_V writes NT copies of'
  write ( *, '(a)' ) '  this data, so that TECPLOT will see a cylindically'
  write ( *, '(a)' ) '  symmetric filled in version of the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It also automatically converts the data from '
  write ( *, '(a)' ) '  (R,Z,T) and (VR,VZ,VT) format to'
  write ( *, '(a)' ) '  (X,Y,Z) and (VX,VY,VZ) format.'

  time_step = 0

  do

    call mydata_read ( iunit1, nr, nz, vr, vz, vt, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

    time_step = time_step + 1
    write ( *, '(a,i6)' ) '  Handling data for time step = ', time_step

    call tec_write_cyl_v ( iunit2, nr, nz, nt, r, z, vr, &
      vz, vt )

  end do

  close ( unit = iunit1 )

  call tec_write_close ( iunit2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of time steps was ', time_step

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 reads data from a file, and writes a 2D version of it for Tecplot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nx = 17
  integer, parameter :: ny = 17

  integer ierror
  character ( len = 255 ) :: input_file = 'mydata.txt'
  integer ios
  integer iunit1
  integer iunit2
  integer, parameter :: nt = 1
  character ( len = 255 ) :: output_file = 'plane_slice.dat'
  integer time_step
  real ( kind = 8 ) vx(nx,ny)
  real ( kind = 8 ) vy(nx,ny)
  real ( kind = 8 ) vz(nx,ny)
  real ( kind = 8 ) x(nx)
  real ( kind = 8 ), parameter :: xlo = 0.0D+00
  real ( kind = 8 ), parameter :: xhi = 1.0D+00
  real ( kind = 8 ) y(ny)
  real ( kind = 8 ), parameter :: ylo = 0.0D+00
  real ( kind = 8 ), parameter :: yhi = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Read my data file.'
  write ( *, '(a)' ) '  Write a Tecplot file with one plane of data.'
!
!  Set X and Y to be dummy values.
!
  write ( *, '(a)' ) '  Create X and Y data.'

  call r8vec_even ( xlo, xhi, nx, x )
  call r8vec_even ( ylo, yhi, ny, y )

  call get_unit ( iunit1 )

  open ( unit = iunit1, file = input_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input data file.'
    return
  end if

  write ( *, '(a)' ) '  Writing the Tecplot file : ' // trim ( output_file )

  call tec_write_open ( output_file, iunit2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the Tecplot header.'

  call tec_write_header ( iunit2, &
    'Cylindrical Flow Velocity Data', &
    '"X","Y","VX","VY","VZ"' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Begin reading sets of data.'

  time_step = 0

  do

    call mydata_read ( iunit1, nx, ny, vx, vy, vz, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

    time_step = time_step + 1
    write ( *, '(a,i6)' ) '  Writing data for time step = ', time_step

    call tec_write_xy_uvw ( iunit2, nx, ny, x, y, vx, vy, vz )

  end do

  close ( unit = iunit1 )

  call tec_write_close ( iunit2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of time steps was ', time_step

  return
end
subroutine mydata_read ( iunit, nr, nz, vr, vz, vt, ierror )

!*****************************************************************************80
!
!! MYDATA_READ reads one set of data from a file.
!
!  Discussion:
!
!    This data file contains two sets of data for a 17 by 17 grid.
!    The grid is a slice through a cylinder.  At each point on the
!    grid, the flow velocity is recorded in the R, Theta and Z directions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN unit number associated with the file.
!
!    Input, integer NR, the number of grid points in the R direction
!    (from the center out).
!
!    Input, integer NZ, the number of grid points in the Z direction (up).
!
!    Output, real ( kind = 8 ) VR(NR,NZ), VZ(NR,NZ), VT(NR,NZ), the radial,
!    vertical, and out-of-plane velocity components at each grid point.
!
  implicit none

  integer nz
  integer nr

  integer i
  integer ierror
  integer ios
  integer iunit
  integer j
  real ( kind = 8 ) vr(nr,nz)
  real ( kind = 8 ) vt(nr,nz)
  real ( kind = 8 ) vz(nr,nz)

  i = 1
  j = 1

  do

    read ( iunit, *, iostat = ios ) vr(i,j), vt(i,j), vz(i,j)

    if ( ios /= 0 ) then
      ierror = 1
      exit
    end if

    i = i + 1
    if ( i > nr ) then
      i = 1
      j = j + 1
      if ( j > nz ) then
        ierror = 0
        exit
      end if
    end if

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 reads data from a file, and writes a 3D slice version for TECPLOT.
!
!  Discussion:
!
!    NT is the number of copies of the data that will be made in the THETA
!    direction.  Actually, to keep Tecplot happy, you want to make one copy
!    at THETA = 0 and one at THETA = 2*PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nr = 17
  integer, parameter :: nz = 17

  integer ierror
  character ( len = 255 ) :: input_file = 'mydata.txt'
  integer ios
  integer iunit1
  integer iunit2
  integer, parameter :: nt = 1
  character ( len = 255 ) :: output_file = 'cylinder_slice.dat'
  real ( kind = 8 ) r(nr)
  real ( kind = 8 ), parameter :: rlo = 0.0D+00
  real ( kind = 8 ), parameter :: rhi = 1.0D+00
  integer time_step
  real ( kind = 8 ) vr(nr,nz)
  real ( kind = 8 ) vt(nr,nz)
  real ( kind = 8 ) vz(nr,nz)
  real ( kind = 8 ) z(nz)
  real ( kind = 8 ), parameter :: zlo = 0.0D+00
  real ( kind = 8 ), parameter :: zhi = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Read my data file.'
  write ( *, '(a)' ) '  Write a Tecplot file with 3D cylinder data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of copies of plane data = ', nt
!
!  Set R and H to be dummy values.
!
  write ( *, '(a)' ) '  Create R and H data.'

  call r8vec_even ( rlo, rhi, nr, r )
  call r8vec_even ( zlo, zhi, nz, z )

  call get_unit ( iunit1 )

  open ( unit = iunit1, file = input_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input data file.'
    return
  end if

  write ( *, '(a)' ) '  Writing the Tecplot file : ' // trim ( output_file )

  call tec_write_open ( output_file, iunit2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the Tecplot header.'

  call tec_write_header ( iunit2, &
    'Cylindrical Flow Velocity Data', &
    '"X","Y","Z","VX","VY", "VZ"' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Begin reading sets of data.'

  time_step = 0

  do

    call mydata_read ( iunit1, nr, nz, vr, vz, vt, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

    time_step = time_step + 1
    write ( *, '(a,i6)' ) '  Writing data for time step = ', time_step

    call tec_write_cyl_v ( iunit2, nr, nz, nt, r, z, vr, vz, vt )

  end do

  close ( unit = iunit1 )

  call tec_write_close ( iunit2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of time steps was ', time_step

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 demonstrates writing PRUV data for a rectangular region.
!
!  Discussion:
!
!    The data exists on a regular grid of points.
!
!    Assume we have flow in a box from (0,0) to (10,3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nx = 11
  integer, parameter :: ny = 4

  integer iunit
  character ( len = 255 ) :: output_file = 'channel.dat'
  real ( kind = 8 ) p(nx,ny)
  real ( kind = 8 ) rho(nx,ny)
  real ( kind = 8 ) u(nx,ny)
  real ( kind = 8 ) v(nx,ny)
  real ( kind = 8 ) x(nx,ny)
  real ( kind = 8 ) xbot(2)
  real ( kind = 8 ) xtop(2)
  real ( kind = 8 ) y(nx,ny)
  real ( kind = 8 ) ybot(2)
  real ( kind = 8 ) ytop(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Assume we have PRUV data on a regular grid '
  write ( *, '(a)' ) '  of XY points.'
  write ( *, '(a)' ) '  Do not assume any underlying finite element structure.'
!
!  Get the data values.
!
  call example_xy_node ( nx, ny, x, y, p, rho, u, v )

  write ( *, '(a)' ) '  Writing the Tecplot file : ' // trim ( output_file )

  call tec_write_open ( output_file, iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the Tecplot header.'

  call tec_write_header ( iunit, &
    'Flow in a Channel', &
    '"X","Y","P","RHO","U", "V"' )

  call tec_write_xy_pruv ( iunit, nx, ny, x, y, p, rho, u, v )

  xtop(1) = 0.0
  ytop(1) = 3.0
  xtop(2) = 10.0
  ytop(2) = 3.0

  call tec_write_xy_line ( iunit, 2, xtop, ytop )

  xbot(1) = 0.0
  ybot(1) = 0.0
  xbot(2) = 10.0
  ybot(2) = 0.0

  call tec_write_xy_line ( iunit, 2, xbot, ybot )
!
!  I really want THIS line to be DASHED, not SOLID, and BLUE, not BLACK.
!  The easiest way to do that is to edit the result TEC file.
!
  xbot(1) = 7.0
  ybot(1) = 0.0
  xbot(2) = 7.0
  ybot(2) = 3.0

  call tec_write_xy_line ( iunit, 2, xbot, ybot )

  call tec_write_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  The TECPLOT file has been written.'

  return
end
subroutine example_xy_node ( nx, ny, x, y, p, rho, u, v )

!*****************************************************************************80
!
!! EXAMPLE_NODE sets up the example node data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the number of nodes in the X and Y directions.
!
!    Output, real ( kind = 8 ) X(NX,NY), Y(NX,NY), the X and Y coordinates
!    of the nodes.
!
!    Output, real ( kind = 8 ) P(NX,NY), the pressure at the nodes.
!
!    Output, real ( kind = 8 ) RHO(NX,NY), the density at the nodes.
!
!    Output, real ( kind = 8 ) U(NX,NY), V(NX,NY), the X and Y velocity
!    components at the nodes.
!
  implicit none

  integer nx
  integer ny

  integer i
  integer j
  real ( kind = 8 ) p(nx,ny)
  real ( kind = 8 ) rho(nx,ny)
  real ( kind = 8 ) u(nx,ny)
  real ( kind = 8 ) v(nx,ny)
  real ( kind = 8 ) x(nx,ny)
  real ( kind = 8 ), parameter:: xlen = 10.0D+00
  real ( kind = 8 ) y(nx,ny)
  real ( kind = 8 ), parameter :: ylen = 3.0D+00

  do i = 1, nx

    do j = 1, ny

      x(i,j) = real ( i - 1, kind = 8 ) * xlen / real ( nx - 1, kind = 8 )
      y(i,j) = real ( j - 1, kind = 8 ) * ylen / real ( ny - 1, kind = 8 )
      rho(i,j) = 1.0D+00 + ( xlen - x(i,j) ) * ( ylen - y(i,j) )
      u(i,j) = y(i,j) * ( ylen - y(i,j) )
      v(i,j) = 0.01D+00 * sin ( 6.28D+00 * x(i,j) / xlen )
      p(i,j) = - 2.0D+00 * x(i,j)

    end do

  end do

  return
end
