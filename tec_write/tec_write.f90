subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen
 
  iunit = 0
 
  do i = 1, 99
 
    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )
 
      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if
 
  end do

  return
end
subroutine r8vec_even ( alo, ahi, n, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  integer i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine tec_write_open ( file_name, iunit )

!*****************************************************************************80
!
!! TEC_WRITE_OPEN opens a TECPLOT ASCII file for output.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the output file to create.
!
!    Output, integer IUNIT, the FORTRAN unit number to be used with the file.
!
  implicit none

  character ( len = * ) file_name
  integer ios
  integer iunit
!
!  Get a free FORTRAN unit and open the output file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, form = 'formatted', &
     access = 'sequential', status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEC_WRITE_OPEN - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  return
end
subroutine tec_write_close ( iunit )

!*****************************************************************************80
!
!! TEC_WRITE_CLOSE closes a TECPLOT ASCII file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN unit number that was used with the file.
!
  implicit none

  integer iunit

  close ( unit = iunit )

  return
end
subroutine tec_write_header ( iunit, title, variables )

!*****************************************************************************80
!
!! TEC_WRITE_HEADER writes the two line header for a TECPLOT ASCII file.
!
!  Discussion:
!
!    You are free to specify any title you want.  The routine will add a
!    "Title" keyword, and delimit the title with double quotes.  For
!    example, if the input value of TITLE is 'FLOW6 data', then what will
!    appear as the title line in the TECPLOT file is:
!
!      Title="FLOW6 data"
!
!    The variables list is significant, since TECPLOT will count
!    the number of names, and assume that a corresponding set of values
!    will follow.  VARIABLES should contain a list of variable names,
!    each in double quotes, separated by commas.  The first two or three
!    variables are usually the spatial coordinates.  A typical value
!    of VARIABLES might be '"X","Y","P","U","V"', in which case this
!    routine will write the line:
!
!      Variables="X","Y","P","U","V"
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
!    Input, integer IUNIT, the FORTRAN unit associated with the file.
!
!    Input, character ( len = * ) TITLE, the title.
!
!    Input, character ( len = * ) VARIABLES, the variable name line.
!
  implicit none

  integer iunit
  character ( len = * ) title
  character ( len = * ) variables

  write ( iunit, '(a)' ) 'Title="' // trim ( title ) // '"'

  write ( iunit, '(a)' ) 'Variables=' // trim ( variables )

  return
end
subroutine tec_write_xy_line ( iunit, np, x, y )

!*****************************************************************************80
!
!! TEC_WRITE_XY_LINE writes out line data in XY geometry for use by TECPLOT.
!
!  Discussion:
!
!    The data is written as a GEOMETRY record.
!
!    The GEOMETRY record arguments, which you might want to adjust, include:
!
!      X = 0.0,
!      Y = 0.0, a fixed offset to be added to the (X,Y) data;
!      T = LINE, specifies that a line is being drawn;
!      CS = GRID, specifies that X and Y are measured in physical units;
!      L = SOLID, chooses a solid line type;
!      LT = 0.005, sets the line thickness, in frame units;
!      C = BLACK, chooses the line color to be black;
!      FC = BLACK, chooses the fill color, to be used if the line forms an area;
!      F = POINT, specifies that data will be (X1,Y1), (X2,Y2), ...(XN,YN);
!      S = GLOBAL, says that this line should appear in all "like" frames;
!
!    I am not completely clear on the filled areas.  In particular, I don't
!    understand whether:
!
!      A) the points are drawn, and if they enclose an area, the area is filled
!         with color FC;
!
!      B) the points are drawn, and the last point connected to the first,
!         and then, if FC is specified, that area is filled.
!
!      C) the points are drawn, and if FC is specified, the last point is
!         connected to the first and the area is filled.
!
!    The ZN parameter can be used to attact the line to a specific zone
!    or XY mapping.
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
!  Reference:
!
!    Section 3.5.2, "Geometry Record",
!    TECPLOT User's Manual, Version 7,
!    AMTEC Engineering, August 1996.
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN unit number associated with the file.
!
!    Input, integer NP, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) Y(NP), the Y coordinates of the nodes.
!
  implicit none

  integer np

  integer i
  integer iunit
  real ( kind = 8 ) x(np)
  real ( kind = 8 ) y(np)

  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'
  write ( iunit, '(a)' ) '  L = SOLID, LT = 0.005, C = BLACK, FC = BLACK, '
  write ( iunit, '(a)' ) '  F = POINT, S = GLOBAL'
  write ( iunit, '(a)' ) '1'
  write ( iunit, '(i6)' ) np

  do i = 1, np
    write ( iunit, '(2g15.6)' ) x(i), y(i)
  end do

  return
end
subroutine tec_write_xy_puv ( iunit, nelem, np, x, y, p, u, v )

!*****************************************************************************80
!
!! TEC_WRITE_XY_PUV writes out PUV data in XY geometry for use by TECPLOT.
!
!  Discussion:
!
!    The data format used is FEDATA, or "finite element" data.
!
!    "PUV" is intended to indicate pressure and X, Y velocity components.
!
!    Before this routine, you should have called
!
!      TEC_WRITE_HEADER
!
!    to write the title and variable names.
!
!    After this routine writes the node data, you should call
!
!      TEC_WRITE_T3 
!
!    or 
!
!      TEC_WRITE_T6
!
!    to record the element data.
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
!    Input, integer IUNIT, the FORTRAN unit number associated with the file.
!
!    Input, integer NELEM, the number of elements.  If you are using 3 node
!    triangles, then NELEM should be this number.  If you are using 6 node
!    triangles, then NELEM should be 4 times that number, since each 6 node
!    triangle you were using will be reported to TECPLOT as 4 triangles 
!    of 3 nodes.
!
!    Input, integer NP, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NP), Y(NP), the X and Y coordinates of
!    the nodes.
!
!    Input, real ( kind = 8 ) P(NP), the pressure.
!
!    Input, real ( kind = 8 ) U(NP), V(NP), the X and Y components of velocity.
!
  implicit none

  integer np

  integer i
  integer iunit
  integer nelem
  real ( kind = 8 ) p(np)
  real ( kind = 8 ) u(np)
  real ( kind = 8 ) v(np)
  real ( kind = 8 ) x(np)
  real ( kind = 8 ) y(np)

  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a,i6,a,i6,a)' ) 'Zone N=', np, ', E=', nelem, &
    ', F=FEPOINT, ET=TRIANGLE'
!
!  Write out the data at each node.
!
  do i = 1, np
    write ( iunit, '(5g15.6)' ) x(i), y(i), p(i), u(i), v(i)
  end do

  return
end
subroutine tec_write_t3 ( iunit, ntri, node )

!*****************************************************************************80
!
!! TEC_WRITE_T3 writes data defining 3 node triangular elements.
!
!  Discussion:
!
!    If you have specified that your data format is FEDATA, then TECPLOT 
!    will be expecting this information, which explains how the nodes are 
!    arranged into triangles to form elements.
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
!    Input, integer IUNIT, the FORTRAN unit number associated with the file.
!
!    Input, integer NTRI, the number of triangles.
!
!    Input, integer NODE(3,NTRI), gives, for each triangle, the indices of
!    the three nodes that make it up.
!
  implicit none

  integer ntri

  integer iunit
  integer j
  integer node(3,ntri)

  do j = 1, ntri
    write ( iunit, '(3i6)' ) node(1,j), node(2,j), node(3,j)
  end do

  return
end
subroutine tec_write_t6 ( iunit, ntri, node )

!*****************************************************************************80
!
!! TEC_WRITE_T6 rewrites 6 node triangular elements as 3 node triangles.
!
!  Discussion:
!
!    If you have specified that your data format is FEDATA, then TECPLOT 
!    will be expecting this information, which explains how the nodes are 
!    arranged into triangles to form elements.
!
!    If the data was computed on 6 node triangles, but TECPLOT can only
!    handle 3 node triangles, so this routine breaks up the 6 node triangle
!    into 4 smaller ones:
!
!    If your numbering of the nodes in the NODE array is different, you
!    will need to adjust the code.
!
!  Diagram:
!
!    6 node    3 node
!    triangle  triangles
!    --------  ---------
!
!    3          3
!    |\         |\
!    | \        | \
!    6  5       6--5
!    |   \      | /|\
!    |    \     |/ | \
!    1--4--2    1--4--2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN unit number associated with the file.
!
!    Input, integer NTRI, the number of triangles.
!
!    Input, integer NODE(6,NTRI), gives, for each triangle, the indices of
!    the six nodes that make it up.  Consult the diagram to see the order
!    that is assumed.
!
  implicit none

  integer ntri

  integer iunit
  integer j
  integer node(6,ntri)

  do j = 1, ntri
    write ( iunit, '(3i6)' ) node(1,j), node(4,j), node(6,j)
    write ( iunit, '(3i6)' ) node(2,j), node(5,j), node(4,j)
    write ( iunit, '(3i6)' ) node(3,j), node(6,j), node(5,j)
    write ( iunit, '(3i6)' ) node(4,j), node(5,j), node(6,j)
  end do

  return
end
subroutine tec_write_xy_uv ( iunit, nx, ny, x, y, u, v )

!*****************************************************************************80
!
!! TEC_WRITE_XY_UV writes out UV data in XY geometry for TECPLOT.
!
!  Discussion:
!
!    It is assumed that a Cartesian XY coordinate system ( X, Y ) is being
!    used, and that 2D velocity data ( U, V ) has been calculated over a grid
!
!      ( X(1:NX), Y(1:NY) )
!
!    This program appends the data to an open TECPLOT file, writing the 
!    data as a single "zone".
!
!    The data format used is POINT data, so no finite element data needs 
!    to be included in the file.
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
!    Input, integer NX, the number of grid points in the X direction.
!
!    Input, integer NY, the number of grid points in the Y direction.
!
!    Input, real ( kind = 8 ) X(NX), the X coordinates of the grid lines.
!
!    Input, real ( kind = 8 ) Y(NY), the Y coordinates of the grid lines.
!
!    Input, real ( kind = 8 ) U(NX,NY), V(NX,NY), the X and Y velocity
!    components.
!
  implicit none

  integer nx
  integer ny

  integer i
  integer iunit
  integer j
  real ( kind = 8 ) u(nx,ny)
  real ( kind = 8 ) v(nx,ny)
  real ( kind = 8 ) x(nx)
  real ( kind = 8 ) y(ny)
!
!  Write the zone header.
!
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a,i6,a,i6,a)' ) 'Zone I=', ny, ', J=', nx, ', F=POINT'
!
!  Write the zone data, one node at a time.
!
  do i = 1, nx
    do j = 1, ny

      write ( iunit, '(2f10.3,2g15.6)' ) x(i), y(j), u(i,j), v(i,j)

    end do
  end do

  return
end
subroutine tec_write_xy_pruv ( iunit, nx, ny, x, y, p, rho, u, v )

!*****************************************************************************80
!
!! TEC_WRITE_XY_PRUV writes out PRUV data in XY geometry for TECPLOT.
!
!  Discussion:
!
!    At a set of nodes scattered over a 2D region, the values of
!    X, Y, pressure, density and velocity are given.
!
!    This program appends the data to an open TECPLOT file, writing the 
!    data as a single "zone".
!
!    The data format used is POINT data, so no finite element data needs 
!    to be included in the file.
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
!    Input, integer IUNIT, the FORTRAN unit number associated with the file.
!
!    Input, integer NX, the number of grid points in the X direction.
!
!    Input, integer NY, the number of grid points in the Y direction.
!
!    Input, real ( kind = 8 ) X(NX,NY), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) Y(NY,NY), the Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) P(NX,NY), the pressure at the nodes.
!
!    Input, real ( kind = 8 ) RHO(NX,NY), the densities at the nodes.
!
!    Input, real ( kind = 8 ) U(NX,NY), V(NX,NY), the X and Y 
!    velocity components.
!
  implicit none

  integer nx
  integer ny

  integer i
  integer iunit
  integer j
  real ( kind = 8 ) p(nx,ny)
  real ( kind = 8 ) rho(nx,ny)
  real ( kind = 8 ) u(nx,ny)
  real ( kind = 8 ) v(nx,ny)
  real ( kind = 8 ) x(nx,ny)
  real ( kind = 8 ) y(nx,ny)
!
!  Write the zone header.
!
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a,i6,a,i6,a)' ) 'Zone I=', ny, ', J=', nx, ', F=POINT'
!
!  Write the zone data, one node at a time.
!
  do i = 1, nx
    do j = 1, ny

      write ( iunit, '(2f10.3,4g15.6)' ) x(i,j), y(i,j), p(i,j), rho(i,j), &
        u(i,j), v(i,j)

    end do
  end do

  return
end
subroutine tec_write_xy_uvw ( iunit, nx, ny, x, y, u, v, w )

!*****************************************************************************80
!
!! TEC_WRITE_XY_UVW writes out UVW data in XY geometry for TECPLOT.
!
!  Discussion:
!
!    Before this routine, you should have called
!
!      TEC_WRITE_HEADER
!
!    to write the title and variable names.
!
!    It is assumed that a Cartesian XY coordinate system ( X, Y ) is being
!    used, and that 3D velocity data ( U, V, W ) has been calculated 
!    over a grid
!
!      ( X(1:NX), Y(1:NY) )
!
!    This program appends the data to an open TECPLOT file, writing the 
!    data as a single "zone".
!
!    The data format used is POINT data, so no finite element data needs 
!    to be included in the file.
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
!    Input, integer NX, the number of grid points in the X direction.
!
!    Input, integer NY, the number of grid points in the Y direction.
!
!    Input, real ( kind = 8 ) X(NX), the X coordinates of the grid lines.
!
!    Input, real ( kind = 8 ) Y(NY), the Y coordinates of the grid lines.
!
!    Input, real ( kind = 8 ) U(NX,NY), V(NX,NY), W(NX,NY), the X, Y and Z 
!    velocity components.
!
  implicit none

  integer nx
  integer ny

  integer i
  integer iunit
  integer j
  real ( kind = 8 ) u(nx,ny)
  real ( kind = 8 ) v(nx,ny)
  real ( kind = 8 ) w(nx,ny)
  real ( kind = 8 ) x(nx)
  real ( kind = 8 ) y(ny)
!
!  Write the zone header.
!
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a,i6,a,i6,a)' ) 'Zone I=', ny, ', J=', nx, ', F=POINT'
!
!  Write the zone data, one node at a time.
!
  do i = 1, nx
    do j = 1, ny

      write ( iunit, '(2f10.3,2g15.6)' ) x(i), y(j), u(i,j), v(i,j), w(i,j)

    end do
  end do

  return
end
subroutine tec_write_cyl_v ( iunit, nr, nz, nt, r, z, vr, vz, vt )

!*****************************************************************************80
!
!! TEC_WRITE_CYL_V writes out V data in cylindrical geometry for TECPLOT.
!
!  Discussion:
!
!    It is assumed that a cylindrical coordinate system ( R, Z, T ) is being 
!    used, and that data has been calculated over a grid 
!
!      ( R(1:NR), Z(1:NZ), 0.0 )
!
!    with the assumption that the data is symmetric in T.
!
!    This program takes the computed 2D data, essentially converts it to
!    3D data by making NT copies of it, and appends the data to an
!    open TECPLOT file, writing the data as a single "zone".
!
!    The data format used is POINT data.
!
!    The data format used is POINT data, so no finite element data needs 
!    to be included in the file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2001
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
!    Input, integer NT, the number of grid points to create in the
!    THETA (out of plane) direction.
!
!    Input, real ( kind = 8 ) R(NR), the R coordinates of the grid lines.
!
!    Input, real ( kind = 8 ) Z(NZ), the Z coordinates of the grid lines.
!
!    Input, real ( kind = 8 ) VR(NR,NZ), VZ(NR,NZ), VT(NR,NZ), the radial,
!    vertical, and out-of-plane velocity components at each grid point.
!
  implicit none

  integer nr
  integer nt
  integer nz

  integer i
  integer ios
  integer iunit
  integer j
  integer k
  real ( kind = 8 ) r(nr)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) vr(nr,nz)
  real ( kind = 8 ) vt(nr,nz)
  real ( kind = 8 ) vx
  real ( kind = 8 ) vy
  real ( kind = 8 ) vz(nr,nz)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z(nz)
!
!  Compute the T grid coordinates.
!
  call tvec_even3 ( nt, t )
!
!  Write the zone header.
!
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a,i6,a,i6,a,i6,a)' ) 'Zone I=', nr, ', J=', nz, 'K=', &
    nt, ', F=POINT'
!
!  Write the zone data, one node at a time.
!
  do k = 1, nt
    do j = 1, nz
      do i = 1, nr

        x = r(i) * cos ( t(k) )
        y = r(i) * sin ( t(k) )

        vx = vr(i,j) * cos ( t(k) ) - vt(i,j) * sin ( t(k) )
        vy = vr(i,j) * sin ( t(k) ) + vt(i,j) * cos ( t(k) )

        write ( iunit, '(3f10.3,3g15.6)' ) x, y, z(j), vx, vy, vz(i,j)

      end do
    end do
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tvec_even3 ( nt, t )

!*****************************************************************************80
!
!! TVEC_EVEN3 computes an evenly spaced set of angles between 0 and 2*PI.
!
!  Discussion:
!
!    The angles begin with 0 and end with 2*PI.
!
!  Example:
!
!    NT = 4
!
!    T = ( 0, 2*PI/3, 4*PI/3 2*PI )
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
!    Input, integer NT, the number of values to compute.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles, in radians.
!
  implicit none

  integer nt

  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t(nt)

  if ( nt == 1 ) then
    t(1) = 0.0D+00
  else
    do i = 1, nt
      t(i) = real ( 2 * ( i - 1 ),  kind = 8 ) * pi &
           / real (      nt - 1,    kind = 8 )
    end do
  end if

  return
end
