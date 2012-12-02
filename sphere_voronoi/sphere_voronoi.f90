program main

!*****************************************************************************80
!
!! MAIN is the main program for SPHERE_VORONOI.
!
!  Discussion:
!
!    SPHERE_VORONOI determines the Voronoi diagram of points on the sphere in 3D.
!
!    The user supplies an XYZ file containing points on the unit sphere.
!    This routine reads that file and calls the appropriate STRIPACK
!    routines to compute and plot the Voronoi diagram. 
!
!    The plot is stored as PostScript files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) lens
  character ( len = 255 ) point_file_name

  ierror = 0

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_VORONOI'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Given a set of XYZ points on the sphere,'
  write ( *, '(a)' ) '  compute the Voronoi diagram.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( arg_num < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input XYZ file name:'
    read ( *, '(a)', iostat = ios ) point_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_VORONOI - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1
    call getarg ( iarg, point_file_name )

  end if

  call stripack_interface ( point_file_name )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_VORONOI'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine stripack_interface ( point_file_name )

!*****************************************************************************80
!
!! STRIPACK_INTERFACE calls STRIPACK routines.
!
!  Discussion:
!
!    10 June 2002: Changed this routine so that it can handle any size
!    problem, by making all arrays allocatable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) POINT_FILE_NAME, the name of the input file.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable, dimension ( : ) :: ds
  real ( kind = 8 ) elat
  real ( kind = 8 ) elon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ), allocatable, dimension ( : ) :: iwk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kt
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: lbtri
  integer ( kind = 4 ), allocatable, dimension ( : ) :: lend
  integer ( kind = 4 ), allocatable, dimension ( : ) :: list
  integer ( kind = 4 ), allocatable, dimension ( : ) :: listc
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ), allocatable, dimension ( : ) :: lptr
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: ltri
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nn
  real ( kind = 8 ) norm
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) ntemp
  logical                numbr
  integer ( kind = 4 ) nv
  real ( kind = 8 ), parameter :: pltsiz = 7.5D+00
  character ( len = * ) point_file_name
  real ( kind = 8 ), allocatable, dimension ( : ) :: rc
  integer ( kind = 4 ) side_max
  real ( kind = 8 ) vlat
  real ( kind = 8 ) vlon
  character ( len = 255 ) voronoi_plot_file_name
  character ( len = 255 ) voronoi_plot_title
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ), allocatable, dimension ( : ) :: xc
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  real ( kind = 8 ), allocatable, dimension ( : ) :: yc
  real ( kind = 8 ), allocatable, dimension ( : ) :: z
  real ( kind = 8 ), allocatable, dimension ( : ) :: zc
!
!  Count the number of lines of (X,Y,Z) data.
!
  call file_row_count ( point_file_name, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  The input file has no data.'
    return
  end if
!
!  Allocate everything.
!
  allocate ( ds(1:n) )
  allocate ( iwk(1:2*n) )
  allocate ( lbtri(1:6,1:n) )
  allocate ( lend(1:n) )
  allocate ( list(1:6*(n-2)) )
  allocate ( listc(1:6*(n-2)) )
  allocate ( lptr(1:6*(n-2)) )
  allocate ( ltri(1:9,1:2*(n-2)) )
  allocate ( rc(1:2*(n-2)) )
  allocate ( x(1:n) )
  allocate ( xc(1:2*(n-2)) )
  allocate ( y(1:n) )
  allocate ( yc(1:2*(n-2)) )
  allocate ( z(1:n) )
  allocate ( zc(1:2*(n-2)) )
!
!  Read the (X,Y,Z) data from a file.
!
  call xyz_read ( point_file_name, n, x, y, z, ierror )
!
!  Make sure the data is on the unit sphere.
!
  do i = 1, n
    norm = sqrt ( x(i)**2 + y(i)**2 + z(i)**2 )
    x(i) = x(i) / norm
    y(i) = y(i) / norm
    z(i) = z(i) / norm
  end do
!
!  Create the triangulation.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ierror )

  if ( ierror == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    stop
  end if

  if ( 0 < ierror ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Create a triangle list.
!
  call trlist ( n, list, lptr, lend, 9, nt, ltri, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRLIST.'
    stop
  end if
!
!  Construct the Voronoi diagram.
!
!  Note that the triangulation data structure is altered if NB > 0.
!
  call crlist ( n, n, x, y, z, list, lend, lptr, lnew, &
    lbtri, listc, nb, xc, yc, zc, rc, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in CRLIST.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Count the number of polygons of each size.
!
  call poly_count_max ( n, lend, lptr, listc, side_max )

  call poly_count ( n, lend, lptr, listc, side_max )
!
!  Plot the portion of the Voronoi diagram contained 
!  in the hemisphere centered at E = (ELAT,ELON), where ELAT and ELON
!  are taken to be the center of the range of
!  the nodal latitudes and longitudes.
!
  elat = 0.0D+00
  elon = 0.0D+00
  a = 90.0D+00
  numbr = ( n <= 200 )
  nt = 2 * n - 4

  voronoi_plot_file_name = point_file_name
  call file_name_ext_swap ( voronoi_plot_file_name, 'eps' )

  voronoi_plot_title = '(' // trim ( point_file_name ) // ')'

  call get_unit ( iunit )

  open ( unit = iunit, file = voronoi_plot_file_name )

  call vrplot ( iunit, pltsiz, elat, elon, a, n, x, y, z, nt, listc, &
    lptr, lend, xc, yc, zc, voronoi_plot_title, numbr, ierror )

  close ( unit = iunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Warning!'
    write ( *, '(a,i8)' ) '  VRPLOT returned error code ', ierror
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VRPLOT created the Voronoi plot file: "' // &
    trim ( voronoi_plot_file_name ) // '".'
!
!  Write the XYZ and XYZF files that define the Voronoi information.
!
  call vr_to_xyzf ( n, xc, yc, zc, listc, lptr, lend, point_file_name )
!
!  Deallocate.
!
  deallocate ( ds )
  deallocate ( iwk )
  deallocate ( lbtri )
  deallocate ( lend )
  deallocate ( list )
  deallocate ( listc )
  deallocate ( lptr )
  deallocate ( ltri )
  deallocate ( rc )
  deallocate ( x )
  deallocate ( xc )
  deallocate ( y )
  deallocate ( yc )
  deallocate ( z )
  deallocate ( zc )

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the numeric value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the numeric value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine file_name_ext_get ( file_name, i, j )

!*****************************************************************************80
!
!! FILE_NAME_EXT_GET determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILE_NAME   I  J
!
!    bob.for     4  7
!    N.B.C.D     6  7
!    Naomi.      6  6
!    Arthur     -1 -1 
!    .com        1  1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and last 
!    characters in the file extension.  
!
!    If no period occurs in FILE_NAME, then
!      I = J = -1;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s_index_last_c

  i = s_index_last_c ( file_name, '.' )

  if ( i == -1 ) then

    j = -1

  else

    j = len_trim ( file_name )

  end if

  return
end
subroutine file_name_ext_swap ( file_name, ext )

!*****************************************************************************80
!
!! FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!  Example:
!
!          Input           Output
!    ================     =========
!    FILE_NAME    EXT     FILE_NAME
!
!    bob.for      obj     bob.obj
!    bob.bob.bob  txt     bob.bob.txt
!    bob          yak     bob.yak
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME, a file name.
!    On output, the extension of the file has been changed.
!
!    Input, character ( len = * ) EXT, the extension to be used on the output
!    copy of FILE_NAME, replacing the current extension if any.
!
  implicit none

  character ( len = * ) ext
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len_max
  integer ( kind = 4 ) len_name

  len_max = len ( file_name )
  len_name = len_trim ( file_name )

  call file_name_ext_get ( file_name, i, j )

  if ( i == -1 ) then

    if ( len_max < len_name + 1 ) then
      return
    end if

    len_name = len_name + 1
    file_name(len_name:len_name) = '.'
    i = len_name + 1

  else

    i = i + 1
    file_name(i:j) = ' '

  end if

  file_name(i:) = ext

  return
end
subroutine file_row_count ( file_name, nline )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of rows in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Blank lines and comment lines, which begin with '#', are not counted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NLINE, the number of lines found in the file.
!    If the file could not be opened, then NLINE is returned as -1.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line
  integer ( kind = 4 ) nline
  logical, parameter :: verbose = .false.

  nline = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then

    nline = -1

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file:'
      write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    end if

    return

  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    nline = nline + 1

  end do

  close ( unit = iunit )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
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
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
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
subroutine poly_count ( n, lend, lptr, listc, side_max )

!*****************************************************************************80
!
!! POLY_COUNT counts the number of polygons of each size in the diagram.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi polygons.
!
!    Input, integer ( kind = 4 ) LEND(N), some kind of pointer.
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), some other kind of pointer.
!
!    Input, integer ( kind = 4 ) LISTC(6*(N-2)), some other kind of pointer.
!
!    Input, integer ( kind = 4 ) SIDE_MAX, the maximum polygonal order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) side_max

  integer ( kind = 4 ) count(0:side_max)
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kv
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) sides
  integer ( kind = 4 ) vertices

  count(0:side_max) = 0

  edges = 0
  vertices = 0

  do n0 = 1, n

    lpl = lend(n0)

    lp = lpl

    sides = 0

    do

      lp = lptr(lp)
      kv = listc(lp)

      vertices = max ( vertices, kv )
      sides = sides + 1
      edges = edges + 1

      if ( lp == lpl ) then
        exit
      end if

    end do

    if ( sides < 3 ) then
      write ( *, * ) '  Polygon ', n0, ' has ', sides, ' sides.'
    end if

    if ( side_max < sides ) then
      write ( *, * ) '  Polygon ', n0, ' has too many sides!'
    end if

    if ( sides <= 0 ) then
      count(0) = count(0) + 1
    else if ( 0 < sides .and. sides <= side_max ) then
      count(sides) = count(sides) + 1
    end if

  end do

  edges = edges / 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLY_COUNT'
  write ( *, '(a)' ) '  Number of polygons of each shape.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Faces =    ', n
  write ( *, '(a,i8)' ) '  Vertices = ', vertices
  write ( *, '(a,i8)' ) '  Edges =    ', edges
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Check Euler''s formula:'
  write ( *, '(a,i8)' ) '  F+V-E-2 =  ', n + vertices - edges - 2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Sides    Number'
  write ( *, '(a)' ) ' '

  do i = 1, side_max
    if ( count(i) /= 0 ) then
      write ( *, '(2x,i8,2x,i8)' ) i, count(i)
    end if
  end do

  return
end
subroutine poly_count_max ( n, lend, lptr, listc, side_max )

!*****************************************************************************80
!
!! POLY_COUNT_MAX determines the highest order polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi polygons.
!
!    Input, integer ( kind = 4 ) LEND(N), some kind of pointer.
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), some other kind of pointer.
!
!    Input, integer ( kind = 4 ) LISTC(6*(N-2)), some other kind of pointer.
!
!    Output, integer ( kind = 4 ) SIDE_MAX, the highest order polgon.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) side_max
  integer ( kind = 4 ) sides

  side_max = - 1

  do n0 = 1, n

    lpl = lend(n0)

    lp = lpl

    sides = 0

    do

      lp = lptr(lp)

      sides = sides + 1

      if ( lp == lpl ) then
        exit
      end if

    end do

    side_max = max ( side_max, sides )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLY_COUNT_MAX'
  write ( *, '(a,i8)' ) '  Highest order polygon is ', side_max

  return
end
function s_index_last_c ( s, c )

!*****************************************************************************80
!
!! S_INDEX_LAST_C finds the LAST occurrence of a given character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character C, the character to search for.
!
!    Output, integer ( kind = 4 ) S_INDEX_LAST_C, the index in S where C occurs
!    last, or -1 if it does not occur.
!
  implicit none

  character              c
  integer ( kind = 4 ) i
  character ( len = * )  s
  integer ( kind = 4 ) s_len
  integer ( kind = 4 ) s_index_last_c

  if ( c == ' ' ) then
    s_len = len ( s )
  else
    s_len = len_trim ( s )
  end if

  do i = s_len, 1, -1

    if ( s(i:i) == c ) then
      s_index_last_c = i
      return
    end if

  end do

  s_index_last_c = -1

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer ( kind = 4 ) part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer ( kind = 4 ) part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( ihave > 1 ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( ihave >= 6 .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. lchar+1 >= nchar ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 )  date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = * )  string
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 )  zone

  call date_and_time ( date, time, zone, values )

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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine vr_to_xyzf ( n, xc, yc, zc, listc, lptr, lend, point_file_name )

!*****************************************************************************80
!
!! VR_TO_XYZF makes an XYZF file of Voronoi diagram data.
!
!  Discussion:
!
!    Actually, this routine makes two files:
!
!    an XYZ file of the Voronoi vertices,
!    an XYZF file of the faces formed by the Voronoi vertices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(2*(N-2)), YC(2*(N-2)), ZC(2*(N-2)),
!    the coordinates  of the Voronoi vertices.
!
!    Input, integer ( kind = 4 ) LISTC(6*(N-2)), LPTR(6*(N-2)), LEND(N), 
!    information defining the triangulation, created by TRMESH.
!
!    Input, character ( len = * ) POINT_FILE_NAME, the name of the
!    file containing the point coordinates.
!
  implicit none

  integer ( kind = 4 ) n

  logical, parameter :: header = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) kv1
  integer ( kind = 4 ) kv2
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  character ( len = *  ) point_file_name
  character ( len = 40 ) string
  character ( len = 255 ) voronoi_face_file_name
  character ( len = 255 ) voronoi_vertex_file_name
  real ( kind = 8 ) xc(2*(n-2))
  real ( kind = 8 ) yc(2*(n-2))
  real ( kind = 8 ) zc(2*(n-2))
!
!  Write the Voronoi diagram vertex file.
!
  voronoi_vertex_file_name = point_file_name

  call file_name_ext_swap ( voronoi_vertex_file_name, 'voronoi.xyz' )

  call get_unit ( iunit )

  open ( unit = iunit, file = voronoi_vertex_file_name )

  if ( header ) then
    write ( iunit, '(a)' ) '# "' // trim ( voronoi_vertex_file_name ) // '".'
    write ( iunit, '(a)' ) '# created by SPHERE_VORONOI.F90.'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) '# This file contains the coordinates of points that'
    write ( iunit, '(a)' ) '# are the vertices of faces in the Voronoi diagram.'
    write ( iunit, '(a)' ) '#'
  end if

  do i = 1, 2 * ( n - 2 )
    write ( iunit, '(3f8.4)' ) xc(i), yc(i), zc(i)
  end do
  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VR_TO_XYZF:'
  write ( *, '(a)' ) '  Wrote the Voronoi XYZ file "' &
    // trim ( voronoi_vertex_file_name ) // '".'
!
!  Write the Voronoi diagram face file.
!
  voronoi_face_file_name = point_file_name

  call file_name_ext_swap ( voronoi_face_file_name, 'voronoi.xyzf' )

  call get_unit ( iunit )

  open ( unit = iunit, file = voronoi_face_file_name )

  if ( header ) then
    write ( iunit, '(a)' ) '# "' // trim ( voronoi_face_file_name ) // '".'
    write ( iunit, '(a)' ) '# created by SPHERE_VORONOI.F90.'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) '# This file contains the indices of Voronoi vertices that'
    write ( iunit, '(a)' ) '# form the faces of the Voronoi diagram.'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) '# The coordinates of these vertices are stored in'
    write ( iunit, '(a)' ) '# "' // trim ( voronoi_vertex_file_name ) // '".'
    write ( iunit, '(a)' ) '#'
  end if

  do i = 1, n

    lpl = lend(i)

    kv2 = listc(lpl) 
    write ( iunit, '(i8)', advance = 'no' ) kv2

    lp = lpl

    do

      lp = lptr(lp)
      kv1 = kv2
      kv2 = listc(lp)

      write ( iunit, '(i8)', advance = 'no' ) kv2

      if ( lp == lpl ) then
        write ( iunit, '(i8)', advance = 'yes' ) -1
        exit
      end if

    end do

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VR_TO_XYZF:'
  write ( *, '(a)' ) '  Wrote the Voronoi XYZF file "'  &
    // trim ( voronoi_face_file_name ) // '".'

  return
end
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Discussion:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!    Also, if there is a trailing comma on the word, it is stripped off.
!    This is to facilitate the reading of lists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!    On input with a fresh string, set DONE to TRUE.
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: lenc = 0
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( ilo > lenc ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( next <= lenc )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do
!
!  Ignore a trailing comma.
!
  if ( s(next-1:next-1) == ',' ) then
    word = s(ilo:next-2)
  else
    word = s(ilo:next-1)
  end if

  return
end
subroutine xyz_read ( point_file_name, n, x, y, z, ierror )

!*****************************************************************************80
!
!! XYZ_READ reads graphics information from an XYZ file.
!
!  Discussion:
!
!    Comment lines begin with '#";
!    The XYZ coordinates of a point are written on a single line.
!
!  Example:
!
!     # cube.xyz
!     #
!     0 0 0
!     0 0 1
!     0 1 0
!     0 1 1
!     1 0 0
!     1 0 1
!     1 1 0
!     1 1 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) XYZ_FILE_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) X(N), Y(N), Z(N), the point coordinates.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lchar
  character ( len = 255 ) line
  integer ( kind = 4 ) n2
  real ( kind = 8 ) temp(3)
  integer ( kind = 4 ) text_num
  character ( len = 255 ) word
  real ( kind = 8 ) x(n)
  character ( len = * ) point_file_name
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  n2 = 0
  ierror = 0
  word = ' '
  text_num = 0

  call get_unit ( iunit )

  open ( unit = iunit, file = point_file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    ierror = 1
    stop
  end if
!
!  Read a line of text from the file.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    text_num = text_num + 1
!
!  If this line begins with '#' , then it's a comment.  Read a new line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  If this line is blank, then record that information.
!
    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  This line records a node's coordinates.
!
    n2 = n2 + 1

    if ( n2 <= n ) then

      done = .true.

      do i = 1, 3

        call word_next_read ( line, word, done )

        call s_to_r8 ( word, temp(i), ierror, lchar )

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'XYZ_READ - Fatal error!'
          write ( *, '(a,i8)' ) '  S_TO_R8 returned IERROR = ', ierror
          write ( *, '(a,i8)' ) '  Reading (X,Y,Z) component ', i
          exit
        end if

      end do

      if ( ierror /= 0 ) then
        exit
      end if

      x(n2) = temp(1)
      y(n2) = temp(2)
      z(n2) = temp(3)

    end if

  end do

  close ( unit = iunit )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) 'XYZ_READ:'
  write ( *, '(a,i8,a)' ) '  Read ', text_num, ' text lines from ' &
    // trim ( point_file_name )
  write ( *, '(a,i8,a)' ) '  Read ', n2, ' sets of (X,Y,Z) coordinates.'

  return
end
