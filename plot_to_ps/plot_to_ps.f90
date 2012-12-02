program main

!*****************************************************************************80
!
!! MAIN is the main program for PLOT_TO_PS.
!
!  Discussion:
!
!    PLOT_TO_PS converts plot commands into a PostScript file.
!
!    The program can be invoked by:
!
!      plot_to_ps  plot_file_name  ps_file_name
!
!    or:
!
!      plot_to_ps  plot_file_name
!
!    in which case the output file name will be constructed by replacing
!    the extension of PLOT_FILE_NAME by ".ps", or:
!
!      plot_to_ps
!
!    in which case the user will be asked to supply both file names.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
  implicit none

  logical, parameter :: debug = .false.
  character ( len = 10 ) ext
  character ( len = 256 ) filein_name
  integer filein_unit
  character ( len = 256 ) fileout_name
  integer iarg
  integer iargc
  integer ierror
  integer ilen
  integer ios
  integer ipxfargc
  integer lens
  integer num_arg
!
!  Initialize.
!
  filein_name = ' '
  fileout_name = ' '
  ierror = 0
!
!  Get the number of command line arguments.
!
!  Old style:
!
  num_arg = iargc ( )
!
!  New style:
!
! num_arg = ipxfargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( num_arg < 1 ) then

    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT_TO_PS:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Convert simple plot descriptions into a PostScript file.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input PLOT file name:'
    read ( *, '(a)', iostat = ios ) filein_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PLOT_TO_PS - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    else if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The input file is ' // trim ( filein_name )
    end if

  else

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, filein_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, filein_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'PLOT_TO_PS - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  end if
!
!  If two command line arguments, the second one is the output file name.
!
  if ( num_arg < 2 ) then

    fileout_name = filein_name
    ext = 'ps'

    call file_name_ext_swap ( fileout_name, ext )

    if ( debug ) then
      write ( *, '(a)' ) '  The output file is ' // trim ( fileout_name )
    end if

  else

    iarg = 2
!
!  Old style:
!
    call getarg ( iarg, fileout_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, fileout_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'PLOT_TO_PS - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  end if
!
!  Open the input file.
!
  call get_unit ( filein_unit )

  open ( unit = filein_unit, file = filein_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT_TO_PS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(a)' ) trim ( filein_name )
    stop
  end if

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT_TO_PS - DEBUG:'
    write ( *, '(a)' ) '  The input file has been opened.'
  end if
!
!  Process the commands.
!
  call process ( filein_unit, fileout_name, ierror )

  close ( unit = filein_unit )

  if ( num_arg < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT_TO_PS:'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  stop
end
subroutine box_clip_line_2d ( xmin, ymin, xmax, ymax, x1, y1, x2, y2, x3, y3, &
  x4, y4, ival )

!*****************************************************************************80
!
!! BOX_CLIP_LINE_2D uses a box to clip a line segment in 2D.
!
!  Discussion:
!
!    The box is assumed to be a rectangle with sides aligned on coordinate
!    axes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, YMIN, XMAX, YMAX, the minimum and maximum
!    X and Y values, which define the box.
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the coordinates of the
!    endpoints of the line segment.
!
!    Output, real ( kind = 8 ) X3, Y3, X4, Y4, the clipped coordinates.
!
!    Output, integer IVAL:
!    -1, no part of the line segment is within the box.
!     0, no clipping was necessary.  The line segment is entirely within
!        the box.
!     1, (X1,Y1) was clipped.
!     2, (X2,Y2) was clipped.
!     3, (X1,Y1) and (X2,Y2) were clipped.
!
  implicit none

  integer ival
  logical l1
  logical l2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  l1 = .false.
  l2 = .false.

  x3 = x1
  y3 = y1
  x4 = x2
  y4 = y2
!
!  Require that XMIN <= X.
!
  if ( x3 < xmin .and. x4 < xmin ) then
    ival = -1
    return
  end if

  if ( x3 < xmin .and. xmin <= x4 ) then
    x = xmin
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x3 = x
    y3 = y
    l1 = .true.
  else if ( xmin <= x3 .and. x4 < xmin ) then
    x = xmin
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x4 = x
    y4 = y
    l2 = .true.
  end if
!
!  Require that X <= XMAX.
!
  if ( xmax < x3 .and. xmax < x4 ) then
    ival = -1
    return
  end if

  if ( xmax < x3 .and. x4 <= xmax ) then
    x = xmax
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x3 = x
    y3 = y
    l1 = .true.
  else if ( x3 <= xmax .and. xmax < x4 ) then
    x = xmax
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x4 = x
    y4 = y
    l2 = .true.
  end if
!
!  Require that YMIN <= Y.
!
  if ( y3 < ymin .and. y4 < ymin ) then
    ival = -1
    return
  end if

  if ( y3 < ymin .and. ymin <= y4 ) then
    y = ymin
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y3 = y
    x3 = x
    l1 = .true.
  else if ( ymin <= y3 .and. y4 < ymin ) then
    y = ymin
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y4 = y
    x4 = x
    l2 = .true.
  end if
!
!  Require that Y <= YMAX.
!
  if ( ymax < y3 .and. ymax < y4 ) then
    ival = -1
    return
  end if

  if ( ymax < y3 .and. y4 <= ymax ) then
    y = ymax
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y3 = y
    x3 = x
    l1 = .true.
  else if ( y3 <= ymax .and. ymax < y4 ) then
    y = ymax
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y4 = y
    x4 = x
    l2 = .true.
  end if

  ival = 0

  if ( l1 ) then
    ival = ival + 1
  end if

  if ( l2 ) then
    ival = ival + 2
  end if

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
  integer itemp

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
!    C_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical C_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
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
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine circle_points ( x0, y0, r, n, x, y )

!*****************************************************************************80
!
!! CIRCLE_POINTS returns N equally spaced points on a circle in 2D.
!
!  Discussion:
!
!    The first point is always ( X0 + R, Y0 ), and subsequent points
!    proceed counterclockwise around the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of
!    the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points
!    on the circle.
!
  implicit none

  integer n

  real ( kind = 8 ) angle
  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) x0
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y(n)

  do i = 1, n
    angle = ( 2.0D+00 * pi * real ( i - 1, kind = 8 ) ) &
      / real ( n, kind = 8 )
    x(i) = x0 + r * cos ( angle )
    y(i) = y0 + r * sin ( angle )
  end do

  return
end
function degrees_to_radians ( angle )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, an angle in degrees.
!
!    Output, real ( kind = 8 ) DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  degrees_to_radians = ( angle / 180.0D+00 ) * pi

  return
end
subroutine ellipse_points ( x0, y0, r1, r2, theta, n, x, y )

!*****************************************************************************80
!
!! ELLIPSE_POINTS returns N points on an ellipse in 2D.
!
!  Discussion:
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter of the ellipse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the ellipse.
!
!    Input, real ( kind = 8 ) R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    Input, real ( kind = 8 ) THETA, the angle that the major axis of the
!    ellipse makes with the X axis.  A value of 0.0 means that the major and
!    minor axes of the ellipse will be the X and Y coordinate axes.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points on
!    the ellipse.
!
  implicit none

  integer n

  real ( kind = 8 ) angle
  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) theta
  real ( kind = 8 ) x0
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y(n)

  do i = 1, n

    angle = ( 2.0D+00 * pi * real ( i - 1, kind = 8 ) ) &
      / real ( n, kind = 8 )

    x(i) = x0 + r1 * cos ( theta ) * cos ( angle ) &
              - r2 * sin ( theta ) * sin ( angle )

    y(i) = y0 + r1 * sin ( theta ) * cos ( angle ) &
              + r2 * cos ( theta ) * sin ( angle )

  end do

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
!    Arthur      0  0
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
!    Output, integer I, J, the indices of the first and last characters
!    in the file extension.
!
!    If no period occurs in FILE_NAME, then
!      I = J = 0;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) file_name
  integer i
  integer j
  integer s_index_last

  i = s_index_last ( file_name, '.' )

  if ( i /= 0 ) then

    j = len_trim ( file_name )

  else

    j = 0

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
  integer i
  integer j
  integer len_max
  integer len_name

  len_max = len ( file_name )
  len_name = len_trim ( file_name )

  call file_name_ext_get ( file_name, i, j )

  if ( i == 0 ) then

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
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end
function point_inside_box_2d ( x1, y1, x2, y2, x, y )

!*****************************************************************************80
!
!! POINT_INSIDE_BOX_2D determines if a point is inside a box in 2D.
!
!  Discussion:
!
!    A "box" is defined by its "left down" corner and its
!    "right up" corner, and all the points between.  It is
!    assumed that the sides of the box align with coordinate directions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the two corners of the box.
!
!    Input, real ( kind = 8 ) X, Y, the point to be checked.
!
!    Output, logical POINT_INSIDE_BOX_2D, is .TRUE. if (X,Y) is inside the
!    box, or on its boundary, and .FALSE. otherwise.
!
  implicit none

  logical point_inside_box_2d
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  if ( x1 <= x .and. x <= x2 .and. &
       y1 <= y .and. y <= y2 ) then
    point_inside_box_2d = .true.
  else
    point_inside_box_2d = .false.
  end if

  return
end
subroutine process ( filein_unit, fileout_name, ierror )

!*****************************************************************************80
!
!! PROCESS processes the user input.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer FILEIN_UNIT, the unit number associated with the input
!    file.  This will be 0 in the case of an interactive session, when
!    input comes from the user, and not a file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output
!    plot file.
!
!    Output, integer IERROR, error flag.
!    0, no error detected,
!    nonzero, an error occurred.
!
  implicit none

  integer, parameter :: maxn = 200

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_min
  real ( kind = 8 ) b
  logical, save :: debug = .false.
  real ( kind = 8 ) degrees_to_radians
  integer filein_badlines
  integer filein_linecount
  integer filein_unit
  character ( len = * ) fileout_name
  integer fileout_unit
  real ( kind = 8 ) g
  real ( kind = 8 ), save :: hist_fat = 1.0D+00
  integer ierror
  integer ios
  integer lchar
  integer lenc
  character ( len = 80 ) line
  integer n
  integer n1
  integer n2
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  logical s_eqi
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  character ( len = 80 ) word
  real ( kind = 8 ) :: xmax = 1.0D+00
  real ( kind = 8 ) :: xmin = 0.0D+00
  real ( kind = 8 ) xvec(maxn)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) :: ymax = 1.0D+00
  real ( kind = 8 ) :: ymin = 0.0D+00
  real ( kind = 8 ) yvec(maxn)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  ierror = 0

  call get_unit ( fileout_unit )
!
!  Read the next line of input.
!
  filein_badlines = 0
  filein_linecount = 0

  do

    if ( 0 < filein_unit ) then
      read ( filein_unit, '(a)', iostat = ios ) line
    else
      read ( *, '(a)', iostat = ios ) line
    end if

    if ( debug ) then
      write ( *, '(a)' ) trim ( line )
    end if

    if ( ios /= 0 ) then
      exit
    end if

    filein_linecount = filein_linecount + 1

    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) <= 0 ) then
      cycle
    end if
!
!  Read the first word of the input line.
!
    call word_extract ( line, word )
!
!  Determine what command the first word represents, read its arguments,
!  and write the appropriate instructions to the file.
!
    if ( s_eqi ( word, 'ARC' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, theta1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, theta2, ierror, lchar )

      call ps_circle_arc ( x1, y1, r1, theta1, theta2 )

    else if ( s_eqi ( word, 'ARROW' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, x2, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y2, ierror, lchar )

      call ps_arrow ( x1, y1, x2, y2 )

    else if ( s_eqi ( word, 'CIRCLE' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call ps_circle ( x1, y1, r1 )

    else if ( s_eqi ( word, 'CIRCLE_FILL' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call ps_circle_fill ( x1, y1, r1 )

    else if ( s_eqi ( word, 'DEBUG' ) ) then

      debug = .not. debug

      if ( debug ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PROCESS:'
        write ( *, '(a)' ) '  Debugging turned on.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PROCESS:'
        write ( *, '(a)' ) '  Debugging turned off.'
      end if

    else if ( s_eqi ( word, 'ELLIPSE' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r2, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, theta1, ierror, lchar )

      n = maxn

      call ellipse_points ( x1, y1, r1, r2, theta1, n, xvec, yvec )

      call ps_line_closed ( n, xvec, yvec )

    else if ( s_eqi ( word, 'ELLIPSE_FILL' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r2, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, theta1, ierror, lchar )

      n = maxn

      call ellipse_points ( x1, y1, r1, r2, theta1, n, xvec, yvec )

      call ps_polygon_fill ( n, xvec, yvec )

    else if ( s_eqi ( word, 'ENDFILE' ) .or. s_eqi ( word, 'END_FILE' ) ) then

      call ps_file_tail

      call ps_file_close ( fileout_unit )

    else if ( s_eqi ( word, 'ENDPAGE' ) .or. s_eqi ( word, 'END_PAGE' ) ) then

      call ps_page_tail

    else if ( s_eqi ( word, 'FILE' ) ) then

      if ( fileout_name == ' ' ) then

        call word_extract ( line, word )

        if ( word /= ' ' ) then
          fileout_name = word
        else
          fileout_name = 'plot_to_ps.ps'
        end if

      end if

      call get_unit ( fileout_unit )

      call ps_file_open ( fileout_name, fileout_unit, ierror )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PROCESS'
        write ( *, '(a,i6)' ) '  File creation error ', ierror
        return
      end if

    else if ( s_eqi ( word, 'FILL_GRAY' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call ps_fill_gray ( x1 )

    else if ( s_eqi ( word, 'FILL_RGB' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, r, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, g, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, b, ierror, lchar )

      call ps_color_fill_set ( r, g, b )

    else if ( s_eqi ( word, 'FONT_SIZE' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call ps_font_size ( x1 )

    else if ( s_eqi ( word, 'GRID' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, x2, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y2, ierror, lchar )
      call word_extract ( line, word )
      call s_to_i4 ( word, n1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_i4 ( word, n2, ierror, lchar )

      call ps_grid_cartesian ( x1, x2, n1, y1, y2, n2 )

    else if ( s_eqi ( word, 'HIST_FAT' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      hist_fat = x1

    else if ( s_eqi ( word, 'HISTOGRAM' ) ) then

      call process_histogram ( filein_linecount, filein_unit, hist_fat, &
        ierror )

    else if ( s_eqi ( word, 'LABEL' ) ) then

      call ps_label ( line )

    else if ( s_eqi ( word, 'LABEL_SLANT' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, angle, ierror, lchar )

      call ps_label_slant ( line, angle )

    else if ( s_eqi ( word, 'LINE' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, x2, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y2, ierror, lchar )

      call ps_moveto ( x1, y1 )
      call ps_lineto ( x2, y2 )

    else if ( s_eqi ( word, 'LINE_GRAY' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, r, ierror, lchar )

      g = r
      b = r

      call ps_color_line_set ( r, g, b )

    else if ( s_eqi ( word, 'LINE_RGB' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, r, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, g, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, b, ierror, lchar )

      call ps_color_line_set ( r, g, b )

    else if ( s_eqi ( word, 'LINE_WIDTH' ) .or. &
              s_eqi ( word, 'LINEWIDTH' ) .or. &
              s_eqi ( word, 'LINE_THICK' ) .or. &
              s_eqi ( word, 'LINETHICK' ) ) then

      call word_extract ( line, word )
      call s_to_i4 ( word, n1, ierror, lchar )

      call ps_line_width ( n1 )

    else if ( s_eqi ( word, 'LINETO' ) .or. &
              s_eqi ( word, 'DRAWTO' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call ps_lineto ( x1, y1 )

    else if ( s_eqi ( word, 'MOVE' ) .or. &
              s_eqi ( word, 'MOVETO' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call ps_moveto ( x1, y1 )

    else if ( s_eqi ( word, 'PAGE' ) ) then

      call ps_page_head ( xmin, ymin, xmax, ymax )

    else if ( s_eqi ( word, 'POINT' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call ps_mark_circle ( x1, y1 )

    else if ( s_eqi ( word, 'POINT_LIST' ) ) then

      do

        if ( 0  < filein_unit ) then
          read ( filein_unit, '(a)', iostat = ios ) line
        else
          read ( *, '(a)', iostat = ios ) line
        end if

        if ( ios /= 0 ) then
          exit
        end if

        filein_linecount = filein_linecount + 1

        if ( line(1:1) == '#' ) then
          cycle
        end if

        if ( len_trim ( line ) <= 0 ) then
          cycle
        end if

        call word_extract ( line, word )

        if ( s_eqi ( word(1:3), 'END' ) ) then
          exit
        end if

        call s_to_r8 ( word, x1, ierror, lchar )
        call word_extract ( line, word )
        call s_to_r8 ( word, y1, ierror, lchar )

        call ps_mark_circle ( x1, y1 )

      end do

    else if ( s_eqi ( word, 'POLYGON' ) ) then

      call process_polygon ( filein_linecount, filein_unit, ierror, MAXN, &
        xvec, yvec )

    else if ( s_eqi ( word, 'POLYGON_FILL' ) ) then

      call process_polygon_fill ( filein_linecount, filein_unit, ierror, MAXN, &
        xvec, yvec )

    else if ( s_eqi ( word, 'RADIUS' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, angle, ierror, lchar )

      call ps_moveto ( x1, y1 )

      x2 = x1 + r1 * cos ( degrees_to_radians ( angle ) )
      y2 = y1 + r1 * sin ( degrees_to_radians ( angle ) )

      call ps_lineto ( x2, y2 )

    else if ( s_eqi ( word, 'SECTOR' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, angle_min, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, angle_max, ierror, lchar )

      call ps_sector ( x1, y1, r1, angle_min, angle_max )

    else if ( s_eqi ( word, 'SECTOR_FILL' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, angle_min, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, angle_max, ierror, lchar )

      call ps_sector_fill ( x1, y1, r1, angle_min, angle_max )

    else if ( s_eqi ( word, 'SPACE' ) ) then

      if ( debug ) then
        write ( *, '(a)' ) 'Begin processing the SPACE command.'
      end if

      call word_extract ( line, word )
      call s_to_r8 ( word, xmin, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, ymin, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, xmax, ierror, lchar )
      call word_extract ( line, word )
      call s_to_r8 ( word, ymax, ierror, lchar )

      call ps_file_head ( fileout_name )

     else if ( s_eqi ( word, 'SQUARE' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call ps_square ( x1, y1, r1 )

    else if ( s_eqi ( word, 'SQUARE_FILL' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call ps_square_fill ( x1, y1, r1 )

    else if ( s_eqi ( word, 'STAR' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call ps_star ( x1, y1, r1 )

    else if ( s_eqi ( word, 'STAR_CIRCLE' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call ps_star_circle ( x1, y1, r1 )

    else if ( s_eqi ( word, 'STAR_DISK' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, r1, ierror, lchar )

      call ps_star_disk ( x1, y1, r1 )

     else if ( s_eqi ( word, 'TRIANGLE' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, x2, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y2, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, x3, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y3, ierror, lchar )

      call ps_triangle ( x1, y1, x2, y2, x3, y3 )

     else if ( s_eqi ( word, 'TRIANGLE_FILL' ) ) then

      call word_extract ( line, word )
      call s_to_r8 ( word, x1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y1, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, x2, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y2, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, x3, ierror, lchar )

      call word_extract ( line, word )
      call s_to_r8 ( word, y3, ierror, lchar )

      call ps_triangle_fill ( x1, y1, x2, y2, x3, y3 )

    else

      if ( .true. ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PROCESS - DEBUG: Bad line:'
        write ( *, '(a)' ) trim ( line )
      end if

      filein_badlines = filein_badlines + 1

    end if

  end do

  if ( 0 < filein_badlines ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PROCESS - Warning:'
    write ( *, '(a,i6)' ) '  Number of input lines:     ', filein_linecount
    write ( *, '(a,i6)' ) '  Number of BAD input lines: ', filein_badlines
  end if

  return
end
subroutine process_histogram ( filein_linecount, filein_unit, hist_fat, &
  ierror )

!*****************************************************************************80
!
!! PROCESS_HISTOGRAM processes the input defining a histogram.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input/output, integer FILEIN_LINECOUNT, the number of lines of
!    input that were read.
!
!    Input, integer FILEIN_UNIT, the unit number associated with the input
!    file.  This will be 0 in the case of an interactive session, when
!    input comes from the user, and not a file.
!
!    Input, real ( kind = 8 ) HIST_FAT, the width of each histogram bar.
!
!    Output, integer IERROR, error flag.
!    0, no error detected,
!    nonzero, an error occurred.
!
  implicit none

  integer, parameter :: maxpoly = 4

  integer filein_linecount
  integer filein_unit
  real ( kind = 8 ) hist_fat
  integer ierror
  integer ios
  integer lchar
  integer lenc
  character ( len = 80 ) line
  integer npoly
  logical s_eqi
  character ( len = 80 ) word
  real ( kind = 8 ) xvec(maxpoly)
  real ( kind = 8 ) x1
  real ( kind = 8 ) yvec(maxpoly)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y1_old

  ierror = 0
  y1_old = 0.0D+00
!
!  Read the next line of input.
!
  do

    if ( 0 < filein_unit ) then
      read ( filein_unit, '(a)', iostat = ios ) line
    else
      read ( *, '(a)', iostat = ios ) line
    end if

    if ( ios /= 0 ) then
      exit
    end if

    filein_linecount = filein_linecount + 1

    if ( line(1:1) == '#' ) then
      cycle
    end if

    lenc = len_trim ( line )
    if ( len_trim ( line ) <= 0 ) then
      cycle
    end if
!
!  Read the first word of the input line.
!
    call word_extract ( line, word )

    if ( s_eqi ( word(1:3), 'END' ) ) then
!
!  If fill color reversed, reverse it again.
!
      return

    end if

    call s_to_r8 ( word, x1, ierror, lchar )
    call word_extract ( line, word )
    call s_to_r8 ( word, y1, ierror, lchar )

    if ( y1 /= 0.0D+00 ) then

      xvec(1) = x1 - 0.5D+00 * hist_fat
      yvec(1) = 0.0D+00

      xvec(2) = x1 - 0.5D+00 * hist_fat
      yvec(2) = y1

      xvec(3) = x1 + 0.5D+00 * hist_fat
      yvec(3) = y1

      xvec(4) = x1 + 0.5D+00 * hist_fat
      yvec(4) = 0.0D+00

      if ( y1_old <= 0.0D+00 .and. 0.0D+00 < y1 ) then
!
!...set fill color
!
      else if ( 0.0D+00 <= y1_old .and. y1 < 0.0D+00 ) then
!
!..set fill color
!
      end if

      npoly = 4
      call ps_polygon_fill ( npoly, xvec, yvec )

      call ps_line_closed ( npoly, xvec, yvec )

      y1_old = y1

    end if

  end do

  ierror = 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROCESS_HISTOGRAM - Fatal error!'
  write ( *, '(a)' ) '  Unexpected end of input.'

  return
end
subroutine process_polygon ( filein_linecount, filein_unit, ierror, MAXN, &
  xvec, yvec )

!*****************************************************************************80
!
!! PROCESS_POLYGON processes the input defining an outlined polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input/output, integer FILEIN_LINECOUNT, the number of lines of
!    input that were read.
!
!    Input, integer FILEIN_UNIT, the unit number associated with the input
!    file.  This will be 0 in the case of an interactive session, when
!    input comes from the user, and not a file.
!
!    Output, integer IERROR, error flag.
!    0, no error detected,
!    nonzero, an error occurred.
!
  implicit none

  integer MAXN

  integer filein_linecount
  integer filein_unit
  integer ierror
  integer ios
  integer lchar
  integer lenc
  character ( len = 80 ) line
  integer npoly
  logical s_eqi
  character ( len = 80 ) word
  real ( kind = 8 ) xvec(MAXN)
  real ( kind = 8 ) x1
  real ( kind = 8 ) yvec(MAXN)
  real ( kind = 8 ) y1

  ierror = 0
  npoly = 0
!
!  Read the next line of input.
!
  do

    if ( 0 < filein_unit ) then
      read ( filein_unit, '(a)', iostat = ios ) line
    else
      read ( *, '(a)', iostat = ios ) line
    end if

    if ( ios /= 0 ) then
      exit
    end if

    filein_linecount = filein_linecount + 1

    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) <= 0 ) then
      cycle
    end if
!
!  Read the first word of the input line.
!
    call word_extract ( line, word )

    if ( s_eqi ( word(1:3), 'END' ) ) then

      npoly = min ( npoly, MAXN )

      if ( 3 <= npoly ) then
        call ps_line_closed ( npoly, xvec, yvec )
      end if

      return

    end if

    call s_to_r8 ( word, x1, ierror, lchar )
    call word_extract ( line, word )
    call s_to_r8 ( word, y1, ierror, lchar )

    npoly = npoly + 1
    if ( npoly <= MAXN ) then
      xvec(npoly) = x1
      yvec(npoly) = y1
    end if

  end do

  ierror = 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROCESS_POLYGON - Fatal error!'
  write ( *, '(a)' ) '  Unexpected end of input.'

  return
end
subroutine process_polygon_fill ( filein_linecount, filein_unit, ierror, MAXN, &
  xvec, yvec )

!*****************************************************************************80
!
!! PROCESS_POLYGON_FILL processes the input defining a filled polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input/output, integer FILEIN_LINECOUNT, the number of lines of
!    input that were read.
!
!    Input, integer FILEIN_UNIT, the unit number associated with the input
!    file.  This will be 0 in the case of an interactive session, when
!    input comes from the user, and not a file.
!
!    Output, integer IERROR, error flag.
!    0, no error detected,
!    nonzero, an error occurred.
!
  implicit none

  integer MAXN

  integer filein_linecount
  integer filein_unit
  integer ierror
  integer ios
  integer lchar
  integer lenc
  character ( len = 80 ) line
  integer npoly
  logical s_eqi
  character ( len = 80 ) word
  real ( kind = 8 ) xvec(MAXN)
  real ( kind = 8 ) x1
  real ( kind = 8 ) yvec(MAXN)
  real ( kind = 8 ) y1

  ierror = 0
  npoly = 0
!
!  Read the next line of input.
!
  do

    if ( 0 < filein_unit ) then
      read ( filein_unit, '(a)', iostat = ios ) line
    else
      read ( *, '(a)', iostat = ios ) line
    end if

    if ( ios /= 0 ) then
      exit
    end if

    filein_linecount = filein_linecount + 1

    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) <= 0 ) then
      cycle
    end if
!
!  Read the first word of the input line.
!
    call word_extract ( line, word )

    if ( s_eqi ( word(1:3), 'END' ) ) then

      npoly = min ( npoly, MAXN )

      if ( 3 <= npoly ) then
        call ps_polygon_fill ( npoly, xvec, yvec )
      end if

      return

    end if

    call s_to_r8 ( word, x1, ierror, lchar )
    call word_extract ( line, word )
    call s_to_r8 ( word, y1, ierror, lchar )

    npoly = npoly + 1
    if ( npoly <= MAXN ) then
      xvec(npoly) = x1
      yvec(npoly) = y1
    end if

  end do

  ierror = 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROCESS_POLYGON_FILL - Fatal error!'
  write ( *, '(a)' ) '  Unexpected end of input.'

  return
end
subroutine ps_arrow ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! PS_ARROW draws an arrow from (X1,Y1) to (X2,Y2).
!
!  Discussion:
!
!    The current point is set to (X2,Y2).
!
!    This routine will clip the line, if necessary, so that the line
!    drawn is entirely within the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, the starting point of the arrow.
!
!    Input, real ( kind = 8 ) X2, Y2, the ending point of the arrow.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha2
  real ( kind = 8 ) alpha3
  real ( kind = 8 ) frac
  integer ival
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  integer state
  integer unit
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5
  real ( kind = 8 ) y6
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  if ( x1 == x2 .and. y1 == y2 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_ARROW - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )
!
!  Clip the line.
!
  call box_clip_line_2d ( xmin, ymin, xmax, ymax, x1, y1, x2, y2, x3, y3, &
    x4, y4, ival )

  if ( ival < 0 ) then
    return
  end if

  r = sqrt ( ( x4 - x3 )**2 + ( y4 - y3 )**2 )

  if ( r == 0.0D+00 ) then
    return
  end if
!
!  FRAC controls the size of the arrow head.  It's specified
!  as a proportion of the length of the line.
!
  frac = 0.1D+00

  r2 = sqrt ( frac**2 + ( 1.0D+00 - frac )**2 ) * r

  alpha2 = atan2 ( y4 - y3, x4 - x3 )
  alpha3 = atan2 ( frac, 1.0D+00 - frac )

  x5 = x3 + r2 * cos ( alpha2 - alpha3 )
  y5 = y3 + r2 * sin ( alpha2 - alpha3 )

  x6 = x3 + r2 * cos ( alpha2 + alpha3 )
  y6 = y3 + r2 * sin ( alpha2 + alpha3 )
!
!  Draw line.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( x3 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y3 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x4 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y4 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw arrow head.
!
  px = plotxmin2 + nint ( alpha * ( x4 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y4 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x5 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y5 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'

  px = plotxmin2 + nint ( alpha * ( x4 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y4 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x6 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y6 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto stroke'

  call ps_setting_real ( 'SET', 'XCUR', x2 )
  call ps_setting_real ( 'SET', 'YCUR', y2 )

  return
end
subroutine ps_circle ( x0, y0, r )

!*****************************************************************************80
!
!! PS_CIRCLE draws a circle.
!
!  Discussion:
!
!    As a side effect, the current point is set to the center of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
  implicit none

  real ( kind = 8 ) alpha
  integer, parameter :: angle_max = 360
  integer, parameter :: angle_min = 0
  integer plotxmin2
  integer plotymin2
  integer pr
  integer pxcen
  integer pycen
  real ( kind = 8 ) r
  integer state
  integer unit
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_CIRCLE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )

  pxcen = plotxmin2 + nint ( alpha * ( x0 - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y0 - ymin ) )
  pr = nint ( alpha * r )

  write ( unit, '(a)' ) 'newpath'
  write ( unit, '(5i6,a)' ) pxcen, pycen, pr, angle_min, angle_max, ' arc'
!
!  Draw the circle.
!
  write ( unit, '(a)' ) 'closepath stroke'

  call ps_setting_real ( 'SET', 'XCUR', x0 )
  call ps_setting_real ( 'SET', 'YCUR', y0 )

  return
end
subroutine ps_circle_arc ( x0, y0, r, theta1, theta2 )

!*****************************************************************************80
!
!! PS_CIRCLE_ARC draws a circular arc.
!
!  Discussion:
!
!    As a side effect, the current point is set to the center of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angular coordinates
!    of the first and last points on the circular arc to be drawn.  These
!    should be ordered in counter-clockwise order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer angle_max
  integer angle_min
  integer plotxmin2
  integer plotymin2
  integer pr
  integer pxcen
  integer pycen
  real ( kind = 8 ) r
  integer state
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  integer unit
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_CIRCLE_ARC - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )

  pxcen = plotxmin2 + nint ( alpha * ( x0 - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y0 - ymin ) )
  pr = nint ( alpha * r )
  angle_min = nint ( theta1 )
  angle_max = nint ( theta2 )

  write ( unit, '(a)' ) 'newpath'
  write ( unit, '(5i6,a)' ) pxcen, pycen, pr, angle_min, angle_max, ' arc'
!
!  Draw.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x0 )
  call ps_setting_real ( 'SET', 'YCUR', y0 )

  return
end
subroutine ps_circle_fill ( x0, y0, r )

!*****************************************************************************80
!
!! PS_CIRCLE_FILL draws a filled circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the disk.
!
!    Input, real ( kind = 8 ) R, the radius of the disk.
!
  implicit none

  integer, parameter :: n = 32

  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  call circle_points ( x0, y0, r, n, x, y )

  call ps_polygon_fill ( n, x, y )

  return
end
subroutine ps_color_fill_set ( r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_FILL_SET sets the fill color.
!
!  Discussion:
!
!    By calling this routine, you guarantee that a check will be made
!    of the current fill color.  If the current and new fill colors are
!    the same, then we skip the extraneous action of setting the color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new fill color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  integer state
  integer unit
!
!  Check the state.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_FILL_SET - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  A PostScript state of 1 or more is required.'
    return
  end if
!
!  Get the current colors.
!
  call ps_setting_real ( 'GET', 'FILL_RED', r_old )
  call ps_setting_real ( 'GET', 'FILL_GREEN', g_old )
  call ps_setting_real ( 'GET', 'FILL_BLUE', b_old )
!
!  If any color has changed, we need to reset them.
!
  if ( r_old /= r .or. g_old /= g .or. b_old /= b ) then

    call ps_setting_int ( 'GET', 'UNIT', unit )

    call ps_comment ( 'Set RGB line color.' )

    write ( unit, '(3f7.4,a)' ) r, g, b, ' setrgbcolor'

    call ps_setting_real ( 'SET', 'FILL_RED', r )
    call ps_setting_real ( 'SET', 'FILL_GREEN', g )
    call ps_setting_real ( 'SET', 'FILL_BLUE', b )

  end if

  return
end
subroutine ps_color_line ( action, r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_LINE handles the line color.
!
!  Discussion:
!
!    By calling this routine, you can temporarily set the line color,
!    draw some lines, and then restore it to whatever it was.
!
!    An earlier version of this routine did not use the SAVE command for
!    the stack arrrays, meaning the stored data was lost.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action.
!    'SET', set the line color to RGB.
!    'GET', set RGB to the current line color.
!    'PUSH', push a value onto the RGB stack.
!    'POP', pop the RGB stack.
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new line color.
!
  implicit none

  integer, parameter :: nstack = 10

  character ( len = * ) action
  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ), save, dimension ( nstack) :: b_stack
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ), save, dimension ( nstack) :: g_stack
  integer, save :: istack = 0
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  real ( kind = 8 ), save, dimension ( nstack) :: r_stack
  logical s_eqi

  if ( s_eqi ( action, 'SET' ) ) then

    call ps_color_line_set ( r, g, b )

  else if ( s_eqi ( action, 'GET' ) ) then

    call ps_setting_real ( 'GET', 'LINE_RED', r )
    call ps_setting_real ( 'GET', 'LINE_GREEN', g )
    call ps_setting_real ( 'GET', 'LINE_BLUE', b )

  else if ( s_eqi ( action, 'POP' ) ) then

    if ( 0 < istack ) then
      r = r_stack(istack)
      g = g_stack(istack)
      b = b_stack(istack)
      istack = istack - 1
    end if

    call ps_color_line_set ( r, g, b )

  else if ( s_eqi ( action, 'PUSH' ) ) then

    call ps_setting_real ( 'GET', 'LINE_RED', r_old )
    call ps_setting_real ( 'GET', 'LINE_GREEN', g_old )
    call ps_setting_real ( 'GET', 'LINE_BLUE', b_old )

    if ( istack <= nstack ) then
      istack = istack + 1
      r_stack(istack) = r_old
      g_stack(istack) = g_old
      b_stack(istack) = b_old
    end if

    call ps_color_line_set ( r, g, b )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_LINE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected ACTION.'
    stop

  end if

  return
end
subroutine ps_color_line_set ( r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_LINE_SET sets the line color.
!
!  Discussion:
!
!    By calling this routine, you guarantee that a check will be made
!    of the current line color.  If the current and new line colors are
!    the same, then we skip the extraneous action of setting the color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new line color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  integer state
  integer unit
!
!  Check the state.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_LINE_SET - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  A PostScript state of at least 1 is required.'
    return
  end if
!
!  Get the current colors.
!
  call ps_setting_real ( 'GET', 'LINE_RED', r_old )
  call ps_setting_real ( 'GET', 'LINE_GREEN', g_old )
  call ps_setting_real ( 'GET', 'LINE_BLUE', b_old )
!
!  If any color has changed, we need to reset them.
!
  if ( r_old /= r .or. g_old /= g .or. b_old /= b ) then

    call ps_setting_int ( 'GET', 'UNIT', unit )

    call ps_comment ( 'Set RGB line color.' )

    write ( unit, '(3f7.4,a)' ) r, g, b, ' setrgbcolor'

    call ps_setting_real ( 'SET', 'LINE_RED', r )
    call ps_setting_real ( 'SET', 'LINE_GREEN', g )
    call ps_setting_real ( 'SET', 'LINE_BLUE', b )

  end if

  return
end
subroutine ps_comment ( string )

!*****************************************************************************80
!
!! PS_COMMENT inserts a comment into the PostScript file.
!
!  Discussion:
!
!    A comment begins with a percent sign in column 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the comment.
!
  implicit none

  character ( len = * ) string
  integer unit
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Write the comment.
!
  if ( len_trim ( string ) == 0 ) then
    write ( unit, '(a)' ) '%'
  else
    write ( unit, '(a)' ) '%'
    write ( unit, '(a2,a)' ) '% ', trim ( string )
    write ( unit, '(a)' ) '%'
  end if

  return
end
subroutine ps_default ( )

!*****************************************************************************80
!
!! PS_DEFAULT sets the internal settings to their default values
!
!  Discussion:
!
!    Certain variables are not reset, including the number of pages,
!    the unit number, the internal state, and variables relating to
!    the size and shape of the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  real ( kind = 8 ) fill_blue
  real ( kind = 8 ) fill_green
  real ( kind = 8 ) fill_red
  real ( kind = 8 ) font_size
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer line_width
  integer marker_size

  line_width = 1
  marker_size = 5

  call ps_setting_int ( 'SET', 'LINE_WIDTH', line_width )
  call ps_setting_int ( 'SET', 'MARKER_SIZE', marker_size )

  fill_blue = 0.7D+00
  fill_green = 0.7D+00
  fill_red = 0.7D+00
  font_size = 0.1D+00
  line_blue = 0.0D+00
  line_green = 0.0D+00
  line_red = 0.0D+00

  call ps_setting_real ( 'SET', 'FILL_BLUE', fill_blue )
  call ps_setting_real ( 'SET', 'FILL_GREEN', fill_green )
  call ps_setting_real ( 'SET', 'FILL_RED', fill_red )
  call ps_setting_real ( 'SET', 'FONT_SIZE', font_size )
  call ps_setting_real ( 'SET', 'LINE_BLUE', line_blue )
  call ps_setting_real ( 'SET', 'LINE_GREEN', line_green )
  call ps_setting_real ( 'SET', 'LINE_RED', line_red )

  return
end
subroutine ps_file_close ( unit )

!*****************************************************************************80
!
!! PS_FILE_CLOSE closes a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer UNIT, the FORTRAN unit to which output was written.
!
  implicit none

  integer state
  integer unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 .or. 4 < state ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_CLOSE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1, 2, 3 or 4 is required.'
    return
  end if

  close ( unit = unit )

  state = 0
  call ps_setting_int ( 'SET', 'STATE', state )

  unit = 0
  call ps_setting_int ( 'SET', 'UNIT', unit )

  return
end
subroutine ps_file_head ( file_name )

!*****************************************************************************80
!
!! PS_FILE_HEAD writes header information to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
  implicit none

  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer margin
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotymax
  integer plotymin
  integer state
  integer unit
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yvec(4)
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1 is required.'
    return
  end if
!
!  Initialization
!
  call ps_default ( )
!
!  Compute the scale factor.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin
!
!  These AD HOC values trim off the unused space.
!
! plotxmin = 36
! plotxmax = 576
! plotymin = 126
! plotymax = 666
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call date_and_time ( date )
!
!  Write the prolog.
!
  write ( unit, '(a)' )     '%!PS-Adobe-1.0'
  write ( unit, '(a)' )     '%%Creator: ps_write.f90'
  write ( unit, '(a)' )     '%%Title: ' // trim ( file_name )
  write ( unit, '(a)' )     '%%CreationDate: ' // trim ( date )
  write ( unit, '(a)' )     '%%Pages: (atend)'
  write ( unit, '(a,4i6)' ) '%%BoundingBox:', plotxmin, plotymin, plotxmax, &
    plotymax
  write ( unit, '(a)' )     '%%Document-Fonts: Times-Roman'
  write ( unit, '(a)' )     '%%LanguageLevel: 1'
  write ( unit, '(a)' )     '%%EndComments'
  write ( unit, '(a)' )     '%%BeginProlog'
  write ( unit, '(a)' )     '/inch {72 mul} def'
  write ( unit, '(a)' )     '%%EndProlog'
!
!  Set the font.
!
  call ps_comment ( 'Set the font:' )

  write ( unit, '(a)' ) '/Times-Roman findfont'
  write ( unit, '(a)' ) '1.00 inch scalefont'
  write ( unit, '(a)' ) 'setfont'
!
!  Set the line color to black.
!
  line_red = 0.0D+00
  line_green = 0.0D+00
  line_blue = 0.0D+00

  call ps_color_line ( 'SET', line_red, line_green, line_blue )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_file_open ( file_name, unit, ierror )

!*****************************************************************************80
!
!! PS_FILE_OPEN opens a new version of a PostScript file with a given name.
!
!  Discussion:
!
!    If a file of the given name already exists, it is deleted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer UNIT, the FORTRAN unit to which output should
!    be written.
!
!    Input, character ( len = 80 ) FILE_NAME, the name of the output file.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    nonzero, the file could not be created.
!
  implicit none

  character ( len = * ) file_name
  integer ierror
  integer ios
  integer state
  integer unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_OPEN - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 0 is required.'
    write ( *, '(a)' ) '  Call PS_FILE_CLOSE first!'
    return
  end if

  ierror = 0
!
!  Now create a new empty file of the given name.
!
  open ( unit = unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    return
  end if

  state = 1
  call ps_setting_int ( 'SET', 'STATE', state )

  call ps_setting_int ( 'SET', 'UNIT', unit )

  return
end
subroutine ps_file_tail ( )

!*****************************************************************************80
!
!! PS_FILE_TAIL writes trailer information to a PostScript file.
!
!  Discussion:
!
!    Looks like that penultimate 'end' line is not wanted, so
!    I commented it out.
!\!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer num_pages
  integer state
  integer unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  A page was open.  It is being forced closed.'
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Retrieve the number of pages.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )
!
!  Write the epilog.
!
  write ( unit, '(a)' ) '%%Trailer'
  write ( unit, '(a,i6)' ) '%%Pages: ', num_pages
! write ( unit, '(a)' ) 'end'
  write ( unit, '(a)' ) '%%EOF'
!
!  Zero out the number of pages.
!
  num_pages = 0

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )
!
!  Reset the state.
!
  state = 4

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_fill_gray ( fill_gray )

!*****************************************************************************80
!
!! PS_FILL_GRAY sets the gray fill for polygons.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FILL_GRAY, the gray fill used to fill polygons.
!    0.0 is black, 1.0 is white, and values in between represent
!    shades of gray.
!
  implicit none

  real ( kind = 8 ) fill_gray
  integer unit

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(f8.4,a)' ) fill_gray, ' setgray'

  call ps_setting_real ( 'SET', 'FILL_BLUE', fill_gray )
  call ps_setting_real ( 'SET', 'FILL_GREEN', fill_gray )
  call ps_setting_real ( 'SET', 'FILL_RED', fill_gray )

  return
end
subroutine ps_font_size ( font_size )

!*****************************************************************************80
!
!! PS_FONT_SIZE sets the font size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FONT_SIZE, the font size, in inches.
!
  implicit none

  real ( kind = 8 ) font_size
  integer state
  integer unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 2 .and. state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FONT_SIZE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 or 3 is required.'
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) '/Times-Roman findfont'
  write ( unit, '(f8.3, a)' ) font_size, ' inch scalefont'
  write ( unit, '(a)' ) 'setfont'

  call ps_setting_real ( 'SET', 'FONT_SIZE', font_size )

  return
end
subroutine ps_grid_cartesian ( xmin, xmax, nx, ymin, ymax, ny )

!*****************************************************************************80
!
!! PS_GRID_CARTESIAN draws a cartesian grid.
!
!  Discussion:
!
!    The current point is not modified by this call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the minimum and maximum values
!    at which X grid lines should be drawn.
!
!    Input, integer NX, the number of X grid lines.
!    If NX is not positive, no X grid lines are drawn.
!    If NX is 1, a single grid line is drawn midway.
!
!    Input, real ( kind = 8 ) YMIN, YMAX, the minimum and maximum values
!    at which Y grid lines should be drawn.
!
!    Input, integer NY, the number of Y grid lines.
!    If NY is not positive, no Y grid lines are drawn.
!    If NY is 1, a single grid line is drawn midway.
!
  implicit none

  real ( kind = 8 ) alpha
  integer i
  integer nx
  integer ny
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
!
!  At least one of NX and NY must be positive.
!
  if ( nx < 1 .and. ny < 1 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_GRID_CARTESIAN - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Get settings.
!
  alpha = 0.0D+00
  xmin2 = 0.0D+00
  ymin2 = 0.0D+00

  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin2 )
  call ps_setting_real ( 'GET', 'YMIN', ymin2 )
!
!  Draw the vertical (X) grid lines.
!
  do i = 1, nx

    if ( 1 < nx ) then

      x = ( real ( nx - i,     kind = 8 ) * xmin   &
          + real (      i - 1, kind = 8 ) * xmax ) &
          / real ( nx     - 1, kind = 8 )

    else if ( nx == 1 ) then

      x = 0.5D+00 * ( xmin + xmax )

    end if

    px = plotxmin2 + nint ( alpha * ( x - xmin2 ) )

    write ( unit, '(a)' ) 'newpath'

    py = plotymin2 + nint ( alpha * ( ymin - ymin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    py = plotymin2 + nint ( alpha * ( ymax - ymin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do
!
!  Draw the horizontal (Y) grid lines.
!
  do i = 1, ny

    if ( 1 < ny ) then

      y = ( real ( ny - i,     kind = 8 ) * ymin   &
          + real (      i - 1, kind = 8 ) * ymax ) &
          / real ( ny     - 1, kind = 8 )

    else if ( ny == 1 ) then

      y = 0.5D+00 * ( ymin + ymax )

    end if

    py = plotymin2 + nint ( alpha * ( y - ymin2 ) )

    write ( unit, '(a)' ) 'newpath'

    px = plotxmin2 + nint ( alpha * ( xmin - xmin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    px = plotxmin2 + nint ( alpha * ( xmax - xmin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do

  return
end
subroutine ps_label ( string )

!*****************************************************************************80
!
!! PS_LABEL prints a label at the current position.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be printed.
!
  implicit none

  character ( len = * ) string
  integer unit

  if ( len_trim ( string ) <= 0 ) then
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) '(' // trim ( string ) // ') show'

  return
end
subroutine ps_label_slant ( string, angle )

!*****************************************************************************80
!
!! PS_LABEL_SLANT prints a slanted label at a given position.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be printed.
!
!    Input, real ( kind = 8 ) ANGLE, the angle of rotation, in degrees.
!
  implicit none

  real ( kind = 8 ) angle
  character ( len = * ) string
  integer unit

  if ( len_trim ( string ) <= 0 ) then
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(f8.4, a)' ) angle, ' rotate'

  write ( unit, '(a)' ) '(' // trim ( string ) // ') show'

  write ( unit, '(f8.4, a)' ) -angle, ' rotate'

  return
end
subroutine ps_line_closed ( npoint, x, y )

!*****************************************************************************80
!
!! PS_LINE_CLOSED adds the graph of a closed line to a PostScript file.
!
!  Discussion:
!
!    A "closed" line is one in which the last point is connected back
!    to the first one.
!
!    The current point is set to the first (and logically last) point
!    in the list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y components
!    of the points.
!
  implicit none

  integer npoint

  real ( kind = 8 ) alpha
  integer i
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(npoint)
  real ( kind = 8 ) ymin
!
!  Refuse to handle fewer than 2 points.
!
  if ( npoint < 2 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE_CLOSED - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Draw lines.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  do i = 2, npoint
    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'
  end do
!
!  Add the final extra segment to the initial point.
!
  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x(1) )
  call ps_setting_real ( 'SET', 'YCUR', y(1) )

  return
end
subroutine ps_line_width ( line_width )

!*****************************************************************************80
!
!! PS_LINE_WIDTH sets the line width.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer LINE_WIDTH, the line width.
!    0 is a valid input, and usually produces the thinnest possible line.
!    1 is a more usual line, 2 is thicker, and so on.
!
  implicit none

  integer line_width
  integer state
  integer unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 2 .and. state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE_WIDTH - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 or 3 is required.'
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(i6,a)' ) line_width, ' setlinewidth'

  call ps_setting_int ( 'SET', 'LINE_WIDTH', line_width )

  return
end
subroutine ps_lineto ( x, y )

!*****************************************************************************80
!
!! PS_LINETO draws a line from the current point to the given point.
!
!  Discussion:
!
!    The current point is updated to the given point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the X and Y components of the new point.
!
  implicit none

  real ( kind = 8 ) alpha
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xcur
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ycur
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XCUR', xcur )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YCUR', ycur )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( xcur - xmin ) )
  py = plotymin2 + nint ( alpha * ( ycur - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x - xmin ) )
  py = plotymin2 + nint ( alpha * ( y - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_mark_circle ( x, y )

!*****************************************************************************80
!
!! PS_MARK_CIRCLE marks a point with a small open circle.
!
!  Discussion:
!
!    The current point is set to the center of the circle.
!
!    The circle is drawn with the current RGB line colors.
!
!    The circle is drawn the current marker size.
!
!    If the point is outside the region, the command is ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point to mark.
!
  implicit none

  real ( kind = 8 ) alpha
  integer marker_size
  integer plotxmin2
  integer plotymin2
  logical point_inside_box_2d
  integer pxcen
  integer pycen
  integer state
  integer unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARK_CIRCLE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )
!
!  If the point is outside the plot box, don't draw it.
!
  if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x, y ) ) then
    return
  end if

  write ( unit, '(a)' ) 'newpath'

  pxcen = plotxmin2 + nint ( alpha * ( x - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y - ymin ) )

  write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
    ' 0 360 arc closepath stroke'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_mark_disk ( x, y )

!*****************************************************************************80
!
!! PS_MARK_DISK marks a point with a small filled disk.
!
!  Discussion:
!
!    The current point is set to the center of the disk.
!
!    The circle is drawn with the current RGB fill colors.
!
!    The circle is drawn the current marker size.
!
!    Points outside the region are not marked.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point to mark.
!
  implicit none

  real ( kind = 8 ) alpha
  integer marker_size
  integer plotxmin2
  integer plotymin2
  logical point_inside_box_2d
  integer pxcen
  integer pycen
  integer state
  integer unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARK_DISK - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )
!
!  If the point is outside the plot box, don't draw it.
!
  if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x, y ) ) then
    return
  end if

  write ( unit, '(a)' ) 'newpath'

  pxcen = plotxmin2 + nint ( alpha * ( x - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y - ymin ) )

  write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
    ' 0 360 arc closepath fill'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_moveto ( x, y )

!*****************************************************************************80
!
!! PS_MOVETO "moves to" a new point, which becomes the current point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the X and Y components of the current point.
!
  implicit none
!
  real ( kind = 8 ) alpha
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MOVETO - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Move to the new point.
!
  px = plotxmin2 + nint ( alpha * ( x - xmin ) )
  py = plotymin2 + nint ( alpha * ( y - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_page_head ( xmin, ymin, xmax, ymax )

!*****************************************************************************80
!
!! PS_PAGE_HEAD writes header information on a new page.
!
!  Discussion:
!
!    I think an earlier version of this code, which wrote
!    "%% Page:" rather than "%%Page:" may have caused problems
!    for some interpreters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, YMIN, XMAX, YMAX, the minimum and maximum X
!    and Y values of the data to be drawn on this page.
!
  implicit none

  real ( kind = 8 ) alpha
  integer num_pages
  integer state
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer margin
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotxmin2
  integer plotymax
  integer plotymin
  integer plotymin2
  integer unit
  real ( kind = 8 ) xcur
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax2
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) ycur
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymax2
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
  real ( kind = 8 ) yvec(4)
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_HEAD - Warning!'
    write ( *, '(a)' ) '  The current open page is forced closed.'
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )

  num_pages = num_pages + 1

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a,i6,i6)' ) '%%Page: ', num_pages, num_pages
  write ( unit, '(a)' ) 'save'
!
!  Reset the state.
!
  state = 3

  call ps_setting_int ( 'SET', 'STATE', state )
!
!  Determine and store parameters.
!
  if ( xmax == xmin ) then
    xmax2 = xmax + 1.0D+00
    xmin2 = xmax - 1.0D+00
  else
    xmax2 = xmax
    xmin2 = xmin
  end if

  if ( ymax == ymin ) then
    ymax2 = ymax + 1.0D+00
    ymin2 = ymax - 1.0D+00
  else
    ymax2 = ymax
    ymin2 = ymin
  end if
!
!  Set the value of "current point".
!
  xcur = xmin
  ycur = ymin
!
!  Set the conversion factors.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin

  alpha = min ( real ( plotxmax - plotxmin, kind = 8 ) / ( xmax2 - xmin2 ), &
                real ( plotymax - plotymin, kind = 8 ) / ( ymax2 - ymin2 ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = nint ( 0.5D+00 * &
    ( real ( plotxmin + plotxmax, kind = 8 ) - alpha * ( xmax2 - xmin2 ) ) )

  plotymin2 = nint ( 0.5D+00 * &
    ( real ( plotymin + plotymax, kind = 8 ) - alpha * ( ymax2 - ymin2 ) ) )
!
!  Store data.
!
  call ps_setting_int ( 'SET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'SET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'SET', 'ALPHA', alpha )
  call ps_setting_real ( 'SET', 'XCUR', xcur )
  call ps_setting_real ( 'SET', 'XMIN', xmin )
  call ps_setting_real ( 'SET', 'XMAX', xmax )
  call ps_setting_real ( 'SET', 'YCUR', ycur )
  call ps_setting_real ( 'SET', 'YMIN', ymin )
  call ps_setting_real ( 'SET', 'YMAX', ymax )
!
!  Make the off-white background.
!
  call ps_comment ( 'Set RGB color to slightly off-white.' )

  line_red = 0.98D+00
  line_green = 0.98D+00
  line_blue = 0.98D+00
  call ps_color_line ( 'PUSH', line_red, line_green, line_blue )
!
!  Paint the background.
!
  call ps_comment ( 'Fill the area with the off-white background color.' )

  xvec(1:4) = (/ xmin, xmax, xmax, xmin /)
  yvec(1:4) = (/ ymin, ymin, ymax, ymax /)

  call ps_polygon_fill ( 4, xvec, yvec )
!
!  Draw a gray border around the page.
!
  line_red = 0.85D+00
  line_green = 0.85D+00
  line_blue = 0.85D+00

  call ps_color_line ( 'PUSH', line_red, line_green, line_blue )

  call ps_comment ( 'Draw a gray border around the page.' )

  xvec(1:4) = (/ xmin, xmax, xmax, xmin /)
  yvec(1:4) = (/ ymin, ymin, ymax, ymax /)

  call ps_line_closed ( 4, xvec, yvec )

  call ps_color_line ( 'POP', line_red, line_green, line_blue )

  return
end
subroutine ps_page_tail ( )

!*****************************************************************************80
!
!! PS_PAGE_TAIL writes tail information at the end of a page.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer state
  integer unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) 'restore showpage'

  call ps_comment ( 'End of page' )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_polygon_fill ( npoint, x, y )

!*****************************************************************************80
!
!! PS_POLYGON_FILL adds a filled polygon to a PostScript file.
!
!  Discussion:
!
!    A closed polygonal path is the sequence of line segments defined
!    by joining consecutive elements of a list of points; the path is
!    closed because the last point is joined to the first.  A filled
!    polygon is the area "inside" a closed polygonal path.  The meaning of
!    the word "inside" can be ambiguous in some cases.
!
!    The polygon fill color should be set before calling this routine.
!
!    The current point is not affected by this call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y components
!    of the points.
!
  implicit none

  integer npoint

  real ( kind = 8 ) alpha
  integer i
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(npoint)
  real ( kind = 8 ) ymin
!
!  Refuse to handle fewer than 2 points.
!
  if ( npoint < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_POLYGON_FILL - Warning!'
    write ( *, '(a)' ) '  Polygon has too few sides.'
    write ( *, '(a,i9)' ) '  NPOINT = ', npoint
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_POLYGON_FILL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Draw lines.
!
  call ps_comment ( 'Draw a polygon' )

  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  do i = 2, npoint
    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'
  end do
!
!  Add the final extra segment to the initial point.
!
  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Fill the polygon.
!
  write ( unit, '(a)' ) 'fill'

  return
end
subroutine ps_sector ( x0, y0, r, angle_min, angle_max )

!*****************************************************************************80
!
!! PS_SECTOR draws a circular sector.
!
!  Discussion:
!
!    As a side effect, the current point is set to the center of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) ANGLE_MIN, ANGLE_MAX, the minimum and
!    maximum angles that define the sector, in degrees.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_min
  integer plotxmin2
  integer plotymin2
  integer pr
  integer pxcen
  integer pycen
  real ( kind = 8 ) r
  integer state
  integer unit
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_SECTOR - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )

  pxcen = plotxmin2 + nint ( alpha * ( x0 - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y0 - ymin ) )
  pr = nint ( alpha * r )

  write ( unit, '(a)' ) 'newpath'
  write ( unit, '(5i6,a)' ) pxcen, pycen, pr, int ( angle_min ), &
    int ( angle_max ), ' arc'
!
!  Draw the circle.
!
  write ( unit, '(a)' ) 'closepath stroke'

  call ps_setting_real ( 'SET', 'XCUR', x0 )
  call ps_setting_real ( 'SET', 'YCUR', y0 )

  return
end
subroutine ps_sector_fill ( x0, y0, r, angle_min, angle_max )

!*****************************************************************************80
!
!! PS_SECTOR_FILL draws a filled circular sector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the disk.
!
!    Input, real ( kind = 8 ) R, the radius of the disk.
!
!    Input, real ( kind = 8 ) ANGLE_MIN, ANGLE_MAX, the minimum and
!    maximum angles that define the sector, in degrees.
!
  implicit none

  integer, parameter :: n = 33

  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_min
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  call sector_points ( x0, y0, r, angle_min, angle_max, n, x, y )

  call ps_comment ( 'Draw a filled polygon.' )

  call ps_polygon_fill ( n, x, y )

  return
end
subroutine ps_setting_int ( action, variable, value )

!*****************************************************************************80
!
!! PS_SETTING_INT sets, gets, or prints integer internal PS_WRITE parameters.
!
!  Discussion:
!
!    Normally, the user does not call this routine.  It is a utility
!    used by the package.
!
!    I'd like a more sophisticated pop and push.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action:
!    'GET' to get the current value of VARIABLE, or
!    'POP' to return the current value and set a new value;
!    'SET' to set a new value of VARIABLE, or
!    'PUSH' to return the current value and set a new value;
!    'PRINT' to print the current value of VARIABLE.
!
!    Input, character ( len = * ) VARIABLE, the variable to get or set:
!    'LINE_WIDTH', the line width.
!      0 is the very thinnest line possible,
!      1 is more usual, 2 is thicker, and so on.
!    'MARKER_SIZE', the size of marker circles and disks, in PostScript points;
!    'NUM_PAGES', the number of pages begun or completed;
!    'PXMIN', the location of the left hand margin of the region
!       in PostScript points;
!    'PYMIN', the location of the lower margin of the region
!       in PostScript points;
!    'STATE', the current internal state,
!      0, file not open,
!      1, file open, no header written, no page open,
!      2, file open, header written, no page open,
!      3, file open, header written, page open.
!      4, file open, header written, trailer written.
!    'UNIT', the FORTRAN output unit associated with the PostScript file.
!
!    Input/output, integer VALUE.
!    If ACTION = 'GET', then VALUE is an output quantity, and is the
!    current internal value of the variable.
!
!    If ACTION = 'SET', then VALUE is an input quantity, and the
!    current internal value of the variable is set to this value.
!
!    If ACTION = 'PRINT', then VALUE is ignored.
!
  implicit none

  character ( len = * ) action
  integer, save :: line_width = 1
  integer, save :: marker_size = 0
  integer, save :: num_pages = 0
  integer, save :: pxmin = 0
  integer, save :: pymin = 0
  integer, save :: state = 0
  integer, save :: unit = 0
  integer value
  character ( len = * ) variable

  if ( variable == 'LINE_WIDTH' ) then

    if ( action == 'GET' ) then
      value = line_width
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Line width, LINE_WIDTH = ', line_width
    else if ( action == 'SET' ) then
      line_width = value
    else if ( action == 'POP' ) then
      call i4_swap ( line_width, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( line_width, value )
    end if

  else if ( variable == 'MARKER_SIZE' ) then

    if ( action == 'GET' ) then
      value = marker_size
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Marker size, MARKER_SIZE = ', marker_size
    else if ( action == 'SET' ) then
      marker_size = value
    else if ( action == 'POP' ) then
      call i4_swap ( marker_size, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( marker_size, value )
    end if

  else if ( variable == 'NUM_PAGES' ) then

    if ( action == 'GET' ) then
      value = num_pages
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Number of pages, NUM_PAGES = ', num_pages
    else if ( action == 'SET' ) then
      num_pages = value
    end if

  else if ( variable == 'PXMIN' ) then

    if ( action == 'GET' ) then
      value = pxmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'PostScript minimum X point, PXMIN = ', pxmin
    else if ( action == 'SET' ) then
      pxmin = value
    else if ( action == 'POP' ) then
      call i4_swap ( pxmin, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( pxmin, value )
    end if

  else if ( variable == 'PYMIN' ) then

    if ( action == 'GET' ) then
      value = pymin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'PostScript minimum Y point, PYMIN = ', pymin
    else if ( action == 'SET' ) then
      pymin = value
    else if ( action == 'POP' ) then
      call i4_swap ( pymin, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( pymin, value )
    end if

  else if ( variable == 'STATE' ) then

    if ( action == 'GET' ) then
      value = state
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Current internal state, STATE = ', state
    else if ( action == 'SET' ) then
      state = value
    else if ( action == 'POP' ) then
      call i4_swap ( state, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( state, value )
    end if

  else if ( variable == 'UNIT' ) then

    if ( action == 'GET' ) then
      value = unit
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Current FORTRAN unit, UNIT = ', unit
    else if ( action == 'SET' ) then
      unit = value
    else if ( action == 'POP' ) then
      call i4_swap ( unit, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( unit, value )
    end if

  end if

  return
end
subroutine ps_setting_real ( action, variable, value )

!*****************************************************************************80
!
!! PS_SETTING_REAL sets, gets, or prints real internal PS_WRITE parameters.
!
!  Discussion:
!
!    I'd like a more sophisticated pop and push.
!
!    This routine has been revised to print an error message and stop
!    if the ACTION or VARIABLE is unexpected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is either:
!    'GET' to get the current value, or
!    'POP' to return the current value and set a new one;
!    'PRINT' to print the current value, or
!    'SET' to set the current value or
!    'PUSH' to set a new value and return the current one.
!
!    Input, character ( len = * ) VARIABLE, the variable to get or set:
!    'ALPHA', the scale factor from XY user space to PostScript points;
!    'FILL_BLUE', the intensity of the blue fill color, between 0.0 and 1.0.
!    'FILL_GREEN', the intensity of the green fill color, between 0.0 and 1.0.
!    'FILL_RED', the intensity of the red fill color, between 0.0 and 1.0.
!    'FONT_SIZE', the font size, in inches.
!    'LINE_BLUE', the blue component of the line color, between 0.0 and 1.0.
!    'LINE_GREEN', the green component of the line color, between 0.0 and 1.0.
!    'LINE_RED', the red component of the line color, between 0.0 and 1.0.
!    'XCUR', the current X location.
!    'XMAX', maximum X value of the data.
!    'XMIN', minimum X value of the data.
!    'YCUR', the current Y location.
!    'YMAX', maximum Y value of the data.
!    'YMIN', minimum Y value of the data.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If ACTION = 'GET', then VALUE is an output quantity, and is the
!    current internal value of the variable.
!
!    If ACTION = 'SET', then VALUE is an input quantity, and the
!    current internal value of the variable is set to this value.
!
!    If ACTION = 'PRINT', then VALUE is ignored.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 0.0D+00
  real ( kind = 8 ), save :: fill_blue = 0.7D+00
  real ( kind = 8 ), save :: fill_green = 0.7D+00
  real ( kind = 8 ), save :: fill_red = 0.7D+00
  real ( kind = 8 ), save :: font_size = 0.1D+00
  real ( kind = 8 ), save :: line_blue = 0.0D+00
  real ( kind = 8 ), save :: line_green = 0.0D+00
  real ( kind = 8 ), save :: line_red = 0.0D+00
  real ( kind = 8 ) value
  character ( len = * ) variable
  real ( kind = 8 ), save :: xcur = 0.0D+00
  real ( kind = 8 ), save :: xmax = 1.0D+00
  real ( kind = 8 ), save :: xmin = 0.0D+00
  real ( kind = 8 ), save :: ycur = 0.0D+00
  real ( kind = 8 ), save :: ymax = 0.0D+00
  real ( kind = 8 ), save :: ymin = 0.0D+00

  if ( variable == 'ALPHA' ) then

    if ( action == 'GET' ) then
      value = alpha
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Scale factor from user to PS, ALPHA = ', alpha
    else if ( action == 'SET' ) then
      alpha = value
    else if ( action == 'POP' ) then
      call r8_swap ( alpha, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( alpha, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_BLUE' ) then

    if ( action == 'GET' ) then
      value = fill_blue
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Blue fill RGB value, FILL_BLUE = ', fill_blue
    else if ( action == 'SET' ) then
      fill_blue = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_blue, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_blue, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_GREEN' ) then

    if ( action == 'GET' ) then
      value = fill_green
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Green fill RGB value, FILL_GREEN = ', fill_green
    else if ( action == 'SET' ) then
      fill_green = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_green, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_green, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_RED' ) then

    if ( action == 'GET' ) then
      value = fill_red
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'RED fill RGB value, FILL_RED = ', fill_red
    else if ( action == 'SET' ) then
      fill_red = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_red, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_red, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FONT_SIZE' ) then

    if ( action == 'GET' ) then
      value = font_size
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Font size, FONT_SIZE = ', font_size
    else if ( action == 'SET' ) then
      font_size = value
    else if ( action == 'POP' ) then
      call r8_swap ( font_size, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( font_size, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_BLUE' ) then

    if ( action == 'GET' ) then
      value = line_blue
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Blue line RGB value, LINE_BLUE = ', line_blue
    else if ( action == 'SET' ) then
      line_blue = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_blue, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_blue, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_GREEN' ) then

    if ( action == 'GET' ) then
      value = line_green
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Green line RGB value, LINE_GREEN = ', line_green
    else if ( action == 'SET' ) then
      line_green = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_green, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_green, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_RED' ) then

    if ( action == 'GET' ) then
      value = line_red
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Red line RGB value, LINE_RED = ', line_red
    else if ( action == 'SET' ) then
      line_red = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_red, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_red, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XCUR' ) then

    if ( action == 'GET' ) then
      value = xcur
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Current X location, XCUR = ', xcur
    else if ( action == 'SET' ) then
      xcur = value
    else if ( action == 'POP' ) then
      call r8_swap ( xcur, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xcur, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XMAX' ) then

    if ( action == 'GET' ) then
      value = xmax
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Maximum X value, XMAX = ', xmax
    else if ( action == 'SET' ) then
      xmax = value
    else if ( action == 'POP' ) then
      call r8_swap ( xmax, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xmax, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XMIN' ) then

    if ( action == 'GET' ) then
      value = xmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Minimum X value, XMIN = ', xmin
    else if ( action == 'SET' ) then
      xmin = value
    else if ( action == 'POP' ) then
      call r8_swap ( xmin, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xmin, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YCUR' ) then

    if ( action == 'GET' ) then
      value = ycur
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Current Y location, YCUR = ', ycur
    else if ( action == 'SET' ) then
      ycur = value
    else if ( action == 'POP' ) then
      call r8_swap ( ycur, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ycur, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YMAX' ) then

    if ( action == 'GET' ) then
      value = ymax
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Maximum Y value, YMAX = ', ymax
    else if ( action == 'SET' ) then
      ymax = value
    else if ( action == 'POP' ) then
      call r8_swap ( ymax, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ymax, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YMIN' ) then

    if ( action == 'GET' ) then
      value = ymin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Minimum Y value, YMIN = ', ymin
    else if ( action == 'SET' ) then
      ymin = value
    else if ( action == 'POP' ) then
      call r8_swap ( ymin, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ymin, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
    write ( *, '(a)' ) '  Unexpected variable!'
    stop

  end if

  return
end
subroutine ps_square ( x0, y0, r )

!*****************************************************************************80
!
!! PS_SQUARE draws a square.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the square.
!
!    Input, real ( kind = 8 ) R, the radius of the square.
!
  implicit none

  integer, parameter :: n = 4

  real ( kind = 8 ) r
  real ( kind = 8 ) xvec(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) yvec(n)
  real ( kind = 8 ) y0

  xvec(1:n) = (/ x0-r, x0+r, x0+r, x0-r /)
  yvec(1:n) = (/ y0-r, y0-r, y0+r, y0+r /)

  call ps_line_closed ( n, xvec, yvec )

  return
end
subroutine ps_square_fill ( x0, y0, r )

!*****************************************************************************80
!
!! PS_SQUARE_FILL draws a filled square.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the box.
!
!    Input, real ( kind = 8 ) R, the radius of the box.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) x(4)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(4)
  real ( kind = 8 ) y0

  x(1:4) = (/ x0-r, x0+r, x0+r, x0-r /)
  y(1:4) = (/ y0-r, y0-r, y0+r, y0+r /)

  call ps_polygon_fill ( 4, x, y )

  return
end
subroutine ps_star ( x0, y0, r )

!*****************************************************************************80
!
!! PS_STAR draws an open star of given radius.
!
!  Discussion:
!
!    The radius refers to the circumscribing circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the star.
!
!    Input, real ( kind = 8 ) R, the radius of the star.
!
  implicit none

  integer, parameter :: n = 5

  real ( kind = 8 ) angle
  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) xvec(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) yvec(n)
  real ( kind = 8 ) y0

  do i = 1, n
    angle = pi * real ( mod ( ( - 7 + i * 8 ), 20 ), kind = 8 ) / 10.0D+00
    xvec(i) = x0 + r * cos ( angle )
    yvec(i) = y0 + r * sin ( angle )
  end do

  call ps_line_closed ( n, xvec, yvec )

  return
end
subroutine ps_star_circle ( x0, y0, r )

!*****************************************************************************80
!
!! PS_STAR_CIRCLE draws open circles at the 10 points on a star of given radius.
!
!  Discussion:
!
!    The radius refers to the circumscribing circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the star.
!
!    Input, real ( kind = 8 ) R, the radius of the star.
!
  implicit none

  integer, parameter :: n = 10

  integer i
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  call ps_star_points ( x0, y0, r, x, y )

  do i = 1, n
    call ps_mark_circle ( x(i), y(i) )
  end do

  return
end
subroutine ps_star_disk ( x0, y0, r )

!*****************************************************************************80
!
!! PS_STAR_DISK draws filled disks at the 10 points on a star of given radius.
!
!  Discussion:
!
!    The radius refers to the circumscribing circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the star.
!
!    Input, real ( kind = 8 ) R, the radius of the star.
!
  implicit none

  integer, parameter :: n = 10

  integer i
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  call ps_star_points ( x0, y0, r, x, y )

  do i = 1, n
    call ps_mark_disk ( x(i), y(i) )
  end do

  return
end
subroutine ps_star_points ( x0, y0, r, x, y )

!*****************************************************************************80
!
!! PS_STAR_POINTS returns 10 points on a star of given radius.
!
!  Discussion:
!
!    The radius refers to the circumscribing circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the star.
!
!    Input, real ( kind = 8 ) R, the radius of the star.
!
!    Output, real ( kind = 8 ) X(10), Y(10), the coordinates of points
!    on the star.
!
  implicit none

  integer, parameter :: n = 10

  real ( kind = 8 ) angle
  integer i
  integer j
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  phi = ( 1.0D+00 + sqrt ( 5.0D+00 ) ) / 2.0D+00

  j = 0
  do i = 1, 5
    j = j + 1
    angle = pi * real ( mod ( ( - 7 + i * 8 ), 20 ), kind = 8 ) / 10.0D+00
    x(j) = x0 + r * cos ( angle )
    y(j) = y0 + r * sin ( angle )
    j = j + 1
    angle = pi * real ( mod ( ( - 5 + i * 8 ), 20 ), kind = 8 ) / 10.0D+00
    x(j) = x0 + r * cos ( angle ) / phi**2
    y(j) = y0 + r * sin ( angle ) / phi**2
  end do

  return
end
subroutine ps_triangle ( x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! PS_TRIANGLE draws an open triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates
!    of the triangle.
!
  implicit none

  integer, parameter :: n = 3

  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xvec(n)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yvec(n)

  xvec(1:n) = (/ x1, x2, x3 /)
  yvec(1:n) = (/ y1, y2, y3 /)

  call ps_line_closed ( n, xvec, yvec )

  return
end
subroutine ps_triangle_fill ( x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! PS_TRIANGLE_FILL draws a filled triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates
!    of the triangle.
!
  implicit none

  integer, parameter :: n = 3

  real ( kind = 8 ) xvec(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) yvec(n)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  xvec(1:n) = (/ x1, x2, x3 /)
  yvec(1:n) = (/ y1, y2, y3 /)

  call ps_polygon_fill ( n, xvec, yvec )

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
function s_index_last ( s, sub )

!*****************************************************************************80
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!  Discussion:
!
!    It returns the location in the string at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either S or SUB is all blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none

  integer i
  integer j
  integer llen1
  integer llen2
  character ( len = * ) s
  integer s_index_last
  character ( len = * ) sub

  s_index_last = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen1 < llen2 ) then
    return
  end if

  do j = 1, llen1+1-llen2

    i = llen1 + 2 - llen2 - j

    if ( s(i:i+llen2-1) == sub ) then
      s_index_last = i
      return
    end if

  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

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
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
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
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
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
      else if ( 1 < ihave ) then
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
      else if ( 6 <= ihave .and. ihave <= 8 ) then
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
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
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
subroutine sector_points ( x0, y0, r, angle_min, angle_max, n, x, y )

!*****************************************************************************80
!
!! SECTOR_POINTS returns N equally spaced points on a circle in 2D.
!
!  Discussion:
!
!    The first point is always ( X0, Y0 ), and the N-1 subsequent points
!    proceed counterclockwise around the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) ANGLE_MIN, ANGLE_MAX, the minimum and
!    maximum angles that define the sector, in degrees.
!
!    Input, integer N, the number of points desired.  N must be at least 3.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points
!    on the circle.
!
  implicit none

  integer n

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_min
  real ( kind = 8 ) degrees_to_radians
  integer i
  real ( kind = 8 ) r
  real ( kind = 8 ) x0
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y(n)

  x(1) = x0
  y(1) = y0

  do i = 2, n
    angle = ( real ( n - i,     kind = 8 ) * angle_min   &
            + real (     i - 2, kind = 8 ) * angle_max ) &
            / real ( n     - 2, kind = 8  )
    angle = degrees_to_radians ( angle )
    x(i) = x0 + r * cos ( angle )
    y(i) = y0 + r * sin ( angle )
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
subroutine word_extract ( s, w )

!*****************************************************************************80
!
!! WORD_EXTRACT extracts the next word from a string.
!
!  Discussion:
!
!    A "word" is a string of characters terminated by a blank or
!    the end of the string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.  On output, the first
!    word has been removed, and the remaining string has been shifted left.
!
!    Output, character ( len = * ) W, the leading word of the string.
!
  implicit none

  integer iget1
  integer iget2
  integer lchar
  character ( len = * ) s
  character ( len = * ) w

  w = ' '

  lchar = len_trim ( s )
!
!  Find the first nonblank.
!
  iget1 = 0

  do

    iget1 = iget1 + 1

    if ( lchar < iget1 ) then
      return
    end if

    if ( s(iget1:iget1) /= ' ' ) then
      exit
    end if

  end do
!
!  Look for the last contiguous nonblank.
!
  iget2 = iget1

  do

    if ( lchar <= iget2 ) then
      exit
    end if

    if ( s(iget2+1:iget2+1) == ' ' ) then
      exit
    end if

    iget2 = iget2 + 1

  end do
!
!  Copy the word.
!
  w = s(iget1:iget2)
!
!  Shift the string.
!
  s(1:iget2) = ' '
  s = adjustl ( s(iget2+1:) )

  return
end
