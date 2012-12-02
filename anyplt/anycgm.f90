subroutine anyplt ( icom )

!*****************************************************************************80
!
!! ANYPLT is an interface routine to a variety of graphics packages.
!
!  Discussion:
!
!    ANYPLT is a subroutine which provides a simple, standard interface
!    between FORTRAN programs and various output devices.  To run a
!    program which calls ANYPLT on a different machine, the program
!    is not modified in any way, but a different version of the ANYPLT
!    program is provided.  Currently, the following versions are available:
!
!    ANYATT - AT&T PC6300 graphics (640 by 400).  Requires ATTPLT.ASM.
!    ANYBUG - Simple debugging output to a file.
!    ANYCAL - CALCOMP file output.  Available on many mainframes.
!    ANYIBM - IBM PC hi resolution (640 by 200).  Requires IBMPLT.ASM.
!    ANYMAC - Macintosh graphics.  Requires auxilliary routine TOOLBX.SUB.
!    ANYNCR - NCAR graphics package.
!    ANYNUL - Does nothing.
!    ANYP10 - PLOT10 interactive graphics. (1024 by 768)
!    ANYTTY - Simple 'typewriter' graphics (80 by 24)
!
!    The symbolic output of characters, numbers and other printable
!    characters was made possible by adaptation of a routine written
!    by Bill Furey of the University of Pittsburgh, Department of
!    Crystallography.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ICOM, specifies the graphics request being made.
!    0, enable graphics.
!    1, disable graphics.
!    2, begin plot.
!    3, define plot size.
!    4, move to a point.
!    5, draw to a point.
!    6, clear screen.
!    7, write string at position.
!    8, use virtual cursor.
!    9, end plot.
!    10, ring bell.
!    11, mark data.
!    12, return screen data.
!    13, return version.
!    14, draw an arrow.
!
  implicit none

  real angle
  character ( len = 80 ) carray
  real csize
  real cwide
  real degrees_to_radians
  character ( len = 10 ), save :: dev
  logical filled
  character flush
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icmax
  integer ( kind = 4 ) icmin
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  integer ( kind = 4 ) itable
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ) lent
  integer ( kind = 4 ) marray
  integer ( kind = 4 ) ndraw
  integer ( kind = 4 ), save :: nplot = 0
  integer ( kind = 4 ) nval
  real pwide
  logical s_eqi
  real x2
  real xdraw(6)
  real xplt1
  real xplt2
  real, save :: xpmax
  real, save :: xpmin
  real xval(2)
  real y2
  real ydraw(6)
  real yplt1
  real yplt2
  real, save :: ypmax
  real, save :: ypmin
  real yval(2)

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2
  common /anychr/ carray
!
!  ICOM = 0  Enable graphics
!
  if ( icom == 0 ) then

    if ( s_eqi ( carray, 'cgm' ) ) then
      dev = 'cgmb'
    else if ( s_eqi ( carray, 'cgmb' ) ) then
      dev = 'cgmb'
    else if ( s_eqi ( carray, 'ps' ) ) then
      dev = 'ps'
    else if ( s_eqi ( carray, 'xws' ) ) then
      dev = 'xws'
    else
      return
    end if

    write ( *, * ) 'Using graphics device ' // trim ( dev )

    call device ( dev )

    if ( dev == 'cgmb' ) then
      call outfil ( 'anyplt.cgm' )
    else if ( dev == 'ps' ) then
      call outfil ( 'anyplt.ps' )
    end if

    icmax = 200
    icmin = 2
    itable = 1
    call color_table_set ( icmin, icmax, itable )

    nplot = 0
!
!  ICOM = 1  Disable graphics
!
  else if ( icom == 1 ) then
    call grfcls
!
!  ICOM = 2  Begin plot
!
  else if ( icom == 2 ) then

    if ( nplot == 0 ) then

      call grfini

      xval(1) = 0.0E+00
      xval(2) = 1.0E+00
      yval(1) = 0.0E+00
      yval(2) = 1.0E+00
      nval = 2
      call setscl ( xval, yval, nval )
      nplot = 1

    else

      call newfrm
      nplot = nplot + 1

    end if
!
!  ICOM = 3  Define plot size
!
  else if ( icom == 3 ) then

    xpmin = xplt1
    xpmax = xplt1 + xplt2
    ypmin = yplt1
    ypmax = yplt1 + yplt2

    xval(1) = xpmin
    xval(2) = xpmax
    yval(1) = ypmin
    yval(2) = ypmax
    nval = 2

    call setscl ( xval, yval, nval )
!
!  ICOM = 4  Move to point
!
  else if ( icom == 4 ) then

    call movcgm ( xplt1, yplt1 )
!
!  ICOM = 5  Draw to point
!
  else if ( icom == 5 ) then

    call drwcgm ( xplt1, yplt1 )
!
!  ICOM = 6  Clear screen
!
  else if ( icom == 6 ) then
!
!    call newfrm
!
!  ICOM = 7,  Write string at position
!
  else if ( icom == 7 ) then

    angle = 0.0E+00
    cwide = 0.025E+00
    pwide = 1.0E+00
    lent = len_trim ( carray )
    flush = 'c'

    call s_plot ( angle, cwide, pwide, trim ( carray ), xplt1, yplt1, flush )
!
!  ICOM = 8  Use virtual cursor
!
  else if ( icom == 8 ) then
!
!  ICOM = 9  End plot
!
  else if ( icom == 9 ) then

    if ( s_eqi ( dev, 'XWS' ) ) then

      call movcgm ( xpmin, ypmin )

      do i = 1, 100
        call drwcgm ( xpmax, ypmin )
        call drwcgm ( xpmax, ypmax )
        call drwcgm ( xpmin, ypmax )
        call drwcgm ( xpmin, ypmin )
      end do

    end if
!
!  ICOM = 10  Ring bell
!
  else if ( icom == 10 ) then
!
!  ICOM = 11  Mark data
!
  else if ( icom == 11 ) then

    filled = .false.
    csize = 0.005E+00
    call circle ( xplt1, yplt1, csize, filled )
!
!  ICOM = 12  Return screen data
!
  else if ( icom == 12 ) then
!
!  ICOM = 13  Return version
!
  else if ( icom == 13 ) then
    carray = 'ANYPLT - Version 1.03  21 November 2000  CGMPLT'
!
!  ICOM = 14, draw an arrow.
!
  else if ( icom == 14 ) then
    x2 = xplt1 + yplt2 * cos ( degrees_to_radians ( xplt2 ) )
    y2 = yplt1 + yplt2 * sin ( degrees_to_radians ( xplt2 ) )
    call arrow ( xplt1, yplt1, x2, y2, xdraw, ydraw )
    call movcgm ( xdraw(1), ydraw(1) )
    do i = 2, 6
      call drwcgm ( xdraw(i), ydraw(i) )
    end do
!
!  Unknown value of ICOM.
!
  else
    write ( *, * ) 'ANYPLT - Fatal error!'
    write ( *, * ) '  Unknown value of ICOM = ',icom
    stop
  end if

  return
end
subroutine arrow ( xstart, ystart, xtip, ytip, xdraw, ydraw )

!*****************************************************************************80
!
!! ARROW returns points that specify an arrow from one point to another.
!
!  Discussion:
!
!    The arrow will stretch between two user specified points.
!
!    The "head" of the arrow may be fatter or thinner than expected
!    if the X and Y scales of the graph are not in the same
!    proportions.
!
!
!                       left(3)
!                        |\
!                        | \
!                        |  \
!    start(1)*-----base(2,6) * tip(4)
!                        |  /
!                        | /
!                        |/
!                       rite(5)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XSTART, YSTART, the starting point for the arrow.
!
!    Input, real XTIP, YTIP, the end point for the arrow.
!
!    Output, real XDRAW(6), YDRAW(6), the X and Y coordinates
!    of the points to connect to draw the arrow.
!
  implicit none

  real alpha
  real del
  real dist
  real, parameter :: pi = 3.141592653589793E+00
  real theta
  real xbase
  real xdraw(6)
  real xleft
  real xrite
  real xstart
  real xtip
  real ybase
  real ydraw(6)
  real yleft
  real yrite
  real ystart
  real ytip

  theta = 0.5E+00 * pi - atan2 ( 2.0E+00, 1.0E+00 )
  dist = sqrt ( ( xtip - xstart )**2 + ( ytip - ystart )**2 )
  alpha = atan2 ( ytip - ystart, xtip - xstart )
  del = sqrt ( 5.0E+00 ) / 3.0E+00

  xbase = ( xstart + 2.0E+00 * xtip ) / 3.0E+00
  ybase = ( ystart + 2.0E+00 * ytip ) / 3.0E+00

  xleft = xstart + del * dist * cos ( alpha - theta )
  yleft = ystart + del * dist * sin ( alpha - theta )

  xrite = xstart + del * dist * cos ( alpha + theta )
  yrite = ystart + del * dist * sin ( alpha + theta )

  xdraw(1:6) = (/ xstart, xbase, xleft, xtip, xrite, xbase /)
  ydraw(1:6) = (/ ystart, ybase, yleft, ytip, yrite, ybase /)

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
subroutine color_table_set ( icmin, icmax, itable )

!*****************************************************************************80
!
!! COLOR_TABLE_SET sets up the color table.
!
!  Discussion:
!
!    This routine replaces the unreliable DRAWCGM routine SETCTB.
!
!    The routine sets the colors between ICMIN and ICMAX, which
!    should typically be 2 and 255.
!
!    It will also set the values of color 0 to white, and
!    color 1 to black.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ICMIN, ICMAX, the minimum and maximum color indices.
!
!    Input, integer ( kind = 4 ) ITABLE, the desired table.
!    1: low black to high white
!    2: low white to high black.
!    3: low blue to high yellow
!    4: low red, high blue, with bands between.
!    5: low red, yellow, green, blue, high white.
!    5: low white, blue, green, yellow, high red.
!    7: low blue to high red
!
  implicit none

  real bhi
  real blo
  real bval
  real ghi
  real glo
  real gval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icmax
  integer ( kind = 4 ) icmin
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) itable
  integer ( kind = 4 ) ival
  real pi
  real rhi
  real rlo
  real rval
  real temp
  real theta

  icmin = max ( icmin, 2 )
  icmax = min ( icmax, 255 )
!
!  1: Low black to high white
!  2: Low white to high black
!  3: Low blue to high yellow.
!  4: Low red, high blue, with bands between.
!  5: Low red, yellow, green, blue, high white.
!  6: Low white, blue, green, yellow, high red.
!  7: Low blue to high red
!
    do i = icmin, icmax

      if ( icmin == icmax ) then
        temp = 0.5E+00
      else
        temp = real ( i - icmin ) / real ( icmax - icmin )
      end if

      if ( itable == 1 ) then
        bval = temp
        gval = temp
        rval = temp
      else if ( itable == 2 ) then
        bval = 1.0E+00 - temp
        gval = 1.0E+00 - temp
        rval = 1.0E+00 - temp
      else if ( itable == 3 ) then
        rval = temp
        gval = temp
        bval = 1.0E+00 - temp
      else if ( itable == 4 ) then
        theta = 0.5E+00 * pi() * temp
        rval = cos ( theta )**2
        bval = sin ( theta )**2
        gval = 0.8E+00 * sin ( 10.0E+00 * theta )**6
      else if ( itable == 5 ) then
        theta = 4.0E+00 * temp
        rval = exp(-(theta-1.0E+00)**2) + exp(-(theta-4.0E+00)**2)
        gval = exp(-(theta-2.0E+00)**2) + exp(-(theta-4.0E+00)**2)
        bval = exp(-(theta-3.0E+00)**2) + exp(-(theta-4.0E+00)**2)
        rval = max ( rval, 0.0E+00 )
        rval = min ( rval, 1.0E+00 )
        gval = max ( gval, 0.0E+00 )
        gval = min ( gval, 1.0E+00 )
        bval = max ( bval, 0.0E+00 )
        bval = min ( bval, 1.0E+00 )
      else if ( itable == 6 ) then
        theta = 4.0E+00 * temp
        rval = exp(-(theta-1.0E+00)**2) + exp(-(theta-4.0E+00)**2)
        gval = exp(-(theta-2.0E+00)**2) + exp(-(theta-4.0E+00)**2)
        bval = exp(-(theta-3.0E+00)**2) + exp(-(theta-4.0E+00)**2)
        rval = max ( rval, 0.0E+00 )
        rval = min ( rval, 1.0E+00 )
        gval = max ( gval, 0.0E+00 )
        gval = min ( gval, 1.0E+00 )
        bval = max ( bval, 0.0E+00 )
        bval = min ( bval, 1.0E+00 )
      else if ( itable == 7 ) then
        rval = temp
        gval = 0.0E+00
        bval = 1.0E+00 - temp
      end if

      ival = i
      call setclr ( ival, rval, gval, bval )

    end do
!
!  Background color 0 is to be white.
!
  ival = 0
  rval = 1.0E+00
  gval = 1.0E+00
  bval = 1.0E+00
  call setclr ( ival, rval, gval, bval )
!
!  Foreground color 1 is to be black.
!
  ival = 1
  rval = 0.0E+00
  gval = 0.0E+00
  bval = 0.0E+00
  call setclr ( ival, rval, gval, bval )

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
!    Input, real ANGLE, an angle in degrees.
!
!    Output, real DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  implicit none

  real angle
  real degrees_to_radians
  real, parameter :: pi = 3.141592653589793E+00

  degrees_to_radians = ( angle / 180.0E+00 ) * pi

  return
end
function s_eqi ( s1, s2 )

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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_plot ( angle, cwide, pwide, s, x, y, flush )

!*****************************************************************************80
!
!! S_PLOT plots a character string onto a graphics image.
!
!  Discussion:
!
!    The string can be at any angle and at any size.
!
!    The plot is assumed to be of size PWIDE by PHITE, although PHITE
!    itself is not input.
!
!    This routine must be modified to work with a particular graphics package.
!    The current code calls two routines:
!      MOVCGM ( X, Y ) moves to a point (X,Y) in the plot;
!      DRWCGM ( X, Y ) draws a line from the current point to (X,Y).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle in degrees at which the
!    string is to be drawn.  0 is typical.  90 degrees would
!    cause the string to be written from top to bottom.
!
!    Input, real CWIDE, the width of the characters.  This
!    is measured in the same units as the plot width PWIDE.
!    For PWIDE = 1, a plot size of 0.025 would be reasonable,
!    since 40 characters would fit, but 2.0 would be nonsense.
!
!    Input, real PWIDE, the width of the plot, in the same
!    units as CWIDE.
!
!    Input, character ( len = * ) S, contains the text to be plotted.
!    Only characters with ASCII codes between 32 and 126 will actually
!    be plotted.  Any other characters are "unprintable", and will be
!    plotted as blanks.
!
!    Input, real X, Y, the coordinates of a point which
!    determines where the string is drawn.  The string will
!    be drawn starting at, centered or, or ending at (X,Y),
!    depending on the value of FLUSH.
!
!    Input, character ( len = * ) FLUSH, a string which specifies how to
!    place the string.  Only the first character of FLUSH is examined, and
!    the case of the character is not important.
!
!    'L' - the string will be drawn flush left.
!    'C' - the string will be centered.
!    'R' - the string will be drawn flush right.
!
  implicit none

  real angle
  real ca
  character c
  real cwide
  real degrees_to_radians
  character ( len = * ) flush
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iascii
  integer ( kind = 4 ) icr
  integer ( kind = 4 ), save, dimension ( 1617 ) :: ifont
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipen
  integer ( kind = 4 ), save, dimension ( 95 ) :: ipoint
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) nvec
  real pwide
  logical rotate
  character ( len = * ) s
  real sa
  real scl2
  real x
  real xb
  real xc
  real xcopy
  real xnew
  real xold
  real xrot
  real xt
  real y
  real yb
  real yc
  real ycopy
  real ynew
  real yold
  real yrot
  real yt
!
!  IPOINT is a pointer array into IFONT.
!
!  IPOINT(I) records where the "strokes" for character I begin
!  in the IFONT array.
!
  data ( ipoint(i), i = 1, 95 ) / &
       1,   3,  26,  45,  66, 102, 130, 156, 166, 186, 206, 222, 233, &
     249, 255, 267, 273, 293, 306, 328, 353, 363, 383, 411, 423, 457, &
     483, 506, 533, 541, 552, 560, 587, 625, 638, 665, 683, 699, 714, &
     727, 754, 770, 786, 805, 818, 826, 838, 848, 868, 884, 909, 930, &
     956, 967, 981, 989,1001,1012,1025,1035,1045,1051,1061,1069,1075, &
    1081,1108,1131,1149,1172,1194,1214,1243,1260,1284,1307,1323,1336, &
    1364,1381,1401,1424,1447,1464,1486,1499,1516,1524,1536,1547,1560, &
    1570,1586,1592,1608 /
!
!  IFONT contains the strokes defining the various symbols.
!
  data ( ifont(i), i = 1, 396 ) / &
     1, 0, 2,10,11, 9,22,10,23,11,22,10,11, 0, 9, 7, 9, 9,11, 9,11, 7, &
     9, 7, 0, 2, 8,17, 7,23, 9,23, 8,17, 0,14,17,13,23,15,23,14,17, 0, &
     4, 9,23, 7, 7, 0,13,23,11, 7, 0, 5,17,15,17, 0, 5,13,15,13, 0, 3, &
    15,19,13,21, 9,21, 7,19, 7,17, 9,15,13,15,15,13,15,11,13, 9, 9, 9, &
     7,11, 0, 9,23, 9, 7, 0,13,23,13, 7, 0, 3, 5,23, 9,23, 9,19, 5,19, &
     5,23, 0,17,23, 5, 7, 0,13, 7,13,11,17,11,17, 7,13, 7, 0, 1,17, 7, &
     7,17, 7,19, 9,21,13,21,15,19,15,17, 5,13, 5,11, 9, 7,13, 7,17,15, &
     0, 1,10,17, 9,23,11,23,10,17, 0, 1,12,23,11,21,10,19, 9,17, 9,15, &
     9,13,10,11,11, 9,12, 7, 0, 1,12,23,13,21,14,19,15,17,15,15,15,13, &
    14,11,13, 9,12, 7, 0, 3, 7,15,15,15, 0,13,19, 9,11, 0, 9,19,13,11, &
     0, 2, 7,15,15,15, 0,11,19,11,11, 0, 1,11, 7, 9, 7, 9, 9,11, 9,11, &
     7,11, 6,10, 4, 0, 1, 7,15,15,15, 0, 1, 9, 7, 9, 9,11, 9,11, 7, 9, &
     7, 0, 1,15,23, 7, 7, 0, 1, 9,23,13,23,15,19,15,11,13, 7, 9, 7, 7, &
    11, 7,19, 9,23, 0, 2, 7,21, 9,23, 9, 7, 0, 7, 7,11, 7, 0, 1, 5,21, &
     9,23,15,23,17,21,17,19,15,17, 7,13, 5,10, 5, 7,17, 7, 0, 2, 5,23, &
    17,23,15,17,13,15, 9,15, 0,13,15,17,13,17,10,14, 7, 8, 7, 5,10, 0, &
     1,13, 7,13,23, 5,13,17,13, 0, 1,17,23, 5,23, 5,17,13,17,17,15,17, &
    11,13, 7, 9, 7, 5,11, 0, 1,17,19,13,23, 9,23, 5,19, 5,13, 9,15,13 /

  data ( ifont(i), i =  397, 792 ) / &
    15,17,13,17,11,13, 7, 9, 7, 5,11, 5,13, 0, 1, 5,19, 5,23,17,23,11, &
    15,11, 7, 0, 1, 8,15, 6,17, 6,21, 8,23,14,23,16,21,16,17,14,15, 8, &
    15, 5,13, 5, 9, 8, 7,14, 7,17, 9,17,13,14,15, 0, 1,17,17,15,15, 7, &
    15, 5,17, 5,21, 7,23,15,23,17,21,17,11,15, 7, 7, 7, 5,11, 0, 2, 9, &
    13, 9,15,11,15,11,13, 9,13, 0, 9, 7, 9, 9,11, 9,11, 7, 9, 7, 0, 2, &
     9,13, 9,15,11,15,11,13, 9,13, 0,11, 7, 9, 7, 9, 9,11, 9,11, 7,11, &
     6,10, 4, 0, 1,17,21, 5,15,17, 9, 0, 2, 7,15,15,15, 0, 7, 9,15, 9, &
     0, 1, 5,21,17,15, 5, 9, 0, 2, 7,21, 9,23,13,23,15,21,15,19,11,15, &
    11,11, 0,10, 7,10, 9,12, 9,12, 7,10, 7, 0, 1,13, 7, 9, 7, 5,11, 5, &
    19, 9,23,13,23,17,19,17,11,15, 9,13,11,12,10,10,10, 9,11, 9,15,10, &
    16,12,16,13,15,13,11, 0, 2, 5, 7,11,23,17, 7, 0, 8,15,14,15, 0, 2, &
     5, 7, 5,23,15,23,17,21,17,17,15,15, 5,15, 0,15,15,17,13,17, 9,15, &
     7, 5, 7, 0, 1,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11, 0, &
     1, 5, 7, 5,23,13,23,17,19,17,11,13, 7, 5, 7, 0, 2,17,23, 5,23, 5, &
     7,17, 7, 0, 5,15,12,15, 0, 2, 5, 7, 5,23,17,23, 0, 5,15,12,15, 0, &
     2,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11,17,15,13,15, 0, &
    17,11,17, 7, 0, 3, 5, 7, 5,23, 0, 5,15,17,15, 0,17,23,17, 7, 0, 3, &
     9,23,13,23, 0,11,23,11, 7, 0, 9, 7,13, 7, 0, 2,15,23,15,11,12, 7 /

  data ( ifont(i), i =  793, 1188 ) / &
     8, 7, 5,11, 5,13, 0,13,23,17,23, 0, 2, 5, 7, 5,23, 0,17,23, 5,15, &
    17, 7, 0, 1, 5,23, 5, 7,17, 7, 0, 1, 5, 7, 5,23,11,11,17,23,17, 7, &
     0, 1, 5, 7, 5,23,17, 7,17,23, 0, 1,17,19,13,23, 9,23, 5,19, 5,11, &
     9, 7,13, 7,17,11,17,19, 0, 1, 5, 7, 5,23,13,23,17,21,17,17,13,15, &
     5,15, 0, 2,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11,17,19, &
     0,13,11,17, 7, 0, 2, 5, 7, 5,23,13,23,17,21,17,17,13,15, 5,15, 0, &
    13,15,17, 7, 0, 1,17,19,13,23, 9,23, 5,20, 5,18, 9,15,13,15,17,12, &
    17,10,13, 7, 9, 7, 5,10, 0, 2, 5,23,17,23, 0,11,23,11, 7, 0, 1, 5, &
    23, 5,10, 8, 7,14, 7,17,10,17,23, 0, 1, 5,23,11, 7,17,23, 0, 1, 5, &
    23, 8, 7,11,17,14, 7,17,23, 0, 2, 5,23,17, 7, 0,17,23, 5, 7, 0, 2, &
     5,23,11,13,17,23, 0,11,13,11, 7, 0, 1, 5,23,17,23, 5, 7,17, 7, 0, &
     1,11,23, 7,23, 7, 7,11, 7, 0, 1, 7,23,15, 7, 0, 1, 7,23,11,23,11, &
     7, 7, 7, 0, 1, 7,21,11,23,15,21, 0, 1, 5, 3,17, 3, 0, 1, 9,23,13, &
    19, 0, 2, 7,14, 9,15,13,15,15,14,15, 7, 0,15,12, 9,12, 7,11, 7, 8, &
     9, 7,13, 7,15, 8, 0, 2, 7,23, 7, 7, 0, 7,13, 9,15,13,15,15,13,15, &
     9,13, 7, 9, 7, 7, 9, 0, 1,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, &
     7,15, 9, 0, 2,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, 7,15, 9, 0, &
    15,23,15, 7, 0, 1, 7,11,15,11,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7 /

  data ( ifont(i), i = 1189, 1584 ) / &
    13, 7,15, 9, 0, 3, 9, 7, 9,23,13,23,13,22, 0, 8,15,12,15, 0, 8, 7, &
    11, 7, 0, 2,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, 7,15, 9, 0,15, &
    13,15, 3,13, 1, 9, 1, 7, 3, 0, 2, 7, 7, 7,23, 0, 7,14, 9,15,13,15, &
    15,14,15, 7, 0, 3, 9,15,11,15,11, 7, 0, 9, 7,13, 7, 0, 9,17, 9,19, &
    11,19,11,17, 9,17, 0, 2, 9,15,11,15,11, 1, 7, 1, 7, 3, 0, 9,17,11, &
    17,11,19, 9,19, 9,17, 0, 3, 7, 7, 7,23, 0,15,15, 7,10, 0, 9,11,15, &
     7, 0, 2, 9,23,11,23,11, 7, 0, 9, 7,13, 7, 0, 3, 7,15, 7, 7, 0, 7, &
    14, 8,15,10,15,11,14,11, 7, 0,11,14,12,15,14,15,15,14,15, 7, 0, 2, &
     7, 7, 7,15, 0, 7,14, 9,15,13,15,15,14,15, 7, 0, 1, 7,13, 9,15,13, &
    15,15,13,15, 9,13, 7, 9, 7, 7, 9, 7,13, 0, 2, 7,13, 9,15,13,15,15, &
    13,15, 9,13, 7, 9, 7, 7, 9, 0, 7,14, 7, 1, 0, 2,15,13,13,15, 9,15, &
     7,13, 7, 9, 9, 7,13, 7,15, 9, 0,15,14,15, 1, 0, 2, 7,15, 9,15, 9, &
     7, 0, 9,13,11,15,13,15,15,13, 0, 1,15,13,13,15, 9,15, 7,13, 9,11, &
    13,11,15, 9,13, 7, 9, 7, 7, 9, 0, 2, 9,23, 9, 7,11, 7, 0, 7,17,11, &
    17, 0, 2, 7,15, 7, 9, 9, 7,13, 7,15, 9, 0,15,15,15, 7, 0, 1, 7,15, &
    11, 7,15,15, 0, 1, 7,15, 9, 7,11,11,13, 7,15,15, 0, 2, 7,15,15, 7, &
     0, 7, 7,15,15, 0, 2, 7,15,11, 7, 0,15,15,10, 5, 7, 1, 0, 1, 7,15, &
    15,15, 7, 7,15, 7, 0, 1,11,23, 7,23, 9,17, 7,15, 9,13, 7, 7,11, 7 /

  data ( ifont(i), i = 1585, 1617 ) / &
     0, 1, 9,23, 9, 7, 0, 1, 7,23,11,23, 9,17,11,15, 9,13,11, 7, 7, 7, &
     0, 1, 5,21, 7,23,15,21,17,23, 0 /
!
  nchar = len_trim ( s )

  if ( pwide <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S_PLOT - Serious error!'
    write ( *, * ) '  The plot width PWIDE is negative!'
    write ( *, * ) '  PWIDE = ', pwide
    return
  end if
!
!  Chop titles that are too long.  To do this, we need to know the
!  width of the plot (PWIDE) in same units as CWIDE.
!
  nmax = int ( pwide / cwide )

  if ( nchar > nmax ) then
    nchar = nmax
  end if
!
!  Shift string if centering or right flush option used.
!
  if ( flush(1:1) == 'l' .or. flush(1:1) == 'L' ) then
    xcopy = x
    ycopy = y
  else if ( flush(1:1) == 'c' .or. flush(1:1) == 'C' ) then
    xcopy = x - 0.5E+00 * nchar * cwide * cos ( degrees_to_radians ( angle ) )
    ycopy = y - 0.5E+00 * nchar * cwide * sin ( degrees_to_radians ( angle ) )
  else if ( flush(1:1) == 'r' .or. flush(1:1) == 'R' ) then
    xcopy = x - nchar * cwide * cos ( degrees_to_radians ( angle ) )
    ycopy = y - nchar * cwide * sin ( degrees_to_radians ( angle ) )
  else
    xcopy = x
    ycopy = y
  end if
!
!  Note that screen coordinates are used.
!  Thus a width of 0.1 is intended to mean 1/10 of screen size.
!
!  Set the scale factor for character height.
!
  scl2 = cwide / 16.0E+00
!
!  Set the starting point for the line of text, the lower left
!  corner of the first character.
!
!  Set the origin about which rotation is performed.
!
  xb = xcopy
  xrot = xcopy
  yb = ycopy
  yrot = ycopy
!
!  Get trig functions if rotation required, converting from
!  degrees to radians.
!
  if ( angle == 0.0E+00 ) then
    rotate = .false.
  else
    ca = cos ( degrees_to_radians ( angle ) )
    sa = sin ( degrees_to_radians ( angle ) )
    rotate = .true.
  end if
!
!  Loop over all characters in the string.
!
  do icr = 1, nchar

    xold = x
    yold = y
    xnew = x
    ynew = y
!
!  Get the ASCII code for the character and shift by 31 so that
!  the first printable character becomes code 1.
!
    c = s(icr:icr)
    iascii = ichar ( c ) - 31
!
!  Replace any nonprintable characters with blanks.
!
    if ( iascii < 1 .or. iascii > 95 ) then
      iascii = 1
    end if
!
!  Get the pointer to this character in font table.
!
    ip = ipoint(iascii)
!
!  Get the number of "vectors" required to draw the character.
!  Here "vectors" means the number of times the pen is lowered, not
!  the number of pen strokes.
!
!  For blanks, this number is 1, due to the way the
!  algorithm is coded.
!
    nvec = ifont(ip)
!
!  Loop over all required pen movements.
!
    do iv = 1, nvec

      ipen = 3
      ip = ip + 1

      do while ( ifont(ip) /= 0 )

        xc = xb + scl2 * ( ifont(ip) - 1 )
        yc = yb + scl2 * ( ifont(ip+1) - 7 )
!
!  Apply rotation if necessary.
!
        if ( rotate ) then
          xt = xc - xrot
          yt = yc - yrot
          xc = xrot + ca * xt - sa * yt
          yc = yrot + sa * xt + ca * yt
        end if
!
!  Plot the pen stroke.
!
        if ( ipen == 3 ) then
          xnew = xc
          ynew = yc
        else
          xold = xnew
          yold = ynew
          xnew = xc
          ynew = yc
          call movcgm ( xold, yold )
          call drwcgm ( xnew, ynew )
        end if

        ipen = 2
        ip = ip + 2

      end do

    end do
!
!  Advance the base to compensate for character just drawn.
!
    xb = xb + cwide

  end do

  return
end
