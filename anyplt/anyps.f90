subroutine anyplt ( icom )

!*****************************************************************************80
!
!! ANYPLT is a generic graphics interface routine.
!
!  Discussion:
!
!    ANYPLT is a subroutine which provides a simple, standard interface
!    between FORTRAN programs and various output devices.  To run a
!    program which calls ANYPLT on a different machine, the program
!    is not modified in any way, but a different version of the ANYPLT
!    program is provided.  
!
!    The following versions are available:
!
!    ANYATT - AT&T PC6300 graphics (640 by 400).  Requires ATTPLT.ASM.
!    ANYBUG - Simple debugging output to a file.  Nominal 1.0 by 1.0 plot.
!    ANYCAL - CALCOMP file output.  Available on many mainframes.
!             8.5 inches by 11.0 inches
!    ANYIBM - IBM PC hi resolution (640 by 200).  Requires IBMPLT.ASM.
!    ANYMAC - Macintosh graphics.  Requires auxilliary routine TOOLBX.SUB.
!             (342 high, 512 wide)
!    ANYNCR - NCAR graphics package.
!    ANYNUL - Does nothing.
!    ANYP10 - PLOT10 interactive graphics. (1024 by 768)
!    ANYTTY - Simple 'typewriter' graphics (80 by 24 "pixels")
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
!    Input, integer ( kind = 4 ) ICOM, the index of the graphics command.
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
!    14, draw an arrow at (XPLT1,YPLT1), of length YPLT2 and angle XPLT2.
!
  implicit none

  character ( len = 80 ) carray
  real degrees_to_radians
  character ( len = 80 ) file_name
  real font_size
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  integer ( kind = 4 ), save :: iunit = 0
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ) marray
  integer ( kind = 4 ), save :: nplot = 0
  real xdraw(6)
  real, save :: xmax = 0.0E+00
  real, save :: xmin = 1.0E+00
  real xplt1
  real xplt2
  real xtip
  real ydraw(6)
  real, save :: ymax = 0.0E+00
  real, save :: ymin = 1.0E+00
  real yplt1
  real yplt2
  real ytip

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2
  common /anychr/ carray
!
!  ICOM = 0  Enable graphics
!
  if ( icom == 0 ) then

    call get_unit ( iunit )

    file_name = 'anyps.ps'

    call ps_file_open ( file_name, iunit, ierror )

    nplot = 0

    call ps_file_head ( file_name )
!
!  ICOM = 1  Disable graphics
!
  else if ( icom == 1 ) then

    call ps_file_tail ( )

    call ps_file_close ( iunit )

    nplot = 0
!
!  ICOM = 2  Begin plot
!
  else if ( icom == 2 ) then

    nplot = nplot + 1
!
!  ICOM = 3  Define plot size
!
  else if ( icom == 3 ) then

    xmin = xplt1
    call ps_setting_real ( 'SET','XMIN', xmin )
    xmax = xplt2
    call ps_setting_real ( 'SET','XMAX', xmax )
    ymin = yplt1
    call ps_setting_real ( 'SET','YMIN', ymin )
    ymax = yplt2
    call ps_setting_real ( 'SET','YMAX', ymax )

    call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  ICOM = 4  Move to point
!
  else if ( icom == 4 ) then

    call ps_moveto ( xplt1, yplt1 )
!
!  ICOM = 5  Draw to point
!
  else if ( icom == 5 ) then

    call ps_lineto ( xplt1, yplt1 )
!
!  ICOM = 6  Clear screen
!
  else if ( icom == 6 ) then
!
!  ICOM = 7,  Write string at position.
!  Need a way to specify the font size.
!
  else if ( icom == 7 ) then

    font_size = 10.0 * xplt2
    call ps_font_size ( font_size )

    call ps_moveto ( xplt1, yplt1 )
    call ps_label ( carray )
!
!  ICOM = 8  Use virtual cursor
!
  else if ( icom == 8 ) then
!
!  ICOM = 9  End plot
!
  else if ( icom == 9 ) then

    call ps_page_tail ( )
!
!  ICOM = 10  Ring bell
!
  else if ( icom == 10 ) then
!
!  ICOM = 11  Mark data
!
  else if ( icom == 11 ) then

    call ps_mark_circle ( xplt1, yplt1 )
!
!  ICOM = 12  Return screen data
!
  else if ( icom == 12 ) then

    call ps_setting_real ( 'GET', 'XMIN', xmin )
    xplt1 = xmin

    call ps_setting_real ( 'GET', 'XMAX', xmax )
    xplt2 = xmax

    call ps_setting_real ( 'GET', 'YMIN', ymin )
    yplt1 = ymin

    call ps_setting_real ( 'GET', 'YMAX', ymax )
    yplt2 = ymax
!
!  ICOM = 13  Return version
!
  else if ( icom == 13 ) then

    carray = 'AnyPlt - Version 1.01  PostScript Graphics  PSPLT'
!
!  ICOM = 14, Draw an arrow.
!
  else if ( icom == 14 ) then

    xtip = xplt1 + yplt2 * cos ( xplt2 )
    ytip = yplt1 + yplt2 * sin ( xplt2 )

    call arrow ( xplt1, yplt1, xtip, ytip, xdraw, ydraw )

    do i = 1, 5
      call ps_moveto ( xdraw(i), ydraw(i) )
      call ps_lineto ( xdraw(i+1), ydraw(i+1) )
    end do    
!
!  Unknown value of ICOM.
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'PSPLT - Fatal error!'
    write ( *, * ) '  Unknown value of ICOM = ', icom
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
  integer ( kind = 4 ) ndraw
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
