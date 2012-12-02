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
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ), parameter :: ldunit = 3
  integer ( kind = 4 ) marray
  real xplt1
  real xplt2
  real yplt1
  real yplt2

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2
  common /anychr/ carray
!
!  ICOM = 0  Enable graphics
!
  if ( icom == 0 ) then

    open ( unit = ldunit, file = 'anyplt.bug', status = 'replace' )

    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 00 - Enable graphics.'
    write ( ldunit, '(a,g14.6,a,g14.6)' ) '  Xmin = ', xplt1, ' Xmax = ', xplt2
    write ( ldunit, '(a,g14.6,a,g14.6)' ) '  Ymin = ', yplt1, ' Ymax = ', yplt2
!
!  ICOM = 1  Disable graphics
!
  else if ( icom == 1 ) then

    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 01 - Disable graphics.'
    close ( unit = ldunit )
!
!  ICOM = 2  Begin plot
!
  else if ( icom == 2 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 02 - Begin plot'
!
!  ICOM = 3  Define plot size
!
  else if ( icom == 3 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 03 - Define plot limits:'
    write ( ldunit, '(a,g14.6,a,g14.6)' ) '  Xmin = ', xplt1, ' Xmax = ', xplt2
    write ( ldunit, '(a,g14.6,a,g14.6)' ) '  Ymin = ', yplt1, ' Ymax = ', yplt2
!
!  ICOM = 4  Move to point
!
  else if ( icom == 4 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 04 - Move to the point:'
    write ( ldunit, '(a,g14.6,a,g14.6,a)' )  '  (', xplt1, ', ', yplt1, ')'
!
!  ICOM = 5  Draw to point
!
  else if ( icom == 5 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 05 - Draw to the point:'
    write ( ldunit, '(a,g14.6,a,g14.6,a)' ) '  (', xplt1, ', ', yplt1, ')'
!
!  ICOM = 6  Clear screen
!
  else if ( icom == 6 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 06 - Clear the screen.'
!
!  ICOM = 7,  Write string at position
!
  else if ( icom == 7 ) then

    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 07 - Write a string.'
    write ( ldunit, '(a,g14.6,a,g14.6,a)' ) '  At the point (', xplt1, ', ', yplt1, ')'
    write ( ldunit, '(a)' ) '  Write the characters: ' // trim ( carray )
!
!  ICOM = 8  Use virtual cursor
!
  else if ( icom == 8 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 08 - Use the virtual cursor.'
!
!  ICOM = 9  End plot
!
  else if ( icom == 9 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 09 - End this plot.'
!
!  ICOM = 10  Ring bell
!
  else if ( icom == 10 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 10 - Ring the bell.'
!
!  ICOM = 11  Mark data
!
  else if ( icom == 11 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 11 - Mark data point.'
    write ( ldunit, '(a,g14.6,a,g14.6,a)' ) '  At the point (', xplt1, ', ', yplt1, ')'
    write ( ldunit, '(a)' ) '  Mark the data with "' // carray(1:1) // '"'
!
!  ICOM = 12  Return screen data
!
  else if ( icom == 12 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 12 - Return maximum screen data.'
!
!  ICOM = 13  Return version
!
  else if ( icom == 13 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'Bugplt - 13 - Return version.'
    carray = 'AnyPlt - Version 1.03  21 November 2000  BugPlt'
    write ( ldunit, '(a)' ) carray
!
!  ICOM = 14, Draw an arrow.
!
  else if ( icom == 14 ) then
    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - 14 - Arrow.'
    write ( ldunit, '(a,g14.6,a,g14.6,a)' ) '  At the point (', xplt1, ', ', yplt1, ')'
    write ( ldunit, * ) '  draw an arrow of length ', yplt2
    write ( ldunit, * ) '  in direction angle ', xplt2
!
!  Unknown value of ICOM.
!
  else

    write ( ldunit, '(a)' ) ' '
    write ( ldunit, '(a)' ) 'BugPlt - Fatal error!'
    write ( ldunit, '(a,i8)' ) '  Unknown value of ICOM = ', icom
    close ( unit = ldunit )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BugPlt - Fatal error!'
    write ( *, '(a,i8)' ) '  Unknown value of ICOM = ', icom
    stop
  end if

  return
end
