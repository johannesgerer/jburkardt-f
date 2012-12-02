program main

!*****************************************************************************80
!
!! DISLIN_EX07 demonstrates demonstrates 3D bar graphs and pie charts.
!
!  Modified:
!
!    09 April 2011
!
!  Reference:
!
!    Helmut Michels,
!    The Data Plotting Software DISLIN - version 10.4,
!    Shaker Media GmbH, January 2010,
!    ISBN13: 978-3-86858-517-9.
!
  use dislin

  implicit none

  character ( len = 80 ) :: cbuf
  integer, dimension (5) :: ic1ray = (/ 50, 150, 100, 200, 175 /)
  integer, dimension (5) :: ic2ray = (/ 50, 150, 100, 200, 175 /)
  real, dimension (5) :: xray  = (/ 2.0, 4.0, 6.0, 8.0, 10.0 /)
  real, dimension (5) :: y1ray = (/ 0.0, 0.0, 0.0, 0.0, 0.0 /)
  real, dimension (5) :: y2ray = (/ 3.2, 1.5, 2.0, 1.0, 3.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX07:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of 3D bar '
  write ( *, '(a)' ) '  and pie charts.'
!
!  Specify the format of the output file.
!
  call metafl ( 'png' )
!
!  Indicate that new data overwrites old data.
!
  call filmod ( 'delete' )
!
!  Specify the name of the output graphics file.
!
  call setfil ( 'dislin_ex07.png' )
!
!  Choose the page size and orientation.
!
  call setpag ( 'usap' )
!
!  For PNG output, reverse the default black background to white.
!
  call scrmod ( 'reverse' )
!
!  Open DISLIN.
!
  call disini ( )
!
!  Plot a border around the page.
!
  call pagera ( )
!
!  Use the standard HARDWARE font.
!
  call hwfont ( )

  call titlin ( '3-d bar graph / 3-d pie chart', 2 )
  call htitle ( 40 )

  call shdpat ( 16)
  call axslen ( 1500, 1000 )
  call axspos ( 300, 1400 )

  call barwth ( 0.5 )
  call bartyp ( '3dvert' )
  call labels ( 'second', 'bars' )
  call labpos ( 'outside', 'bars' )
  call labclr ( 255, 'bars' )
  call graf ( 0.0, 12.0, 0.0, 2.0, 0.0, 5.0, 0.0, 1.0 )
  call title ( )
  call color( 'red' )
  call bars ( xray, y1ray, y2ray, 5 )
  call endgrf ( )

  call shdpat ( 16 )
  call labels ( 'data', 'pie' )
  call labclr ( 255, 'pie' )
  call chnpie ( 'none' )
  call pieclr ( ic1ray, ic2ray, 5 )
  call pietyp ( '3d' )
  call axspos ( 300, 2700 )
  call piegrf ( cbuf, 0, y2ray, 5 )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX07:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end

