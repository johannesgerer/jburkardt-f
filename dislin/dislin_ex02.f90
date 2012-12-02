program main

!*****************************************************************************80
!
!! DISLIN_EX02 demonstrates the use of POLAR to plot (R,Theta) data.
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

  integer, parameter :: m = 10
  integer, parameter :: n = 300

  real a
  integer i
  real, parameter :: pi = 3.1415927
  real step
  real x2(m)
  real xray(n)
  real y2(m)
  real yray(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX02:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the use of POLAR to plot '
  write ( *, '(a)' ) '  (R,Theta) data.'

  step = 360.0 / real ( n - 1 )

  do i = 1, n
    a = real ( i - 1 ) * step
    a = a * pi / 180.0
    yray(i) = a
    xray(i) = sin ( 5.0 * a )
  end do

  do i = 1, m
    x2(i) = real ( i )
    y2(i) = real ( i )
  end do
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
  call setfil ( 'dislin_ex02.png' )
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

  call titlin ( 'Polar Plots', 2 )
  call ticks ( 3, 'Y' )
  call axends ( 'NOENDS', 'X' )
  call labdig ( -1, 'Y' )
  call axslen ( 1000, 1000 )
  call axsorg ( 1050, 900 )

  call polar ( 1.0, 0.0, 0.2, 0.0, 30.0 )
  call curve ( xray, yray, n )
  call htitle ( 50 )
  call title ( )
  call endgrf ( )

  call labdig ( -1, 'X' )
  call axsorg ( 1050, 2250 )
  call labtyp ( 'VERT', 'Y' )
  call barwth ( 5.0 )
  
  call polar ( 10.0, 0.0, 2.0, 0.0, 30.0 )
  call barwth ( -5.0 )
  call polcrv ( 'FBARS' )
  call curve ( x2, y2, m )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX02:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
