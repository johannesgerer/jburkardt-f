program main

!*****************************************************************************80
!
!! DISLIN_EX01 demonstrates the use of CURVE to plot (X,Y) data.
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

  integer, parameter :: n=100

  real :: fpi
  integer :: i
  real, parameter :: pi = 3.1415926  
  real :: x
  real, dimension  (n) :: xray
  real, dimension ( n ) :: y1ray
  real, dimension ( n ) :: y2ray

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX01:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the use of CURVE to plot '
  write ( *, '(a)' ) '  (X,Y) data.'
!
!  Set up the X and Y data for the plot.
!
  do i = 1, n
    xray(i) = real (  i - 1 ) * 360.0 / real ( n - 1 )
  end do

  y1ray(1:n) = sin ( pi * xray(1:n) / 180.0 )
  y2ray(1:n) = cos ( pi * xray(1:n) / 180.0 )
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
  call setfil ( 'dislin_ex01.png' )
!
!  Choose the page size and orientation.
!
  call setpag ( 'usal' )
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
!  Use the COMPLEX font.
!
  call complx ( )
  call axspos ( 450, 1800 )
  call axslen ( 2200, 1200 )

  call name ( 'X-axis', 'X' )
  call name ( 'Y-axis', 'Y' )

  call labdig ( -1, 'x' )
  call ticks ( 10, 'xy' )

  call titlin ( 'Demonstration of CURVE', 1 )
  call titlin ( 'sin(x), cos(x)', 3 )

  call graf ( 0.0, 360.0, 0.0, 90.0, -1.0, 1.0, -1.0, 0.5 )
  call title ( )
!
!  Draw XRAY versus Y1RAY in red.
!
  call color ( 'red' )
  call curve ( xray, y1ray, n )
!
!  Draw XRAY versus Y2RAY in green.
!
  call color ( 'green' )
  call curve ( xray, y2ray, n )

  call color ( 'fore' ) 
  call dash ( )
  call xaxgit ( )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX01:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
