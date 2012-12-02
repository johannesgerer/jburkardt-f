program main

!*****************************************************************************80
!
!! QUICKPLOT_CURVE demonstrates the DISLIN quickplot command QPLOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2011
!
!  Author:
!
!    John Burkardt
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

  integer, parameter :: n = 100

  integer i
  real, parameter :: pi = 3.1415926  
  real xray(n)
  real yray(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_CURVE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command QPLOT'
  write ( *, '(a)' ) '  to plot a curve.'
!
!  Set up the X and Y data for the plot.
!
  do i = 1, n
    xray(i) = real (  i - 1 ) * 360.0 / real ( n - 1 )
  end do

  yray(1:n) = sin ( pi * xray(1:n) / 180.0 )
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
  call setfil ( 'quickplot_curve.png' )
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
!  Label the axes and the plot.
!
  call name ( '<-- Angle in Degrees -->', 'X' )
  call name ( '<-- Sine (angle) -->', 'Y' )
  call titlin ( 'Quick plot by QPLOT', 2 )
!
!  Draw the curve.
!
  call qplot ( xray, yray, n )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_CURVE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
