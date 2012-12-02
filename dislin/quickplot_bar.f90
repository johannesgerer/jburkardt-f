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

  integer, parameter :: n = 14

  integer i
  real, dimension ( n ) :: xray = (/ &
     1.0, 15.0, 38.0, 22.0, 16.0, &
    16.0, 26.0, 55.0, 50.0, 40.0, &
    16.0,  3.0,  0.0,  1.0 /)
    
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_BAR:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command QPLBAR'
  write ( *, '(a)' ) '  to plot a bar chart.'
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
  call setfil ( 'quickplot_bar.png' )
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
  call name ( '<-- Minutes -->', 'X' )
  call name ( '<-- Frequency -->', 'Y' )
  call titlin ( 'Quick plot by QPLBAR', 2 )
!
!  Draw the curve.
!
  call qplbar ( xray, n )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_BAR:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
