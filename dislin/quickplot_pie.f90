program main

!*****************************************************************************80
!
!! QUICKPLOT_PIE demonstrates the DISLIN quickplot command QPLPIE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2012
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

  integer, parameter :: n = 5

  integer i
  real, dimension (5) :: xray = (/ 10.0, 20.0, 15.0, 5.0, 50.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_PIE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command QPLPIE'
  write ( *, '(a)' ) '  to plot a pie chart.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we plot 10 percent luck, 20 percent skill,'
  write ( *, '(a)' ) '  15 percent concentrated power of will, 5 percent pleasure,'
  write ( *, '(a)' ) '  50 percent pain.'
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
  call setfil ( 'quickplot_pie.png' )
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
!  Label the plot.
!
  call titlin ( 'Quick plot by QPLPIE', 2 )
!
!  Draw the pie chart.
!
  call qplpie ( xray, n )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_PIE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
