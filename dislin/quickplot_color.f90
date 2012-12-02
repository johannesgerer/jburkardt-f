program main

!*****************************************************************************80
!
!! QUICKPLOT_COLOR demonstrates the DISLIN quickplot command QPLCLR.
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

  integer, parameter :: n = 100

  real fpi
  integer i
  integer j
  real, parameter :: pi = 3.1415926  
  real step
  real x
  real y
  real z(n,n)

  fpi = pi / 180.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_COLOR:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command QPLCLR'
  write ( *, '(a)' ) '  to make a color plot of a matrix of data.'
!
!  Set up the X and Y data for the plot.
!
  step = 360.0 / real ( n - 1 )
  do i = 1, n
    x = real ( i - 1 ) * step
    do j = 1, n
      y = real ( j - 1 ) * step
      z(i,j) = 2.0 * sin ( x * fpi ) * sin ( y * fpi )
    end do
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
  call setfil ( 'quickplot_color.png' )
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
  call name ( "X-axis", "x" );
  call name ( "Y-axis", "y" );
  call titlin ( 'Quick plot by QPLCLR', 2 )
!
!  Draw the curve.
!
  call qplclr ( z, n, n )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_COLOR:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
