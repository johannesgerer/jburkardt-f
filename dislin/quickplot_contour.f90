program main

!*****************************************************************************80
!
!! QUICKPLOT_CONTOUR demonstrates the DISLIN quickplot command QPLCON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2012
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

  integer, parameter :: m = 100
  integer, parameter :: n = 100

  integer i
  integer j
  integer levels
  real x
  real y
  real zmat(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_CONTOUR:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command QPLCON'
  write ( *, '(a)' ) '  to make a contour plot of data stored as a matrix.'
!
!  Set up the data.
!
  do i = 1, m
    x = 1.6 * real ( i - 1 ) / real ( m - 1 )
    do j = 1, n
      y = 1.6 * real ( j - 1 ) / real ( n - 1 )
      zmat(i,j) = ( x * x - 1.0 )**2 + ( y * y - 1.0 )**2
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
  call setfil ( 'quickplot_contour.png' )
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
  call name ( '<-- X -->', 'X' )
  call name ( '<-- Y -->', 'Y' )
  call titlin ( 'Quick plot by QPLCON', 2 )
!
!  Draw the curve.
!
  levels = 20
  call qplcon ( zmat, m, n, levels )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_CONTOUR:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
