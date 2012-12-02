program main

!*****************************************************************************80
!
!! QUICKPLOT_SURFACE demonstrates the DISLIN quickplot command QPLSUR.
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
  real, parameter :: pi = 3.1415926
  real x
  real y
  real zmat(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_SURFACE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command QPLSUR'
  write ( *, '(a)' ) '  to plot a surface Z(X,Y) stored as a matrix.'
!
!  Set up the data.
!
  do i = 1, m
    x = 2.0 * pi * real ( i - 1 ) / real ( m - 1 )
    do j = 1, n
      y = 2.0 * pi * real ( j - 1 ) / real ( n - 1 )
      zmat(i,j) = 2.0 * sin ( x ) * sin ( y )
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
  call setfil ( 'quickplot_surface.png' )
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
  call titlin ( 'Quick plot by QPLSUR', 2 )
!
!  Draw the curve.
!
  call qplsur ( zmat, m, n )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_SURFACE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
