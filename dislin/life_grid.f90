program main

!*****************************************************************************80
!
!! LIFE_GRID uses DISLIN to draw a grid for the game of Life.
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

  character ( len = 60 ) ctit
  integer i
  integer j
  integer nr
  integer nx
  integer ny
  integer pat
  real r
  real x
  real xvec(2)
  real y
  real yvec(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LIFE_GRID:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Use DISLIN routines to plot a grid for Life.'
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
  call setfil ( 'life_grid.png' )
!
!  Choose the page size and orientation.
!  'USA' is 2160 plot units wide and 2790 plot units high.
!  'P' requests PROFILE rather than LANDSCAPE orientation.
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
!  Use the COMPLEX font.
!
  call complx ( )
!
!  Use a color table, which is required if we want to do color graphics.
!  For this color table, in particular,
!  1 = black,
!  2 = red,
!  3 = green,
!  4 = blue.
!
  call setvlt ( 'small' )
!
!  Define the X and Y sizes of the axis system in plot units.
!
  call axslen ( 1000, 1000 )
!
!  Specify how the lower X, left Y, upper X and right Y axes are labeled.
!
  call setgrf ( 'line', 'line', 'line', 'line' )
!
!  Set the axis origin 500 plot units to the right, and 1500 plot units DOWN.
!
  call axspos ( 500, 1500 )
!
!  Relate the physical coordinates to the axes.
!
  call graf ( 0.0, 100.0, 0.0, 0.5, 0.0, 100.0, 0.0, 0.5 )
!
!  Draw 21 horizontal lines.
!
  do j = 0, 100, 5
    y = real ( j )
    xvec(1) = 0.0
    xvec(2) = 100.0
    yvec(1) = y
    yvec(2) = y
    call curve ( xvec, yvec, 2 )
  end do
!
!  Draw 21 vertical lines.
!
  do i = 0, 100, 5
    x = real ( i )
    xvec(1) = x
    xvec(2) = x
    yvec(1) = 0.0
    yvec(2) = 100.0
    call curve ( xvec, yvec, 2 )
  end do
!
!  Select the shading pattern.
!
  pat = 16
  call shdpat ( pat )
!
!  Select color 3 (green) from the color table.
!
  call setclr ( 3 )
!
!  Draw one circle near the origin.
!
  x = 2.5
  y = 2.5
  r = 2.0
  call rlcirc ( x, y, r )
!
!  Select color 2 (red).
!
  call setclr ( 2 )
!
!  Draw a glider.
!
  x = 7.5
  y = 37.5
  r = 2.0
  call rlcirc ( x, y, r )

  x = 12.5
  y = 37.5
  r = 2.0
  call rlcirc ( x, y, r )

  x = 12.5
  y = 47.5
  r = 2.0
  call rlcirc ( x, y, r )

  x = 17.5
  y = 37.5
  r = 2.0
  call rlcirc ( x, y, r )

  x = 17.5
  y = 42.5
  r = 2.0
  call rlcirc ( x, y, r )
!
!  Select color 4 (blue)
!
  call setclr ( 4 )
!
!  Select open shading pattern.
!
  pat = 0
  call shdpat ( pat )
!
!  Draw three open circles.
!
  x = 62.5
  y = 62.5
  r = 2.0
  call rlcirc ( x, y, r )

  x = 67.5
  y = 57.5
  r = 2.0
  call rlcirc ( x, y, r )

  x = 72.5
  y = 52.5
  r = 2.0
  call rlcirc ( x, y, r )
!
!  Select character height in plot units.
!
  call height ( 50 )
!
!  Select color 1 (black) from the color table.
!
  call setclr ( 1 )
!
!  Define axis system titles.
!
  ctit = 'Grid for Game of Life'
  call titlin ( ctit, 1 )
!
!  Draw the title.
!
  call title ( )
!
!  End this plot.
!
  call endgrf ( )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LIFE_GRID:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
