program main

!*****************************************************************************80
!
!! SCATTER_PLOT uses DISLIN to draw a scatter plot of X Y data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2011
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

  integer, parameter :: n = 500

  integer i
  integer j
  integer nr
  integer nx
  integer ny
  integer pat
  real r
  real r4_uniform_01
  real s
  integer seed
  real x
  real xvec(n)
  real y
  real yvec(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SCATTER_PLOT:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Use DISLIN routines to make a scatterplot.'
!
!  Generate the data.  
!  We average 4 random values to get data that tends to cluster
!  near (0.5,0.5).
!
  seed = 123456789

  do i = 1, n
    s = 0.0E+00
    do j = 1, 4
      s = s + r4_uniform_01 ( seed )
    end do
    xvec(i) = s / 4.0E+00
  end do

  do i = 1, n
    s = 0.0E+00
    do j = 1, 4
      s = s + r4_uniform_01 ( seed )
    end do
    yvec(i) = s / 4.0E+00
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
  call setfil ( 'scatter_plot.png' )
!
!  Choose the page size and orientation.
!  'USA' is 2160 plot units wide and 2790 plot units high.
!  'P' requests PROFILE rather than LANDSCAPE orientation.
!
  call setpag ( 'usap' )
!
!  For PNG output, use reverse the default black background to white.
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
!  Define the X and Y sizes of the axis system in plot units.
!
  call axslen ( 1800, 1800 )
!
!  Specify how the lower X, left Y, upper X and right Y axes are labeled.
!
  call setgrf ( 'line', 'line', 'line', 'line' )
!
!  Set the axis origin 180 plot units to the right, and 2610 plot units DOWN.
!
  call axspos ( 180, 2610 )
!
!  Relate the physical coordinates to the axes.
!
  call graf ( 0.0, 1.0, 0.0, 0.1, 0.0, 1.0, 0.0, 0.1 )
!
!  Add a grid, with one grid line for every tick mark in the X and Y axes.
!
  call grid ( 1, 1 )
!
!  Select the shading pattern.
!
  pat = 16
  call shdpat ( pat )
!
!  Set the color to blue.
!
  call color ( "blue" )
!
!  At every data point, draw a circle of radius 0.01.
!
  do i = 1, n
    call rlcirc ( xvec(i), yvec(i), 0.01 )
  end do
!
!  Select character height in plot units.
!
  call height ( 50 )
!
!  We choose "white" for the title, which is actually black
!  because we reversed black and white earlier!
!
  call color ( "white" )
!
!  Define axis system titles.
!
  call titlin ( 'Scatter Plot', 1 )
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
  write ( *, '(a)' ) 'SCATTER_PLOT:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
