program main

!*****************************************************************************80
!
!! DISLIN_EX04 demonstrates various interpolation methods for data.
!
!  Discussion:
!
!    Create a plot containing six subplots.  Each subplot shows the same
!    data, interpolated in a different way.
!
!  Modified:
!
!    12 April 2011
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
  integer :: i
  integer :: nx
  integer :: ny
  integer :: nya = 2700
  character ( len = 6 ), dimension (6) :: cpol = (/ &
    'Spline', 'Stem  ', 'Bars  ', 'Step  ', 'Stairs', 'Linear' /)

  real, dimension (16) :: x = (/ &
     0.0,  1.0,  3.0,  4.5,  6.0,  8.0, 9.0, 11.0, 12.0, 12.5, &
    13.0, 15.0, 16.0, 17.0, 19.0, 20.0 /)
  real, dimension (16) :: y = (/ &
     2.0,  4.0,  4.5,  3.0,  1.0,  7.0, 2.0,  3.0,  5.0,  2.0, &
     2.5,  2.0,  4.0,  6.0,  5.5,  4.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX04:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the various interpolation'
  write ( *, '(a)' ) '  methods available for (X,Y) data.'
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
  call setfil ( 'dislin_ex04.png' )
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
!  Use the COMPLEX font.
!
  call complx ( )
!
!  Mark data points on the curve, incrementing by 1.
!
  call incmrk ( 1 )
!
!  Use a plot symbol height measured in plot coordinates.
!
  call hsymbl ( 25 )
!
!  Define axis system titles.
!
  ctit = 'Interpolation methods'
  call titlin ( ctit, 1 )
!
!  Define the X and Y sizes of the axis system.
!
  call axslen ( 1500, 350 )
!
!  Specify how the lower X, left Y, upper X and right Y axes are labeled.
!
  call setgrf ( 'line', 'line', 'line', 'line' )
!
!  The subplots are drawn from bottom to top.
!
  do i = 1, 6

    call axspos ( 350, nya-(i-1)*350 )
    call polcrv ( cpol(i) )
    call marker ( 0 )

    call graf ( 0.0, 20.0, 0.0, 5.0, 0.0, 10.0, 0.0, 5.0 )
    nx = nxposn ( 1.0 )
    ny = nyposn ( 8.0 )
    call messag ( cpol(i), nx, ny )
    call curve ( x, y, 16 )

    if ( i == 6 ) then
      call height ( 50 )
      call title ( )
    end if

    call endgrf ( )

  end do
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX04:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
