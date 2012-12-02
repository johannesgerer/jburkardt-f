program main

!*****************************************************************************80
!
!! DISLIN_EX09 demonstrates 3D color plots.
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

  integer, parameter :: n = 100

  real :: fpi
  integer :: i
  integer :: j
  real :: step
  real :: x
  real :: y
  real, dimension (n,n) :: zmat

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX09:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of a '
  write ( *, '(a)' ) '  3D color plot.'

  fpi = 3.1415927 / 180.0
  step = 360.0 / real ( n - 1 )

  do i = 1, n
    x = real ( i - 1 ) * step
    do j = 1, n
      y = real ( j - 1 ) * step
      zmat(i,j) = 2.0 * sin ( x * fpi ) * sin ( y * fpi )
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
  call setfil ( 'dislin_ex09.png' )
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
!  Use the standard HARDWARE font.
!
  call hwfont ( )

  call titlin ( '3-d color plot of the function', 2 )
  call titlin ( 'f(x,y) = 2 * sin(x) * sin(y)', 4 )

  call name ( 'x-axis', 'x' )
  call name ( 'y-axis', 'y' )
  call name ( 'z-axis', 'z' )

  call intax ( )
  call autres ( n, n )
  call axspos ( 300, 1850 )
  call ax3len ( 2200, 1400, 1400 )

  call graf3 ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0, -2.0, 2.0, -2.0, 1.0 )
  call crvmat ( zmat, n, n, 1, 1 )

  call height ( 50 )
  call title ( )
  call mpaepl ( 3 )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX09:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
