program main

!*****************************************************************************80
!
!! DISLIN_EX10 demonstrates a 3D surface plot.
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

  integer, parameter :: n = 50

  character ( len = 60 ) :: ctit1 = 'surface plot (surmat)'
  character ( len = 60 ) :: ctit2 = 'f(x,y) = 2*sin(x)*sin(y)'
  real :: fpi
  integer :: i
  integer :: j
  real :: step
  real :: x
  real :: y
  real, dimension (n,n) :: zmat

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX10:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of '
  write ( *, '(a)' ) '  a surface plot.'

  fpi = 3.14159 / 180.0
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
  call setfil ( 'dislin_ex10.png' )
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
  call axspos ( 200, 2600 )
  call axslen ( 1800, 1800 )

  call name ( 'x-axis', 'x' )
  call name ( 'y-axis', 'y' )
  call name ( 'z-axis', 'z' )

  call titlin ( ctit1, 2 )
  call titlin ( ctit2, 4 )

  call view3d ( -5.0, -5.0, 4.0, 'abs' )
  call graf3d ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0, -3.0, 3.0, -3.0, 1.0 )
  call height ( 50 )
  call title ( )

  call color ( 'green' )
  call surmat ( zmat, n, n, 1, 1 )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX10:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
