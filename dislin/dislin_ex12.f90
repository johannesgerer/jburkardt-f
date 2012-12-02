program main

!*****************************************************************************80
!
!! DISLIN_EX12 demonstrates a shaded contour plot.
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

  integer, parameter :: m = 50
  integer, parameter :: n = 50

  integer :: i
  integer :: j
  real :: x
  real, dimension (m) :: xray
  real :: y
  real, dimension (n) :: yray
  real, dimension (12) :: zlev
  real, dimension (m,n) :: zmat

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX12:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of '
  write ( *, '(a)' ) '  a shaded contour plot.'

  do i = 1, m
    x = 1.6 * real ( i - 1 ) / real ( m - 1 )
    xray(i) = x
    do j = 1, n
      y = 1.6 * real ( j - 1 ) / real ( n - 1 )
      yray(j) = y
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
  call setfil ( 'dislin_ex12.png' )
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

  call mixalf ( )
  call titlin ( 'Shaded contour plot', 1 )
  call titlin ( 'f(x,y) = (x[2$ - 1)[2$ + (y[2$ - 1)[2$', 3 )
  call name ( 'x-axis', 'x' )
  call name ( 'y-axis', 'y' )

  call shdmod ( 'poly', 'contur' )
  call axspos ( 450, 2670 )
  call graf ( 0.0, 1.6, 0.0, 0.2, 0.0, 1.6, 0.0, 0.2 )

  do i = 1, 12
    zlev(13-i) = 0.1 + real( i - 1 ) * 0.1
  end do

  call conshd ( xray, m, yray, n, zmat, zlev, 12 )

  call height ( 50 )
  call title ( )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX12:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end

