program main

!*****************************************************************************80
!
!! DISLIN_EX11 demonstrates a contour plot.
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

  real :: fpi
  integer :: i
  integer :: j
  real :: step
  real, dimension (n) :: xray
  real, dimension (n) :: yray
  real :: zlev
  real, dimension (n,n) :: zmat

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX11:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of '
  write ( *, '(a)' ) '  a contour plot.'

  fpi = 3.14159 / 180.0
  step = 360.0 / real ( n - 1 )

  do i = 1, n
    xray(i) = real ( i - 1 ) * step
    yray(i) = real ( i - 1 ) * step
  end do

  do i = 1, n
    do j = 1, n
      zmat(i,j) = 2.0 * sin ( xray(i) * fpi ) * sin ( yray(j) * fpi )
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
  call setfil ( 'dislin_ex11.png' )
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

  call titlin ( 'Contour plot', 1 )
  call titlin ( 'f(x,y) = 2 * sin(x) * sin(y)', 3 )

  call name ( 'x-axis', 'x' )
  call name ( 'y-axis', 'y' )

  call intax ( )
  call axspos ( 450, 2670 )
  call graf ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0 )

  call height ( 30 )

  do i = 1, 9

    zlev = -2.0 + real ( i - 1 ) * 0.5
    call setclr ( i * 25 )

    if ( i == 5 ) then
      call labels ( 'none', 'contur' )
    else
      call labels ( 'float', 'contur' )
    end if

    call contur ( xray, n, yray, n, zmat, zlev )

  end do

  call height ( 50 )
  call color ( 'fore' )
  call title ( )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX11:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
