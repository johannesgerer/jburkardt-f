program main

!*****************************************************************************80
!
!! DISLIN_EX07B demonstrates demonstrates 3D bar graphs.
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

  integer, parameter :: n = 18

  character ( len = 80 ) :: cbuf
  integer :: i
  integer, dimension (n) :: icray = (/ &
     30,  30,  30,  30,  30,  30, 100, 100, 100, 100, &
    100, 100, 170, 170, 170, 170, 170, 170 /)
  real, dimension (n) :: xray  = (/ &
    1.0, 3.0, 8.0, 1.5, 9.0, 6.3, 5.8, 2.3, 8.1, 3.5, &
    2.2, 8.7, 9.2, 4.8, 3.4, 6.9, 7.5, 3.8 /)
  real, dimension (n) :: xwray
  real, dimension (n) :: yray  = (/ &
    5.0, 8.0, 3.5, 2.0, 7.0, 1.0, 4.3, 7.2, 6.0, 8.5, &
    4.1, 5.0, 7.3, 2.8, 1.6, 8.9, 9.5, 3.2 /)
  real, dimension (n) :: ywray
  real, dimension (n) :: z1ray = (/ &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  real, dimension (n) :: z2ray = (/ &
    4.0, 5.0, 3.0, 2.0, 3.5, 4.5, 2.0, 1.6, 3.8, 4.7, &
    2.1, 3.5, 1.9, 4.2, 4.9, 2.8, 3.6, 4.3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX07B:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of 3D bar graphs.'

  do i = 1, n
    xwray(i) = 0.5
    ywray(i) = 0.5
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
  call setfil ( 'dislin_ex07b.png' )
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
!  Use the standard HARDWARE font.
!
  call hwfont ( )
  call axspos ( 200, 2600 )
  call axslen ( 1800, 1800 )

  call name ( 'x-axis', 'x' )
  call name ( 'y-axis', 'y' )
  call name ( 'z-axis', 'z' )

  call titlin ( '3-d bars / bars3d', 3 )

  call labl3d ( 'hori' )
  call graf3d ( 0.0, 10.0, 0.0, 2.0, 0.0, 10.0, 0.0, 2.0, 0.0, 5.0, 0.0, 1.0 )
  call grid3d ( 1, 1, 'bottom' )

  call bars3d ( xray, yray, z1ray, z2ray, xwray, ywray, icray, n )

  call legini ( cbuf, 3, 20 )
  call legtit ( ' ' )
  call legpos ( 1300, 1100 )
  call leglin ( cbuf, 'first', 1 )
  call leglin ( cbuf, 'second', 2 )
  call leglin ( cbuf, 'third', 3 )
  call legend ( cbuf, 3 )

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
  write ( *, '(a)' ) 'DISLIN_EX07B:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
