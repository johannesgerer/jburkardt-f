program main

!*****************************************************************************80
!
!! DISLIN_EX05 demonstrates the creation of bar graphs.
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

  character ( len = 24 ) :: cbuf
  character ( len = 60 ) :: ctit = 'Bar graphs (BARS)'
  integer i
  integer, parameter :: nya = 2700
  real, dimension (9) :: x  = (/ &
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 /)
  real, dimension ( 9 ) :: y = (/ &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  real, dimension ( 9 ) :: y1 = (/ &
    1.0, 1.5, 2.5, 1.3, 2.0, 1.2, 0.7, 1.4, 1.1 /)
  real, dimension ( 9 ) :: y2 = (/ &
    2.0, 2.7, 3.5, 2.1, 3.2, 1.9, 2.0, 2.3, 1.8 /)
  real, dimension ( 9 ) :: y3 = (/ &
    4.0, 3.5, 4.5, 3.7, 4.0, 2.9, 3.0, 3.2, 2.6 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX05:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of bar graphs.'
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
  call setfil ( 'dislin_ex05.png' )
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
  call ticks ( 1, 'x' )
  call intax ( )
  call axslen ( 1600, 700 )
  call titlin ( ctit, 3 )

  call legini ( cbuf, 3, 8 )
  call leglin ( cbuf, 'First', 1 )
  call leglin ( cbuf, 'Second', 2 )
  call leglin ( cbuf, 'Third', 3 )
  call legtit ( ' ')
  call shdpat ( 5 )

  do i = 1, 3

    if ( 1 < i ) then
      call labels ( 'none', 'x' )
    end if

    call axspos ( 300, nya-(i-1)*800 )

    call graf ( 0.0, 10.0, 0.0, 1.0, 0.0, 5.0, 0.0, 1.0 )

    if ( i == 1 ) then
      call bargrp ( 3, 0.15 )
      call color ( 'red' )
      call bars ( x, y, y1, 9 )
      call color ( 'green' )
      call bars ( x, y, y2, 9 )
      call color ( 'blue' )
      call bars ( x, y, y3, 9 )
      call color ( 'fore' )
      call reset ( 'bargrp' )
    else if ( i == 2 ) then
      call height ( 30 )
      call labels ( 'delta', 'bars' )
      call labpos ( 'center', 'bars' )
      call color ( 'red' )
      call bars ( x, y, y1, 9 )
      call color ( 'green' )
      call bars ( x, y1, y2, 9 )
      call color ( 'blue' )
      call bars ( x, y2, y3, 9 )
      call color ( 'fore' )
      call reset ( 'height' )
    else if ( i == 3 ) then
      call labels ( 'Second', 'bars' )
      call labpos ( 'outside', 'bars' )
      call color ( 'red' )
      call bars ( x, y, y1, 9 )
      call color ( 'fore' )
    end if

    if ( i < 3 ) then
      call legend ( cbuf, 7 )
    else if ( i == 3 ) then
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
  write ( *, '(a)' ) 'DISLIN_EX05:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
