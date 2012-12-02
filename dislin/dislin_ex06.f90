program main

!*****************************************************************************80
!
!! DISLIN_EX06 demonstrates the creation of pie charts.
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

  character ( len = 40 ) :: cbuf
  character ( len = 60 ) :: ctit = 'Pie charts (PIEGRF)'
  integer :: i
  integer :: nya = 2800
  real, dimension (5) :: xray = (/ 1.0, 2.5, 2.0, 2.7, 1.8 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX06:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of pie charts.'
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
  call setfil ( 'dislin_ex06.png' )
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
  call axslen ( 1600, 1000 )
  call titlin ( ctit, 2 )
  call chnpie ( 'both' )

  call legini ( cbuf, 5, 8 )
  call leglin ( cbuf, 'First', 1 )
  call leglin ( cbuf, 'Second', 2 )
  call leglin ( cbuf, 'Third', 3 )
  call leglin ( cbuf, 'Fourth', 4 )
  call leglin ( cbuf, 'Fifth', 5 )

  call patcyc ( 1, 7 )
  call patcyc ( 2, 4 )
  call patcyc ( 3, 13 )
  call patcyc ( 4, 3 )
  call patcyc ( 5, 5 )

  do i = 1, 2

    call axspos ( 250, nya-(i-1)*1200 )

    if ( i == 2 ) then
      call labels ( 'data', 'pie' )
      call labpos ( 'external', 'pie' )
    end if

    call piegrf ( cbuf, 1, xray, 5 )

    if ( i == 2 ) then
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
  write ( *, '(a)' ) 'DISLIN_EX06:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
