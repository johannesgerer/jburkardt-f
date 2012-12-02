program main

!*****************************************************************************80
!
!! DISLIN_EX03 demonstrates the creation and display of special symbols.
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

  character ( len = 2 ) :: cstr
  character ( len = 20 ) :: ctit = 'Symbols'
  integer :: i
  integer :: nl
  integer :: nxp
  integer :: ny

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX03:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the use of the SYMBOL routine'
  write ( *, '(a)' ) '  to create and display special symbols.'
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
  call setfil ( 'dislin_ex03.png' )
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

  call height ( 60 )
  nl = nlmess ( ctit )
  call messag ( ctit, (2100-nl)/2, 200 )

  call height ( 50 )
  call hsymbl ( 120 )

  ny = 150

  do i = 0, 21

    if ( mod ( i, 4 ) == 0 ) then
      ny = ny + 400
      nxp = 550
    else
      nxp = nxp + 350
    end if

    if ( i < 10 ) then
      write ( cstr, '(i1)' ) i
    else
      write ( cstr, '(i2)' ) i
    end if

    nl = nlmess ( cstr ) / 2
    call messag ( cstr, nxp-nl, ny+150 )
    call symbol ( i, nxp, ny )

  end do
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX03:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
