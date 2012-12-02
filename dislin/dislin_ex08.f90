program main

!*****************************************************************************80
!
!! DISLIN_EX08 demonstrates various shading patterns.
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

  character ( len = 2 ) cstr
  character ( len = 60 ) ctit
  integer i
  integer iclr
  integer ii
  integer, dimension ( 4 ) :: ix = (/ 0, 300, 300, 0 /)
  integer ixp(4)
  integer, dimension ( 4 ) :: iy = (/ 0, 0, 400, 400 /)
  integer iyp(4)
  integer j
  integer k
  integer nl
  integer nx
  integer nx0
  integer ny
  integer ny0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX08:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the use of shading patterns.'

  ctit = 'Shading Patterns (AREAF)'
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
  call setfil ( 'dislin_ex08.png' )
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
!  Use the COMPLEX font.
!
  call complx ( )
  call setvlt ( 'SMALL' )

  call height ( 50 )
  nl = nlmess ( ctit )
  nx = ( 2970 - nl ) / 2
  call messag ( ctit, nx, 200 )

  nx0 = 335
  ny0 = 350

  do i = 1, 3

    ny = ny0 + ( i - 1 ) * 600

    do j = 1, 6

      iclr = ( i - 1 ) * 6 + j - 1
      iclr = mod ( iclr, 15 )

      if ( iclr == 0 ) then
        iclr = 15
      end if

      call setclr ( iclr )

      nx = nx0 + ( j - 1 ) * 400
      ii = ( i - 1 ) * 6 + j - 1
      call shdpat ( ii )
      write ( cstr, '(i2)' ) ii

      do k = 1, 4
        ixp(k) = ix(k) + nx
        iyp(k) = iy(k) + ny
      end do

      call areaf ( ixp, iyp, 4 )

      nl = nlmess ( cstr )
      nx = nx + ( 300 - nl ) / 2
      call messag ( cstr, nx, ny+460 )

    end do

  end do
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX08:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
