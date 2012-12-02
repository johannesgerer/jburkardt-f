program main

!*****************************************************************************80
!
!! DISLIN_EX13 demonstrates a map plot.
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISLIN_EX13:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the creation of a map plot.'
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
  call setfil ( 'dislin_ex13.png' )
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

  call frame ( 3 )
  call axspos ( 400, 1850 )
  call axslen ( 2400, 1400 )

  call name ( 'longitude', 'x' )
  call name ( 'latitude', 'y' )
  call titlin ( 'World coastlines and lakes', 3 )

  call labels ( 'map', 'xy' )
  call grafmp ( -180.0, 180.0, -180.0, 90.0, -90.0, 90.0, -90.0, 30.0 ) 

  call gridmp ( 1, 1 )
  call color ( 'green' )
  call world ( )
  call color ( 'fore' )

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
  write ( *, '(a)' ) 'DISLIN_EX13:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
