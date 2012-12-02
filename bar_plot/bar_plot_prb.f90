program main

!*****************************************************************************80
!
!! MAIN is the main program for BAR_PLOT_PRB.
!
!  Discussion:
!
!    BAR_PLOT_PRB tests the BAR_PLOT routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAR_PLOT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BAR_PLOT library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAR_PLOT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BAR_DATA_EXAMPLE, BAR_PLOT_PPMA
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ncol = 800
  integer, parameter :: nrow = 200
  integer, parameter :: nx = 100
  integer, parameter :: ny = 6

  real ( kind = 8 ) b(nrow,ncol)
  character ( len = 80 ) file_name
  real ( kind = 8 ) g(nrow,ncol)
  integer ierror
  real ( kind = 8 ) r(nrow,ncol)
  real ( kind = 8 ) ytab(nx,ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  BAR_DATA_EXAMPLE sets up some sample data.'
  write ( *, '(a)' ) '  BAR_PLOT_PPMA plots it.'
!
!  Get some sample data.
!
  call bar_data_example ( nx, ny, ytab )
!
!  Convert the data to RGB arrays.
!
  call bar_data_to_rgb ( nx, ny, ytab, nrow, ncol, r, g, b )
!
!  Write the RGB arrays to a PPMA file.
!
  file_name = 'sines.ppma'

  call ppma_write ( file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'PPMA_WRITE returns IERROR = ', ierror
  end if

  return
end
