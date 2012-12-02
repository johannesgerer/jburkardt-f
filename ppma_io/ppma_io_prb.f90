program main

!*****************************************************************************80
!
!! MAIN is the main program for PPMA_IO_PRB.
!
!  Discussion:
!
!    PPMA_IO_PRB calls the PPMA_IO test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PPMA_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PPMA_IO library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PPMA_IO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests PPM_EXAMPLE and PPMA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: ncol = 300
  integer   ( kind = 4 ), parameter :: nrow = 300

  integer   ( kind = 4 ) b(nrow,ncol)
  character ( len = 80 ) :: file_name = 'test01.ascii.ppm'
  integer   ( kind = 4 ) g(nrow,ncol)
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) r(nrow,ncol)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  PPMA_EXAMPLE sets up sample PPMA data.'
  write ( *, '(a)' ) '  PPMA_WRITE writes an ASCII PPMA file.'

  call ppma_example ( nrow, ncol, r, g, b )

  call ppma_write ( file_name, nrow, ncol, r, g, b, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a,i6)' ) 'PPMA_WRITE returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the header and data for "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i6)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i6)' ) '  Number of columns of data = ', ncol

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests PPMA_READ_DATA and PPMA_READ_HEADER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: b
  character ( len = 80 ) :: file_name = 'test02.ascii.ppm'
  integer   ( kind = 4 ) file_unit
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: g
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) k
  integer   ( kind = 4 ) ncol
  integer   ( kind = 4 ) nrow
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: r
  integer   ( kind = 4 ) rgb_max

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  PPMA_READ_HEADER reads the header.'
  write ( *, '(a)' ) '  PPMA_READ_HEADER reads the data.'

  call ppma_write_test ( file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PPMA_WRITE_TEST created some data.'

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    return
  end if

  call ppma_read_header ( file_unit, nrow, ncol, rgb_max, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Error while reading the header.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PPMA_READ_HEADER read the header.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i6)' ) '  Number of columns of data = ', ncol
  write ( *, '(a,i6)' ) '  Maximum RGB value =         ', rgb_max

  allocate ( r(nrow,ncol) )
  allocate ( g(nrow,ncol) )
  allocate ( b(nrow,ncol) )

  call ppma_read_data ( file_unit, nrow, ncol, r, g, b, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  PPMA_READ_DATA failed.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PPMA_READ_DATA read the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Ten sample values:'
  write ( *, '(a)' ) ' '
  do k = 1, 10
    i = ( ( 10 - k ) * 1 + ( k - 1 ) * nrow ) / ( 10 - 1 )
    j = ( ( 10 - k ) * 1 + ( k - 1 ) * ncol ) / ( 10 - 1 )
    write ( *, '(2i4,4x,3i6)' ) i, j, r(i,j), g(i,j), b(i,j)
  end do

  call ppma_check_data ( nrow, ncol, rgb_max, r, g, b, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Error!'
    write ( *, '(a,i6)' ) '  The data was not accepted by PPMA_CHECK_DATA.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was accepted by PPMA_CHECK_DATA.'

  deallocate ( r )
  deallocate ( g )
  deallocate ( b )

  return
end
