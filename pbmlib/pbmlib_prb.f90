program main

!*****************************************************************************80
!
!! MAIN is the main program for PBMLIB_PRB.
!
!  Discussion:
!
!    PBMLIB_PRB calls the PBMLIB test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_char = 100
  integer ( kind = 4 ), parameter :: max_col = 5
  integer ( kind = 4 ), parameter :: max_row = 7

  integer ( kind = 4 ) bits(max_row,max_col,max_char)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ipoint(0:255)
  integer ( kind = 4 ) nchar

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PBMLIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PBMLIB Library.'

  call test01 ( 'pbmlib.pbma' )
  call test02 ( 'pbmlib.pbma' )

  call test03 ( 'pbmlib.pgma' )
  call test04 ( 'pbmlib.pgma' )

  call test05 ( 'pbmlib.ppma' )
  call test06 ( 'pbmlib.ppma' )

  call test07 ( 'pbmlib.pbmb' )
  call test08 ( 'pbmlib.pbmb' )

  call test09 ( 'pbmlib.pgmb' )
  call test10 ( 'pbmlib.pgmb' )

  call test11 ( 'pbmlib.ppmb' )
  call test12 ( 'pbmlib.ppmb' )

  call test13 ( bits, ipoint, max_char, max_col, max_row, nchar, ierror )

  if ( ierror == 0 ) then
    call test14 ( bits, ipoint, max_char, max_col, max_row )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Skipping TEST14 because of error in TEST13!'
  end if

  call test15 ( bits, ipoint, max_char, max_col, max_row, nchar )

  call test16 ( )

  call test17 ( )
  call test18 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PBMLIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( file_name )

!*****************************************************************************80
!
!! TEST01 tests PBM_EXAMPLE, PBMA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 200
  integer ( kind = 4 ), parameter :: nrow = 200

  integer ( kind = 4 ) b(nrow,ncol)
  character ( len = * ) file_name
  integer ( kind = 4 ) ierror

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  PBMA_EXAMPLE sets some example PBM data.'
  write ( *, '(a)' ) '  PBMA_WRITE writes an ASCII PBM file.'

  call pbm_example ( nrow, ncol, b )

  call pbma_write ( file_name, ierror, nrow, ncol, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PBMA_WRITE returns IERROR = ', ierror
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Wrote the header and data for "' &
      // trim ( file_name ) //'".'
    write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
    write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  end if

  return
end
subroutine test02 ( file_name )

!*****************************************************************************80
!
!! TEST02 tests PBMA_READ_HEADER, PBMA_READ_DATA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension (:,:) :: b
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For an ASCII PBM file,'
  write ( *, '(a)' ) '  PBMA_READ_HEADER reads the header;'
  write ( *, '(a)' ) '  PBMA_READ_DATA reads the data.'

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    return
  end if
!
!  Read the header.
!
  call pbma_read_header ( file_unit, nrow, ncol, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not read the file header.'
    return
  end if
!
!  Allocate memory.
!
  allocate ( b(1:nrow,1:ncol) )
!
!  Read the data.
!
  call pbma_read_data ( file_unit, nrow, ncol, b, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not read the file data.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Read the header and data from "' &
      // trim ( file_name ) //'".'
    write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
    write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  end if

  deallocate ( b )

  return
end
subroutine test03 ( file_name )

!*****************************************************************************80
!
!! TEST03 tests PGM_EXAMPLE, PGMA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nrow = 200
  integer ( kind = 4 ), parameter :: ncol = 600

  character ( len = * ) file_name
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) ierror

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  PGMA_EXAMPLE sets up sample PGM data.'
  write ( *, '(a)' ) '  PGMA_WRITE writes an ASCII PGM file.'

  call pgm_example ( nrow, ncol, g )

  call pgma_write ( file_name, ierror, nrow, ncol, g )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PGMA_WRITE returns IERROR = ', ierror
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Wrote the header and data for "' &
      // trim ( file_name ) //'".'
    write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
    write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  end if

  return
end
subroutine test04 ( file_name )

!*****************************************************************************80
!
!! TEST04 tests PGMA_READ_HEADER, PGMA_READ_DATA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension (:,:) :: g
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) maxgray
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For an ASCII PGM file,'
  write ( *, '(a)' ) '  PGMA_READ_HEADER reads the header;'
  write ( *, '(a)' ) '  PGMA_READ_DATA reads the data.'

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    return
  end if
!
!  Read the header.
!
  call pgma_read_header ( file_unit, nrow, ncol, maxgray, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Fatal error!'
    write ( *, '(a)' ) '  Could not read the file header.'
    return
  end if
!
!  Allocate memory.
!
  allocate ( g(1:nrow,1:ncol) )
!
!  Read the data.
!
  call pgma_read_data ( file_unit, nrow, ncol, g, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Fatal error!'
    write ( *, '(a)' ) '  Could not read the file data.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header and data from "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

  deallocate ( g )

  return
end
subroutine test05 ( file_name )

!*****************************************************************************80
!
!! TEST05 tests PPM_EXAMPLE, PPMA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 300
  integer ( kind = 4 ), parameter :: nrow = 300

  integer ( kind = 4 ) b(nrow,ncol)
  character ( len = * ) file_name
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) r(nrow,ncol)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  PPM_EXAMPLE sets up sample PPM data.'
  write ( *, '(a)' ) '  PPMA_WRITE writes an ASCII PPM file.'

  call ppm_example ( nrow, ncol, r, g, b )

  call ppma_write ( file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PPMA_WRITE returns IERROR = ', ierror
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Wrote the header and data for "' &
      // trim ( file_name ) //'".'
    write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
    write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  end if

  return
end
subroutine test06 ( file_name )

!*****************************************************************************80
!
!! TEST06 tests PPMA_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_p = 300*300

  integer ( kind = 4 ) b(max_p)
  character ( len = * ) file_name
  integer ( kind = 4 ) g(max_p)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) r(max_p)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  PPMA_READ reads an ASCII PPM file.'

  call ppma_read ( file_name, ierror, maxcol, max_p, nrow, ncol, r, g, b  )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06'
    write ( *, '(a,i8)' ) '  PPMA_READ returns IERROR = ', ierror
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header and data from "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

  return
end
subroutine test07 ( file_name )

!*****************************************************************************80
!
!! TEST07 tests PBM_EXAMPLE, PBMB_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 200
  integer ( kind = 4 ), parameter :: nrow = 200

  integer ( kind = 4 ) b(nrow,ncol)
  character ( len = * ) file_name
  integer ( kind = 4 ) ierror

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  PBM_EXAMPLE sets some sample PBM data.'
  write ( *, '(a)' ) '  PBMB_WRITE writes a binary PBM file.'

  call pbm_example ( nrow, ncol, b )
!
!  Reverse black and white.
!
  b(1:nrow,1:ncol) = 1 - b(1:nrow,1:ncol)

  call pbmb_write ( file_name, ierror, nrow, ncol, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PBMB_WRITE returns IERROR = ', ierror
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Wrote the header and data for "' &
      // trim ( file_name ) //'".'
    write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
    write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  end if

  return
end
subroutine test08 ( file_name )

!*****************************************************************************80
!
!! TEST08 tests PBMB_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_g = 200*200

  integer ( kind = 4 ) b(max_g)
  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  PBMB_READ reads a binary PBM file.'

  call pbmb_read ( file_name, ierror, max_g, nrow, ncol, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08'
    write ( *, '(a,i8)' ) '  PBMB_READ returns IERROR = ', ierror
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header and data from "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

  return
end
subroutine test09 ( file_name )

!*****************************************************************************80
!
!! TEST09 tests PGM_EXAMPLE, PGMB_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 200
  integer ( kind = 4 ), parameter :: nrow = 200

  character ( len = * ) file_name
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) ierror

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  PGMB_WRITE writes a binary PGM file.'

  call pgm_example ( nrow, ncol, g )

  call pgmb_write ( file_name, ierror, nrow, ncol, g )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PGMB_WRITE returns IERROR = ', ierror
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Wrote the header and data for "' &
      // trim ( file_name ) //'".'
    write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
    write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  end if

  return
end
subroutine test10 ( file_name )

!*****************************************************************************80
!
!! TEST10 tests PGMB_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_g = 200*200

  character ( len = * ) file_name
  integer ( kind = 4 ) g(max_g)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  PGMB_READ reads a binary PGM file.'

  call pgmb_read ( file_name, ierror, maxcol, max_g, nrow, ncol, g )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PGMB_READ returns IERROR = ', ierror
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header and data from "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

  return
end
subroutine test11 ( file_name )

!*****************************************************************************80
!
!! TEST11 tests PPM_EXAMPLE,s PPMB_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 300
  integer ( kind = 4 ), parameter :: nrow = 300

  integer ( kind = 4 ) b(nrow,ncol)
  character ( len = * ) file_name
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) r(nrow,ncol)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  PPM_EXAMPLE sets up sample PPM data.'
  write ( *, '(a)' ) '  PPMB_WRITE writes a binary PPM file.'

  call ppm_example ( nrow, ncol, r, g, b )

  call ppmb_write ( file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PPMB_WRITE returns IERROR = ', ierror
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Wrote the header and data for "' &
      // trim ( file_name ) //'".'
    write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
    write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  end if

  return
end
subroutine test12 ( file_name )

!*****************************************************************************80
!
!! TEST12 tests PPMB_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_g = 300*300

  integer ( kind = 4 ) b(max_g)
  character ( len = * ) file_name
  integer ( kind = 4 ) g(max_g)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) r(max_g)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  PPMB_READ reads a binary PPM file.'

  call ppmb_read ( file_name, ierror, maxcol, max_g, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST12'
    write ( *, '(a,i8)' ) 'PPMB_READ returns IERROR = ', ierror
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header and data from "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

  return
end
subroutine test13 ( bits, ipoint, maxchar, maxcol, maxrow, nchar, ierror )

!*****************************************************************************80
!
!! TEST13 tests FONT_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxchar
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) maxrow

  integer ( kind = 4 ) bits(maxrow,maxcol,maxchar)
  character ( len = 80 ) :: file_name = 'alphabits.txt'
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) inunit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipoint(0:255)
  integer ( kind = 4 ) nchar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  FONT_READ reads in a simple bit map font.'

  ierror = 0
  inunit = 1

  open ( unit = inunit, file = file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST13'
    write ( *, '(a)' ) '  Could not open the font file!'
    return
  end if

  call font_read ( bits, ierror, inunit, ipoint, maxchar, maxcol, maxrow, &
    nchar )

  close ( unit = inunit )

  return
end
subroutine test14 ( bits, ipoint, maxchar, maxcol, maxrow )

!*****************************************************************************80
!
!! TEST14 tests FONT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxchar
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) maxrow

  integer ( kind = 4 ) bits(maxrow,maxcol,maxchar)
  integer ( kind = 4 ) ipoint(0:255)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  FONT_PRINT prints the font we just read in.'

  call font_print ( bits, ipoint, maxchar, maxcol, maxrow )

  return
end
subroutine test15 ( bits, ipoint, maxchar, maxcol, maxrow, nchar )

!*****************************************************************************80
!
!! TEST15 tests FONT_DATA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxchar
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) maxrow

  integer ( kind = 4 ) bits(maxrow,maxcol,maxchar)
  integer ( kind = 4 ) ipoint(0:255)
  integer ( kind = 4 ) nchar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  FONT_DATA prints out a font data statement.'

  call font_data ( bits, ipoint, maxchar, maxcol, maxrow, nchar )

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests BITCHR75.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character, parameter, dimension ( 0:1 ) :: bit = (/ ' ', '*' /)
  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) pattern(7,5)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  BITCHR75 returns a 7 x 5 representation of'
  write ( *, '(a)' ) '    SOME characters.'

  do k = 32, 96

    c = char ( k )

    call bitchr75 ( c, pattern )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) 'Character ', k, ' is "' // c // '"'
    write ( *, '(a)' ) ' '

    do i = 1, 7
      write ( *, '(5a1)' ) ( bit ( pattern(i,j) ), j = 1, 5 )
    end do

  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests PPMB_WRITE.
!
!  Discussion:
!
!    This example uses color to record which starting values converge to
!    which roots in a Newton iteration in the complex plane.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 1152
  integer ( kind = 4 ), parameter :: nrow = 864

  integer ( kind = 4 ) b(nrow,ncol)
  character ( len = 80 ) :: file_name = 'newton.ppmb'
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 20
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r(nrow,ncol)
  real ( kind = 8 ) radius
  complex ( kind = 8 ) root1
  complex ( kind = 8 ) root2
  complex ( kind = 8 ) root3
  real ( kind = 8 ) theta
  real ( kind = 8 ), parameter :: tol = 0.01D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xmax = 1.4D+00
  real ( kind = 8 ), parameter :: xmin = -1.4D+00
  real ( kind = 8 ) y
  real ( kind = 8 ), parameter :: ymax = 1.2D+00
  real ( kind = 8 ), parameter :: ymin = -1.2D+00
  complex ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  PPMB_WRITE writes a PPMB file.'
  write ( *, '(a)' ) '  This displays the basins of attraction for'
  write ( *, '(a)' ) '  Newton''s method applied to a particular'
  write ( *, '(a)' ) '  nonlinear equation in the complex plane.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The equation is Z^3 = -1, which has three roots.'
  write ( *, '(a)' ) ' '

  radius = 1.0D+00
  theta = pi / 3.0D+00
  call rt_to_xy ( radius, theta, x, y )
  root1 = cmplx ( x, y, kind = 8 )
  write ( *, '(a,2g14.6)' ) '  ROOT1 = ', root1

  theta = pi
  call rt_to_xy ( radius, theta, x, y )
  root2 = cmplx ( x, y, kind = 8 )
  write ( *, '(a,2g14.6)' ) '  ROOT2 = ', root2

  theta = 5.0D+00 * pi / 3.0D+00
  call rt_to_xy ( radius, theta, x, y )
  root3 = cmplx ( x, y, kind = 8 )
  write ( *, '(a,2g14.6)' ) '  ROOT3 = ', root3

  do i = 1, nrow

    y = ( real ( nrow - i,     kind = 8 ) * ymax   &
        + real (        i - 1, kind = 8 ) * ymin ) &
        / real ( nrow     - 1, kind = 8 )

    do j = 1, ncol

      x = ( real ( ncol - j,     kind = 8 ) * xmin   &
          + real (        j - 1, kind = 8 ) * xmax ) &
          / real ( ncol     - 1, kind = 8 )

      z = cmplx ( x, y, kind = 8 )

      r(i,j) = 100
      g(i,j) = 100
      b(i,j) = 100

      do it = 1, it_max
!
!  This is one step of Newton's iteration for f(z) = z^3 + 1.
!
        z = z - ( z * z * z + 1.0D+00 ) / ( 3.0D+00 * z * z )
!
!  This is one step of Sterling's iteration for f(z) = z^3 + 1.
!  We have commented this one out.
!
!       z = z - ( z * z * z + 1.0D+00 ) / ( 3.0D+00 * ( z - z**3 - 1.0D+00 )**2 )
!
!  It the current iterate is close enough to a root, exit.
!  But mark the starting point with a shade of red, green or blue
!  depending on which root it converged to.
!
!  To make visualization easier, we gray out the red, green or blue
!  depending on the number of iterations taken.
!
!  Also, every third iteration we set the color to white, to make the
!  bands easier to see.
!
        if ( abs ( z - root1 ) < tol ) then

          if ( mod ( it, 3 ) /= 2 ) then

            r(i,j) = ( ( it_max - it + 1 ) * 255   &
                     + (          it - 1 ) * 127 ) &
                     /   it_max

            g(i,j) = ( ( it_max - it + 1 ) *   0   &
                     + (          it - 1 ) * 127 ) &
                     /   it_max

            b(i,j) = ( ( it_max - it + 1 ) *   0   &
                     + (          it - 1 ) * 127 ) &
                     /   it_max

          else

            r(i,j) = 255
            g(i,j) = 255
            b(i,j) = 255
          end if

          exit

        else if ( abs ( z - root2 ) < tol ) then

          if ( mod ( it, 3 ) /= 2 ) then
            r(i,j) = ( ( it_max + 1 - it ) *   0 + ( it - 1 ) * 127 ) / it_max
            g(i,j) = ( ( it_max + 1 - it ) * 255 + ( it - 1 ) * 127 ) / it_max
            b(i,j) = ( ( it_max + 1 - it ) *   0 + ( it - 1 ) * 127 ) / it_max
          else
            r(i,j) = 255
            g(i,j) = 255
            b(i,j) = 255
          end if

          exit

        else if ( abs ( z - root3 ) < tol ) then

          if ( mod ( it, 3 ) /= 2 ) then
            r(i,j) = ( ( it_max + 1 - it ) *   0 + ( it - 1 ) * 127 ) / it_max
            g(i,j) = ( ( it_max + 1 - it ) *   0 + ( it - 1 ) * 127 ) / it_max
            b(i,j) = ( ( it_max + 1 - it ) * 255 + ( it - 1 ) * 127 ) / it_max
          else
            r(i,j) = 255
            g(i,j) = 255
            b(i,j) = 255
          end if

          exit

        end if

      end do

    end do

  end do
!
!  Now fill with gray the three roots.
!
  do i = 1, nrow

    y = ( real ( nrow - i,     kind = 8 ) * ymax   &
        + real (        i - 1, kind = 8 ) * ymin ) &
        / real ( nrow     - 1, kind = 8 )

    do j = 1, ncol

      x = ( real ( ncol - j,     kind = 8 ) * xmin   &
          + real (        j - 1, kind = 8 ) * xmax ) &
          / real ( ncol     - 1, kind = 8 )

      z = cmplx ( x, y, kind = 8 )

      if ( abs ( z - root1 ) < 4.0D+00 * tol ) then

        r(i,j) = 100
        g(i,j) = 100
        b(i,j) = 100

      else if ( abs ( z - root2 ) < 4.0D+00 * tol ) then

        r(i,j) = 100
        g(i,j) = 100
        b(i,j) = 100

      else if ( abs ( z - root3 ) < 4.0D+00 * tol ) then

        r(i,j) = 150
        g(i,j) = 150
        b(i,j) = 150

      end if

    end do

  end do
!
!  Now write all this data to a file.
!
  call ppmb_write ( file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PPMB_WRITE returns IERROR = ', ierror
  end if

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests PPMA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 300
  integer ( kind = 4 ), parameter :: nrow = 300

  real ( kind = 8 ) angle
  integer ( kind = 4 ) b(nrow,ncol)
  character ( len = 80 ) :: file_name = 'hexcol.ppma'
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r(nrow,ncol)
  real ( kind = 8 ) rb
  real ( kind = 8 ) rg
  real ( kind = 8 ) rr
  real ( kind = 8 ) scale
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  PPMA_WRITE writes an ASCII PPM file.'

  do i = 1, nrow

    y = real ( 2 * i - nrow, kind = 8 )

    do j = 1, ncol

      x = real ( 2 * j - ncol, kind = 8 )

      angle = 180.0D+00 * ( atan2 ( y, x ) / pi )

      call hexcol ( angle, rr, rg, rb )

      r(i,j) = int ( 255.0D+00 * rr )
      g(i,j) = int ( 255.0D+00 * rg )
      b(i,j) = int ( 255.0D+00 * rb )

    end do
  end do

  call ppma_write ( file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PPMA_WRITE returns IERROR = ', ierror
  end if

  return
end
