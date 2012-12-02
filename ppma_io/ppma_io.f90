subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine getint ( done, ierror, inunit, ival, string )

!*****************************************************************************80
!
!! GETINT reads an integer from a file.
!
!  Discussion:
!
!    The file, or at least the part read by GETINT, is assumed to
!    contain nothing but integers.  These integers may be separated
!    by spaces, or appear on separate lines.  Comments, which begin
!    with "#" and extend to the end of the line, may appear anywhere.
!
!    Each time GETINT is called, it tries to read the next integer
!    it can find.  It remembers where it was in the current line
!    of text.
!
!    The user should open a text file on FORTRAN unit INUNIT,
!    set STRING = ' ' and DONE = TRUE.  The GETINT routine will take
!    care of reading in a new STRING as necessary, and extracting
!    as many integers as possible from the line of text before
!    reading in the next line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, logical DONE.
!    On input, if this is the first call, or the user has changed
!    STRING, then set DONE = TRUE.
!    On output, if there is no more data to be read from STRING,
!    then DONE is TRUE.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred while trying to read the integer.
!
!    Input, integer ( kind = 4 ) INUNIT, the FORTRAN unit from which to read.
!
!    Output, integer ( kind = 4 ) IVAL, the integer that was read.
!
!    Input/output, character ( len = * ) STRING, the text of the most recently
!    read line of the file.
!
  implicit none

  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) inunit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) string
  character ( len = 80 ) word

  do

    call word_next_rd ( string, word, done )

    if ( .not. done ) then
      exit
    end if

    read ( inunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

    i = index ( string, '#' )
    if ( i /= 0 ) then
      string(i:) = ' '
    end if

  end do

  call s_to_i4 ( word, ival, ierror, last )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GETINT - Fatal error!'
    write ( *, '(a)' ) '  Error trying to convert string to integer.'
    return
  end if

  return
end
subroutine ppma_check_data ( row_num, col_num, rgb_max, r, g, b, ierror )

!*****************************************************************************80
!
!! PPMA_CHECK_DATA checks pixel data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COL_NUM, ROW_NUM, the number of rows
!    and columns of data.
!
!    Input, integer ( kind = 4 ) RGB_MAX, the maximum RGB value.
!
!    Input, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
!    B(ROW_NUM,COL_NUM), contains the RGB pixel data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, the data is illegal.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) r(row_num,col_num)
  integer ( kind = 4 ) rgb_max

  ierror = 0
!
!  Make sure no color is negative.
!
  if ( minval ( r(1:row_num,1:col_num) ) < 0 .or. &
       minval ( g(1:row_num,1:col_num) ) < 0 .or. &
       minval ( b(1:row_num,1:col_num) ) < 0 ) then
    ierror = 1
    return
  end if
!
!  Make sure no color is greater than RGB_MAX.
!
  if ( rgb_max < maxval ( r(1:row_num,1:col_num) ) .or. &
       rgb_max < maxval ( g(1:row_num,1:col_num) ) .or. &
       rgb_max < maxval ( b(1:row_num,1:col_num) ) ) then
    ierror = 1
    return
  end if

  return
end
subroutine ppma_example ( row_num, col_num, r, g, b )

!*****************************************************************************80
!
!! PPMA_EXAMPLE sets up sample RGB data suitable for a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and
!    columns of data.  A reasonable value is 200 for ROW_NUM and 600 for COL_NUM.
!
!    Output, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
!    B(ROW_NUM,COL_NUM), the RGB data.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) r(row_num,col_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do i = 1, row_num

    y = real ( row_num - i, kind = 8 ) / real ( row_num - 1, kind = 8 )

    do j = 1, col_num

      x = real ( j - 1, kind = 8 ) / real ( col_num - 1, kind = 8 )

      f1 = 4.0D+00 * ( x - 0.5D+00 )**2
      f2 = sin ( 3.14159265D+00 * x )
      f3 = x

      if ( y <= f1 ) then
        r(i,j) = int ( 255.0D+00 * f1 )
      else
        r(i,j) = 50
      end if

      if ( y <= f2 ) then
        g(i,j) = int ( 255.0D+00 * f2 )
      else
        g(i,j) = 150
      end if

      if ( y <= f3 ) then
        b(i,j) = int ( 255.0D+00 * f3 )
      else
        b(i,j) = 250
      end if

    end do
  end do

  return
end
subroutine ppma_read_data ( file_in_unit, row_num, col_num, r, g, b, ierror )

!*****************************************************************************80
!
!! PPMA_READ_DATA reads the data in a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_IN_UNIT, the unit number of the file.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and
!    columns of data.
!
!    Output, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
!    B(ROW_NUM,COL_NUM),
!    the RGB data.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  logical done
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) r(row_num,col_num)
  character ( len = 80 ) string

  ierror = 0
  done = .true.
  string = ' '

  do i = 1, row_num
    do j = 1, col_num

      call getint ( done, ierror, file_in_unit, r(i,j), string )

      if ( ierror /= 0 ) then
        ierror = 5
        close ( unit = file_in_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_READ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Problem reading R data.'
        return
      end if

      call getint ( done, ierror, file_in_unit, g(i,j), string )

      if ( ierror /= 0 ) then
        ierror = 5
        close ( unit = file_in_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_READ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Problem reading G data.'
        return
      end if

      call getint ( done, ierror, file_in_unit, b(i,j), string )

      if ( ierror /= 0 ) then
        ierror = 5
        close ( unit = file_in_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_READ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Problem reading B data.'
        return
      end if

    end do
  end do

  return
end
subroutine ppma_read_header ( file_in_unit, row_num, col_num, rgb_max, ierror )

!*****************************************************************************80
!
!! PPMA_READ_HEADER reads the header of a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_IN_UNIT, the unit number of the file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows
!    and columns of data.
!
!    Output, integer ( kind = 4 ) RGB_MAX, the maximum RGB value.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  logical done
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 2 ) magic
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ) rgb_max
  logical s_eqi
  character ( len = 80 ) string
!
!  Read the first line of data, which must begin with the magic number.
!
  read ( file_in_unit, '(a)', iostat = ios ) magic

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  End or error while reading file.'
    ierror = 2
    return
  end if

  if ( .not. s_eqi ( magic, 'P3' ) ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_HEADER - Fatal error.'
    write ( *, '(a)' ) '  First two bytes are not magic number "P3".'
    write ( *, '(a)' ) '  First two bytes are: "' // magic // '".'
    return
  end if
!
!  Now search for COL_NUM, ROW_NUM, and RGB_MAX.
!
  done = .true.
  string = ' '

  call getint ( done, ierror, file_in_unit, col_num, string )

  if ( ierror /= 0 ) then
    close ( unit = file_in_unit )
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading COL_NUM.'
    return
  end if

  call getint ( done, ierror, file_in_unit, row_num, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_in_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading ROW_NUM.'
    return
  end if

  call getint ( done, ierror, file_in_unit, rgb_max, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_in_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading RGB_MAX.'
    return
  end if

  return
end
subroutine ppma_read_test ( file_in_name, ierror )

!*****************************************************************************
!
!! PPMA_READ_TEST tests the ASCII portable pixel map read routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the file
!    containing the ASCII portable pixel map data.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag which is nonzero if
!    there was an error.
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: b
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: g
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: r
  integer ( kind = 4 ) rgb_max

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_TEST - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    return
  end if
!
!  Read the header.
!
  call ppma_read_header ( file_in_unit, row_num, col_num, rgb_max, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_TEST - Fatal error!'
    write ( *, '(a)' ) '  PPMA_READ_HEADER failed.'
    return
  end if
!
!  Allocate the data.
!
  allocate ( r(row_num,col_num) )
  allocate ( g(row_num,col_num) )
  allocate ( b(row_num,col_num) )
!
!  Read the data.
!
  call ppma_read_data ( file_in_unit, row_num, col_num, r, g, b, ierror )

  if ( ierror /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_TEST - Fatal error!'
    write ( *, '(a)' ) '  PPMA_READ_HEADER failed.'

    deallocate ( r )
    deallocate ( g )
    deallocate ( b )

    return
  end if

  close ( unit = file_in_unit )
!
!  Check the data.
!
  call ppma_check_data ( row_num, col_num, rgb_max, r, g, b, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_TEST - Fatal error!'
    write ( *, '(a)' ) '  PPMA_CHECK_DATA did not approve the data.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ_TEST:'
    write ( *, '(a)' ) '  PPMA_CHECK_DATA has approved the data from the file.'
  end if

  deallocate ( r )
  deallocate ( g )
  deallocate ( b )

  return
end
subroutine ppma_write ( file_out_name, row_num, col_num, r, g, b, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE writes an ASCII portable pixel map file.
!
!  Discussion:
!
!    PPM files can be viewed by XV.
!
!    Programs to convert files to this format include:
!
!      GIFTOPPM  - GIF file
!      PGMTOPPM  - Portable Gray Map file
!      PICTTOPPM - Macintosh PICT file
!      XPMTOPPM  - X11 pixmap file
!
!    Various programs can convert other formats to PPM format, including:
!
!      BMPTOPPM - Microsoft Windows BMP file.
!
!    A PPM file can also be converted to other formats, by programs:
!
!      PPMTOACAD - AutoCAD file
!      PPMTOGIF  - GIF file
!      PPMTOPGM  - Portable Gray Map file
!      PPMTOPICT - Macintosh PICT file
!      PPMTOPUZZ - X11 puzzle file
!      PPMTORGB3 - 3 Portable Gray Map files
!      PPMTOXPM  - X11 pixmap file
!      PPMTOYUV  - Abekas YUV file
!
!  Example:
!
!    P3
!    # feep.ppma created by PBMLIB(PPMA_WRITE).
!    4 4
!    15
!     0  0  0    0  0  0    0  0  0   15  0 15
!     0  0  0    0 15  7    0  0  0    0  0  0
!     0  0  0    0  0  0    0 15  7    0  0  0
!    15  0 15    0  0  0    0  0  0    0  0  0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows
!    and columns of data.
!
!    Input, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
!    B(ROW_NUM,COL_NUM), the red, green and blue values of each pixel.  These
!    should be positive.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  logical, parameter :: debug = .false.
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) r(row_num,col_num)
  integer ( kind = 4 ) rgb_max

  ierror = 0
!
!  Compute the maximum color value.
!
  rgb_max = max ( &
    maxval ( r(1:row_num,1:col_num) ), &
    maxval ( g(1:row_num,1:col_num) ), &
    maxval ( b(1:row_num,1:col_num) ) )
!
!  Open the file.
!
  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    return
  end if
!
!  Write the header.
!
  call ppma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
    rgb_max, ierror )
!
!  Write the data.
!
  call ppma_write_data ( file_out_unit, row_num, col_num, r, g, b, ierror )
!
!  Close the file.
!
  close ( unit = file_out_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Note:'
    write ( *, '(a)' ) '  The data was checked and written.'
    write ( *, '(a,i8)' ) '  Number of data rows ROW_NUM =    ', row_num
    write ( *, '(a,i8)' ) '  Number of data columns COL_NUM = ', col_num
    write ( *, '(a,i8)' ) '  Maximum RGB value RGB_MAX =      ', rgb_max
  end if

  return
end
subroutine ppma_write_data ( file_out_unit, row_num, col_num, r, g, b, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE_DATA writes the data of a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows
!    and columns of data.
!
!    Input, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
!    B(ROW_NUM,COL_NUM), the red, green and blue values of each pixel.  These
!    should be positive.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) r(row_num,col_num)

  ierror = 0
!
!  Write the header.
!
  do i = 1, row_num
    do jlo = 1, col_num, 4
      jhi = min ( jlo + 3, col_num )
      write ( file_out_unit, '(12i5)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
    end do
  end do

  return
end
subroutine ppma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
  rgb_max, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE_HEADER writes the header of a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and
!    columns of data.
!
!    Input, integer ( kind = 4 ) RGB_MAX, the maximum value of any data component.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  character ( len = * )  file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ierror
  character ( len = 2 ) :: magic = 'P3'
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ) rgb_max

  ierror = 0
!
!  Write the header.
!
  write ( file_out_unit, '(a2)' ) magic
  write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
    // ' created by PPMA_IO::PPMA_WRITE.F90.'
  write ( file_out_unit, '(i5,2x,i5)' ) col_num, row_num
  write ( file_out_unit, '(i5)' ) rgb_max

  return
end
subroutine ppma_write_test ( file_out_name )

!*****************************************************************************80
!
!! PPMA_WRITE_TEST tests the ASCII portable pixel map write routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file
!    to contain the ASCII portable pixel map data.
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: b
  character ( len = * ) file_out_name
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: g
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: r

  row_num = 300
  col_num = 300
!
!  Allocate memory.
!
  allocate ( r(row_num,col_num) )
  allocate ( g(row_num,col_num) )
  allocate ( b(row_num,col_num) )
!
!  Set the data.
!
  call ppma_example ( row_num, col_num, r, g, b )
!
!  Write the data to the file.
!
  call ppma_write ( file_out_name, row_num, col_num, r, g, b, ierror )

  deallocate ( r );
  deallocate ( g );
  deallocate ( b );

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE_TEST - Fatal error!'
    write ( *, '(a)' ) '  PPMA_WRITE failed.'
  end if

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the value read from the string.
!    If blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character that was
!    part of the representation of IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lens
  character ( len = * ) s

  ierror = 0
  istate = 0

  isgn = 1
  ival = 0

  lens = len ( s )

  i = 0

  do

    i = i + 1

    c = s(i:i)

    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        exit
      end if

    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        exit
      end if

    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        istate = 3
      end if

    end if
!
!  Continue or exit?
!
    if ( istate == 3 ) then
      ival = isgn * ival
      last = i - 1
      exit
    else if ( lens <= i ) then
      if ( istate == 2 ) then
        ival = isgn * ival
        last = lens
      else
        ierror = 1
        last = 0
      end if
      exit
    end if

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine word_next_rd ( line, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_RD "reads" words from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    words separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE,
!      WORD contains the "next" word read from LINE.
!    Else
!      WORD is blank.
!
!    Input/output, logical DONE.
!    On input, on the first call, or with a fresh value of LINE,
!      set DONE to TRUE.
!    Else
!      leave it at the output value of the previous call.
!    On output, if a new nonblank word was extracted from LINE
!      DONE is FALSE
!    ELSE
!      DONE is TRUE.
!    If DONE is TRUE, then you need to provide a new LINE of data.
!
!  Local Parameters:
!
!    NEXT is the next location in LINE that should be searched.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lenl
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1
  character ( len = 1 ), parameter :: TAB = char(9)
  character ( len = * ) word

  lenl = len_trim ( line )

  if ( done ) then
    next = 1
    done = .false.
  end if
!
!  Beginning at index NEXT, search LINE for the next nonblank.
!
  ilo = next

  do
!
!  ...LINE(NEXT:LENL) is blank.  Return with WORD=' ', and DONE=TRUE.
!
    if ( lenl < ilo ) then
      word = ' '
      done = .true.
      next = lenl + 1
      return
    end if
!
!  ...If the current character is blank, skip to the next one.
!
    if ( line(ilo:ilo) /= ' ' .and. line(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  To get here, ILO must be the index of the nonblank starting
!  character of the next word.
!
!  Now search for the LAST nonblank character.
!
  next = ilo + 1

  do

    if ( lenl < next ) then
      word = line(ilo:next-1)
      return
    end if

    if ( line(next:next) == ' ' .or. line(next:next) == TAB ) then
      exit
    end if

    next = next + 1

  end do

  word = line(ilo:next-1)

  return
end
