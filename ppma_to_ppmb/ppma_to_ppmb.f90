program main

!*****************************************************************************80
!
!! MAIN is the main program of PPMA_TO_PPMB.
!
!  Discussion:
!
!    PPMA_TO_PPMB converts an ASCII PPM file to binary PPM format.
!
!    PPMA_TO_PPMB is a sample application of the PBMLIB library.
!
!    It calls on the PBMLIB PPMA_READ routine to open and read a
!    user-specified file in the ASCII PPM format.
!
!    It then calls on the PBMLIB PPMB_WRITE routine to create and write
!    a copy of the data in a file in the binary PPM format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Usage:
!
!    ppma_to_ppmb file.ppma file.ppmb
!
!  Parameters:
!
!    Character FILE.PPMA, is the name of the input ASCII PPM file to read.
!
!    Character FILE.PPMB, is the name of the output binary PPM file to write.
!
!
!  MAXP is the maximum number of pixels that the program can handle.
!  Simply increase the value of MAXP to handle larger files.
!
  integer ( kind = 4 ), parameter :: maxp = 500000

  integer ( kind = 4 ) b(maxp)
  integer ( kind = 4 ) g(maxp)
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  character ( len = 255 ) input_file_name
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) numarg
  character ( len = 255 ) output_file_name
  integer ( kind = 4 ) r(maxp)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PPMA_TO_PPMB'
  write ( *, '(a)' ) '  FORTRAN90 version'

  ierror = 0

  numarg = iargc ( )

  if ( numarg < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_TO_PPMB - Error!'
    write ( *, '(a)' ) '  Usage:  ppma_to_ppmb file.ppma file.ppmb'
    stop
  end if
!
!  Read the input file.
!
  iarg = 1

  call getarg ( iarg, input_file_name )

  call ppma_read ( input_file_name, ierror, maxcol, maxp, nrow, ncol, r, g, b  )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_TO_PPMB - Error!'
    write ( *, '(a,i6)' ) '  PPMA_READ returns IERROR = ', ierror
  end if
!
!  Write the output file.
!
  iarg = 2

  call getarg ( iarg, output_file_name )

  call ppmb_write ( output_file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_TO_PPMB - Error!'
    write ( *, '(a,i6)' ) '  PPMB_WRITE returns IERROR = ', ierror
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PPMA_TO_PPMB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
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
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
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
!
!    On input, if this is the first call, or the user has changed
!    STRING, then set DONE = TRUE.
!
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
subroutine ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

!*****************************************************************************80
!
!! PPM_CHECK_DATA checks pixel data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contains the
!    RGB pixel data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, the data is illegal.
!
!    Input, integer ( kind = 4 ) MAXCOL, the maximum value.
!
!    Input, integer ( kind = 4 ) NCOL, NROW, the number of rows and columns of data.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) r(nrow,ncol)

  ierror = 0
!
!  Make sure no color is negative.
!
  if ( minval ( r(1:nrow,1:ncol) ) < 0 .or. &
       minval ( g(1:nrow,1:ncol) ) < 0 .or. &
       minval ( b(1:nrow,1:ncol) ) < 0 ) then
    ierror = 1
    return
  end if
!
!  Make sure no color is greater than MAXCOL.
!
  if ( maxcol < maxval ( r(1:nrow,1:ncol) ) .or. &
       maxcol < maxval ( g(1:nrow,1:ncol) ) .or. &
       maxcol < maxval ( b(1:nrow,1:ncol) ) ) then
    ierror = 1
    return
  end if

  return
end
subroutine ppma_read ( file_name, ierror, maxcol, maxp, nrow, ncol, r, g, b )

!*****************************************************************************80
!
!! PPMA_READ reads an ASCII portable pixel map file.
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
!    # feep.ppm
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
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file from which
!    the data should be read.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error occurred.
!    1, the file could not be opened.
!    2, end or error while reading file.
!    3, bad magic number (the first two bytes must be 'P3').
!    4, trouble reading NROW or NCOL or MAXCOL.
!    5, trouble reading one of the pixel values.
!    6, at least one pixel value was less than 0 or greater than MAXCOL.
!    7, NROW*NCOL exceeds MAXP.
!
!    Output, integer ( kind = 4 ) MAXCOL, the maximum pixel color value.
!
!    Input, integer ( kind = 4 ) MAXP, the number of entries available in R, G and B.
!    If MAXP is smaller than NROW*NCOL, then the data will not be read.
!
!    Output, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns of data.
!
!    Output, integer ( kind = 4 ) R(MAXP), contains the NROW by NCOL red data.
!    The (I,J) entry is in R( (J-1)*NROW + I ), the usual
!    FORTRAN indexing method.
!
!    Output, integer ( kind = 4 ) G(MAXP), contains the NROW by NCOL green data.
!    The (I,J) entry is in G( (J-1)*NROW + I ), the usual
!    FORTRAN indexing method.
!
!    Output, integer ( kind = 4 ) B(MAXP), contains the NROW by NCOL blue data.
!    The (I,J) entry is in B( (J-1)*NROW + I ), the usual
!    FORTRAN indexing method.
!
  implicit none

  integer ( kind = 4 ) maxp

  integer ( kind = 4 ) b(maxp)
  logical, parameter :: debug = .false.
  logical done
  character ( len = * ) file_name
  integer ( kind = 4 ) g(maxp)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 2 ) magic
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) r(maxp)
  logical s_eqi
  character ( len = 80 ) string

  ierror = 0
  maxcol = 0
  ncol = 0
  nrow = 0
!
!  Open the file.
!
  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 1
    return
  end if
!
!  Read the first line of data, which must begin with the magic number.
!
  read ( file_unit, '(a)', iostat = ios ) magic

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
    write ( *, '(a)' ) '  End or error while reading file.'
    ierror = 2
    return
  end if

  if ( .not. s_eqi ( magic, 'P3' ) ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Fatal error.'
    write ( *, '(a)' ) '  First two bytes are not magic number "P3".'
    write ( *, '(a)' ) '  First two bytes are: ' // magic
    return
  end if
!
!  Now search for NCOL, NROW, and MAXCOL.
!
  done = .TRUE.
  string = ' '

  call getint ( done, ierror, file_unit, ncol, string )

  if ( ierror /= 0 ) then
    close ( unit = file_unit )
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
    write ( *, '(a)' ) '  Problem reading NCOL.'
    return
  end if

  call getint ( done, ierror, file_unit, nrow, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
    write ( *, '(a)' ) '  Problem reading NROW.'
    return
  end if

  call getint ( done, ierror, file_unit, maxcol, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
    write ( *, '(a)' ) '  Problem reading MAXCOL.'
    return
  end if
!
!  Check that there is enough room.
!
  if ( maxp < nrow * ncol ) then
    ierror = 7
    close ( unit = file_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
    write ( *, '(a,i12)' ) '  Needed NROW*NCOL = ', nrow * ncol
    write ( *, '(a,i12)' ) '  Available MAXP = ', maxp
    return
  end if
!
!  Now read the pixel data.
!
  k = 0
  do i = 1, nrow
    do j = 1, ncol

      k = ( j - 1 ) * nrow + i

      call getint ( done, ierror, file_unit, r(k), string )

      if ( ierror /= 0 ) then
        ierror = 5
        close ( unit = file_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
        write ( *, '(a)' ) '  Problem reading R data.'
        return
      end if

      call getint ( done, ierror, file_unit, g(k), string )

      if ( ierror /= 0 ) then
        ierror = 5
        close ( unit = file_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
        write ( *, '(a)' ) '  Problem reading G data.'
        return
      end if

      call getint ( done, ierror, file_unit, b(k), string )

      if ( ierror /= 0 ) then
        ierror = 5
        close ( unit = file_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
        write ( *, '(a)' ) '  Problem reading B data.'
        return
      end if

    end do
  end do
!
!  Close the file.
!
  close ( unit = file_unit )
!
!  Check the data.
!
  call ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Fatal error!'
    write ( *, '(a)' ) '  Bad data detected by PPM_CHECK_DATA.'
    ierror = 6
    return
  end if
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_READ - Note:'
    write ( *, '(a)' ) '  The file was read and checked.'
    write ( *, '(a,i6)' ) '  Number of data rows NROW =    ', nrow
    write ( *, '(a,i6)' ) '  Number of data columns NCOL = ', ncol
    write ( *, '(a,i6)' ) '  Maximum color value MAXCOL =  ', maxcol
  end if

  return
end
subroutine ppmb_write ( file_name, ierror, nrow, ncol, r, g, b )

!*****************************************************************************80
!
!! PPMB_WRITE writes a binary portable pixel map file.
!
!  Discussion:
!
!    PPM files can be viewed by XV.
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
!      PPMTORGB3 - 3 Portable Gray Map files
!      PPMTOXPM  - X11 pixmap file
!      PPMTOYUV  - Abekas YUV file
!
!    DIRECT ACCESS is used for the output file just so that we can
!    avoid the internal carriage returns and things that FORTRAN
!    seems to want to add.
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
!    Input, character ( len = * ) FILE_NAME, the name of the file to which
!    the data should be written.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contain
!    the red, green and blue values of each pixel.  These should all
!    be values between 0 and 255.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  logical, parameter :: debug = .false.
  character ( len = * ) file_name
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) istring(17)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) nval
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) r(nrow,ncol)
  integer ( kind = 4 ) record
  integer ( kind = 4 ) record_length
  logical, parameter :: reverse = .false.
  character ( len = 68 ) string

  ierror = 0
!
!  Compute the maximum color value.
!
  maxcol = max ( &
    maxval ( r(1:nrow,1:ncol) ), &
    maxval ( g(1:nrow,1:ncol) ), &
    maxval ( b(1:nrow,1:ncol) ) )
!
!  Check that no color data exceeds 255.
!
  if ( 255 < maxcol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMB_WRITE - Fatal error!'
    write ( *, '(a)' ) '  The color data exceeds 255!'
    write ( *, '(a,i12)' ) '  MAXCOL = ', maxcol
    ierror = 1
    return
  end if
!
!  Check the data.
!
  call ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMB_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Bad data detected by PPM_CHECK_DATA!'
    ierror = 1
    return
  end if
!
!  Open the file.
!
!  The smallest amount of information we can write at a time is
!  1 word = 4 bytes = 32 bits.
!
  call get_unit ( output_unit )
!
!  For the SGI:
!
  record_length = 4
!
!  For the DEC Alpha:
!
! record_length = 1

  open ( unit = output_unit, file = file_name, status = 'replace', &
    form = 'unformatted', access = 'direct', recl = record_length, &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMB_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    return
  end if

  record = 0
!
!  Write the data.
!
  string = ' '

  string(1:4) = 'P6  '

  write ( string(5:9), '(i5)' ) ncol

  string(10:11) = ' '

  write ( string(12:16), '(i5)' ) nrow

  string(17:18) = ' '

  write ( string(19:23), '(i5)' ) maxcol

  string(24:24) = ' '

  nchar = 24

  do i = 1, nrow
    do j = 1, ncol

      if ( nchar == 68 ) then
        call s_to_i4vec ( string(1:nchar), nval, istring, reverse )
        do l = 1, nval
          record = record + 1
          write ( output_unit, rec = record ) istring(l)
        end do
        string = ' '
        nchar = 0
      end if

      nchar = nchar + 1
      string(nchar:nchar) = char ( r(i,j) )

      if ( nchar == 68 ) then
        call s_to_i4vec ( string(1:nchar), nval, istring, reverse )
        do l = 1, nval
          record = record + 1
          write ( output_unit, rec = record ) istring(l)
        end do
        string = ' '
        nchar = 0
      end if

      nchar = nchar + 1
      string(nchar:nchar) = char ( g(i,j) )

      if ( nchar == 68 ) then
        call s_to_i4vec ( string(1:nchar), nval, istring, reverse )
        do l = 1, nval
          record = record + 1
          write ( output_unit, rec = record ) istring(l)
        end do
        string = ' '
        nchar = 0
      end if

      nchar = nchar + 1
      string(nchar:nchar) = char ( b(i,j) )


    end do
  end do

  if ( 0 < nchar ) then
    call s_to_i4vec ( string(1:nchar), nval, istring, reverse )
    do l = 1, nval
      record = record + 1
      write ( output_unit, rec = record ) istring(l)
    end do
    string = ' '
    nchar = 0
  end if
!
!  Close the file.
!
  close ( unit = output_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMB_WRITE - Note:'
    write ( *, '(a)' ) '  The data was checked and written.'
    write ( *, '(a,i6)' ) '  Number of words =             ', record
    write ( *, '(a,i6)' ) '  Number of data rows NROW =    ', nrow
    write ( *, '(a,i6)' ) '  Number of data columns NCOL = ', ncol
    write ( *, '(a,i6)' ) '  Maximum color value MAXCOL =  ', maxcol
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
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
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
subroutine s_to_i4vec ( s, n, ivec, reverse )

!*****************************************************************************80
!
!! S_TO_I4VEC converts an string of characters into an I4VEC.
!
!  Discussion:
!
!    This routine can be useful when trying to write character data to an
!    unformatted direct access file.
!
!    Depending on the internal byte ordering used on a particular machine,
!    the parameter REVERSE_ORDER may need to be set TRUE or FALSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string of characters.
!    Each set of 4 characters is assumed to represent an integer.
!
!    Output, integer ( kind = 4 ) N, the number of integers read from the string.
!
!    Output, integer ( kind = 4 ) IVEC(*), an array of N integers which contains
!    the information from S.
!
!    Input, logical REVERSE, is TRUE if the bytes in a word need to be
!    reversed.
!
  implicit none

  integer ( kind = 4 ) from
  integer ( kind = 4 ) frompos
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ivec(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: length = 8
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nchar
  logical reverse
  character ( len = * ) s
  integer ( kind = 4 ) to
  integer ( kind = 4 ) topos

  nchar = len ( s )
  n = 0
  frompos = 0

  do ilo = 1, nchar, 4

    n = n + 1
    ihi = min ( ilo + 3, nchar )
    to = 0

    do j = ilo, ihi

      from = ichar ( s(j:j) )

      if ( reverse ) then
        topos = length * ( j - ilo )
      else
        topos = length * ( ilo + 3 - j )
      end if

      call mvbits ( from, frompos, length, to, topos )

    end do

    ivec(n) = to

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
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

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
!    06 January 2009
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
!
!    If DONE is FALSE,
!      WORD contains the "next" word read from LINE.
!    Else
!      WORD is blank.
!
!    Input/output, logical DONE.
!
!    On input, on the first call, or with a fresh value of LINE,
!      set DONE to TRUE.
!    Else
!      leave it at the output value of the previous call.
!
!    On output, if a new nonblank word was extracted from LINE
!      DONE is FALSE
!    ELSE
!      DONE is TRUE.
!
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
