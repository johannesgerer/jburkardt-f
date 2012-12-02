program main

!*****************************************************************************80
!
!! MAIN is the main program for QUOTE.
!
!  Discussion:
!
!    QUOTE extracts a random quote from a file.
!
!    This version of QUOTE is descended from a program written
!    by David Moses.
!
!    The quote file contains a series of quotes, separated by a blank line.
!    When printing out a quote from the file, any line that ends with a
!    quotation mark is followed by an extra blank line.
!
!  Usage:
!
!    quote
!
!      extracts a quote from the default quote file.
!
!    quote FILE
!
!      extracts a quote from the file FILE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  character ( len = 255 ) line
  integer ( kind = 4 ) line_index
  integer ( kind = 4 ) line_num
  character ( len = 255 ) list_file_name
  integer ( kind = 4 ) numarg
  character ( len = 255 ) quote_file_name
  character ( len = 255 ) quote_file_title
  integer ( kind = 4 ) quote_index
  integer ( kind = 4 ) quote_num
  integer ( kind = 4 ) seed

  if ( DEBUG ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUOTE:'
    write ( *, '(a)' ) '  FORTRAN90 version'
  end if

  write ( *, '(a)' ) ' '

  call timestamp ( )
!
!  Get a seed for the random number generator based on the current time.
!
  call get_seed ( seed )
!
!  Count the command line arguments.
!
  numarg = iargc ( )
!
!  If a quote file was not specified, open QUOTE_FILES.TXT and
!  pick a file at random.
!
  if ( numarg == 0 ) then

!   list_file_name = &
!     '/a/fs.csit.fsu.edu/u8/users/burkardt/public_html/f_src/quote/quote_files.txt'
!   list_file_name = &
!     '/Users/burkardt/public_html/f_src/quote/quote_files.txt'
!   list_file_name = &
!     '/Users/jburkardt/public_html/f_src/quote/quote_files.txt'
    list_file_name = &
      '/panfs/panasas1/users/jburkardt/public_html/f_src/quote/quote_files.txt'

    call file_line_count ( list_file_name, line_num )

    if ( line_num < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUOTE - Fatal error.'
      write ( *, '(a)' ) '  The list of quote files is missing!'
      stop
    end if

    if ( DEBUG ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUOTE: Debug'
      write ( *, '(a,i8)' ) &
        '  The number of quote file listing entries is ', line_num
    end if

    line_index = i4_uniform ( 1, line_num, seed )

    if ( DEBUG ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUOTE: Debug'
      write ( *, '(a,i8)' ) '  Choosing file number ', line_index
    end if

    call file_line_get ( list_file_name, line_index, line )

    call word_next2 ( line, quote_file_name, quote_file_title )

    if ( DEBUG ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUOTE: Debug'
      write ( *, '(a)' ) trim ( quote_file_title )
    end if

  else

    iarg = 1
    call getarg ( iarg, quote_file_name )

    quote_file_title = 'One of my favorite quotations.'

  end if
!
!  Compute QUOTE_NUM, the number of quotes in the file.
!
  call file_para_count ( quote_file_name, quote_num )

  if ( quote_num <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUOTE - Fatal error!'
    write ( *, '(a)' ) '  FILE_PARA_COUNT could not count the quotes.'
    stop
  end if
!
!  Set QUOTE_INDEX, the number of the quote to be displayed.
!
  quote_index = i4_uniform ( 1, quote_num, seed )
 
  if ( quote_index <= 0 ) then
    quote_index = 1
  end if
 
  if ( quote_num < quote_index ) then
    quote_index = quote_num
  end if

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUOTE: Debug'
    write ( *, '(a,i8)' ) '  Extract quote number ', quote_index
  end if
!
!  Read from the quote file until you reach the desired quote.
!  QUOTE_INDEX tells us which quote we are currently reading.
!  Quotes are presumed to be separated by a single blank line.
!
  call quote_file_print ( quote_file_name, quote_file_title, quote_index, &
    ierror )
!
!  Terminate.
!
  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUOTE:'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  stop
end
subroutine file_line_count ( file_name, line_num )

!*****************************************************************************80
!
!! FILE_LINE_COUNT counts the number of lines in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer LINE_NUM, the number of lines found in the file.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) line_num

  line_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    line_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not find a free FORTRAN unit.'
    return
  end if

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    line_num = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' // trim ( file_name ) // '".'
    return
  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios )

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

  end do

  close ( unit = iunit )

  return
end
subroutine file_line_get ( file_name, line_index, line )

!*****************************************************************************80
!
!! FILE_LINE_GET gets a particular line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input, integer LINE_INDEX, the line number to be read.
!
!    Output, character ( len = * ) LINE, the text of the line.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = * ) line
  integer ( kind = 4 ) line_index
!
!  Open the file.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    line = ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_GET - Fatal error!'
    write ( *, '(a)' ) '  Could not find a free FORTRAN unit.'
    return
  end if

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    line = ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_GET - Fatal error!'
    write ( *, '(a)' ) '  Could not open the the file "' // trim ( file_name ) // '".'
    return
  end if
!
!  Count the lines.
!
  do i = 1, line_index

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      line = ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_LINE_GET - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file.'
      return
    end if

  end do

  close ( unit = iunit )

  return
end
subroutine file_para_count ( file_name, para_num )

!*****************************************************************************80
!
!! FILE_PARA_COUNT counts the number of paragraphs in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  A paragraph is
!    a sequence of nonblank lines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer PARA_NUM, the number of paragraphs found in the file.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lenc_old
  character ( len = 255 ) line
  integer ( kind = 4 ) para_num

  para_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    para_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_PARA_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not find a free FORTRAN unit.'
    return
  end if

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    para_num = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_PARA_COUNT - Fatal error!'
    write ( *, '(a)') '  Could not open the the file "' // trim ( file_name ) // '".'
    return
  end if
!
!  Count the paragraphs.
!
  lenc = 0

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    lenc_old = lenc
    lenc = len_trim ( line )

    if ( 0 < lenc .and. lenc_old <= 0 ) then
      para_num = para_num + 1
    end if

  end do

  close ( unit = iunit )

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) / 11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) / 30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) / 23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp / 6.0D+00
!
!  Force 0 < TEMP <= 1.
!
  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == huge ( 1 ) ) then
    seed = seed - 1
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
!    Output, integer IUNIT, the free unit number.
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
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine quote_file_print ( quote_file_name, quote_file_title, quote_index, &
  ierror )

!*****************************************************************************80
!
!! QUOTE_FILE_PRINT prints a given quote from a quote file.
!
!  Discussion:
!
!    The index of the desired quote is specified.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) QUOTE_FILE_NAME, the name of the quote file.
!
!    Input, character ( len = * ) QUOTE_FILE_TITLE, the title of the quote file.
!
!    Input, integer ( kind = 4 ) QUOTE_INDEX, the index of the quote to be printed.
!    The first quote has QUOTE_INDEX = 1, and so on.  If there are fewer
!    that QUOTE_INDEX quotes in the file, then nothing will be printed.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error was encountered.
!    nonzero, an error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) jquote
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lencold
  character ( len = 255 ) line
  character ( len = * ) quote_file_name
  character ( len = * ) quote_file_title
  integer ( kind = 4 ) quote_file_unit
  integer ( kind = 4 ) quote_index
!
!  Open the quote file.
!
  call get_unit ( quote_file_unit )

  if ( quote_file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUOTE_FILE_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Could not find a free FORTRAN unit.'
    return
  end if

  open ( unit = quote_file_unit, file = quote_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUOTE_FILE_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the quote the file "' &
      // trim ( quote_file_name ) // '".'
    return
  end if

  jquote = 1
  lenc = 0
 
  write ( *, '(a)' ) ' '
 
  do
 
    read ( quote_file_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then

      if ( quote_index /= jquote ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUOTE_FILE_PRINT - Fatal error!'
        write ( *, '(a)' ) '  Could not find the desired quote.'
      end if
      close ( unit = quote_file_unit )
      return
    end if

    lencold = lenc
    lenc = len_trim ( line )
 
    if ( 0 < lenc .and. lencold <= 0 ) then
      jquote = jquote + 1
    end if

    if ( jquote < quote_index ) then

    else if ( jquote == quote_index ) then

      write ( *, '(2x,a)' ) trim ( line )
      if ( 0 < lenc ) then
        if ( line(lenc:lenc) == '"' ) then
          write ( *, '(a)' ) ' '
        end if
      end if
 
    else if ( quote_index < jquote ) then

      exit

    end if
 
  end do
 
  write ( *, '(2x,a,i6)' ) 'Quotation ', jquote
  write ( *, '(2x,a)' ) trim ( quote_file_title )
  write ( *, '(a)' ) ' '

  close ( unit = quote_file_unit )

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
subroutine word_next2 ( s, first, last )

!*****************************************************************************80
!
!! WORD_NEXT2 returns the first word in a string.
!
!  Discussion:
!
!    "Words" are any string of characters, separated by commas or blanks.
!
!    The routine returns:
!    * FIRST, the first string of nonblank, noncomma characters;
!    * LAST, the characters of the string that occur after FIRST and
!      the commas and blanks.
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
!    Input/output, character ( len = * ) S, the string of words to be analyzed.
!
!    Output, character ( len = * ) FIRST, the next word in the string.
!
!    Output, character ( len = * ) LAST, the remaining string.
!
  implicit none

  character c
  character ( len = * ) first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) ilast
  character ( len = * ) last
  integer ( kind = 4 ) lenf
  integer ( kind = 4 ) lenl
  integer ( kind = 4 ) lens
  character ( len = * ) s

  first = ' '
  last = ' '

  ifirst = 0
  ilast = 0

  lens = len_trim ( s )
  lenf = len ( first )
  lenl = len ( last )

  ido = 0

  do i = 1, lens

    c = s(i:i)

    if ( ido == 0 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ido = 1
      end if
    end if

    if ( ido == 1 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ifirst = ifirst + 1
        if ( ifirst <= lenf ) then
          first(ifirst:ifirst) = c
        end if
      else
        ido = 2
      end if
    end if

    if ( ido == 2 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ido = 3
      end if
    end if

    if ( ido == 3 ) then
      ilast = ilast + 1
      if ( ilast <= lenl ) then
        last(ilast:ilast) = c
      end if
    end if

  end do

  return
end
