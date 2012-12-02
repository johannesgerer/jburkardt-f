program main

!*****************************************************************************80
!
!! MAIN is the main program for MOVIE_DATA_REFORMAT.
!
!  Discussion:
!
!    This is an interactive program which allows a user to reformat
!    movie data copied from THE NUMBERS web page.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) agross
  integer ( kind = 4 ) arg_num
  character ( len = 80 ) distributor
  character ( len = 80 ) genre
  character ( len = 80 ) gross
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) input_count
  character ( len = 255 ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 80 ) movie
  character ( len = 80 ) mpaa
  integer ( kind = 4 ) output_count
  character ( len = 255 ) output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) rank
  character ( len = 80 ) release
  character ( len = 80 ) tickets

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MOVIE_DATA_REFORMAT:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Convert a text file of movie data into a CSV file.'
  write ( *, '(a)' ) '  The original data comes from the THE NUMBERS web page:'
  write ( *, '(a)' ) '  http://www.the-numbers.com.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A browser can copy the tabular movie data for one year'
  write ( *, '(a)' ) '  but stores the copied data as a list, one item per line.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program puts all the data for a movie on one line,'
  write ( *, '(a)' ) '  separating data items by commas, quoting string data,'
  write ( *, '(a)' ) '  removing dollar signs and commas from numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Empty string data is set to "?".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It also removes control characters.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, input_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the name of the input file.'

    read ( *, '(a)' ) input_filename

  end if
!
!  Create the output filename.
!
  output_filename = input_filename
  call file_name_ext_swap ( output_filename, 'csv' )
!
!  Open the files.
!
  call get_unit ( input_unit )
  open ( unit = input_unit, file = input_filename, status = 'old' )

  call get_unit ( output_unit )
  open ( unit = output_unit, file = output_filename )
!
!  Read 9 data items from input, write one data record to output.
!
  write ( output_unit, '(a)' ) &
    '"Rank", "Title", "Release", "Distributor", "Genre", "MPAA", "Gross", "Tickets", "Adjusted Gross"'

  input_count = 0
  output_count = 0

  do

    read ( input_unit, *, iostat = ios ) rank

    if ( ios /= 0 ) then
      exit
    end if

    read ( input_unit, '(a80)' ) movie

    read ( input_unit, '(a80)' ) release

    read ( input_unit, '(a80)' ) distributor
    call s_control_blank ( distributor )
    if ( len_trim ( distributor ) == 0 ) then
      distributor = '?'
    end if

    read ( input_unit, '(a80)' ) genre
    call s_control_blank ( genre )
    if ( len_trim ( genre ) == 0 ) then
      genre = '?'
    end if

    read ( input_unit, '(a80)' ) mpaa
    call s_control_blank ( mpaa )
    if ( len_trim ( mpaa ) == 0 ) then
      mpaa = '?'
    end if

    read ( input_unit, '(a80)' ) gross
    call s_ch_delete ( gross, '$' )
    call s_ch_delete ( gross , ',' )

    read ( input_unit, '(a80)' ) tickets
    call s_ch_delete ( tickets , ',' )

    read ( input_unit, '(a80)' ) agross
    call s_ch_delete ( agross, '$' )
    call s_ch_delete ( agross , ',' )

    input_count = input_count + 9

    write ( output_unit, '(i3,a,a,a,a,a,a,a,a,a)' ) &
      rank, ',', &
      '"' // trim ( movie ) // '",', &
      '"' // trim ( release ) // '",', &
      '"' // trim ( distributor ) // '",', &
      '"' // trim ( genre ) // '",', &
      '"' // trim ( mpaa ) // '",', &
      trim ( gross ) // ',', &
      trim ( tickets ) // ',' , &
      trim ( agross )

    output_count = output_count + 1

  end do

  close ( unit = input_unit )
  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(i6,a)' ) input_count, ' movie data items read from "' &
    // trim ( input_filename ) // '".'
  write ( *, '(a,i6,a)' ) '1 header line and ', output_count, &
    ' output movie records were written to "' // trim ( output_filename ) &
    // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MOVIE_DATA_REFORMAT:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function ch_index_last ( s, ch )

!*****************************************************************************80
!
!! CH_INDEX_LAST is the last occurrence of a character in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character CH, the character to be searched for.
!
!    Output, integer ( kind = 4 ) CH_INDEX_LAST, the location of the last 
!    occurrence of the character in the string, or -1 if it does not occur.
!
  implicit none

  character ch
  integer ( kind = 4 ) ch_index_last
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  ch_index_last = -1
  s_length = len_trim ( s )

  do i = s_length, 1, -1

    if ( s(i:i) == ch ) then
      ch_index_last = i
      return
    end if

  end do
 
  return
end
function ch_is_control ( ch )

!*****************************************************************************80
!
!! CH_IS_CONTROL is TRUE if a character is a control character.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
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
!    Input, character CH, the character to be tested.
!
!    Output, logical CH_IS_CONTROL, TRUE if the character is a control
!    character, and FALSE otherwise.
!
  implicit none

  character ch
  logical ch_is_control

  if ( iachar ( ch ) <= 31 .or. 127 <= iachar ( ch ) ) then
    ch_is_control = .true.
  else
    ch_is_control = .false.
  end if

  return
end
subroutine file_name_ext_get ( file_name, i, j )

!*****************************************************************************80
!
!! FILE_NAME_EXT_GET determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILE_NAME   I  J
!
!    bob.for     4  7
!    N.B.C.D     6  7
!    Naomi.      6  6
!    Arthur     -1 -1 
!    .com        1  1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and 
!    last characters in the file extension.  
!    If no period occurs in FILE_NAME, then
!      I = J = -1;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ch_index_last

  i = ch_index_last ( file_name, '.' )

  if ( i == -1 ) then

    j = -1

  else

    j = len_trim ( file_name )

  end if

  return
end
subroutine file_name_ext_swap ( file_name, ext )

!*****************************************************************************80
!
!! FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!  Example:
!
!          Input           Output
!    ================     =========
!    FILE_NAME    EXT     FILE_NAME
!
!    bob.for      obj     bob.obj
!    bob.bob.bob  txt     bob.bob.txt
!    bob          yak     bob.yak
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME, a file name.
!    On output, the extension of the file has been changed.
!
!    Input, character ( len = * ) EXT, the extension to be used on the output
!    copy of FILE_NAME, replacing the current extension if any.
!
  implicit none

  character ( len = * ) ext
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len_max
  integer ( kind = 4 ) len_name

  len_max = len ( file_name )
  len_name = len_trim ( file_name )

  call file_name_ext_get ( file_name, i, j )

  if ( i == -1 ) then

    if ( len_max < len_name + 1 ) then
      return
    end if

    len_name = len_name + 1
    file_name(len_name:len_name) = '.'
    i = len_name + 1

  else

    i = i + 1
    file_name(i:j) = ' '

  end if

  file_name(i:) = ext

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
!    26 October 2008
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
subroutine s_ch_delete ( s, ch )

!*****************************************************************************80
!
!! S_CH_DELETE removes all occurrences of a character from a string.
!
!  Discussion:
!
!    Each time the given character is found in the string, the characters
!    to the right of the string are shifted over one position.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
!    Input, character CH, the character to be removed.
!
  implicit none

  character ch
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

  put = 1

  do get = 1, s_length

    if ( s(get:get) == ch ) then

    else if ( put == get ) then
      put = put + 1
    else
      s(put:put) = s(get:get)
      put = put + 1
    end if

  end do

  s(put:s_length) = ' '

  return
end
subroutine s_control_blank ( s )

!*****************************************************************************80
!
!! S_CONTROL_BLANK replaces control characters with blanks.
!
!  Discussion:
!
!    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
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
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  logical ch_is_control
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length
    if ( ch_is_control ( s(i:i) ) ) then
      s(i:i) = ' '
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
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2008
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
