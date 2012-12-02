program main

!*****************************************************************************80
!
!! MAIN is the main program for INCLUDE_FILES.
!
!  Discussion:
!
!    INCLUDE_FILES reads a FORTRAN file containing INCLUDE statements,
!    and makes a copy in which the indicated files have been included.
!
!    The format of an INCLUDE statement is presumed to be
!
!      include 'file_name'
!
!    or
!
!      include "file_name"
!
!    If the indicated file cannot be found, the copied file retains
!    the INCLUDE statement.  Otherwise, the INCLUDE statement is replaced
!    by the text of the indicated file.
!
!  Usage:
!
!    include_files file1.f file2.f
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) include_bad_num
  integer ( kind = 4 ) include_num
  character ( len = 255 ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) num_arg
  character ( len = 255 ) output_file
  integer ( kind = 4 ) output_unit

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INCLUDE_FILES'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a FORTRAN file with INCLUDE statements;'
  write ( *, '(a)' ) '  Make a copy in which those files are included.'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input file name:'
    read ( *, '(a)', iostat = ios ) input_file

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INCLUDE_FILES - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1
    call getarg ( iarg, input_file )

  end if
!
!  If two command line arguments, the second one is the output file name.
!
  if ( num_arg < 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the output file name:'
    read ( *, '(a)', iostat = ios ) output_file

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INCLUDE_FILES - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 2
    call getarg ( iarg, output_file )

  end if
!
!  Now we know what to do.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INCLUDE_FILES'
  write ( *, '(a)' ) '  Read file:  "' // trim ( input_file ) // '".'
  write ( *, '(a)' ) '  Write file: "' // trim ( output_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INCLUDE_FILES - Fatal error!'
    write ( *, '(a)' ) '  There was an error while trying to open the'
    write ( *, '(a)' ) '  input file "' // trim ( input_file ) // '".'
    stop
  end if

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace' )

  call include ( input_unit, output_unit, include_num, include_bad_num )

  close ( unit = input_unit )
  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(i6,a)' ) include_num, ' INCLUDE statements were encountered.'
  write ( *, '(a)' ) ' '
  write ( *, '(i6,a)' ) include_bad_num, ' times, the necessary INCLUDE file'
  write ( *, '(a)' ) '  could not be found.'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INCLUDE_FILES'
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
!    A "free" FORTRAN unit number is an integer ( kind = 4 ) between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer ( kind = 4 ) between 1 and 99, representing a
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
subroutine include ( input_unit, output_unit, include_num, &
  include_bad_num )

!*****************************************************************************80
!
!! INCLUDE makes a copy of a FORTRAN file with INCLUDE files included.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, OUTPUT_UNIT, the I/O units 
!    associated with the input and output files respectively.
!
!    Output, integer ( kind = 4 ) INCLUDE_NUM, the number of times an 
!    INCLUDE statement was encountered.
!
!    Output, integer ( kind = 4 ) INCLUDE_BAD_NUM, the number of times an 
!    INCLUDE statement was encountered, but the INCLUDE file could not be found.
!
  implicit none

  integer ( kind = 4 ) include_bad_num
  character ( len = 255 ) include_file
  integer ( kind = 4 ) include_num
  integer ( kind = 4 ) include_status
  integer ( kind = 4 ) include_unit
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = 255 ) line2
  integer ( kind = 4 ) n
  logical s_eqi
  integer ( kind = 4 ) output_unit

  include_bad_num = 0
  include_num = 0

  do
!
!  Read the next FORTRAN line.
!
    read ( input_unit, '(a)', iostat = input_status ) line
    
    if ( input_status /= 0 ) then
      exit
    end if

    line2 = line

    call s_blank_delete ( line2 )

    if ( .not. s_eqi ( line2(1:8), 'include"' ) .and. &
         .not. s_eqi ( line2(1:8), 'include''' ) ) then
      write ( output_unit, '(a)' ) trim ( line )
      cycle
    end if

    include_num = include_num + 1

    include_file = line2(9:)
    n = len_trim ( include_file )
    include_file(n:n) = ' '

    call get_unit ( include_unit )

    open ( unit = include_unit, file = include_file, &
      status = 'old', iostat = include_status )

    if ( include_status /= 0 ) then
      include_bad_num = include_bad_num + 1
      write ( output_unit, '(a)' ) trim ( line )
      cycle
    end if

    do

      read ( include_unit, '(a)', iostat = include_status ) line

      if ( include_status /= 0 ) then
        close ( unit = include_unit )
        exit
      end if

      write ( output_unit, '(a)' ) trim ( line )

    end do

  end do

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1998
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

  character c
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)
 
    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if
 
  end do
 
  s(put+1:nchar) = ' '
  
  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
    s3 = ' '
  else if ( s1 == ' ' ) then
    s3 = s2
  else if ( s2 == ' ' ) then
    s3 = s1
  else
    s3 = trim ( s1 ) // trim ( s2 )
  end if

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.  
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )
 
  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do
 
  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do
 
  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do
 
  s_eqi = .true.
 
  return
end
function s_index_last ( s, sub )

!*****************************************************************************80
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!  Discussion:
!
!    It returns the location in the string at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either S or SUB is all blanks.
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
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) llen1
  integer ( kind = 4 ) llen2
  character ( len = * ) s
  integer ( kind = 4 ) s_index_last
  character ( len = * ) sub

  s_index_last = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen2 > llen1 ) then
    return
  end if

  do j = 1, llen1+1-llen2

    i = llen1 + 2 - llen2 - j

    if ( s(i:i+llen2-1) == sub ) then
      s_index_last = i
      return
    end if

  end do

  return
end
subroutine s_tab_blanks ( s )

!*****************************************************************************80
!
!! S_TAB_BLANKS replaces TAB characters by 6 spaces.
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
!    Input/output, character ( len = * ) S, the string to be modified.  On
!    output, some significant characters at the end of S may have
!    been lost.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) ntab
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
!  If no TAB's occur in the line, there is nothing to do.
!
  if ( index ( s, TAB ) == 0 ) then
    return
  end if
!
!  Otherwise, find out how long the string is.
!
  lenc = len_trim ( s )
  lens = len ( s )
!
!  Count the TAB's.
!
  ntab = 0
  do i = 1, lenc
    if ( s(i:i) == TAB ) then
      ntab = ntab + 1
    end if
  end do
!
!  Now copy the string onto itself, going backwards.
!  As soon as we've processed the first TAB, we're done.
!
  iput = lenc + 5 * ntab

  do iget = lenc, 1, - 1

    if ( s(iget:iget) /= TAB ) then

      if ( iput <= lens ) then
        s(iput:iput) = s(iget:iget)
      end if

      iput = iput - 1

    else

      do i = iput, iput - 5, -1
        if ( i <= lens ) then
          s(i:i) = ' '
        end if
      end do

      iput = iput - 6
      ntab = ntab - 1

      if ( ntab == 0 ) then
        return
      end if

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
!    31 May 2001
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
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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
