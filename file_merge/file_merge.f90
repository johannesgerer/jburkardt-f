program main

!*****************************************************************************80
!
!! MAIN is the main program for FILE_MERGE.
!
!  Discussion:
!
!    FILE_MERGE merges two sorted files into a third.
!
!  Usage:
!
!    file_merge input1 input2 output
!
!  Discussion:
!
!    No check is made that the file is already sorted, and the program
!    will become seriously confused if it is not.
!
!    Duplicate entries are removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  character ( len = 255 ) dict_file_name_1
  character ( len = 255 ) dict_file_name_2
  character ( len = 255 ) dict_file_name_3
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  call timestamp ( )
!
!  Count the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the arguments from the command line or the user.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, dict_file_name_1 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter input file #1:'
    read ( *, '(a)' ) dict_file_name_1
  end if

  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, dict_file_name_2 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter input file #2:'
    read ( *, '(a)' ) dict_file_name_2
  end if

  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, dict_file_name_3 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter output file nme:'
    read ( *, '(a)' ) dict_file_name_3
  end if

  call merge_sorted_files ( dict_file_name_1, dict_file_name_2, &
    dict_file_name_3, n1, n2, n3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_MERGE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Words in file1: ', n1
  write ( *, '(a,i8)' ) '  Words in file2: ', n2
  write ( *, '(a,i8)' ) '  Words in file3: ', n3
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Duplicates:	', n1 + n2 - n3
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_MERGE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_low ( c )

!*****************************************************************************80
!
!! CH_LOW lowercases a single character.
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
!    Input/output, character C, the character to be lowercased.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 65 <= itemp .and. itemp <= 90 ) then
    c = char ( itemp + 32 )
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
subroutine merge_sorted_files ( file_name_1, file_name_2, file_name_3, &
  n1, n2, n3 )

!*****************************************************************************80
!
!! MERGE_SORTED_FILES merges two sorted files into a third.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME_1, FILE_NAME_2, the names of the
!    two input files to be merged.
!
!    Input, character ( len = * ) FILE_NAME_3, the name of the output file to
!    be created.
!
!    Output, integer ( kind = 4 ) N1, N2, N3, the number of lines of text in the
!    two input files and the output file.
!
  implicit none

  character ( len = * ) file_name_1
  character ( len = * ) file_name_2
  character ( len = * ) file_name_3
  integer ( kind = 4 ) file_unit_1
  integer ( kind = 4 ) file_unit_2
  integer ( kind = 4 ) file_unit_3
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  character ( len = 255 ) word1
  character ( len = 255 ) word2
  character ( len = 255 ) word3

  n1 = 0
  n2 = 0
  n3 = 0

  call get_unit ( file_unit_1 )

  open ( unit = file_unit_1, file = file_name_1, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MERGE_SORTED_FILES - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // trim ( file_name_1 ) // '".'
    return
  end if

  call get_unit ( file_unit_2 )

  open ( unit = file_unit_2, file = file_name_2, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MERGE_SORTED_FILES - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // trim ( file_name_2 ) // '".'
    return
  end if

  call get_unit ( file_unit_3 )

  open ( unit = file_unit_3, file = file_name_3, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MERGE_SORTED_FILES - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' // trim ( file_name_3 ) // '".'
    return
  end if

  word1 = ' '
  word2 = ' '

  do

    if ( word1 == ' ' ) then

      read ( file_unit_1, '(a)', iostat = ios ) word1

      if ( ios == 0 ) then
        call s_low ( word1 )
        n1 = n1 + 1
      end if

      if ( word1 == ' ' ) then
        word1 = '_END_'
        close ( unit = file_unit_1 )
      end if

    end if

    if ( word2 == ' ' ) then

      read ( file_unit_2, '(a)', iostat = ios ) word2

      if ( ios == 0 ) then
        call s_low ( word2 )
        n2 = n2 + 1
      end if

      if ( word2 == ' ' ) then
        word2 = '_END_'
        close ( unit = file_unit_2 )
      end if

    end if

    if ( word1 == '_END_' .and. word2 == '_END_' ) then

      exit

    else if ( word1 /= '_END_' .and. word2 == '_END_' ) then

      word3 = word1
      word1 = ' '

    else if ( word1 == '_END_' .and. word2 /= '_END_' ) then

      word3 = word2
      word2 = ' '

    else

      if ( llt ( word1, word2 ) ) then

        word3 = word1
        word1 = ' '

      else if ( word1 == word2 ) then

        word3 = word1
        word1 = ' '
        word2 = ' '

      else

        word3 = word2
        word2 = ' '

      end if

    end if

    write ( file_unit_3, '(a)' ) trim ( word3 )
    n3 = n3 + 1

  end do

  endfile ( unit = file_unit_3 )

  close ( unit = file_unit_3 )

  return
end
subroutine s_low ( s )

!*****************************************************************************80
!
!! S_LOW replaces all uppercase letters by lowercase ones.
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
!    Input/output, character ( len = * ) S, the string to be
!    transformed.  On output, the string is all lowercase.
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = * ) s

  do i = 1, len_trim ( s )
    call ch_low ( s(i:i) )
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
