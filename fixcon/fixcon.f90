program main

!*****************************************************************************80
!
!! MAIN is the main program for FIXCON.
!
!  Discussion:
!
!    FIXCON gets information from a user about the FORTRAN77 file to convert.
!
!    The program's task is simple but tedious.  FORTRAN77 imposed a maximum
!    line length of 72 characters (and also reserved the first 6 columns for
!    labels and continuation marks).  So a FORTRAN77 programmer, when writing
!    a long statement, would have to break the statement up in a way that
!    fitted within the column limitations, while signaling to the compiler
!    that these several statements were really one.  FORTRAN77 used a system
!    of "continuation" characters in this situation.  If a statement was
!    continued over several lines, then a special character was put in
!    column 6 of all the lines which represented continuations of the previous
!    line.  The continuation character could be any legal FORTRAN77 character
!    EXCEPT a blank or the character '0'.
!
!    Thus, a long FORTRAN77 statement, when broken up using continuation,
!    might look like this:
!
!      x = 1 + 2 + 3 + 4 + 
!     &    5 + 6 + 7 + 8 +
!     &    9 + 10
!
!    where we are here using "&" for the continuation character.
!
!    In FORTRAN90, a long line is continued by putting an ampersand (this is
!    the ONLY continuation character FORTRAN90 allows) at the END of each 
!    line for which there is more coming.  So the same long statement
!    in FORTRAN90 might look like
!
!      x = 1 + 2 + 3 + 4 + &
!          5 + 6 + 7 + 8 + &
!          9 + 10
!
!    This program reads in a file that uses the FORTRAN77 continuation scheme,
!    and writes out a copy that uses the FORTRAN90 continuation scheme.
!
!    FORTRAN77 ignored all data beyond column 72, and so people naturally
!    put stuff there (to number the lines, for instance).  This program
!    will try to simulate this fact by truncating the input line at column 72.
!    (Of course, doing this raises issues about the meaning and treatment
!    of TAB characters but I'm not the one who sold my soul to the 
!    TAB devil!)
!
!    FORTRAN77 comment lines (which begin variously with 'C', 'c',
!    'D', 'd', '*' and '!') are rewritten with an initial '!', the
!    FORTRAN90 comment initiator.
!
!    It is the author's belief that TAB characters are a blight.  Having
!    been bitten numerous times by TAB characters while using FIXCON,
!    the program has been altered to carry out the following procedure
!    in the presence of TAB characters:
!
!    * Upon detection of the first TAB character, a warning message is
!      issued.
!
!    * Every TAB character is replaced by 6 blanks, whether this was the
!      user's intention or not.
!
!    * Continuation replacement proceeds thereafter.
!
!  Usage:
!
!      fixcon file.f77 file.f90
!
!    or
!
!      fixcon file.f77
!
!    in which case the output file will automatically be called
!    "file.f90"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  character ( len = 255 ) input_file
  integer ( kind = 4 ), parameter :: input_unit = 1
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) num_arg
  character ( len = 255 ) output_file
  integer ( kind = 4 ), parameter :: output_unit = 2

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FIXCON'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a file with FORTRAN77 style continuation lines;'
  write ( *, '(a)' ) '  Write the information with FORTRAN90 continuation.'
!
!  Get the number of command line arguments.
!
!  Old style:
!
  num_arg = iargc ( )
!
!  New style:
!
! num_arg = ipxfargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input file name:'
    read ( *, '(a)', iostat = ios ) input_file

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FIXCON - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, input_file )
!
!  New style:
!
!   call pxfgetarg ( iarg, input_file, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'FIXCON - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  end if
!
!  If two command line arguments, the second one is the output file name.
!
  if ( num_arg < 2 ) then

    output_file = input_file
    call file_name_ext_swap ( output_file, 'f90' )

  else

    iarg = 2
!
!  Old style:
!
    call getarg ( iarg, output_file )
!
!  New style:
!
!   call pxfgetarg ( iarg, output_file, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'FIXCON - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  end if
!
!  Now we know what to do.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FIXCON'
  write ( *, '(a)' ) '  Read F77 file:  "' // trim ( input_file ) // '".'
  write ( *, '(a)' ) '  Write F90 file: "' // trim ( output_file ) // '".'

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FIXCON - Fatal error!'
    write ( *, '(a)' ) '  There was an error while trying to open the'
    write ( *, '(a)' ) '  input file "' // trim ( input_file ) // '".'
    stop
  end if

  open ( unit = output_unit, file = output_file, status = 'replace' )

  call fix_continuation ( input_unit, output_unit )

  close ( unit = input_unit )
  close ( unit = output_unit )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FIXCON'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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

  if ( i == 0 ) then

    if ( len_name + 1 > len_max ) then
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
!    Arthur      0  0
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
!
!    If no period occurs in FILE_NAME, then
!      I = J = 0;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s_index_last

  i = s_index_last ( file_name, '.' )

  if ( i /= 0 ) then

    j = len_trim ( file_name )

  else

    j = 0

  end if

  return
end
subroutine fix_continuation ( input, output )

!*****************************************************************************80
!
!! FIX_CONTINUATION copies a FORTRAN77 file, using FORTRAN90 continuation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, OUTPUT, the I/O units associated 
!    with the input (F77) and output (F90) files respectively.
!
  implicit none

  logical, save :: found_tab = .false.
  logical have_line
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 255 ) new_line
  character ( len = 255 ) old_line
  integer ( kind = 4 ) output
  character, parameter :: TAB = char ( 9 )

  have_line = .false.
  old_line = ' '

  do
!
!  Read the next line.
!
    read ( input, '(a)', iostat = ios ) new_line

    if ( ios /= 0 ) then
      exit
    end if
!
!  If the line contains any TAB characters, then
!
    if ( index ( new_line, TAB ) /= 0 ) then

      if ( .not. found_tab ) then
        found_tab = .true.
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FIXCON - Warning!'
        write ( *, '(a)' ) '  This file contains TAB characters.'
        write ( *, '(a)' ) '  They will be replaced by 6 blanks.'
        write ( *, '(a)' ) '  The results may not be what you want.'
      end if

      call s_tab_blanks ( new_line )

    end if
!
!  If it's an F77 comment line, write it out immediately...
!  ...but change to F90 format.
!
    if ( &
      new_line(1:1) == '!' .or. &
      new_line(1:1) == '*' .or. &
      new_line(1:1) == 'c' .or. &
      new_line(1:1) == 'C' .or. &
      new_line(1:1) == 'd' .or. &
      new_line(1:1) == 'D' ) then

      if ( have_line ) then
        write ( output, '(a)' ) trim ( old_line )
        have_line = .false.
      end if

      new_line(1:1) = '!'
      write ( output, '(a)' ) trim ( new_line )
      cycle

    end if
!
!  Now is the time to truncate the input line to 72 columns.
!
    new_line = new_line(1:72)
!
!  If there's no line in storage, cycle.
!
    if ( .not. have_line ) then
      have_line = .true.
      old_line = new_line
      cycle
    end if
!
!  If the new line is a continuation of the old line, then...
!
    if ( new_line(6:6) /= ' ' ) then
      call s_cat ( old_line, ' &', old_line )
      new_line(6:6) = ' '
    end if

    write ( output, '(a)' ) trim ( old_line )

    old_line = new_line

  end do
!
!  When we jump out of the loop, we may have had one last line unwritten.
!
  if ( have_line ) then
    write ( output, '(a)' ) trim ( old_line )
  end if

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
  character ( len = 10 ) time
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
