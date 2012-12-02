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
  integer itemp

  itemp = ichar ( c )
 
  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if
 
  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.  
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
function ch_is_alpha ( c )

!*****************************************************************************80
!
!! CH_IS_ALPHA is TRUE if C is an alphabetic character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, a character to check.
!
!    Output, logical CH_IS_ALPHA is TRUE if C is an alphabetic character.
!
  implicit none

  character c
  logical ch_is_alpha

  if ( ( lle ( 'a', c ) .and. lle ( c, 'z' ) ) .or. &
       ( lle ( 'A', c ) .and. lle ( c, 'Z' ) ) ) then
    ch_is_alpha = .true.
  else
    ch_is_alpha = .false.
  end if

  return
end
subroutine dtable_data_write ( output_unit, m, n, table )

!*****************************************************************************80
!
!! DTABLE_DATA_WRITE writes data to a double precision table file.
!
!  Discussion:
!
!    This routine writes a single line of output for each point,
!    containing its spatial coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the output unit.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer m
  integer n

  integer output_unit
  integer j
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Create the format string.
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
  call s_blank_delete ( string )

  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do

  return
end
subroutine dtable_header_write ( output_file_name, output_unit, m, n )

!*****************************************************************************80
!
!! DTABLE_HEADER_WRITE writes the header to a double precision table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer OUTPUT_UNIT, the output unit.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
  implicit none

  integer m
  integer n
  character ( len = * ) output_file_name
  integer output_unit
  character ( len = 40 ) string
  real ( kind = 8 ), parameter :: x = 1.0D+00

  call timestring ( string )

  write ( output_unit, '(a)'       ) '#  ' // trim ( output_file_name )
  write ( output_unit, '(a)'       ) '#  created by TABLE_IO.F90'
  write ( output_unit, '(a)'       ) '#  at ' // trim ( string )
  write ( output_unit, '(a)'       ) '#'
  write ( output_unit, '(a,i8)'    ) '#  Spatial dimension M = ', m
  write ( output_unit, '(a,i8)'    ) '#  Number of points N = ', n
  write ( output_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
    epsilon ( x )
  write ( output_unit, '(a)'       ) '#'

  return
end
subroutine dtable_write ( output_file_name, m, n, table, header )

!*****************************************************************************80
!
!! DTABLE_WRITE writes a double precision table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
!    Input, logical HEADER, is TRUE if the header is to be included.
!
  implicit none

  integer m
  integer n

  logical header
  integer ios
  character ( len = * ) output_file_name
  integer output_unit
  real ( kind = 8 ) table(m,n)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_file_name ) // '" on unit ', output_unit
    stop
  end if

  if ( header ) then
    call dtable_header_write ( output_file_name, output_unit, m, n )
  end if

  call dtable_data_write ( output_unit, m, n, table )

  close ( unit = output_unit )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a9to99.txt'     'a0to00.txt'
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer change
  integer digit
  character ( len = * ) file_name
  integer i
  integer lens

  lens = len_trim ( file_name )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_INC - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  if ( change == 0 ) then
    file_name = ' '
    return
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

  integer i
  integer ios
  integer iunit
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
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)') i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a8)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine itable_data_write ( output_unit, m, n, table )

!*****************************************************************************80
!
!! ITABLE_DATA_WRITE writes data to an integer table file.
!
!  Discussion:
!
!    This routine writes a single line of output for each point,
!    containing its spatial coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the output unit.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer TABLE(M,N), the table data.
!
  implicit none

  integer m
  integer n

  integer output_unit
  integer j
  character ( len = 30 ) string
  integer table(m,n)
!
!  Create the format string.
!
  write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
  call s_blank_delete ( string )

  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do

  return
end
subroutine itable_header_write ( output_file_name, output_unit, m, n )

!*****************************************************************************80
!
!! ITABLE_HEADER_WRITE writes the header to an integer table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer OUTPUT_UNIT, the output unit.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
  implicit none

  integer m
  integer n
  character ( len = * ) output_file_name
  integer output_unit
  character ( len = 40 ) string

  call timestring ( string )

  write ( output_unit, '(a)'       ) '#  ' // trim ( output_file_name )
  write ( output_unit, '(a)'       ) '#  created by TABLE_IO.F90'
  write ( output_unit, '(a)'       ) '#  at ' // trim ( string )
  write ( output_unit, '(a)'       ) '#'
  write ( output_unit, '(a,i8)'    ) '#  Spatial dimension M = ', m
  write ( output_unit, '(a,i8)'    ) '#  Number of points N = ', n
  write ( output_unit, '(a)'       ) '#'

  return
end
subroutine itable_write ( output_file_name, m, n, table, header )

!*****************************************************************************80
!
!! ITABLE_WRITE writes an integer table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer TABLE(M,N), the table data.
!
!    Input, logical HEADER, is TRUE if the header is to be included.
!
  implicit none

  integer m
  integer n

  logical header
  integer ios
  character ( len = * ) output_file_name
  integer output_unit
  integer table(m,n)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_file_name ) // '" on unit ', output_unit
    stop
  end if

  if ( header ) then
    call itable_header_write ( output_file_name, output_unit, m, n )
  end if

  call itable_data_write ( output_unit, m, n, table )

  close ( unit = output_unit )

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer m
  integer n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)') i
    end do

    write ( *, '(''        Row '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '        Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(2x,i8,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine s_alpha_last ( s, iloc )

!*****************************************************************************80
!
!! S_ALPHA_LAST returns the location of the last alphabetic character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Output, integer ILOC, the location of the last alphabetic
!    character in the string.  If there are no alphabetic
!    characters, ILOC is returned as 0.
!
  implicit none

  logical ch_is_alpha
  integer i
  integer iloc
  character ( len = * ) s

  do i = len ( s ), 1, -1
    if ( ch_is_alpha ( s(i:i) ) ) then
      iloc = i
      return
    end if
  end do
 
  iloc = 0
 
  return
end
function s_begin ( s1, s2 )

!*****************************************************************************80
!
!! S_BEGIN is TRUE if one string matches the beginning of the other.
!
!  Discussion:
!
!    The strings are compared, ignoring blanks, spaces and capitalization.
!
!  Example:
!
!     S1              S2      S_BEGIN
!
!    'Bob'          'BOB'     TRUE
!    '  B  o b '    ' bo b'   TRUE
!    'Bob'          'Bobby'   TRUE
!    'Bobo'         'Bobb'    FALSE
!    ' '            'Bob'     FALSE    (Do not allow a blank to match
!                                       anything but another blank string.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to be compared.
!
!    Output, logical S_BEGIN, is TRUE if the strings match up to
!    the end of the shorter string, ignoring case.
!
  implicit none

  logical ch_eqi
  integer i1
  integer i2
  integer len1
  integer len2
  logical s_begin
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len_trim ( s1 )
  len2 = len_trim ( s2 )
!
!  If either string is blank, then both must be blank to match.
!  Otherwise, a blank string matches anything, which is not 
!  what most people want.
!
  if ( len1 == 0 .or. len2 == 0 ) then

    if ( len1 == 0 .and. len2 == 0 ) then
      s_begin = .true.
    else
      s_begin = .false.
    end if

    return

  end if

  i1 = 0
  i2 = 0
!
!  Find the next nonblank in S1.
!
  do

    do

      i1 = i1 + 1

      if ( len1 < i1 ) then
        s_begin = .true.
        return
      end if

      if ( s1(i1:i1) /= ' ' ) then
        exit
      end if

    end do
!
!  Find the next nonblank in S2.
!
    do

      i2 = i2 + 1
  
      if ( len2 < i2 ) then
        s_begin = .true.
        return
      end if

      if ( s2(i2:i2) /= ' ' ) then
        exit
      end if

    end do
!
!  If the characters match, get the next pair.
!
    if ( .not. ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
      exit
    end if

  end do

  s_begin = .false.

  return
end
subroutine s_behead_substring ( s, sub )

!*****************************************************************************80
!
!! S_BEHEAD_SUBSTRING "beheads" a string, removing a given substring.
!
!  Discussion:
!
!    Initial blanks in the string are removed first.
!
!    Then, if the initial part of the string matches the substring,
!    that part is removed and the remainder shifted left.
!
!    Initial blanks in the substring are NOT ignored.
!
!    Capitalization is ignored.
!
!    If the substring is equal to the string, then the resultant
!    string is returned as a single blank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
!    Input, character ( len = * ) SUB, the substring to be removed from
!    the beginning of the string.
!
  implicit none

  character ( len = * ) s
  logical s_eqi
  integer s_len
  character ( len = * ) sub
  integer sub_len
!
!  Remove leading blanks from the string.
!
  s = adjustl ( s )
!
!  Get lengths.
!
  s_len = len_trim ( s )
  sub_len = len_trim ( sub )

  if ( s_len < sub_len ) then
    return
  end if
!
!  If the string begins with the substring, chop it off.
!
  if ( s_eqi ( s(1:sub_len), sub(1:sub_len) ) ) then

    if ( sub_len < s_len ) then
      s = s(sub_len+1:s_len)
      s = adjustl ( s )
    else
      s = ' '
    end if

  end if

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.  
!
!  Examples:
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
  integer i
  integer len1
  integer len2
  integer lenc
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
subroutine s_inc ( s, ierror )

!*****************************************************************************80
!
!! S_INC "increments" a string.  
!
!  Discussion:
!
!    The routine tries to produce the next string, in dictionary order, 
!    following the input value of a string.  Digits, spaces, and other 
!    nonalphabetic characters are ignored.  Case is respected; in other 
!    words, the case of every alphabetic character on input will be the 
!    same on output.
!
!    The following error conditions can occur:
!
!      There are no alphabetic characters in the string.  No
!      incrementing is possible.
!
!      All alphabetic characters are equal to 'Z' or 'z'.  In this
!      case, an error value is returned, but the string is also "wrapped
!      around" so that all alphabetic characters are "A" or "a".
!
!    If the word "Tax" were input, the successive outputs would be 
!    "Tay", "Taz", "Tba", "Tbb", ...  If the input word "January 4, 1989" 
!    were input, the output would be "Januarz 4, 1989".
!
!    This routine could be useful when trying to create a unique file
!    name or variable name at run time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string whose
!    alphabetic successor is desired.  On output, if IERROR = 0,
!    S has been replaced by its successor.  If IERROR = 2,
!    S will be "wrapped around" so that all alphabetic
!    characters equal "A" or "a".
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, no alphabetic characters occur in the string.
!    2, all alphabetic characters are "Z" or "z".  S is wrapped around so 
!       that all alphabetic characters are "A" or "a".
!
  implicit none

  integer ierror
  integer ihi
  integer ilo
  integer iloc
  character ( len = * ) s

  ierror = 0
  ilo = 1
  ihi = len ( s )
!
!  Find the last alphabetic character in the string.
!
  do
 
    call s_alpha_last ( s(ilo:ihi), iloc )
!
!  If there is no alphabetic character, we can't help.
!
    if ( iloc == 0 ) then
      ierror = 1
      exit
    end if
 
    if ( s(iloc:iloc) == char ( 122 ) ) then

      s(iloc:iloc) = char ( 97 )
      ihi = iloc - 1

      if ( ihi <= 0 ) then
        ierror = 2
        exit
      end if

    else if ( s(iloc:iloc) == char ( 90 ) ) then

      s(iloc:iloc) = char ( 65 )
      ihi = iloc - 1

      if ( ihi <= 0 ) then
        ierror = 2
        exit
      end if

    else

      s(iloc:iloc) = char ( ichar ( s(iloc:iloc) ) + 1 )
      exit

    end if
 
  end do

  return
end
subroutine s_replace_ch ( s, c1, c2 )

!*****************************************************************************80
!
!! S_REPLACE_CH replaces all occurrences of one character by another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.
!
!    Input, character C1, C2, the character to be replaced, and the
!    replacement character.
!
  implicit none

  character c1
  character c2
  integer i
  character ( len = * ) s

  do i = 1, len ( s )
    if ( s(i:i) == c1 ) then
      s(i:i) = c2
    end if
  end do

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

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
!    13 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LENGTH, the number of characters of S used to make 
!    the integer.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer length
  character ( len = * ) s
  integer state
  integer value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + ichar ( c ) - ichar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

  end if

  return
end
subroutine s_word_count ( s, word_num )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer WORD_NUM, the number of "words" in the string.  
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer i
  character ( len = * ) s
  integer s_len
  integer word_num

  word_num = 0
  s_len = len ( s )

  if ( s_len <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, s_len

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      word_num = word_num + 1
      blank = .false. 
    end if

  end do

  return
end
subroutine s_word_extract ( s, w )

!*****************************************************************************80
!
!! S_WORD_EXTRACT extracts the next word from a string.
!
!  Discussion:
!
!    A "word" is a string of characters terminated by a blank or
!    the end of the string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.  On output, the first
!    word has been removed, and the remaining string has been shifted left.
!
!    Output, character ( len = * ) W, the leading word of the string.
!
  implicit none

  integer get1
  integer get2
  character ( len = * ) s
  integer s_len
  character ( len = * ) w

  w = ' '

  s_len = len_trim ( s )

  if ( s_len < 1 ) then
    return
  end if
!
!  Find the first nonblank.
!
  get1 = 0

  do

    get1 = get1 + 1

    if ( s_len < get1 ) then
      return
    end if

    if ( s(get1:get1) /= ' ' ) then
      exit
    end if

  end do
!
!  Look for the last contiguous nonblank.
!
  get2 = get1

  do

    if ( s_len <= get2 ) then
      exit
    end if

    if ( s(get2+1:get2+1) == ' ' ) then
      exit
    end if

    get2 = get2 + 1

  end do
!
!  Copy the word.
!
  w = s(get1:get2)
!
!  Shift the string.
!
  s(1:get2) = ' '
  s = adjustl ( s(get2+1:) )

  return
end
subroutine tec_data_read ( tec_file_name, tec_file_unit, dim_num, &
  node_num, element_num, element_order, node_data_num, node_coord, &
  element_node, node_data )

!*****************************************************************************80
!
!! TEC_DATA_READ reads the data from a TEC file.
!
!  Discussion:
!
!    This routine assumes that the TEC file has already been opened,
!    and that the optional TITLE record, VARIABLES record and ZONE
!    record have been read, so that the file is positioned at the
!    next record (the first data record).
!
!    After this call, the user may close the file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the name of the file.
!
!    Input, integer TEC_FILE_UNIT, the unit associated with the file.
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer NODE_DATA_NUM, the number of data items per node.
!
!    Output, real ( kind = 8 ) NODE_COORD(DIM_NUM,NODE_NUM), the coordinates 
!    of nodes.
!
!    Output, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM); 
!    the global index of local node I in element J.
!
!    Output, real ( kind = 8 ) NODE_DATA(NODE_DATA_NUM,NODE_NUM), the 
!    data values associated with each node.
!
  implicit none

  integer dim_num
  integer element_num
  integer element_order
  integer node_data_num
  integer node_num

  integer element
  integer element_node(element_order,element_num)
  integer node
  real ( kind = 8 ) node_coord(dim_num,node_num)
  real ( kind = 8 ) node_data(node_data_num,node_num)
  character ( len = * ) tec_file_name
  integer tec_file_unit
!
!  Read the node coordinates and node data.
!
  do node = 1, node_num
    read ( tec_file_unit, * ) &
      node_coord(1:dim_num,node), node_data(1:node_data_num,node)
  end do
!
!  Read the element-node connectivity.
!
  do element = 1, element_num
    read ( tec_file_unit, * ) element_node(1:element_order,element)
  end do

  return
end
subroutine tec_data_write ( tec_file_name, tec_file_unit, dim_num, &
  node_num, element_num, element_order, node_data_num, node_coord, &
  element_node, node_data )

!*****************************************************************************80
!
!! TEC_DATA_WRITE writes the data to a TEC file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the name of the file.
!
!    Input, integer TEC_FILE_UNIT, the unit associated with the file.
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer NODE_DATA_NUM, the number of data items per node.
!
!    Input, real ( kind = 8 ) NODE_COORD(DIM_NUM,NODE_NUM), the coordinates 
!    of nodes.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM); 
!    the global index of local node I in element J.
!
!    Input, real ( kind = 8 ) NODE_DATA(NODE_DATA_NUM,NODE_NUM), the 
!    data values associated with each node.
!
  implicit none

  integer dim_num
  integer element_num
  integer element_order
  integer node_data_num
  integer node_num

  integer dim
  integer element
  integer element_node(element_order,element_num)
  character ( len = 40 ) format_string
  integer node
  real ( kind = 8 ) node_coord(dim_num,node_num)
  real ( kind = 8 ) node_data(node_data_num,node_num)
  character ( len = * ) tec_file_name
  integer tec_file_unit
!
!  Write the node coordinates and node data.
!
  write ( format_string, '(a,i2,a)' ) &
    '(', dim_num + node_data_num, '(2x,g14.6))'

  do node = 1, node_num
    write ( tec_file_unit, format_string ) &
      node_coord(1:dim_num,node), node_data(1:node_data_num,node)
  end do
!
!  Write the element-node connectivity.
!
  write ( format_string, '(a,i2,a)' ) '(', element_order, '(2x,i8))'

  do element = 1, element_num
    write ( tec_file_unit, format_string ) element_node(1:element_order,element)
  end do

  return
end
subroutine tec_header_print ( dim_num, node_num, element_num, &
  element_order, node_data_num )

!*****************************************************************************80
!
!! TEC_HEADER_PRINT prints the header to a TEC file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer NODE_DATA_NUM, the number of data items per node.
!
  implicit none

  integer dim_num
  integer element_num
  integer element_order
  integer node_data_num
  integer node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension         = ', dim_num
  write ( *, '(a,i8)' ) '  Number of nodes           = ', node_num
  write ( *, '(a,i8)' ) '  Number of elements        = ', element_num
  write ( *, '(a,i8)' ) '  Element order             = ', element_order
  write ( *, '(a,i8)' ) '  Number of node data items = ', node_data_num

  return
end
subroutine tec_header_read ( tec_file_name, tec_file_unit, dim_num, node_num, &
  element_num, element_order, node_data_num )

!*****************************************************************************80
!
!! TEC_HEADER_READ reads the header from a TEC file.
!
!  Discussion:
!
!    This routine assumes that the TEC file has already been opened on
!    unit TEC_FILE_UNIT, and that it contains finite element data.
!
!    The routine reads the optional TITLE record, the VARIABLES line
!    and the ZONE line.  It leaves the file open, positioned at the next
!    record, which begins the data section.  The user may either close
!    the file, or call TEC_DATA_READ to read the data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character TEC_FILE_NAME(*), the name of the TEC file.
!
!    Input, integer TEC_FILE_UNIT, the unit number associated with the TEC file.
!
!    Output, integer DIM_NUM, the spatial dimension, inferred from the
!    names of the variables.
!
!    Output, integer NODE_NUM, the number of nodes, determined by the 
!    "N=" argument.
!
!    Output, integer ELEMENT_NUM, the number of elements, inferred from the
!    "E=" argument.
!
!    Output, integer ELEMENT_ORDER, the order of the elements, inferred from
!    the "ZONETYPE=" argument.
!
!    Output, integer NODE_DATA_NUM, the number of data items per node,
!    inferred from the the number of node data items, minus those which are
!    inferred to be spatial coordinates.
!
!    Output, real NODE_COORD(DIM_NUM,NODE_NUM), the coordinates of nodes.
!
  implicit none

  integer begin
  logical ch_eqi
  integer dim_num
  integer element_num
  integer element_order
  character ( len = 80 ) element_type
  character ( len = 255 ) line
  character ( len = 20 ) name
  integer name_len
  integer node_data_num
  integer node_num
  logical s_begin
  logical s_eqi
  character ( len = * ) tec_file_name
  integer tec_file_status
  integer tec_file_unit
  integer variable
  character ( len = 255 ) variable_name
  integer, allocatable, dimension ( : ) :: variable_name_length
  integer variable_num
!
!  Read and parse the TITLE line.
!  But it is optional, so you may have just read the VARIABLES line instead!
!
  line = ' '

  do while ( len_trim ( line ) <= 0 )

    read ( tec_file_unit, '(a)', iostat = tec_file_status ) line

    if ( tec_file_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEC_FILE_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading the file,'
      write ( *, '(a)' ) '  searching for TITLE line.'
      stop
    end if

  end do
!
!  Read the VARIABLES line.
!
!  Because the TITLE line is apparently optional, we may have already
!  read the VARIABLES line!
!
  if ( .not. s_begin ( line, 'VARIABLES=' ) ) then
    line = ' '
    do while ( len_trim ( line ) == 0 )
      read ( tec_file_unit, '(a)', iostat = tec_file_status ) line

      if ( tec_file_status /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEC_FILE_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Error while reading the file,'
        write ( *, '(a)' ) '  searching for VARIABLES line.'
        stop
      end if

    end do
  end if

  if ( .not. s_begin ( line, 'VARIABLES=' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEC_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  The VARIABLES = line is missing in the file.'
    stop
  end if
!
!  Parse the VARIABLES line.
!  VARIABLES = name1 name2 name3...
!  The names may be quoted, and are separated by quotes, commas or spaces.
!
!  Remove the initial "VARIABLES="
!
  call s_behead_substring ( line, 'VARIABLES' )
  call s_behead_substring ( line, '=' )
!
!  Replace single quotes, double quotes, commas and periods by blanks.
!
  call s_replace_ch ( line, '''', ' ' )
  call s_replace_ch ( line, '"', ' ' )
  call s_replace_ch ( line, ',', ' ' )
  call s_replace_ch ( line, '.', ' ' )
!
!  Count the words.
!
  call s_word_count ( line, variable_num )

  allocate ( variable_name_length(variable_num) )
!
!  Extract the words.
!
  begin = 0

  do variable = 1, variable_num
    call s_word_extract ( line, name )
    name_len = len_trim ( name )
    variable_name_length(variable) = name_len
    variable_name(begin+1:begin+name_len) = name(1:name_len)
    begin = begin + name_len
  end do
!
!  Based on the variable names, determine the spatial dimension and the number
!  of node data items.
!
!  For now, we SIMPLY ASSUME that the spatial coordinates are listed first.
!  Hence, when we read the node data, we assume that the first DIM_NUM values
!  represent X, Y and possibly Z.
!
  dim_num = 0
  node_data_num = variable_num

  begin = 0

  do variable = 1, variable_num

    if ( variable_name_length(variable) == 1 ) then
      name = variable_name(begin+1:begin+1)
      if ( ch_eqi ( name, 'X' ) .or. & 
          ch_eqi ( name, 'Y' ) .or. &
          ch_eqi ( name, 'Z' ) ) then
        dim_num = dim_num + 1
        node_data_num = node_data_num - 1
      end if
    end if

    begin = begin + variable_name_length(variable)

  end do
!
!  Read and parse the ZONE line.
!
  line = ' '
  do while ( len_trim ( line ) == 0 )
    read ( tec_file_unit, '(a)' ) line
  end do
    
  if ( .not. s_begin ( line, 'ZONE' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEC_READ - Fatal error!'
    write ( *, '(a)' ) '  The ZONE = line is missing in the file.'
    stop
  end if

  call tec_zone_line_parse ( line, node_num, element_num, element_type )
!
!  Based on ELEMENT_TYPE, determine the element order.
!
  if ( s_eqi ( element_type, 'FETRIANGLE' ) ) then
    element_order = 3
  elseif ( s_eqi ( element_type, 'FEQUADRILATERAL' ) ) then
    element_order = 4
  elseif ( s_eqi ( element_type, 'FETETRAHEDRON' ) ) then
    element_order = 4
  elseif ( s_eqi ( element_type, 'FEBRICK' ) ) then
    element_order = 8
  else
    element_order = -1
  end if

  deallocate ( variable_name_length )

 return
end
subroutine tec_header_write ( tec_file_name, tec_file_unit, dim_num, &
  node_num, element_num, element_order, node_data_num )

!*****************************************************************************80
!
!! TEC_HEADER_WRITE writes the header to a TEC file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the name of the file.
!
!    Input, integer TEC_FILE_UNIT, the unit associated with the file.
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer NODE_DATA_NUM, the number of data items per node.
!
  implicit none

  integer dim
  integer dim_num
  integer element
  integer element_num
  integer element_order
  integer ierror
  character ( len = 8 ) name
  integer name_num
  character ( len = 80 ) name_string
  integer node
  integer node_data_num
  integer node_num
  character ( len = * ) tec_file_name
  integer tec_file_unit
  character ( len = 15 ) zonetype
!
!  Write the title.
!
  write ( tec_file_unit, '(a)' ) 'TITLE = "' // trim ( tec_file_name ) // '"'
!
!  Write the variable names.
!
  name_string = 'VARIABLES = "'

  name = 'X'
  name_num = 0
  do dim = 1, dim_num
    name_num = name_num + 1
    if ( 1 < name_num ) then
      name_string = trim ( name_string ) // '", "'
    end if
    name_string = trim ( name_string ) // trim ( name )
    call s_inc ( name, ierror )
  end do

  name = 'data_001'
  do dim = 1, node_data_num
    name_num = name_num + 1
    if ( 1 < name_num ) then
      name_string = trim ( name_string ) // '", "'
    end if
    name_string = trim ( name_string ) // trim ( name )
    call file_name_inc ( name )
  end do

  write ( tec_file_unit, '(a)' ) trim ( name_string )
!
!  Write the ZONE record.
!
  if ( dim_num == 2 .and. element_order == 3 ) then
    zonetype = 'FETRIANGLE'
  elseif ( dim_num == 2 .and. element_order == 4 ) then
    zonetype = 'FEQUADRILATERAL'
  elseif ( dim_num == 3 .and. element_order == 4 ) then
    zonetype = 'FETETRAHEDRON'
  elseif ( dim_num == 3 .and. element_order == 8 ) then
    zonetype = 'FEBRICK'
  else
    zonetype = 'FEUNKNOWN'
  end if

  write ( tec_file_unit, '(a,i8,a,i8,a)' ) 'ZONE  N = ', node_num, ',  E = ', &
    element_num, ',  DATAPACKING = POINT,  ZONETYPE = ' // trim ( zonetype )

  return
end
subroutine tec_open_read ( tec_file_name, tec_file_unit )

!*****************************************************************************80
!
!! TEC_OPEN_READ opens a TEC file for reading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the name of the file.
!
!    Output, integer TEC_FILE_UNIT, the unit on which the file has
!    been opened.  If the file could not be opened, then TEC_FILE_UNIT
!    is returned with the value of -1.
!
  implicit none

  character ( len = * ) tec_file_name
  integer tec_file_status
  integer tec_file_unit

  call get_unit ( tec_file_unit )

  open ( unit = tec_file_unit, file = tec_file_name, status = 'old', &
    iostat = tec_file_status )

  if ( tec_file_status /= 0 ) then
    tec_file_unit = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEC_OPEN_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( tec_file_name ) // '".'
    stop
  end if

  return
end
subroutine tec_open_write ( tec_file_name, tec_file_unit )

!*****************************************************************************80
!
!! TEC_OPEN_WRITE opens a TEC file for writing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the name of the file.
!
!    Output, integer TEC_FILE_UNIT, the unit on which the file has
!    been opened.  If the file could not be opened, then TEC_FILE_UNIT
!    is returned with the value of -1.
!
  implicit none

  character ( len = * ) tec_file_name
  integer tec_file_status
  integer tec_file_unit

  call get_unit ( tec_file_unit )

  open ( unit = tec_file_unit, file = tec_file_name, status = 'replace', &
    iostat = tec_file_status )

  if ( tec_file_status /= 0 ) then
    tec_file_unit = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEC_OPEN_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( tec_file_name ) // '".'
    stop
  end if

  return
end
subroutine tec_write ( tec_file_name, dim_num, node_num, element_num, &
  element_order, node_data_num, node_coord, element_node, node_data )

!*****************************************************************************80
!
!! TEC_WRITE writes finite element data to a TEC file.
!
!  Discussion:
!
!    This routine writes the node, element and data files that define
!    a finite element geometry and data based on that geometry:
!    * a set of nodes, 
!    * a set of elements based on those nodes, 
!    * a set of data values associated with each node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the name of the file.
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer NODE_DATA_NUM, the number of data items per node.
!
!    Input, real ( kind = 8 ) NODE_COORD(DIM_NUM,NODE_NUM), the coordinates 
!    of nodes.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM); 
!    the global index of local node I in element J.
!
!    Input, real ( kind = 8 ) NODE_DATA(NODE_DATA_NUM,NODE_NUM), the 
!    data values associated with each node.
!
  implicit none

  integer dim_num
  integer element_num
  integer element_order
  integer node_data_num
  integer node_num

  integer dim
  integer element
  integer element_node(element_order,element_num)
  character ( len = 40 ) format_string
  character ( len = 8 ) name
  integer name_num
  character ( len = 80 ) name_string
  integer node
  real ( kind = 8 ) node_coord(dim_num,node_num)
  real ( kind = 8 ) node_data(node_data_num,node_num)
  character ( len = * ) tec_file_name
  integer tec_file_unit
  character ( len = 15 ) zonetype
!
!  Open the file.
!
  call tec_open_write ( tec_file_name, tec_file_unit )

  if ( tec_file_unit == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEC_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if
!
!  Write the header.
!
  call tec_header_write ( tec_file_name, tec_file_unit, dim_num, &
    node_num, element_num, element_order, node_data_num )
!
!  Write the node coordinates and node data.
!
  call tec_data_write ( tec_file_name, tec_file_unit, dim_num, &
    node_num, element_num, element_order, node_data_num, node_coord, &
    element_node, node_data )
!
!  Close the file.
!
  close ( unit = tec_file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_WRITE wrote all data to "' &
    // trim ( tec_file_name ) // '".'

  return
end
subroutine tec_zone_line_parse ( line, node_num, element_num, element_type )

!*****************************************************************************80
!
!! TEC_ZONE_LINE_PARSE parses the "ZONE" line of a TEC file.
!
!  Discussion:
!
!    The string begins with the substring "ZONE" and is followed by
!    a sequence of keywords and values joined by an equals sign.
!
!    We expect the following, but in arbitrary order, separated 
!    by spaces or commas:
!
!      N = number of nodes
!      E = number of elements
!      T = optional zone title (we can't handle this right now)
!      PACKING = POINT
!      ZONETYPE = FETRIANGLE or FEQUADRILATERAL or FETETRAHEDRON or FEBRICK.
!
!    Other arguments that may appear on this line will be ignore.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string of characters, representing the
!    "VARIABLES=" line of the file.
!
!    Output, integer NODE_NUM, the number of nodes.
!
!    Output, integer ELEMENT_NUM, the number of elements.
!
!    Output, character ( len = * ) ELEMENT_TYPE, the element type: 
!    FETRIANGLE or FEQUADRILATERAL or FETETRAHEDRON or FEBRICK.
!
  implicit none

  logical ch_eqi
  integer element_num
  character ( len = * ) element_type
  integer found_num
  integer ierror
  integer length
  character ( len = * ) line
  character ( len = 80 ) name
  integer node_num
  logical s_eqi
  character ( len = 80 ) value
!
!  Remove the initial "ZONE"
!
  call s_behead_substring ( line, 'ZONE' )
!
!  Replace each EQUALS sign by a space.
!  Also get rid of commas and periods.
!  Do single and double quotes have to go, also?
!
  call s_replace_ch ( line, '=', ' ' )
  call s_replace_ch ( line, ',', ' ' )
  call s_replace_ch ( line, '.', ' ' )
!
!  Now each pair of words represents a name and a value.
!
  node_num = -1
  element_num = -1
  element_type = ' '

  found_num = 0

  do

    call s_word_extract ( line, name )

    if ( len_trim ( name ) == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEC_ZONE_LINE_PARSE - Fatal error!'
      write ( *, '(a)' ) '  Unexpected End of input.'
      stop
    end if

    call s_word_extract ( line, value )

    if ( len_trim ( value ) == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEC_ZONE_LINE_PARSE - Fatal error!'
      write ( *, '(a)' ) '  Unexpected End of input.'
      stop
    end if

    if ( ch_eqi ( name(1:1), 'N' ) .and. node_num == -1 ) then

      call s_to_i4 ( value, node_num, ierror, length )
      found_num = found_num + 1

    elseif ( ch_eqi ( name(1:1), 'E' ) .and. element_num == -1 ) then

      call s_to_i4 ( value, element_num, ierror, length )
      found_num = found_num + 1

    elseif ( s_eqi ( name, 'DATAPACKING' ) ) then

      if ( .not. s_eqi ( value, 'POINT' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEC_ZONE_LINE_PARSE - Fatal error!'
        write ( *, '(a)' ) '  Value of DATAPACKING argument must be POINT.'
        stop
      end if

    elseif ( s_eqi ( name, 'ZONETYPE' ) .and. &
      len_trim ( element_type ) == 0 ) then

      found_num = found_num + 1
      element_type = value

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Ignoring "' // trim ( name ) &
        // '" = "' // trim ( value ) // '".'

    end if

    if ( found_num == 3 ) then
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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  character ( len = 10 )  time
  integer values(8)
  integer y
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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
