program main

!*****************************************************************************80
!
!! MAIN is the main program for TEC_TO_FEM.
!
!  Discussion:
!
!    TEC_TO_FEM reads a TECPLOT FEPOINT file and writes finite element files.
!
!  Usage:
!
!    tec_to_fem file.tec
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) element_file_name
  integer iarg
  integer iargc
  integer ierror
  integer ilen
  integer ios
  integer ipxfargc
  character ( len = 255 ) node_coord_file_name
  character ( len = 255 ) node_data_file_name
  integer num_arg
  character ( len = 255 ) tec_file_name

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_TO_FEM'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a TECPLOT FEPOINT file;'
  write ( *, '(a)' ) '  Write finite element files.'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input file name:'
    read ( *, '(a)', iostat = ios ) tec_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEC_TO_FEM - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, tec_file_name )

  end if
!
!  Create the output file names from the input file name.
!
  node_coord_file_name = tec_file_name
  call file_name_ext_swap ( node_coord_file_name, 'nodes.txt' )

  element_file_name = tec_file_name
  call file_name_ext_swap ( element_file_name, 'elements.txt' )

  node_data_file_name = tec_file_name
  call file_name_ext_swap ( node_data_file_name, 'values.txt' )
!
!  Now we know what to do.
!
  call tec_to_fem_handle ( tec_file_name, node_coord_file_name, &
    element_file_name, node_data_file_name )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_TO_FEM'
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
!    CH_EQI ( 'A', 'a' ) is .TRUE.
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

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

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
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

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

  write ( output_unit, '(a)'       ) '#  ' // trim ( output_file_name )
  write ( output_unit, '(a)'       ) '#  created by DTABLE_HEADER_WRITE.F90'
  write ( output_unit, '(a)'       ) '#'
  write ( output_unit, '(a,i8)'    ) '#  Spatial dimension M = ', m
  write ( output_unit, '(a,i8)'    ) '#  Number of points N = ', n
  write ( output_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
    epsilon ( x )
  write ( output_unit, '(a)'       ) '#'

  return
end
subroutine fem_write ( node_coord_file_name, element_file_name, &
  node_data_file_name, dim_num, node_num, element_num, element_order, &
  node_data_num, node_coord, element_node, node_data )

!*****************************************************************************80
!
!! FEM_WRITE writes data files associated with a finite element solution.
!
!  Discussion:
!
!    This program writes the node, element and data files that define
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
!    01 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_COORD_FILE_NAME, the name of the node
!    coordinate file.  If this argument is empty, no node coordinate file will
!    be written.
!
!    Input, character ( len = * ) ELEMENT_FILE_NAME, the name of the element
!    file.  If this argument is empty, no element file will be written.
!
!    Input, character ( len = * ) NODE_DATA_FILE_NAME, the name of the node
!    data file.  If this argument is empty, no node data file will be written.
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
!    Input, real ( kind = 8 ) NODE_DATA(NODE_DATA_NUM,NODE_NUM), the data
!    values associated with each node.
!
  implicit none

  integer dim_num
  integer element_num
  integer element_order
  integer node_data_num
  integer node_num

  character ( len = * ) element_file_name
  integer element_file_status
  integer element_file_unit
  integer element_node(element_order,element_num)
  real ( kind = 8 ) node_coord(dim_num,node_num)
  character ( len = * ) node_coord_file_name
  integer node_coord_file_status
  integer node_coord_file_unit
  real ( kind = 8 ) node_data(node_data_num,node_num)
  character ( len = * ) node_data_file_name
  integer node_data_file_status
  integer node_data_file_unit
!
!  Write the node coordinate file.
!
  if ( 0 < len_trim ( node_coord_file_name ) ) then

    call get_unit ( node_coord_file_unit )

    open ( unit = node_coord_file_unit, file = node_coord_file_name, &
      status = 'replace', iostat = node_coord_file_status )

    if ( node_coord_file_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM_WRITE - Error!'
      write ( *, '(a)' ) '  Could not open the node coordinate file.'
      stop
    end if

    call dtable_header_write ( node_coord_file_name, node_coord_file_unit, &
      dim_num, node_num )

    call dtable_data_write ( node_coord_file_unit, dim_num, node_num, &
      node_coord );

    close ( unit = node_coord_file_unit )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM_WRITE wrote node coordinates to "' &
      // trim ( node_coord_file_name ) // '".'

  end if
!
!  Write the element file.
!
  if ( 0 < len_trim ( element_file_name ) ) then

    call get_unit ( element_file_unit )

    open ( unit = element_file_unit, file = element_file_name, &
      status = 'replace', iostat = element_file_status )

    if ( element_file_status /= 0 )  then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM_WRITE - Error!'
      write ( *, '(a)' ) '  Could not open the element file.'
      stop
    end if

    call itable_header_write ( element_file_name, element_file_unit, &
      element_order, element_num )

    call itable_data_write ( element_file_unit, &
      element_order, element_num, element_node )

    close ( unit = element_file_unit )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM_WRITE wrote element data to "' // &
      trim ( element_file_name ) // '".'

  end if
!
!  Write the node data file.
!
  if ( 0 < len_trim ( node_data_file_name ) ) then

    call get_unit ( node_data_file_unit )

    open ( unit = node_data_file_unit, file = node_data_file_name, &
      status = 'replace', iostat = node_data_file_status )

    if ( node_data_file_status /= 0 )  then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM_WRITE - Error!'
      write ( *, '(a)' ) '  Could not open the node data file.'
      stop
    end if

    call dtable_header_write ( node_data_file_name, node_data_file_unit, &
      node_data_num, node_num )

    call dtable_data_write ( node_data_file_unit, &
      node_data_num, node_num, node_data )

    close ( unit = node_data_file_unit )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM_WRITE wrote node data to "' &
      // trim ( node_data_file_name ) // '".'

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
!    Output, integer I, J, the indices of the first and last characters
!    in the file extension.
!
!    If no period occurs in FILE_NAME, then
!      I = J = 0;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) file_name
  integer i
  integer j
  integer s_index_last

  i = s_index_last ( file_name, '.' )

  if ( i /= 0 ) then

    j = len_trim ( file_name )

  else

    j = 0

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
  integer i
  integer j
  integer len_max
  integer len_name

  len_max = len ( file_name )
  len_name = len_trim ( file_name )

  call file_name_ext_get ( file_name, i, j )

  if ( i == 0 ) then

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
  write ( output_unit, '(a)'       ) '#  created by ITABLE_HEADER_WRITE.F90'
  write ( output_unit, '(a)'       ) '#  at ' // trim ( string )
  write ( output_unit, '(a)'       ) '#'
  write ( output_unit, '(a,i8)'    ) '#  Spatial dimension M = ', m
  write ( output_unit, '(a,i8)'    ) '#  Number of points N = ', n
  write ( output_unit, '(a)'       ) '#'

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
  integer get
  integer put
  integer nchar
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
!    Output, integer S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none

  integer i
  integer j
  integer llen1
  integer llen2
  character ( len = * ) s
  integer s_index_last
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
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
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
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( ihave > 1 ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( ihave >= 6 .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. lchar+1 >= nchar ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

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
subroutine tec_to_fem_handle ( tec_file_name, node_coord_file_name, &
  element_file_name, node_data_file_name )

!*****************************************************************************80
!
!! TEC_TO_FEM_HANDLE reads data from a TECPLOT file and writes FEM files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the input TEC file name.
!
!    Input, character ( len = * ) NODE_COORD_FILE_NAME, the output
!    node coordinate file name.
!
!    Input, character ( len = * ) ELEMENT_FILE_NAME, the input
!    element file name.
!
!    Input, character ( len = * ) NODE_DATA_FILE_NAME, the input
!    node data file name.
!
  implicit none

  integer dim_num
  character ( len = * ) element_file_name
  integer, allocatable, dimension ( :, : ) :: element_node
  integer element_num
  integer element_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_coord
  character ( len = * ) node_coord_file_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_data
  character ( len = * ) node_data_file_name
  integer node_data_num
  integer node_num
  character ( len = * ) tec_file_name
  integer tec_file_status
  integer tec_file_unit

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading TEC file "' // trim ( tec_file_name ) // '".'

  call get_unit ( tec_file_unit )

  open ( unit = tec_file_unit, file = tec_file_name, status = 'old', &
    iostat = tec_file_status )

  if ( tec_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEC_TO_FEM'
    write ( *, '(a)' ) '  Unable to open TEC file "' &
      // trim ( tec_file_name ) // '".'
    stop
  end if

  call tec_header_read ( tec_file_name, tec_file_unit, dim_num, node_num, &
    element_num, element_order, node_data_num )

  call tec_header_print ( dim_num, node_num, element_num, &
    element_order, node_data_num )
!
!  Allocate space for the data, and read the data.
!
  allocate ( node_coord(1:dim_num,1:node_num) )
  allocate ( node_data(1:node_data_num,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )

  call tec_data_read ( tec_file_name, tec_file_unit, dim_num, &
    node_num, element_num, element_order, node_data_num, node_coord, &
    element_node, node_data )

  close ( unit = tec_file_unit )
!
!  Write the FEM files.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing FEM files:"'
  write ( *, '(a)' ) '    "' // trim ( node_coord_file_name ) // '".'
  write ( *, '(a)' ) '    "' // trim ( element_file_name ) // '".'
  write ( *, '(a)' ) '    "' // trim ( node_data_file_name ) // '".'

  call fem_write ( node_coord_file_name, element_file_name, &
    node_data_file_name, dim_num, node_num, element_num, element_order, &
    node_data_num, node_coord, element_node, node_data )

  deallocate ( node_coord )
  deallocate ( node_data )
  deallocate ( element_node )

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
!    STRING = '31 May 2001   9:45:54.872 AM'
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
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8 ) ampm
  integer d
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
  integer values(8)
  integer y

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

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
