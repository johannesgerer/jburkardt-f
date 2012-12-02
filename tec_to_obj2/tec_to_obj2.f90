program main

!*****************************************************************************80
!
!! MAIN is the main program for TEC_TO_OBJ2.
!
!  Discussion:
!
!    TEC_TO_OBJ2 reads a TECPLOT 3D surface file and writes an OBJ file.
!
!  Usage:
!
!    tec_to_obj2 file.dat
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2006
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
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  character ( len = 255 ) obj_file_name
  integer ( kind = 4 ) num_arg
  character ( len = 255 ) tec_file_name

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_TO_OBJ2'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a TECPLOT FEPOINT file describing a 3D surface;'
  write ( *, '(a)' ) '  Write an OBJ file.'
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
      write ( *, '(a)' ) 'TEC_TO_OBJ2 - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, tec_file_name )

  end if
!
!  Create the output file name from the input file name.
!
  obj_file_name = tec_file_name
  call file_name_ext_swap ( obj_file_name, 'obj' )
!
!  Now we know what to do.
!
  call tec_to_obj_handle ( tec_file_name, obj_file_name )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_TO_OBJ2'
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine digit_to_ch ( digit, ch )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Discussion:
!
!    Instead of CHAR, we now use the ACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!    DIGIT   CH
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character CH, the corresponding character.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    ch = achar ( digit + 48 )

  else

    ch = '*'

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
!    Output, integer ( kind = 4 ) I, J, the indices of the first and last 
!    characters in the file extension.
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
subroutine i4_to_s_zero ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_ZERO converts an I4 to a string, with zero padding.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ).
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  000001
!        -1  -00001
!         0  000000
!      1952  001952
!    123456  123456
!   1234567  ******  <-- Not enough room!
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
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right justified, and zero padded.  
!    If there is not enough space, the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = -ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  Working from right to left, strip off the digits of the integer
!  and place them into S(ILO:IHI).
!
  ipos = ihi

  do while ( ival /= 0 .or. ipos == ihi )

    idig = mod ( ival, i4_10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

  end do
!
!  Fill the empties with zeroes.
!
  do i = ilo, ipos
    s(i:i) = '0'
  end do
 
  return
end
subroutine obj_header_write ( obj_file_name, obj_file_unit, tec_file_name, &
  obj_text_num )

!*****************************************************************************80
!
!! OBJ_HEADER_WRITE writes the header for an OBJ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OBJ_FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) OBJ_FILE_UNIT, the unit number of the file.
!
!    Input, character ( len = * ) TEC_FILE_NAME, the name of the input TEC file.
!
!    Output, integer ( kind = 4 ) OBJ_TEXT_NUM, the number of lines in the OBJ file.
!
  implicit none

  integer   ( kind = 4 ) obj_file_unit
  character ( len = *  ) obj_file_name
  integer   ( kind = 4 ) obj_text_num
  character ( len = *  ) tec_file_name

  write ( obj_file_unit, '(a)' ) '# "' // trim ( obj_file_name ) // '"'
  write ( obj_file_unit, '(a)' ) '# created by TEC_TO_OBJ.F90'
  write ( obj_file_unit, '(a)' ) '# from data extracted from "' &
    // trim ( tec_file_name ) // '".'

  obj_text_num = 3

  return
end
subroutine obj_data_write ( obj_file_unit, node_num, face_num, normal_num, &
  order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal, &
  element_offset, node_offset, obj_group_num, obj_text_num )

!*****************************************************************************80
!
!! OBJ_DATA_WRITE writes graphics information to an Alias OBJ file.
!
!  Discussion:
!
!    If no normal vectors are supplied (NORMAL_NUM <= 0) then
!    a simple format is used for the "F" records.  Otherwise,
!    the "v//vn" format is used.
!
!  Example:
!
!    #  no_normals.obj
!
!    g Group002
!
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!
!    f 8 9 11 10
!    f 12 13 15 14
!    ...
!    f 788 806 774
!
!    #  normals_supplied.obj
!
!    g Group001
!
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!
!    vn 0.0 1.0 0.0
!    vn 1.0 0.0 0.0
!    ...
!    vn 0.0 0.0 1.0
!
!    f 8//1 9//2 11//3 10//4
!    f 12//5 13//6 15//7 14//8
!    ...
!    f 788//800 806//803 774//807
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OBJ_FILE_UNIT, the unit number of the output file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) NORMAL_NUM, the number of normal vectors.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates of points.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) FACE_NODE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input, real ( kind = 8 ) NORMAL_VECTOR(3,NORMAL_NUM), normal vectors.
!
!    Input, integer ( kind = 4 ) VERTEX_NORMAL(ORDER_MAX,FACE_NUM), the indices of normal
!    vectors per vertex.
!
!    Input/output, integer ( kind = 4 ) OBJ_GROUP_NUM, the number of groups.
!
!    Input/output, integer ( kind = 4 ) OBJ_TEXT_NUM, the number of lines of text
!    in the OBJ file.
!
  implicit none

  integer   ( kind = 4 ) face_num
  integer   ( kind = 4 ) node_num
  integer   ( kind = 4 ) normal_num
  integer   ( kind = 4 ) order_max

  integer   ( kind = 4 ) element_offset
  integer   ( kind = 4 ) face
  integer   ( kind = 4 ) face_node(order_max,face_num)
  integer   ( kind = 4 ) face_order(face_num)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) node
  integer   ( kind = 4 ) node_offset
  real      ( kind = 8 ) node_xyz(3,node_num)
  integer   ( kind = 4 ) normal
  real      ( kind = 8 ) normal_vector(3,normal_num)
  integer   ( kind = 4 ) obj_file_unit
  integer   ( kind = 4 ) obj_file_status
  integer   ( kind = 4 ) obj_group_num
  integer   ( kind = 4 ) obj_text_num
  character ( len = 3  ) s3
  character ( len = 256 ) text
  character ( len = 256 ) text2
  integer   ( kind = 4 ) vertex
  integer   ( kind = 4 ) vertex_normal(order_max,face_num)
  real      ( kind = 8 ) w

  obj_group_num = obj_group_num + 1

  write ( obj_file_unit, '(a)' ) ' '
  call i4_to_s_zero ( obj_group_num, s3 )
  write ( obj_file_unit, '(a)' ) 'g Group' // s3

  obj_text_num = obj_text_num + 2
!
!  V: vertex coordinates.
!  For some reason, a fourth "coordinate" may be recommended.
!  What is its meaning?
!
  if ( 0 < node_num ) then
    write ( obj_file_unit, '(a)' ) ' '
    obj_text_num = obj_text_num + 1
  end if

  w = 1.0D+00
  do node = 1, node_num
    write ( text, '(a1,2x,4g14.6)' ) 'v', node_xyz(1:3,node), w
    call s_blanks_delete ( text )
    write ( obj_file_unit, '(a)' ) trim ( text )
    obj_text_num = obj_text_num + 1
  end do
!
!  VN: normal vectors.
!
  if ( 0 < normal_num ) then

    write ( obj_file_unit, '(a)' ) ' '
    obj_text_num = obj_text_num + 1

    do normal = 1, normal_num

      write ( text, '(a2,2x,3f7.3)' ) 'vn', normal_vector(1:3,normal)
      call s_blanks_delete ( text )
      write ( obj_file_unit, '(a)' ) trim ( text )
      obj_text_num = obj_text_num + 1

    end do

  end if
!
!  F: Faces, specified as a list of triples, one triple for each vertex:
!    vertex index/vertex texture index/vertex normal index
!
  if ( 0 < face_num ) then
    write ( obj_file_unit, '(a)' ) ' '
    obj_text_num = obj_text_num + 1
  end if

  do face = 1, face_num

    text = 'f'

    if ( normal_num <= 0 ) then

      do vertex = 1, face_order(face)
        text2 = ' '
        write ( text2(2:), '(i8)' ) face_node(vertex,face) + node_offset
        call s_blank_delete ( text2(2:) )
        call s_cat ( text, text2, text )
      end do

    else

      do vertex = 1, face_order(face)
        text2 = ' '
        write ( text2(2:), '(i8, ''//'', i8 )' ) face_node(vertex,face) + node_offset, &
          vertex_normal(vertex,face)
        call s_blank_delete ( text2(2:) )
        call s_cat ( text, text2, text )
      end do

    end if

    write ( obj_file_unit, '(a)' ) trim ( text )
    obj_text_num = obj_text_num + 1

  end do

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
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
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
  integer ( kind = 4 ) s_len
  character ( len = * ) sub
  integer ( kind = 4 ) sub_len
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
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    Thanks to Bill Richmond for pointing out a programming flaw which
!    meant that, as characters were slid to the left through multiple
!    blanks, their original images were not blanked out.  This problem
!    is easiest resolved by using a copy of the string.
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 September 2004
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character newchr
  character oldchr
  character ( len = * ) s
  character ( len = len ( s ) ) s_copy
  integer ( kind = 4 ) s_length
  character, parameter :: TAB = char ( 9 )

  s_length = len ( s )

  j = 0
  s_copy(1:s_length) = s(1:s_length)
  s(1:s_length) = ' '

  newchr = ' '

  do i = 1, s_length

    oldchr = newchr
    newchr = s_copy(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

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
  integer ( kind = 4 ) i
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
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
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
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real    ( kind = 8 ) r
  real    ( kind = 8 ) rbot
  real    ( kind = 8 ) rexp
  real    ( kind = 8 ) rtop
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
!    Output, integer ( kind = 4 ) WORD_NUM, the number of "words" in the string.  
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_len
  integer ( kind = 4 ) word_num

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

  integer ( kind = 4 ) get1
  integer ( kind = 4 ) get2
  character ( len = * ) s
  integer ( kind = 4 ) s_len
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
subroutine tec_zone_data_read ( tec_file_unit, line, dim_num, &
  node_num, element_num, element_order_uniform, node_data_num, node_coord, &
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
!    Input, integer ( kind = 4 ) TEC_FILE_UNIT, the unit associated with the file.
!
!    Input/output, character ( len = * ) LINE.  On input, the first line
!    of the data records.  On output, blank.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER_UNIFORM, the order of the elements.
!
!    Input, integer ( kind = 4 ) NODE_DATA_NUM, the number of data items per node.
!
!    Output, real ( kind = 8 ) NODE_COORD(DIM_NUM,NODE_NUM), the coordinates 
!    of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER_UNIFORM,ELEMENT_NUM); 
!    the global index of local node I in element J.
!
!    Output, real ( kind = 8 ) NODE_DATA(NODE_DATA_NUM,NODE_NUM), the 
!    data values associated with each node.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order_uniform
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order_uniform,element_num)
  character ( len = * ) line
  integer ( kind = 4 ) node
  real    ( kind = 8 ) node_coord(dim_num,node_num)
  real    ( kind = 8 ) node_data(node_data_num,node_num)
  integer ( kind = 4 ) tec_file_unit
!
!  Read the node coordinates and node data.
!
  do node = 1, node_num
    if ( 1 < node ) then
      read ( tec_file_unit, '(a)' ) line
    end if
    read ( line, * ) &
      node_coord(1:dim_num,node), node_data(1:node_data_num,node)

  end do
!
!  Read the element-node connectivity.
!
  do element = 1, element_num
    read ( tec_file_unit, '(a)' ) line
    read ( line, * ) element_node(1:element_order_uniform,element)
  end do

  line = ' '

  return
end
subroutine tec_header_read ( tec_file_name, title, variable_num, dim_num, &
  node_data_num, zone_num, node_total, element_total, &
  element_order_uniform )

!*****************************************************************************80
!
!! TEC_HEADER_READ reads all the header information from a TEC file.
!
!  Discussion:
!
!    This routine opens the TEC file and reads it, ignoring most of the
!    data, but counting certain items that will be useful in allocating
!    memory and processing the data.
!
!    In particular, the routine expects to see:
!
!    An optional "TITLE=" record.
!
!    A "VARIABLES=" record, which may extend over several lines, containingg
!    items of the form "NAME=VALUE"; the special names "X", "Y" and "Z" are
!    assumed to indicate spatial coordinates, and their presence or absence
!    is in indicator of the spatial dimension of the data.
!
!    One or more ZONE records.  A ZONE record consists of a ZONE header and
!    ZONE data.  A ZONE header begins with the word "ZONE",
!    may extend over several lines,  and contains items of the form 
!    "NAME=VALUE"; the special name "N" or "NODES" counts
!    the nodes in this zone; "E" or "ELEMENTS" counts the elements.
!    The ZONE header is followed by the ZONE data, which is a line of numeric 
!    data for each node.  The number of items in the line is the same as the 
!    number of variables in the VARIABLES record.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character TEC_FILE_NAME(*), the name of the TEC file.
!
!    Output, character ( len = * ) TITLE, the title.
!
!    Output, integer ( kind = 4 ) VARIABLE_NUM, the number of variables.
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension, inferred 
!    from the names of the variables.
!
!    Output, integer ( kind = 4 ) NODE_DATA_NUM, the number of data items per 
!    node, inferred from the the number of node data items, minus those which 
!    are inferred to be spatial coordinates.
!
!    Output, integer ( kind = 4 ) ZONE_NUM, the number of zones.
!
!    Output, integer ( kind = 4 ) NODE_TOTAL, the total number of nodes, over 
!    all zones.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the total number of elements,
!    over all zones.
!
!    Output, integer ( kind = 4 ) ELEMENT_ORDER_UNIFORM, the element order
!    that is presumably used by all elements.
!
  implicit none

  integer   ( kind = 4 )  begin
  logical                 ch_eqi
  integer   ( kind = 4 )  dim_num
  integer   ( kind = 4 )  element_num
  integer   ( kind = 4 )  element_order
  integer   ( kind = 4 )  element_order_uniform
  integer   ( kind = 4 )  element_total
  character ( len = 40 )  element_type
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  length
  character ( len = 255 ) line
  integer   ( kind = 4 )  line_num
  character ( len = 20 )  name
  integer   ( kind = 4 )  node_data_num
  integer   ( kind = 4 )  node_num
  integer   ( kind = 4 )  node_total
  integer   ( kind = 4 )  read_status
  logical                 s_begin
  logical                 s_eqi
  character ( len = * )   tec_file_name
  integer   ( kind = 4 )  tec_file_status
  integer   ( kind = 4 )  tec_file_unit
  character ( len = 255 ) title
  character ( len = 40  ) value
  integer   ( kind = 4 )  variable
  character ( len = 255 ) variable_name
  integer   ( kind = 4 )  variable_inc
  integer   ( kind = 4 )  variable_num
  integer   ( kind = 4 )  zone_num

  line_num = 0

  call get_unit ( tec_file_unit )

  open ( unit = tec_file_unit, file = tec_file_name, iostat = tec_file_status, &
    status = 'old' )
!
!  The TITLE = line is optional.
!
  read ( tec_file_unit, '(a)', iostat = tec_file_status ) line

  if ( tec_file_status /= 0 ) then
    return
  end if

  line_num = line_num + 1

  if ( s_begin ( line, 'TITLE=' ) ) then
    call s_behead_substring ( line, 'TITLE' )
    call s_behead_substring ( line, '=' )
    title = line
    line = ' '
  else
    title = '(No title supplied)'
  end if
!
!  Read and parse the 'VARIABLES =' lines.
!  But it is optional, so you may have just read the VARIABLES line instead!
!
  read_status = 0

  variable_num = 0
  dim_num = 0
  node_data_num = 0

  do

    read ( tec_file_unit, '(a)', iostat = tec_file_status ) line

    if ( tec_file_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEC_FILE_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading the file,'
      write ( *, '(a)' ) '  searching for TITLE line.'
      stop
    end if

    line_num = line_num + 1

    if ( read_status == 0 ) then

      if ( s_begin ( line, 'VARIABLES=' ) ) then

        read_status = 1

        call s_behead_substring ( line, 'VARIABLES' )
        call s_behead_substring ( line, '=' )
!
!  Blank lines are OK.
!
      else if ( len_trim ( line ) == 0 ) then

!
!  Anything else unacceptable at this point.
!
      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEC_FILE_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Expecting "VARIABLES =" or blank line,'
        write ( *, '(a)' ) '  Encountered "' // trim ( line ) // '" instead.'
        stop

      end if
!
!  Otherwise, if the first character is not a quote, we're done with variables.
!
    else if ( read_status == 1 ) then

      if ( s_begin ( line, '"' ) .or. &
           s_begin ( line, '''' ) .or. &
           len_trim ( line ) == 0 ) then
      else
        read_status = 2
        exit
      end if

    end if
!
!  Parse one or more variable names.
!  They are single or double quoted and may be separated by commas or spaces.
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
    call s_word_count ( line, variable_inc )
    variable_num = variable_num + variable_inc
!
!  Extract the words just to check whether they represent coordinate
!  directions or other node data.
!
    do variable = 1, variable_inc
      call s_word_extract ( line, name )
      if ( ch_eqi ( name, 'X' ) .or. & 
          ch_eqi ( name, 'Y' ) .or. &
          ch_eqi ( name, 'Z' ) ) then
        dim_num = dim_num + 1
      else
        node_data_num = node_data_num + 1
      end if
    end do

  end do
!
!  We got here because LINE contains text that was rejected by the VARIABLES
!  processor.  So it probably already contains the beginning of our ZONE record.
!
!  READ_STATUS = 0: We haven't recognized a ZONE record yet.
!  READ_STATUS = 1: We've read "ZONE", and haven't seen the end of the zone stuff.
!  READ_STATUS = 2: We're reading (and ignoring) numeric data for a zone.
!
  zone_num = 0
  node_total = 0
  element_total = 0
  element_order_uniform = 0

  read_status = 0

  do
!
!  Ignore blank lines.
!
    if ( 0 < len_trim ( line ) ) then
!
!  READ_STATUS = 0 expects a ZONE record.
!
      if ( read_status == 0 ) then

        if ( s_begin ( line, 'ZONE' ) ) then
          call s_behead_substring ( line, 'ZONE' )
          read_status = 1
          zone_num = zone_num + 1
          node_num = -1
          element_num = -1
          element_type = ' '
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEC_HEADER_READ:'
          write ( *, '(a)' ) '  Puzzled by line "' // trim ( line ) // '".'
          stop
        end if
!
!  READ_STATUS = 2 expects more numbers, or ZONE.
!
      else if ( read_status == 2 ) then

        if ( s_begin ( line, 'ZONE' ) ) then
          call s_behead_substring ( line, 'ZONE' )
          read_status = 1
          zone_num = zone_num + 1
          node_num = -1
          element_num = -1
          element_type = ' '
        end if

      end if
!
!  If we are in READ_STATUS = 1, but the line contains no equal signs,
!  we probably need to move to READ_STATUS = 2.
!
      if ( read_status == 1 ) then
        if ( 0 < len_trim ( line ) ) then
          if ( index ( line, '=' ) <= 0 ) then
            read_status = 2
          end if
        end if
      end if

      if ( read_status == 1 ) then
!
!  This part of code is only reached for READ_STATUS = 1.
!  LINE may be blank, or contain an unknown number of pairs of
!  "Name=Value" items,
!
!  We are particularly interested in node, element and element order information.
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
      do

        call s_word_extract ( line, name )

        if ( len_trim ( name ) <= 0 ) then
          exit
        end if

        call s_word_extract ( line, value )

        if ( len_trim ( value ) == 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEC_HEADER_READ - Fatal error!'
          write ( *, '(a)' ) '  Unexpected End of input line.'
          write ( *, '(a,i8)' ) '  Reading line ', line_num
          write ( *, '(a)' ) '  Pretend nothing is wrong...'
          value = '0'
        end if

        if ( ( ch_eqi ( name, 'N' ) .or. ch_eqi ( name, 'NODES' ) ) .and. &
             node_num == -1 ) then

          call s_to_i4 ( value, node_num, ierror, length )
          node_total = node_total + node_num

        elseif ( ( ch_eqi ( name, 'E' ) .or. ch_eqi ( name, 'ELEMENTS' ) ) .and. &
          element_num == -1 ) then

          call s_to_i4 ( value, element_num, ierror, length )
          element_total = element_total + element_num

        elseif ( s_eqi ( name, 'DATAPACKING' ) ) then

          if ( .not. s_eqi ( value, 'POINT' ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TEC_HEADER_READ - Fatal error!'
            write ( *, '(a)' ) '  Unacceptable DATAPACKING value.'
            write ( *, '(a)' ) '  Only "DATAPACKING = POINT" is supported.'
            stop
          end if

        elseif ( s_eqi ( name, 'ZONETYPE' ) .and. &
          len_trim ( element_type ) == 0 ) then

          element_type = value

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

          if ( element_order_uniform == 0 ) then
            element_order_uniform = element_order
          else if ( element_order_uniform == element_order ) then

          else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TEC_HEADER_READ - Fatal error!'
            write ( *, '(a)' ) '  More than one element order in file.'
            stop
          end if

        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Ignoring "' // trim ( name ) &
            // '" = "' // trim ( value ) // '".'

        end if

      end do

      end if

    end if
!
!  Read the next line.
!
    read ( tec_file_unit, '(a)', iostat = tec_file_status ) line

    if ( tec_file_status /= 0 ) then
      exit
    end if

    line_num = line_num + 1

  end do

  close ( unit = tec_file_unit )

 return
end
subroutine tec_header_skip ( tec_file_unit, line )

!*****************************************************************************80
!
!! TEC_HEADER_SKIP skips the TITLE and VARIABLES records in a TEC file.
!
!  Discussion:
!
!    This routine assumes the TEC file has just been opened, positioned
!    to read the first record.
!
!    It reads the optional "TITLE=" record, and the multiline
!    "VARIABLES=" record.  
!
!    As soon as it reads the beginning of the first ZONE record, it
!    stops reading, and returns the text of the first line of the first
!    ZONE record in LINE.
!
!    That's all.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEC_FILE_UNIT, the unit number of the Tec file.
!
!    Output, character ( len = * ) LINE, the last line of text read by this
!    routine, which presumably contains the beginning of a ZONE record.
!
  implicit none

  integer   ( kind = 4 )  begin
  logical                 ch_eqi
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  length
  character ( len = 255 ) line
  integer   ( kind = 4 )  line_num
  character ( len = 20 )  name
  integer   ( kind = 4 )  node_data_num
  integer   ( kind = 4 )  read_status
  logical                 s_begin
  logical                 s_eqi
  integer   ( kind = 4 )  tec_file_status
  integer   ( kind = 4 )  tec_file_unit
  character ( len = 40  ) value
  integer   ( kind = 4 )  variable
  character ( len = 255 ) variable_name
  integer   ( kind = 4 )  variable_inc
  integer   ( kind = 4 )  variable_num

  line_num = 0
!
!  The TITLE = line is optional.
!
  read ( tec_file_unit, '(a)', iostat = tec_file_status ) line

  if ( tec_file_status /= 0 ) then
    return
  end if

  line_num = line_num + 1

  if ( s_begin ( line, 'TITLE=' ) ) then
    call s_behead_substring ( line, 'TITLE' )
    call s_behead_substring ( line, '=' )
    line = ' '
  end if
!
!  Read and parse the 'VARIABLES =' lines.
!  But it is optional, so you may have just read the VARIABLES line instead!
!
  read_status = 0

  do

    read ( tec_file_unit, '(a)', iostat = tec_file_status ) line

    if ( tec_file_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEC_FILE_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading the file,'
      write ( *, '(a)' ) '  searching for TITLE line.'
      stop
    end if

    line_num = line_num + 1

    if ( read_status == 0 ) then

      if ( s_begin ( line, 'VARIABLES=' ) ) then

        read_status = 1

        call s_behead_substring ( line, 'VARIABLES' )
        call s_behead_substring ( line, '=' )
!
!  Blank lines are OK.
!
      else if ( len_trim ( line ) == 0 ) then
!
!  Anything else unacceptable at this point.
!
      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEC_FILE_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Expecting "VARIABLES =" or blank line,'
        write ( *, '(a)' ) '  Encountered "' // trim ( line ) // '" instead.'
        stop

      end if
!
!  Otherwise, if the first character is not a quote, we're done with variables.
!
    else if ( read_status == 1 ) then

      if ( s_begin ( line, '"' ) .or. &
           s_begin ( line, '''' ) .or. &
           len_trim ( line ) == 0 ) then
      else
        read_status = 2
        exit
      end if

    end if

  end do

 return
end
subroutine tec_to_obj_handle ( tec_file_name, obj_file_name )

!*****************************************************************************80
!
!! TEC_TO_OBJ_HANDLE reads data from a TECPLOT file and writes an OBJ file.
!
!  Discussion:
!
!    There are MANY kinds of TECPLOT file.
!
!    This routine is only intended to work in the very specific case
!    where the TECPLOT file describes a triangulated 3D surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the input TEC file name. 
!
!    Input, character ( len = * ) OBJ_FILE_NAME, the output 
!    OBJ file name. 
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DIM_NUM, the spatial dimension, inferred
!    by whether the variables include just "X", or "X" and "Y", or
!    "X", "Y" and "Z".
!
!    Local, integer ( kind = 4 ) ELEMENT_ORDER_UNIFORM, the order (number of
!    vertices or nodes) of all elements.  This value is assumed to be 
!    uniform across all zones.
!
!    Local, integer ( kind = 4 ) ELEMENT_TOTAL, the total number of elements
!    (over all zones) in the TEC file.
!
!    Local, integer ( kind = 4 ) NODE_DATA_NUM, the number of variables
!    which are not spatial coordinates.
!
!    Local, integer ( kind = 4 ) NODE_TOTAL, the total number (over all
!    zones) of nodes in the TEC file.
!
!    Local, character ( len = * ) TITLE, the TEC file title.
!
!    Local, integer ( kind = 4 ) VARIABLE_NUM, the number of variables
!    declared by the TEC file.  This includes geometric coordinates
!    (such as X, Y, Z) as well as associated data.
!
!    Local, integer ( kind = 4 ) ZONE_NUM, the number of zones in the
!    TEC file.
!
  implicit none

  integer   ( kind = 4  ) dim_num
  integer   ( kind = 4  ), allocatable, dimension ( :, : ) :: element_node
  integer   ( kind = 4  ) element_num
  integer   ( kind = 4  ) element_offset
  integer   ( kind = 4  ), allocatable, dimension ( : ) :: element_order
  integer   ( kind = 4  ) element_order_uniform
  integer   ( kind = 4  ) element_total
  character ( len = 80  ) element_type
  character ( len = 255 ) line
  real      ( kind = 8  ), allocatable, dimension ( :, : ) :: node_coord
  real      ( kind = 8  ), allocatable, dimension ( :, : ) :: node_data
  integer   ( kind = 4  ) node_data_num
  integer   ( kind = 4  ) node_num
  integer   ( kind = 4  ) node_offset
  integer   ( kind = 4  ) node_total
  integer   ( kind = 4  ) normal_num
  real      ( kind = 8  ) normal_vector(1,1)
  character ( len = *   ) obj_file_name
  integer   ( kind = 4  ) obj_file_status
  integer   ( kind = 4  ) obj_file_unit
  integer   ( kind = 4  ) obj_group_num
  integer   ( kind = 4  ) obj_text_num
  logical                 s_begin
  logical                 s_eqi
  character ( len = *   ) tec_file_name
  integer   ( kind = 4  ) tec_file_status
  integer   ( kind = 4  ) tec_file_unit
  integer   ( kind = 4  ) tec_zone_num
  integer   ( kind = 4  ) text_num
  character ( len = 255 ) title
  integer   ( kind = 4  ) variable
  character ( len = 20  ), allocatable, dimension ( : ) :: variable_name
  integer   ( kind = 4  ) variable_num
  integer   ( kind = 4  ) vertex_normal(1,1)
  integer   ( kind = 4  ) zone
  integer   ( kind = 4  ) zone_num
!
!  #1: Read the TEC file once just to get "header" information.
!
  call tec_header_read ( tec_file_name, title, variable_num, dim_num, &
    node_data_num, zone_num, node_total, element_total, &
    element_order_uniform )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_HEADER_READ:'
  write ( *, '(a)' ) '  TITLE = "' // trim ( title ) // '"'
  write ( *, '(a,i8)' ) '  Number of variables = ', variable_num
  write ( *, '(a,i8)' ) '    Spatial variables = ', dim_num
  write ( *, '(a,i8)' ) '    Other variables   = ', node_data_num
  write ( *, '(a,i8)' ) '  Number of zones =    ', zone_num
  write ( *, '(a,i8)' ) '  Number of nodes =    ', node_total
  write ( *, '(a,i8)' ) '  Number of elements = ', element_total
  write ( *, '(a,i8)' ) '  Order of elements =  ', element_order_uniform
!
!  #2: Advance the TEC file up to the first ZONE record.
!
  call get_unit ( tec_file_unit )

  open ( unit = tec_file_unit, file = tec_file_name, iostat = tec_file_status, &
    status = 'old' )

  call tec_header_skip ( tec_file_unit, line )
!
!  #3: Write the beginning of the OBJ file.
!
  call get_unit ( obj_file_unit )

  open ( unit = obj_file_unit, file = obj_file_name, iostat = obj_file_status, &
    status = 'replace' )

  call obj_header_write ( obj_file_name, obj_file_unit, tec_file_name, obj_text_num )
!
!  #4: Read each zone, and write it to the OBJ file.
!
  element_offset = 0
  node_offset = 0
  obj_group_num = 0
  
  do zone = 1, zone_num

    call tec_zone_header_read ( tec_file_unit, line, node_num, element_num )

    allocate ( node_coord(dim_num,node_num) )
    allocate ( node_data(1:node_data_num,node_num) )
    allocate ( element_node(element_order_uniform,element_num) )
    allocate ( element_order(element_num) )

    call tec_zone_data_read ( tec_file_unit, line, dim_num, node_num, &
      element_num, element_order_uniform, node_data_num, node_coord, element_node, &
      node_data )

    element_order(1:element_num) = element_order_uniform

    normal_num = 0

    call obj_data_write ( obj_file_unit, node_num, element_num, normal_num, &
      element_order_uniform, node_coord, element_order, element_node, &
      normal_vector, vertex_normal, element_offset, node_offset, &
      obj_group_num, obj_text_num )

    element_offset = element_offset + element_num
    node_offset = node_offset + node_num

    deallocate ( node_coord )
    deallocate ( node_data )
    deallocate ( element_node )
    deallocate ( element_order )

  end do
!
!  #5: All done.  Close the files.
!
  close ( unit = tec_file_unit )
  close ( unit = obj_file_unit )

  return
end
subroutine tec_zone_header_read ( tec_file_unit, line, node_num, element_num )

!*****************************************************************************80
!
!! TEC_ZONE_HEADER_READ reads the header information for one zone.
!
!  Discussion:
!
!    On input, LINE contains the first record of the zone header, and
!    there may be more in the file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEC_FILE_UNIT, the unit number for the open TEC file.
!
!    Input/output, character ( len = * ) LINE.  On output, the text of
!    a line of the file which indicates the beginning of a ZONE header.
!    On output, the first line after the ZONE header.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer   ( kind = 4 )  begin
  logical                 ch_eqi
  integer   ( kind = 4 )  dim_num
  integer   ( kind = 4 )  element_num
  integer   ( kind = 4 )  element_order
  character ( len = 40 )  element_type
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  length
  character ( len = 255 ) line
  integer   ( kind = 4 )  line_num
  character ( len = 20 )  name
  integer   ( kind = 4 )  node_num
  integer   ( kind = 4 )  read_status
  logical                 s_begin
  logical                 s_eqi
  integer   ( kind = 4 )  tec_file_status
  integer   ( kind = 4 )  tec_file_unit
  character ( len = 255 ) title
  character ( len = 40  ) value
!
!  We got here because LINE contains text that was rejected by the VARIABLES
!  processor.  So it probably already contains the beginning of our ZONE record.
!
!  READ_STATUS = 0: We haven't recognized a ZONE record yet.
!  READ_STATUS = 1: We've read "ZONE", and haven't seen the end of the zone stuff.
!  READ_STATUS = 2: We're reading (and ignoring) numeric data for a zone.
!
  read_status = 0

  do
!
!  Ignore blank lines.
!
    if ( 0 < len_trim ( line ) ) then
!
!  READ_STATUS = 0 expects a ZONE record.
!
      if ( read_status == 0 ) then

        if ( s_begin ( line, 'ZONE' ) ) then
          call s_behead_substring ( line, 'ZONE' )
          read_status = 1
          node_num = -1
          element_num = -1
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEC_ZONE_HEADER_READ:'
          write ( *, '(a)' ) '  Puzzled by line "' // trim ( line ) // '".'
          stop
        end if
!
!  READ_STATUS = 2 expects more numbers, or ZONE.
!
      else if ( read_status == 2 ) then

        if ( s_begin ( line, 'ZONE' ) ) then
          call s_behead_substring ( line, 'ZONE' )
          read_status = 1
          node_num = -1
          element_num = -1
        end if

      end if
!
!  If we are in READ_STATUS = 1, but the line contains no equal signs,
!  we probably need to move to READ_STATUS = 2 and RETURN.
!
      if ( read_status == 1 ) then
        if ( 0 < len_trim ( line ) ) then
          if ( index ( line, '=' ) <= 0 ) then
            read_status = 2
            return
          end if
        end if
      end if

      if ( read_status == 1 ) then
!
!  This part of code is only reached for READ_STATUS = 1.
!  LINE may be blank, or contain an unknown number of pairs of
!  "Name=Value" items,
!
!  We are particularly interested in node, element and element order information.
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
      do

        call s_word_extract ( line, name )

        if ( len_trim ( name ) <= 0 ) then
          exit
        end if

        call s_word_extract ( line, value )

        if ( len_trim ( value ) == 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEC_ZONE_HEADER_READ - Fatal error!'
          write ( *, '(a)' ) '  Unexpected End of input.'
          write ( *, '(a)' ) '  Pretend nothing wrong...'
          value = '0'
        end if

        if ( ( ch_eqi ( name, 'N' ) .or. ch_eqi ( name, 'NODES' ) ) .and. &
             node_num == -1 ) then

          call s_to_i4 ( value, node_num, ierror, length )

        elseif ( ( ch_eqi ( name, 'E' ) .or. ch_eqi ( name, 'ELEMENTS' ) ) .and. &
          element_num == -1 ) then

          call s_to_i4 ( value, element_num, ierror, length )

        elseif ( s_eqi ( name, 'DATAPACKING' ) ) then

          if ( .not. s_eqi ( value, 'POINT' ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TEC_ZONE_HEADER_READ - Fatal error!'
            write ( *, '(a)' ) '  Unacceptable DATAPACKING value.'
            write ( *, '(a)' ) '  Only "DATAPACKING = POINT" is supported.'
            stop
          end if

        elseif ( s_eqi ( name, 'ZONETYPE' ) .and. &
          len_trim ( element_type ) == 0 ) then

        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Ignoring "' // trim ( name ) &
            // '" = "' // trim ( value ) // '".'

        end if

      end do

      end if

    end if
!
!  Read the next line.
!
    read ( tec_file_unit, '(a)', iostat = tec_file_status ) line

    if ( tec_file_status /= 0 ) then
      exit
    end if

  end do

 return
end
subroutine tec_zone_line_parse ( line, node_num, element_num, element_type, &
  element_order )

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
!    by spaces or commas and even carriage returns:
!
!      N = number of nodes 
!   or NODES = number of nodes
!
!      E = number of elements
!   or ELEMENTS = number of elements
!
!      T = optional zone title (we can't handle this right now)
!
!      PACKING = 'POINT'
!
!      ZONETYPE = 'FETRIANGLE' 
!   or ZONETYPE = 'FEQUADRILATERAL'
!   or ZONETYPE = 'FETETRAHEDRON' 
!   or ZONETYPE = 'FEBRICK'.
!
!    Other arguments that may appear will be ignored.
!
!    The end of the ZONE record is only implicitly indicated, by the 
!    occurrence of a line of numeric data!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string of characters, representing 
!    the "VARIABLES=" line of the file.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, character ( len = * ) ELEMENT_TYPE, the element type: 
!    'FETRIANGLE' or 'FEQUADRILATERAL' or 'FETETRAHEDRON' or 'FEBRICK'.
!
  implicit none

  logical ch_eqi
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  character ( len = * ) element_type
  integer ( kind = 4 ) found_num
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  character ( len = * ) line
  character ( len = 80 ) name
  integer ( kind = 4 ) node_num
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

    if ( ( ch_eqi ( name, 'N' ) .or. ch_eqi ( name, 'NODES' ) ) .and. &
         node_num == -1 ) then

      call s_to_i4 ( value, node_num, ierror, length )
      found_num = found_num + 1

    elseif ( ( ch_eqi ( name, 'E' ) .or. ch_eqi ( name, 'ELEMENTS' ) ) .and. &
      element_num == -1 ) then

      call s_to_i4 ( value, element_num, ierror, length )
      found_num = found_num + 1

    elseif ( s_eqi ( name, 'DATAPACKING' ) ) then

      if ( .not. s_eqi ( value, 'POINT' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEC_ZONE_LINE_PARSE - Fatal error!'
        write ( *, '(a)' ) '  Unacceptable DATAPACKING value.'
        write ( *, '(a)' ) '  Only "DATAPACKING = POINT" is supported.'
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
