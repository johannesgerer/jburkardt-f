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
subroutine poly_header_read ( poly_file_name, node_num, edge_num, hole_num, &
  region_num )

!*****************************************************************************80
!
!! POLY_HEADER_READ reads header data from a POLY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) POLY_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
!    Output, integer ( kind = 4 ) REGION_NUM, the number of regions.
!
  implicit none

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) hole
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  character ( len = 255 ) line
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  character ( len = * )   poly_file_name
  integer ( kind = 4 )  poly_file_status
  integer ( kind = 4 )  poly_file_unit
  integer ( kind = 4 )  region
  integer ( kind = 4 )  region_num
  character ( len = 255 ) word

  line_num = 0

  node = 0
  node_num = -1
  edge = 0
  edge_num = -1
  hole = 0
  hole_num = -1
  region = 0
  region_num = -1

  call get_unit ( poly_file_unit )

  open ( unit = poly_file_unit, file = poly_file_name, status = 'old', &
    iostat = poly_file_status )

  if ( poly_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLY_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input POLY file:'
    write ( *, '(a)' ) '  "' //  trim ( poly_file_name ) // '".'
    stop
  end if
!
!  Read lines, skip comments, search for header values.
!
  do

    read ( poly_file_unit, '(a)', iostat = poly_file_status ) line

    if ( poly_file_status /= 0 ) then

      if ( region_num /= -1 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLY_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading POLY file:'
      write ( *, '(a)' ) '  "' //  trim ( poly_file_name ) // '".'
      write ( *, '(a,i8)' ) '  Line number of file was ', line_num
      stop

    end if

    line_num = line_num + 1
!
!  Skip comment lines.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  If NODE_NUM = -1, then you expect the line with node numbers.
!
    if ( node_num == -1 ) then

      call s_word_extract ( line, word )

      if ( len_trim ( word ) == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLY_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if

      call s_to_i4 ( word, node_num, ierror, length )
!
!  If 0 < NODE_NUM and NODE < NODE_NUM,
!  you expect a line with node coordinates.
!  Without checking, we just assume this is a node definition line.
!
    else if ( node < node_num ) then

      node = node + 1
!
!  If 0 <= NODE_NUM and NODE == NODE_NUM, but
!  EDGE_NUM = -1, then you expect the line with edge numbers.
!
    else if ( edge_num == -1 ) then

      call s_word_extract ( line, word )

      if ( len_trim ( word ) == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLY_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if

      call s_to_i4 ( word, edge_num, ierror, length )
!
!  If 0 < EDGE_NUM, and EDGE < EDGE_NUM,
!  you expect a line with edge information.
!
    else if ( edge < edge_num ) then

      edge = edge + 1
!
!  If 0 <= NODE_NUM and NODE = NODE_NUM,
!  and 0 <= EDGE_NUM and EDGE = EDGE_NUM
!  and HOLE_NUM = -1,
!  you expect the line with hole numbers.
!
    else if ( hole_num == -1 ) then

      call s_word_extract ( line, word )

      if ( len_trim ( word ) == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLY_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if

      call s_to_i4 ( word, hole_num, ierror, length )
!
!  If 0 <= NODE_NUM and NODE = NODE_NUM,
!  and 0 <= EDGE_NUM and EDGE = EDGE_NUM
!  and 0 < HOLE_NUM,
!  you expect the lines with hole points.
!
    else if ( hole < hole_num ) then

      hole = hole + 1
!
!  If 0 <= NODE_NUM and NODE = NODE_NUM,
!  and 0 <= EDGE_NUM and EDGE = EDGE_NUM
!  and 0 <= HOLE_NUM and HOLE = HOLE_NUM,
!  you expect the line with region numbers.
!
    else if ( region_num <= -1 ) then

      call s_word_extract ( line, word )

      if ( len_trim ( word ) == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLY_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if

      call s_to_i4 ( word, region_num, ierror, length )
!
!  Unexpected lines.
!
    else

    end if

  end do

  close ( unit = poly_file_unit )

  return
end
subroutine poly_write ( poly_file_name, node_num, segment, edge_num, &
  edge_nodes, hole_num, hole_point )

!*****************************************************************************80
!
!! POLY_WRITE writes data to a POLY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) POLY_FILE_NAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) SEGMENT(2,NODE_NUM), the nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_NODES(2,EDGE_NUM), the nodes that
!    form each edge.
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
!    Input, real ( kind = 8 ) HOLE_POINT(2,HOLE_NUM), a point in each hole.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num

  character ( len = 8 ) date
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_nodes(2,edge_num)
  integer ( kind = 4 ) hole
  real ( kind = 8 ) hole_point(2,hole_num)
  integer ( kind = 4 ) node
  character ( len = * ) poly_file_name
  integer ( kind = 4 ) poly_file_status
  integer ( kind = 4 ) poly_file_unit
  integer ( kind = 4 ), parameter :: region_num = 0
  real ( kind = 8 ) segment(2,node_num)
  character ( len = 40 ) string

  call get_unit ( poly_file_unit )

  open ( unit = poly_file_unit, file = poly_file_name, status = 'replace', &
   iostat = poly_file_status )

  if ( poly_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLY_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output POLY file.'
    stop
  end if

  write ( poly_file_unit, '(a)' ) '#  ' // trim ( poly_file_name )
  write ( poly_file_unit, '(a)' ) '#  Created by poly_write.f90'
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(a)' ) '#                               Boundary'
  write ( poly_file_unit, '(a)' ) '#   Vertex Dimension Attribute    Marker'
  write ( poly_file_unit, '(a)' ) '#    Count     Count     Count     Count'
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(2x,i8,a)' ) &
    node_num, '         2         0         0'
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(a)' ) &
    '#   Vertex      X           Y       Attributes  Marker'
  write ( poly_file_unit, '(a)' ) '#    Index'
  write ( poly_file_unit, '(a)' ) '#'
  do node = 1, node_num
    write ( poly_file_unit, '(2x,i8,2x,f10.6,2x,f10.6)' ) &
      node, segment(1,node), segment(2,node)
  end do
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(a)' ) '#           Boundary'
  write ( poly_file_unit, '(a)' ) '#  Segment    Marker'
  write ( poly_file_unit, '(a)' ) '#    Count     Count'
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(2x,i8,a)' ) edge_num, '         0'
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(a)' ) '#  Segment     Node1     Node2    Marker'
  write ( poly_file_unit, '(a)' ) '#    Index'
  write ( poly_file_unit, '(a)' ) '#'
  do edge = 1, edge_num
    write ( poly_file_unit, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) &
      edge, edge_nodes(1,edge), edge_nodes(2,edge), 0
  end do
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(a)' ) '#     Hole'
  write ( poly_file_unit, '(a)' ) '#    Count'
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(2x,i8)' ) hole_num
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(a)' ) '#     Hole      X           Y'
  write ( poly_file_unit, '(a)' ) '#    index'
  write ( poly_file_unit, '(a)' ) '#'
  do hole = 1, hole_num
    write ( poly_file_unit, '(2x,i8,2x,f10.6,2x,f10.6)' ) &
      hole, hole_point(1,hole), hole_point(2,hole)
  end do
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(a)' ) '#   Region'
  write ( poly_file_unit, '(a)' ) '#    Count'
  write ( poly_file_unit, '(a)' ) '#'
  write ( poly_file_unit, '(2x,i8)' ) region_num

  close ( unit = poly_file_unit )

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
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S
!    used to make the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
  integer ( kind = 4 ) value

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
