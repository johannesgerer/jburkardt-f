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
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer ( kind = 4 ) value of a base 10 digit.
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding value.  If C was
!    'illegal', then DIGIT is -1.
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
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

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
subroutine grf_data_print ( node_num, edge_num, edge_pointer, edge_data, xy )

!*****************************************************************************80
!
!! GRF_DATA_PRINT prints the data of a GRF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Skiena,
!    Implementing Discrete Mathematics,
!    Combinatorics and Graph Theory with Mathematica,
!    Addison-Wesley, 1990.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_POINTER(NODE_NUM+1), pointers to 
!    the beginning of edge data for each node.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(EDGE_NUM), the edge data.
!
!    Input, real ( kind = 8 ) XY(2,NODE_NUM), the node coordinates.
!
  implicit none

  integer ( kind = 4 )  edge_num
  integer ( kind = 4 )  node_num

  integer ( kind = 4 )  edge
  integer ( kind = 4 )  edge_data(edge_num)
  integer ( kind = 4 )  edge_pointer(node_num+1)
  integer ( kind = 4 )  node
  real ( kind = 8 )  xy(2,node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Edge pointers:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node     First      Last'
  write ( *, '(a)' ) ' '
  do node = 1, node_num
    write ( *, '(2x,i4,2x,i8,2x,i8)' ) &
      node, edge_pointer(node), edge_pointer(node+1)-1
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Edge data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node     Adjacent nodes'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i4)', advance = 'NO' ) node
    do edge = edge_pointer(node), edge_pointer(node+1) - 1
      write ( *, '(2x,i8)', advance = 'NO' ) edge_data(edge)
    end do
    write ( *, '(1x)', advance = 'YES' )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node        X          Y'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i4,2x,f10.6,2x,f10.6)' ) node, xy(1:2,node)
  end do

  return
end
subroutine grf_data_read ( input_filename, node_num, edge_num, edge_pointer, &
  edge_data, xy )

!*****************************************************************************80
!
!! GRF_DATA_READ reads the data of a GRF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Skiena,
!    Implementing Discrete Mathematics,
!    Combinatorics and Graph Theory with Mathematica,
!    Addison-Wesley, 1990.
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Output, integer ( kind = 4 ) EDGE_POINTER(NODE_NUM+1), pointers to 
!    the beginning of edge data for each node.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(EDGE_NUM), the edge data.
!
!    Output, real ( kind = 8 ) XY(2,NODE_NUM), the node coordinates.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(edge_num)
  integer ( kind = 4 ) edge_pointer(node_num+1)
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) node_i
  integer ( kind = 4 ) node_j
  character ( len = 255 ) text
  integer ( kind = 4 ) text_pos
  real ( kind = 8 ) xy(2,node_num)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval

  edge_data(1:edge_num) = -1
  edge_pointer(1:node_num+1) = -1

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) &
      '  Could not open the input file "' // trim ( input_filename ) // '".'
    stop
  end if
!
!  Read information about each node.
!
  edge = 0
  edge_pointer(1) = 1

  do

    read ( input_unit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( text ) <= 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if

    text_pos = 1
!
!  Extract the node index, NODE.
!
    call s_to_i4 ( text(text_pos:), node_i, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable node index value.'
      stop
    end if

    text_pos = text_pos + lchar

    edge_pointer(node_i+1) = edge_pointer(node_i)
!
!  Extract the X, Y coordinates of the node.
!
    call s_to_r8 ( text(text_pos:), xval, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable X coordinate for node.'
      stop
    end if

    xy(1,node_i) = xval

    text_pos = text_pos + lchar

    call s_to_r8 ( text(text_pos:), yval, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable Y coordinate for node.'
      stop
    end if

    xy(2,node_i) = yval

    text_pos = text_pos + lchar
!
!  Read the indices of the nodes to which the node is connected.
!
    do

      call s_to_i4 ( text(text_pos:), node_j, ierror, lchar )

      if ( ierror /= 0 .and. ierror /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unreadable node neighbor value.'
        stop
      end if

      text_pos = text_pos + lchar

      if ( lchar <= 0 ) then
        exit
      end if

      if ( node_j <= 0 ) then
        cycle
      end if

      edge = edge + 1
      edge_data(edge) = node_j
      edge_pointer(node_i+1) = edge_pointer(node_i+1) + 1

    end do

  end do

  close ( unit = input_unit )

  return
end
subroutine grf_data_write ( output_unit, node_num, edge_num, edge_pointer, &
  edge_data, xy )

!*****************************************************************************80
!
!! GRF_DATA_WRITE prints the data to a GRF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Skiena,
!    Implementing Discrete Mathematics,
!    Combinatorics and Graph Theory with Mathematica,
!    Addison-Wesley, 1990.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_POINTER(NODE_NUM+1), pointers to 
!    the beginning of edge data for each node.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(EDGE_NUM), the edge data.
!
!    Input, real ( kind = 8 ) XY(2,NODE_NUM), the node coordinates.
!
  implicit none

  integer ( kind = 4 )  edge_num
  integer ( kind = 4 )  node_num

  integer ( kind = 4 )  edge
  integer ( kind = 4 )  edge_data(edge_num)
  integer ( kind = 4 )  edge_pointer(node_num+1)
  integer ( kind = 4 )  node
  integer ( kind = 4 )  output_unit
  real ( kind = 8 )  xy(2,node_num)

  do node = 1, node_num
    write ( output_unit, '(i4,2x,f10.6,2x,f10.6)', advance = 'no' ) &
      node, xy(1:2,node)
    do edge = edge_pointer(node), edge_pointer(node+1) - 1
      write ( output_unit, '(2x,i4)', advance = 'NO' ) edge_data(edge)
    end do
    write ( output_unit, '(1x)', advance = 'YES' )
  end do

  return
end
subroutine grf_example ( node_num, edge_num, edge_pointer, edge_data, xy )

!*****************************************************************************80
!
!! GRF_EXAMPLE sets up a GRF example.
!
!  Discussion:
!
!    The example is known as the Coxeter graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Skiena,
!    Implementing Discrete Mathematics,
!    Combinatorics and Graph Theory with Mathematica,
!    Addison-Wesley, 1990.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Output, integer ( kind = 4 ) EDGE_POINTER(NODE_NUM+1), pointers to 
!    the beginning of edge data for each node.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(EDGE_NUM), the edge data.
!
!    Output, real ( kind = 8 ) XY(2,NODE_NUM), the node coordinates.
!
  implicit none

  integer ( kind = 4 )  edge_num
  integer ( kind = 4 )  node_num

  integer ( kind = 4 )  edge
  integer ( kind = 4 )  edge_data(edge_num)
  integer ( kind = 4 )  edge_pointer(node_num+1)
  integer ( kind = 4 )  node
  real ( kind = 8 )  xy(2,node_num)

  edge_pointer = (/ &
    1,  4,  7, 10, 13, 16, 19, 22, 25, 28, &
   31, 34, 37, 40, 43, 46, 49, 52, 55, 58, &
   61, 64, 67, 70, 73, 76, 79, 82, 85 /)

  edge_data = (/ &
     8,   2,   3, &
    14,   1,   5, &
     9,   4,   1, &
    10,   7,   3, &
    13,   2,   6, &
    12,   5,   7, &
    11,   6,   4, &
    25,  20,   1, &
    24,  21,   3, &
    23,  15,   4, &
    22,  16,   7, &
    28,  17,   6, &
    27,  18,   5, &
    26,  19,   2, &
    10,  18,  19, &
    11,  19,  20, &
    12,  21,  20, &
    13,  15,  21, &
    14,  16,  15, &
     8,  17,  16, &
     9,  18,  17, &
    11,  27,  24, &
    10,  28,  25, &
     9,  26,  22, &
     8,  23,  27, &
    14,  24,  28, & 
    13,  25,  22, &  
    12,  26,  23 /)

  xy = reshape ( (/ &
    0.412,   0.984, &
    0.494,   0.984, &
    0.366,   0.926, &
    0.388,   0.862, &
    0.546,   0.926, &
    0.518,   0.860, &
    0.458,   0.818, &
    0.152,   0.684, &
    0.264,   0.682, &
    0.354,   0.680, &
    0.458,   0.670, &
    0.554,   0.672, &
    0.658,   0.668, &
    0.774,   0.692, &
    0.164,   0.450, &
    0.228,   0.448, &
    0.274,   0.390, & 
    0.242,   0.330, &
    0.194,   0.278, &
    0.146,   0.328, & 
    0.102,   0.390, &
    0.668,   0.472, &
    0.638,   0.416, & 
    0.656,   0.334, &   
    0.714,   0.270, &  
    0.798,   0.326, &  
    0.830,   0.408, & 
    0.754,   0.466 /), (/ 2, node_num/) )

  return
end
subroutine grf_example_size ( node_num, edge_num )

!*****************************************************************************80
!
!! GRF_EXAMPLE_SIZE sizes a GRF example.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Skiena,
!    Implementing Discrete Mathematics,
!    Combinatorics and Graph Theory with Mathematica,
!    Addison-Wesley, 1990.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) node_num

  node_num = 28
  edge_num = 84

  return
end
subroutine grf_header_print ( node_num, edge_num )

!*****************************************************************************80
!
!! GRF_HEADER_PRINT prints the header of a GRF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Skiena,
!    Implementing Discrete Mathematics,
!    Combinatorics and Graph Theory with Mathematica,
!    Addison-Wesley, 1990.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes NODE_NUM = ', node_num
  write ( *, '(a,i8)' ) '  The number of edges EDGE_NUM = ', edge_num

  return
end
subroutine grf_header_read ( input_filename, node_num, edge_num )

!*****************************************************************************80
!
!! GRF_HEADER_READ reads the header of a GRF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Skiena,
!    Implementing Discrete Mathematics,
!    Combinatorics and Graph Theory with Mathematica,
!    Addison-Wesley, 1990.
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
  implicit none

  integer ( kind = 4 )  edge_num
  integer ( kind = 4 )  ierror
  character ( len = * )   input_filename
  integer ( kind = 4 )  input_unit
  integer ( kind = 4 )  ios
  integer ( kind = 4 )  lchar
  integer ( kind = 4 )  node_i
  integer ( kind = 4 )  node_j
  integer ( kind = 4 )  node_num
  character ( len = 255 ) text
  integer ( kind = 4 )  text_pos
  real ( kind = 8 )  xval
  real ( kind = 8 )  yval

  edge_num = -1
  node_num = -1

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) &
      '  Could not open the input file "' // trim ( input_filename ) // '".'
    stop
  end if
!
!  Read information about each node.
!
  node_num = 0
  edge_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( text ) <= 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if

    text_pos = 1
!
!  Extract the node index, NODE.
!
    call s_to_i4 ( text(text_pos:), node_i, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable node index value.'
      stop
    end if

    node_num = node_num + 1

    text_pos = text_pos + lchar

    if ( node_i /= node_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Record ', node_num, ' is labeled ', node_i
      write ( *, '(a)' ) '  but these values should be equal.'
      stop
    end if
!
!  Extract the X, Y coordinates of the node.
!
    call s_to_r8 ( text(text_pos:), xval, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable X coordinate for node.'
      stop
    end if

    text_pos = text_pos + lchar

    call s_to_r8 ( text(text_pos:), yval, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable Y coordinate for node.'
      stop
    end if

    text_pos = text_pos + lchar
!
!  Read the indices of the nodes to which the node is connected.
!
    do

      call s_to_i4 ( text(text_pos:), node_j, ierror, lchar )

      if ( ierror /= 0 .and. ierror /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unreadable node neighbor value.'
        stop
      end if

      text_pos = text_pos + lchar

      if ( lchar <= 0 ) then
        exit
      end if

      if ( node_j <= 0 ) then
        cycle
      end if

      edge_num = edge_num + 1

    end do

  end do

  close ( unit = input_unit )

  return
end
subroutine grf_header_write ( output_filename, output_unit, node_num, &
  edge_num )

!*****************************************************************************80
!
!! GRF_HEADER_WRITE writes the header of a GRF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) node_num
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
!
!  Write the header.
!
  write ( output_unit, '(a)'    ) '#  ' // trim ( output_filename )
  write ( output_unit, '(a)'    ) '#  created by grf_io::grf_header_write.f90'
  write ( output_unit, '(a)'    ) '#'
  write ( output_unit, '(a,i8)' ) '#  Number of nodes  =     ', node_num
  write ( output_unit, '(a,i8)' ) '#  Number of edges =      ', edge_num
  write ( output_unit, '(a)'    ) '#'

  return
end
subroutine grf_write ( output_filename, node_num, edge_num, edge_pointer, &
  edge_data, xy )

!*****************************************************************************80
!
!! GRF_WRITE writes a GRF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the file
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_POINTER(NODE_NUM+1), pointers to the
!    first edge item for each node.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(EDGE_NUM), indices of adjacent nodes.
!
!    Input, real ( kind = 8 ) XY(2,NODE_NUM), the node coordinates.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) edge_data(edge_num)
  integer ( kind = 4 ) edge_pointer(node_num+1)
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) xy(2,node_num)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRF_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '".'
    stop
  end if
!
!  Write the header.
!
  if ( .false. ) then
    call grf_header_write ( output_filename, output_unit, node_num, edge_num )
  end if
!
!  Write the data.
!
  call grf_data_write ( output_unit, node_num, edge_num, edge_pointer, &
    edge_data, xy )
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2009
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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters 
!    of S used to make the integer.
!
  implicit none

  character              c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * )  s
  integer ( kind = 4 ) state
  character              :: TAB = achar ( 9 )
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

      if ( c == ' ' .or. c == TAB ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

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
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
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
!    12 January 2009
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
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character              c
  logical                ch_eqi
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * )  s
  integer ( kind = 4 ) s_length
  character           :: TAB = achar ( 9 )

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
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
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
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
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

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
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

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
