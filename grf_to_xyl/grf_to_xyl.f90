program main

!*****************************************************************************80
!
!! MAIN is the main program for GRF_TO_XYL.
!
!  Discussion:
!
!    GRF_TO_XYL takes a GRF file and creates XY and XYL files that may
!    be used to make a plot of the graph.
!
!  Usage:
!
!    grf_to_xyl file.grf
!
!    creates the files "file.xy" and "file.xyl".
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
  implicit none

  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  character ( len = 255 ) grf_filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) num_arg
  character ( len = 255 ) xy_filename
  character ( len = 255 ) xyl_filename

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRF_TO_XYL'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a GRF file;'
  write ( *, '(a)' ) '  Create corresponding XY (point) and XYL (line) files.'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input file name:'
    read ( *, '(a)', iostat = ios ) grf_filename

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_TO_XYL - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, grf_filename )

  end if
!
!  Construct the output file names.
!
  xy_filename = grf_filename
  call file_name_ext_swap ( xy_filename, 'xy' )

  xyl_filename = grf_filename
  call file_name_ext_swap ( xyl_filename, 'xyl' )
!
!  Now we know what to do.
!
  call grf_to_xyl ( grf_filename, xy_filename, xyl_filename )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XY  information written to "' // trim ( xy_filename ) // '",'
  write ( *, '(a)' ) '  XYL information written to "' // trim ( xyl_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRF_TO_XYL'
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
!! CH_TO_DIGIT returns the value of a base 10 digit.
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

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(edge_num)
  integer ( kind = 4 ) edge_pointer(node_num+1)
  integer ( kind = 4 ) node
  real ( kind = 8 ) xy(2,node_num)

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

  open ( unit = input_unit, file = input_filename, status = 'old', iostat = ios )

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

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) node_i
  integer ( kind = 4 ) node_j
  integer ( kind = 4 ) node_num
  character ( len = 255 ) text
  integer ( kind = 4 ) text_pos
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval

  edge_num = -1
  node_num = -1

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', iostat = ios )

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
subroutine grf_to_xyl ( grf_filename, xy_filename, xyl_filename )

!*****************************************************************************80
!
!! GRF_TO_XYL converts GRF information to XYL information.
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
!    Input, character ( len = * ) GRF_FILENAME, the name of the GRF file.
!
!    Input, character ( len = * ) XY_FILENAME, XYL_FILENAME, the names of the
!    XY and XYL files to be created.
!
  implicit none

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ), allocatable :: edge_data(:)
  integer ( kind = 4 ), allocatable :: edge_pointer(:)
  character ( len = * )  grf_filename
  integer ( kind = 4 ) line
  integer ( kind = 4 ), allocatable :: line_data(:)
  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ), allocatable :: line_pointer(:)
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  real ( kind = 8 ), allocatable :: xy(:,:)
  character ( len = * )  xy_filename
  character ( len = * )  xyl_filename

  call grf_header_read ( grf_filename, node_num, edge_num )

  if ( .true. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  GRF data:'
    call grf_header_print ( node_num, edge_num )
  end if

  allocate ( edge_pointer(1:node_num+1) )
  allocate ( edge_data(1:edge_num) )
  allocate ( xy(1:2,1:node_num) )

  call grf_data_read ( grf_filename, node_num, edge_num, edge_pointer, &
    edge_data, xy )

  if ( .false. ) then
    call grf_data_print ( node_num, edge_num, edge_pointer, edge_data, xy )
  end if

  line_num = edge_num
  line_data_num = 2 * edge_num

  allocate ( line_pointer(1:line_num+1) )
  allocate ( line_data(1:line_data_num) )

  line = 0
  line_data_num = 0
  line_pointer(1) = 1

  do node1 = 1, node_num
    do edge = edge_pointer(node1), edge_pointer(node1+1) - 1
      node2 = edge_data(edge)
      line = line + 1
      line_data_num = line_data_num + 1
      line_data(line_data_num) = node1
      line_data_num = line_data_num + 1
      line_data(line_data_num) = node2
      line_pointer(line+1) = line_data_num + 1    
    end do
  end do

  if ( .false. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  XY data:'
    call xy_header_print ( node_num )
  end if

  if ( .false. ) then
    call xy_data_print ( node_num, xy )
  end if

  call xy_write ( xy_filename, node_num, xy )

  if ( .true. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  XYL data:'
    call xyl_header_print ( node_num, line_num, line_data_num )
  end if

  if ( .false. ) then
    call xyl_data_print ( node_num, line_num, line_data_num, &
      line_pointer, line_data )
  end if

  call xyl_write ( xyl_filename, node_num, line_num, line_data_num, &
    line_pointer, line_data )

  deallocate ( edge_data )
  deallocate ( edge_pointer )
  deallocate ( line_data )
  deallocate ( line_pointer )
  deallocate ( xy )

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
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
  character ( len = * )  s
  integer ( kind = 4 ) s_index_last
  character ( len = * )  sub

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

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * )  s
  integer ( kind = 4 ) state
  character :: TAB = achar ( 9 )
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

  character c
  logical ch_eqi
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
  character :: TAB = achar ( 9 )

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
  character ( len = 8 )  date
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
  character ( len = 5 )  zone

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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = *  ) string
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

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine xy_data_print ( point_num, xy )

!*****************************************************************************80
!
!! XY_DATA_PRINT prints the data of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) j
  real    ( kind = 8 ) xy(2,point_num)

  do j = 1, point_num
    write ( *, '(2x,f14.6,2x,f14.6)' ) xy(1:2,j)
  end do

  return
end
subroutine xy_data_write ( output_unit, point_num, xy )

!*****************************************************************************80
!
!! XY_DATA_WRITE writes the data of an XY file.
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
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) j
  integer ( kind = 4 ) output_unit
  real    ( kind = 8 ) xy(2,point_num)

  do j = 1, point_num
    write ( output_unit, '(2x,f14.6,2x,f14.6)' ) xy(1:2,j)
  end do

  return
end
subroutine xy_header_print ( point_num )

!*****************************************************************************80
!
!! XY_HEADER_PRINT prints the header of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
  implicit none

  integer ( kind = 4 ) point_num

  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  return
end
subroutine xy_header_write ( output_filename, output_unit, point_num )

!*****************************************************************************80
!
!! XY_HEADER_WRITE writes the header of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2003
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
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
  implicit none

  character ( len = * )  output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) point_num
  character ( len = 40 ) string

  call timestring ( string )
!
!  Write the header.
!
  write ( output_unit, '(a)'    ) '#  ' // trim ( output_filename )
  write ( output_unit, '(a)'    ) '#  created by xy_io::xy_header_write.f90'
  write ( output_unit, '(a)'    ) '#  on ' // trim ( string )
  write ( output_unit, '(a)'    ) '#'
  write ( output_unit, '(a,i8)' ) '#  Number of points = ', point_num
  write ( output_unit, '(a)'    ) '#'

  return
end
subroutine xy_write ( output_filename, point_num, xy )

!*****************************************************************************80
!
!! XY_WRITE writes an XY file.
!
!  Example:
!
!    # my_file.xy
!    # created by XY_IO::XY_WRITE.
!    #
!    #  Number of points = 5
!    #
!    0.0  0.0
!    1.0  2.0
!    3.0  5.0
!    2.0  1.0
!    8.0  7.5
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
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) point_num

  logical, parameter ::  debug = .false.
  integer ( kind = 4 ) ios
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) xy(2,point_num)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '".'
    stop
  end if
!
!  Write the header.
!
  call xy_header_write ( output_filename, output_unit, point_num )
!
!  Write the data.
!
  call xy_data_write ( output_unit, point_num, xy )
!
!  Close the file.
!
  close ( unit = output_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_WRITE - Note:'
    write ( *, '(a)' ) '  The data was written.'
    write ( *, '(a,i8)' ) '  Number of points =    ', point_num
  end if

  return
end
subroutine xyl_data_print ( point_num, line_num, line_data_num, &
  line_pointer, line_data )

!*****************************************************************************80
!
!! XYL_DATA_PRINT prints the data of an XYL file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
!    Input, integer ( kind = 4 ) LINE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each line.
!
!    Input, integer ( kind = 4 ) LINE_DATA(LINE_DATA_NUM), indices
!    of points that form lines.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) line
  integer ( kind = 4 ) line_data(line_data_num)
  integer ( kind = 4 ) line_pointer(line_num+1)

  do line = 1, line_num
    write ( *, '(2x,i4,2x,i8,2x,i8)' ) &
      line, line_pointer(line), line_pointer(line+1)-1
  end do

  write ( *, '(a)' ) ' '

  do line = 1, line_num
    do i = line_pointer(line), line_pointer(line+1) - 1
      write ( *, '(2x,i8)', advance = 'NO' ) line_data(i)
    end do
    write ( *, '(1x)', advance = 'YES' )
  end do

  return
end
subroutine xyl_data_write ( output_unit, point_num, line_num, line_data_num, &
  line_pointer, line_data )

!*****************************************************************************80
!
!! XYL_DATA_WRITE writes the data of an XYL file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
!    Input, integer ( kind = 4 ) LINE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each line.
!
!    Input, integer ( kind = 4 ) LINE_DATA(LINE_DATA_NUM), indices
!    of points that form lines.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) line
  integer ( kind = 4 ) line_data(line_data_num)
  integer ( kind = 4 ) line_pointer(line_num+1)
  integer ( kind = 4 ) output_unit

  do line = 1, line_num
    do i = line_pointer(line), line_pointer(line+1) - 1
      write ( output_unit, '(2x,i8)', advance = 'NO' ) line_data(i)
    end do
    write ( output_unit, '(1x)', advance = 'YES' )
  end do

  return
end
subroutine xyl_header_print ( point_num, line_num, line_data_num )

!*****************************************************************************80
!
!! XYL_HEADER_PRINT prints the header of an XYL file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of points =     ', point_num
  write ( *, '(a,i8)' ) '  Number of lines =      ', line_num
  write ( *, '(a,i8)' ) '  Number of line items = ', line_data_num

  return
end
subroutine xyl_header_write ( output_filename, output_unit, point_num, &
  line_num, line_data_num )

!*****************************************************************************80
!
!! XYL_HEADER_WRITE writes the header of an XYL file.
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
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) point_num
  character ( len = 40 ) string

  call timestring ( string )
!
!  Write the header.
!
  write ( output_unit, '(a)'    ) '#  ' // trim ( output_filename )
  write ( output_unit, '(a)'    ) '#  created by xy_io::xyl_header_write.f90'
  write ( output_unit, '(a)'    ) '#  on ' // trim ( string )
  write ( output_unit, '(a)'    ) '#'
  write ( output_unit, '(a,i8)' ) '#  Number of points =     ', point_num
  write ( output_unit, '(a,i8)' ) '#  Number of lines =      ', line_num
  write ( output_unit, '(a,i8)' ) '#  Number of line items = ', line_data_num
  write ( output_unit, '(a)'    ) '#'

  return
end
subroutine xyl_write ( output_filename, point_num, line_num, line_data_num, &
  line_pointer, line_data )

!*****************************************************************************80
!
!! XYL_WRITE writes an XYL file.
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
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
!    Input, integer ( kind = 4 ) LINE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each line.
!
!    Input, integer ( kind = 4 ) LINE_DATA(LINE_DATA_NUM), indices
!    of points that form lines.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) line_data(line_data_num)
  integer ( kind = 4 ) line_pointer(line_num+1)
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_unit
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYL_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '".'
    stop
  end if
!
!  Write the header.
!
  call xyl_header_write ( output_filename, output_unit, point_num, line_num, &
    line_data_num )
!
!  Write the data.
!
  call xyl_data_write ( output_unit, point_num, line_num, line_data_num, &
    line_pointer, line_data )
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
