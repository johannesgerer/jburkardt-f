program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM_TO_TEC.
!
!  Discussion:
!
!    FEM_TO_TEC reads a set of FEM files and writes a corresponding TEC file.
!
!  Usage:
!
!    fem_to_tec coord.txt element.txt data.txt file.tec
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) element_file_name
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  character ( len = 255 ) node_coord_file_name
  character ( len = 255 ) node_data_file_name
  integer ( kind = 4 ) num_arg
  character ( len = 255 ) tec_file_name

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_TO_TEC'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a set of FEM files;'
  write ( *, '(a)' ) '  Write a corresponding TEC file.'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input coordinate file name:'
    read ( *, '(a)', iostat = ios ) node_coord_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM_TO_TEC - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, node_coord_file_name )

  end if

  if ( num_arg < 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the element file name:'
    read ( *, '(a)', iostat = ios ) element_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM_TO_TEC - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 2

    call getarg ( iarg, element_file_name )

  end if

  if ( num_arg < 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the node data file name:'
    read ( *, '(a)', iostat = ios ) node_data_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM_TO_TEC - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 3

    call getarg ( iarg, node_data_file_name )

  end if

  if ( num_arg < 4 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the output TECPLOT file name:'
    read ( *, '(a)', iostat = ios ) tec_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM_TO_TEC - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 4

    call getarg ( iarg, tec_file_name )

  end if
!
!  Now we know what to do.
!
  call fem_to_tec_handle ( node_coord_file_name, element_file_name, &
    node_data_file_name, tec_file_name )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_TO_TEC'
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

  if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine fem_data_read ( node_coord_file_name, element_file_name, &
  node_data_file_name, dim_num, node_num, element_num, element_order, &
  node_data_num, node_coord, element_node, node_data )

!*****************************************************************************80
!
!! FEM_DATA_READ reads data from a set of FEM files.
!
!  Discussion:
!
!    This program reads the node, element and data files that define
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
!    02 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_COORD_FILE_NAME, the name of the node
!    coordinate file.  If this argument is not supplied, it will be requested.
!    If the interactive response is blank, or otherwise defective, then the
!    program terminates.
!
!    Input, character ( len = * ) ELEMENT_FILE_NAME, the name of the element
!    file.  If this argument is not supplied, it will be requested.  If the
!    interactive response is blank, then the program will assume that no
!    element information is to be supplied.  (But the node coordinates must
!    be available and may be plotted.  And if a node data file is supplied,
!    then the data can be plotted against the node coordinates without using
!    any finite element structure.)
!
!    Input, character ( len = * ) NODE_DATA_FILE_NAME, the name of the node
!    data file.  If this argument is not supplied, it will be requested.  If
!    the interactive response is blank, then the program will assume that
!    no node data information is to be supplied.  (But the node coordinates
!    will be available and may be plotted.
!    And if an element file is supplied, then the elements can also be
!    displayed.)
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) NODE_DATA_NUM, the number of data items per node.
!
!    Output, real ( kind = 8 ) NODE_COORD(DIM_NUM,NODE_NUM), the coordinates
!    of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    the global index of local node I in element J.
!
!    Output, real ( kind = 8 ) NODE_DATA(NODE_DATA_NUM,NODE_NUM), the data
!    values associated with each node.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num

  character ( len = * ) element_file_name
  integer ( kind = 4 ), dimension (element_order,element_num) :: element_node
  real ( kind = 8 ), dimension (dim_num,node_num) :: node_coord
  character ( len = * ) node_coord_file_name
  real ( kind = 8 ), dimension (node_data_num,node_num) :: node_data
  character ( len = * ) node_data_file_name

  if ( len_trim ( node_coord_file_name ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM_DATA_READ:'
    write ( *, '(a)' ) '  No node coordinate file name was supplied!'
    write ( *, '(a)' ) '  NO DATA WILL BE READ!'
    return
  end if

  call r8mat_data_read ( node_coord_file_name, dim_num, node_num, node_coord )

  if ( 0 < len_trim ( element_file_name ) ) then
    call i4mat_data_read ( element_file_name, element_order, &
      element_num, element_node )
  end if

  if ( 0 < len_trim ( node_data_file_name ) ) then
    call r8mat_data_read ( node_data_file_name, node_data_num, node_num, &
      node_data )
  end if

  return
end
subroutine fem_header_print ( dim_num, node_num, element_num, &
  element_order, node_data_num )

!*****************************************************************************80
!
!! FEM_HEADER_PRINT prints the header to set of FEM files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) NODE_DATA_NUM, the number of data items per node.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension         = ', dim_num
  write ( *, '(a,i8)' ) '  Number of nodes           = ', node_num
  write ( *, '(a,i8)' ) '  Number of elements        = ', element_num
  write ( *, '(a,i8)' ) '  Element order             = ', element_order
  write ( *, '(a,i8)' ) '  Number of node data items = ', node_data_num

  return
end
subroutine fem_header_read ( node_coord_file_name, element_file_name, &
  node_data_file_name, dim_num, node_num, element_num, element_order, &
  node_data_num )

!*****************************************************************************80
!
!! FEM_HEADER_READ reads the sizes of arrays in a set of FEM files.
!
!  Discussion:
!
!    This program reads the node, element and data files that define
!    a finite element geometry and data based on that geometry:
!    * a set of nodes,
!    * a set of elements based on those nodes,
!    * a set of data values associated with each node.
!    and returns the sizes DIM_NUM, NODE_NUM, ELEMENT_NUM, ELEMENT_ORDER,
!    and NODE_DATA_NUM required to allocate space for these arrays.
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
!    coordinate file.  If this argument is not supplied, it will be requested.
!    If the interactive response is blank, or otherwise defective, then the
!    program terminates.
!
!    Input, character ( len = * ) ELEMENT_FILE_NAME, the name of the element
!    file.  If this argument is not supplied, it will be requested.  If the
!    interactive response is blank, then the program will assume that no
!    element information is to be supplied.  (But the node coordinates must
!    be available and may be plotted.  And if a node data file is supplied,
!    then the data can be plotted against the node coordinates without using
!    any finite element structure.)
!
!    Input, character ( len = * ) NODE_DATA_FILE_NAME, the name of the node
!    data file.  If this argument is not supplied, it will be requested.  If
!    the interactive response is blank, then the program will assume that
!    no node data information is to be supplied.  (But the node coordinates
!    will be available and may be plotted.
!    And if an element file is supplied, then the elements can also be
!    displayed.)
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension, inferred from the
!    "shape" of the data in the node file.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes, inferred from the
!    number of lines of data in the node coordinate file.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements, inferred from the
!    number of lines of data in the element file.
!
!    Output, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements, inferred from
!    the number of items in the first line of the element file.
!
!    Output, integer ( kind = 4 ) NODE_DATA_NUM, the number of data items per node,
!    inferred from the number of items in the first line of the node data file.
!
  implicit none

  integer ( kind = 4 ) dim_num
  character ( len = * ) element_file_name
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  character ( len = * ) node_coord_file_name
  character ( len = * ) node_data_file_name
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) node_num2

  if ( len_trim ( node_coord_file_name ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM_HEADER_READ - Error!'
    write ( *, '(a)' ) '  No node coordinate file name was supplied!'
    write ( *, '(a)' ) '  NO DATA WILL BE READ!'
    return
  end if

  if ( len_trim ( element_file_name ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM_HEADER_READ:'
    write ( *, '(a)' ) '  No element file name was supplied.'
    write ( *, '(a)' ) '  Therefore, no element data will be returned.'
  end if

  if ( len_trim ( node_data_file_name ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM_HEADER_READ:'
    write ( *, '(a)' ) '  No node data file name was supplied!'
    write ( *, '(a)' ) '  Therefore, no node data will be returned.'
  end if
!
!  Read the node coordinate file.
!
  call r8mat_header_read ( node_coord_file_name, dim_num, node_num )

  if ( 0 < len_trim ( element_file_name ) ) then
    call i4mat_header_read ( element_file_name, element_order, element_num )
  else
    element_order = 0
    element_num = 0
  end if

  if ( 0 < len_trim ( node_data_file_name ) ) then

    call r8mat_header_read ( node_data_file_name, node_data_num, node_num2 )

    if ( node_num2 /= node_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The number of nodes in the node coordinate '
      write ( *, '(a)' ) '  file is ', node_num, ' but the number of nodes'
      write ( *, '(a)' ) '  in the node data file is ', node_num2
      write ( *, '(a)' ) '  Because of this, no node data will be stored.'
      node_data_num = 0;
    end if

  else
    node_data_num = 0;
  end if

  return
end
subroutine fem_to_tec_handle ( node_coord_file_name, element_file_name, &
  node_data_file_name, tec_file_name )

!*****************************************************************************80
!
!! FEM_TO_TEC_HANDLE reads data from a set of FEM files and writes a TEC file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_COORD_FILE_NAME, the input
!    node coordinate file name.
!
!    Input, character ( len = * ) ELEMENT_FILE_NAME, the input
!    element file name.
!
!    Input, character ( len = * ) NODE_DATA_FILE_NAME, the input
!    node data file name.
!
!    Input, character ( len = * ) TEC_FILE_NAME, the output TEC file name.
!
  implicit none

  integer ( kind = 4 ) dim_num
  character ( len = * ) element_file_name
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_coord
  character ( len = * ) node_coord_file_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_data
  character ( len = * ) node_data_file_name
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num
  character ( len = * ) tec_file_name

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading FEM files:'
  write ( *, '(a)' ) '    "' // trim ( node_coord_file_name ) // '".'
  write ( *, '(a)' ) '    "' // trim ( element_file_name ) // '".'
  write ( *, '(a)' ) '    "' // trim ( node_data_file_name ) // '".'

  call fem_header_read ( node_coord_file_name, element_file_name, &
    node_data_file_name, dim_num, node_num, element_num, element_order, &
    node_data_num )

  call fem_header_print ( dim_num, node_num, element_num, &
    element_order, node_data_num )
!
!  Allocate space for the data, and read the data.
!
  allocate ( node_coord(1:dim_num,1:node_num) )
  allocate ( node_data(1:node_data_num,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )

  call fem_data_read ( node_coord_file_name, element_file_name, &
    node_data_file_name, dim_num, node_num, element_num, element_order, &
    node_data_num, node_coord, element_node, node_data )
!
!  Write the TEC file.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing TEC file: "' // trim ( tec_file_name ) // '".'

  call tec_write ( tec_file_name, dim_num, node_num, element_num, &
  element_order, node_data_num, node_coord, element_node, node_data )

  deallocate ( node_coord )
  deallocate ( node_data )
  deallocate ( element_node )

  return
end
subroutine file_column_count ( file_name, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file are presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines, which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  character ( len = * ) file_name
  logical got_one
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 256 ) line
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( iunit )

    do

      read ( iunit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = iunit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

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
  integer ( kind = 4 ) change
  integer ( kind = 4 ) digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

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
subroutine file_row_count ( file_name, line_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of rows in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Blank lines and comment lines, which begin with '#', are not counted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines found in the file.
!    If the file could not be opened, then LINE_NUM is returned as -1.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line
  integer ( kind = 4 ) line_num
  logical, parameter :: verbose = .false.

  line_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then

    line_num = -1

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file:'
      write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    end if

    return

  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    line_num = line_num + 1

  end do

  close ( unit = iunit )

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
subroutine i4mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_DATA_READ reads data from an I4MAT file.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!    The file may contain more than N points, but this routine
!    will return after reading N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  integer ( kind = 4 ) table(m,n)
  integer ( kind = 4 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_i4vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine i4mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! I4MAT_HEADER_READ reads the header from an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

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
!    Output, integer ( kind = 4 ) ILOC, the location of the last alphabetic
!    character in the string.  If there are no alphabetic
!    characters, ILOC is returned as 0.
!
  implicit none

  logical ch_is_alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iloc
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
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, no alphabetic characters occur in the string.
!    2, all alphabetic characters are "Z" or "z".  S is wrapped around so
!       that all alphabetic characters are "A" or "a".
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) iloc
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
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer ( kind = 4 ) value from a string.
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
!    Output, integer ( kind = 4 ) VALUE, the value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters
!    used to make the value.
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
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an integer ( kind = 4 ) vector from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, integer ( kind = 4 ) IVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
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
!       3 integer ( kind = 4 ) part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer ( kind = 4 ) part,
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
!    07 September 2004
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
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

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

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

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
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
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
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

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
!    Output, integer ( kind = 4 ) WORD_NUM, the number of "words" in the
!    string.  Words are presumed to be separated by one or more blanks.
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
!    Input, integer ( kind = 4 ) TEC_FILE_UNIT, the unit associated with the file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) NODE_DATA_NUM, the number of data items per node.
!
!    Input, real ( kind = 8 ) NODE_COORD(DIM_NUM,NODE_NUM), the coordinates
!    of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    the global index of local node I in element J.
!
!    Input, real ( kind = 8 ) NODE_DATA(NODE_DATA_NUM,NODE_NUM), the
!    data values associated with each node.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  character ( len = 40 ) format_string
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_coord(dim_num,node_num)
  real ( kind = 8 ) node_data(node_data_num,node_num)
  character ( len = * ) tec_file_name
  integer ( kind = 4 ) tec_file_unit
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
!    Input, integer ( kind = 4 ) TEC_FILE_UNIT, the unit associated with the file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) NODE_DATA_NUM, the number of data items per node.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) ierror
  character ( len = 8 ) name
  integer ( kind = 4 ) name_num
  character ( len = 80 ) name_string
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num
  character ( len = * ) tec_file_name
  integer ( kind = 4 ) tec_file_unit
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

  name_string = trim ( name_string ) // '".'

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
!    Output, integer ( kind = 4 ) TEC_FILE_UNIT, the unit on which the file has
!    been opened.  If the file could not be opened, then TEC_FILE_UNIT
!    is returned with the value of -1.
!
  implicit none

  character ( len = * ) tec_file_name
  integer ( kind = 4 ) tec_file_status
  integer ( kind = 4 ) tec_file_unit

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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) NODE_DATA_NUM, the number of data items per node.
!
!    Input, real ( kind = 8 ) NODE_COORD(DIM_NUM,NODE_NUM), the coordinates
!    of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    the global index of local node I in element J.
!
!    Input, real ( kind = 8 ) NODE_DATA(NODE_DATA_NUM,NODE_NUM), the
!    data values associated with each node.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  character ( len = 40 ) format_string
  character ( len = 8 ) name
  integer ( kind = 4 ) name_num
  character ( len = 80 ) name_string
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_coord(dim_num,node_num)
  real ( kind = 8 ) node_data(node_data_num,node_num)
  character ( len = * ) tec_file_name
  integer ( kind = 4 ) tec_file_unit
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
