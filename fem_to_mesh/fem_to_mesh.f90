program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM_TO_MESH.
!
!  Discussion:
!
!    FEM_TO_MESH converts data from a FEM format to a MESH format.
!
!    The FEM format defines "node", "element", and "boundary_node_mask",
!    files for a mesh.  A typical set of such files might have the names
!    "suv_nodes.txt", "suv_elements.txt" and "suv_boundary_node_mask.txt".
!
!    This program reads these files and creates a MESH file, whose
!    name might be "suv.mesh".
!
!  Usage:
!
!    fem_to_mesh prefix
!
!    reads the FEM files
!      "prefix"_nodes.txt
!      "prefix"_elements.txt and
!      "prefix"_boundary_node_mask.txt
!    and creates the MESH file
!      "prefix".mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) boundary_node_mask_filename
  character ( len = 255 ) element_filename
  integer   ( kind = 4 )  iarg
  integer   ( kind = 4 )  iargc
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  ilen
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  ipxfargc
  character ( len = 255 ) mesh_filename
  character ( len = 255 ) node_filename
  integer   ( kind = 4 )  num_arg
  character ( len = 255 ) prefix

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_TO_MESH'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a set of FEM files;'
  write ( *, '(a)' ) '  Write a corresponding MESH file.'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the filename prefix:'
    read ( *, '(a)', iostat = ios ) prefix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM_TO_MESH - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, prefix )

  end if
!
!  Set the filenames.
!
  node_filename    = trim ( prefix ) // '_nodes.txt'
  element_filename = trim ( prefix ) // '_elements.txt'
  boundary_node_mask_filename    = trim ( prefix ) // '_boundary_node_mask.txt'
  mesh_filename    = trim ( prefix ) // '.mesh'
!
!  Now we know what to do.
!
  call fem_to_mesh_handle ( node_filename, element_filename, &
    boundary_node_mask_filename, mesh_filename )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_TO_MESH'
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
  integer   ( kind = 4 ) itemp

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
  logical   ch_eqi

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
  integer   ( kind = 4 ) digit

  if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine fem_data_read ( node_filename, element_filename, &
  boundary_node_mask_filename, dim_num, node_num, element_num, element_order, &
  node_coord, element_node, boundary_node_mask )

!*****************************************************************************80
!
!! FEM_DATA_READ reads data from a set of FEM files.
!
!  Discussion:
!
!    This program reads the node, element and mask files that define
!    a finite element geometry and data based on that geometry:
!    * a set of nodes,
!    * a set of elements based on those nodes,
!    * a mask that indicates which nodes are on the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILENAME, the name of the node
!    coordinate file.  If this argument is not supplied, it will be requested.
!    If the interactive response is blank, or otherwise defective, then the
!    program terminates.
!
!    Input, character ( len = * ) ELEMENT_FILENAME, the name of the element
!    file.  If this argument is not supplied, it will be requested.  If the
!    interactive response is blank, then the program will assume that no
!    element information is to be supplied.  (But the node coordinates must
!    be available and may be plotted.  And if a node data file is supplied,
!    then the data can be plotted against the node coordinates without using
!    any finite element structure.)
!
!    Input, character ( len = * ) BOUNDARY_NODE_MASK_FILENAME, the name of the
!    boundary node mask file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Output, real ( kind = 8 ) NODE_COORD(DIM_NUM,NODE_NUM), the coordinates
!    of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    the global index of local node I in element J.
!
!    Output, integer ( kind = 4 ) BOUNDARY_NODE_MASK(NODE_NUM), is 0 for
!    each interior node and 1 for each boundary node.
!
  implicit none

  integer   ( kind = 4 ) dim_num
  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) element_order
  integer   ( kind = 4 ) node_num

  integer   ( kind = 4 ) boundary_node_mask(node_num)
  character ( len = * )  boundary_node_mask_filename
  character ( len = * )  element_filename
  integer   ( kind = 4 ) element_node(element_order,element_num)
  real ( kind = 8 ) node_coord(dim_num,node_num)
  character ( len = * )  node_filename

  call r8mat_data_read ( node_filename, dim_num, node_num, node_coord )

  call i4mat_data_read ( element_filename, element_order, &
    element_num, element_node )

  call i4mat_data_read ( boundary_node_mask_filename, 1, node_num, &
    boundary_node_mask )

  return
end
subroutine fem_header_print ( dim_num, node_num, element_num, element_order )

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
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension         = ', dim_num
  write ( *, '(a,i8)' ) '  Number of nodes           = ', node_num
  write ( *, '(a,i8)' ) '  Number of elements        = ', element_num
  write ( *, '(a,i8)' ) '  Element order             = ', element_order

  return
end
subroutine fem_header_read ( node_filename, element_filename, dim_num, &
  node_num, element_num, element_order )

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
!    and returns the sizes DIM_NUM, NODE_NUM, ELEMENT_NUM, ELEMENT_ORDER,
!    required to allocate space for these arrays.
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
!    Input, character ( len = * ) NODE_FILENAME, the name of the node
!    coordinate file.  If this argument is not supplied, it will be requested.
!    If the interactive response is blank, or otherwise defective, then the
!    program terminates.
!
!    Input, character ( len = * ) ELEMENT_FILENAME, the name of the element
!    file.  If this argument is not supplied, it will be requested.  If the
!    interactive response is blank, then the program will assume that no
!    element information is to be supplied.  (But the node coordinates must
!    be available and may be plotted.  And if a node data file is supplied,
!    then the data can be plotted against the node coordinates without using
!    any finite element structure.)
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension, inferred from
!    the "shape" of the data in the node file.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes, inferred from
!    the number of lines of data in the node coordinate file.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements, inferred
!    from the number of lines of data in the element file.
!
!    Output, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements,
!    inferred from the number of items in the first line of the element file.
!
  implicit none

  integer   ( kind = 4 ) dim_num
  character ( len = * ) element_filename
  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) element_order
  character ( len = * ) node_filename
  integer   ( kind = 4 ) node_num
  integer   ( kind = 4 ) node_num2

  call r8mat_header_read ( node_filename, dim_num, node_num )

  call i4mat_header_read ( element_filename, element_order, element_num )

  return
end
subroutine fem_to_mesh_handle ( node_filename, element_filename, &
  boundary_node_mask_filename, mesh_filename )

!*****************************************************************************80
!
!! FEM_TO_MESH_HANDLE copies data from a FEM dataset to a MESH data set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILENAME, the FEM node filename.
!
!    Input, character ( len = * ) ELEMENT_FILENAME, the FEM element filename.
!
!    Input, character ( len = * ) BOUNDARY_NODE_MASK_FILENAME, the FEM
!    boundary node mask filename.
!
!    Input, character ( len = * ) MESH_FILENAME, the MESH filename.
!
  implicit none

  integer   ( kind = 4 ), allocatable :: boundary_node_mask(:)
  character ( len = * ) boundary_node_mask_filename
  integer   ( kind = 4 ) dim
  integer   ( kind = 4 ) dim_num
  integer   ( kind = 4 ), allocatable :: edge_label(:)
  integer   ( kind = 4 ), allocatable :: edge_vertex(:,:)
  integer   ( kind = 4 ) edges
  character ( len = * ) element_filename
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) element_order
  integer   ( kind = 4 ), allocatable :: hexahedron_label(:)
  integer   ( kind = 4 ), allocatable :: hexahedron_vertex(:,:)
  integer   ( kind = 4 ) hexahedrons
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  character ( len = * )  mesh_filename
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_coord
  character ( len = * ) node_filename
  integer   ( kind = 4 ) node_num
  integer   ( kind = 4 ), allocatable :: quadrilateral_label(:)
  integer   ( kind = 4 ), allocatable :: quadrilateral_vertex(:,:)
  integer   ( kind = 4 ) quadrilaterals
  integer   ( kind = 4 ), allocatable :: tetrahedron_label(:)
  integer   ( kind = 4 ), allocatable :: tetrahedron_vertex(:,:)
  integer   ( kind = 4 ) tetrahedrons
  integer   ( kind = 4 ), allocatable :: triangle_label(:)
  integer   ( kind = 4 ), allocatable :: triangle_vertex(:,:)
  integer   ( kind = 4 ) triangles
  real ( kind = 8 ), allocatable :: vertex_coordinate(:,:)
  integer   ( kind = 4 ), allocatable :: vertex_label(:)
  integer   ( kind = 4 ) vertices

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading FEM files:'
  write ( *, '(a)' ) '    "' // trim ( node_filename ) // '".'
  write ( *, '(a)' ) '    "' // trim ( element_filename ) // '".'
  write ( *, '(a)' ) '    "' // trim ( boundary_node_mask_filename ) // '".'

  call fem_header_read ( node_filename, element_filename, &
    dim_num, node_num, element_num, element_order )

  call fem_header_print ( dim_num, node_num, element_num, element_order )
!
!  Allocate space for the data, and read the data.
!
  allocate ( node_coord(1:dim_num,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( boundary_node_mask(1:node_num) )

  call fem_data_read ( node_filename, element_filename, &
    boundary_node_mask_filename, dim_num, node_num, element_num, &
    element_order, node_coord, element_node, boundary_node_mask )
!
!  Set up the MESH data.
!
  dim = dim_num
  vertices = node_num
  edges = 0
  if ( dim == 2 ) then
    triangles = element_num
  else
    triangles = 0
  end if
  quadrilaterals = 0
  if ( dim == 3 ) then
    tetrahedrons = element_num
  else
    tetrahedrons = 0
  end if
  hexahedrons = 0

  allocate ( edge_label(edges) )
  allocate ( edge_vertex(2,edges) )
  allocate ( hexahedron_label(hexahedrons) )
  allocate ( hexahedron_vertex(8,hexahedrons) )
  allocate ( quadrilateral_label(quadrilaterals) )
  allocate ( quadrilateral_vertex(4,quadrilaterals) )
  allocate ( tetrahedron_label(tetrahedrons) )
  allocate ( tetrahedron_vertex(4,tetrahedrons) )
  allocate ( triangle_label(triangles) )
  allocate ( triangle_vertex(3,triangles) )
  allocate ( vertex_coordinate(dim,vertices) )
  allocate ( vertex_label(vertices) )

  triangle_vertex(1:3,1:triangles) = element_node(1:3,1:triangles)
  triangle_label(1:triangles) = 0

  tetrahedron_vertex(1:4,1:tetrahedrons) = element_node(1:4,1:tetrahedrons)
  tetrahedron_label(1:tetrahedrons) = 0

  vertex_coordinate(1:dim,1:vertices) = node_coord(1:dim,1:vertices)
  vertex_label(1:vertices) = boundary_node_mask(1:vertices)
!
!  Write the MESH file.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing MESH file: "' // trim ( mesh_filename ) // '".'

  call mesh_write ( mesh_filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )
!
!  Deallocate memory.
!
  deallocate ( edge_label )
  deallocate ( edge_vertex )
  deallocate ( hexahedron_label )
  deallocate ( hexahedron_vertex )
  deallocate ( quadrilateral_label )
  deallocate ( quadrilateral_vertex )
  deallocate ( tetrahedron_label )
  deallocate ( tetrahedron_vertex )
  deallocate ( triangle_label )
  deallocate ( triangle_vertex )
  deallocate ( vertex_coordinate )
  deallocate ( vertex_label )

  deallocate ( boundary_node_mask )
  deallocate ( node_coord )
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

  integer   ( kind = 4 )  column_num
  character ( len = * )   file_name
  logical                 got_one
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  iunit
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

  character              c
  integer   ( kind = 4 ) change
  integer   ( kind = 4 ) digit
  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lens

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
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines found in the
!    file.  If the file could not be opened, then LINE_NUM is returned as -1.
!
  implicit none

  character ( len = * )   file_name
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  iunit
  character ( len = 256 ) line
  integer   ( kind = 4 )  line_num
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
  logical              lopen

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

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 255 )  line
  integer   ( kind = 4 )   table(m,n)
  integer   ( kind = 4 )   x(m)

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

  character ( len = * )  input_filename
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

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
subroutine mesh_write ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
  vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! MESH_WRITE writes sizes and data to a MESH file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to be created.
!    Ordinarily, the name should include the extension ".mesh".
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should
!    be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, real ( kind = 8 ) VERTEX_COORDINATE(DIM,VERTICES), the coordinate
!    of each vertex.
!
!    Input, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Input, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Input, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for each
!    triangle.
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
!    that form each hexahedron.
!
!    Input, integer ( kind = 4 ) HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges)
  character ( len = * ) filename
  integer ( kind = 4 ) fileunit
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  real ( kind = 8 ) vertex_coordinate(dim,vertices)
  integer ( kind = 4 ) vertex_label(vertices)
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if

  write ( fileunit, '(a)' ) 'MeshVersionFormatted 1'
  write ( fileunit, '(a)' ) '#  Created by mesh_write.f90'
!
!  Dimension information.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'Dimension'
  write ( fileunit, '(i8)' ) dim
!
!  Vertices.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'Vertices'
  write ( fileunit, '(i8)' ) vertices
  if ( dim == 2 ) then
    do j = 1, vertices
      write ( fileunit, '(2(2x,f10.6),2x,i8)' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  else if ( dim == 3 ) then
    do j = 1, vertices
      write ( fileunit, '(3(2x,f10.6),2x,i8)' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  end if
!
!  Edges.
!
  if ( 0 < edges ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Edges'
    write ( fileunit, '(i8)' ) edges
    do j = 1, edges
      write ( fileunit, '(2(2x,i8),2x,i8)' ) &
        edge_vertex(1:2,j), edge_label(j)
    end do
  end if
!
!  Triangles.
!
  if ( 0 < triangles ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Triangles'
    write ( fileunit, '(i8)' ) triangles
    do j = 1, triangles
      write ( fileunit, '(3(2x,i8),2x,i8)' ) &
        triangle_vertex(1:3,j), triangle_label(j)
    end do
  end if
!
!  Quadrilaterals.
!
  if ( 0 < quadrilaterals ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Quadrilaterals'
    write ( fileunit, '(i8)' ) quadrilaterals
    do j = 1, quadrilaterals
      write ( fileunit, '(4(2x,i8),2x,i8)' ) &
        quadrilateral_vertex(1:4,j), quadrilateral_label(j)
    end do
  end if
!
!  Tetrahedron.
!
  if ( 0 < tetrahedrons ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Tetrahedra'
    write ( fileunit, '(i8)' ) tetrahedrons
    do j = 1, tetrahedrons
      write ( fileunit, '(4(2x,i8),2x,i8)' ) &
        tetrahedron_vertex(1:4,j), tetrahedron_label(j)
    end do
  end if
!
!  Hexahedron.
!
  if ( 0 < hexahedrons ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Hexahedra'
    write ( fileunit, '(i8)' ) hexahedrons
    do j = 1, hexahedrons
      write ( fileunit, '(8(2x,i8),2x,i8)' ) &
        hexahedron_vertex(1:8,j), hexahedron_label(j)
    end do
  end if
!
!  End.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'End'

  close ( unit = fileunit )

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

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 255 )  line
  real ( kind = 8 )   table(m,n)
  real ( kind = 8 )   x(m)

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

  character ( len = * )  input_filename
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

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

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) length
  character ( len = * )  s
  integer   ( kind = 4 ) state
  integer   ( kind = 4 ) value

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

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) ivec(n)
  integer   ( kind = 4 ) length
  character ( len = * )  s

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

  character              c
  logical                ch_eqi
  real ( kind = 8 ) dval
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) jbot
  integer   ( kind = 4 ) jsgn
  integer   ( kind = 4 ) jtop
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) nchar
  integer   ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * )  s

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

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * )  s

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

  logical                blank
  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_len
  integer   ( kind = 4 ) word_num

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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
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
