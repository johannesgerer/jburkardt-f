program main

!*****************************************************************************80
!
!! MESH_TO_ICE reads "ICE" data from a MESH file and writes it to a NETCDF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, string PREFIX, the filename prefix.  The input file is
!    assumed to be "prefix.nc" and the output file will be "prefix.mesh".
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) dim
  integer ( kind = 4 ), allocatable :: edge_label(:)
  integer ( kind = 4 ), allocatable :: edge_vertex(:,:)
  integer ( kind = 4 ) edges
  character ( len = 255 ) filename_mesh
  character ( len = 255 ) filename_nc
  integer ( kind = 4 ), allocatable :: hexahedron_label(:)
  integer ( kind = 4 ), allocatable :: hexahedron_vertex(:,:)
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) iarg
  character ( len = 255 ) prefix
  integer ( kind = 4 ), allocatable :: quadrilateral_label(:)
  integer ( kind = 4 ), allocatable :: quadrilateral_vertex(:,:)
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ), allocatable :: tetrahedron_label(:)
  integer ( kind = 4 ), allocatable :: tetrahedron_vertex(:,:)
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ), allocatable :: triangle_label(:)
  integer ( kind = 4 ), allocatable :: triangle_vertex(:,:)
  integer ( kind = 4 ) triangles
  real ( kind = 8 ), allocatable :: vertex_coordinate(:,:)
  integer ( kind = 4 ), allocatable :: vertex_label(:)
  integer ( kind = 4 ) vertices

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MESH_TO_ICE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read ICE data from a MESH file, write to a NETCDF file.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Check the input argument.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, prefix )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the filename prefix:'
    read ( *, '(a)' ) prefix

  end if
!
!  Create the file names.
!
  filename_nc  = trim ( prefix ) // '.nc'
  filename_mesh = trim ( prefix ) // '.mesh'
!
!  Read sizes;
!
  call mesh_size_read ( filename_mesh, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons )
!
!  Print sizes.
!
  call mesh_size_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons )
!
!  Allocate memory.
!
  allocate ( vertex_coordinate ( dim, vertices ) )
  allocate ( vertex_label ( vertices ) )
  allocate ( edge_vertex ( 2, edges ) )
  allocate ( edge_label ( edges ) )
  allocate ( triangle_vertex ( 3, triangles ) )
  allocate ( triangle_label ( triangles ) )
  allocate ( quadrilateral_vertex ( 4, quadrilaterals ) )
  allocate ( quadrilateral_label ( quadrilaterals ) )
  allocate ( tetrahedron_vertex ( 4, tetrahedrons ) )
  allocate ( tetrahedron_label ( tetrahedrons ) )
  allocate ( hexahedron_vertex ( 8, hexahedrons ) )
  allocate ( hexahedron_label ( hexahedrons ) )
!
!  Read the data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading "' // trim ( filename_mesh ) // '".'

  call mesh_data_read ( filename_mesh, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, &
    triangle_label, quadrilateral_vertex, quadrilateral_label, &
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, &
    hexahedron_label )
!
!  Print the data.
!
  if ( vertices < 250 ) then

    call mesh_data_print ( dim, vertices, edges, triangles, quadrilaterals, &
      tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,&
      edge_vertex, edge_label, triangle_vertex, triangle_label, &
      quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
      tetrahedron_label, hexahedron_vertex, hexahedron_label )

  end if
!
!  Write the data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing "' // trim ( filename_nc ) // '".'

  call ice_write ( filename_nc, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label,  triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )

  write ( *, '(a)' ) '  Conversion completed.'
!
!  Free memory.
!
  deallocate ( vertex_coordinate )
  deallocate ( vertex_label )
  deallocate ( edge_vertex )
  deallocate ( edge_label )
  deallocate ( triangle_vertex )
  deallocate ( triangle_label )
  deallocate ( quadrilateral_vertex )
  deallocate ( quadrilateral_label )
  deallocate ( tetrahedron_vertex )
  deallocate ( tetrahedron_label )
  deallocate ( hexahedron_vertex )
  deallocate ( hexahedron_label )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MESH_TO_ICE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character              ch
  integer   ( kind = 4 ) itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Discussion:
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
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
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
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character              ch
  integer   ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    digit = iachar ( ch ) - 48

  else if ( ch == ' ' ) then

    digit = 0

  else

    digit = -1

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
!    26 October 2008
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
subroutine ice_write ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, &
  edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! ICE_WRITE writes 3D ICE sizes and data to a NETCDF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2010
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
!    Ordinarily, the name should include the extension '.nc'.
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Input, real VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Input, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Input, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Input, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Input, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Input, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for
!    each triangle.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Input, integer ( kind = 4 ) HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
!    that form each hexahedron.
!
!    Input, integer ( kind = 4 ) HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  use netcdf

  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  integer ( kind = 4 ) dim_dimension
  integer ( kind = 4 ) dim_edges
  integer ( kind = 4 ) dim_eight
  integer ( kind = 4 ) dim_four
  integer ( kind = 4 ) dim_hexahedrons
  integer ( kind = 4 ) dim_quadrilaterals
  integer ( kind = 4 ) dim_tetrahedrons
  integer ( kind = 4 ) dim_three
  integer ( kind = 4 ) dim_triangles
  integer ( kind = 4 ) dim_two
  integer ( kind = 4 ) dim_vertices
  integer ( kind = 4 ) dimids(2)
  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges)
  character ( len = * ) filename
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) ncid
  integer ( kind = 4 ) ndims
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) status
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  integer ( kind = 4 ) var_edge_label
  integer ( kind = 4 ) var_edge_vertex
  integer ( kind = 4 ) var_hexahedron_label
  integer ( kind = 4 ) var_hexahedron_vertex
  integer ( kind = 4 ) var_quadrilateral_label
  integer ( kind = 4 ) var_quadrilateral_vertex
  integer ( kind = 4 ) var_tetrahedron_label
  integer ( kind = 4 ) var_tetrahedron_vertex
  integer ( kind = 4 ) var_triangle_label
  integer ( kind = 4 ) var_triangle_vertex
  integer ( kind = 4 ) var_vertex_coordinate
  integer ( kind = 4 ) var_vertex_label
  real ( kind = 8 ) vertex_coordinate(dim,vertices)
  integer ( kind = 4 ) vertex_label(vertices)
  integer ( kind = 4 ) xtype
!
!  Create the file.  This automatically 'opens' it as well.
!
  mode = NF90_CLOBBER
  status = nf90_create ( filename, mode, ncid )
!
!  Put NETCDF into 'define' mode.
!
  status = nf90_redef ( ncid )
!
!  Dimension information.
!
!  If a dimension has length 0, it seems to be taken to be the unlimited
!  dimension (not what you want) and then if you have two such dimensions,
!  you get a complaint that you have tried to define the unlimited dimension
!  twice.  The fix requires the programmer not to write anything whose
!  dimension is zero.
!
  status = nf90_def_dim ( ncid, 'Dimension', dim, dim_dimension )

  status = nf90_def_dim ( ncid, 'Vertices', vertices, dim_vertices )

  if ( 0 < edges ) then
    status = nf90_def_dim ( ncid, 'Edges', edges, dim_edges )
  end if

  if ( 0 < triangles ) then
    status = nf90_def_dim ( ncid, 'Triangles', triangles, dim_triangles )
  end if

  if ( 0 < quadrilaterals ) then
    status = nf90_def_dim ( ncid, 'Quadrilaterals', quadrilaterals, &
      dim_quadrilaterals )
  end if

  if ( 0 < tetrahedrons ) then
    status = nf90_def_dim ( ncid, 'Tetrahedrons', tetrahedrons, &
      dim_tetrahedrons )
  end if

  if ( 0 < hexahedrons ) then
    status = nf90_def_dim ( ncid, 'Hexahedrons', hexahedrons, dim_hexahedrons )
  end if

  status = nf90_def_dim ( ncid, 'Two', 2, dim_two )
  status = nf90_def_dim ( ncid, 'Three', 3, dim_three )
  status = nf90_def_dim ( ncid, 'Four', 4, dim_four )
  status = nf90_def_dim ( ncid, 'Eight', 8, dim_eight )
!
!  Define variables.
!
  ndims = 2
  if ( dim == 2 ) then
    dimids(1) = dim_two
  else if ( dim == 3 ) then
    dimids(1) = dim_three
  end if
  dimids(2) = dim_vertices
  status = nf90_def_var ( ncid, 'Vertex_Coordinate', NF90_DOUBLE, &
    dimids, var_vertex_coordinate )

  ndims = 1
  dimids(1) = dim_vertices
  status = nf90_def_var ( ncid, 'Vertex_Label', NF90_INT, dimids, &
    var_vertex_label )

  if ( 0 < edges ) then
    ndims = 2
    dimids(1) = dim_two
    dimids(2) = dim_edges
    status = nf90_def_var ( ncid, 'Edge_Vertex', NF90_INT, dimids, &
      var_edge_vertex )

    ndims = 1
    dimids(1) = dim_edges
    status = nf90_def_var ( ncid, 'Edge_Label', NF90_INT, dimids, &
      var_edge_label )
  end if

  if ( 0 < triangles ) then
    ndims = 2
    dimids(1) = dim_three
    dimids(2) = dim_triangles
    status = nf90_def_var ( ncid, 'Triangle_Vertex', NF90_INT, dimids, &
      var_triangle_vertex )

    ndims = 1
    dimids(1) = dim_triangles
    status = nf90_def_var ( ncid, 'Triangle_Label', NF90_INT, dimids, &
      var_triangle_label )
  end if

  if ( 0 < quadrilaterals ) then
    ndims = 2
    dimids(1) = dim_four
    dimids(2) = dim_quadrilaterals
    status = nf90_def_var ( ncid, 'Quadrilateral_Vertex', NF90_INT, &
      dimids, var_quadrilateral_vertex )

    ndims = 1
    dimids(1) = dim_quadrilaterals
    status = nf90_def_var ( ncid, 'Quadrilateral_Label', NF90_INT, &
      dimids, var_quadrilateral_label )
  end if

  if ( 0 < tetrahedrons ) then
    ndims = 2
    dimids(1) = dim_four
    dimids(2) = dim_tetrahedrons
    status = nf90_def_var ( ncid, 'Tetrahedron_Vertex', NF90_INT, &
      dimids, var_tetrahedron_vertex )

    ndims = 1
    dimids(1) = dim_tetrahedrons
    status = nf90_def_var ( ncid, 'Tetrahedron_Label', NF90_INT, &
      dimids, var_tetrahedron_label )
  end if

  if ( 0 < hexahedrons ) then
    ndims = 2
    dimids(1) = dim_eight
    dimids(2) = dim_hexahedrons
    status = nf90_def_var ( ncid, 'Hexahedron_Vertex', NF90_INT, &
      dimids, var_hexahedron_vertex )

    ndims = 1
    dimids(1) = dim_hexahedrons
    status = nf90_def_var ( ncid, 'Hexahedron_Label', NF90_INT, dimids, &
      var_hexahedron_label )
  end if
!
!  Terminate the definition phase.
!
  status = nf90_enddef ( ncid )
!
!  Write the data.
!
  status = nf90_put_var ( ncid, var_vertex_coordinate, vertex_coordinate )
  status = nf90_put_var ( ncid, var_vertex_label, vertex_label )

  if ( 0 < edges ) then
    status = nf90_put_var ( ncid, var_edge_vertex, edge_vertex )
    status = nf90_put_var ( ncid, var_edge_label, edge_label )
  end if

  if ( 0 < triangles ) then
    status = nf90_put_var ( ncid, var_triangle_vertex, triangle_vertex )
    status = nf90_put_var ( ncid, var_triangle_label, triangle_label )
  end if

  if ( 0 < quadrilaterals ) then
    status = nf90_put_var ( ncid, var_quadrilateral_vertex, &
      quadrilateral_vertex )
    status = nf90_put_var ( ncid, var_quadrilateral_label, quadrilateral_label )
  end if

  if ( 0 < tetrahedrons ) then
    status = nf90_put_var ( ncid, var_tetrahedron_vertex, tetrahedron_vertex )
    status = nf90_put_var ( ncid, var_tetrahedron_label, tetrahedron_label )
  end if

  if ( 0 < hexahedrons ) then
    status = nf90_put_var ( ncid, var_hexahedron_vertex, hexahedron_vertex )
    status = nf90_put_var ( ncid, var_hexahedron_label, hexahedron_label )
  end if
!
!  Close the file.
!
  status = nf90_close ( ncid )

  return
end
subroutine mesh_data_print ( dim, vertices, edges, triangles, quadrilaterals, &
  tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,  edge_vertex, &
  edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
  quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
  hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! MESH_DATA_PRINT prints mesh data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2010
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
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Input, real VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Input, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Input, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Input, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Input, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Input, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for each
!    triangle.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
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
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  real    ( kind = 8 ) vertex_coordinate(dim,vertices)
  integer ( kind = 4 ) vertex_label(vertices)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices:'
  write ( *, '(a)' ) ' '
  if ( dim == 2 ) then
    do j = 1, vertices
      write ( *, '(2(2x,f10.4),2x,''('',i4,'')'')' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  else if ( dim == 3 ) then
    do j = 1, vertices
      write ( *, '(3(2x,f10.4),2x,''('',i4,'')'')' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  end if

  if ( 0 < edges ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Edges:'
    write ( *, '(a)' ) ' '
    do j = 1, edges
      write ( *, '(2(2x,i8),2x,''('',i4,'')'')' ) &
        edge_vertex(1:2,j), edge_label(j)
    end do
  end if

  if ( 0 < triangles ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Triangles:'
    write ( *, '(a)' ) ' '
    do j = 1, triangles
      write ( *, '(3(2x,i8),2x,''('',i4,'')'')' ) &
        triangle_vertex(1:3,j), triangle_label(j)
    end do
  end if

  if ( 0 < quadrilaterals ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Quadrilaterals:'
    write ( *, '(a)' ) ' '
    do j = 1, quadrilaterals
      write ( *, '(4(2x,i8),2x,''('',i4,'')'')' ) &
        quadrilateral_vertex(1:4,j), quadrilateral_label(j)
    end do
  end if

  if ( 0 < tetrahedrons ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Tetrahedrons:'
    write ( *, '(a)' ) ' '
    do j = 1, tetrahedrons
      write ( *, '(4(2x,i8),2x,''('',i4,'')'')' ) &
        tetrahedron_vertex(1:4,j), tetrahedron_label(j)
    end do
  end if

  if ( 0 < hexahedrons ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Hexahedrons:'
    write ( *, '(a)' ) ' '
    do j = 1, hexahedrons
      write ( *, '(8(2x,i8),2x,''('',i4,'')'')' ) &
        hexahedron_vertex(1:8,j), hexahedron_label(j)
    end do
  end if

  return
end
subroutine mesh_data_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! MESH_DATA_READ reads data from a MESH file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2010
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
!    Input, character ( len = * ) FILENAME, the name of the MESH file.
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRAONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Output, real VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Output, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Output, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Output, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Output, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Output, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for each
!    triangle.
!
!    Output, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Output, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Output, integer ( kind = 4 ) HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
!    that form each hexahedron.
!
!    Output, integer ( kind = 4 ) HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges)
  character ( len = * ) filename
  integer ( kind = 4 ) fileunit
  integer ( kind = 4 ) hexahedron
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) i4vec(9)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  character ( len = 80 ) keyword
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) quadrilateral
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  real    ( kind = 8 ) r8vec(9)
  logical s_begin
  logical s_eqi
  integer ( kind = 4 ) tetrahedron
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  character ( len = 255 ) text
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  integer ( kind = 4 ) vertex
  real    ( kind = 8 ) vertex_coordinate(dim,vertices)
  integer ( kind = 4 ) vertex_label(vertices)
!
!  Initialize everything to nothing.
!
  vertex_coordinate(1:dim,1:vertices) = 0.0D+00
  vertex_label(1:vertices) = 0
  edge_vertex(1:2,1:edges) = 0
  edge_label(1:edges) = 0
  triangle_vertex(1:3,1:triangles) = 0
  triangle_label(1:triangles) = 0
  quadrilateral_vertex(1:4,1:quadrilaterals) = 0
  quadrilateral_label(1:quadrilaterals) = 0
  tetrahedron_vertex(1:4,1:tetrahedrons) = 0
  tetrahedron_label(1:tetrahedrons) = 0
  hexahedron_vertex(1:8,1:hexahedrons) = 0
  hexahedron_label(1:hexahedrons) = 0
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if
!
!  Read lines til you get alphanumerics and determine a "mode"
!
  line_num = 0
  keyword = 'NONE'

  do

    read ( fileunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( len_trim ( text ) == 0 ) then
      keyword = 'NONE'
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if
!
!  Remove initial blanks.
!
    text = adjustl ( text )
!
!  Expecting a keyword.
!
        if ( s_eqi ( text, 'CORNERS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'DIMENSION' ) ) then

      keyword = 'DIMENSION'

    else if ( s_eqi ( text, 'EDGES' ) ) then

      keyword = 'EDGES'

    else if ( s_eqi ( text, 'END' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  END statement encountered.'
      exit

    else if ( s_eqi ( text, 'HEXAHEDRA' ) .or. &
              s_eqi ( text, 'HEXAHEDRONS' ) ) then

      keyword = 'HEXAHEDRONS'

    else if ( s_begin ( text, 'MESHVERSIONFORMATTED' ) ) then

    else if ( s_eqi ( text, 'NORMALATQUADRILATERALVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATTRIANGLEVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'QUADRILATERALS' ) ) then

      keyword = 'QUADRILATERALS'

    else if ( s_eqi ( text, 'REQUIREDEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'REQUIREDVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'RIDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTATEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TETRAHEDRA' ) .or. &
              s_eqi ( text, 'TETRAHEDRONS' ) ) then

      keyword = 'TETRAHEDRONS'

    else if ( s_eqi ( text, 'TRIANGLES' ) ) then

      keyword = 'TRIANGLES'

    else if ( s_eqi ( text, 'VERTICES' ) ) then

      keyword = 'VERTICES'
!
!  Presumably, numeric data to be processed by keyword.
!
    else if ( s_eqi ( keyword, 'DIMENSION' ) ) then

      call s_to_i4 ( text, dim, ierror, length )

      keyword = 'NONE'

    else if ( s_eqi ( keyword, 'EDGES' ) ) then

      call s_to_i4 ( text, edges, ierror, length )

      keyword = 'EDGE_VERTEX'
      edge = 0

    else if ( s_eqi ( keyword, 'EDGE_VERTEX' ) ) then

      call s_to_i4vec ( text, 3, i4vec, ierror )
      edge = edge + 1
      edge_vertex(1:2,edge) = i4vec(1:2)
      edge_label(edge) = i4vec(3)

    else if ( s_eqi ( keyword, 'HEXAHEDRONS' ) ) then

      call s_to_i4 ( text, hexahedrons, ierror, length )

      keyword = 'HEXAHEDRON_VERTEX'
      hexahedron = 0

    else if ( s_eqi ( keyword, 'HEXAHEDRON_VERTEX' ) ) then

      call s_to_i4vec ( text, 9, i4vec, ierror )
      hexahedron = hexahedron + 1
      hexahedron_vertex(1:8,hexahedron) = i4vec(1:8)
      hexahedron_label(hexahedron) = i4vec(9)

    else if ( s_eqi ( keyword, 'QUADRILATERALS' ) ) then

      call s_to_i4 ( text, quadrilaterals, ierror, length )

      keyword = 'QUADRILATERAL_VERTEX'
      quadrilateral = 0

    else if ( s_eqi ( keyword, 'QUADRILATERAL_VERTEX' ) ) then

      call s_to_i4vec ( text, 5, i4vec, ierror )
      quadrilateral = quadrilateral + 1
      quadrilateral_vertex(1:4,quadrilateral) = i4vec(1:4)
      quadrilateral_label(quadrilateral) = i4vec(5)

    else if ( s_eqi ( keyword, 'TETRAHEDRONS' ) ) then

      call s_to_i4 ( text, tetrahedrons, ierror, length )

      keyword = 'TETRAHEDRON_VERTEX'
      tetrahedron = 0

    else if ( s_eqi ( keyword, 'TETRAHEDRON_VERTEX' ) ) then

      call s_to_i4vec ( text, 5, i4vec, ierror )
      tetrahedron = tetrahedron + 1
      tetrahedron_vertex(1:4,tetrahedron) = i4vec(1:4)
      tetrahedron_label(tetrahedron) = i4vec(5)

    else if ( s_eqi ( keyword, 'TRIANGLES' ) ) then

      call s_to_i4 ( text, triangles, ierror, length )

      keyword = 'TRIANGLE_VERTEX'
      triangle = 0

    else if ( s_eqi ( keyword, 'TRIANGLE_VERTEX' ) ) then

      call s_to_i4vec ( text, 4, i4vec, ierror )
      triangle = triangle + 1
      triangle_vertex(1:3,triangle) = i4vec(1:3)
      triangle_label(triangle) = i4vec(4)

    else if ( s_eqi ( keyword, 'VERTICES' ) ) then

      call s_to_i4 ( text, vertices, ierror, length )

      keyword = 'VERTEX_COORDINATE'
      vertex = 0

    else if ( s_eqi ( keyword, 'VERTEX_COORDINATE' ) ) then

      call s_to_r8vec ( text, dim + 1, r8vec, ierror )
      vertex = vertex + 1
      vertex_coordinate(1:dim,vertex) = r8vec(1:dim)
      vertex_label(vertex) = int ( r8vec(dim+1) )

    else if ( s_eqi ( keyword, 'SKIP' ) ) then

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_DATA_READ - Fatal error!'
      write ( *, '(a,i8)' ) &
        '  Could not find keyword while reading line ', line_num
      write ( *, '(a)' ) '"' // trim ( text ) // '".'
      stop

    end if

  end do
!
!  Close the file.
!
  close ( unit = fileunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Read ', line_num, &
    ' lines from "' // trim ( filename ) // '".'

  return
end
subroutine mesh_size_print ( dim, vertices, edges, triangles, quadrilaterals, &
  tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! MESH_SIZE_PRINT prints mesh sizes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2010
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
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of dimensions = ', dim
  write ( *, '(a,i8)' ) '  Number of vertices = ', vertices
  write ( *, '(a,i8)' ) '  Number of edges = ', edges
  write ( *, '(a,i8)' ) '  Number of triangles = ', triangles
  write ( *, '(a,i8)' ) '  Number of quadrilaterals = ', quadrilaterals
  write ( *, '(a,i8)' ) '  Number of tetrahedrons = ', tetrahedrons
  write ( *, '(a,i8)' ) '  Number of hexahedrons = ', hexahedrons

  return
end
subroutine mesh_size_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! MESH_SIZE_READ reads sizes from a MESH file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2010
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
!    Input, character ( len = * ) FILENAME, the name of the MESH file.
!
!    Output, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Output, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Output, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Output, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Output, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Output, integer ( kind = 4 ) TETRAHEDRAONS, the number of tetrahedrons
!    (may be 0).
!
!    Output, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  character ( len = * ) filename
  integer ( kind = 4 ) fileunit
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 80 ) keyword
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) quadrilaterals
  logical s_begin
  logical s_eqi
  integer ( kind = 4 ) tetrahedrons
  character ( len = 255 ) text
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices
!
!  Initialize everything to nothing.
!
  dim = 0
  vertices = 0
  edges = 0
  triangles = 0
  quadrilaterals = 0
  tetrahedrons = 0
  hexahedrons = 0
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_SIZE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if
!
!  Read lines til you get alphanumerics and determine a "mode"
!
  line_num = 0
  keyword = 'NONE'

  do

    read ( fileunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( len_trim ( text ) == 0 ) then
      keyword = 'NONE'
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if
!
!  Remove initial blanks.
!
    text = adjustl ( text )
!
!  Expecting a keyword.
!
        if ( s_eqi ( text, 'CORNERS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'DIMENSION' ) ) then

      keyword = 'DIMENSION'

    else if ( s_eqi ( text, 'EDGES' ) ) then

      keyword = 'EDGES'

    else if ( s_eqi ( text, 'END' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  END statement encountered.'
      exit

    else if ( s_eqi ( text, 'HEXAHEDRA' ) .or. &
              s_eqi ( text, 'HEXAHEDRONS' ) ) then

      keyword = 'HEXAHEDRONS'

    else if ( s_begin ( text, 'MESHVERSIONFORMATTED' ) ) then

    else if ( s_eqi ( text, 'NORMALATQUADRILATERALVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATTRIANGLEVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'QUADRILATERALS' ) ) then

      keyword = 'QUADRILATERALS'

    else if ( s_eqi ( text, 'REQUIREDEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'REQUIREDVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'RIDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTATEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TETRAHEDRA' ) .or. &
              s_eqi ( text, 'TETRAHEDRONS' ) ) then

      keyword = 'TETRAHEDRONS'

    else if ( s_eqi ( text, 'TRIANGLES' ) ) then

      keyword = 'TRIANGLES'

    else if ( s_eqi ( text, 'VERTICES' ) ) then

      keyword = 'VERTICES'
!
!  Presumably, numeric data to be processed by keyword.
!
    else if ( s_eqi ( keyword, 'DIMENSION' ) ) then

      call s_to_i4 ( text, dim, ierror, length )

      keyword = 'NONE'

    else if ( s_eqi ( keyword, 'EDGES' ) ) then

      call s_to_i4 ( text, edges, ierror, length )

      keyword = 'EDGE_VERTEX'

    else if ( s_eqi ( keyword, 'EDGE_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'HEXAHEDRONS' ) ) then

      call s_to_i4 ( text, hexahedrons, ierror, length )

      keyword = 'HEXAHEDRON_VERTEX'

    else if ( s_eqi ( keyword, 'HEXAHEDRON_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'QUADRILATERALS' ) ) then

      call s_to_i4 ( text, quadrilaterals, ierror, length )

      keyword = 'QUADRILATERAL_VERTEX'

    else if ( s_eqi ( keyword, 'QUADRILATERAL_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'TETRAHEDRONS' ) ) then

      call s_to_i4 ( text, tetrahedrons, ierror, length )

      keyword = 'TETRAHEDRON_VERTEX'

    else if ( s_eqi ( keyword, 'TETRAHEDRON_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'TRIANGLES' ) ) then

      call s_to_i4 ( text, triangles, ierror, length )

      keyword = 'TRIANGLE_VERTEX'

    else if ( s_eqi ( keyword, 'TRIANGLE_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'VERTICES' ) ) then

      call s_to_i4 ( text, vertices, ierror, length )

      keyword = 'VERTEX_COORDINATE'

    else if ( s_eqi ( keyword, 'VERTEX_COORDINATE' ) ) then

    else if ( s_eqi ( keyword, 'SKIP' ) ) then

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_SIZE_READ - Fatal error!'
      write ( *, '(a,i8)' ) &
        '  Could not find keyword while reading line ', line_num
      write ( *, '(a)' ) '"' // trim ( text ) // '".'
      stop

    end if

  end do
!
!  Close the file.
!
  close ( unit = fileunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Read ', line_num, &
    ' lines from "' // trim ( filename ) // '".'

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

  logical                ch_eqi
  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) i2
  logical                s_begin
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  integer   ( kind = 4 ) s2_length

  s1_length = len_trim ( s1 )
  s2_length = len_trim ( s2 )
!
!  If either string is blank, then both must be blank to match.
!  Otherwise, a blank string matches anything, which is not
!  what most people want.
!
  if ( s1_length == 0 .or. s2_length == 0 ) then

    if ( s1_length == 0 .and. s2_length == 0 ) then
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

      if ( s1_length < i1 ) then
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

      if ( s2_length < i2 ) then
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
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Discussion:
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

  character              c1
  character              c2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lenc
  logical                s_eqi
  character ( len = *  ) s1
  integer   ( kind = 4 ) s1_length
  character ( len = *  ) s2
  integer   ( kind = 4 ) s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )

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

  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

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
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) length
  character ( len = * )  s
  integer   ( kind = 4 ) state
  character              :: TAB = achar ( 9 )
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
subroutine s_to_i4vec ( s, n, i4vec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an integer vector from a string.
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
!    Output, integer ( kind = 4 ) I4VEC(N), the values read from the string.
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
  integer   ( kind = 4 ) i4vec(n)
  integer   ( kind = 4 ) length
  character ( len = * )  s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), i4vec(i), ierror, length )

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
  real      ( kind = 8 ) dval
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) jbot
  integer   ( kind = 4 ) jsgn
  integer   ( kind = 4 ) jtop
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) ndig
  real      ( kind = 8 ) rbot
  real      ( kind = 8 ) rexp
  real      ( kind = 8 ) rtop
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
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

    if ( s_length < length + 1 ) then
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
  if ( iterm /= 1 .and. length + 1 == s_length ) then
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
subroutine s_to_r8vec ( s, n, r8vec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Discussion:
!
!    An R8VEC is a vector of real values, of type "real ( kind = 8 )".
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
!    Output, real ( kind = 8 ) R8VEC(N), the values read from the string.
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
  real      ( kind = 8 ) r8vec(n)
  character ( len = * )  s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), r8vec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
