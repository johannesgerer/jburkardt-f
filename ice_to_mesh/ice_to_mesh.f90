program main

!*****************************************************************************80
!
!! ICE_TO_MESH reads ICE data from a NETCDF file and writes to a MESH file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  write ( *, '(a)' ) 'ICE_TO_MESH:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read ICE data from NETCDF file, write to MESH file.'
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
  call size_read ( filename_nc, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons )
!
!  Print sizes.
!
  call size_print ( dim, vertices, edges, triangles, quadrilaterals, &
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
  write ( *, '(a)' ) '  Reading "' // trim ( filename_nc ) // '".'

  call data_read ( filename_nc, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, &
    triangle_label, quadrilateral_vertex, quadrilateral_label, &
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, &
    hexahedron_label )
!
!  Print the data.
!
  if ( vertices < 250 ) then

    call data_print ( dim, vertices, edges, triangles, quadrilaterals, &
      tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,&
      edge_vertex, edge_label, triangle_vertex, triangle_label, &
      quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
      tetrahedron_label, hexahedron_vertex, hexahedron_label )

  end if
!
!  Write the data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing "' // trim ( filename_mesh ) // '".'

  call mesh_write ( filename_mesh, dim, vertices, edges, triangles, &
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
  write ( *, '(a)' ) 'ICE_TO_MESH:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine data_print ( dim, vertices, edges, triangles, quadrilaterals, &
  tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, &
  edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
  quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
  hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! DATA_PRINT prints the data of an ICE grid dataset.
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
!    Input, double VERTEX_COORDINATE(DIM,VERTICES), the coordinates
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
!    Input, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS),
!    the vertices that form each quadrilateral.
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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  real ( kind = 8 ) vertex_coordinate(dim,vertices)
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
subroutine data_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
  vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! DATA_READ reads ICE data from a NETCDF file.
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
!    Output, real VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Output, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
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
  integer ( kind = 4 ) dimid
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
  integer ( kind = 4 ) varid
  real ( kind = 8 ) vertex_coordinate(3,vertices)
  integer ( kind = 4 ) vertex_label(vertices)
  integer ( kind = 4 ) xtype
!
!  Open the file.
!
  mode = NF90_NOCLOBBER
  status = nf90_open ( filename, mode, ncid )
!
!  Vertices.
!
  status = nf90_inq_varid ( ncid, 'Vertex_Coordinate', varid )
  status = nf90_get_var ( ncid, varid, vertex_coordinate )

  status = nf90_inq_varid ( ncid, 'Vertex_Label', varid )
  status = nf90_get_var ( ncid, varid, vertex_label )
!
!  Edges.
!
  if ( 0 < edges ) then
    status = nf90_inq_varid ( ncid, 'Edge_Vertex', varid )
    status = nf90_get_var ( ncid, varid, edge_vertex )

    status = nf90_inq_varid ( ncid, 'Edge_Label', varid )
    status = nf90_get_var ( ncid, varid, edge_label )
  end if
!
!  Triangles.
!
  if ( 0 < triangles ) then
    status = nf90_inq_varid ( ncid, 'Triangle_Vertex', varid )
    status = nf90_get_var ( ncid, varid, triangle_vertex )

    status = nf90_inq_varid ( ncid, 'Triangle_Label', varid )
    status = nf90_get_var ( ncid, varid, triangle_label )
  end if
!
!  Quadrilaterals.
!
  if ( 0 < quadrilaterals ) then
    status = nf90_inq_varid ( ncid, 'Quadrilateral_Vertex', varid )
    status = nf90_get_var ( ncid, varid, quadrilateral_vertex )

    status = nf90_inq_varid ( ncid, 'Quadrilateral_Label', varid )
    status = nf90_get_var ( ncid, varid, quadrilateral_label )
  end if
!
!  Tetrahedrons.
!
  if ( 0 < tetrahedrons ) then
    status = nf90_inq_varid ( ncid, 'Tetrahedron_Vertex', varid )
    status = nf90_get_var ( ncid, varid, tetrahedron_vertex )

    status = nf90_inq_varid ( ncid, 'Tetrahedron_Label', varid )
    status = nf90_get_var ( ncid, varid, tetrahedron_label )
  end if
!
!  Hexahedrons.
!
  if ( 0 < hexahedrons ) then
    status = nf90_inq_varid ( ncid, 'Hexahedron_Vertex', varid )
    status = nf90_get_var ( ncid, varid, hexahedron_vertex )

    status = nf90_inq_varid ( ncid, 'Hexahedron_Label', varid )
    status = nf90_get_var ( ncid, varid, hexahedron_label )
  end if
!
!  Close the file.
!
  status = nf90_close ( ncid )

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
!    Input, real ( kind = 8 ) VERTEX_COORDINATE(DIM,VERTICES), the coordinates
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
subroutine size_print ( dim, vertices, edges, triangles, quadrilaterals, &
  tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! SIZE_PRINT prints the sizes of an ICE grid dataset.
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
subroutine size_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! SIZE_READ reads ICE sizes from a NETCDF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2010
!
!  Author:
!
!   John Burkardt
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
!    Input, string FILENAME, the name of the file to be read.
!    Ordinarily, the name should include the extension '.nc'.
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
!    Output, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Output, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
  use netcdf

  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dimid
  integer ( kind = 4 ) edges
  character ( len = * ) filename
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) ncid
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) status
  integer ( kind = 4 ) tetrahedrons
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
  mode = NF90_NOWRITE
  status = nf90_open ( filename, mode, ncid )
!
!  Get the dimension information.
!
!  In an act of perplexing effect, the F90 NETCDF does not
!  include a "NF90_INQ_DIMLEN" function to allow you to determine
!  the length of a dimension from its ID.  Instead, NF90_INQUIRE_DIMENSION
!  returns both name and length.  Sadly, then, I must use the name
!  to get the ID, and then the ID to get the dimension...and the name again.
!  Oh, I see, this lets them use OPTIONAL arguments.  Whoopee!
!
  if ( nf90_inq_dimid ( ncid, 'Dimension', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = dim )
  end if

  if ( nf90_inq_dimid ( ncid, 'Vertices', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = vertices )
  end if

  if ( nf90_inq_dimid ( ncid, 'Edges', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = edges )
  end if

  if ( nf90_inq_dimid ( ncid, 'Triangles', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = triangles )
  end if

  if ( nf90_inq_dimid ( ncid, 'Quadrilaterals', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = quadrilaterals )
  end if

  if ( nf90_inq_dimid ( ncid, 'Tetrahedra', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = tetrahedrons )
  end if

  if ( nf90_inq_dimid ( ncid, 'Tetrahedrons', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = tetrahedrons )
  end if

  if ( nf90_inq_dimid ( ncid, 'Hexahedra', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = hexahedrons )
  end if

  if ( nf90_inq_dimid ( ncid, 'Hexahedrons', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = hexahedrons )
  end if
!
!  Close the file.
!
  status = nf90_close ( ncid )

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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
