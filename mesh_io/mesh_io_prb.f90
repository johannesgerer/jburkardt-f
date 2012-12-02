program main

!*****************************************************************************80
!
!! MESH_IO_PRB tests the MESH_IO library.
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
  implicit none

  character ( len = 255 ) filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MESH_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the MESH_IO library.'
!
!  Create the file hexahexa_2x2x2.mesh
!
  call test01 ( )
!
!  Read and print the file hexahexa_2x2x2.mesh.
!
  filename = 'hexahexa_2x2x2.mesh'
  call test03 ( filename )
!
!  Create the file cyl248.mesh
!
  call test02 ( )
!
!  Read and print the sizes of file cyl248.mesh.
!
  filename = 'cyl248.mesh'
  call test03 ( filename )
!
!  Read and print the data in file cyl248.mesh.
!
  filename = 'cyl248.mesh'
  call test04 ( filename )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MESH_IO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 creates a MESH dataset and writes it to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ), allocatable :: edge_label(:)
  integer ( kind = 4 ), allocatable :: edge_vertex(:,:)
  integer ( kind = 4 ) edges
  character ( len = 255 ) filename
  integer ( kind = 4 ), allocatable :: hexahedron_label(:)
  integer ( kind = 4 ), allocatable :: hexahedron_vertex(:,:)
  integer ( kind = 4 ) hexahedrons
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Create a hexahedral mesh and write it to a file.'
!
!  Get sizes.
!
  call hexahexa_2x2x2_size ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons )
!
!  Allocate memory.
!
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
  allocate ( vertex_coordinate(3,vertices) )
  allocate ( vertex_label(vertices) )
!
!  Get the data.
!
  call hexahexa_2x2x2_data ( dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, &
    edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )
!
!  Write the data.
!
  filename = 'hexahexa_2x2x2.mesh';

  call mesh_write ( filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the file "' // trim ( filename ) // '".'
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

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 creates a MESH dataset and writes it to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ), allocatable :: edge_label(:)
  integer ( kind = 4 ), allocatable :: edge_vertex(:,:)
  integer ( kind = 4 ) edges
  character ( len = 255 ) filename
  integer ( kind = 4 ), allocatable :: hexahedron_label(:)
  integer ( kind = 4 ), allocatable :: hexahedron_vertex(:,:)
  integer ( kind = 4 ) hexahedrons
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Create a tetrahedral mesh and write it to a file.'
!
!  Get sizes.
!
  call cyl248_size ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons )
!
!  Allocate memory.
!
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
  allocate ( vertex_coordinate(3,vertices) )
  allocate ( vertex_label(vertices) )
!
!  Get the data.
!
  call cyl248_data ( dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, &
    edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )
!
!  Write the data.
!
  filename = 'cyl248.mesh';

  call mesh_write ( filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the file "' // trim ( filename ) // '".'
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

  return
end
subroutine test03 ( filename )

!*****************************************************************************80
!
!! MESH_IO_TEST03 reads and prints the sizes in a MESH dataset.
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
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  character ( len = * ) filename
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MESH_IO_TEST03'
  write ( *, '(a)' ) '  Read a mesh file and print its sizes.'
!
!  Read sizes.
!
  call mesh_size_read ( filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons )
!
!  Print sizes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Header information for "' // trim ( filename ) // '".'

  call mesh_size_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons )

  return
end
subroutine test04 ( filename )

!*****************************************************************************80
!
!! MESH_IO_TEST04 reads a MESH dataset and prints its data.
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
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ), allocatable :: edge_label(:)
  integer ( kind = 4 ), allocatable :: edge_vertex(:,:)
  integer ( kind = 4 ) edges
  character ( len = * ) filename
  integer ( kind = 4 ), allocatable :: hexahedron_label(:)
  integer ( kind = 4 ), allocatable :: hexahedron_vertex(:,:)
  integer ( kind = 4 ) hexahedrons
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MESH_IO_TEST04'
  write ( *, '(a)' ) '  Read a mesh file and print its data.'
!
!  Read sizes.
!
  call mesh_size_read ( filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons )
!
!  Allocate memory.
!
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
  allocate ( vertex_coordinate(3,vertices) )
  allocate ( vertex_label(vertices) )
!
!  Read the data.
!
  call mesh_data_read ( filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )
!
!  Print the data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data for file "' // trim ( filename ) // '".'

  call mesh_data_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,  edge_vertex, &
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
    hexahedron_vertex, hexahedron_label )
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

  return
end
