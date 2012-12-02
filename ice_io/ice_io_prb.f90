program main

!*****************************************************************************80
!
!! ICE_IO_PRB tests the ICE_IO library.
!
!  Discussion:
!
!    We begin by creating a file.
!
!    The FORTRAN90 version of NETCDF was so unpleasant to install
!    (modules, what a concept!) that this file was finished a month
!    after the C version.
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
!    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
!    The NETCDF User"s Guide,
!    Unidata Program Center, March 2009.
!
  implicit none

  character ( len = 255 ) filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ICE_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ICE_IO library.'
!
!  Create "hexahexa_2x2x2.nc"
!
  filename = 'hexahexa_2x2x2.nc'
  call test01 ( filename )
!
!  Read "hexahexa_2x2x2.nc"
!
   call test02 ( filename )
!
!  Create "cyl248.nc"
!
   filename = 'cyl248.nc'
   call test03 ( filename )
!
!  Read "cyl248.nc"
!
   call test02 ( filename )
!
! Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ICE_IO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( filename )

!*****************************************************************************80
!
!! TEST01 creates the HEXAHEXA_2X2X2 dataset and writes it to NETCDF.
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
!  Reference:
!
!    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
!    The NETCDF User"s Guide,
!    Unidata Program Center, March 2009.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to be created.
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
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Create an ICE grid dataset, print it,'
  write ( *, '(a)' ) '  and write it to an NETCDF file.'
!
!  Get sizes.
!
  call hexahexa_2x2x2_size ( dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons )
!
!  Print sizes;
!
  call size_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons )
!
!  Allocate memory.
!
  allocate ( vertex_coordinate(3,vertices) )
  allocate ( vertex_label(vertices) )
  allocate ( edge_vertex(2,edges) )
  allocate ( edge_label(edges) )
  allocate ( triangle_vertex(3,triangles) )
  allocate ( triangle_label(triangles) )
  allocate ( quadrilateral_vertex(4,quadrilaterals) )
  allocate ( quadrilateral_label(quadrilaterals) )
  allocate ( tetrahedron_vertex(4,tetrahedrons) )
  allocate ( tetrahedron_label(tetrahedrons) )
  allocate ( hexahedron_vertex(8,hexahedrons) )
  allocate ( hexahedron_label(hexahedrons) )
!
!  Get data.
!
  call hexahexa_2x2x2_data ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, &
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
    hexahedron_vertex, hexahedron_label )
!
!  Print the data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data to be written to "' // trim ( filename ) // '".'

  call data_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, &
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
    hexahedron_vertex, hexahedron_label )
!
!  Create the file.
!
  call ice_write ( filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the file "' // trim ( filename ) // '".'
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

  return
end
subroutine test02 ( filename )

!*****************************************************************************80
!
!! TEST02 reads an ICE grid dataset from a NETCDF file.
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
!    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
!    The NETCDF User"s Guide,
!    Unidata Program Center, March 2009.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to be read.
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
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Read an ICE grid dataset from a NETCDF file,'
  write ( *, '(a)' ) '  and print the data.'
!
!  Read sizes;
!
  call size_read ( filename, dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons )
!
!  Print sizes;
!
  call size_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons )
!
!  Allocate memory.
!
  allocate ( vertex_coordinate ( 3, vertices ) )
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
!  Read the file
!
  call data_read ( filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )
!
!  Print the data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data from file ' // trim ( filename ) // '".'

  call data_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, &
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
    hexahedron_vertex, hexahedron_label )
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

  return
end
subroutine test03 ( filename )

!*****************************************************************************80
!
!! TEST03 creates the HEXAHEXA_2X2X2 dataset and writes it to NETCDF.
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
!    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
!    The NETCDF User"s Guide,
!    Unidata Program Center, March 2009.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to be created.
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
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  Create an ICE grid dataset, print it,'
  write ( *, '(a)' ) '  and write it to an NETCDF file.'
!
!  Get sizes.
!
  call cyl248_size ( dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons )
!
!  Print sizes;
!
  call size_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons )
!
!  Allocate memory.
!
  allocate ( vertex_coordinate(3,vertices) )
  allocate ( vertex_label(vertices) )
  allocate ( edge_vertex(2,edges) )
  allocate ( edge_label(edges) )
  allocate ( triangle_vertex(3,triangles) )
  allocate ( triangle_label(triangles) )
  allocate ( quadrilateral_vertex(4,quadrilaterals) )
  allocate ( quadrilateral_label(quadrilaterals) )
  allocate ( tetrahedron_vertex(4,tetrahedrons) )
  allocate ( tetrahedron_label(tetrahedrons) )
  allocate ( hexahedron_vertex(8,hexahedrons) )
  allocate ( hexahedron_label(hexahedrons) )
!
!  Get data.
!
  call cyl248_data ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, &
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
    hexahedron_vertex, hexahedron_label )
!
!  Print the data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data to be written to "' // trim ( filename ) // '".'

  call data_print ( dim, vertices, edges, triangles, quadrilaterals, &
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, &
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
    hexahedron_vertex, hexahedron_label )
!
!  Create the file.
!
  call ice_write ( filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the file "' // trim ( filename ) // '".'
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

  return
end
