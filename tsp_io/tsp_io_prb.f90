program main

!*****************************************************************************80
!
!! TSP_IO_PRB tests the TSP_IO library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TSP_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TSP_IO library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TSP_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 reports the data in 'P01.TSP'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 )    tsp_capacity
  character ( len = 255 ) tsp_comment
  integer ( kind = 4 )    tsp_dimension
  character ( len = 255 ) tsp_display_data_type
  character ( len = 255 ) tsp_edge_data_format
  character ( len = 255 ) tsp_edge_weight_format
  character ( len = 255 ) tsp_edge_weight_type
  character ( len = 255 ) tsp_filename
  character ( len = 255 ) tsp_name
  character ( len = 255 ) tsp_node_coord_type
  character ( len = 255 ) tsp_type

  tsp_filename = 'p01.tsp'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Read the header from "' // trim ( tsp_filename ) // '"'
  write ( *, '(a)' ) '  and print the values.'

  call tsp_header_read ( tsp_filename, tsp_name, tsp_type, tsp_comment, &
    tsp_dimension, tsp_capacity, tsp_edge_weight_type, tsp_edge_weight_format, &
    tsp_edge_data_format, tsp_node_coord_type, tsp_display_data_type )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NAME: "' // trim ( tsp_name ) // '".'
  write ( *, '(a)' ) '  TYPE: "' // trim ( tsp_type ) // '".'
  write ( *, '(a)' ) '  COMMENT: "' // trim ( tsp_comment ) // '".'
  write ( *, '(a,i8)' ) '  DIMENSION: ', tsp_dimension
  write ( *, '(a,i8)' ) '  CAPACITY: ', tsp_capacity
  write ( *, '(a)' ) '  EDGE_WEIGHT_TYPE: "' // trim ( tsp_edge_weight_type ) // '".'
  write ( *, '(a)' ) '  EDGE_WEIGHT_FORMAT: "' // trim ( tsp_edge_weight_format ) // '".'
  write ( *, '(a)' ) '  EDGE_DATA_FORMAT: "' // trim ( tsp_edge_data_format ) // '".'
  write ( *, '(a)' ) '  NODE_COORD_TYPE: "' // trim ( tsp_node_coord_type ) // '".'
  write ( *, '(a)' ) '  DISPLAY_DATA_TYPE: "' // trim ( tsp_display_data_type ) // '".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 retrieves the edge weight data in 'P01.TSP'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 )    tsp_capacity
  character ( len = 255 ) tsp_comment
  integer ( kind = 4 )    tsp_dimension
  character ( len = 255 ) tsp_display_data_type
  character ( len = 255 ) tsp_edge_data_format
  integer ( kind = 4 ), allocatable :: tsp_edge_weight(:,:)
  character ( len = 255 ) tsp_edge_weight_format
  character ( len = 255 ) tsp_edge_weight_type
  character ( len = 255 ) tsp_filename
  character ( len = 255 ) tsp_name
  character ( len = 255 ) tsp_node_coord_type
  character ( len = 255 ) tsp_type

  tsp_filename = 'p01.tsp'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Read the edge weights from "' // trim ( tsp_filename ) // '"'
  write ( *, '(a)' ) '  and print them.'

  call tsp_header_read ( tsp_filename, tsp_name, tsp_type, tsp_comment, &
    tsp_dimension, tsp_capacity, tsp_edge_weight_type, tsp_edge_weight_format, &
    tsp_edge_data_format, tsp_node_coord_type, tsp_display_data_type )

  allocate ( tsp_edge_weight(tsp_dimension,tsp_dimension) )

  call tsp_edge_weight_read ( tsp_filename, tsp_dimension, &
    tsp_edge_weight_type, tsp_edge_weight_format, tsp_edge_weight )

  call i4mat_print ( tsp_dimension, tsp_dimension, tsp_edge_weight, &
    '  Edge Weight matrix (distances)' )

  deallocate ( tsp_edge_weight )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 retrieves the edge weight data in 'P01.TSP' and writes it to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) output_filename
  integer ( kind = 4 )    tsp_capacity
  character ( len = 255 ) tsp_comment
  integer ( kind = 4 )    tsp_dimension
  character ( len = 255 ) tsp_display_data_type
  character ( len = 255 ) tsp_edge_data_format
  integer ( kind = 4 ), allocatable :: tsp_edge_weight(:,:)
  character ( len = 255 ) tsp_edge_weight_format
  character ( len = 255 ) tsp_edge_weight_type
  character ( len = 255 ) tsp_filename
  character ( len = 255 ) tsp_name
  character ( len = 255 ) tsp_node_coord_type
  character ( len = 255 ) tsp_type

  tsp_filename = 'p01.tsp'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Read the edge weights from "' // trim ( tsp_filename ) // '"'
  write ( *, '(a)' ) '  and write them to a file.'

  call tsp_header_read ( tsp_filename, tsp_name, tsp_type, tsp_comment, &
    tsp_dimension, tsp_capacity, tsp_edge_weight_type, tsp_edge_weight_format, &
    tsp_edge_data_format, tsp_node_coord_type, tsp_display_data_type )

  allocate ( tsp_edge_weight(tsp_dimension,tsp_dimension) )

  call tsp_edge_weight_read ( tsp_filename, tsp_dimension, &
    tsp_edge_weight_type, tsp_edge_weight_format, tsp_edge_weight )

  output_filename = 'p01_d.txt'
  call i4mat_write ( output_filename, tsp_dimension, tsp_dimension, &
    tsp_edge_weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Distance table written to "' &
    // trim ( output_filename ) // '".'

  deallocate ( tsp_edge_weight )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 retrieves the edge weight data in 'ATT48.TSP' and writes it to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) output_filename
  integer ( kind = 4 )    tsp_capacity
  character ( len = 255 ) tsp_comment
  integer ( kind = 4 )    tsp_dimension
  character ( len = 255 ) tsp_display_data_type
  character ( len = 255 ) tsp_edge_data_format
  integer ( kind = 4 ), allocatable :: tsp_edge_weight(:,:)
  character ( len = 255 ) tsp_edge_weight_format
  character ( len = 255 ) tsp_edge_weight_type
  character ( len = 255 ) tsp_filename
  character ( len = 255 ) tsp_name
  character ( len = 255 ) tsp_node_coord_type
  character ( len = 255 ) tsp_type

  tsp_filename = 'att48.tsp'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Read the edge weights from "' // trim ( tsp_filename ) // '"'
  write ( *, '(a)' ) '  and write them to a file.'

  call tsp_header_read ( tsp_filename, tsp_name, tsp_type, tsp_comment, &
    tsp_dimension, tsp_capacity, tsp_edge_weight_type, tsp_edge_weight_format, &
    tsp_edge_data_format, tsp_node_coord_type, tsp_display_data_type )

  allocate ( tsp_edge_weight(tsp_dimension,tsp_dimension) )

  call tsp_edge_weight_read ( tsp_filename, tsp_dimension, &
    tsp_edge_weight_type, tsp_edge_weight_format, tsp_edge_weight )

  output_filename = 'att48_d.txt'
  call i4mat_write ( output_filename, tsp_dimension, tsp_dimension, &
    tsp_edge_weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Distance table written to "' &
    // trim ( output_filename ) // '".'

  deallocate ( tsp_edge_weight )

  return
end
