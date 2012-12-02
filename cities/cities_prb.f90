program main

!*****************************************************************************80
!
!! MAIN is the main program for CITIES_PRB.
!
!  Discussion:
!
!    CITIES_PRB tests routines from the CITIES library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CITIES_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CITIES library.'

  call test01 ( 'wg22' )
  call test02 ( 'usca312' )
  call test03 ( 'usca312' )
  call test04 ( 'spaeth2_09' )
  call test05 ( 'uscap' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CITIES_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( prefix )

!*****************************************************************************80
!
!! TEST01 tests POINT_TO_DIST_TABLE.
!
!  Discussion:
!
!    Get the XY coordinates of a set of cities, and compute
!    the city-to-city distance table.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) PREFIX, the common file prefix.
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: dist
  character ( len = 255 ) :: dist_filename
  character ( len = 255 ) :: main_filename
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename
  integer ( kind = 4 ) point_num
  character ( len = * )  prefix

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  POINT_TO_DIST_TABLE computes a distance table from a'
  write ( *, '(a)' ) '  list of point locations.'

  call s_cat ( prefix, '_main.txt', main_filename )
  call s_cat ( prefix, '_xy.txt', point_filename )
  call s_cat ( prefix, '_dist.txt', dist_filename )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The main filename is "' // trim ( main_filename ) // '"'
  write ( *, '(a)' ) '  The point filename is "' // trim ( point_filename ) // '"'
  write ( *, '(a)' ) '  The distance table filename will be "' // trim ( dist_filename ) // '"'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension is ', dim_num
  write ( *, '(a,i8)' ) '  The number of points is  ', point_num

  allocate ( point(1:dim_num,1:point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )

  call r8mat_transpose_print ( dim_num, point_num, point, '  The points:' )

  allocate ( dist(1:point_num,1:point_num) )

  call point_to_dist_table ( dim_num, point_num, point, dist )

  dist(1:point_num,1:point_num) = aint ( dist(1:point_num,1:point_num) )

  call r8mat_print_some ( point_num, point_num, dist, 1, 1, 5, 5, &
    '  Initial 5x5 distance subtable:' )

  call r8mat_write ( dist_filename, point_num, point_num, dist )

  deallocate ( dist )
  deallocate ( point )

  return
end
subroutine test02 ( prefix )

!*****************************************************************************80
!
!! TEST02 tests MAIN_READ_SIZE, MAIN_READ_DMS, MAIN_READ_NAME.
!
!  Discussion:
!
!    Get the DMS coordinates of a set of cities, and compute
!    the city-to-city distance table, using distances on a sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) PREFIX, the common file prefix.
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: dist
  character ( len = 255 ) :: dist_filename
  character ( len = 255 ) :: dms_filename
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: lat_dms
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: long_dms
  character ( len = 255 ) :: main_filename
  integer ( kind = 4 ) n
  character ( len = * ) prefix

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Get the DMS coordinates of a set of cities.'
  write ( *, '(a)' ) '  Compute the city-to-city distance table,'
  write ( *, '(a)' ) '  assuming the cities lie on a sphere (the earth).'

  call s_cat ( prefix, '_main.txt', main_filename )
  call s_cat ( prefix, '_dms.txt', dms_filename )
  call s_cat ( prefix, '_dist.txt', dist_filename )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The main filename is "' // trim ( main_filename ) // '"'
  write ( *, '(a)' ) '  The dms filename is "' // trim ( dms_filename ) // '"'
  write ( *, '(a)' ) '  The distance filename will be "' // trim ( dist_filename ) // '"'

  call main_read_size ( main_filename, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of data items is ', n

  allocate ( dist(1:n,1:n) )
  allocate ( lat_dms(1:4,1:n) )
  allocate ( long_dms(1:4,1:n) )

  call dms_read ( dms_filename, n, lat_dms, long_dms )

  call dms_print ( n, lat_dms, long_dms, '  The longitude/latitude data:' )

  call dms_to_dist ( n, lat_dms, long_dms, dist )

  dist(1:n,1:n) = aint ( dist(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Distance from Atlanta to Boston  = ', dist(12,34)
  write ( *, '(a)' ) '  Road distance is 1037'
  write ( *, '(a,g14.6)' ) '  Distance from Boston to Chicago  = ', dist(34,58)
  write ( *, '(a)' ) '  Road distance is  963'
  write ( *, '(a,g14.6)' ) '  Distance from Chicago to Atlanta = ', dist(58,12)
  write ( *, '(a)' ) '  Road distance is  674'

  call r8mat_write ( dist_filename, n, n, dist )

  deallocate ( dist )
  deallocate ( lat_dms )
  deallocate ( long_dms )

  return
end
subroutine test03 ( prefix )

!*****************************************************************************80
!
!! TEST03 tests DMS_TO_XY.
!
!  Discussion:
!
!    Get the DMS coordinates of a set of cities, and compute
!    the XY coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) PREFIX, the common file prefix.
!
  implicit none

  character ( len = 255 ) :: dms_filename
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: lat_dms
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: long_dms
  character ( len = 255 ) :: main_filename
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_xy
  character ( len = 255 ) :: point_xy_filename
  character ( len = * ) prefix

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  DMS_TO_XY takes latitude and longitude'
  write ( *, '(a)' ) '  information, and assigns pseudo XY coordinates.'

  call s_cat ( prefix, '_main.txt', main_filename )
  call s_cat ( prefix, '_dms.txt', dms_filename )
  call s_cat ( prefix, '_xy.txt', point_xy_filename )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The main filename is "' // trim ( main_filename ) // '"'
  write ( *, '(a)' ) '  The dms filename is "' // trim ( dms_filename ) // '"'
  write ( *, '(a)' ) '  The point XY filename will be "' // trim ( point_xy_filename ) // '"'

  call main_read_size ( main_filename, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of data items is ', n

  allocate ( lat_dms(1:4,1:n) )
  allocate ( long_dms(1:4,1:n) )
  allocate ( point_xy(1:2,1:n) )

  call dms_read ( dms_filename, n, lat_dms, long_dms )

  call dms_print ( n, lat_dms, long_dms, '  The longitude/latidude data:' )

  call dms_to_xy ( n, lat_dms, long_dms, point_xy )

  call r8mat_transpose_print ( 2, n, point_xy, '  The computed point values:' )

  call r8mat_write ( point_xy_filename, 2, n, point_xy )

  deallocate ( lat_dms )
  deallocate ( long_dms )
  deallocate ( point_xy )

  return
end
subroutine test04 ( prefix )

!*****************************************************************************80
!
!! TEST04 tests DIST_TABLE_CHECK.
!
!  Discussion:
!
!    Read a distance matrix and check it for consistency.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) PREFIX, the common file prefix.
!
  implicit none

  integer ( kind = 4 ) check
  real ( kind = 8 ), allocatable :: dist_table(:,:)
  character ( len = 255 ) dist_table_filename
  integer ( kind = 4 ) dist_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  character ( len = * )  prefix

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DIST_TABLE_CHECK checks a distance table.'

  call s_cat ( prefix, '_dist.txt', dist_table_filename )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The distance table filename is "' // trim ( dist_table_filename ) // '"'

  call r8mat_header_read ( dist_table_filename, n1, n2 )

  if ( n1 /= n2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a)' ) '  The distance table is not square.'
    return
  end if

  n = n1
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is ', n

  allocate ( dist_table(n,n) )

  call r8mat_data_read ( dist_table_filename, n, n, dist_table )

  call dist_table_check ( n, dist_table, check )

  write ( *, '(a)' ) ' '
  if ( check == 0 ) then
    write ( *, '(a)' ) '  0: The distance table passed all checks.'
  else if ( check == 1 ) then
    write ( *, '(a)' ) '  1: The table failed the nonnegativity check.'
  else if ( check == 2 ) then
    write ( *, '(a)' ) '  2: The table failed the zero self-distance check.'
  else if ( check == 3 ) then
    write ( *, '(a)' ) '  3: The table failed the symmetry check.'
  else if ( check == 4 ) then
    write ( *, '(a)' ) '  4: The table failed the triangle check.'
  end if

  deallocate ( dist_table )

  return
end
subroutine test05 ( prefix )

!*****************************************************************************80
!
!! TEST05 tests LL_DEGREES_TO_XY.
!
!  Discussion:
!
!    Get the LL coordinates of a set of cities, and compute
!    the XY coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) PREFIX, the common file prefix.
!
  implicit none

  character ( len = 255 ) :: ll_filename
  real ( kind = 8 ), allocatable, dimension ( : ) :: lat
  real ( kind = 8 ), allocatable, dimension ( : ) :: long
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  character ( len = 255 ) :: xy_filename
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  character ( len = * ) prefix

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  LL_DEGREES_TO_XY takes latitude and longitude'
  write ( *, '(a)' ) '  information, and assigns pseudo XY coordinates.'

  call s_cat ( prefix, '_ll.txt', ll_filename )
  call s_cat ( prefix, '_xy.txt', xy_filename )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The LL filename is "' // trim ( ll_filename ) // '"'
  write ( *, '(a)' ) '  The XY filename is "' // trim ( xy_filename ) // '"'

  call r8vec2_header_read ( ll_filename, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of data items is ', n

  allocate ( lat(1:n) )
  allocate ( long(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  call r8vec2_data_read ( ll_filename, n, lat, long )

  call r8vec2_print ( n, lat, long, '  The longitude/latidude data:' )

  call ll_degrees_to_xy ( n, lat, long, x, y )

  call r8vec2_print ( n, x, y, '  The computed point values:' )

  call r8vec2_write ( xy_filename, n, x, y )

  deallocate ( lat )
  deallocate ( long )
  deallocate ( x )
  deallocate ( y )

  return
end
