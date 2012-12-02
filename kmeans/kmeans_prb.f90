program main

!*****************************************************************************80
!
!! MAIN is the main program for KMEANS_PRB.
!
!  Discussion:
!
!    KMEANS_PRB tests various KMEANS programs.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KMEANS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the KMEANS library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KMEANS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tries out the HMEANS_01 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test01_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test01_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test the HMEANS_01 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call hmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tries out the HMEANS_02 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test02_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test02_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Test the HMEANS_02 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call hmeans_02 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tries out the KMEANS_01 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test03_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test03_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Test the KMEANS_01 algorithm.'
  write ( *, '(a)' ) '  (Applied Statistics Algorithm #58)'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call kmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of KMEANS_01 iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tries out the KMEANS_02 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test04_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test04_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test the KMEANS_02 algorithm.'
  write ( *, '(a)' ) '  (Applied Statistics Algorithm #136)'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the centers.
!
  call cluster_initialize_1 ( dim_num, point_num, cluster_num, point, &
    cluster_center )

  cluster(1:cluster_num) = 0

  call kmeans_02 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tries out the KMEANS_03 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test05_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test05_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Test the KMEANS_03 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the centers.
!
  call cluster_initialize_1 ( dim_num, point_num, cluster_num, point, &
    cluster_center )

  cluster(1:cluster_num) = 0

  call kmeans_03 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tries out the HMEANS_01 + KMEANS_01 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test06_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test06_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_max_h
  integer ( kind = 4 ) it_max_k
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Test the HMEANS_01 + KMEANS_01 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789
  it_max_h = 3
  it_max_k = it_max

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of HMEANS_01 iterations allowed is ', it_max_h
  write ( *, '(a,i8)' ) '  Number of KMEANS_01 iterations allowed is ', it_max_k
!
!  Initialize the centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of HMEANS_01 iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call kmeans_01 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of KMEANS_01 iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tries out the HMEANS_01 + KMEANS_02 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test07_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test07_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_max_h
  integer ( kind = 4 ) it_max_k
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Test the HMEANS_01 + KMEANS_02 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  it_max_h = 3
  it_max_k = it_max

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of HMEANS_01 iterations allowed is ', it_max_h
  write ( *, '(a,i8)' ) '  Number of KMEANS_02 iterations allowed is ', it_max_k
!
!  Initialize the centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of HMEANS_01 iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call kmeans_02 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of KMEANS_02 iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tries out the HMEANS_01 + KMEANS_03 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test08_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test08_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_max_h
  integer ( kind = 4 ) it_max_k
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Test the HMEANS_01 + KMEANS_03 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789
  it_max_h = 3
  it_max_k = it_max

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initialize by using a few steps of HMEANS_02:'
  write ( *, '(a,i8)' ) '  Number of HMEANS_01 iterations allowed is ', it_max_h
  write ( *, '(a,i8)' ) '  Number of KMEANS_03 iterations allowed is ', it_max_k
!
!  Initialize the centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )
!
!  Initialize the clusters.
!
  cluster(1:cluster_num) = 0

  call hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of HMEANS_01 iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call kmeans_03 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, &
    cluster, cluster_center, cluster_population, cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of KMEANS_03 iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tries out the HMEANS_W_01 routine.
!
!  Discussion:
!
!    The weights are all equal, so the results should
!    be identical to those for HMEANS_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test09_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test09_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: weight(:)
  integer ( kind = 4 ) weight_dim
  integer ( kind = 4 ) weight_num
  character ( len = 255 ) :: weight_filename = 'weights_equal_100.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  Test the HMEANS_W_01 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Read the weights.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights will be read from "' // trim ( weight_filename ) // '".'

  call r8mat_header_read ( weight_filename, weight_dim, weight_num )

  if ( weight_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of weight array is not 1.'
    stop
  end if

  if ( weight_num /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Number of weights not equal to number of points.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( weight_filename, weight_dim, weight_num, weight )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the cluster centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call hmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, &
    cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tries out the HMEANS_W_02 routine.
!
!  Discussion:
!
!    The weights are all equal, so the results should
!    be identical to those for HMEANS_02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test10_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test10_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: weight(:)
  integer ( kind = 4 ) weight_dim
  integer ( kind = 4 ) weight_num
  character ( len = 255 ) :: weight_filename = 'weights_equal_100.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  Test the HMEANS_W_02 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Read the weights.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights will be read from "' // trim ( weight_filename ) // '".'

  call r8mat_header_read ( weight_filename, weight_dim, weight_num )

  if ( weight_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of weight array is not 1.'
    stop
  end if

  if ( weight_num /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Number of weights not equal to number of points.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( weight_filename, weight_dim, weight_num, weight )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the cluster centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call hmeans_w_02 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, &
    cluster_energy, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tries out the KMEANS_W_01 routine.
!
!  Discussion:
!
!   The weights are all equal, so the results should
!    be identical to those for KMEANS_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test11_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test11_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: weight(:)
  integer ( kind = 4 ) weight_dim
  integer ( kind = 4 ) weight_num
  character ( len = 255 ) :: weight_filename = 'weights_equal_100.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Test the KMEANS_W_01 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Read the weights.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights will be read from "' // trim ( weight_filename ) // '".'

  call r8mat_header_read ( weight_filename, weight_dim, weight_num )

  if ( weight_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of weight array is not 1.'
    stop
  end if

  if ( weight_num /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Number of weights not equal to number of points.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( weight_filename, weight_dim, weight_num, weight )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the cluster centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call kmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, &
    cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tries out the KMEANS_W_03 routine.
!
!  Discussion:
!
!    The weights are all equal, so the results should
!    be identical to those for KMEANS_03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test12_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test12_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: weight(:)
  integer ( kind = 4 ) weight_dim
  integer ( kind = 4 ) weight_num
  character ( len = 255 ) :: weight_filename = 'weights_equal_100.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  Test the KMEANS_W_03 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Read the weights.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights will be read from "' // trim ( weight_filename ) // '".'

  call r8mat_header_read ( weight_filename, weight_dim, weight_num )

  if ( weight_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of weight array is not 1.'
    stop
  end if

  if ( weight_num /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Number of weights not equal to number of points.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( weight_filename, weight_dim, weight_num, weight )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the cluster centers.
!
  call cluster_initialize_1 ( dim_num, point_num, cluster_num, point, &
    cluster_center )

  cluster(1:cluster_num) = 0

  call kmeans_w_03 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, &
    cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tries out the HMEANS_W_01 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test13_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test13_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: weight(:)
  integer ( kind = 4 ) weight_dim
  integer ( kind = 4 ) weight_num
  character ( len = 255 ) :: weight_filename = 'weights_unequal_100.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  Test the HMEANS_W_01 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Read the weights.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights will be read from "' // trim ( weight_filename ) // '".'

  call r8mat_header_read ( weight_filename, weight_dim, weight_num )

  if ( weight_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of weight array is not 1.'
    stop
  end if

  if ( weight_num /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Number of weights not equal to number of points.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( weight_filename, weight_dim, weight_num, weight )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the cluster centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call hmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, &
    cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tries out the HMEANS_W_02 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test14_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test14_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: weight(:)
  integer ( kind = 4 ) weight_dim
  integer ( kind = 4 ) weight_num
  character ( len = 255 ) :: weight_filename = 'weights_unequal_100.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  Test the HMEANS_W_02 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Read the weights.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights will be read from "' // trim ( weight_filename ) // '".'

  call r8mat_header_read ( weight_filename, weight_dim, weight_num )

  if ( weight_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of weight array is not 1.'
    stop
  end if

  if ( weight_num /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Number of weights not equal to number of points.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( weight_filename, weight_dim, weight_num, weight )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the cluster centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call hmeans_w_02 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, &
    cluster_energy, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tries out the KMEANS_W_01 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test15_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test15_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: weight(:)
  integer ( kind = 4 ) weight_dim
  integer ( kind = 4 ) weight_num
  character ( len = 255 ) :: weight_filename = 'weights_unequal_100.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  Test the KMEANS_W_01 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Read the weights.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights will be read from "' // trim ( weight_filename ) // '".'

  call r8mat_header_read ( weight_filename, weight_dim, weight_num )

  if ( weight_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of weight array is not 1.'
    stop
  end if

  if ( weight_num /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Number of weights not equal to number of points.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( weight_filename, weight_dim, weight_num, weight )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the cluster centers.
!
  call cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

  cluster(1:cluster_num) = 0

  call kmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, &
    cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tries out the KMEANS_W_03 routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  character ( len = 255 ) :: center_filename = 'test16_centers.txt'
  integer ( kind = 4 ), allocatable :: cluster(:)
  real ( kind = 8 ), allocatable :: cluster_center(:,:)
  real ( kind = 8 ), allocatable :: cluster_energy(:)
  character ( len = 255 ) :: cluster_filename = 'test16_clusters.txt'
  integer ( kind = 4 ), allocatable :: cluster_population(:)
  real ( kind = 8 ), allocatable :: cluster_variance(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  character ( len = 255 ) :: point_filename = 'points_100.txt'
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: weight(:)
  integer ( kind = 4 ) weight_dim
  integer ( kind = 4 ) weight_num
  character ( len = 255 ) :: weight_filename = 'weights_unequal_100.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  Test the KMEANS_W_03 algorithm.'
!
!  Read the data points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points will be read from "' // trim ( point_filename ) // '".'

  call r8mat_header_read ( point_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Point spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( point(dim_num,point_num) )

  call r8mat_data_read ( point_filename, dim_num, point_num, point )
!
!  Read the weights.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weights will be read from "' // trim ( weight_filename ) // '".'

  call r8mat_header_read ( weight_filename, weight_dim, weight_num )

  if ( weight_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of weight array is not 1.'
    stop
  end if

  if ( weight_num /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Number of weights not equal to number of points.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( weight_filename, weight_dim, weight_num, weight )
!
!  Clustering parameters.
!
  cluster_num = 5
  it_max = 20
  seed = 123456789

  allocate ( cluster(point_num) )
  allocate ( cluster_center(dim_num,cluster_num) )
  allocate ( cluster_energy(cluster_num) )
  allocate ( cluster_population(cluster_num) )
  allocate ( cluster_variance(cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations allowed is ', it_max
!
!  Initialize the cluster centers.
!
  call cluster_initialize_1 ( dim_num, point_num, cluster_num, point, &
    cluster_center )

  cluster(1:cluster_num) = 0

  call kmeans_w_03 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, &
    cluster_energy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken is ', it_num

  call cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

  call cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

  call r8mat_write ( center_filename, dim_num, cluster_num, &
    cluster_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Cluster centers written to "' // trim ( center_filename ) // '".'

  call i4mat_write ( cluster_filename, 1, point_num, cluster )

  write ( *, '(a)' ) &
    '  Cluster assignments written to "' // trim ( cluster_filename ) // '".'

  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( cluster_energy )
  deallocate ( cluster_population )
  deallocate ( cluster_variance )
  deallocate ( point )
  deallocate ( weight )

  return
end
