!! MPI_STUBS_F90.H is the include file for MPI_STUBS.F90.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2012
!
!  Author:
!
!    John Burkardt
!
  integer mpi_comm_world
  parameter ( mpi_comm_world = 0 )
!
!  Return values
!
  integer mpi_failure
  parameter ( mpi_failure = 1 )
  integer mpi_success
  parameter ( mpi_success = 0 )
!
!  recv message status
!
  integer mpi_status_size
  parameter ( mpi_status_size = 3 )
  integer mpi_source
  parameter ( mpi_source = 1 )
  integer mpi_tag
  parameter ( mpi_tag = 2 )
  integer mpi_count
  parameter ( mpi_count = 3 )
!
!  recv flags
!
  integer mpi_any_source
  parameter ( mpi_any_source = -1 )
  integer mpi_any_tag
  parameter ( mpi_any_tag = -1 )
!
!  data types and sizes
!
  integer mpi_integer
  parameter ( mpi_integer = 1 )
  integer mpi_real
  parameter ( mpi_real = 2 )
  integer mpi_double_precision
  parameter ( mpi_double_precision = 3 )
  integer mpi_logical
  parameter ( mpi_logical = 4 )
  integer mpi_character
  parameter ( mpi_character = 5 )
!
!  allreduce operations
!
  integer mpi_sum
  parameter ( mpi_sum = 1 )
  integer mpi_max
  parameter ( mpi_max = 2 )
  integer mpi_min
  parameter ( mpi_min = 3 )
  integer mpi_product
  parameter ( mpi_product = 4 )
!
!  timer
!
  double precision mpi_wtime
