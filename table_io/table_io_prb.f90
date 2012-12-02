program main

!*****************************************************************************80
!
!! MAIN is the main program for TABLE_IO_PRB.
!
!  Discussion:
!
!    TABLE_IO_PRB calls the TABLE_IO test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the TABLE_IO library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_IO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests R8MAT_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: n = 20
  integer   ( kind = 4 ), parameter :: m = 5

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  character ( len = 80 ) :: output_filename = 'r8mat_05_00020.txt'
  real      ( kind = 8 ) table(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) &
    '  R8MAT_WRITE0 writes an R8MAT file.'

  call r8mat_indicator ( m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i8)' ) '  Number of points N  = ', n

  call r8mat_print_some ( m, n, table, 1, 1, 5, 5, &
    '  5x5 portion of the data written to file:' )

  call r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, &
    '  5x5 portion of the TRANSPOSED data:' )

  call r8mat_write ( output_filename, m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the file "' &
    // trim ( output_filename ) //'".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R8MAT_HEADER_READ, R8MAT_DATA_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: input_filename = 'r8mat_05_00020.txt'
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: table

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For an R8MAT file,'
  write ( *, '(a)' ) '  R8MAT_HEADER_READ reads the header information'
  write ( *, '(a)' ) '  (about the dimensions of the data);'
  write ( *, '(a)' ) '  R8MAT_DATA_READ reads the data.'

  call r8mat_header_read (  input_filename, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' // trim ( input_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i8)' ) '  Number of points N  = ', n

  allocate ( table(1:m,1:n) )

  call r8mat_data_read (  input_filename, m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' // trim ( input_filename ) //'".'

  call r8mat_print_some ( m, n, table, 1, 1, 5, 5, &
    '  5x5 portion of data read from file:' )

  deallocate ( table )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests I4MAT_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: n = 20
  integer   ( kind = 4 ), parameter :: m = 5

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  character ( len = 80 ) :: output_filename = 'i4mat_05_00020.txt'
  integer   ( kind = 4 ) table(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  I4MAT_WRITE writes an I4MAT file.'

  call i4mat_indicator ( m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i8)' ) '  Number of points N  = ', n

  call i4mat_print_some ( m, n, table, 1, 1, 5, 5, &
    '  5 x 5 portion of data written to file:' )

  call i4mat_write ( output_filename, m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the file "' &
    // trim ( output_filename ) //'".'

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests I4MAT_HEADER_READ, I4MAT_DATA_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: input_filename = 'i4mat_05_00020.txt'
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: table

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For an I4MAT file,'
  write ( *, '(a)' ) '  I4MAT_HEADER_READ reads the header information'
  write ( *, '(a)' ) '  (about the dimensions of the data);'
  write ( *, '(a)' ) '  I4MAT_DATA_READ reads the data.'

  call i4mat_header_read (  input_filename, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' // trim ( input_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i8)' ) '  Number of points N  = ', n

  allocate ( table(1:m,1:n) )

  call i4mat_data_read (  input_filename, m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' // trim ( input_filename ) //'".'

  call i4mat_print_some ( m, n, table, 1, 1, 5, 5, &
    '  5 x 5 portion of data read from file:' )

  deallocate ( table )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests R8MAT_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ), dimension (m,n) :: table

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  R8MAT_UNIFORM_01 sets a random R8MAT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i8)' ) '  Number of points N  = ', n

  call r8mat_uniform_01 ( m, n, seed, table )

  call r8mat_print_some ( m, n, table, 1, 1, 5, 10, &
    '  5x10 portion of random real table dataset:' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests I4MAT_BORDER_ADD, I4MAT_BORDER_CUT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ), dimension (m,n) :: table
  integer ( kind = 4 ), dimension (m-2,n-2) :: table2
  integer ( kind = 4 ), dimension (m,n) :: table3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  I4MAT_BORDER_CUT cuts off the border;'
  write ( *, '(a)' ) '  I4MAT_BORDER_ADD adds a zero border.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i8)' ) '  Number of points N  = ', n

  call i4mat_indicator ( m, n, table )

  call i4mat_print ( m, n, table, '  Initial dataset:' )

  call i4mat_border_cut ( m, n, table, table2 )

  call i4mat_print ( m-2, n-2, table2, '  "Cut" dataset:' )

  call i4mat_border_add ( m-2, n-2, table2, table3 )

  call i4mat_print ( m, n, table3, '  "Added" dataset:' )

  return
end
