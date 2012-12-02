program main

!*****************************************************************************80
!
!! LATINIZE_PRB tests the LATINIZE routines.
!
!  Discussion:
!
!    The dataset is presumed to be an M by N array of real numbers,
!    where M is the spatial dimension, and N is the number of sample points.
!
!    The dataset is presumed to be stored in a file, with N records,
!    one per each sample point.  (Comment records may be included,
!    which begin with '#'.)
!
!    The program reads the data file, "latinizes" the data, and writes
!    the latinized data to a new file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATINIZE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LATINIZE library.'

  call test01 ( 'cvt_02_00010.txt' )
  call test01 ( 'cvt_03_00007.txt' )
  call test01 ( 'cvt_03_00056.txt' )
  call test01 ( 'cvt_07_00100.txt' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATINIZE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( input_filename )

!*****************************************************************************80
!
!! TEST01 reads a datafile, latinizes it, and writes the new data out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = * ) :: input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = 255 ) :: output_filename = ' '
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: table
!
!  Get the row and column dimensions of the dataset.
!
  call r8mat_header_read (  input_filename, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' // trim ( input_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i6)' ) '  Number of points N  = ', n
!
!  Set up space for the array.
!
  allocate ( table(1:m,1:n) )
!
!  Read the array.
!
  call r8mat_data_read ( input_filename, m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' // trim ( input_filename ) //'".'

  call r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, &
    '  Small portion of data read from file:' )
!
!  Latinize the array.
!
  call r8mat_latinize ( m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Latinized the data.'
!
!  Print out a small sample of the array.
!
  call r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, &
    '  Small portion of Latinized data:' )
!
!  Make up a name for the output file that is likely to be related
!  to the input file name, and unique.
!
  output_filename = input_filename
  call file_name_ext_swap ( output_filename, 'latin.txt' )
!
!  Write the latinized array to a file.
!
  call r8mat_write ( output_filename, m, n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the latinized data to "' &
    // trim ( output_filename ) //'".'

  deallocate ( table )

  return
end
