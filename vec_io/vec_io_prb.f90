program main

!*****************************************************************************80
!
!! MAIN is the main program for VEC_IO_PRB.
!
!  Discussion:
!
!    VEC_IO_PRB is tests the VEC_IO library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VEC_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the VEC_IO library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VEC_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests the R8UD_IO routine.
!
!  Discussion:
!
!    It is assumed that the user wants to store and retrieve a number
!    of vectors.  All the vectors are of the same size, and the user
!    always specifies an index or record number when writing or retrieving
!    a particular vector.  At any time, the user can write a new vector,
!    or retrieve any vector that has been written to the file earlier.
!
!    The data type is "R8" (double precision / real ( kind = 8 ) );
!    The storage format is "U" ( unformatted or binary );
!    and the access method is "D" ( direct access ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 100

  character action
  character ( len = 80 ) :: file_name = 'r8ud_io_test01.dat'
  integer file_unit
  integer i
  integer record
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  R8UD_IO stores and retrieves data vectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data type is "double precision", (R8)'
  write ( *, '(a)' ) '  The storage format is "unformatted", (U)'
  write ( *, '(a)' ) '  The access mode is "direct" (D).'
!
!  Create the file.
!
  action = 'c'
  call r8ud_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  File "' // trim ( file_name ) &
    // '" has been assigned I/O unit ', file_unit
!
!  Write the records in order.
!
!  (We could actually write them in any order,
!  we could skip some records, and write some records several times...)
!
  action = 'w'

  write ( *, '(a)' ) ' '

  do record = 1, 10

    call random_number ( harvest = x(1:n) )
    call r8ud_io ( action, file_name, file_unit, record, n, x )

    write ( *, '(a,i6,a,g14.6)' ) &
      '  Wrote record ', record, ' for which X(1) = ', x(1)

  end do
!
!  Read the records in random order.
!
  write ( *, '(a)' ) ' '

  action = 'r'
  do i = 1, 10
    call random_number ( harvest = t )
    record = 1 + int ( 10.0D+00 * t )
    call r8ud_io ( action, file_name, file_unit, record, n, x )
    write ( *, '(a,i6,a,g14.6)' ) &
      '  Read record  ', record, ' for which X(1) = ', x(1)
  end do
!
!  Reads and writes can happen in any order.
!  To prove this, we will:
!
!    Read records 6 and 7 and print 1 entry.
!
!    Write a new version of record 7.
!
!    Read records 6 (unchanged) and record 7 (now changed).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate that you can change some records'
  write ( *, '(a)' ) '  and others will be left alone.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read records 6 and 7,'
  write ( *, '(a)' ) ' '

  action = 'r'
  record = 6
  call r8ud_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Read record  ', record, ' for which X(1) = ', x(1)

  action = 'r'
  record = 7
  call r8ud_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Read record  ', record, ' for which X(1) = ', x(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Write a new version of record 7,'
  write ( *, '(a)' ) ' '

  call random_number ( harvest = x(1:n) )

  action = 'w'
  record = 7
  call r8ud_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Wrote record ', record, ' for which X(1) = ', x(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read records 6 and 7 again.'
  write ( *, '(a)' ) '  Record 6 should be the same, record 7 changed.'
  write ( *, '(a)' ) ' '

  action = 'r'
  record = 6
  call r8ud_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Read record  ', record, ' for which X(1) = ', x(1)

  action = 'r'
  record = 7
  call r8ud_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Read record  ', record, ' for which X(1) = ', x(1)
!
!  Get statistics.
!
  action = 's'
  call r8ud_io ( action, file_name, file_unit, record, n, x )
!
!  Delete the file.
!
  action = 'd'
  call r8ud_io ( action, file_name, file_unit, record, n, x )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests the R8UD_IO routine.
!
!  Discussion:
!
!    It is assumed that the user wants to store and retrieve a number
!    of vectors.  All the vectors are of the same size, and the user
!    always specifies an index or record number when writing or retrieving
!    a particular vector.  At any time, the user can write a new vector,
!    or retrieve any vector that has been written to the file earlier.
!
!    The data type is "R8" (double precision / real ( kind = 8 ) );
!    The storage format is "U" ( unformatted or binary );
!    and the access method is "D" ( direct access ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n1 = 10
  integer, parameter :: n2 = 7
  integer, parameter :: n3 = 12

  character action
  character ( len = 80 ) :: file_name1 = 'r8ud_io_test02_1.dat'
  character ( len = 80 ) :: file_name2 = 'r8ud_io_test02_2.dat'
  character ( len = 80 ) :: file_name3 = 'r8ud_io_test02_3.dat'
  integer file_unit1
  integer file_unit2
  integer file_unit3
  integer i
  integer record
  real ( kind = 8 ) t
  real ( kind = 8 ) x1(n1)
  real ( kind = 8 ) x2(n2)
  real ( kind = 8 ) x3(n3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R8UD_IO stores and retrieves data vectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data type is "double precision", (R8)'
  write ( *, '(a)' ) '  The storage format is "unformatted", (U)'
  write ( *, '(a)' ) '  The access mode is "direct" (D).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we create and manipulate'
  write ( *, '(a)' ) '  three files at the same time.'
!
!  Create file 1.
!
  action = 'c'
  call r8ud_io ( action, file_name1, file_unit1, record, n1, x1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  File "' // trim ( file_name1 ) &
    // '" has been assigned I/O unit ', file_unit1
!
!  Create file 2.
!
  action = 'c'
  call r8ud_io ( action, file_name2, file_unit2, record, n2, x2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  File "' // trim ( file_name2 ) &
    // '" has been assigned I/O unit ', file_unit2
!
!  Write the records in order.
!
!  (We could actually write them in any order,
!  we could skip some records, and write some records several times...)
!
  action = 'w'

  write ( *, '(a)' ) ' '

  do record = 1, 10

    call random_number ( harvest = x1(1:n1) )
    call r8ud_io ( action, file_name1, file_unit1, record, n1, x1 )

    write ( *, '(a,i6,a,i6,a,g14.6)' ) &
      '  Wrote record ', record,      &
      ' to unit ', file_unit1,        &
      ' for which X1(1) = ', x1(1)

  end do
!
!  Read the records in random order.
!
  write ( *, '(a)' ) ' '

  action = 'r'
  do i = 1, 10
    call random_number ( harvest = t )
    record = 1 + int ( 10.0D+00 * t )
    call r8ud_io ( action, file_name1, file_unit1, record, n1, x1 )
    write ( *, '(a,i6,a,g14.6)' ) &
      '  Read record  ', record, ' for which X1(1) = ', x1(1)
  end do
!
!  Create file 3.
!
  action = 'c'
  call r8ud_io ( action, file_name3, file_unit3, record, n3, x3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  File "' // trim ( file_name3 ) &
    // '" has been assigned I/O unit ', file_unit3
!
!  Create and write a record to each file.
!
  write ( *, '(a)' ) ' '

  action = 'w'

  record = 1
  call random_number ( harvest = x1(1:n1) )
  call r8ud_io ( action, file_name1, file_unit1, record, n1, x1 )

  write ( *, '(a,i6,a,i6,a,g14.6)' ) &
    '  Wrote record ', record,      &
    ' to unit ', file_unit1,        &
    ' for which X1(1) = ', x1(1)

  record = 1
  call random_number ( harvest = x2(1:n2) )
  call r8ud_io ( action, file_name2, file_unit2, record, n2, x2 )

  write ( *, '(a,i6,a,i6,a,g14.6)' ) &
    '  Wrote record ', record,      &
    ' to unit ', file_unit2,        &
    ' for which X2(1) = ', x2(1)

  record = 2
  call random_number ( harvest = x1(1:n1) )
  call r8ud_io ( action, file_name1, file_unit1, record, n1, x1 )

  write ( *, '(a,i6,a,i6,a,g14.6)' ) &
    '  Wrote record ', record,      &
    ' to unit ', file_unit1,        &
    ' for which X1(1) = ', x1(1)

  record = 1
  call random_number ( harvest = x3(1:n3) )
  call r8ud_io ( action, file_name3, file_unit3, record, n3, x3 )

  write ( *, '(a,i6,a,i6,a,g14.6)' ) &
    '  Wrote record ', record,      &
    ' to unit ', file_unit3,        &
    ' for which X3(1) = ', x3(1)

  record = 2
  call random_number ( harvest = x2(1:n2) )
  call r8ud_io ( action, file_name2, file_unit2, record, n2, x2 )

  write ( *, '(a,i6,a,i6,a,g14.6)' ) &
    '  Wrote record ', record,      &
    ' to unit ', file_unit2,        &
    ' for which X2(1) = ', x2(1)

  record = 3
  call random_number ( harvest = x1(1:n1) )
  call r8ud_io ( action, file_name1, file_unit1, record, n1, x1 )

  write ( *, '(a,i6,a,i6,a,g14.6)' ) &
    '  Wrote record ', record,      &
    ' to unit ', file_unit1,        &
    ' for which X1(1) = ', x1(1)
!
!  Recover and verify the data we just wrote to the different files.
!
  action = 'r'

  do record = 1, 3

    call r8ud_io ( action, file_name1, file_unit1, record, n1, x1 )

    write ( *, '(a,i6,a,i6,a,g14.6)' ) &
      '  Read record ', record,      &
      ' from unit ', file_unit1,        &
      ' for which X1(1) = ', x1(1)

  end do

  do record = 1, 2

    call r8ud_io ( action, file_name2, file_unit2, record, n2, x2 )

    write ( *, '(a,i6,a,i6,a,g14.6)' ) &
      '  Read record ', record,      &
      ' from unit ', file_unit2,        &
      ' for which X2(1) = ', x2(1)

  end do

  do record = 1, 1

    call r8ud_io ( action, file_name3, file_unit3, record, n3, x3 )

    write ( *, '(a,i6,a,i6,a,g14.6)' ) &
      '  Read record ', record,      &
      ' from unit ', file_unit3,        &
      ' for which X3(1) = ', x3(1)

  end do
!
!  Delete the files.
!
  action = 'd'
  call r8ud_io ( action, file_name1, file_unit1, record, n1, x1 )
  call r8ud_io ( action, file_name2, file_unit2, record, n2, x2 )
  call r8ud_io ( action, file_name3, file_unit3, record, n3, x3 )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests the R8AS_IO routine.
!
!  Discussion:
!
!    It is assumed that the user wants to store and retrieve a number
!    of vectors.  All the vectors are of the same size, and the user
!    always specifies an index or record number when writing or retrieving
!    a particular vector.  At any time, the user can write a new vector,
!    or retrieve any vector that has been written to the file earlier.
!
!    The data type is "R8" (double precision / real ( kind = 8 ) );
!    The storage format is "A" ( formatted or ASCII );
!    and the access method is "D" ( direct access ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 100

  character action
  character ( len = 80 ) :: file_name = 'r8ad_io_test03.dat'
  integer file_unit
  integer i
  integer record
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  R8AD_IO stores and retrieves data vectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data type is "double precision", (R8)'
  write ( *, '(a)' ) '  The storage format is "ASCII", (A)'
  write ( *, '(a)' ) '  The access mode is "direct" (D).'
!
!  Create the file.
!
  action = 'c'
  call r8ad_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  File "' // trim ( file_name ) &
    // '" has been assigned I/O unit ', file_unit
!
!  Write the records in order.
!
!  (We could actually write them in any order,
!  we could skip some records, and write some records several times...)
!
  action = 'w'

  write ( *, '(a)' ) ' '

  do record = 1, 10

    call random_number ( harvest = x(1:n) )
    call r8ad_io ( action, file_name, file_unit, record, n, x )

    write ( *, '(a,i6,a,g14.6)' ) &
      '  Wrote record ', record, ' for which X(1) = ', x(1)

  end do
!
!  Read the records in random order.
!
  write ( *, '(a)' ) ' '

  action = 'r'
  do i = 1, 10
    call random_number ( harvest = t )
    record = 1 + int ( 10.0D+00 * t )
    call r8ad_io ( action, file_name, file_unit, record, n, x )
    write ( *, '(a,i6,a,g14.6)' ) &
      '  Read record  ', record, ' for which X(1) = ', x(1)
  end do
!
!  Reads and writes can happen in any order.
!  To prove this, we will:
!
!    Read records 6 and 7 and print 1 entry.
!
!    Write a new version of record 7.
!
!    Read records 6 (unchanged) and record 7 (now changed).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate that you can change some records'
  write ( *, '(a)' ) '  and others will be left alone.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read records 6 and 7,'
  write ( *, '(a)' ) ' '

  action = 'r'
  record = 6
  call r8ad_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Read record  ', record, ' for which X(1) = ', x(1)

  action = 'r'
  record = 7
  call r8ad_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Read record  ', record, ' for which X(1) = ', x(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Write a new version of record 7,'
  write ( *, '(a)' ) ' '

  call random_number ( harvest = x(1:n) )

  action = 'w'
  record = 7
  call r8ad_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Wrote record ', record, ' for which X(1) = ', x(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read records 6 and 7 again.'
  write ( *, '(a)' ) '  Record 6 should be the same, record 7 changed.'
  write ( *, '(a)' ) ' '

  action = 'r'
  record = 6
  call r8ad_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Read record  ', record, ' for which X(1) = ', x(1)

  action = 'r'
  record = 7
  call r8ad_io ( action, file_name, file_unit, record, n, x )

  write ( *, '(a,i6,a,g14.6)' ) &
    '  Read record  ', record, ' for which X(1) = ', x(1)
!
!  Get statistics.
!
  action = 's'
  call r8ad_io ( action, file_name, file_unit, record, n, x )
!
!  Delete the file.
!
  action = 'd'
  call r8ad_io ( action, file_name, file_unit, record, n, x )

  return
end
