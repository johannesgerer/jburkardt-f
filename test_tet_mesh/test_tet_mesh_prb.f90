program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_TET_MESH_PRB.
!
!  Discussion:
!
!    TEST_TET_MESH_PRB runs the TEST_TET_MESH tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TET_MESH_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_TET_MESH library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TET_MESH_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests P00_TEST_NUM, P00_TITLE, and P00_HEADER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  P00_TEST_NUM reports the number of problems.'
  write ( *, '(a)' ) '  P00_TITLE returns a title for each problem.'
  write ( *, '(a)' ) '  P00_HEADER prints some information about each problem.'

  call p00_test_num ( test_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of tests available = ', test_num

  do test = 1, test_num

    call p00_title ( test, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Test number:        ', test
    write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

    call p00_header ( test )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests P00_TEST_NUM and P00_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_start = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  P00_TEST_NUM reports the number of problems.'
  write ( *, '(a)' ) '  P00_SAMPLE returns sample points from the region.'

  call p00_test_num ( test_num )

  do test = 1, test_num

    seed = seed_start

    call p00_title ( test, title )

    n = 20

    write ( *, '(a)'     ) ' '
    write ( *, '(a,i6)'  ) '  Test number       =    ', test
    write ( *, '(a, a)'  ) '  Title:            =   "', trim ( title ) // '"'
    write ( *, '(a,i6)'  ) '  Dimension DIM_NUM =    ', dim_num
    write ( *, '(a,i6)'  ) '  Number of samples N =  ', n
    write ( *, '(a,i12)' ) '  Initial SEED:       =  ', seed
    write ( *, '(a)'     ) ' '

    allocate ( point(1:dim_num,1:n) )

    call p00_sample ( test, n, seed, point )

    call r8mat_transpose_print ( dim_num, n, point, '  The sample points:' )

    deallocate ( point )

  end do

  return
end
