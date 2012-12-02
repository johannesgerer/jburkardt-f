program main

!*****************************************************************************80
!
!! TEST_MATRIX_EXPONENTIAL_TEST tests some matrix exponential algorithms.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( );
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MATRIX_EXPONENTIAL_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_MATRIX_EXPONENTIAL library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'

  call test_matrix_exponential_test01 ( );
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MATRIX_EXPONENTIAL_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( );

  return
end
subroutine test_matrix_exponential_test01 ( )

!*****************************************************************************80
!
!! TEST_MATRIX_EXPONENTIAL_TEST01 retrieves the test data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: expa(:,:)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MATRIX_EXPONENTIAL_TEST01:'
  write ( *, '(a)' ) '  Retrieve the data for each matrix exponential test.'

  call mexp_test_num ( test_num )

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Test #', test

    call mexp_n ( test, n )

    call mexp_story ( test );

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Matrix order N = ', n

    allocate ( a(1:n,1:n) )

    call mexp_a ( test, n, a )

    call r8mat_print ( n, n, a, '  Matrix A:' );

    allocate ( expa(1:n,1:n) )

    call mexp_expa ( test, n, expa )
    call r8mat_print ( n, n, expa, '  Exact Exponential exp(A):' )

    deallocate ( a )
    deallocate ( expa )

  end do

  return
end
