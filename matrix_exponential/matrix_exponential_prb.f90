program main

!*****************************************************************************80
!
!! MATRIX_EXPONENTIAL_TEST tests some matrix exponential algorithms.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the MATRIX_EXPONENTIAL library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  This test needs the TEST_MATRIX_EXPONENTIAL library.'

  call matrix_exponential_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine matrix_exponential_test01 ( )

!*****************************************************************************80
!
!! MATRIX_EXPONENTIAL_TEST01 compares matrix exponential algorithms.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: a_exp(:,:)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST01:'
  write ( *, '(a)' ) '  EXPM is MATLAB''s matrix exponential function'
  write ( *, '(a)' ) '  EXPM1 is an M-file equivalent to EXPM'
  write ( *, '(a)' ) '  EXPM2 uses a Taylor series approach'
  write ( *, '(a)' ) '  EXPM3 relies on an eigenvalue calculation.'

  call mexp_test_num ( test_num )

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Test #', test

    call mexp_story ( test )

    call mexp_n ( test, n )

    write ( *, '(a,i4)' ) '  Matrix order N = ', n

    allocate ( a(1:n,1:n) )

    call mexp_a ( test, n, a )

    call r8mat_print ( n, n, a, '  Matrix:' )

    allocate ( a_exp(1:n,1:n) )

    call expm1 ( n, a, a_exp )
    call r8mat_print ( n, n, a_exp, '  EXPM1(A):' )

    call expm2 ( n, a, a_exp )
    call r8mat_print ( n, n, a_exp, '  EXPM2(A):' )

!   call expm3 ( n, a, a_exp )
!   call r8mat_print ( n, n, a_exp, '  EXPM3(A):' )

    call mexp_expa ( test, n, a_exp )
    call r8mat_print ( n, n, a_exp, '  Exact Exponential:' )

    deallocate ( a )
    deallocate ( a_exp )

  end do

  return
end
