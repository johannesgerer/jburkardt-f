program main

!*****************************************************************************80
!
!! LAPACK_PRB shows how to solve A*x=b with the LAPACK routine SGESV.
!
!  Discussion:
!
!    The exact solution to the given system is
!
!      x = ( 1, 0, 0, 0, ..., 0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2002
!
!  Author:
!
!    John Burkardt
!
  integer, parameter :: n = 5
  integer, parameter :: nrhs = 1

  integer, parameter :: lda = n
  integer, parameter :: ldb = n

  real a(lda,n)
  real b(ldb)
  integer i
  integer info
  integer ipiv(n)
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAPACK_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the use of LAPACK to'
  write ( *, '(a)' ) '  solve a simple linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The correct solution is '
  write ( *, '(a)' ) '  ( 1, 0, 0, ... 0 ).'

  b(1:n) = 1.0E+00

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( min ( i, j ) )
    end do
  end do

  call sgesv ( n, nrhs, a, n, ipiv, b, ldb, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SGESV returns solution vector X:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(g14.6)' ) b(i)
    end do
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  SGESV returns error flag INFO = ', info
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end

