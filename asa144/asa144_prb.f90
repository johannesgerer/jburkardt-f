program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA144_PRB.
!
!  Discussion:
!
!    ASA144_PRB calls a set of problems for ASA144.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA144_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA144 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA144_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests RCONT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nrow = 5
  integer ( kind = 4 ), parameter :: ncol = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  logical key
  integer ( kind = 4 ) matrix(nrow,ncol)
  integer ( kind = 4 ), dimension ( ncol ) :: ncolt = (/ &
    2, 2, 2, 2, 1 /)
  integer ( kind = 4 ), dimension ( nrow ) :: nrowt = (/ &
    3, 2, 2, 1, 1 /)
  integer ( kind = 4 ) nsubt(ncol)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  RCONT constructs a random matrix with'
  write ( *, '(a)' ) '  given row and column sums.'

  call i4vec_print ( nrow, nrowt, '  The rowsum vector:' )
  call i4vec_print ( ncol, ncolt, '  The columnsum vector: ' )

  key = .false.

  do test = 1, test_num

    call rcont ( nrow, ncol, nrowt, ncolt, nsubt, matrix, key, ifault )

    if ( ifault /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  RCONT returned IFAULT = ', ifault
      return
    end if

    call i4mat_print ( nrow, ncol, matrix, '  The rowcolsum matrix:' )

  end do

  return
end
