program main

!*****************************************************************************80
!
!! SPARSEKIT_PRB14 demonstrates CSR to NCF conversion.
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB14'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This test demonstrates the conversion of a'
  write ( *, '(a)' ) '  sparse matrix from CSR to NCF formats.'

  call test01

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB14'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 sets up a CSR matrix and converts it to NCF form.
!
!  Discussion:
!
!    The matrix is
!
!    11  12   0  14   0
!    21  22  23   0  25
!     0  32  33  34   0
!    41   0  43  44  45
!     0  52   0  54  55
!
!    The CSR format is
!
!    A = ( 11, 12, 14, 21, 22, 23, 25, 32, 33, 34, 41, 43, 44, 45, 52, 54, 55 )
!    
!    IA = ( 1,          4,              8,         11,             15,
!          18 )
!
!    JA = ( 1,  2,  4,  1,  2,  3,  5,  2,  3,  4,  1,  3,  4,  5,  2, 4,  5 )
!
!    The NCF format is
!
!    COEF =   ( 11, 22, 33, 44, 55, 12, 14, 21, 23, 25, 32, 34, 41, 43, 45, 52, 54 )
!    JCOEF' = ( 1,  2,  3,  4,  5,  1,  1,  2,  2,  2,  3,  3,  4,  4,  4,  5,  5  )
!             (  1,  2,  3,  4,  5,  2,  4,  1,  3,  5,  2,  4,  1,  3,  5,  2,  4 )
!
  implicit none

  integer, parameter :: n = 5
  integer, parameter :: nzmax = 17

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) coef(nzmax)
  integer i
  integer ia(n+1)
  integer ierr
  integer ja(nzmax)
  integer jcoef(nzmax,2)
  integer nz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Set up a CSR matrix, and call'
  write ( *, '(a)' ) '  CSRNCF to convert it to NCF format.'

  a(1:nzmax) = (/ &
    11.0D+00, 12.0D+00,           14.0D+00,           &
    21.0D+00, 22.0D+00, 23.0D+00,           25.0D+00, &
              32.0D+00, 33.0D+00, 34.0D+00,           &
    41.0D+00,           43.0D+00, 44.0D+00, 45.0D+00, &
              52.0D+00,           54.0D+00, 55.0D+00 /)

  ia(1:n+1) = (/ 1, 4, 8, 11, 15, 18 /)

  ja(1:nzmax) = (/ &
    1, 2,    4,    &
    1, 2, 3,    5, &
       2, 3, 4,    &
    1,    3, 4, 5, & 
       2, 4, 5 /)

  call csrncf ( n, a, ja, ia, nzmax, nz, coef, jcoef, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  CSRNCF returned IERR = ', ierr
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I       J        A(I,J)'
  write ( *, '(a)' ) ' '

  do i = 1, nz
    write ( *, '(2x,i6,2x,i6,2x,g14.6)' ) jcoef(i,1), jcoef(i,2), coef(i)
  end do

  return
end
