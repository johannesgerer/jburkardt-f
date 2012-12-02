program main

!*****************************************************************************80
!
!! MAIN is the main program for SPHERE_MONTE_CARLO_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_MONTE_CARLO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPHERE_MONTE_CARLO library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_MONTE_CARLO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses SPHERE_SAMPLE_01 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e(3)
  integer ( kind = 4 ) :: e_test(3,7) = reshape ( (/ &
    0, 0, 0, &
    2, 0, 0, &
    0, 2, 0, &
    0, 0, 2, &
    4, 0, 0, &
    2, 2, 0, &
    0, 0, 4 /), (/ 3, 7 /) )
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ) result(7)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X^2             Y^2' // &
    '             Z^2             X^4           X^2Y^2           Z^4'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:3,1:n) )

    call sphere01_sample ( n, seed, x )

    do j = 1, 7

      e(1:3) = e_test(1:3,j)

      call monomial_value ( 3, n, x, e, value )

      result(j) = 4.0D+00 * pi * sum ( value(1:n) ) / real ( n, kind = 8 )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 7

    e(1:3) = e_test(1:3,j)

    call sphere01_monomial_integral ( e, result(j) )

  end do

  write ( *, '(2x,a8,7(2x,g14.6))' ) '   Exact', result(1:7)

  return
end
