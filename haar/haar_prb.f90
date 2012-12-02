program main

!*****************************************************************************80
!
!! MAIN is the main program for HAAR_PRB.
!
!  Discussion:
!
!    HAAR_PRB tests the HAAR library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HAAR_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HAAR library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HAAR_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests HAAR_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: w(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  HAAR_1D computes the Haar transform of a vector.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )
  v(1:n) = u(1:n)

  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, u(i), v(i), w(i)
  end do

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )
!
!  Constant signal.
!
  n = 8
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  u(1:n) = 1.0D+00
  v(1:n) = u(1:n)

  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, u(i), v(i), w(i)
  end do

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )
!
!  Linear signal.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  a_first = 1.0D+00
  a_last = dble ( n )
  call r8vec_linspace ( n, a_first, a_last, u )
  v(1:n) = u(1:n)

  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, u(i), v(i), w(i)
  end do

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )
!
!  Quadratic data.
!
  n = 8
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  u(1) = 25.0D+00
  u(2) = 16.0D+00
  u(3) = 9.0D+00
  u(4) = 4.0D+00
  u(5) = 1.0D+00
  u(6) = 0.0D+00
  u(7) = 1.0D+00
  u(8) = 4.0D+00
  v(1:n) = u(1:n)

  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, u(i), v(i), w(i)
  end do

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HAAR_2D and HAAR_2D_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 16
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) seed
  real ( kind = 8 ) u(m,n)
  real ( kind = 8 ) v(m,n)
  real ( kind = 8 ) w(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HAAR_2D computes the Haar transform of an array.'
  write ( *, '(a)' ) '  HAAR_2D_INVERSE inverts the transform.'
!
!  Demonstrate successful inversion.
!
  seed = 123456789
  call r8mat_uniform_01 ( m, n, seed, u )

  call r8mat_print ( m, n, u, '  Input array U:' )

  v(1:m,1:n) = u(1:m,1:n)
  call haar_2d ( m, n, v )

  call r8mat_print ( m, n, v, '  Transformed array V:' )

  w(1:m,1:n) = v(1:m,1:n)
  call haar_2d_inverse ( m, n, w )

  call r8mat_print ( m, n, w, '  Recovered array W:' )

  return
end
