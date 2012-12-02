program main

!*****************************************************************************80
!
!! MAIN is the main program for WAVELET_PRB.
!
!  Discussion:
!
!    WAVELET_PRB calls the WAVELET library.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WAVELET_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the WAVELET library.'
!
!  Test transforms.
!
  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WAVELET_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DAUB2_TRANSFORM and DAUB2_TRANSFORM_INVERSE.
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  DAUB2_TRANSFORM computes the DAUB2 transform of a vector.'
  write ( *, '(a)' ) '  DAUB2_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub2_transform ( n, u, v )

  call daub2_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D2(U)    D2inv(D2(U))'
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

  call daub2_transform ( n, u, v )

  call daub2_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D2(U)    D2inv(D2(U))'
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

  call daub2_transform ( n, u, v )

  call daub2_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D2(U)    D2inv(D2(U))'
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

  call daub2_transform ( n, u, v )

  call daub2_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D2(U)    D2inv(D2(U))'
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
!! TEST02 tests DAUB4_TRANSFORM and DAUB4_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  DAUB4_TRANSFORM computes the DAUB4 transform of a vector.'
  write ( *, '(a)' ) '  DAUB4_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub4_transform ( n, u, v )

  call daub4_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D4(U)    D4inv(D4(U))'
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

  call daub4_transform ( n, u, v )

  call daub4_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D4(U)    D4inv(D4(U))'
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

  call daub4_transform ( n, u, v )

  call daub4_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D4(U)    D4inv(D4(U))'
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

  call daub4_transform ( n, u, v )

  call daub4_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D4(U)    D4inv(D4(U))'
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
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DAUB6_TRANSFORM and DAUB6_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  DAUB6_TRANSFORM computes the DAUB6 transform of a vector.'
  write ( *, '(a)' ) '  DAUB6_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub6_transform ( n, u, v )

  call daub6_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D6(U)    D6inv(D6(U))'
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

  call daub6_transform ( n, u, v )

  call daub6_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D6(U)    D6inv(D6(U))'
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

  call daub6_transform ( n, u, v )

  call daub6_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D6(U)    D6inv(D6(U))'
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

  call daub6_transform ( n, u, v )

  call daub6_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D6(U)    D6inv(D6(U))'
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
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DAUB8_TRANSFORM and DAUB8_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DAUB8_TRANSFORM computes the DAUB8 transform of a vector.'
  write ( *, '(a)' ) '  DAUB8_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub8_transform ( n, u, v )

  call daub8_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D8(U)    D8inv(D8(U))'
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

  call daub8_transform ( n, u, v )

  call daub8_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D8(U)    D8inv(D8(U))'
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

  call daub8_transform ( n, u, v )

  call daub8_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D8(U)    D8inv(D8(U))'
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

  call daub8_transform ( n, u, v )

  call daub8_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U          D8(U)    D8inv(D8(U))'
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
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests DAUB10_TRANSFORM and DAUB10_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  DAUB10_TRANSFORM computes the DAUB10 transform of a vector.'
  write ( *, '(a)' ) '  DAUB10_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub10_transform ( n, u, v )

  call daub10_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D10(U)   D10inv(D10(U))'
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

  call daub10_transform ( n, u, v )

  call daub10_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D10(U)   D10inv(D10(U))'
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

  call daub10_transform ( n, u, v )

  call daub10_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D10(U)   D10inv(D10(U))'
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

  call daub10_transform ( n, u, v )

  call daub10_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D10(U)   D10inv(D10(U))'
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
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests DAUB12_TRANSFORM and DAUB12_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  DAUB12_TRANSFORM computes the DAUB12 transform of a vector.'
  write ( *, '(a)' ) '  DAUB12_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub12_transform ( n, u, v )

  call daub12_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D12(U)   D12inv(D12(U))'
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

  call daub12_transform ( n, u, v )

  call daub12_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D12(U)   D12inv(D12(U))'
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

  call daub12_transform ( n, u, v )

  call daub12_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D12(U)   D12inv(D12(U))'
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

  call daub12_transform ( n, u, v )

  call daub12_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D12(U)   D12inv(D12(U))'
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
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests DAUB14_TRANSFORM and DAUB14_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  DAUB14_TRANSFORM computes the DAUB14 transform of a vector.'
  write ( *, '(a)' ) '  DAUB14_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub14_transform ( n, u, v )

  call daub14_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D14(U)   D14inv(D14(U))'
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

  call daub14_transform ( n, u, v )

  call daub14_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D14(U)   D14inv(D14(U))'
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

  call daub14_transform ( n, u, v )

  call daub14_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D14(U)   D14inv(D14(U))'
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

  call daub14_transform ( n, u, v )

  call daub14_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D14(U)   D14inv(D14(U))'
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
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests DAUB16_TRANSFORM and DAUB16_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  DAUB16_TRANSFORM computes the DAUB16 transform of a vector.'
  write ( *, '(a)' ) '  DAUB16_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub16_transform ( n, u, v )

  call daub16_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D16(U)   D16inv(D16(U))'
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

  call daub16_transform ( n, u, v )

  call daub16_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D16(U)   D16inv(D16(U))'
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

  call daub16_transform ( n, u, v )

  call daub16_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D16(U)   D16inv(D16(U))'
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

  call daub16_transform ( n, u, v )

  call daub16_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D16(U)   D16inv(D16(U))'
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
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests DAUB18_TRANSFORM and DAUB18_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  DAUB18_TRANSFORM computes the DAUB18 transform of a vector.'
  write ( *, '(a)' ) '  DAUB18_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub18_transform ( n, u, v )

  call daub18_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D18(U)   D18inv(D18(U))'
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

  call daub18_transform ( n, u, v )

  call daub18_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D18(U)   D18inv(D18(U))'
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

  call daub18_transform ( n, u, v )

  call daub18_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D18(U)   D18inv(D18(U))'
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

  call daub18_transform ( n, u, v )

  call daub18_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D18(U)   D18inv(D18(U))'
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
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests DAUB20_TRANSFORM and DAUB20_TRANSFORM_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2011
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  DAUB20_TRANSFORM computes the DAUB20 transform of a vector.'
  write ( *, '(a)' ) '  DAUB20_TRANSFORM_INVERSE inverts it.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, u )

  call daub20_transform ( n, u, v )

  call daub20_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D20(U)   D20inv(D20(U))'
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

  call daub20_transform ( n, u, v )

  call daub20_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D20(U)   D20inv(D20(U))'
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

  call daub20_transform ( n, u, v )

  call daub20_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D20(U)   D20inv(D20(U))'
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

  call daub20_transform ( n, u, v )

  call daub20_transform_inverse ( n, v, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U         D20(U)   D20inv(D20(U))'
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
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests DAUB*_MATRIX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: b(:,:)
  real ( kind = 8 ), allocatable :: c(:,:)
  real ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  DAUB*_MATRIX computes the DAUB* matrix.'
  write ( *, '(a)' ) '  Verify that each matrix is orthogonal.'

  n = 20

  allocate ( a(n,n) )
  allocate ( b(n,n) )
  allocate ( c(n,n) )

  call daub2_matrix ( n, a )

  b = transpose ( a )
  c = matmul ( a, b )
  call r8mat_is_identity ( n, c, error_frobenius )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a,g14.6)' ) '  DAUB2,  N = ', n, '  || A*A'' - I|| = ', error_frobenius

  call daub4_matrix ( n, a )

  b = transpose ( a )
  c = matmul ( a, b )
  call r8mat_is_identity ( n, c, error_frobenius )

  write ( *, '(a,i4,a,g14.6)' ) '  DAUB4,  N = ', n, '  || A*A'' - I|| = ', error_frobenius

  call daub6_matrix ( n, a )

  b = transpose ( a )
  c = matmul ( a, b )
  call r8mat_is_identity ( n, c, error_frobenius )

  write ( *, '(a,i4,a,g14.6)' ) '  DAUB6,  N = ', n, '  || A*A'' - I|| = ', error_frobenius

  call daub8_matrix ( n, a )

  b = transpose ( a )
  c = matmul ( a, b )
  call r8mat_is_identity ( n, c, error_frobenius )

  write ( *, '(a,i4,a,g14.6)' ) '  DAUB8,  N = ', n, '  || A*A'' - I|| = ', error_frobenius

  call daub10_matrix ( n, a )

  b = transpose ( a )
  c = matmul ( a, b )
  call r8mat_is_identity ( n, c, error_frobenius )

  write ( *, '(a,i4,a,g14.6)' ) '  DAUB10, N = ', n, '  || A*A'' - I|| = ', error_frobenius

  call daub12_matrix ( n, a )

  b = transpose ( a )
  c = matmul ( a, b )
  call r8mat_is_identity ( n, c, error_frobenius )

  write ( *, '(a,i4,a,g14.6)' ) '  DAUB12, N = ', n, '  || A*A'' - I|| = ', error_frobenius

  deallocate ( a )
  deallocate ( b )
  deallocate ( c )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests DAUB*_SCALE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) daub2_scale
  real ( kind = 8 ) daub4_scale
  real ( kind = 8 ) daub6_scale
  real ( kind = 8 ) daub8_scale
  real ( kind = 8 ) daub10_scale
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ) x_save
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  DAUB*_SCALE uses recursion to estimate a scaling function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N        X       D2       D4       D6       D8      D10'

  seed = 123456789

  do test = 1, 5
    write ( *, '(a)' ) ' '
    x_save = r8_uniform_01 ( seed )
    do n = 1, 6
      x = x_save
      y1 = daub2_scale ( n, x )
      x = x_save
      y2 = daub4_scale ( n, x )
      x = x_save
      y3 = daub6_scale ( n, x )
      x = x_save
      y4 = daub8_scale ( n, x )
      x = x_save
      y5 = daub10_scale ( n, x )
      write ( *, '(2x,i4,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4)' ) &
        n, x, y1, y2, y3, y4, y5
    end do
  end do

  return
end
