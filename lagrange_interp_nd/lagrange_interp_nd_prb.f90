program main

!*****************************************************************************80
!
!! TEST tests the LAGRANGE_INTERP_ND library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LAGRANGE_INTERP_ND library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
!
!  Use the interface that passes in the orders directly.
!
  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Repeat tests 1, 2, 3, and 4,
!  using the interface that passes in the orders indirectly,
!  based on the "level".
!
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Experiment with anisotropic orders.
!
  call test09 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 interpolates in 1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: n_1d(:)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Interpolate in 1D, using orders.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_GRID sets the interpolant.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE evaluates it.'

  m = 1

  allocate ( n_1d(1:m) )
  n_1d(1:m) = 5

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  call lagrange_interp_nd_size ( m, n_1d, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )

  call lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
  call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = 5
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )

  allocate ( zi(1:ni) )
  call lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j) )
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 interpolates in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: n_1d(:)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Interpolate in 2D, using orders.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_GRID sets the interpolant.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE evaluates it.'

  m = 2

  allocate ( n_1d(1:m) )
  n_1d(1:m) = 5

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  call lagrange_interp_nd_size ( m, n_1d, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )

  call lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
  call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = 5
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )

  allocate ( zi(1:ni) )
  call lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j) ) 
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 interpolates in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: n_1d(:)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  Interpolate in 3D, using orders.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_GRID sets the interpolant.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE evaluates it.'

  m = 3

  allocate ( n_1d(1:m) )
  n_1d(1:m) = 5

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  call lagrange_interp_nd_size ( m, n_1d, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )

  call lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
  call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = 5
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )

  allocate ( zi(1:ni) )
  call lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j) ) 
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 interpolates in 3D, using increasing resolution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ) e
  external f_sinr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: n_1d(:)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) order
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  Interpolate in 3D, using orders.'
  write ( *, '(a)' ) '  Use a sequence of increasing orders.'

  m = 3

  allocate ( n_1d(1:m) )

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  ni = 20
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( zi(1:ni) )

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Level     Order   Average Error'
  write ( *, '(a)' ) ''

  do l = 0, 5

    if ( l == 0 ) then
      order = 1
    else
      order = 2 ** l + 1
    end if
    n_1d(1:m) = order
  
    call lagrange_interp_nd_size ( m, n_1d, nd )

    allocate ( xd(1:m,1:nd) )
    allocate ( zd(1:nd) )

    call lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
    call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
    call lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi, zi )

    e = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = 8 )

    write ( *, '(2x,i5,2x,i8,2x,e10.2)' ) l, nd, e

    deallocate ( xd )
    deallocate ( zd )
    
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 repeats test 1 using levels.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  Repeat test #1, using levels.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE2 evaluates it.'

  m = 1

  allocate ( ind(1:m) )
  ind(1:m) = 2

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  call lagrange_interp_nd_size2 ( m, ind, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )

  call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
  call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = 5
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )

  allocate ( zi(1:ni) )
  call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j) ) 
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 repeats test 2 using levels.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  Repeat test #2, using levels.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE2 evaluates it.'

  m = 2

  allocate ( ind(1:m) )
  ind(1:m) = 2

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  call lagrange_interp_nd_size2 ( m, ind, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )

  call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
  call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = 5
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )

  allocate ( zi(1:ni) )
  call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j) ) 
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 repeats test 3 using levels.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  Repeat test #3,  using levels.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE2 evaluates it.'

  m = 3

  allocate ( ind(1:m) )
  ind(1:m) = 2

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  call lagrange_interp_nd_size2 ( m, ind, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )

  call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
  call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = 5
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )

  allocate ( zi(1:ni) )
  call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j) ) 
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 repeats test 4 using levels.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ) e
  external f_sinr
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) order
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  Interpolate in 3D, using levels.'
  write ( *, '(a)' ) '  Use a sequence of increasing levels.'

  m = 3

  allocate ( ind(1:m) )

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  ni = 20
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( zi(1:ni) )

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Level     Order   Average Error'
  write ( *, '(a)' ) ''

  do l = 0, 5

    ind(1:m) = l
  
    call lagrange_interp_nd_size2 ( m, ind, nd )

    allocate ( xd(1:m,1:nd) )
    allocate ( zd(1:nd) )

    call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
    call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
    call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

    e = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = 8 )

    write ( *, '(2x,i5,2x,i8,2x,e10.2)' ) l, nd, e

    deallocate ( xd )
    deallocate ( zd )

  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 interpolates in 3D, using anisotropic resolution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ) e
  external f_poly352
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: n_1d(:)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) order
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  Interpolate in 3D, using orders.'
  write ( *, '(a)' ) '  Use a sequence of increasing orders.'
  write ( *, '(a)' ) '  Use anisotropic resolution.'
  write ( *, '(a)' ) '  The interpoland is a polynomial of degrees 3, 5, 2'
  write ( *, '(a)' ) '  so our orders need to be at least 4, 6, 3 to match it.'

  m = 3

  allocate ( n_1d(1:m) )

  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  ni = 20
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )

  allocate ( zi(1:ni) )

  allocate ( ze(1:ni) )
  call f_poly352 ( m, ni, xi, ze )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Level     Orders   Average Error'
  write ( *, '(a)' ) ''

  do l = 0, 10

    if ( l == 0 ) then
      n_1d(1:m) = (/ 1, 1, 1 /)
    else if ( l == 1 ) then
      n_1d(1:m) = (/ 2, 1, 1 /)
    else if ( l == 2 ) then
      n_1d(1:m) = (/ 1, 2, 1 /)
    else if ( l == 3 ) then
      n_1d(1:m) = (/ 1, 1, 2 /)
    else if ( l == 4 ) then
      n_1d(1:m) = (/ 4, 2, 2 /)
    else if ( l == 5 ) then
      n_1d(1:m) = (/ 2, 4, 2 /)
    else if ( l == 6 ) then
      n_1d(1:m) = (/ 2, 2, 4 /)
    else if ( l == 8 ) then
      n_1d(1:m) = (/ 6, 4, 4 /)
    else if ( l == 9 ) then
      n_1d(1:m) = (/ 4, 6, 4 /)
    else if ( l == 10 ) then
      n_1d(1:m) = (/ 4, 4, 6 /)
    end if
  
    call lagrange_interp_nd_size ( m, n_1d, nd )

    allocate ( xd(1:m,1:nd) )
    allocate ( zd(1:nd) )

    call lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
    call f_poly352 ( m, nd, xd, zd )
!
!  Evaluate.
!
    call lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi, zi )

    e = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = 8 )

    write ( *, '(2x,i5,2x,i5,2x,i5,2x,i5,2x,e10.2)' ) l, n_1d(1:3), e

    deallocate ( xd )
    deallocate ( zd )

  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine f_sinr ( m, n, x, z )

!*****************************************************************************80
!
!! F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the points.
!
!    Output, real ( kind = 8 ) Z(N), the value of the function at each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)

  r(1:n) = sqrt ( sum ( x(1:m,1:n)**2, dim = 1 ) )

  z(1:n) = sin ( r(1:n) )

  return
end
subroutine f_poly352 ( m, n, x, z )

!*****************************************************************************80
!
!! F_POLY253 is a scalar function of a 3-dimensional argument, to be interpolated.
!
!  Discussion:
!
!    The polynomial is of maximum degrees 3, 5, and 2, in the three variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the points.
!
!    Output, real ( kind = 8 ) Z(N,1), the value of the function at each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)

  z(1:n) = 1.0 + x(1,1:n) ** 2 * x(2,1:n) ** 5 * x(3,1:n) ** 2 &
    + x(1,1:n) * x(2,1:n) ** 2 * x(3,1:n) + x(1,1:n) ** 3

  return
end
