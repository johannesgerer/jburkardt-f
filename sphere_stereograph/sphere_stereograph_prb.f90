program main

!*****************************************************************************80
!
!! MAIN is the main program for SPHERE_STEREOGRAPH_PRB.
!
!  Discussion:
!
!    SPHERE_STEREOGRAPH_PRB calls the SPHERE_STEREOGRAPH tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPHERE_STEREOGRAPH_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPHERE_STEREOGRAPH library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPHERE_STEREOGRAPH_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 checks that the two functions are inverses.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) dif
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ), allocatable :: p1(:,:)
  real ( kind = 8 ), allocatable :: p2(:,:)
  real ( kind = 8 ), allocatable :: q(:,:)
  real ( kind = 8 ), allocatable :: q1(:,:)
  real ( kind = 8 ), allocatable :: q2(:,:)
  real ( kind = 8 ) r8mat_norm_fro_affine
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SPHERE_STEREOGRAPH maps from sphere to plane.'
  write ( *, '(a)' ) '  SPHERE_STEREOGRAPH_INVERSE is the inverse map.'
  write ( *, '(a)' ) '  Check that these two functions are inverses.'

  m = 3
  n = 100

  allocate ( p1(m,n) )
  allocate ( q(m,n) )
  allocate ( p2(m,n) )
!
!  Check #1.
!
  seed = 123456789

  call uniform_on_sphere01_map ( m, n, seed, p1 )
  call sphere_stereograph ( m, n, p1, q )
  call sphere_stereograph_inverse ( m, n, q, p2 )

  dif = r8mat_norm_fro_affine ( m, n, p1, p2 )
  
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Map points from sphere to plane to sphere.'
  write ( *, '(a,i6,a,g14.6)' ) '  Frobenius difference for ', n, ' points was ', dif
!
!  Check #2.
!
  allocate ( q1(m,n) )
  allocate ( p(m,n) )
  allocate ( q2(m,n) )

  call r8mat_uniform_01 ( m, n, seed, q1 )
  q1(m,1:n) = 1.0D+00

  call sphere_stereograph_inverse ( m, n, q1, p )
  call sphere_stereograph ( m, n, p, q2 )

  dif = r8mat_norm_fro_affine ( m, n, q1, q2 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Map points from plane to sphere to plane.'
  write ( *, '(a,i6,a,g14.6)' ) '  Frobenius difference for ', n, ' points was ', dif

  deallocate ( p )
  deallocate ( p1 )
  deallocate ( p2 )
  deallocate ( q )
  deallocate ( q1 )
  deallocate ( q2 )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 checks the generalized mapping.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: center(:)
  real ( kind = 8 ) dif
  real ( kind = 8 ), allocatable :: focus(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: p1(:,:)
  real ( kind = 8 ), allocatable :: p2(:,:)
  real ( kind = 8 ), allocatable :: q1(:,:)
  real ( kind = 8 ), allocatable :: q2(:,:)
  real ( kind = 8 ) r8mat_norm_fro_affine
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  SPHERE_STEREOGRAPH standard mapping from sphere to plane.'
  write ( *, '(a)' ) '  SPHERE_STEREOGRAPH2 generalized mapping:'
  write ( *, '(a)' ) '  (focus and center are arbitrary)'
  write ( *, '(a)' ) '  Check that these two functions can agree.'

  m = 3
  n = 100

  allocate ( focus(1:m) )
  focus(1:m-1) = 0.0D+00
  focus(m) = -1.0D+00

  allocate ( center(1:m) )
  center(1:m) = 0.0D+00

  allocate ( p1(1:m,1:n) )
  allocate ( p2(1:m,1:n) )
  allocate ( q1(1:m,1:n) )
  allocate ( q2(1:m,1:n) )
!
!  Check #1.
!
  seed = 123456789
  call uniform_on_sphere01_map ( m, n, seed, p1 )

  call sphere_stereograph ( m, n, p1, q1 )

  call sphere_stereograph2 ( m, n, p1, focus, center, q2 )

  dif = r8mat_norm_fro_affine ( m, n, q1, q2 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Map points from sphere to plane.'
  write ( *, '(a,i6,a,g14.6)' ) '  Frobenius difference for ', n, ' points was ', dif
!
!  Check #2.
!
  call r8mat_uniform_01 ( m, n, seed, q1 )
  do j = 1, n
    q1(m,j) = 1.0D+00
  end do

  call sphere_stereograph_inverse ( m, n, q1, p1 )

  call sphere_stereograph2_inverse ( m, n, q1, focus, center, p2 )

  dif = r8mat_norm_fro_affine ( m, n, p1, p2 )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Map points from plane to sphere.'
  write ( *, '(a,i6,a,g14.6)' ) '  Frobenius difference for ', n, ' points was ', dif

  deallocate ( p1 )
  deallocate ( p2 )
  deallocate ( q1 )
  deallocate ( q2 )

  deallocate ( center )
  deallocate ( focus )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 checks that the two functions are inverses.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: alpha(:)
  real ( kind = 8 ), allocatable :: beta(:)
  real ( kind = 8 ), allocatable :: center(:)
  real ( kind = 8 ) dif
  real ( kind = 8 ), allocatable :: focus(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: normal(:)
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ), allocatable :: p1(:,:)
  real ( kind = 8 ), allocatable :: p2(:,:)
  real ( kind = 8 ), allocatable :: pq(:)
  real ( kind = 8 ), allocatable :: pr(:)
  real ( kind = 8 ), allocatable :: q(:,:)
  real ( kind = 8 ), allocatable :: q1(:,:)
  real ( kind = 8 ), allocatable :: q2(:,:)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8mat_norm_fro_affine
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: tang(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  SPHERE_STEREOGRAPH2 maps from sphere to plane.'
  write ( *, '(a)' ) '  SPHERE_STEREOGRAPH2_INVERSE is the inverse map.'
  write ( *, '(a)' ) '  Check that these two functions are inverses.'

  m = 3
  n = 100
  seed = 123456789

  allocate ( focus(1:m) )
  allocate ( center(1:m) )

  call r8vec_uniform_01 ( m, seed, focus )
  call r8vec_uniform_01 ( m, seed, center )
  r = r8vec_norm_affine ( m, focus, center )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Using radius = ', r
  call r8vec_transpose_print ( m, center, '  Center:' )
  call r8vec_transpose_print ( m, focus, '  Focus:' )
!
!  Check #1.
!
  allocate ( p1(1:m,1:n) )
  allocate ( p2(1:m,1:n) )
  allocate ( q(1:m,1:n) )

  call uniform_on_sphere01_map ( m, n, seed, p1 )

  do j = 1, n
    p1(1:m,j) = center(1:m) + r * p1(1:m,j)
  end do

  call sphere_stereograph2 ( m, n, p1, focus, center, q )

  call sphere_stereograph2_inverse ( m, n, q, focus, center, p2 )

  dif = r8mat_norm_fro_affine ( m, n, p1, p2 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Map points from sphere to plane to sphere.'
  write ( *, '(a,i6,a,g14.6)' ) '  Frobenius difference for ', n, ' points was ', dif

  deallocate ( p1 )
  deallocate ( p2 )
  deallocate ( q )
!
!  Check #2.
!  We have to work hard to get random points on the plane, since
!  all we know to begin with is the point of tangency and the normal.
!
  allocate ( alpha(1:n) )
  allocate ( beta(1:n) )
  allocate ( normal(1:m) )
  allocate ( p(1:m,1:n) )
  allocate ( pr(1:m) )
  allocate ( pq(1:m) )
  allocate ( q1(1:m,1:n) )
  allocate ( q2(1:m,1:n) )
  allocate ( tang(1:m) )

  tang(1:m) = 2.0 * center(1:m) - focus(1:m)

  normal(1:m) = center(1:m) - focus(1:m)

  call plane_normal_basis_3d ( tang, normal, pr, pq )

  call r8vec_uniform_01 ( n, seed, alpha )
  call r8vec_uniform_01 ( n, seed, beta )

  do j = 1, n
    q1(1:m,j) = tang(1:m) + pr(1:m) * alpha(j) + pq(1:m) * beta(j)
  end do
  call sphere_stereograph2_inverse ( m, n, q1, focus, center, p )
  call sphere_stereograph2 ( m, n, p, focus, center, q2 )

  dif = r8mat_norm_fro_affine ( m, n, q1, q2 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Map points from plane to sphere to plane.'
  write ( *, '(a,g14.6)' ) '  Frobenius difference for ', n, ' points was ', dif

  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( normal )
  deallocate ( p )
  deallocate ( pq )
  deallocate ( pr )
  deallocate ( q1 )
  deallocate ( q2 )
  deallocate ( tang )

  return
end
