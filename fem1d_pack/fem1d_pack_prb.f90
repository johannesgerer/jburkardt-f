program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM1D_PACK_PRB.
!
!  Discussion:
!
!    FEM1D_PACK_PRB calls the various FEM1D_PACK tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_PACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FEM1D_PACK library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_PACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 verifies LOCAL_BASIS_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) node_x(node_num)
  real ( kind = 8 ) phi(node_num)
  real ( kind = 8 ) phi_matrix(node_num,node_num)
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  LOCAL_BASIS_1D evaluates the local basis functions'
  write ( *, '(a)' ) '  for a 1D element.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test that the basis functions, evaluated at the nodes,'
  write ( *, '(a)' ) '  form the identity matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  node_x(1:node_num) = (/ 1.0D+00, 2.0D+00, 4.0D+00, 4.5D+00 /);

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node coordinates:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, node_x(j)
  end do

  do j = 1, node_num
    x = node_x(j)
    call local_basis_1d ( node_num, node_x, x, phi_matrix(1,j) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(I,J) = PHI(I) at node (J):'
  write ( *, '(a)' ) ' '
  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi_matrix(i,1:node_num)
  end do

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The PHI functions should sum to 1 at random X values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X        Sum ( PHI(:)(X) )'
  write ( *, '(a)' ) ' '

  do j = 1, 5
    x = r8_uniform ( 1.0D+00, 4.5D+00, seed )
    call local_basis_1d ( node_num, node_x, x, phi )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, sum ( phi(1:node_num) )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 verifies LOCAL_BASIS_PRIME_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 4

  real ( kind = 8 ) dphidx(node_num)
  real ( kind = 8 ) dphidx_matrix(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) node_x(node_num)
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  LOCAL_BASIS_PRIME_1D evaluates the local basis function'
  write ( *, '(a)' ) '  derivatives for a 1D element.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  node_x(1:node_num) = (/ 1.0D+00, 2.0D+00, 4.0D+00, 4.5D+00 /);

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node coordinates:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, node_x(j)
  end do

  do j = 1, node_num
    x = node_x(j)
    call local_basis_prime_1d ( node_num, node_x, x, dphidx_matrix(1,j) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(I,J) = dPHIdx(I) at node(J):'
  write ( *, '(a)' ) '  The diagonal should be 0.'
  write ( *, '(a)' ) '  Columns should sum to 0.'
  write ( *, '(a)' ) ' '
  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) dphidx_matrix(i,1:node_num)
  end do

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dPHIdx functions should sum to 0 at random X values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X        Sum ( dPHIdx(:)(X) )'
  write ( *, '(a)' ) ' '

  do j = 1, 5
    x = r8_uniform ( 1.0D+00, 4.5D+00, seed )
    call local_basis_prime_1d ( node_num, node_x, x, dphidx )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, sum ( dphidx(1:node_num) )
  end do

  return
end
