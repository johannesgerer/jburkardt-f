subroutine sge_slt ( b, x )

!*******************************************************************************
!
!! SGE_SLT solves a transposed system factored by SGE_FA.
!
!  Discussion:
!
!    In order to use this routine, the user must first have
!    specified the values of N and A, and passed them to SGE_CREATE,
!    and then called SGE_FA.
!
!    SGE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real B(N), the right hand side vector.
!
!    Output, real X(N), the solution vector.
!
  use sge_module

  implicit none

  real b(n_sge)
  integer k
  integer l
  real x(n_sge)

  if ( .not. a_created ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been defined yet.'
    write ( *, '(a)' ) '  Call SGE_CREATE first!'
    stop
  end if

  if ( .not. lu_computed ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been factored yet.'
    write ( *, '(a)' ) '  Call SGE_FA first!'
    stop
  end if

  if ( info_sge /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix is singular.'
    stop
  end if

  x(1:n_sge) = b(1:n_sge)
!
!  Solve U' * Y = B.
!
  do k = 1, n_sge
    x(k) = ( x(k) - dot_product ( x(1:k-1), lu_sge(1:k-1,k) ) ) / lu_sge(k,k)
  end do
!
!  Solve ( PL )' * X = Y.
!
  do k = n_sge-1, 1, -1

    x(k) = x(k) + dot_product ( x(k+1:n_sge), lu_sge(k+1:n_sge,k) )

    l = pivot_sge(k)

    if ( l /= k ) then
      call r_swap ( x(l), x(k) )
    end if

  end do

  return
end
