subroutine sge_create ( n, a )

!*******************************************************************************
!
!! SGE_CREATE stores information about a single SGE matrix.
!
!  Discussion:
!
!    The user calls this routine with input values of N and A.
!    A copy of N is stored.  A copy of A is stored in the LU factor
!    array, but the LU factors are not computed by this routine.
!    Any old information in the module is cancelled.
!
!  Modified:
!
!    03 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the array.
!
!    Input, real A(N,N), the array.
!
  use sge_module

  implicit none

  integer n

  real a(n,n)
!
!  Make space for the LU factors and pivot vector.
!
  allocate ( lu_sge(1:n,1:n) )
  allocate ( pivot_sge(1:n) )
!
!  "Remember" the value of N.
!
  n_sge = n
!
!  Note that we have not computed the LU factor yet.
!
  info_sge = -1
!
!  Initialize the LU factor.
!
  lu_sge(1:n,1:n) = a(1:n,1:n)
!
!  Note what we have and have not done.
!
  a_created = .true.
  lu_computed = .false.
  inv_computed = .false.
  det_computed = .false.

  return
end
