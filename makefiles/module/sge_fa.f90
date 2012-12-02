subroutine sge_fa ( )

!*******************************************************************************
!
!! SGE_FA factors a general matrix.
!
!  Discussion:
!
!    In order to use this routine, the user must first have
!    specified the values of N and A, and passed them to SGE_CREATE.
!
!    SGE_FA is a simplified version of the LINPACK routine SGEFA.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
  use sge_module

  implicit none

  integer i
  integer j
  integer k
  integer l

  if ( .not. a_created ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_FA - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been defined yet.'
    write ( *, '(a)' ) '  Call SGE_CREATE first!'
    stop
  end if

  info_sge = 0

  do k = 1, n_sge-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n_sge
      if ( abs ( lu_sge(l,k) ) < abs ( lu_sge(i,k) ) ) then
        l = i
      end if
    end do

    pivot_sge(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( lu_sge(l,k) == 0.0E+00 ) then
      info_sge = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGE_FA - Fatal error!'
      write ( *, '(a,i6)' ) '  Zero pivot on step ', info_sge
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      call r_swap ( lu_sge(l,k), lu_sge(k,k) )
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    lu_sge(k+1:n_sge,k) = -lu_sge(k+1:n_sge,k) / lu_sge(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n_sge

      if ( l /= k ) then
        call r_swap ( lu_sge(l,j), lu_sge(k,j) )
      end if

      lu_sge(k+1:n_sge,j) = lu_sge(k+1:n_sge,j) &
        + lu_sge(k+1:n_sge,k) * lu_sge(k,j)

    end do

  end do

  pivot_sge(n_sge) = n_sge

  if ( lu_sge(n_sge,n_sge) == 0.0E+00 ) then
    info_sge = n_sge
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_FA - Fatal error!'
    write ( *, '(a,i6)' ) '  Zero pivot on step ', info_sge
  end if

  lu_computed = .true.

  return
end
