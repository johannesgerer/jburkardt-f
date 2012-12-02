function sge_det ( )

!*******************************************************************************
!
!! SGE_DET computes the determinant of the current matrix.
!
!  Discussion:
!
!    In order to use this routine, the user must first have
!    specified the values of N and A, and passed them to SGE_CREATE,
!    and then called SGE_FA.
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real SGE_DET, the determinant of the matrix.
!
  use sge_module

  implicit none

  real sge_det
  integer i

  if ( .not. a_created ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_DET - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been defined yet.'
    write ( *, '(a)' ) '  Call SGE_CREATE first!'
    stop
  end if

  if ( .not. lu_computed ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_DET - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been factored yet.'
    write ( *, '(a)' ) '  Call SGE_FA first!'
    stop
  end if

  det_sge = 1.0E+00

  do i = 1, n_sge
    det_sge = det_sge * lu_sge(i,i)
  end do

  do i = 1, n_sge
    if ( pivot_sge(i) /= i ) then
      det_sge = - det_sge
    end if
  end do

  det_computed = .true.

  sge_det = det_sge

  return
end
