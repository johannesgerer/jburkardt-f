subroutine sge_delete ( )

!*******************************************************************************
!
!! SGE_DELETE deletes information about a single SGE matrix.
!
  use sge_module

  implicit none

  deallocate ( lu_sge )
  deallocate ( pivot_sge )

  n_sge = -1
  info_sge = -1

  a_created = .false.
  lu_computed = .false.
  inv_computed = .false.
  det_computed = .false.

  return
end
