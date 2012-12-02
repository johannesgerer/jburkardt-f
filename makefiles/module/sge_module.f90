module sge_module

!*******************************************************************************
!
!! SGE_MODULE is a module of information about a single SGE matrix.
!
!  Discussion:
!
!    Using a module, it is possible to set up an allocatable array
!    that is shared by various routines.  In particular, a user's 
!    program, by including the statement "use sge_module", will be
!    able to examine or alter any of the variables in this module.
!
  implicit none

  real, save :: det_sge = 0.0E+00
  logical, save :: a_created = .false.
  logical, save :: lu_computed = .false.
  logical, save :: inv_computed = .false.
  logical, save :: det_computed = .false.
  integer, save :: n_sge = -2
  integer, save :: info_sge = -2
  integer, save, allocatable, dimension ( : ) :: pivot_sge
  real, save, allocatable, dimension ( :, : ) :: lu_sge

end
