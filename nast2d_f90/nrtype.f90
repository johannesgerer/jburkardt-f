module nrtype

  integer, parameter :: i1b = selected_int_kind(2)
  integer, parameter :: i2b = selected_int_kind(4)
  integer, parameter :: i4b = selected_int_kind(9)

  integer, parameter :: sp = kind ( 1.0e+00 )
  integer, parameter :: dp = kind ( 1.0d+00 )
!
!  Choose 
!    single precision by "RP = SP",
!    double precision by "RP = DP" 
!    
  integer, parameter :: rp = dp

  type :: particle
    real ( rp ) :: x 
    real ( rp ) :: y 
    type ( particle ), pointer :: next 
  end type

  type :: particleline
    integer :: length 
    type ( particle ) :: particles 
  end type

end module nrtype
