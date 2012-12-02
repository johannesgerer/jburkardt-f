subroutine r_swap ( x, y )

!*******************************************************************************
!
!! R_SWAP swaps two real values.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real x
  real y
  real z

  z = x
  x = y
  y = z

  return
end
