subroutine midpoint ( a, b, f, int_num, quad )

!*******************************************************************************
!
!! MIDPOINT approximates an integral using the composite midpoint rule.
!
!  Modified:
!
!    03 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the endpoints of the interval of integration.
!
!    Input, real external F, the name of the function to be integrated.
!
!    Input, integer INT_NUM, the number of intervals to be used.
!
!    Output, real QUAD, the approximate value of the integral.
!
  implicit none

  real a
  real b
  real, external :: f
  integer i
  integer int_num
  real int_width
  real quad
  real t

  quad = 0.0

  do i = 1, int_num
    
    t = ( real ( 2 * int_num - 2 * i + 1 ) * a   &
        + real (               2 * i - 1 ) * b ) &
        / real ( 2 * int_num )

    quad = quad + f ( t )

  end do

  int_width = ( b - a ) / real ( int_num )

  quad = quad * int_width

  return
end
