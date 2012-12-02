function f ( t )

!*******************************************************************************
!
!! F evaluates the power consumption function F(T).
!
!  Modified:
!
!    03 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Rubin Landau,
!    A First Course in Scientific Computing,
!    Princeton, 2005,
!    ISBN: 0-691-12183-4
!    LC: Q183.9.L36
!
!  Parameters:
!
!    Input, real T, the argument of the function.
!
!    Output, real F, the value of the function.
!
  implicit none

  real f
  real, parameter :: pi = 3.14159265
  real t

  f = ( 4.0 + t / 365.0 + 0.5 * sin ( pi * t / 91.0 ) ) * &
      ( 2.0 + exp ( - sin ( 2.0 * pi * t ) ) )

  return
end
