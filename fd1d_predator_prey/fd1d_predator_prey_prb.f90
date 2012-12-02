subroutine u_init ( n, x, u0 )

!*****************************************************************************80
!
!! U_INIT supplies the initial value of U at each node.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = 8 ) X(N), the X coordinates of the nodes.
!
!    Output, real ( kind = 8 ) U0(N), the initial value of U at the nodes.
!
  implicit none

  integer n

  real ( kind = 8 ) u0(n)
  real ( kind = 8 ) x(n)

  u0(1:n) = exp ( - ( x(1:n) - real ( 100.0, kind = 8 ) )**2 ) &
    / real ( 5.0, kind = 8 )

  return
end
subroutine v_init ( n, x, v0 )

!*****************************************************************************80
!
!! V_INIT supplies the initial value of V at each node.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = 8 ) X(N), the X coordinates of the nodes.
!
!    Output, real ( kind = 8 ) U0(N), the initial value of U at the nodes.
!
  implicit none

  integer n

  real ( kind = 8 ) v0(n)
  real ( kind = 8 ) x(n)

  v0(1:n) = real ( 2.0, kind = 8 ) / real ( 5.0, kind = 8 )

  return
end
