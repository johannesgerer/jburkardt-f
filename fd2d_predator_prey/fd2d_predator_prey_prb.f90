subroutine u_init ( alpha, beta, gamma, delta, m, n, x, y, u0 )

!*****************************************************************************80
!
!! U_INIT supplies the initial value of U at each node.
!
!  Modified:
!
!    01 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA, GAMMA, DELTA, the parameters.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) X(M,N), Y(M,N), the X and Y coordinates
!    of the nodes.
!
!    Output, real ( kind = 8 ) U0(M,N), the initial value of U at the nodes.
!
  implicit none

  integer m
  integer n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) delta
  real ( kind = 8 ) gamma
  real ( kind = 8 ) u0(m,n)
  real ( kind = 8 ) ustar
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) y(m,n)

  ustar = gamma * alpha / ( beta - gamma )

  u0(1:m,1:n) = ustar - 2.0D-07 &
    * ( x(1:m,1:n) - 0.1D+00 * y(1:m,1:n) - 225.0D+00 ) &
    * ( x(1:m,1:n) - 0.1D+00 * y(1:m,1:n) - 675.0D+00 )

  return
end
subroutine v_init ( alpha, beta, gamma, delta, m, n, x, y, v0 )

!*****************************************************************************80
!
!! V_INIT supplies the initial value of U at each node.
!
!  Modified:
!
!    30 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA, GAMMA, DELTA, the parameters.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) X(M,N), Y(M,N), the X and Y coordinates
!    of the nodes.
!
!    Output, real ( kind = 8 ) V0(M,N), the initial value of U at the nodes.
!
  implicit none

  integer m
  integer n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) delta
  real ( kind = 8 ) gamma
  real ( kind = 8 ) ustar
  real ( kind = 8 ) v0(m,n)
  real ( kind = 8 ) vstar
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) y(m,n)

  ustar = gamma * alpha / ( beta - gamma )
  vstar = ( 1.0D+00 - ustar ) * ( alpha + ustar )

  v0(1:m,1:n) = vstar - 3.0D-05 * ( x(1:m,1:n) - 450.0D+00 ) &
    - 1.2D-04 * ( y(1:m,1:n) - 150.0D+00 )

  return
end
