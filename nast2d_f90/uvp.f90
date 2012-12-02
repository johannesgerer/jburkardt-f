subroutine comp_temp ( u, v, temp, flag, work, imax, jmax, delt, &
  delx, dely, gamma, re, pr )

!*******************************************************************************
!
!! COMP_TEMP computes the temperature field.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    the velocity field.
!
!    Input/output, real TEMP(0:IMAX+1,0:JMAX+1), the temperature field.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Workspace, real WORK(0:IMAX+1,0:JMAX+1).
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELT, the time step.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
!    Input, real GAMMA, the upwind differencing factor.
!
!    Input, real RE, the Reynolds number.
!
!    Input, real PR, the Prandtl number.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ) :: dutdx
  real ( rp ) :: dvtdy
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), intent ( in ) :: gamma
  integer i
  integer, intent ( in ) :: imax
  real ( rp ) :: indelx2
  real ( rp ) :: indely2
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ) :: laplt
  real ( rp ), intent ( in ) :: pr
  real ( rp ), intent ( in ) :: re
  real ( rp ), dimension(0:,0:), intent ( inout ) :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v
  real ( rp ), dimension(0:,0:), intent ( out ) :: work

  include 'defs.h'
      
  indelx2 = 1.0 / ( delx * delx )
  indely2 = 1.0 / ( dely * dely )

  do i = 1, imax
    do j = 1, jmax

      if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then

        laplt = ( temp(i+1,j) - 2.0 * temp(i,j) + temp(i-1,j) ) * indelx2 &
              + ( temp(i,j+1) - 2.0 * temp(i,j) + temp(i,j-1) ) * indely2

        dutdx = ( ( u(i,j)   * 0.5 * ( temp(i,j)   + temp(i+1,j) )   &
                  - u(i-1,j) * 0.5 * ( temp(i-1,j) + temp(i,j)   ) ) &
          + gamma * ( abs ( u(i,j)   ) * 0.5 * ( temp(i,j)   - temp(i+1,j) ) &
                    - abs ( u(i-1,j) ) * 0.5 * ( temp(i-1,j) - temp(i,j) ) ) &
                 ) / delx

        dvtdy = ( ( v(i,j)   * 0.5 * ( temp(i,j)   + temp(i,j+1) ) &
              -     v(i,j-1) * 0.5 * ( temp(i,j-1) + temp(i,j)   ))&
              + gamma * ( abs ( v(i,j)   )   * 0.5 * (temp(i,j)  -temp(i,j+1)) &
              -           abs ( v(i,j-1) ) * 0.5 * (temp(i,j-1)-temp(i,j)))&
              ) / dely

        work(i,j) = temp(i,j) + delt * ( laplt / re / pr - dutdx - dvtdy )

      end if

    end do
  end do

  do i = 1, imax
    do j = 1, jmax

      if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
        temp(i,j) = work(i,j)
      end if

    end do
  end do

  return
end
subroutine comp_fg ( u, v, temp, f, g, flag, imax, jmax, delt, delx, dely, &
  gx, gy, gamma, re, beta )

!*******************************************************************************
!
!! COMP_FG computes the F and G fields.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    TEMP(0:IMAX+1,0:JMAX+1), the velocity and temperature fields.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELT, the time step.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
!    Input, real GX, GY, the X and Y components of a volume force.
!
!    Input, real GAMMA, the upwind differencing factor.
!
!    Input, real RE, the Reynolds number.
!
!    Input, real BETA, the coefficient of volume expansion.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: beta
  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ) :: du2dx
  real ( rp ) :: duvdx
  real ( rp ) :: duvdy
  real ( rp ) :: dv2dy
  real ( rp ), dimension(0:,0:), intent ( out ) :: f
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( out ) :: g
  real ( rp ), intent ( in ) :: gamma
  real ( rp ), intent ( in ) :: gx
  real ( rp ), intent ( in ) :: gy
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ) :: laplu
  real ( rp ) :: laplv
  real ( rp ), intent ( in ) :: re
  real ( rp ), dimension(0:,0:), intent ( in ) :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v

  include 'defs.h'
!
!  Compute flux field F
! 
!  Only if both adjacent cells are fluid cells.
!
  do i = 1, imax-1
    do j = 1, jmax

      if ((( iand ( flag(i,j),  c_f) /= 0) .and. flag(i,j)   < c_e ) .and.&
          (( iand ( flag(i+1,j),c_f) /= 0) .and. flag(i+1,j) < c_e ) ) then

        du2dx = ((u(i,j)+u(i+1,j))*(u(i,j)+u(i+1,j)) &
                     + gamma*abs ( u(i,j)+u(i+1,j))*(u(i,j)-u(i+1,j)) &
                     - (u(i-1,j)+u(i,j))*(u(i-1,j)+u(i,j)) &
                     - gamma*abs ( u(i-1,j)+u(i,j))*(u(i-1,j)-u(i,j))) &
                     / (4.0*delx)

        duvdy = ((v(i,j)+v(i+1,j))*(u(i,j)+u(i,j+1)) &
                     + gamma*abs ( v(i,j)+v(i+1,j))*(u(i,j)-u(i,j+1)) &
                     - (v(i,j-1)+v(i+1,j-1))*(u(i,j-1)+u(i,j)) &
                     - gamma*abs ( v(i,j-1)+v(i+1,j-1))*(u(i,j-1)-u(i,j))) &
                     / ( 4.0 * dely )

        laplu = ( u(i+1,j) - 2.0 * u(i,j) + u(i-1,j) ) / delx / delx &
                     + ( u(i,j+1) - 2.0 * u(i,j) + u(i,j-1) ) / dely / dely
   
        f(i,j) = u(i,j) + delt * ( laplu / re - du2dx - duvdy + gx ) &
          - delt * beta * gx * ( temp(i,j) + temp(i+1,j) ) / 2.0

      else

        f(i,j) = u(i,j)

      end if

    end do
  end do
!
!  Compute flux field G
! 
!  Only if both adjacent cells are fluid cells
!
  do i = 1, imax
    do j = 1, jmax-1

            if ((( iand ( flag(i,j),  c_f)/=0).and.(flag(i,j)<c_e)).and.&
                (( iand ( flag(i,j+1),c_f)/=0).and.(flag(i,j+1)<c_e))) then

               duvdx = ((u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) &
                     + gamma*abs ( u(i,j)+u(i,j+1))*(v(i,j)-v(i+1,j)) &
                     - (u(i-1,j)+u(i-1,j+1))*(v(i-1,j)+v(i,j)) &
                     - gamma*abs ( u(i-1,j)+u(i-1,j+1))*(v(i-1,j)-v(i,j))) &
	             / ( 4.0 * delx )

               dv2dy = ((v(i,j)+v(i,j+1))*(v(i,j)+v(i,j+1)) &
                     + gamma*abs ( v(i,j)+v(i,j+1))*(v(i,j)-v(i,j+1)) &
                     - (v(i,j-1)+v(i,j))*(v(i,j-1)+v(i,j)) &
                     - gamma*abs ( v(i,j-1)+V(i,j))*(v(i,j-1)-v(i,j))) &
	             / ( 4.0 * dely )

               laplv = ( v(i+1,j) - 2.0 * v(i,j) + v(i-1,j) ) / delx / delx &
                     + ( v(i,j+1) - 2.0 * v(i,j) + v(i,j-1) ) / dely / dely

               g(i,j) = v(i,j) + delt * ( laplv / re - duvdx - dv2dy + gy ) &
                      - delt * beta * gy * ( temp(i,j) + temp(i,j+1) ) / 2.0    

            else
 
               g(i,j) = v(i,j)

            end if

         end do
      end do
!
!  F and G at external boundary
!
  do j = 1, jmax
    f(0,j)    = u(0,j)
    f(imax,j) = u(imax,j)
  end do
 
  do i = 1, imax
    g(i,0)    = v(i,0)
    g(i,jmax) = v(i,jmax)
  end do

  return
end
subroutine comp_rhs ( f, g, rhs, flag, imax, jmax, delt, delx, dely )

!*******************************************************************************
!
!! COMP_RHS computes the righthand side field.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, real F(0:IMAX+1,0:JMAX+1), G(0:IMAX+1,0:JMAX+1), ?
!
!    Output, real RHS(0:IMAX+1,0:JMAX+1), ?
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELT, the time step.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ), dimension(0:,0:), intent ( in ) :: f
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( in ) :: g
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( out ) :: rhs

  include 'defs.h'

  do i = 1, imax
    do j = 1, jmax
!
!  Only for fluid and non-surface cells.
!
      if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_o ) then  
        rhs(i,j) = ( ( f(i,j) - f(i-1,j) ) / delx &
                 +   ( g(i,j) - g(i,j-1) ) / dely ) / delt
      end if

    end do
  end do

  return
end
subroutine poisson ( p, rhs, flag, imax, jmax, delx, dely, &
  eps, iter, itermax, omega, res, ifull, p_bound )

!*******************************************************************************
!
!! POISSON carries out an iterative solution of Poisson's equation.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, real P(0:IMAX+1,0:JMAX+1), the pressure field.
!
!    Input, real RHS(0:IMAX+1,0:JMAX+1), ?
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
!    Input, real EPS, the stopping tolerance for the pressure iteration.
!
!    Output, integer ITER, the number of iterations taken.
!
!    Input, integer ITERMAX, the maximum number of pressure iterations
!    in one time step.
!
!    Input, real OMEGA, the SOR relaxation parameter.
!
!    Output, real RES, the residual.
!
!    Input, integer IFULL, ?
!
!    Input, integer P_BOUND, flag for the treatment of the pressure boundary.
!
  use nrtype
     
  implicit none

  real ( rp ) :: add
  real ( rp ) :: beta_2
  real ( rp ) :: beta_mod
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ), intent ( in ) :: eps
  integer eps_e
  integer eps_n
  integer eps_s
  integer eps_w
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer i
  integer, intent ( in ) :: ifull
  integer, intent ( in ) :: imax
  integer, intent ( inout ) :: iter
  integer, intent ( in ) :: itermax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), intent ( in ) :: omega
  real ( rp ), dimension(0:,0:), intent ( inout ) :: p
  integer, intent ( in ) :: p_bound
  real ( rp ) :: p0
  real ( rp ) :: rdx2
  real ( rp ) :: rdy2
  real ( rp ), intent ( out ) :: res
  real ( rp ), dimension(0:,0:), intent ( in ) :: rhs

  include 'defs.h'

  rdx2 = 1.0 / delx / delx
  rdy2 = 1.0 / dely / dely
  beta_2 = - omega / ( 2.0 * ( rdx2 + rdy2 ) )

  p0 = 0.0

  do i = 1, imax
    do j = 1, jmax

      if ( iand ( flag(i,j), c_f ) /= 0 ) then
        p0 = p0 + p(i,j) * p(i,j)
      end if

    end do
  end do

  p0 = sqrt ( p0 / ifull )

  if ( p0 < 0.0001 ) then
    p0 = 1.0
  end if
!
!  SOR iteration
!
  do iter = 1, itermax
!
!  P_BOUND is the boundary-condition type for pressure
!
!  Modify the equation at the boundary.
!
    if ( p_bound == 1 ) then
!
!  Relaxation for fluid cells.
!
      do i = 1, imax
        do j = 1, jmax 
!
!  Five point star for interior fluid cells.
!
          if ( flag(i,j) == c_a ) then 

            p(i,j) = ( 1.0 - omega ) * p(i,j) & 
              - beta_2 * ( ( p(i+1,j) + p(i-1,j) ) * rdx2 &
              +            ( p(i,j+1) + p(i,j-1) ) * rdy2 - rhs(i,j) )
!
!  Modified star near boundary.
!
          else if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_o ) then

            if ( flag(i+1,j) < c_f ) then
              eps_e = 0
            else
              eps_e = 1
            end if

            if ( flag(i-1,j) < c_f ) then
              eps_w = 0
            else
              eps_w = 1
            end if

            if ( flag(i,j+1) < c_f ) then
              eps_n = 0
            else
              eps_n = 1
            end if

            if ( flag(i,j-1) < c_f ) then
              eps_s = 0
            else
              eps_s = 1
            end if

            beta_mod = -omega / ( + ( eps_e + eps_w ) * rdx2 &
              + ( eps_n + eps_s ) * rdy2 )

            p(i,j) = ( 1.0 - omega ) * p(i,j) - beta_mod * ( &
                ( eps_e * p(i+1,j) + eps_w * p(i-1,j) ) * rdx2 &
              + ( eps_n * p(i,j+1) + eps_s * p(i,j-1) ) * rdy2 - rhs(i,j) )

          end if

        end do
      end do
!
!  Compute the residual.
!
      res = 0.0

      do i = 1, imax
        do j = 1, jmax
!
!  Only fluid cells
!
          if ( flag(i+1,j) < c_f ) then
            eps_e = 0
          else
            eps_e = 1
          end if

          if ( flag(i-1,j) < c_f ) then
            eps_w = 0
          else
            eps_w = 1
          end if

          if ( flag(i,j+1) < c_f ) then
            eps_n = 0
          else
            eps_n = 1
          end if

          if ( flag(i,j-1) < c_f ) then
            eps_s = 0
          else
            eps_s = 1
          end if

          if ( ( iand ( flag(i,j), c_f ) /= 0 ) .and. (flag(i,j) < c_o ) ) then
            add = ( eps_e * ( p(i+1,j) - p(i,j) ) &
                -   eps_w * ( p(i,j)   - p(i-1,j) ) ) * rdx2 &
                + ( eps_n * ( p(i,j+1) - p(i,j) ) &
                -   eps_s * ( p(i,j)   - p(i,j-1) ) ) * rdy2 -  rhs(i,j)
            res = res + add * add
 	  end if 

        end do
      end do

      res = sqrt ( res / ifull ) / p0
 
      if ( res < eps ) then
        return
      end if

    else if ( p_bound == 2 ) then
!
!  Copy values at external boundary...  
!
      do i = 1, imax
        p(i,0)      = p(i,1)
	p(i,jmax+1) = p(i,jmax)
      end do

      do j = 1, jmax
	p(0,j)      = p(1,j)
	p(imax+1,j) = p(imax,j)
      end do
!
!  and at interior boundary cells.
!
      do i = 1, imax
        do j = 1, jmax

	  if ( b_n <= flag(i,j) .and. flag(i,j) <= b_se ) then

	        if ( flag(i,j) == b_n ) then
                  p(i,j) = p(i,j+1)
                else if ( flag(i,j) == b_e ) then
                  p(i,j) = p(i+1,j)
                else if ( flag(i,j) == b_s ) then
                  p(i,j) = p(i,j-1)
                else if ( flag(i,j) == b_w ) then
                  p(i,j) = p(i-1,j)
                else if ( flag(i,j) == b_ne ) then
                  p(i,j) = 0.5 * ( p(i,j+1) + p(i+1,j) )
                else if ( flag(i,j) == b_se ) then
                  p(i,j) = 0.5 * ( p(i,j-1) + p(i+1,j) )
                else if ( flag(i,j) == b_sw ) then
                  p(i,j) = 0.5 * ( p(i,j-1) + p(i-1,j) )
                else if ( flag(i,j) == b_nw ) then
                  p(i,j) = 0.5 * ( p(i,j+1) + p(i-1,j) )
                end if
              end if

            end do
          end do
!
!  Relaxation method for fluid cells.
!
          do i = 1,imax 
            do j = 1,jmax

	      if ( iand ( flag(i,j),c_f) /= 0 .and. flag(i,j) < c_o ) then
	        p(i,j) = ( 1.0 - omega ) * p(i,j) &
                       - beta_2 * ( (p(i+1,j) + p(i-1,j)) * rdx2 &
                       +            (p(i,j+1) + p(i,j-1)) * rdy2 &
                       - rhs(i,j) )
              end if
       
            end do
          end do
!
!  Computation of the residual.
!
          res = 0.0

          do i = 1,imax
            do j = 1,jmax
!
!  Only fluid cells
!
	  if ( ( iand ( flag(i,j),c_f) /= 0 ) .and. flag(i,j) < c_o ) then
            add = ( p(i+1,j) - 2.0 * p(i,j) + p(i-1,j)) * rdx2 &
              + ( p(i,j+1) - 2.0 * p(i,j) + p(i,j-1)) * rdy2 - rhs(i,j)
            res = res + add * add
          end if

        end do
      end do
	
      res = sqrt ( res / ifull ) / p0

      if ( res < eps ) then
        return
      end if
  
    end if
      
  end do

  return
end
subroutine adap_uv ( u, v, f, g, p, flag, imax, jmax, delt, delx, dely )

!*******************************************************************************
!
!! ADAP_UV computes the U and V fields if adjacent cells are fluid cells.
!
!  Discussion:
!
!    We only update U or V when both adjacent cells are fluid cells.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input/output, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    the velocity field.
!
!    Input, real F(0:IMAX+1,0:JMAX+1), G(0:IMAX+1,0:JMAX+1), ?
!
!    Input, real P(0:IMAX+1,0:JMAX+1), the pressure field.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELT, the time step.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ), dimension(0:,0:), intent ( in ) :: f
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( in ) :: g
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( in ) :: p
  real ( rp ), dimension(0:,0:), intent ( out ) :: u
  real ( rp ), dimension(0:,0:), intent ( out ) :: v

  include 'defs.h'

  do i = 1, imax - 1
    do j = 1, jmax

      if ((( iand ( flag(i,j),   c_f ) > 0 ) .and. (flag(i,  j) < c_e ) ) .and. &
          (( iand ( flag(i+1,j), c_f ) > 0 ) .and. (flag(i+1,j) < c_e ) )) then

        u(i,j) = f(i,j) - ( p(i+1,j) - p(i,j) ) * delt / delx

      end if

    end do
  end do

  do i = 1, imax
    do j = 1, jmax - 1

      if ((( iand ( flag(i,j),  c_f ) > 0 ) .and. (flag(i,  j) < c_e ) ) .and.&
          (( iand ( flag(i,j+1),c_f ) > 0 ) .and. (flag(i,j+1) < c_e ) )) then

        v(i,j) = g(i,j) - ( p(i,j+1) - p(i,j) ) * delt / dely

      end if

    end do
  end do

  return
end
subroutine comp_delt ( delt, t, imax, jmax, delx, dely, u, v, re, pr, &
  tau, iwrite, del_trace, del_inj, del_streak, del_vec )

!*******************************************************************************
!
!! COMP_DELT computes the adaptive time stepsize.
!
!  Discussion:
!
!    The adaptive time stepsize must satisfy the CFL stability criteria.
!
!    The routine also sets the flag "write" if some data
!    has to be written into a file. 
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Output, real DELT, the time step.
!
!    Input, real T, the current time.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
!    Input, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    the velocity field.
!
!    Input, real RE, the Reynolds number.
!
!    Input, real PR, the Prandtl number.
!
!    Input, real TAU, the safety factor for timestep control.
!
!    Input, integer IWRITE, ?
!
!    Input, real DEL_TRACE, ?
!
!    Input, real DEL_INJ, ?
!
!    Input, real DEL_STREAK, ?
!
!    Input, real DEL_VEC, ?
!
  use nrtype

  implicit none
 
  real ( rp ), intent ( in ) :: del_inj
  real ( rp ), intent ( in ) :: del_streak
  real ( rp ), intent ( in ) :: del_trace
  real ( rp ), intent ( in ) :: del_vec
  real ( rp ) :: deltrepr
  real ( rp ), intent ( out ) :: delt
  real ( rp ) :: deltu
  real ( rp ) :: deltv
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer i
  integer, intent ( in ) :: imax
  integer, intent ( out ) :: iwrite
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), intent ( in ) :: pr
  real ( rp ), intent ( in ) :: re
  real ( rp ), intent ( in ) :: t
  real ( rp ) :: t_new
  real ( rp ), intent ( in ) :: tau
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ) :: umax
  real ( rp ), dimension(0:,0:), intent ( in ) :: v
  real ( rp ) :: vmax
!
!  Satisfy the CFL conditions.
!
!  If a very small TAU, then no time stepsize control.
!
  if ( 1.0e-10 <= tau ) then

    umax = 1.0e-10
    vmax = 1.0e-10 

    do i = 0, imax+1
      do j = 1, jmax+1
        if ( umax < abs ( u(i,j) ) ) then
          umax = abs ( u(i,j) )
        end if
      end do
    end do

    do i = 1, imax+1
      do j = 0, jmax+1
        if ( vmax < abs ( v(i,j) ) ) then
          vmax = abs ( v(i,j) )
        end if
      end do
    end do

    deltu = delx / umax
    deltv = dely / vmax 

    if ( pr < 1 ) then
      deltrepr = 1.0 / ( 1.0 / ( delx * delx ) &
                       + 1.0 / ( dely * dely ) ) * re * pr / 2.0
    else
      deltrepr = 1.0 / ( 1.0 / ( delx * delx ) &
                       + 1.0 / ( dely * dely ) ) * re / 2.0
    end if

    if ( deltu < deltv ) then
      if ( deltu < deltrepr ) then
        delt = deltu
      else
        delt = deltrepr
      end if
    else
      if ( deltv < deltrepr ) then
        delt = deltv
      else
        delt = deltrepr
      end if
    end if
!
!  Multiply by safety factor.
!
    delt = tau * delt   

  end if
!
!  Look if some data has to be written to a file in the next time step.
!
  iwrite = 0
  t_new = t + delt

  if ( int ( t / del_trace ) /= int ( t_new / del_trace ) ) then
    iwrite = iwrite + 1
  end if

  if ( int ( t / del_inj ) /= int ( t_new / del_inj ) ) then
    iwrite = iwrite + 2
  end if

  if ( int ( t / del_streak ) /= int ( t_new / del_streak ) ) then
    iwrite = iwrite + 4
  end if

  if ( int ( t / del_vec ) /= int ( t_new / del_vec ) ) then
    iwrite = iwrite + 8
  end if

  return
end
