subroutine setbcond ( u, v, p, temp, flag, imax, jmax, ww, we, wn, ws )

!*******************************************************************************
!
!! SETBCOND sets the boundary conditions at the boundary strip.
!
!  Discussion:
!          
!    Moreover, no-slip conditions are set at internal obstacle cells
!    by default.                                                    
!
!    For temperature, adiabatic boundary conditions are set.        
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM, 1998.
!
!  Parameters:
!
!    Input/output, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    P(0:IMAX+1,0:JMAX+1), TEMP(0:IMAX+1,0:JMAX+1), the velocity, 
!    pressure, and temperature fields.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, integer WW, WE, WN, WS, specify the boundary condition
!    to be applied on the west, east, north and south walls.
!    1 = slip
!    2 = no-slip                             
!    3 = outflow
!    4 = periodic  
!
  use nrtype

  implicit none

  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( inout ) :: p
  real ( rp ), dimension(0:,0:), intent ( inout ) :: temp
  real ( rp ), dimension(0:,0:), intent ( inout ) :: u
  real ( rp ), dimension(0:,0:), intent ( inout ) :: v
  integer, intent ( in ) :: we
  integer, intent ( in ) :: wn
  integer, intent ( in ) :: ws
  integer, intent ( in ) :: ww

  include 'defs.h'

  do j = 0, jmax+1 
!
!  Western and eastern boundary.
!
!  Free slip, U = 0, dVdN = 0.
!
    if ( ww == 1 ) then

      u(0,j) = 0.0
      v(0,j) = v(1,j)
!
!  No slip, U = 0, V = 0 at the boundary by averaging.
!
    else if ( ww == 2 ) then

      u(0,j) = 0.0
      v(0,j) = (-1.0) * v(1,j)
!
!  Outflow
!
    else if ( ww == 3 ) then

      u(0,j) = u(1,j) 
      v(0,j) = v(1,j) 
!
!  Periodic, left and right cells overlap.
!
    else if ( ww == 4 ) then

      u(0,j) = u(imax-1,j) 
      v(0,j) = v(imax-1,j)
      v(1,j) = v(imax,j)
      p(1,j) = p(imax,j)  

    end if
!
!  dTdN = 0
!
    temp(0,j) = temp(1,j)
!
!  Free slip
!
    if ( we == 1 ) then

      u(imax,j) = 0.0          
      v(imax+1,j) = v(imax,j)   
!
!  No slip
!
    else if ( we == 2 ) then

      u(imax,j) = 0.0
      v(imax+1,j) = -v(imax,j)
!
!  Outflow
!
    else if ( we == 3 ) then

      u(imax,j) = u(imax-1,j)
      v(imax+1,j) = v(imax,j)
!
!  Periodic
!
    else if ( we == 4 ) then

      u(imax,j) = u(1,j)
      v(imax+1,j) = v(2,j)
    end if

    temp(imax+1,j) = temp(imax,j)

  end do
!
!  Northern and southern boundary
!
  do i = 0, imax+1

    if ( wn == 1 ) then
      v(i,jmax) = 0.0
      u(i,jmax+1) = u(i,jmax)
    else if ( wn == 2 ) then 
      v(i,jmax) = 0.0
      u(i,jmax+1) = -u(i,jmax)
    else if ( wn == 3 ) then
      v(i,jmax) = v(i,jmax-1)
      u(i,jmax+1) = u(i,jmax)
    else if ( wn == 4 ) then
      v(i,jmax) = v(i,1)
      u(i,jmax+1) = u(i,2)
    end if

    temp(i,0) = temp(i,1) 

    if ( ws == 1 ) then
      v(i,0) = 0.0
      u(i,0) = u(i,1)
    else if ( ws == 2 ) then
      v(i,0) = 0.0
      u(i,0) = -u(i,1)
    else if ( ws == 3 ) then
      v(i,0) = v(i,1)
      u(i,0) = u(i,1)
    else if ( ws == 4 ) then
      v(i,0) = v(i,jmax-1)
      u(i,0) = u(i,jmax-1)
      u(i,1) = u(i,jmax)
      p(i,1) = p(i,jmax)
    end if

    temp(i,jmax+1) = temp(i,jmax) 

  end do
!
!  Set the boundary values at inner obstacle cells (only no-slip).
!
  do i = 1, imax
    do j = 1, jmax
!
!  Mask C_X = 000f filters the obstacle cells adjacent to fluid cells.
!
      if ( iand ( flag(i,j), c_x ) /= 0 ) then

        select case ( flag(i,j) )

          case (b_n)
	    v(i,j)   = 0.0
            u(i,j)   = -u(i,j+1)
            u(i-1,j) = -u(i-1,j+1)
            temp(i,j) = temp(i,j+1)
          case (b_e)
	    u(i,j)   = 0.0
            v(i,j)   = -v(i+1,j)
            v(i,j-1) = -v(i+1,j-1)
            temp(i,j) = temp(i+1,j)
          case (b_s)
	    v(i,j-1) = 0.0
            u(i,j)   = -u(i,j-1)
            u(i-1,j) = -u(i-1,j-1)
            temp(i,j) = temp(i,j-1)
          case (b_w)
	    u(i-1,j) = 0.0
            v(i,j)   = -v(i-1,j)
            v(i,j-1) = -v(i-1,j-1)
            temp(i,j) = temp(i-1,j)
          case (b_ne)
	    v(i,j)   = 0.0
            u(i,j)   = 0.0
            v(i,j-1) = -v(i+1,j-1)
            u(i-1,j) = -u(i-1,j+1)
            temp(i,j) = 0.5 * ( temp(i,j+1) + temp(i+1,j) )
          case (b_se)
	    v(i,j-1) = 0.0
            u(i,j)   = 0.0
            v(i,j)   = -v(i+1,j)
            u(i-1,j) = -u(i-1,j-1)
            temp(i,j) = 0.5 * ( temp(i,j-1) + temp(i+1,j) )
          case (b_sw)
	    v(i,j-1) = 0.0
            u(i-1,j) = 0.0
            v(i,j)   = -v(i-1,j)
            u(i,j)   = -u(i,j-1)
            temp(i,j) = 0.5 * ( temp(i,j-1) + temp(i-1,j) )
          case (b_nw)
	    v(i,j)   = 0.0
            u(i-1,j) = 0.0
            v(i,j-1) = -v(i-1,j-1)
            u(i,j)   = -u(i,j+1)
            temp(i,j) = 0.5 * ( temp(i,j+1) + temp(i-1,j) )
	  case default

        end select

      end if

    end do
  end do

  return
end
subroutine setspecbcond ( problem, u, v, temp, imax, jmax, ui, vi )

!*******************************************************************************
!
!! SETSPECBCOND sets problem-specific boundary conditions.
!
!  Modified:
!
!    17 February 2004
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM, 1998.
!
!  Parameters:
!
!    Input, character ( len = * ) PROBLEM, indicates the problem.
!
!    Input/output, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    TEMP(0:IMAX+1,0:JMAX+1), the velocity and temperature fields.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real UI, VI, scalars that may be used to initialize the flow.
!
  use nrtype

  implicit none

  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  character ( len = * ), intent ( in ) :: problem
  real ( rp ), dimension(0:,0:), intent(inout) :: temp
  real ( rp ), dimension(0:,0:), intent(inout) :: u
  real ( rp ), intent ( in ) :: ui
  real ( rp ), dimension(0:,0:), intent(inout) :: v
  real ( rp ), intent ( in ) :: vi
!
!  BACKSTEP
!  U = 1.0 at the left boundary                   
!
  if ( problem == "backstep" ) then
       
    do j = jmax/2+1, jmax
      u(0,j) = 1.0
    end do
!
!  CIRCLE
!  U = 1.0 at left boundary              
!    
  else if ( problem == "circle" ) then

    v(0,0) = 2.0 * vi - v(1,0)
    do j = 1, jmax
      u(0,j) = ui
      v(0,j) = 2.0 * vi - v(1,j)
    end do
!
!  CONVECTION: 
!
!  left T = 0.5 
!  right T = -0.5     
!  upper and lower wall adiabatic          
!
  else if ( problem == "convection" ) then

    do j = 0, jmax+1
      temp(0,j) = 2.0 * ( 0.5 ) - temp(1,j)
      temp(imax+1,j) = 2.0 * ( -0.5 ) - temp(imax,j)
    end do
       
    do i = 0, imax+1
      temp(i,0) = temp(i,1)
      temp(i,jmax+1) = temp(i,jmax)
    end do
!
!  DAM:
!
  else if ( problem == "dam" ) then
!
!  DCAVITY: 
!
!  U = 1.0 at the upper boundary              
!
  else if ( problem == "dcavity" ) then

    do i = 0, imax
      u(i,jmax+1) = 2.0 - u(i,jmax)   
    end do
!
!  DROP:
!
  else if ( problem == "drop" ) then
!
!  FLUIDTRAP:
! 
!  left T = 0.5 
!  right T = -0.5     
!  upper and lower wall adiabatic          
!
  else if ( problem == "fluidtrap" ) then

    do j = 0, jmax+1
      temp(0,j) = 2.0 * (0.5) - temp(1,j) 
      temp(imax+1,j) = 2.0 * (-0.5) - temp(imax,j)
    end do
       
    do i = 0, imax+1
      temp(i,0) = temp(i,1)
      temp(i,jmax+1) = temp(i,jmax)
    end do
!
!  MOLDING:
!
!  U = 1.0 in the mid of left boundary   
!
  else if ( problem == "molding" ) then

    do j = int (0.4*jmax) + 1, int (0.6*jmax)
      u(0,j) = 1.0       
    end do
!
!  PLATE:
!
!  U = 1.0 at left boundary              
!    
  else if ( problem == "plate" ) then

    v(0,0) = 2.0 * vi - v(1,0)
    do j = 1, jmax
      u(0,j) = ui
      v(0,j) = 2.0 * vi - v(1,j)
    end do
!
!  RAYLEIGH: 
!
!  top T = -0.5 
!  bottom T = 0.5  
!  left and right adiabatic     
!
  else if ( problem == "rayleigh" ) then

    do j = 0, jmax+1
      temp(0,j) = temp(1,j)         
      temp(imax+1,j) = temp(imax,j)
    end do
  
    do i = 0, imax+1
      temp(i,0) = 2.0 * (0.5) - temp(i,1)
      temp(i,jmax+1) = 2.0 * (-0.5) - temp(i,jmax)
    end do
!
!  WAVE:
!
!  U = 1.0 at the left boundary                   
!  
  else if ( problem == "wave" ) then
       
    do j = jmax/2+1, jmax
      u(0,j) = 1.0
    end do
!
!  Unrecognized problem.
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETSPECBCOND - Fatal error!'
    write ( *, '(a)' ) '  Value of problem specification not recognized.'
    write ( *, '(a)' ) '  PROBLEM = "' // trim ( problem ) // '".'
    stop

  end if

  return
end
