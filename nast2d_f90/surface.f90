subroutine init_particles ( n, imax, jmax, delx, dely, ppc, problem, &
  u, v, partlines )

!*******************************************************************************
!
!! INIT_PARTICLES initializes particles for free boundary problems.
!
!  Modified:
!
!    08 September 2005
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input/output, integer N, the number of particles.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
!    Input, integer PPC, the number of particles per cell.
!
!    Input, character ( len = * ) PROBLEM, indicates the problem.
!
!    Input/output, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    the velocity field.
!
!    Output, particle PARTLINES(N), ?
!
  use nrtype

  implicit none

  integer, intent ( inout ) :: n

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ) :: height
  integer i
  integer, intent ( in ) :: imax
  integer ip
  integer j
  integer, intent ( in ) :: jmax
  integer jp
  real ( rp ) :: mpx
  real ( rp ) :: mpy
  type ( particleline ), dimension(n), intent ( inout ) :: partlines
  integer, intent ( in ) :: ppc
  character ( len = * ), intent ( in ) :: problem
  real ( rp ) :: rad
  real ( rp ), dimension(0:,0:), intent ( inout ) :: u
  real ( rp ), dimension(0:,0:), intent ( inout ) :: v
  real ( rp ) :: vstart
  real ( rp ) :: x
  real ( rp ) :: y

  interface
    subroutine set_part ( partline, x, y )

    use nrtype

    implicit none
 
    type ( particle ), pointer :: part
    type ( particleline ), intent ( inout ), target :: partline
    real ( rp ), intent ( in ) :: x
    real ( rp ), intent ( in ) :: y

    end subroutine set_part
  end interface

! include 'interfaces.h'
!
!  Initialization of some parameters.
!
  if ( problem == "dam" ) then

    n = 1
!
!  HEIGHT = height of the basin.
!  RAD = radius of the drop.
!  (MPX,MPY) is the midpoint of the drop.
!  VSTART is the initial velocity of the drop.
!
  else if ( problem == "drop" ) then

    n = 2
    height = 1.0 / 2.0 * jmax * dely
    rad    = 0.1 * jmax * dely
    mpx    = 0.5 * imax * delx
    mpy    = 2.0 / 3.0 * jmax * dely
    vstart = -2.0

  end if

  do i = 1, n
    partlines(i)%length = 0
    partlines(i)%particles%x = -1.0
    partlines(i)%particles%y = -1.0
  end do
!
!  Set the particles 
!
  do i = 1, imax
    do j = 1, jmax

      do ip = 1, ppc

         x = ( i - 1 ) * delx + ( ip - 0.5 ) / real ( ppc ) * delx

         do jp = 1, ppc

           y = ( j - 1 ) * dely + ( jp - 0.5 ) / real ( ppc ) * dely
    
           if ( problem == "dam" ) then

             if ( x < 0.2 * imax * delx ) then
               call set_part ( partlines(1), x, y )
             end if

           else if ( problem == "drop" ) then

             if ( y < height ) then
               call set_part ( partlines(1), x, y )
             else if ( ( x - mpx )**2 + ( y - mpy )**2 <= rad**2 ) then
               call set_part ( partlines(2), x, y )
               v(i,j) = vstart
             end if

           end if

         end do
       end do

     end do
   end do

  return
end
subroutine set_part ( partline, x, y )

!*******************************************************************************
!
!! SET_PART adds a particle to "Partline" at (x,y).                             
!
!  Modified:
!
!    08 September 2005
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input/output, particleline PARTLINE, the line of particles.
!    On output, one new particle has been added to the line.
!
!    Input, real X, Y, the coordinates of the new particle.
! 
  use nrtype

  implicit none
 
  type ( particle ), pointer :: part
  type ( particleline ), intent ( inout ), target :: partline
  real ( rp ), intent ( in ) :: x
  real ( rp ), intent ( in ) :: y
!
!  Create a new particle.
!
  allocate ( part )
!
!  Set the coordinates of the particle.
!
  part%x = x
  part%y = y
!
!  Add the particle to "Partline" in the first position after the dummy.
!
  part%next => partline%particles%next
  partline%particles%next => part
!
!  Note that PARTLINE has increased in length by one new particle.
!
  partline%length = partline%length + 1

  return
end
subroutine mark_cells ( flag, imax, jmax, delx, dely, ifull, isurf, &
  n, partlines )

!*******************************************************************************
!
!! MARK_CELLS marks the cells of the fluid domain                            
!
!  Modified:
!
!    12 March 2004
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input/output, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the 
!    type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
!    Output, integer IFULL, ISURF, the number of interior and surface cells.
!
!    Input, integer N, the number of particles.
!
!    Input/output, particle PARTLINES(N), a set of N particle lines.
!
  use nrtype

  implicit none

  integer, intent ( in ) :: n

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( inout ) :: flag
  type ( particle ), pointer :: help
  integer i
  integer, intent ( out ) :: ifull
  integer, intent ( in ) :: imax
  integer, intent ( out ) :: isurf
  integer j
  integer, intent ( in ) :: jmax
  integer k
  type ( particle ), pointer :: part
  type ( particleline ), dimension(n), intent ( inout ), target :: partlines
  type ( particle ), pointer :: temp
  real ( rp ) :: x
  real ( rp ) :: y

  include 'defs.h'
!
!  Set all cells which are not obstacle cells to empty cells.
!
  do i = 0, imax+1
    do j = 0, jmax+1
      if ( c_f <= flag(i,j) ) then
        flag(i,j) =  iand ( ior ( flag(i,j), c_e ), not ( c_nswo ) )
      end if
    end do
  end do
!
!  Mark cells containing particles as fluid cells.
!  Loop over the particles.
!
  do k = 1, n

    if ( partlines(k)%length <= 0 ) then
      cycle
    end if

    part => partlines(k)%particles

    do

      if ( .not. ( associated ( part%next ) ) ) then 
        exit
      end if

      temp => part%next
      x = temp%x
      y = temp%y
      i = int ( x / delx ) + 1;
      j = int ( y / dely ) + 1;
!
!  Delete particles that have moved into obstacle cells.  Note that 
!  the predecessor pointer doesn't advance if a deletion is involved.
!
      if ( flag(i,j) < c_f ) then
        help => temp%next
        deallocate ( part%next )
        part%next => help
        partlines(k)%length = partlines(k)%length - 1
      else
        flag(i,j) =  iand ( flag(i,j), not ( c_e ) )
      end if

      part => part%next
           
    end do

  end do
!
!  Mark the surface cells.
!
  ifull = 0
  isurf = 0

  do j = 1, jmax
    do i = 1, imax

      if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
         
        if ( iand ( flag(i-1,j), c_e ) /= 0 ) then
          flag(i,j) = ior ( flag(i,j), c_w )
        end if

        if ( iand ( flag(i+1,j), c_e ) /= 0 ) then
          flag(i,j) = ior ( flag(i,j), c_o )
        end if

        if ( iand ( flag(i,j-1), c_e ) /= 0 ) then
          flag(i,j) = ior ( flag(i,j), c_s )
        end if

        if ( iand ( flag(i,j+1), c_e ) /= 0 ) then
          flag(i,j) = ior ( flag(i,j), c_n )
        end if

        if ( flag(i,j) < c_o ) then  
          ifull = ifull + 1
        else
          isurf = isurf + 1
        end if

      end if
      
    end do

  end do
!
!  DIAGNOSTIC:  Output geometry of the fluid domain
!
!     WRITE (6,*) ' '
!     WRITE (6,*) ' Geometry of the fluid domain'
!     WRITE (6,*) ' '
!     do j = jmax+1,0,-1
!        WRITE (6,'(1x,200i5)') (flag(i,j), i=0,imax+1)
!     end do

  return
end
subroutine set_uvp_surface ( u, v, p, flag, gx, gy, imax, jmax, &
  re, delx, dely, delt )

!*******************************************************************************
!
!! SET_UVP_SURFACE sets boundary values at a free surface.                    
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
!    P(0:IMAX+1,0:JMAX+1), the velocity and pressure fields.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, real GX, GY, the X and Y components of a volume force.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real RE, the Reynolds number.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
!    Input, real DELT, the time step.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), intent ( in ) :: gx
  real ( rp ), intent ( in ) :: gy
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( inout ) :: p
  real ( rp ), intent ( in ) :: re
  real ( rp ), dimension(0:,0:), intent ( inout ) :: u
  real ( rp ), dimension(0:,0:), intent ( inout ) :: v

  include 'defs.h'
!
!  Set velocity values in empty cells to zero 
!
  do j = 1, jmax
    do i = 1, imax-1
      if ( iand ( flag(i,j),   c_e ) /= 0 .and. &
           iand ( flag(i+1,j), c_e ) /= 0 ) then
        u(i,j) = 0.0
      end if
    end do
  end do

  do j = 1, jmax-1 
    do i = 1, imax
      if ( iand ( flag(i,j),   c_e ) /= 0 .and. &
           iand ( flag(i,j+1), c_e ) /= 0 ) then
        v(i,j) = 0.0
      end if
    end do
  end do
!
!  Treat only surface cells 
!
!  mask NSWO_E=0x0f00 filters surface cells.
!
  do j = 1, jmax
    do i = 1, imax

      if ( iand ( flag(i,j), c_e ) == 0 .or. flag(i,j) < c_o ) then
 
        select case ( iand ( flag(i,j), c_nswo ) ) 

          case (C_N)   
            v(i,j) = v(i,j-1) - dely / delx * ( u(i,j) - u(i-1,j) )
            if ( iand ( flag(i-1,j+1), c_e ) /= 0 ) then
              u(i-1,j+1) = u(i-1,j) - dely / delx * ( v(i,j) - v(i-1,j) )
            end if

          case (C_S)  
            v(i,j-1) = v(i,j) + dely / delx * ( u(i,j) - u(i-1,j) )
            if ( iand ( flag(i-1,j-1), c_e ) /= 0 ) then 
              u(i-1,j-1) = u(i-1,j) + dely / delx * ( v(i,j-1) - v(i-1,j-1) )
            end if
 
	  case (C_O)
            u(i,j) = u(i-1,j) - delx / dely * ( v(i,j) - v(i,j-1) )
            if ( iand ( flag(i+1,j-1),c_e) /= 0 ) then 
              v(i+1,j-1) = v(i,j-1) - delx / dely * ( u(i,j) - u(i,j-1) )
            end if

	  case (C_W) 
            u(i-1,j) = u(i,j) + delx / dely * ( v(i,j) - v(i,j-1) )
            if ( iand ( flag(i-1,j-1), c_e ) /= 0 ) then 
              v(i-1,j-1) = v(i,j-1) + delx / dely * ( u(i-1,j) - u(i-1,j-1) )
            end if

          case (C_NO) 
            u(i,j) = u(i-1,j)
            v(i,j) = v(i,j-1)
            if ( iand ( flag(i-1,j+1),c_e) /= 0 ) then
              u(i-1,j+1) = u(i-1,j)-dely/delx*(v(i,j)-v(i-1,j))
            end if
            if ( iand ( flag(i+1,j+1),c_e) /= 0 ) then
              u(i,j+1) = u(i,j)
              v(i+1,j) = v(i,j)
            end if
            if ( iand ( flag(i+1,j-1),c_e) /= 0 ) then
              v(i+1,j-1) = v(i,j-1)-delx/dely*(u(i,j)-u(i,j-1))
            end if

          case (C_NW)
            u(i-1,j) = u(i,j)
            v(i,j)   = v(i,j-1)
            if ( iand ( flag(i-1,j+1),c_e) /= 0 ) then
              u(i-1,j+1) = u(i-1,j)
              v(i-1,j)   = v(i,j)
            end if
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              v(i-1,j-1) = v(i,j-1)+delx/dely*(u(i-1,j)-u(i-1,j-1))
            end if

          case (C_SW)
            u(i-1,j) = u(i,j)
            v(i,j-1) = v(i,j)
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              u(i-1,j-1) = u(i-1,j)
              v(i-1,j-1) = v(i,j-1)
            end if

          case (C_SO)
            u(i,j)   = u(i-1,j)
            v(i,j-1) = v(i,j)
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              u(i-1,j-1) = u(i-1,j)+dely/delx*(v(i,j-1)-v(i-1,j-1))
            end if
            if ( iand ( flag(i+1,j-1),c_e) /= 0 ) then
              u(i,j-1)   = u(i,j)
              v(i+1,j-1) = v(i,j-1)
            end if

          case (C_WO)
            u(i,j)   = u(i,j)   + delt * gx
            u(i-1,j) = u(i-1,j) + delt * gx
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              v(i-1,j-1) = v(i,j-1)+delx/dely*(u(i-1,j)-u(i-1,j-1))
            end if
            if ( iand ( flag(i+1,j-1),c_e) /= 0 ) then
              v(i+1,j-1) = v(i,j-1) - delx / dely * (u(i,j)-u(i,j-1))
            end if

          case (C_NS)
            v(i,j)   = v(i,j)   + delt * gy
            v(i,j-1) = v(i,j-1) + delt * gy
            if ( iand ( flag(i-1,j+1),c_e) /= 0 ) then
              u(i-1,j+1) = u(i-1,j)-dely/delx*(v(i,j)-v(i-1,j))
            end if
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              u(i-1,j-1) = u(i-1,j)+dely/delx*(v(i,j-1)-v(i-1,j-1))
            end if

          case (C_NWO)
            v(i,j) = v(i,j-1)-dely / delx * (u(i,j)-u(i-1,j))
            u(i,j) = u(i,j) + delt * gx
            u(i-1,j) = u(i-1,j) + delt * gx
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              v(i-1,j-1) = v(i,j-1)+delx/dely*(u(i-1,j)-u(i-1,j-1))
            end if
            if ( iand ( flag(i+1,j-1),c_e) /= 0 ) then
              v(i+1,j-1) = v(i,j-1) - delx / dely * (u(i,j)-u(i,j-1))
            end if
            if ( iand ( flag(i-1,j+1),c_e) /= 0 ) then
              v(i-1,j)   = v(i,j)
              u(i-1,j+1) = u(i-1,j)
            end if
            if ( iand ( flag(i+1,j+1),c_e) /= 0 ) then
              v(i+1,j) = v(i,j)
              u(i,j+1) = u(i,j)    
            end if

          case (C_NSW)
            u(i-1,j) = u(i,j) + delx / dely*(v(i,j)-v(i,j-1))
            v(i,j)   = v(i,j) + delt * gy
            v(i,j-1) = v(i,j-1) + delt * gy
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              v(i-1,j-1)  = v(i,j-1)
              u(i-1,j-1)  = u(i-1,j)
            end if
            if ( iand ( flag(i-1,j+1),c_e) /= 0 ) then
              v(i-1,j)   = v(i,j)
              u(i-1,j+1) = u(i-1,j)
            end if

          case (C_SWO)
            v(i,j-1) = v(i,j) + dely / delx * (u(i,j)-u(i-1,j))
            u(i,j)   = u(i,j) + delt * gx
            u(i-1,j) = u(i-1,j) + delt * gx
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              u(i-1,j-1) = u(i-1,j)
              v(i-1,j-1) = v(i,j-1)
            end if
            if ( iand ( flag(i+1,j-1),c_e) /= 0 ) then
              u(i,j-1)    = u(i,j)
              v(i+1,j-1)  = v(i,j-1)
            end if

          case (C_NSO)
            u(i,j)   = u(i-1,j)-delx / dely * (v(i,j)-v(i,j-1))
            v(i,j)   = v(i,j) + delt * gy
            v(i-1,j) = v(i-1,j) + delt * gy
            if ( iand ( flag(i-1,j+1),c_e) /= 0 ) then
              u(i-1,j+1) = u(i-1,j) - dely / delx * (v(i,j)-v(i-1,j))
            end if
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              u(i-1,j-1) = u(i-1,j) + dely / delx * (v(i,j-1)-v(i-1,j-1))
            end if
            if ( iand ( flag(i+1,j-1),c_e) /= 0 ) then
              u(i,j-1)   = u(i,j)
              v(i+1,j-1) = v(i,j-1)
            end if
            if ( iand ( flag(i+1,j+1),c_e) /= 0 ) then
              u(i,j+1)    = u(i,j)
              v(i+1,j)    = v(i,j)
            end if

          case (C_NSWO)
            u(i,j)   = u(i,j)   + delt * gx
            u(i-1,j) = u(i-1,j) + delt * gx
            v(i,j)   = v(i,j)   + delt * gy
            v(i,j-1) = v(i,j-1) + delt * gy
            if ( iand ( flag(i-1,j+1), c_e ) /= 0 ) then
              u(i-1,j+1) = u(i-1,j)
              v(i-1,j)   = v(i,j)
            end if
            if ( iand ( flag(i+1,j+1),c_e) /= 0 ) then
              u(i,j+1) = u(i,j)
              v(i+1,j) = v(i,j)
            end if
            if ( iand ( flag(i-1,j-1),c_e) /= 0 ) then
              u(i-1,j-1)  = u(i-1,j)
              v(i-1,j-1)  = v(i,j-1)
            end if
            if ( iand ( flag(i+1,j-1),c_e) /= 0 ) then
              u(i,j-1)    = u(i,j)
              v(i+1,j-1)  = v(i,j-1)
            end if

	  case default

        end select

      end if

    end do

  end do
! 
!  Second loop do pressure boundary values 
!
  do j = 1, jmax 
    do i = 1, imax
       
      if ( .not. ( iand ( flag(i,j),c_e) /= 0 .or. flag(i,j) < c_o ) ) then

        select case ( iand ( flag(i,j), C_NSWO ) )  

	  case (C_N)
            p(i,j) = 2.0 / re / dely * (v(i,j)-v(i,j-1))
	  case (C_S)
            p(i,j) = 2.0 / re / dely * (v(i,j)-v(i,j-1))  
	  case (C_O)
            p(i,j) = 2.0 / re / delx * (u(i,j)-u(i-1,j)) 
	  case (C_W)
            p(i,j) = 2.0 / re / delx * (u(i,j)-u(i-1,j)) 
          case (C_NO)
            p(i,j) = 1.0 / re / 2.0 * &
               ( (u(i,j)+u(i-1,j)-u(i,j-1)-u(i-1,j-1)) / dely &
             +   (v(i,j)+v(i,j-1)-v(i-1,j)-v(i-1,j-1)) / delx )
          case (C_NW)
            p(i,j) = -1.0 / re / 2.0 * &
               ( (u(i,j)+u(i-1,j)-u(i,j-1)-u(i-1,j-1)) / dely &
             +   (v(i+1,j)+v(i+1,j-1)-v(i,j)-v(i,j-1)) / delx )
          case (C_SW)
            p(i,j) = 1.0 / re / 2.0 * &
               ( (u(i,j+1)+u(i-1,j+1)-u(i,j)-u(i-1,j)) / dely &
             +   (v(i+1,j)+v(i+1,j-1)-v(i,j)-v(i,j-1)) / delx )
          case (C_SO)
            p(i,j) = -1.0 / re / 2.0 * &
               ( (u(i,j+1)+u(i-1,j+1)-u(i,j)-u(i-1,j)) / dely &
             +   (v(i,j)+v(i,j-1)-v(i-1,j)-v(i-1,j-1)) / delx )
          case (C_WO)
            p(i,j) = 0.0
          case (C_NS)
            p(i,j) = 0.0
          case (C_NWO)
            p(i,j) = 0.0
          case (C_NSW)
            p(i,j) = 0.0
          case (C_SWO)
            p(i,j) = 0.0
          case (C_NSO)
            p(i,j) = 0.0
          case (C_NSWO)
            p(i,j) = 0.0
          case default

        end select

      end if

    end do
  end do

  return
end
