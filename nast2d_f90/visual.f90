subroutine outputvec_bin ( u, v, p, temp, psi, zeta, heat, flag, work1, &
  xlength, ylength, imax, jmax, vecfile )

!*******************************************************************************
!
!! OUTPUTVEC_BIN writes state variables to a binary file for visualization.
!
!  Modified:
!
!    03 March 2004
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
!    P(0:IMAX+1,0:JMAX+1), TEMP(0:IMAX+1,0:JMAX+1), PSI(0:IMAX+1,0:JMAX+1),
!    ZETA(1:IMAX,1:JMAX), HEAT(0:IMAX+1,0:JMAX+1), the velocity, pressure,
!    temperature, stream function, vorticity and heat fields.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Workspace, real WORK1(0:IMAX+1,0:JMAX+1).
!
!    Input, real XLENGTH, YLENGTH, the width and height of the flow region.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, character ( len = * ) VECFILE, the name of the file to which
!    the data is written.
!
  use nrtype

  implicit none

  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( in ) :: heat
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( in ) :: p
  real ( rp ), dimension(0:,0:), intent ( in ) :: psi
  real ( rp ), dimension(0:,0:), intent ( in ) :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v
  character ( len = * ), intent ( in ) :: vecfile
  real ( rp ), dimension(1:,1:), intent ( out ) :: work1
  real ( rp ), intent ( in ) :: xlength
  real ( rp ), intent ( in ) :: ylength
  real ( rp ), dimension(1:,1:), intent ( in ) :: zeta

  include 'defs.h'

  open ( unit = 1, file = vecfile, status = 'replace', form = 'unformatted' ) 

  write ( unit = 1 ) xlength, ylength
  write ( unit = 1 ) imax, jmax

  do j = 1, jmax
    do i = 1, imax
      if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
        work1(i,j) = ( u(i,j) + u(i-1,j) ) / 2.0
      else
        work1(i,j) = 0.0
      end if
    end do
  end do

  write ( unit = 1 ) work1

  do j = 1, jmax
    do i = 1, imax
      if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
        work1(i,j) = ( v(i,j) + v(i,j-1) ) / 2.0
      else
        work1(i,j) = 0.0
      end if
    end do
  end do

  write ( unit = 1 ) work1

  do j = 1, jmax
    do i = 1, imax
      if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
        work1(i,j) = p(i,j)
      else
        work1(i,j) = 0.0
      end if
    end do
  end do

  write ( unit = 1 ) work1

  do j = 1, jmax
    do i = 1, imax
      if ( iand ( flag(i,j), c_f ) /= 0 .and. flag(i,j) < c_e ) then
        work1(i,j) = temp(i,j)
      else
        work1(i,j) = -0.5
      end if
    end do
  end do

  write ( unit = 1 ) work1

  write ( unit = 1 ) zeta
  write ( unit = 1 ) psi
  write ( unit = 1 ) heat

  close ( unit = 1 )

  return
end
subroutine comp_psi_zeta ( u, v, psi, zeta, flag, imax, jmax, delx, dely )

!*******************************************************************************
!
!! COMP_PSI_ZETA computes the stream function and vorticity.
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
!    Output, real PSI(0:IMAX+1,0:JMAX+1), the stream function.
!
!    Output, real ZETA(1:IMAX,1:JMAX), the vorticity.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of a cell.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( out ) :: psi
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v
  real ( rp ), dimension(1:,1:), intent ( out ) :: zeta

  include 'defs.h'
!
!  Computation of the vorticity ZETA at the upper right corner 
!  of cell (i,j) (only if the corner is surrounded by fluid cells)
!
  do i = 1, imax-1   
    do j = 1, jmax-1

      if (((iand(flag(i,  j),c_f) /= 0 ) .and. flag(i,j) < c_e ) .and.& 
          ((iand(flag(i+1,j),c_f) /= 0 ) .and. flag(i+1,j) < c_e ) .and.&
          ((iand(flag(i,j+1),c_f) /= 0 ) .and. flag(i,j+1) < c_e ) .and.&
          ((iand(flag(i+1,j+1),c_f) /= 0 ) .and. flag(i+1,j+1) < c_e ))&
      then

        zeta(i,j) = ( u(i,j+1) - u(i,j) ) / dely &
                   -( v(i+1,j) - v(i,j) ) / delx
 
      else

        zeta(i,j) = 0.0

      end if

    end do
  end do
!
!  Computation of the stream function at the upper right corner 
!  of cell (I,J), but only if both lower cells are fluid cells.
!
  do i = 0, imax
    psi(i,0) = 0.0
    do j = 1, jmax

      if (((iand(flag(i,j),  c_f) /= 0 ) .and. (flag(i,  j)<c_e)) .or. &
          ((iand(flag(i+1,j),c_f) /= 0 ) .and. (flag(i+1,j)<c_e))) then
        psi(i,j) = psi(i,j-1) + u(i,j) * dely
      else
        psi(i,j) = psi(i,j-1)
      end if

    end do
  end do

  return
end
subroutine comp_heat ( u, v, temp, heat, flag, re, pr, imax, jmax, delx, dely )

!*******************************************************************************
!
!! COMP_HEAT computes the heat function 
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
!    Output, real HEAT(0:IMAX+1,0:JMAX+1), the heat field.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, real RE, the Reynolds number.
!
!    Input, real PR, the Prandtl number.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of a cell.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( out ) :: heat
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  real ( rp ), intent ( in ) :: pr
  real ( rp ), intent ( in ) :: re
  real ( rp ), dimension(0:,0:), intent ( in ) :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v

  include 'defs.h'

  do i = 0, imax
    heat(i,0) = 0.0
    do j = 1, jmax
     if  (((iand(flag(i,j),  c_f) /= 0 ) .and. (flag(i,j)   < c_e)) .or. &
          ((iand(flag(i+1,j),c_f) /= 0 ) .and. (flag(i+1,j) < c_e))) then
       heat(i,j) = heat(i,j-1) &
         + dely * (u(i,j) * 0.5 * ( 1.0 + temp(i+1,j) + temp(i,j)) * re * pr &
         -                      ( temp(i+1,j) - temp(i,j) ) / delx )
      end if       
    end do
  end do

  return
end
subroutine set_particles ( n, pos1x, pos1y, pos2x, pos2y, partlines )

!*******************************************************************************
!
!! SET_PARTICLES sets initial coordinates where particles are injected. 
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, integer N, the number of particle lines.
!
!    Input, real POS1X, POS1Y, POS2X, POS2Y, the X and Y coordinates of
!    the first and last particles in the particle line.
!
!    Output, particleline PARTLINES(N), ?
!
  use nrtype

  implicit none

  integer, intent ( in ) :: n

  real ( rp ) :: hx
  real ( rp ) :: hy
  integer i
  type ( particleline ), dimension(n), intent ( out ) :: partlines
  real ( rp ), intent ( in ) :: pos1x
  real ( rp ), intent ( in ) :: pos1y
  real ( rp ), intent ( in ) :: pos2x
  real ( rp ), intent ( in ) :: pos2y
  real ( rp ) :: x
  real ( rp ) :: y

  if ( 2 <= n ) then
    hx  = ( pos2x - pos1x ) / real ( n - 1 )
    hy  = ( pos2y - pos1y ) / real ( n - 1 )
  else
    hx = 0.0
    hy = 0.0
  end if

  do i = 1, n
    x = pos1x + hx * real ( i - 1 )
    y = pos1y + hy * real ( i - 1 )
    partlines(i)%particles%x = x
    partlines(i)%particles%y = y
    partlines(i)%length = 0
    write ( *, '(a,i6,a,2g14.6)' ) '  Particle I = ', i, '  (X,Y) = ', x, y
  end do

  return
end
subroutine advance_particles ( imax, jmax, delx, dely, delt, u, v, flag, &
  n, partlines )

!*******************************************************************************
!
!! ADVANCE_PARTICLES advances particles by the Euler method.
!
!  Modified:
!
!    09 September 2005
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of a cell.
!
!    Input, real DELT, the time step.
!
!    Input/output, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    the velocity field.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer N, the number of particle lines.
!
!    ?, particleline PARTLINES(N), ?
!
  use nrtype

  implicit none

  integer, intent ( in ) :: n

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  type ( particle ), pointer :: help
  integer i
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  integer k
  type ( particle ), pointer :: part
  type ( particleline ), dimension(n), intent ( inout ), target :: partlines
  type ( particle ), pointer :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ) :: uu
  real ( rp ), dimension(0:,0:), intent ( in ) :: v
  real ( rp ) :: vv
  real ( rp ) :: x
  real ( rp ) :: x1
  real ( rp ) :: x2
  real ( rp ) :: y
  real ( rp ) :: y1
  real ( rp ) :: y2

  include 'defs.h'

  interface
    subroutine advance_at_bound ( i, j, x, y, uu, vv, u, v, flag, delx, dely, &
      delt )

    use nrtype
 
    implicit none

    real ( rp ), intent ( in ) :: delt
    real ( rp ), intent ( in ) :: delx
    real ( rp ), intent ( in ) :: dely
    integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
    integer, intent ( in ) :: i
    integer, intent ( in ) :: j
    real ( rp ), dimension(0:,0:), intent ( in ) :: u
    real ( rp ), intent ( inout ) :: uu
    real ( rp ), dimension(0:,0:), intent ( in ) :: v   
    real ( rp ), intent ( inout ) :: vv
    real ( rp ), intent ( inout ) :: x
    real ( rp ), intent ( inout ) :: y

    end subroutine advance_at_bound
  end interface

  do k = 1, n

    part => partlines(k)%particles 

    do

      if ( .not. ( associated ( part%next ) ) ) then
        exit
      end if

      temp => part%next
!
!  Advance all elements except the first element, which contains the injection
!  point for new particles.
!
      x = temp%x
      y = temp%y
!
!  Computation of new X coordinates by discretizing dX/dT=U.
!
      i = int ( x / delx ) + 1
      j = int ( ( y + 0.5 * dely ) / dely ) + 1

      x1 = ( i - 1 ) * delx
      y1 = ( ( j - 1 ) - 0.5 ) * dely
      x2 = i * delx
      y2 = ( j - 0.5 ) * dely
!
!  Bilinear interpolation.
!
      uu = ( (x2-x) * (y2-y) * u(i-1,j-1) &
        +    (x-x1) * (y2-y) * u(i,j-1)   &
        +    (x2-x) * (y-y1) * u(i-1,j)   &
        +    (x-x1) * (y-y1) * u(i,j)) / delx / dely
!
!  Computation of new Y coordinates by discretizing dY/dT=V.
!
      i = int ( ( x + 0.5 * delx ) / delx ) + 1
      j = int ( y / dely ) + 1
 
      x1 = ( real ( i - 1 ) - 0.5 ) * delx
      y1 = real ( j - 1 ) * dely
      x2 = ( real ( i ) - 0.5 ) * delx
      y2 = real ( j ) * dely
!
!  Bilinear interpolation.
!
      vv = ( (x2-x) * (y2-y) * v(i-1,j-1) &
        +    (x-x1) * (y2-y) * v(i,j-1)   &
        +    (x2-x) * (y-y1) * v(i-1,j)   &
        +    (x-x1) * (y-y1) * v(i,j)) / delx / dely
!
!  Velocity updates.
!
      x = x + delt * uu
      y = y + delt * vv
!
!  Determine new cell for the particle.
!
      i = int ( x / delx ) + 1
      j = int ( y / dely ) + 1
!
!  if the particle exits the fluid domain, delete it.
!
      if ( imax * delx <= x .or. &
           jmax * dely <= y .or. &
           x <= 0.0 .or. &
           y <= 0.0 ) then
!
!  Is TEMP the last particle?
!
        if ( .not. ( associated(temp%next) ) ) then
          deallocate(temp)
          nullify(part%next)
!
!  TEMP is NOT the last particle.
!
        else
          help => temp%next
          deallocate(temp)
          part%next => help
        end if

        partlines(k)%length = partlines(k)%length - 1

      else
!
!  Special treatment if particle would be in an inner obstacle cell.
!
        if ( flag(i,j) < c_f ) then

          call advance_at_bound ( i, j, x, y, uu, vv, u, v, flag, &
            delx, dely, delt )
        end if

        temp%x = x
        temp%y = y

      end if

      part => part%next

    end do

  end do

  return
end
subroutine advance_at_bound ( i, j, x, y, uu, vv, u, v, flag, delx, dely, delt )

!*******************************************************************************
!
!! ADVANCE_AT_BOUND computes new particle locations near a no-slip wall.
!
!  Discussion:
!
!    The computation guarantees that the new position is not in the obstacle
!    cell.  Here a modified interpolation algorithm is applied, using the
!    fact that at no-skip walls, the velocity is not only given at the
!    midpoint of the edge but on the whole edge. 
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, integer I, J, the coordinates of the cell.
!
!    Input/output, real X, Y.  On input, the particle position.
!    On output, the particle position has been updated.
!
!    Input, real UU, VV, the velocity components of the particle.
!
!    Input, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    the velocity field.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, real DELX, DELY, the width and height of a cell.
!
!    Input, real DELT, the time step.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer, intent ( in ) :: i
  integer iold
  integer, intent ( in ) :: j
  integer  jold
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ) :: ul
  real ( rp ) :: ur
  real ( rp ), intent ( inout ) :: uu
  real ( rp ), dimension(0:,0:), intent ( in ) :: v
  real ( rp ) :: vo
  real ( rp ) :: vu     
  real ( rp ), intent ( inout ) :: vv
  real ( rp ), intent ( inout ) :: x
  real ( rp ) :: x1
  real ( rp ) :: x2
  real ( rp ) :: xold
  real ( rp ), intent ( inout ) :: y
  real ( rp ) :: y1
  real ( rp ) :: y2
  real ( rp ) :: yold

  include 'defs.h'
!
!  Get the old particle position.
!
  xold = x - delt * uu
  yold = y - delt * vv
  iold = int ( xold / delx ) + 1
  jold = int ( yold / dely ) + 1
!
!  Compute new (X,Y).
!
  if ( i /= iold ) then
!
!  Compute new X.
!
    if ( flag(iold+1,jold) < c_f ) then

      ur = 0.0

    else

      if ( (jold-0.5)*dely <= yold ) then

        if ( flag(iold+1,jold+1) < c_f ) then
          y2 = jold *dely
        else
          y1 = (jold-0.5)*dely
          y2 = (jold+0.5)*dely
          ur = (u(iold,jold)  *(y2-yold) &
             +  u(iold,jold+1)*(yold-y1)) / dely
        end if

      else   
   
        if ( flag(iold+1,jold-1) < c_f ) then
          y1 = (jold-1.0)*dely
          ur = u(iold,jold)*(yold-y1)*2.0 / dely
        else
          y1 = (jold-1.5) * dely
          y2 = (jold-0.5) * dely
        end if

      end if

    end if

    if ( flag(iold-1,jold) < c_f ) then

      ul = 0.0  

    else

      if ( ( jold - 0.5 ) * dely <= yold ) then

        if ( flag(iold-1,jold+1) < c_f ) then
          y2 = jold *dely
          ul = u(iold-1,jold) * ( y2 - yold ) * 2.0 / dely
        else   
          y1 = ( jold - 0.5 ) * dely
          y2 = ( jold + 0.5 ) * dely
          ul = ( u(iold-1,jold)  * ( y2 - yold ) &
             +  u(iold-1,jold+1) * ( yold - y1 ) ) / dely
	end if 

      else   
    
        if ( flag(iold-1,jold-1) < c_f ) then
          y1 = ( jold - 1.0 ) * dely
          ul = u(iold-1,jold) * ( yold - y1 ) * 2.0 / dely
        else
          y1 = ( jold - 1.5 ) * dely
          y2 = ( jold - 0.5 ) * dely
          ul = ( u(iold-1,jold-1) * ( y2 - yold ) &
             +  u(iold-1,jold)  * ( yold - y1 ) ) / dely
        end if

      end if 

    end if

    uu = ( ul * ( iold * delx - xold ) &
         + ur * ( xold - ( iold - 1 ) * delx ) ) / delx
    x = xold + uu * delt

  end if  
!
!  New X is finished.
!
  if ( j /= jold) then
!
!  Compute new Y.
!
    if ( flag(iold,jold+1) < c_f ) then

      vo = 0.0  
 
    else

      if ( (iold-0.5) * delx <= xold ) then

        if ( flag(iold+1,jold+1) < c_f ) then
          x2 = iold*delx
          vo = v(iold,jold) * (x2-xold) * 2.0 / delx
        else  
          x1 = (iold-0.5) * delx
          x2 = (iold+0.5) * delx
          vo = (v(iold,jold)  *(x2-xold) &
             +  v(iold+1,jold)*(xold-x1)) / delx
        end if

      else    
  
        if ( flag(iold-1,jold+1) < c_f ) then
          x1 = (iold-1.0) * delx
          vo = v(iold,jold) * (xold-x1) * 2.0 / delx
        else
          x1 = (iold-1.5)*delx
          x2 = (iold-0.5)*delx
          vo = (v(iold-1,jold) * (x2-xold) &
             +  v(iold,jold)  * (xold-x1))/delx
        end if

      end if

    end if

    if ( flag(iold,jold-1) < c_f ) then

      vu = 0.0  

    else

      if ( (iold-0.5)*delx <= xold ) then

        if (flag(iold+1,jold-1) < c_f) then
          x2 = iold * delx
          vu = v(iold,jold-1) * (x2-xold) * 2.0 / delx
        else   
          x1 = (iold-0.5) * delx
          x2 = (iold+0.5) * delx
          vu = (v(iold,jold-1)   * (x2-xold) &
             +  v(iold+1,jold-1) * (xold-x1)) / delx
        end if

      else 
      
        if ( flag(iold-1,jold-1) < c_f ) then
          x1 = (iold-1.0) * delx
          vu = v(iold,jold-1) * (xold-x1)*2.0 / delx
        else 
          x1 = (iold-1.5) * delx
          x2 = (iold-0.5) * delx
          vu = (v(iold-1,jold-1) * (x2-xold) &
              + v(iold,jold-1) * (xold-x1)) / delx
        end if

      end if

    end if

    vv = ( vu * ( jold * dely - yold ) &
         + vo * ( yold - ( jold - 1 ) * dely ) ) / dely

    y = yold + vv * delt
  
  end if  
!
!  new y is finished
!
  return
end
subroutine write_particles_ascii ( tracefile, itype, t, n, partlines )

!*******************************************************************************
!
!! WRITE_PARTICLES_ASCII appends particle positions to an ASCII file.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, character ( len = * ) TRACEFILE, the file containing particle data.
!
!    Input, integer ITYPE, ?
!
!    Input, real T, the current time.
!
!    Input, integer N, the number of particle lines.
!
!    ?, particleline PARTLINES(N), ?
!
  use nrtype

  implicit none

  integer, intent ( in ) :: n

  integer count
  integer, intent ( in ) :: itype
  integer k
  integer length
  type ( particle ), pointer :: part
  type ( particleline ), dimension(n), intent ( in ), target :: partlines
  real ( rp ), intent ( in ) :: t
  character ( len = * ), intent ( in ) :: tracefile

  open ( unit = 1, file = tracefile, position = 'append', &
    status = 'old', form = 'formatted' )
!
!  Write time stamp only if a streakfile (not particle trace).
!
  if ( itype == 1 ) then
    write ( 1, * ) t, ' -1.0    -1    -1' 
  end if

  do k = 1, n

    count = 0
    part => partlines(k)%particles
    length = partlines(k)%length

    do while ( 0  < length ) 

      count = count + 1
      part => part%next
      write ( 1, '(1x,g12.6,1x,g12.6,1x,i4,1x,i4)' ) part%x, part%y, &
        count, length

      if ( .not. ( associated ( part%next ) ) ) then
        exit
      end if

    end do

  end do

  close ( unit = 1 )

  return
end
subroutine particle_tracing ( tracefile, t, imax, jmax, delx, dely, delt, &
  u, v, flag, n, partlines, iwrite )

!*******************************************************************************
!
!! PARTICLE_TRACING moves particles and append them to a file if desired.
!
!  Modified:
!
!    09 September 2005
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, character ( len = * ) TRACEFILE, the file containing particle data.
!
!    Input, real T, the current time.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of a cell.
!
!    Input, real DELT, the time step.
!
!    Input, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    the velocity field.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer N, the number of particle lines.
!
!    ?, particleline PARTLINES(N), ?
!
!    Input, integer IWRITE, ?
!
  use nrtype

  implicit none

  integer, intent ( in ) :: n

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: iwrite
  integer, intent ( in ) :: jmax
  type ( particleline ), dimension(n), intent ( inout ) :: partlines
  real ( rp ), intent ( in ) :: t
  character ( len = * ), intent ( in ) :: tracefile
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v

  interface
    subroutine advance_particles ( imax, jmax, delx, dely, delt, u, v, flag, &
      n, partlines )

    use nrtype

    implicit none

    integer, intent ( in ) :: n

    real ( rp ), intent ( in ) :: delt
    real ( rp ), intent ( in ) :: delx
    real ( rp ), intent ( in ) :: dely
    integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
    integer, intent ( in ) :: imax
    integer, intent ( in ) :: jmax
    type ( particleline ), dimension(n), intent ( inout ), target :: partlines
    real ( rp ), dimension(0:,0:), intent ( in ) :: u
    real ( rp ), dimension(0:,0:), intent ( in ) :: v

    end subroutine advance_particles
  end interface

  interface
    subroutine write_particles_ascii ( tracefile, itype, t, n, partlines )

    use nrtype

    implicit none

    integer, intent ( in ) :: n

    integer, intent ( in ) :: itype
    type ( particleline ), dimension(n), intent ( in ), target :: partlines
    real ( rp ), intent ( in ) :: t
    character ( len = * ), intent ( in ) :: tracefile

    end subroutine write_particles_ascii
  end interface

  if ( t <= 0 ) then
    open ( unit = 1, file = tracefile, status = 'replace', form = 'formatted' ) 
    close ( unit = 1 )
    call write_particles_ascii ( tracefile, 0, t, n, partlines )
  end if

  call advance_particles ( imax, jmax, delx, dely, delt, &
    U, V, flag, n, partlines )

  if ( iand ( iwrite, 1 ) /= 0 ) then
    call write_particles_ascii ( tracefile, 0, t+delt, n, partlines )
  end if
 
  return    
end
subroutine inject_particles ( n, partlines )

!*******************************************************************************
!
!! INJECT_PARTICES injects new particles for streaklines
!
!  Discussion:
!
!    Note that the most recent particle is placed at the head of list, 
!    just after the injection point.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, integer N, the number of particle lines.
!
!    ?, particleline PARTLINES(N), ?
!
  use nrtype

  implicit none

  integer, intent ( in ) :: n

  integer k
  type ( particle ), pointer :: part
  type ( particleline ), dimension(n), intent ( inout ), target :: partlines
  type ( particle ), pointer :: predecessor
  integer status_var

  do k = 1, n

    allocate ( part, stat = status_var )

    if ( status_var /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INJECT_PARTICLES: memory allocation failure.'
      STOP 'inject'
    end if

    part%x = partlines(k)%particles%x
    part%y = partlines(k)%particles%y
    predecessor => partlines(k)%particles
    part%next => predecessor%next
    predecessor%next => part
    partlines(k)%length = partlines(k)%length + 1  

  end do

  return
end
subroutine streaklines ( streakfile, iwrite, imax, jmax, delx, dely, delt, &
  t, u, v, flag, n, partlines )

!*******************************************************************************
!
!! STREAKLINES moves particles for streaklines, inject and write particle positions
!
!  Modified:
!
!    09 September 2005
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Parameters:
!
!    Input, character ( len = * ) STREAKFILE, ?
!
!    Input, integer IWRITE, ?
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of a cell.
!
!    Input, real DELT, the time step.
!
!    Input, real T, the current time.
!
!    Input, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    the velocity field.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer N, the number of particle lines.
!
!    ?, particleline PARTLINES(N), ?
!
  use nrtype

  implicit none
 
  integer, intent ( in ) :: n

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: iwrite
  integer, intent ( in ) :: jmax
  type ( particleline ), dimension(n), intent ( out ) :: partlines
  character ( len = * ), intent ( in ) :: streakfile
  real ( rp ), intent ( in ) :: t
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v

  interface
    subroutine advance_particles ( imax, jmax, delx, dely, delt, u, v, flag, &
      n, partlines )

    use nrtype

    implicit none

    integer, intent ( in ) :: n

    real ( rp ), intent ( in ) :: delt
    real ( rp ), intent ( in ) :: delx
    real ( rp ), intent ( in ) :: dely
    integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
    integer, intent ( in ) :: imax
    integer, intent ( in ) :: jmax
    type ( particleline ), dimension(n), intent ( inout ), target :: partlines
    real ( rp ), dimension(0:,0:), intent ( in ) :: u
    real ( rp ), dimension(0:,0:), intent ( in ) :: v

    end subroutine advance_particles
  end interface

  interface
    subroutine inject_particles ( n, partlines )

    use nrtype

    implicit none

    integer, intent ( in ) :: n

    type ( particleline ), dimension(n), intent ( inout ), target :: partlines

    end subroutine inject_particles
  end interface

  interface
    subroutine write_particles_ascii ( tracefile, itype, t, n, partlines )

    use nrtype

    implicit none

    integer, intent ( in ) :: n

    integer, intent ( in ) :: itype
    type ( particleline ), dimension(n), intent ( in ), target :: partlines
    real ( rp ), intent ( in ) :: t
    character ( len = * ), intent ( in ) :: tracefile

    end subroutine write_particles_ascii
  end interface

  if ( t <= 0 ) then
    open ( unit = 1, file = streakfile, status = 'replace', form = 'formatted' )
    close ( unit = 1 )
    call write_particles_ascii ( streakfile, 1, t, n, partlines )
  end if

  call advance_particles ( imax, jmax, delx, dely, delt, u, v, flag, &
    n, partlines )

  if ( iand ( iwrite, 2 ) /= 0 ) then
    call inject_particles ( n, partlines )
  end if

  if ( iand ( iwrite, 4 ) /= 0 ) then
    call write_particles_ascii ( streakfile, 1, t+delt, n, partlines )
  end if

  return
end
