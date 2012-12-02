subroutine read_parameters ( problem, inputfile, infile, outfile, xlength, &
  ylength, imax, jmax, t_end, delt, tau, delx, dely, del_vec, del_trace, &
  del_streak, del_inj, vecfile, tracefile, streakfile, n, pos1x, pos1y, &
  pos2x, pos2y, itermax, eps, omega, gamma, p_bound, re, pr, beta, gx, &
  gy, ui, vi, ti, ww, we, wn, ws )

!*******************************************************************************
!
!! READ_PARAMETERS reads parameter data and returns to main program
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM, 1998.
!
!  Parameters:
!
!    Output, character ( len = * ) PROBLEM, indicates the problem.
!
!    Output, character ( len = * ) INPUTFILE, ?
! 
!    Output, character ( len = * ) INFILE, ?
!
!    Output, character ( len = * ) OUTFILE, ?
!
!    Output, real XLENGTH, YLENGTH, the width and height of the region.
!
!    Output, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Output, real T_END, the final time.
!
!    Output, real DELT, the time step.
!
!    Output, real TAU, the safety factor for time step control.
!
!    Output, real DELX, DELY, the width and height of one cell.
!
!    Output, real DEL_VEC, ?
!
!    Output, real DEL_TRACE, ?
!
!    Output, real DEL_STREAK, ?
!
!    Output, real DEL_INJ, ?
!    
!    Output, character ( len = * ) VECFILE, ?
!
!    Output, character ( len = * ) TRACEFILE, ?
!
!    Output, character ( len = * ) STREAKFILE, ?
!
!    Output, integer N, ?
!
!    Output, real POS1X, POS1Y, POS2X, POS2Y, ?
!
!    Output, integer ITERMAX, the maximum number of pressure iterations.
!
!    Output, real EPS, the stopping tolerance for the pressure iteration.
!
!    Output, real OMEGA, the relaxation parameter for the SOR iteration.
!
!    Output, real GAMMA, the upwind differencing factor.
!
!    Output, real P_BOUND, the boundary condition type for pressure.
!
!    Output, real RE, the Reynolds number.
!
!    Output, real PR, the Prandtl number.
!
!    Output, real BETA, coefficient of thermal expansion.
!
!    Output, real GX, GY, the X and Y components of a volume force.
!
!    Output, real UI, VI, TI, scalars that may be used to initialize 
!    the flow and temperature.
!
!    Output, integer WW, WE, WN, WS, specify the boundary condition
!    to be applied on the west, east, north and south walls.
!    1 = slip
!    2 = no-slip                             
!    3 = outflow
!    4 = periodic  
!
  use nrtype

  implicit none

  real ( rp ), intent ( out ) :: beta
  real ( rp ), intent ( out ) :: del_inj
  real ( rp ), intent ( out ) :: del_streak
  real ( rp ), intent ( out ) :: del_trace
  real ( rp ), intent ( out ) :: del_vec
  real ( rp ), intent ( out ) :: delt
  real ( rp ), intent ( out ) :: delx
  real ( rp ), intent ( out ) :: dely
  real ( rp ), intent ( out ) :: eps
  real ( rp ), intent ( out ) :: gamma
  real ( rp ), intent ( out ) :: gx
  real ( rp ), intent ( out ) :: gy
  integer, intent ( out ) :: imax
  character ( len = * ), intent ( out ) :: infile
  character ( len = * ), intent ( in ) :: inputfile
  integer, intent ( out ) :: itermax
  integer, intent ( out ) :: jmax
  integer, intent ( out ) :: n
  real ( rp ), intent ( out ) :: omega
  character ( len = * ), intent ( out ) :: outfile
  integer, intent ( out ) :: p_bound
  real ( rp ), intent ( out ) :: pos1x
  real ( rp ), intent ( out ) :: pos1y
  real ( rp ), intent ( out ) :: pos2x
  real ( rp ), intent ( out ) :: pos2y
  real ( rp ), intent ( out ) :: pr
  character ( len = * ), intent ( out ) :: problem
  real ( rp ), intent ( out ) :: re
  character ( len = * ), intent ( out ) :: streakfile
  real ( rp ), intent ( out ) :: t_end
  real ( rp ), intent ( out ) :: tau
  real ( rp ), intent ( out ) :: ti
  character ( len = * ), intent ( out ) :: tracefile
  real ( rp ), intent ( out ) :: ui
  character ( len = * ), intent ( out ) :: vecfile
  real ( rp ), intent ( out ) :: vi
  integer, intent ( out ) :: we
  integer, intent ( out ) :: wn
  integer, intent ( out ) :: ws
  integer, intent ( out ) :: ww
  real ( rp ), intent ( out ) :: xlength
  real ( rp ), intent ( out ) :: ylength
!
!  Read the data.
!
  read (5,*) problem
  read (5,'(1a1)') 

  read (5,*) infile
  read (5,*) outfile
  read (5,'(1a1)') 

  read (5,*) xlength
  read (5,*) ylength
  read (5,*) imax
  read (5,*) jmax

  delx = xlength / imax
  dely = ylength / jmax

  read (5,'(1a1)') 

  read (5,*) t_end
  read (5,*) delt
  read (5,*) tau
  read (5,'(1a1)') 

  read (5,*) del_trace
  read (5,*) del_inj
  read (5,*) del_streak
  read (5,*) del_vec
  read (5,'(1a1)') 
  
  read (5,*) vecfile
  read (5,*) tracefile
  read (5,*) streakfile
  read (5,'(1a1)') 

  read (5,*) n
  read (5,*) pos1x
  read (5,*) pos1y
  read (5,*) pos2x
  read (5,*) pos2y
  read (5,'(1a1)') 

  read (5,*) itermax
  read (5,*) eps 
  read (5,*) omega
  read (5,*) gamma
  read (5,*) p_bound
  read (5,'(1a1)') 
  read (5,'(1a1)') 
  read (5,'(1a1)') 

  read (5,*) re
  read (5,*) pr
  read (5,*) beta
  read (5,*) gx
  read (5,*) gy
  read (5,*) ui
  read (5,*) vi
  read (5,*) ti
  read (5,'(1a1)') 

  read (5,*) ww
  read (5,*) we
  read (5,*) wn
  read (5,*) ws
!
!  Print the data.
!
  write (6,*) " problem: ", problem

  write (6,*) " xlength = ", xlength
  write (6,*) " ylength = ", ylength

  write (6,*) " imax = ", imax
  write (6,*) " jmax = ", jmax

  write (6,*) " delx = ", delx
  write (6,*) " dely = ", dely

  write (6,*) " delt = ", delt
  write (6,*) " t_end = ", t_end
  write (6,*) " tau = ", tau

  write (6,*) " del_trace = ", del_trace
  write (6,*) " del_inj = ", del_inj
  write (6,*) " del_streak = ", del_streak
  write (6,*) " del_vec = ", del_vec

  write (6,*) " vecfile: ", vecfile
  write (6,*) " tracefile: ", tracefile
  write (6,*) " streakfile: ", streakfile
  write (6,*) " infile: ", infile
  write (6,*) " outfile: ", outfile

  write (6,*) " n = ", n
  write (6,*) " pos1x = ", pos1x
  write (6,*) " pos1y = ", pos1y
  write (6,*) " pos2x = ", pos2x
  write (6,*) " pos2y = ", pos2y

  write (6,*) " itermax = ", itermax
  write (6,*) " eps = ", eps
  write (6,*) " omega = ", omega
  write (6,*) " gamma = ", gamma

  if ( p_bound /= 1  .and. p_bound /= 2 ) then
    write (6,*) " p_bound must be 1 or 2"
    stop ' init'
  else
    write (6,*) " p_bound = ", p_bound
  end if

  write (6,*) " re = ", re
  write (6,*) " pr = ", pr
  write (6,*) " beta = ", beta
  write (6,*) " gx = ", gx
  write (6,*) " gy = ", gy
  write (6,*) " ui = ", ui
  write (6,*) " vi = ", vi
  write (6,*) " ti = ", ti
!
!  Perform some simple checks on the data.
!
  if ( ww < 1 .or. 4 < ww ) then
    write (6,*) "ww must be 1,2,3, or 4"
    stop ' ww'
  else
    write (6,*) " ww = ", ww
  end if

  if ( we < 1 .or. 4 < we ) then
    write (6,*) "we must be 1,2,3, or 4"
    stop ' we'
  else
    write (6,*) " we = ", we
  end if

  if ( wn < 1 .or. 4 < wn ) then
    write (6,*) "wn must be 1,2,3, or 4"
    stop ' wn'
  else
    write (6,*) " wn = ", wn
  end if

  if ( ws < 1 .or. 4 < ws ) then
    write (6,*) "ws must be 1,2,3, or 4"
    stop ' ws'
  else
    write (6,*) " ws = ", ws
  end if

  if (((ww == 4).and.(we /= 4)) .or. (ww /= 4).and.(we == 4)) then
    write (6,*) "Periodic boundary conditions need ww=we=4"
    stop ' ww or we'
  end if

  if (((ws == 4).and.(wn /= 4)) .OR. (ws /= 4).and.(wn == 4)) then
    write (6,*) "Periodic boundary conditions need ws=wn=4"
    stop ' wn or ws'
  end if

  return
end
subroutine init_uvp ( problem, u, v, p, temp, imax, jmax, ui, vi, ti )

!*******************************************************************************
!
!! INIT_UVP initializes the U, V, P, and TEMP fields.
!
!  Modified:
!
!    27 February 2004
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
!    P(0:IMAX+1,0:JMAX+1), TEMP(0:IMAX+1,0:JMAX+1), the velocity, 
!    pressure, and temperature fields.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real UI, VI, TI, scalars that may be used to initialize 
!    the flow and temperature.
!
  use nrtype

  implicit none

  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( inout ) :: p
  character ( len = * ), intent ( in ) :: problem
  real ( rp ), dimension(0:,0:), intent ( inout ) :: temp
  real ( rp ), intent ( in ) :: ti
  real ( rp ), dimension(0:,0:), intent ( inout ) :: u
  real ( rp ), intent ( in ) :: ui
  real ( rp ), dimension(0:,0:), intent ( inout ) :: v
  real ( rp ), intent ( in ) :: vi

  u(0:imax+1,0:jmax+1) = ui
  v(0:imax+1,0:jmax+1) = vi
  p(0:imax+1,0:jmax+1) = 0.0
  temp(0:imax+1,0:jmax+1) = ti

  if ( problem == "backstep" ) then
    u(0:imax+1,0:jmax/2) = 0.0
  end if

  return
end
subroutine read_bin ( filename, u, v, p, temp, flag, imax, jmax, check )

!*******************************************************************************
!
!! READ_BIN reads binary files of fields U, V, P, and TEMP for use in restarts
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM, 1998.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to read.
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
!    Output, integer CHECK, is 0 if the operation was successful.
!
  use nrtype

  implicit none

  integer, intent ( out ) :: check
  character ( len = * ), intent ( in ) :: filename
  integer ( i2b ), dimension(0:,0:), intent ( out ) :: flag
  integer i_in
  integer, intent ( in ) :: imax
  integer j_in
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( out ) :: p
  real ( rp ), dimension(0:,0:), intent ( out ) :: temp
  real ( rp ), dimension(0:,0:), intent ( out ) :: u
  real ( rp ), dimension(0:,0:), intent ( out ) :: v

  open ( unit = 1, file = filename, status = 'old', form = 'unformatted', &
    iostat = check )

  if ( check == 0 ) then

    read ( 1, iostat = check ) i_in, j_in

    if ( i_in /= imax .or. j_in /= jmax ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'READ_BIN - Fatal error!'
      write ( *, '(a)' ) '  Problem dimension mismatch.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  i_in = ', i_in
      write ( *, '(a,i6)' ) '  imax = ', imax
      write ( *, '(a,i6)' ) '  j_in = ', j_in
      write ( *, '(a,i6)' ) '  jmax = ', jmax
      stop
    end if

    read ( 1, iostat = check ) u
    read ( 1, iostat = check ) v
    read ( 1, iostat = check ) p
    read ( 1, iostat = check ) temp
    read ( 1, iostat = check ) flag

    if ( check /= 0 ) then
      check = 2
    end if
 
  else

    check = 1

  end if

  close ( unit = 1 )

  return
end
subroutine write_bin ( filename, u, v, p, temp, flag, imax, jmax )

!*******************************************************************************
!
!! WRITE_BIN writes binary restart files of U, V, P, and TEMP.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM, 1998.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to write.
!
!    Input, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1), 
!    P(0:IMAX+1,0:JMAX+1), TEMP(0:IMAX+1,0:JMAX+1), the velocity, 
!    pressure, and temperature fields.
!
!    Input, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
  use nrtype

  implicit none

  character ( len = * ), intent ( in ) :: filename
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( in ) :: p
  real ( rp ), dimension(0:,0:), intent ( in ) :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v

  open ( unit = 1, file = filename, status = 'replace', form = 'unformatted' )

  write ( unit = 1 ) imax, jmax

  write ( unit = 1 ) u
  write ( unit = 1 ) v
  write ( unit = 1 ) p
  write ( unit = 1 ) temp
  write ( unit = 1 ) flag

  close ( unit = 1 )

  return
end
subroutine init_flag ( problem, flag, imax, jmax, delx, dely, ibound )

!*******************************************************************************
!
!! INIT_FLAG defines and prints the geometry for the selected case.
!
!  Modified:
!
!    12 March 2004
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
!    Output, integer FLAG(0:IMAX+1,0:JMAX+1), indicates the type of each cell.
!
!    Input, integer IMAX, JMAX, the index of the last computational
!    row and column of the grid.
!
!    Input, real DELX, DELY, the width and height of one cell.
!
!    Output, integer IBOUND, the number of boundary cells.
!
  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( inout ) :: flag
  integer i
  integer, intent ( out ) :: ibound
  integer ihi
  integer ilo
  integer, intent ( in ) :: imax
  integer j
  integer, intent ( in ) :: jmax
  integer :: low
  real ( rp ) :: mx
  real ( rp ) :: my
  character, dimension(0:200) :: one_line
  character ( len = * ), intent ( in ) :: problem
  real ( rp ) :: rad1
  integer :: up
  real ( rp ) :: x
  real ( rp ) :: y

  include 'defs.h'

  if ( 200 < imax ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INIT_FLAG - Fatal error!'
    write ( *, '(a)' ) '  IMAX is too big.'
    stop
  end if
!
!  Initialize the boundary cells.
!
  do i = 0, imax+1
    flag(i,0) = c_b
    flag(i,jmax+1) = c_b
  end do

  do j = 1, jmax
    flag(0,j) = c_b
    flag(imax+1,j) = c_b
  end do
!
!  Initialize the fluid cells.
!
  do i = 1, imax
    do j = 1, jmax
      flag(i,j) = c_f
    end do
  end do
!
!  Specialize to the problem
!
  if ( problem == "backstep" ) then

    do i = 1, jmax
      do j = 1, jmax/2
        flag(i,j) = c_b
      end do
    end do

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "circle" ) then

    mx = 20.0 / 41.0 * jmax * dely
    my = mx
    rad1 = 5.0 / 41.0 * jmax * dely

    do i = 1, imax
      do j = 1, jmax
        x = (i-0.5) * delx
        y = (j-0.5) * dely
        if ((x-mx)*(x-mx)+(y-my)*(y-my) <= rad1 * rad1 ) then
          flag(i,j) = c_b
        end if
      end do
    end do

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "convection" ) then

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "dam" ) then

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "dcavity" ) then

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "drop" ) then

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "fluidtrap" ) then

    do i = 9*imax/22+1, 13*imax/22
      do j = 1, 4*jmax/11
        flag(i,j) = c_b
      end do
      do j = 8*jmax/11+1, jmax
        flag(i,j) = c_b
      end do
    end do

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "molding" ) then

    mx = jmax * dely / 2
    my = jmax * dely / 2
    rad1 = jmax * dely / 6

    do i = 1, imax
      do j = 1, jmax
        x = (i-0.5)*delx
        y = (j-0.5)*dely
        if ((x-mx)*(x-mx)+(y-my)*(y-my) <= rad1 * rad1 ) then
          flag(i,j) = c_b
        end if
      end do
    end do

    write ( *, '(a)' ) '  problem = ', problem
!
!  Flow past an inclined plate
!
  else if ( problem == "plate" ) then

    low = 2*jmax/5
    up = 3*jmax/5
    flag(low,low) = c_b
    flag(low,low+1) = c_b
    flag(up,up-1) = c_b
    flag(up,up) = c_b

    do i = low+1, up-1
      do j = i-1, i+1
        flag(i,j) = c_b
      end do
    end do

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "rayleigh" ) then

    write ( *, '(a)' ) '  problem = ', problem

  else if ( problem == "wave" ) then

    do i = 1, jmax
      do j = 1, jmax/2
        flag(i,j) = c_b
      end do
    end do

    write ( *, '(a)' ) '  problem = ', problem

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INIT_FLAG:'
    write ( *, '(a)' ) '  Unrecognized problem = "' // trim ( problem ) // '".'
    stop

  end if
!
!  Graphical output of geometry.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INIT_FLAG:'
  write ( *, '(a)' ) '  Geometry of the fluid domain: '
  write ( *, '(a)' ) ' '

  do ilo = 0, imax+1, 70

    ihi = min ( ilo + 70 - 1, imax + 1 )

    write ( *, '(a,i6,a,i6)' ) '  Columns ', ilo, ' to ', ihi
    write ( *, '(a)' ) ' '

    do j = jmax+1, 0, -1
      do i = ilo, ihi
        if ( iand ( flag(i,j), c_f ) /= c_f ) then
          one_line(i-ilo) = '*'
        else
          one_line(i-ilo) = '.'
        end if
      end do
      write ( *, '(201a1)' ) one_line(0:ihi-ilo)
    end do
  end do

  write ( *, '(a)' ) ' '
!
!  Flags for boundary cells
!
  ibound = 0

  do i = 1, imax
    do j = 1, jmax

      if ( iand ( flag(i,j), c_f ) /= c_f ) then

        ibound = ibound + 1

      end if

      flag(i,j) = flag(i,j) + ( iand ( flag(i-1,j), c_f ) * b_w  &
                            +   iand ( flag(i+1,j), c_f ) * b_e  &
                            +   iand ( flag(i,j-1), c_f ) * b_s  &
                            +   iand ( flag(i,j+1), c_f ) * b_n ) / c_f

    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INIT_FLAG:'
  write ( *, '(a,i6)' ) ' IBOUND = ', ibound
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FLAG array:'
  write ( *, '(a)' ) ' '

  do ilo = 0, imax+1, 35

    ihi = min ( ilo + 35 - 1, imax + 1 )

    write ( *, '(a,i6,a,i6)' ) '  Columns ', ilo, ' to ', ihi
    write ( *, '(a)' ) ' '

    do j = jmax+1, 0, -1
      write ( *, '(2x,200i2)' ) flag(ilo:ihi,j)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '

  end do

  return
end
