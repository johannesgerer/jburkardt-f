program main

!*****************************************************************************80
!
!! MAIN is the main program for MHD_CONTROL.
!
!  Discussion:
!
!    MHD_CONTROL is the main routine for the MHD control program.
!
!    This program solves a control problem for an electrically conducting
!    fluid on a rectangular domain.  The control is applied to the right hand
!    side of each equation.  The user supplies the target velocities and
!    magnetic field.
!
!    A gradient type iteration is set up where a right hand side is guessed,
!    then the state equations are solved, the functional evaluated, the adjoint
!    problem solved and its solution used to update control.  See the paper
!    by Trenchea et al for a description.
!
!    The equations are time dependent.  A backward Euler approximation is used
!    for time derivatives of the state variables; the adjoint problem is
!    backward in time  with the analogous time derivative.
!
!    Along the boundary, the horizontal and vertical velocities are
!    specified to be zero.  The normal component of the magnetic field
!    is specified to be zero, while the other condition on the magnetic
!    field is the natural condition.
!
!    Note that the program includes a routine that can specify an inflow
!    or outflow velocity boundary condition; however this capability is
!    not currently being used.
!
!    If a parabolic inflow is to be specified at the bottom portion (lower 1/4)
!    of the left boundary and outflow at the portion (upper 1/4) of right
!    boundary where grad u dot n = 0, then you must choose the number
!    of intervals in the Y direction to be a multiple of 4.
!
!    The MHD problem is formulated in terms of primitive variables, the
!    velocity components (U,V), the pressure P, and the magnetic field (B1,B2).
!
!    The target can either be given as an explicit function or input as nodal
!    values from the file 'soln_target.dat'.
!
!    This version uses finite element techniques with piecewise linear
!    functions on triangles to approximate the pressure and quadratics on
!    triangles for the velocity and magnetic field.  This is the Taylor-Hood
!    element.
!
!    This program writes the target solution, adjoint solution, state solution
!    and old and new controls at all timesteps to a file because of storage
!    requirements.
!
!    The term
!      Integral ( div B * div test function )
!    is added to the equation for B since otherwise B is not guaranteed to
!    be exactly weakly divergent free.
!
!  List of routines:
!
!  ADJOINT  solves the adjoint equations
!
!  BANDED (banded solver) & files to increase file name & files to write
!  off arrays.
!
!  BASIS_LINEAR   evaluates linear basis function centered at given node at (x,y)
!
!  BASIS_QUAD  evaluates basis function and x and y derivatives at given point
!
!  GEOMETRY  sets up the geometry for the problem giving x,y coordinates of each point,
!  the node array and the index array
!
!  COST_FNL  evaluates the cost functional given the state and target
!
!  DECOMPOSE takes a 1-d array of unknowns indexed by unknown number and puts into a 2-D array
!  with row being node numnber and column u,v,B1 or B2
!
!  ERROR calculates the L^2 error  between state and target normalized by L^2 error of target
!
!  EVAL  evaluates the velocity and magnetic field at a point (x,y) given the nodal values in a 1-d
!   indexed by unknown number
!
!  EVAL_2d evaluates solution whose nodal values are given in terms of a 2-d array
!
!  L2NORM computes the L^2 norm of a vector (such as control)
!
!  MAG_NSTOKE  solves the Navier Stokes equations coupled with magnetic field equations
!  uses a combination of Simple and Newton iteration
!
!  OUTPUT routine outputs the velocity and magnetic field for plotting
!
!  POST PROCESS - makes pressure satisfy zero mean
!
!  UBDRY, a user-supplied routine which sets boundry values for velocity.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Reference:
!
!    Max Gunzburger, Catalin Trenchea,
!    Analysis and Discretization of an Optimal Control Problem
!    for the Time-Periodic MHD Equations,
!    Journal of Mathematical Analysis and Applications,
!    Volume 308, Number 2, 2005, pages 440-466.
!
  implicit none
!
!  PARAMETERS
!
  real ( kind = 8 ), parameter :: coupling_par = 10.0
  real ( kind = 8 ), parameter :: reynold = 1.0
  real ( kind = 8 ), parameter :: rey_mag = 1.0
  real ( kind = 8 ), parameter :: visc = 1./ reynold

!
!  geometry
   integer ( kind = 4 ), parameter  :: nx =9 , ny =9        ! number of points in x,y directions (excluding midpo
   integer ( kind = 4 ), parameter  :: max_elements =  2*( nx-1)* (ny-1)    ! max number of elements for dimensio
   integer ( kind = 4 ), parameter  :: max_nodes =  (2*nx-1) * (2*ny - 1)   ! max number of nodes for dimensionin
   integer ( kind = 4 ), parameter  :: n_local = 6        ! number of local  nodes per element
   integer ( kind = 4 ), parameter  :: n_quad = 3          ! number of quadrature points for assembly
   integer ( kind = 4 ), parameter  :: nuk = 5             ! number of unknowns (u,v,b1,b2,p)

  real ( kind = 8 ), parameter :: xlength = 1., ylength = 1.      ! size of domain
!
!  control  parameters
!
  integer ( kind = 4 ), parameter  ::  flag_target = 1      !  = 0 => read in nodal values of target; = 1 => targ
  integer ( kind = 4 ), parameter  :: max_steps_control = 40        !  maximum steps in iteration for control
  real ( kind = 8 ) :: alpha1  =  10.        !  coefficient in cost functional of  velocity - target
  real ( kind = 8 ) :: alpha2  =  10.0        !  coefficient in cost functional of   magnetic field -
  real ( kind = 8 ) :: beta1  =  0.01        !  coefficient in cost functional of control (rhs) for v
  real ( kind = 8 ) :: beta2  =  0.01        !  coefficient in cost functional of control (rhs) for m
  real ( kind = 8 ) :: lambda_max = 20.      !  maximum value for stepsize in gradient method
  real ( kind = 8 ) :: tol_grad = 1.e-3    !  convergence tolerance for gradient method
!
!  time stepping parameters
  integer ( kind = 4 ), parameter  ::  steps_time = 10         ! total  number of time steps
  real ( kind = 8 ), parameter ::   dt = 0.005                !  fixed time step

!  nonlinear solver
  integer ( kind = 4 ), parameter  :: max_newton =  15   ! max number of steps for nonlinear solver
  integer ( kind = 4 ), parameter  :: n_simple = 2        ! number of quadrature points for assembly
   real ( kind = 8 ), parameter :: tol_newton = 0.5e-8
!
!  output parameters
  integer ( kind = 4 ), parameter  :: flag_write = 2       ! indicates how much output you want  (0 - minimum, 5
  integer ( kind = 4 )  ::  plot_steps = steps_time/2 + 1
!
!  GEOMETRY VARIABLES
!
   real ( kind = 8 ), dimension(max_elements)  ::  area       ! area of element
   integer ( kind = 4 ), dimension(max_nodes,nuk)  ::  index                ! unknown numbering array (ordering:
   integer ( kind = 4 ) :: n_lband                                       ! lower bandwidth
   integer ( kind = 4 ) :: n_band                                       ! total bandwidth
   integer ( kind = 4 ) :: n_nodes                    ! number of nodes
   integer ( kind = 4 ) :: n_elements                                   ! number of triangles
   integer ( kind = 4 ) :: n_unknowns                  ! number of unknowns
   integer ( kind = 4 ), dimension(max_elements, n_local)  :: node            ! array associating local and globa
   real ( kind = 8 ), dimension(max_nodes)  ::  xc,yc         ! x and y coordinates of nodes
   real ( kind = 8 ), dimension(max_elements, n_quad)  ::  xq, yq         ! x and y quadrature arrays
  real ( kind = 8 )  ::  y_inflow       !    nonzero inflow condition from 0 to y_inflow
  real ( kind = 8 )  ::  y_outflow       !  nonzero outflow condition from   y_outflow to ylength
  integer ( kind = 4 )  ::  y_in_node      !    node number where inflow begins
  integer ( kind = 4 )  ::  y_out_node       !  node number where outflow begins
!
!  CONTROL VARIABLES
!
  real ( kind = 8 ), dimension(max_nodes, 4  )  ::  control_new    !  value of current control (rhs)
  real ( kind = 8 ), dimension(max_nodes, 4  )  ::  control_old    !  value of old control (rhs) at e
  real ( kind = 8 )  :: cost        !    cost functional  at given time step
  real ( kind = 8 )  :: cost_new        ! accumulated cost functional at new gradient iteration
  real ( kind = 8 )  :: cost_old      !  accumulated (at time t^n) cost functional at old gradient it
  real ( kind = 8 )  :: lambda       ! parameter in gradient iteration which -> 0 from 1
  real ( kind = 8 ), dimension (max_nodes, 4 ) :: soln_adj   ! adjoint solution at each node at each
  real ( kind = 8 ), dimension (max_nodes, 4 ) :: soln_state   ! state solution at each node at each


!
!  STATE SOLUTION VARIABLES
!
  real ( kind = 8 ), allocatable, dimension (:) :: f_new   ! array for current Newton iterate of stat
  real ( kind = 8 ), allocatable, dimension (:) :: f_old   ! array for old Newton iterate of state so
  real ( kind = 8 ), allocatable, dimension (:) :: f_prev   ! array for solution at previous time ste
  real ( kind = 8 ), allocatable, dimension (:,:) :: a   ! coefficient matrix stored here using bande
  real ( kind = 8 ), allocatable, dimension (:) :: work  ! work array for solver
  integer ( kind = 4 ), allocatable, dimension (:)  :: ipivot         ! pivot array for solver
!
!  TARGET VARIABLES
!
    real ( kind = 8 ), dimension(max_nodes, 4  )  ::   soln_target       !  stores target velocities
!  real ( kind = 8 ), dimension(max_nodes, 4, steps_time+1 )  ::  opt_control       !  stores control


  character*20 uv_file_name
  character*20 mag_file_name
  character*20 uvinit_file_name
  character*20 maginit_file_name
  character*20 uvfinal_file_name
  character*20 magfinal_file_name
  character*20 uv_tar_file_name
  character*20 mag_tar_file_name
  character action
  character ( len = 80 ) :: control_old_name = 'control_old.dat'
  character ( len = 80 ) :: control_new_name = 'control_new.dat'
  character ( len = 80 ) :: soln_state_name = 'soln_state.dat'
  character ( len = 80 ) :: soln_adj_name = 'soln_adj.dat'
  character ( len = 80 ) :: soln_target_name = 'soln_target.dat'
!
!  LOCAL VARIABLES
!
  integer ( kind = 4 ) :: control_unit_txt
  integer ( kind = 4 )  ::  file_unit_control_old     ! unit number to write on for control
  integer ( kind = 4 )  ::  file_unit_control_new     ! unit number to write on for control
  integer ( kind = 4 )  ::  file_unit_state     ! unit number to write on for state
  integer ( kind = 4 )  ::  file_unit_adj     ! unit numfile_unit_adber to write on for adjoint
  integer ( kind = 4 )  ::  file_unit_target     ! unit number to write on for target
  integer ( kind = 4 )  ::  flag_eqns   ! flag to indicate whether solving state or adjoint; used to indicate spe
  integer ( kind = 4 ) ::  i, idum, ip, i_time, iuk, iuk_b1,  iuk_b2,  k
  integer ( kind = 4 )  ::  n_control    ! counter for number of iterations in gradient algorithm for control
  integer ( kind = 4 )  ::  n_storage   ! length of vector to be written off at each time step  = 4 * number of n
  integer ( kind = 4 )  ::  n_time    ! counter for time steps
  integer ( kind = 4 )  ::  n_time_p1    ! n_time + 1
  real ( kind = 8 ) :: cost_diff,  cost_test
  real ( kind = 8 ) :: dum

  real ( kind = 8 ), dimension(4)   ::  norm_control     ! little l2 norm of optimal control
  real ( kind = 8 ) ::  time_cur
  real ( kind = 8 ) ::  x, xx,  y, yy, zero

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MHD_CONTROL'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Solution of incompressible Navier-Stokes equations'
  write ( *, '(a)' ) '  coupled with magnetic field.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Time dependent equations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Domain:'
  write ( *, '(a,g14.6)' ) '    X length: ', xlength
  write ( *, '(a,g14.6)' ) '    Y length: ', ylength
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Boundary conditions:'
  write ( *, '(a)' ) '    u = v = 0, '
  write ( *, '(a)' ) '    B dot n = 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Target solution given by explicit function '
  write ( *, '(a)' ) '  in subroutine.'

  i = ( ny + ny - 2  ) / 4
  y_in_node = i + 1
  y_out_node = ( 2 * nx - 1 ) * ( 2 * ny - 1 ) - i
!
!  Set up files to output data for plotting.
!
  uvinit_file_name = 'vel_init.txt'
  maginit_file_name = 'mag_init.txt'
  uv_tar_file_name = 'vel_tar.txt'
  mag_tar_file_name = 'mag_tar.txt'
  uvfinal_file_name = 'vel_final000.txt'
  magfinal_file_name = 'mag_final000.txt'

  open ( unit = 7, file = 'mhd_output.txt' )

  open ( unit = 8, file = 'cost_func.txt' )
!
!  SET UP GEOMETRY FOR PROBLEM, WRITE OUT ARRAYS & ALLOCATE ARRAYS
!
  call geometry (flag_write, max_elements, max_nodes, n_local,  &
    n_quad, nuk, nx, ny, xlength, ylength, y_in_node, y_out_node,&
    area, index,  n_lband, n_nodes, node, n_elements, n_unknowns,   &
    xc, yc, xq, yq )
!
!  Write the XY coordinates to a text file.
!
  call xy_write ( 'xy.txt', n_nodes, xc, yc )
!
!  Set position for inflow and outflow boundary locations in y-direction
!  set to 1/4 * ylength at bottom for inflow and 1/4 from top at outflow
!  note:  not used in this version of code
!
  y_inflow = yc ( y_in_node)
  y_outflow = yc(y_out_node )

!  write(*,*)  ' inflow boundary conditions:  u  quadratic from 0 to ', y_inflow
!  write(*,*)  ' outflow boundary conditions: grad u dot n=0  from ', y_outflow, ' to ', ylength
!
  write(*,*)  '                   U = 0 everywhere  '
  write(*,*)  '                   V = 0 everywhere '
  write(*,*)  '             B dot N = 0 (essential) '
  write(*,*)  '            curl B cross N = 0 (natural) '
  write(*,*)  '  '
!
!  Write out geometry information
!
  write(*,*) ' '
  write(*,*) '  Number of nodes is ', n_nodes
  write(*,*) '  Number of elements is ', n_elements
  write(*,*) '  Number of unknowns is ', n_unknowns
  n_band = 2 * n_lband + 1
  write(*,*) '  Total bandwidth is  ', n_band
  write(*,*) ' '
!
!  allocate arrays for assembly and solver
!
  allocate ( a( n_unknowns,  n_band) )   !  allocate space for coefficient matrix
  allocate ( f_new(n_unknowns ) )   !  allocate space for Newton iterate
  allocate ( f_old(n_unknowns ) )   !  allocate space for old Newton iterate
  allocate ( f_prev(n_unknowns ) )   !  allocate space for solution at previous time step
  allocate ( ipivot(n_unknowns ) )   !  allocate space for pivoting array
  allocate (work(n_unknowns ) )   !  allocate space for work array for solver

  n_storage = 4 * n_nodes  !  length of vector to be written off for out of core storage
!
!  Prepare for CONTROL PROBLEM.
!
!  set up files and unit numbers to be used for out of core storage
!
  action = 'c'
  call r8ud_io ( action, control_old_name, file_unit_control_old, idum, &
    n_storage, f_prev )

  call r8ud_io ( action, control_new_name, file_unit_control_new, idum, &
    n_storage, f_prev )

  call r8ud_io ( action, soln_adj_name, file_unit_adj, idum, &
    n_storage, f_prev )

  call r8ud_io ( action, soln_state_name, file_unit_state, idum, &
    n_storage, f_prev )

  call r8ud_io ( action, soln_target_name, file_unit_target, idum, &
    n_storage, f_prev )
!
!  Note:  data is stored starting at t=0 and going up to t=T;  so the value
!  of the variable at final time is in location steps_time + 1.
!
!  Read in the target and optimal controls (rhs) which generated the
!  target solution.  The optimal control is only used for comparison.
!
  if ( flag_target == 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The target solution is already available'
    write ( *, '(a)' ) '  in "soln_target.dat".'

  else

    action = 'w'

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Setting the target solution.'

    do n_time = 1, steps_time + 1

      do i = 1, n_nodes
        x = xc(i)
        y = yc(i)
        call target (x, y, soln_target(i, 1 ), soln_target(i, 2 ),  &
          soln_target(i, 3 ),  soln_target(i, 4 ) )
      end do

      call r8ud_io ( action, soln_target_name, file_unit_target, n_time, &
        n_storage, soln_target )
!
!  Plot the target solution at one time step
!
      if ( n_time == plot_steps+1 ) then
        call output ( flag_write, index, mag_tar_file_name, max_nodes,  &
          n_nodes, soln_target(1:max_nodes,1:4), uv_tar_file_name, xc, yc )
      end if

    end do

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CONTROL PARAMETERS:'
  write(*,*)  '   Alpha1      Alpha2     Beta1       Beta2   '
  write(*,'(4f12.4)' )  alpha1, alpha2, beta1, beta2
  write(*,*)

  write(7,*)  ' CONTROL PARAMETERS '
  write(7,*)  '   alpha1      alpha2     beta1       beta2   '
  write(7,'(4f12.4)' )  alpha1, alpha2, beta1, beta2
  write(7,*)
!
!  BEGIN CONTROL PROBLEM
!
  write(*,*)
  write(*,*) ' '
  write(*,*) '             BEGIN CONTROL PROBLEM'
  write(*,*) ' '
!
!  Choose an initial guess for the controls.
!
  action = 'w'

  do n_time = 1, steps_time + 1
    control_old(1:n_nodes,  1:2  ) = 0.0
    control_old(1:n_nodes,  3:4  ) = 0.0
    call r8ud_io ( action, control_old_name, file_unit_control_old, &
      n_time, n_storage, control_old )
  end do

  write(*,*)
  write(*,*)  ' STEP 0  -   INITIALIZE '  ;   write(*,*)
  write(*,*)  ' Set initial guess for controls to be constant value = 0.0 '
  write(*,*) ' Solve the state equations for the initial control at each time '
  write(*,*)
  write(7,*)
  write(7,*)  ' STEP 0  -   INITIALIZE '  ;   write(*,*)
  write(7,*)  ' Set initial guess for controls to be constant value = 0.0 '
  write(7,*) ' Solve the state equations for the initial control at each time '
  write(*,*)
!
!  STEP 0  -  Solve state equations for initial guess of controls
!        at each time step:
!  1.  solve state equations with initial guess for control
!     (store in f_state and then write to file)
!  2.  accumulate the cost functional  (cost_old )
!
  time_cur = 0.0
!
!  Set the initial condition for the state solution.
!
!  (fprev is used by state solver and soln_state saves nodal values
!  for plotting)
!
  call set_initial ( index, max_nodes,  n_nodes, n_unknowns,  &
    xc, yc,  f_prev, soln_state(1:max_nodes, 1:4  )    )

  action = 'w'
  n_time = 1
  call r8ud_io ( action, soln_state_name, file_unit_state,n_time, &
    n_storage, soln_state )
!
!  Initialize the Newton iterate.
!
  f_old(1:n_unknowns) =  f_prev(1:n_unknowns)

  cost_old = 0.0

  do n_time = 1, steps_time

    n_time_p1 = n_time + 1
    time_cur = time_cur + dt
    write(*,*) ' '
    write(*,*) ' Solving the state equations  at t = ', time_cur; write(*,*)
    write(7,*) ' '
    write(7,*) ' Solving the state equations  at t = ', time_cur; write(*,*)
!
!  Solve state equations using Newton's method at particular time level
!
    action = 'r'    ! get controls from file
    call r8ud_io ( action, control_old_name, file_unit_control_old, &
      n_time_p1, n_storage, control_old )

    call mag_nstoke (  area, control_old(1:n_nodes,  1:4 ),   &
      coupling_par,  dt, flag_write,   f_old,  f_prev, &
      index, max_elements, max_newton, max_nodes,  &
      n_band, n_elements,  n_lband, n_local, n_nodes, node, &
      n_quad, n_simple, nuk, n_unknowns, rey_mag,   time_cur, &
      tol_newton, visc, xc, xq, yc, yq, ylength, y_inflow, y_outflow,  &
      a,  f_new, ipivot, work )

    if ( flag_write >4 ) then
      write(*,*) ' state solution', ( i, f_new(i), i=1,n_unknowns)
    end if
!
!  Store state solution in soln_state by node, unknown and time level
!  Break 1-D array f_new into a 2D-array (soln_state) which gives
!  nodal values of each unknown
!
    flag_eqns = 1      ! indicates solving state eqns

    call decompose (flag_eqns, f_new,  index,  max_nodes, n_nodes,   &
      xc, yc, ylength, y_inflow, y_outflow,         &
      soln_state(1:max_nodes, 1:4  )   )

    action = 'w'    ! write state to file
    write(*,*) '   write state solution to file'
    call r8ud_io ( action, soln_state_name,file_unit_state, &
      n_time_p1, n_storage, soln_state )
!
!  Calculate L2 norm of difference in state solution and target solution
!
!  Read the target from the file.
!
    action = 'r'

    call r8ud_io ( action, soln_target_name, file_unit_target, n_time_p1, &
      n_storage, soln_target )

    call error ( area, flag_write, index, max_elements, max_nodes,  &
      n_elements,   n_local, n_nodes, node, n_quad, n_unknowns, &
      soln_state(1:max_nodes, 1:4 ), soln_target(1:max_nodes, 1:4  ) ,    &
      xc, xq, yc, yq  )
!
!  Write off initial state solution for plotting.
!
    if ( n_time == plot_steps ) then
      call output ( flag_write, index, maginit_file_name, max_nodes, &
        n_nodes, soln_state(1:max_nodes,1:4), uvinit_file_name, xc, yc )
    end if
!
!  Evaluate the cost functional at current time and accumulate.
!
    call cost_fnl ( alpha1, alpha2, area, beta1, beta2,   &
      coupling_par, dt, flag_write, control_old(1:max_nodes,1:4), &
      soln_state(1:max_nodes,1:4), soln_target(1:max_nodes,1:4),   &
      index, max_elements, max_nodes,  &
      n_elements, n_local, n_nodes, node,  &
      n_quad, n_unknowns,   &
      xc, xq, yc, yq, cost )

    cost_old = cost_old  + cost
    if ( n_time == 2 ) then
      cost_test = cost_old
    end if
!
!  Get ready for the next time step.
!
    f_prev(1:n_unknowns) =  f_new(1:n_unknowns)
    f_old(1:n_unknowns) =  f_new(1:n_unknowns)

  end do    !  end loop over timesteps for state equations

  write(8, * ) 0, cost_old

  write(*,*)
  write(*,*)  ' value of cost functional is ',  cost_old
  write(*,*)
  write(7,*)
  write(7,*)  ' value of cost functional is ',  cost_old
  write(*,*)

  write(*,*)
  write(*,*)  ' COMPLETED INITIALIZATION '
  write(*,*)
  write(*,*) ' '
  write(*,*)
  write(*,*)  ' BEGIN GRADIENT ITERATION'
  write(*,*) ' '
!
!  MAIN LOOP OVER ITERATIONS IN GRADIENT ROUTINE
!
  lambda = 0.5
  write(*,*) ' initial value of stepsize parameter is ', lambda
  write(7,*) ' initial value of stepsize parameter is ', lambda

  do n_control = 1, max_steps_control

    write(*,*) ' '
    write(*,*) ' '
    write(*,*) ' Gradient Iteration Number ', n_control
    write(*,*) ' '
    write(*,*) ' '
    write(7,*) ' '
    write(7,*) ' '
    write(7,*) ' Gradient Iteration Number ', n_control
    write(7,*) ' '
    write(7,*) ' '
!
!  STEP 1  -  Solve linear adjoint equations backward
!  at each time step:
!  1.  solve adjoint equations with state solution generated at
!  previous iteration
!  2.  update control
!
    write(*,*)
    write(*,*) ' STEP 1 :  Solve the adjoint equations at each time'
    write(*,*) ' level and update controls '
    write(*,*)
!
!  Set final condition for adjoint equation and update controls there.
!
    f_prev(1:n_unknowns) = 0.0   ! final condition for adjoint equation
    soln_adj(1:max_nodes, 1:4   )  = 0.0
     action = 'w'    ! write adjoint to file at t=T
    call r8ud_io ( action, soln_adj_name, file_unit_adj,steps_time+1, &
      n_storage, soln_adj )

    action = 'r'    ! get old controls from file at final time and update
    call r8ud_io ( action, control_old_name, file_unit_control_old, &
      steps_time+1, n_storage, control_old )
!
!  JVB: Try no control on velocity.
!
!   control_new(1:n_nodes,1:2) = 0.0
!
    control_new(1:n_nodes,1:2) = ( 1.0 - lambda * beta1 ) &
      * control_old(1:n_nodes,1:2)

    control_new(1:n_nodes,3:4) = ( 1.0 - lambda * beta2 ) &
      * control_old(1:n_nodes,3:4)

    action = 'w'    ! write new controls to file at final time
    call r8ud_io ( action, control_new_name, file_unit_control_new, &
                  steps_time+1, n_storage, control_new )
    time_cur = real ( steps_time ) * dt

    do  n_time = 1,  steps_time

      time_cur = time_cur - dt
      i_time = steps_time - n_time +1
!
!  solve linear adjoint equations at time=time_cur
!
      action = 'r'    !read state solution from file
      call r8ud_io ( action, soln_state_name, file_unit_state, i_time, &
        n_storage, soln_state )
      action = 'r'    !read target solution from file
      call r8ud_io ( action, soln_target_name, file_unit_target, i_time , &
        n_storage, soln_target)

      call adjoint (alpha1, alpha2,  area, beta1, beta2,    &
        coupling_par,  dt, flag_write,   f_prev, &
        soln_state(1:max_nodes, 1:4  ), soln_target(1:max_nodes, 1:4  ),  &
        index, max_elements,   max_nodes,  &
        n_band, n_elements,  n_lband, n_local, n_nodes, node, &
        n_quad,   nuk, n_unknowns, rey_mag,  time_cur, &
        visc, xc, xq, yc, yq, ylength,  y_inflow,  y_outflow,    &
        a,  f_new, ipivot, work )

      write(*,*)  'Solved the linear adjoint equations at time = ', time_cur
      write(7,*)  'Solved the linear adjoint equations at time = ', time_cur
!
!  Store adjoint solution in soln_adj by node, unknown and time level
!  Break 1-D array f_new into a 2-array which gives nodal values
!  of each unknown
!
      flag_eqns = 2      ! indicates solving adjoint eqns
      call decompose (flag_eqns, f_new,  index,  max_nodes, n_nodes,   &
        xc, yc, ylength, y_inflow, y_outflow,         &
        soln_adj(1:max_nodes, 1:4   )   )

      action = 'w'    ! write adjoint solution to file
      call r8ud_io ( action, soln_adj_name, file_unit_adj, i_time , &
        n_storage, soln_adj )
!
!  Update the control at current time
!
      action = 'r'    ! get old controls from file at final time and update
      call r8ud_io ( action, control_old_name,  file_unit_control_old, &
        i_time, n_storage, control_old )
!
!  JVB: Try no control on velocity.
!
!   control_new(1:n_nodes,1:2) = 0.0
!
      control_new(1:n_nodes, 1:2 )  =  (1.0 - lambda * beta1) &
        * control_old(1:n_nodes, 1:2 ) - lambda * soln_adj(1:n_nodes, 1:2 )

      control_new(1:n_nodes, 3:4 )  =  (1.0 - lambda * beta2) &
        * control_old(1:n_nodes, 3:4 ) - lambda * soln_adj(1:n_nodes, 3:4 )

      action = 'w'    ! write new controls to file at final time
      call r8ud_io ( action, control_new_name, file_unit_control_new, &
        i_time, n_storage, control_new )
      write(*,*)  'Updated controls '
      write(7,*)  'Updated controls '

      f_prev(1:n_unknowns) = f_new(1:n_unknowns)

    end do

 100   continue
!
!  STEP 2  -  Solve state equations with updated controls
!  at each time step:
!  1.  solve state equations with current control (store in f_state)
!  2.  accumulate the cost functional  (cost_new_acc)
!
    write(*,*)
    write(*,*) ' STEP 2 :  Solve the state equations at each time level'
    write(*,*) ' using new controls '
    write(*,*)
!
!  Set the initial condition for the state solution.
!
    time_cur = 0.0
    call set_initial ( index, max_nodes,  n_nodes, n_unknowns,  &
      xc, yc,  f_prev, soln_state(1:max_nodes, 1:4  )    )

    f_old(1:n_unknowns) =  f_prev(1:n_unknowns)
    cost_new= 0.0

    do n_time = 1,  steps_time

      n_time_p1 = n_time + 1
      time_cur = time_cur + dt
!
!  Solve state equations using Newton's method at particular time level
!
      action = 'r'    ! get controls from file
      call r8ud_io ( action, control_new_name, file_unit_control_new, &
        n_time_p1, n_storage, control_new )
      write(*,*)
      write(*,*)  'Solving state equations at time = ', time_cur ;  write(*,*)

      call mag_nstoke (  area, control_new(1:n_nodes,  1:4 ),   &
        coupling_par,  dt, flag_write,   f_old,  f_prev, &
        index, max_elements, max_newton, max_nodes,  &
        n_band, n_elements,  n_lband, n_local, n_nodes, node, &
        n_quad, n_simple, nuk, n_unknowns, rey_mag,   time_cur, &
        tol_newton, visc, xc, xq, yc, yq, ylength, y_inflow, y_outflow,  &
        a,  f_new, ipivot, work )
!
!  Store state solution in soln_state by node, unknown and time level
!  Break 1-D array f_new into a 2-array which gives nodal values of each unknown
!
      flag_eqns = 1      ! indicates solving state eqns
      call decompose (flag_eqns, f_new,  index,  max_nodes, n_nodes,   &
        xc, yc, ylength, y_inflow, y_outflow,         &
        soln_state(1:max_nodes, 1:4  )   )
      action = 'w'    ! write state solution to file
      call r8ud_io ( action, soln_state_name, file_unit_state,  &
        n_time_p1, n_storage, soln_state )
!
!  Calculate L2 norm of difference in state solution and target solution
!
      action = 'r'    !read target solution from file
      call r8ud_io ( action, soln_target_name, file_unit_target,   &
        n_time_p1 , n_storage, soln_target)

      call error (    area,   flag_write,    &
        index, max_elements,   max_nodes,  &
        n_elements,   n_local, n_nodes, node, n_quad,   n_unknowns, &
        soln_state(1:max_nodes, 1:4  ), soln_target(1:max_nodes, 1:4  ) ,    &
        xc, xq, yc, yq  )
!
!  Evaluate the cost functional at current time and accumulate
!
      call cost_fnl (alpha1, alpha2, area, beta1, beta2,   &
        coupling_par, dt,   flag_write, control_new(1:max_nodes, 1:4 ), &
        soln_state(1:max_nodes, 1:4 ), soln_target(1:max_nodes, 1:4 ),   &
        index, max_elements,   max_nodes,  &
        n_elements,  n_local, n_nodes,  node,  &
        n_quad,    n_unknowns,   &
        xc, xq, yc, yq,    cost )

      cost_new = cost_new  + cost
!
!  Get ready for next time step.
!
      f_prev(1:n_unknowns) =  f_new(1:n_unknowns)

    end do

    write(*,*) ' '
    write(*,*)  ' value of cost functional is ',  cost_new
    write(*,*) ' '
    write(*,*)
    write(7,*) ' '
    write(7,*)  ' value of cost functional is ',  cost_new
    write(7,*) ' '
!
!  STEP 3  -  Check for convergence
!  1.  if cost functional increasing,  don't accept step;
!  reduce lambda by half and recompute controls & state
!  2.  if normalized difference in cost functional < tolerance STOP
!  3.  if convergence criteria not met, but cost functional decreasing,
!  accept iteration and reduce lambda
!
    write(*,*)
    write(*,*)  ' STEP 3 :  Decide whether to accept step  '
    write(*,*)  '     If cost functional decreases, increase stepsize by 50%'
    write(*,*)  '     If cost functional increases, decrease stepsize by half'
    write(*,*)
!
!  Case 1 - cost functional increasing.
!
    if ( cost_new  >= cost_old  ) then  ! cost functional increasing
      lambda = lambda / 2.
      write(*,*) ' cost functional increasing so don''t accept iteration '
      write(*,*) ' cost increased though cost for first 2 steps didn''t'
      write(*,*) '   cut lambda in half to ', lambda
      write(*,*) ' recompute controls  '
      write(7,*) '  cost functional increasing so don''t accept iteration '
      write(7,*) '   cut lambda in half to ', lambda
      write(7,*) ' recompute controls  '
!
!  Recompute  the control from t=0 to t=T with new lambda
!
      do n_time = 1, steps_time +1

        action = 'r'    ! get adjoints
        call r8ud_io ( action, soln_adj_name, file_unit_adj, n_time, &
          n_storage, soln_adj )
        action = 'r'    ! get old controls from file at final time and update
        call r8ud_io ( action, control_old_name,file_unit_control_old, &
          n_time, n_storage, control_old )
!
!  JVB: Try no control on velocity.
!
!       control_new(1:n_nodes,1:2) = 0.0
!
        control_new(1:n_nodes, 1:2 )  =  (1.0 - lambda * beta1) &
          * control_old(1:n_nodes, 1:2 ) - lambda * soln_adj(1:n_nodes, 1:2 )

        control_new(1:n_nodes, 3:4 )  =  (1.0 - lambda * beta2) &
          * control_old(1:n_nodes, 3:4 ) - lambda * soln_adj(1:n_nodes, 3:4 )

        action = 'w'    ! write new controls to file at each time
        call r8ud_io ( action, control_new_name, file_unit_control_new, &
          n_time, n_storage, control_new )
      end do

      write(*,*)  ' recompute  state and try again'
      go to 100  !  recomputed controls so recompute the corresponding state
!
!  Cases 2,3 - cost functional decreasing.
!
    else

      cost_diff = abs ( (cost_new  - cost_old ) / cost_new )
      write(*,*) ' normalized difference in successive costs is ', cost_diff
!
!  Case 3   -    Convergence achieved
!
      if ( cost_diff  <= tol_grad ) then

        write(*,*) 'gradient iteration has converged'
        write(7,*) 'gradient iteration has converged'
        write(8,*) n_control, cost_new
!
!  Write solution at all timesteps for plotting.
!
        do n_time = 1, steps_time + 1

          call file_name_inc ( uvfinal_file_name )
          call file_name_inc ( magfinal_file_name )

          action = 'r'

          call r8ud_io ( action, soln_state_name, file_unit_state, &
            n_time, n_storage, soln_state )

          call output ( flag_write, index, magfinal_file_name, &
            max_nodes, n_nodes, soln_state(1:max_nodes,1:4),   &
            uvfinal_file_name, xc, yc )

        end do
!
!  Calculate and write off L2 norm of control at each timestep when converged.
!
        call get_unit ( control_unit_txt )

        open ( unit = control_unit_txt, file = 'control.txt' )

        do n_time = 1, steps_time + 1

          action = 'r'
          call r8ud_io ( action, control_new_name, file_unit_control_new, &
            n_time, n_storage, control_new )
          call l2norm ( area, control_new, flag_write, index, &
            max_elements, max_nodes, n_elements, n_local, n_nodes, node, &
            n_quad,   n_unknowns,  xc, xq, yc, yq,  norm_control  )

          write(control_unit_txt,'(5(2x,g14.6))' ) &
            real(n_time-1)*dt, norm_control(1:4)

        end do

        close ( unit = control_unit_txt )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MHD_CONTROL:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )

        stop     ! gradient iteration converged

      else

        write(*,*)
        write(*,*) '  cost functional decreasing so accept step, set lambda to'
        lambda = 1.5 * lambda
        if ( lambda >= lambda_max  ) lambda = lambda_max
        write(*,*)   lambda
        write(*,*)
        write(7,*)
        write(7,*) '  cost functional decreasing so accept step, set lambda to'
        write(7,*)   lambda
        write(7,*)
        write(8,*) n_control, cost_new

        do n_time = 1, steps_time +1

          action = 'r'    ! get new controls from file
          call r8ud_io ( action, control_new_name, file_unit_control_new,  &
            n_time, n_storage, control_new )

          action = 'w'    ! write new controls to file for old controls
          call r8ud_io ( action, control_old_name, file_unit_control_old, &
            n_time, n_storage, control_new )
        end do

        cost_old = cost_new

      end if

    end if

  end do    !   end loop over gradient iteration for control
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MHD_CONTROL'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine adjoint ( alpha1, alpha2, area, beta1, beta2, coupling_par, dt,   &
  flag_write, f_prev, f_state, f_target, index, max_elements, max_nodes,  &
  n_band, n_elements, n_lband, n_local, n_nodes, node, n_quad, nuk, &
  n_unknowns, rey_mag, time_cur, visc, xc, xq, yc, yq, ylength, y_inflow, &
  y_outflow, a, f_adj, ipivot, work )

!*****************************************************************************80
!
!! ADJOINT solves the linear adjoint equations for the control problem.
!
!  Discussion:
!
!    The equations are linear, using Taylor-Hood elements for velocity/pressure
!    and quadratics on triangles for the magnetic field B.
!
!    The f_prev  array contains the solution at previous time step (t^n+1)
!
!    The f_adj array contains the right hand side initially and then the
!    adjoint solution.
!
!    f_state contains the state solution from the current value of controls
!    - 2 d array
!
!    f_target contains the target solution  - 2 d array
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Reference:
!
!    Max Gunzburger, Catalin Trenchea,
!    Analysis and Discretization of an Optimal Control Problem
!    for the Time-Periodic MHD Equations,
!    Journal of Mathematical Analysis and Applications,
!    Volume 308, Number 2, 2005, pages 440-466.
!
  implicit none
!
!  INPUT
!
   integer ( kind = 4 )  :: flag_write          !  flag for determining how much is written out
   integer ( kind = 4 )  :: max_nodes          ! max number of nodes
   integer ( kind = 4 ), dimension(max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
   integer ( kind = 4 )  :: max_elements       ! max number of elements
   integer ( kind = 4 ) :: n_band               ! total bandwidth and leading dimension of a
   integer ( kind = 4 ) :: n_elements                 ! number of triangles
   integer ( kind = 4 )  ::  n_local          ! number of local  nodes per element
   integer ( kind = 4 ) :: n_lband                             ! lower bandwidth
   integer ( kind = 4 )  :: n_nodes           ! number of global nodes
   integer ( kind = 4 ), dimension(max_elements,*)  :: node        ! array associating local and global nodes
   integer ( kind = 4 )  :: n_quad            ! number of quadrature points for assembly
   integer ( kind = 4 ) :: n_unknowns                  ! number of unknowns
   integer ( kind = 4 )  :: nuk               ! number of unknowns - (3) u,v,p

   real ( kind = 8 ) :: alpha1, alpha2, beta1, beta2   ! parameters in cost functional
   real ( kind = 8 ), dimension(*)  ::  area       ! area of element
   real ( kind = 8 )  ::   coupling_par               ! coupling parameter
   real ( kind = 8 )  ::   dt                  ! delta t
   real ( kind = 8 )  ::  f_prev   ! solution of adjoint equations at time t^(n+1)
   real ( kind = 8 ), dimension(max_nodes, *) :: f_state          ! state solution at given time leve
   real ( kind = 8 ), dimension(max_nodes, *) :: f_target         ! target solution at given time lev
   real ( kind = 8 ), dimension(2) :: mag_prev
   real ( kind = 8 )  ::  rey_mag                  ! magnetic Reynolds number
    real ( kind = 8 ) :: time_cur     ! current time
   real ( kind = 8 ), dimension(2)  :: vel_prev
   real ( kind = 8 ) ::  visc     ! viscosity = 1/Reynolds number
   real ( kind = 8 ), dimension(*)  ::  xc,yc      ! x and y coordinates of nodes
   real ( kind = 8 ), dimension(max_elements, *)  ::  xq, yq         ! x and y quadrature arrays
   real ( kind = 8 ) :: ylength,  y_inflow,  y_outflow  ! not really used but needed for calling argu

!
!  OUTPUT
!
  real ( kind = 8 ), dimension(n_unknowns, *) ::  a      ! coefficient matrix
  real ( kind = 8 ), dimension(*) :: f_adj           ! adjoint solution (and rhs initially)
  integer ( kind = 4 ), dimension(*)  ::  ipivot              ! pivot array for solver
  real ( kind = 8 ), dimension(*) :: work           ! work array for solver
!
!  FUNCTIONS

   real ( kind = 8 ) ::  basis_linear

!
!  LOCAL
!
   integer ( kind = 4 ) :: info   ! flag for solver
   integer ( kind = 4 ) :: iter   ! counter for nonlinear iteration
   integer ( kind = 4 ) ::  i,  id, ip, ipp, iq, iqq, iquad, it, iuk, iukk, iuse,j , nrhs
   real ( kind = 8 )  :: par_sim    ! parameter for simple iteration ( = 0 => simple iteration, =0 =>
    real ( kind = 8 ) ::  aij,  ar, x, y,  yy
  real ( kind = 8 ) :: curl_b   !  curl of magnetic field
   real ( kind = 8 ), dimension(4) :: dum    ! dummy array
   real ( kind = 8 ) ::  test, testx, testy, trial, trialx, trialy, test_p, trial_p
   real ( kind = 8 ) ::  diff, normf
   real ( kind = 8 ) ::  rhs
   real ( kind = 8 )  ::  temp

   real ( kind = 8 ), dimension(4) :: state    ! value of  state unknowns u,v,B1,B2 at quadrature poi
   real ( kind = 8 ), dimension(4) :: statex    ! value of x derivative of state unknowns  at quadrat
   real ( kind = 8 ), dimension(4) :: statey    ! value of y derivative of  state unknowns t quadratu
   real ( kind = 8 ), dimension(4) :: target    ! value of target unknowns u,v,B1,B2 at quadrature po
!
!  Initialize.
!
  a(1:n_unknowns,1:n_band) = 0.0
  f_adj(1:n_unknowns) = 0.0
!
!  Indicate that we are evaluating the adjoint equations.
!
  id = 2
!
!  Matrix assembly triangle by triangle
!
  do it = 1, n_elements

    ar = area(it) / 3.0D+00

    do iquad = 1, n_quad    ! loop over quadrature points
      y = yq(it,iquad)
      x = xq(it,iquad)
!
!  Evaluate adjoint solution at time t^(n+1) where we are solving at time t^n
!
      call eval  (f_prev, id,  index, it, max_elements, max_nodes,  &
        n_local, node, x, y, xc, yc, ylength, y_inflow, y_outflow, &
        vel_prev, temp, temp, mag_prev, temp, temp )!
!
!  Evaluate state solution at quadrature point
!
      call eval_2d  (f_state,  it, max_elements, max_nodes,  &
        n_local, node, x, y, xc, yc, state, statex, statey )

      curl_b = statex(4) - statey(3)
!
!  Evaluate target solution at quadrature point
!
      call eval_2d  (f_target,  it, max_elements, max_nodes,  &
        n_local, node, x, y, xc, yc, target, dum, dum )

      do iq = 1, n_local  ! loop over local nodes for test functions

        ip = node(it,iq)
        call basis_quad (it, iq, max_elements, node, x, y, xc, yc, &
          test, testx, testy)

        if ( iq <= 3) then
          test_p = basis_linear( it, iq, 1, max_elements, node, x, y, xc, yc )
        end if

        do iuk = 1, nuk

          i =  index(ip,iuk)

          if ( i > 0 .and. i /= n_unknowns) then

            if ( iuk == 1)   &
              f_adj(i) = f_adj(i) &
                + ar * ( test * alpha1 * (state(1) - target(1) )   &
                + vel_prev(1)*test   / dt    )


            if ( iuk == 2)   &
              f_adj(i) = f_adj(i) &
                + ar* (  test *  alpha1 * (state(2) - target(2) )   &
                + vel_prev(2) * test  / dt  )

            if ( iuk == 3)   &
              f_adj(i) = f_adj(i) &
                +  ar * (  test *  alpha2 * (state(3) - target(3) )    &
                +  mag_prev(1) * test  / dt       )

            if ( iuk == 4)   &
              f_adj(i) = f_adj(i) &
                +  ar * (     test * alpha2 * (state(4) - target(4) ) &
                +  mag_prev(2) * test / dt   )

            do iqq = 1, n_local   ! loop over local nodes for trial functions

              ipp = node(it,iqq)
              yy = yc(ipp)

              call basis_quad ( it, iqq, max_elements, node, x, y, xc, &
                yc, trial, trialx, trialy)

              if ( iqq <= 3) then
                trial_p = basis_linear( it, iqq, 1, max_elements, node, x, &
                  y, xc, yc )
              end if

              do iukk = 1, nuk

                j =  index(ipp,iukk)

                if ( j > 0 ) then

                  aij = 0.0

                  if  (iuk == 1 ) then   ! equation for w

                     if (iukk == 1)  then   ! unknown for w1

                       aij = test * trial / dt  &
                         +  visc * (testx * trialx + testy * trialy )  &
                         + test *   trial  * statex(1)    &
                         +  trial * ( state(1) * testx + state(2) * testy )

                     else if (iukk == 2) then   ! unknown for w2

                       aij =   test * trial * statex(2)

                     else if (iukk == 3) then   ! unknown for D1

                       aij =  - coupling_par * trial &
                         * ( testy * state(4) + test * statey(4) )

                     else if (iukk == 4) then   ! unknown for D2

                       aij = coupling_par *  trial &
                         * ( testx * state(4) + test * statex(4) )

                     else    ! unknown for p

                       aij = -testx * trial_p

                     end if   ! over iukk

                   else if (iuk == 2) then   ! equation for w2

                     if (iukk == 1) then

                       aij =     test * trial   * statey(1)  ! unknown for w1

                     else if (iukk.eq.2) then   ! unknown for w2

                       aij =  test * trial / dt &
                         + visc * ( testx * trialx + testy * trialy) &
                         +  trial * (  state(1) * testx + state(2) * testy) &
                         + test * trial * statey(2)

                     else if (iukk == 3) then   ! unknown for D1

                       aij = coupling_par * trial * ( testy * state(3) &
                         + test * statey(3) )

                     else if (iukk == 4) then   ! unknown for D2

                       aij =  -coupling_par * trial * ( testx * state(3) &
                         + test * statex(3) )

                     else                ! unknown for r

                       aij = -testy * trial_p

                     end if

                        else if (iuk == 3) then   ! equation for D1

                          if (iukk == 1) then
                              aij =   -state(4) * testy * trial  ! unknown for w1

                            else if (iukk.eq.2) then   ! unknown for w1

                              aij =   - test  * curl_b * trial - state(3)*testy * trial

                            else if (iukk == 3) then   ! unknown for D1
                              aij =   test * trial / dt  +  (1./rey_mag) * testy * trialy &
                                      + statey(2)*test * trial + state(2) * testy * trial  &
                                      +(1./rey_mag) * testx * trialx    !  div div term
!
!
                            else if (iukk == 4) then   ! unknown for D2
                              aij = - (1./rey_mag) *  testy  *trialx   &
                                     -statex(2) * test * trial - state(2) * testx * trial   &
                                    +(1./rey_mag) * testx * trialy    !  div div term


                            else                ! unknown for r
                              aij = 0.0
                          end if    ! over iukk

                          else if (iuk == 4) then   ! equation for D2

                          if (iukk == 1) then     !unknown for w1
                              aij =   test* trial * curl_b +    state(4) * testx * trial

                            else if (iukk.eq.2) then   ! unknown for w2
                              aij =  - state(3) * testx * trial

                            else if (iukk == 3) then   ! unknown for D1
                              aij =  -(1./ rey_mag) * testx *  trialy   &
                                     - trial * ( statey(1) * test + state(1) * testy ) &
                                   +(1./rey_mag) * testy * trialx    !  div div term

                            else if (iukk == 4) then   ! unknown for D2
                              aij =  test * trial / dt + (1./ rey_mag) * testx * trialx   &
                                      - trial * ( statex(1) * test + state(1) * testx )    &
                                     +(1./rey_mag) * testy * trialy    !  div div term

                            else                ! unknown for p
                              aij =  0.0
                          end if    ! over iukk

                        else                      ! equation for pressure
                          if (iukk == 1  ) aij = -trialx * test_p
                          if (iukk == 2  ) aij = - trialy * test_p

                                  end if   !  end if over iuk=1,2,3,4,5

                      if ( j > 0 ) then
                          iuse = j - i + n_lband + 1
                          a(i, iuse) = a(i, iuse) + aij*ar    !  add entries into matrix
!
                        else if  ( j < 0 ) then     !  add terms to rhs side due to inhomogeneous bc

                       end if

                  endif   ! end if over j>= 0

              end do    ! end of loop over unknowns for trial
            end do      ! end of loop over local nodes for trial

            end if          ! end if for i > 0 conditional

          end do        ! end of loop over unknowns for test
        end do          ! end of loop over local nodes for test
      end do            ! end of loop over quadrature points
    end do              ! end of loop over triangles
!
!  To avoid singularity of the pressure system, the last pressure
!  is simply assigned a value of 1.
!
    f_adj( n_unknowns) = 1.
    a (  n_unknowns, n_band ) = 1.
!
!  Factor the matrix and solve the linear system.
!
  call banded ( a, f_adj, work, ipivot, n_unknowns, n_lband, n_lband, &
    n_unknowns )
!
!  Post process the pressure by computing the mean and subtracting
!  from nodal values of p
!
  call post_process ( area, flag_write, f_adj, index, max_elements,    &
    max_nodes,   n_elements, n_nodes, node )

  return
end
subroutine banded(a,b,wk,is,n,m1,m2,maxu)

!*****************************************************************************80
!
!! BANDED is the banded equation solver.
!
!  Discussion:
!
!    BANDED ELIMINATION CODE - USES PARTIAL PIVOTING, BUT
!    DOES NOT REQUIRE ANY ADDITIONAL STORAGE - KNOWS THAT
!    THE LAST PRESSURE UNKNOWN MAY BE SET ARBITRARILY
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Parameters:
!
!    Input/output, a    = banded matrix, with non-zero bands stored as columns;
!              the (i,j) entry in dense mode is stored in
!              the (i,m1+1+j-i) entry of banded storage, where
!              m1=lower half-bandwidth; thus the main diagonal
!              is stored in the column number (m1+1)
!
!    Input/output, B(N).  On input, the right hand side vector of the
!    linear system.  On output, the solution vector for the linear system.
!
!    Workspace, real ( kind = 8 ) WK(*).
!
!    Workspace, integer IS(*).
!
!    Input, integer N, the number of equations.
!
!    Input, m1   = lower half-bandwidth
!
!    Input, m2   = upper half-bandwidth; number of nonzero diagonals is
!              thus (m1+m2+1)
!
!    Input, integer MAXU, the leading dimension of A.
!
  implicit none

  integer ( kind = 4 ) ::  maxu
  integer ( kind = 4 ) :: n

  real ( kind = 8 ), dimension (maxu, *) ::   a
  real ( kind = 8 ), dimension(*)  ::  b
  real ( kind = 8 ), dimension(*)  ::  wk      ! work array

  integer ( kind = 4 ), dimension (*) ::   is     ! work array for pivots
  integer ( kind = 4 )  ::  m1        ! lower bandwidth
  integer ( kind = 4 )  ::  m2        ! upper bandwidth
!
!  Local variables
  real ( kind = 8 ) ::  p, q,  tols, x
  integer ( kind = 4 )  ::  i, ii, im1, j, je, js,  k, kp1, l, m, nm1, np1

  tols=.00001d0
  m=m1+m2+1

  do   i=1,n

    p=abs(a(i,1))

    do  j=2,m

          q=abs(a(i,j))
          if(q.gt.p) p=q

    end do

        wk(i)=p

  end do

  do  i=1,m1

    je=m2+i
    js=m1+1-i
    a(i,1:je) = a(i, 1+js: js +je)
    je=je+1
    a(i, je:m ) = 0.0
  end do

    nm1=n-1
    l=m1
    is(1:n) = 0

    do   k=1,nm1
      kp1=k+1
      if(l.lt.n) l=l+1
      i=k
      p=abs(a(k,1))/wk(k)

      do   j=kp1,l
        q=abs(a(j,1))/wk(j)

        if(q.gt.p) then
          p=q
          i=j
        end if
      end do

      if( abs(p) > tols) then

          if( i /= k) then

              do   j=1,m
                x=a(k,j)
                a(k,j)=a(i,j)
                a(i,j)=x
              end do

              x=b(k)
              b(k)=b(i)
              b(i)=x

          end if

          do  i=kp1,l
            x=a(i,1)/a(k,1)

            do   j=2,m
              a(i,j-1)=a(i,j)-x*a(k,j)
            end do

            a(i,m)=0.d0
            b(i)=b(i)-x*b(k)
          end do

        else

          is(k)=1

          do   i=kp1,l
            a(i, 1:m-1) = a(i,2:m)
            a(i,m)=0.d0
          end do

      end if

    end do

      if(dabs(a(n,1)).le.tols) is(n)=1
      np1=n+1
      l=1

      do   ii=1,n
        i=np1-ii

        if(is(i) == 1) then
        b(i)=1.
      else
        im1=i-1
        x=b(i)

        if(i < n) then
            do k=2,l
              x=x-a(i,k)*b(im1+k)
            end do
        end if

        b(i)=x/a(i,1)

        end if

        if(l.lt.m) l=l+1

  end do

  return
end
function basis_linear (it, iq, id, max_elements, node, x, y, xc, yc )

!*****************************************************************************80
!
!! BASIS_LINEAR evaluates the linear basis functions associated with pressure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
!  INPUT
  real ( kind = 8 ) :: x, y    ! point at which basis function is evaluated
  real ( kind = 8 ), dimension(*) :: xc, yc
  integer ( kind = 4 ) :: it   ! triangle basis function is in
  integer ( kind = 4 ) :: iq   ! local node basis function is centered at
  integer ( kind = 4 ) :: id   ! flag for determining basis function or its derivative
  integer ( kind = 4 )  :: max_elements          ! max number of nodes
  integer ( kind = 4 ), dimension(max_elements,*)  :: node        ! array associating local and global nodes
!
!  OUTPUT
  real ( kind = 8 ) basis_linear
!
!  Local variables
  integer ( kind = 4 ) :: iq1, iq2, iq3, i1, i2, i3
  real ( kind = 8 ) :: d

  iq1 = iq
  iq2 = mod(iq,3)+1
  iq3 = mod(iq+1,3)+1
  i1 = node(it,iq1)
  i2 = node(it,iq2)
  i3 = node(it,iq3)
  d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1) )

  if ( id == 1) then  ! evaluates value of basis function

     basis_linear = 1.0+((yc(i2)-yc(i3))*(x -xc(i1))+(xc(i3)-xc(i2))*( y -yc(i1)))/d

    else if (id == 2) then  ! evaluates value of x derivative of basis function

      basis_linear = (yc(i2)-yc(i3))/d

    else if (id == 3) then

      basis_linear = (xc(i3)-xc(i2))/d

    else

      write (*,*) 'bsp - fatal error!'
      write (*,*) 'unknown value of id=',id
      stop
  end if

  return
end
subroutine basis_quad ( it, in, max_elements, node, x, y, xc, yc,  bb, bx, by )

!*****************************************************************************80
!
!! BASIS_QUAD evaluates a quadratic basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
!  INPUT
  real ( kind = 8 )  :: x, y    ! point at which we are evaluating basis function
  real ( kind = 8 ), dimension(*) :: xc, yc    ! coordinate arrays
  integer ( kind = 4 ) :: it     ! triangle point is in
  integer ( kind = 4 ) :: in     ! local node the basis function is centered at
  integer ( kind = 4 )  :: max_elements       ! max number of elements
  integer ( kind = 4 ), dimension(max_elements,*)  :: node        ! array associating local and global nodes
!
!  OUTPUT
  real ( kind = 8 ) :: bb, bx, by   ! value of basis function and x and y derivatives at pt x,y
!
!  Local variables
  integer ( kind = 4 ) ::   inn, in1, in2, in3, i1, i2, i3, j1, j2, j3
  real ( kind = 8 ) :: d, t, c, s

  if ( in <= 3) then   ! compute basis function centered at vertex
        in1 = in
        in2 = mod(in,3)+1
        in3 = mod(in+1,3)+1
        i1 = node(it,in1)
        i2 = node(it,in2)
        i3 = node(it,in3)

        d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

        t = 1.+((yc(i2)-yc(i3))*(x-xc(i1))+(xc(i3)-xc(i2))*(y-yc(i1)))/d

        bb = t*(2.*t-1.)
        bx = (yc(i2)-yc(i3))*(4.*t-1.)/d
        by = (xc(i3)-xc(i2))*(4.*t-1.)/d

    else        ! compute basis function centered at midpoint of side

        inn = in-3
        in1 = inn
        in2 = mod(inn,3)+1
        in3 = mod(inn+1,3)+1
        i1 = node(it,in1)
        i2 = node(it,in2)
        i3 = node(it,in3)
        j1 = i2
        j2 = i3
        j3 = i1
        d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

        c = (xc(j2)-xc(j1))*(yc(j3)-yc(j1))-(xc(j3)-xc(j1))*(yc(j2)- yc(j1))

        t = 1.+((yc(i2)-yc(i3))*(x-xc(i1))+(xc(i3)-xc(i2))*(y-yc(i1)))/d

        s = 1.+((yc(j2)-yc(j3))*(x-xc(j1))+(xc(j3)-xc(j2))*(y-yc(j1)))/c

        bb = 4.*s*t
        bx = 4.*(t*(yc(j2)-yc(j3))/c+s*(yc(i2)-yc(i3))/d)
        by = 4.*(t*(xc(j3)-xc(j2))/c+s*(xc(i3)-xc(i2))/d)

  end if

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c .eq. ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine cost_fnl (alpha1, alpha2, area, beta1, beta2, coupling_par, dt, &
  flag_write, f_control, f_state, f_target, index, max_elements, max_nodes, &
  n_elements, n_local, n_nodes, node, n_quad, n_unknowns, xc, xq, yc, yq, cost )

!*****************************************************************************80
!
!! COST_FNL evaluates the cost functional at a given time for a specific target
!  and a given state
!
!  The f_control contains nodal values of  control at given time level
!  The f_state contains nodal values of  state solution at given time level
!  The f_target contains nodal values of target solution at given time level
!
!  cost is the output functional
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Reference:
!
!    Max Gunzburger, Catalin Trenchea,
!    Analysis and Discretization of an Optimal Control Problem
!    for the Time-Periodic MHD Equations,
!    Journal of Mathematical Analysis and Applications,
!    Volume 308, Number 2, 2005, pages 440-466.
!
  implicit none
!
!  INPUT
!
   integer ( kind = 4 )  :: flag_write          !  flag for determining how much is written out
   integer ( kind = 4 )  :: max_nodes          ! max number of nodes
   integer ( kind = 4 ), dimension(max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
   integer ( kind = 4 )  :: max_elements       ! max number of elements
   integer ( kind = 4 ) :: n_elements                 ! number of triangles
   integer ( kind = 4 )  ::  n_local          ! number of local  nodes per element
   integer ( kind = 4 )  :: n_nodes           ! number of global nodes
   integer ( kind = 4 ), dimension(max_elements,*)  :: node        ! array associating local and global nodes
   integer ( kind = 4 )  :: n_quad            ! number of quadrature points for assembly
   integer ( kind = 4 ) :: n_unknowns                  ! number of unknowns

   real ( kind = 8 ) ::  alpha1, alpha2, beta1, beta2   !  parameters in cost functional
   real ( kind = 8 ), dimension(*)  ::  area       ! area of element
   real ( kind = 8 )  ::   coupling_par               ! coupling parameter
   real ( kind = 8 )  ::   dt                  ! delta t
   real ( kind = 8 ), dimension(max_nodes,*) :: f_control          ! nodal values of control
   real ( kind = 8 ), dimension(max_nodes,*) :: f_state          ! nodal values of state
   real ( kind = 8 ), dimension(max_nodes,*) :: f_target          ! nodal values of target


   real ( kind = 8 ), dimension(*)  ::  xc,yc      ! x and y coordinates of nodes
    real ( kind = 8 ), dimension(max_elements, *)  ::  xq, yq         ! x and y quadrature arrays

!
!  OUTPUT
!
  real ( kind = 8 )  ::  cost      ! cost functional at given time level

!
!
!  LOCAL
!
   integer ( kind = 4 ) ::  iquad, it, iuk
   real ( kind = 8 ) ::  ar, x, y
   real ( kind = 8 ), dimension(4)  :: dum   ! dummy array
   real ( kind = 8 ), dimension(4)  ::  norm_control   !  l2 norm squared of control
   real ( kind = 8 ), dimension(4)  ::  norm_state      !l2 norm squared of state - target
   real ( kind = 8 ), dimension(4) :: state    ! value of  state unknowns u,v,B1,B2 at quadrature poi
   real ( kind = 8 ), dimension(4) :: target    ! value of target unknowns u,v,B1,B2 at quadrature po
   real ( kind = 8 ), dimension(4) :: control  ! value of  controls   at quadrature point
!
! initialize
!
  norm_state(1:4) = 0.0
  norm_control(1:4) = 0.0
!
!  calculate L2 norm over each element and combine
!
  do it = 1, n_elements

    ar = area(it) / 3.0D+00

    do  iquad = 1, n_quad

      y = yq(it,iquad)
      x = xq(it,iquad)
!
!  Evaluate state solution at quadrature point
!
      call eval_2d  (f_state,  it, max_elements, max_nodes,  &
        n_local, node, x, y, xc, yc, state, dum, dum )
!
!  Evaluate target solution at quadrature point
!
      call eval_2d  (f_target,  it, max_elements, max_nodes,  &
        n_local, node, x, y, xc, yc, target, dum, dum )
!
!  Evaluate controls at quadrature point
!
      call eval_2d  (f_control,  it, max_elements, max_nodes,  &
        n_local, node, x, y, xc, yc, control, dum, dum )

      do iuk = 1, 4

        norm_state(iuk) = norm_state(iuk) + ar * ( state(iuk) - target(iuk) )**2
        norm_control(iuk) = norm_control(iuk) + ar * ( control(iuk) )**2

      end do

    end do

  end do

  cost = 0.5 * alpha1                * dt * ( norm_state(1)   + norm_state(2) )   &
       + 0.5 * alpha2 * coupling_par * dt * ( norm_state(3)   + norm_state(4) )   &
       + 0.5 * beta1                 * dt * ( norm_control(1) + norm_control(2) ) &
       + 0.5 * beta2  * coupling_par * dt * ( norm_control(3) + norm_control(4) )

  return
end
subroutine decompose ( flag, f_1d, index, max_nodes, n_nodes, xc, yc, &
  ylength, y_inflow, y_outflow, f_2d )

!*****************************************************************************80
!
!! DECOMPOSE breaks a 1-d vector into a 2-d vector by unknown and node
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
!  INPUT
   integer ( kind = 4 )  :: flag    !  if flag= 1 then we are solving state eqns;  =2 adjoint (needed for bc)
   real ( kind = 8 ), dimension(*) :: f_1d     !  1-d array of unknowns

  integer ( kind = 4 ) :: max_nodes  ! max number of nodes

  integer ( kind = 4 ), dimension(max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
  integer ( kind = 4 ) :: n_nodes  !   number of nodes
  real ( kind = 8 ), dimension(*) :: xc, yc   !x,  y coordinates
  real ( kind = 8 )  :: ylength,  y_inflow, y_outflow
!
!  OUTPUT
   real ( kind = 8 ), dimension(max_nodes, *) :: f_2d     !  2-d array of unknowns
!
!  Local variables
  integer ( kind = 4 ) :: i, ip,   iuk
  real ( kind = 8 ) ::  yy, ubc

!
!  FUNCTIONS
  real ( kind = 8 ) :: ubdry

  f_2d(1:max_nodes, 1:4) = 0.0

  do ip = 1, n_nodes   !  loop over nodes

    do iuk = 1, 4      !  loop over unknowns

      i = index(ip, iuk)

      if ( i > 0 ) then
          f_2d( ip, iuk) = f_1d (i)
        else if ( i < 0 .and. flag == 1 ) then  ! using state solution, not adjoint
          yy = yc(ip)
          ubc = ubdry ( i, yy, ylength, y_inflow, y_outflow  )
          f_2d( ip, iuk ) = ubc
      end if

    end do

  end do

  return
end
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit .eq. -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit .eq. 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 .le. digit .and. digit .le. 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine error ( area, flag_write, index, max_elements, max_nodes, &
  n_elements, n_local, n_nodes, node, n_quad, n_unknowns, f_state, f_target, &
  xc, xq, yc, yq )

!*****************************************************************************80
!
!! ERROR computes the L^2 error in the state and target solutions.
!
!  INPUT:
!      f_state   -  state solution in 2-D array dimensioned (node #, unknown)
!      f_target   -  target solution in 2-D array dimensioned (node #, unknown)
!      geometry info
!
!  OUTPUT:
!      writes norms to screen and output file
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
!  INPUT
!
   integer ( kind = 4 )  :: flag_write          !  flag for determining how much is written out
   integer ( kind = 4 )  :: max_nodes          ! max number of nodes
   integer ( kind = 4 ), dimension(max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
   integer ( kind = 4 )  :: max_elements       ! max number of elements
   integer ( kind = 4 ) :: n_elements                 ! number of triangles
   integer ( kind = 4 )  ::  n_local          ! number of local  nodes per element
   integer ( kind = 4 )  :: n_nodes           ! number of global nodes
   integer ( kind = 4 ), dimension(max_elements,*)  :: node        ! array associating local and global nodes
  integer ( kind = 4 )  :: n_quad            ! number of quadrature points for assembly
   integer ( kind = 4 ) :: n_unknowns                  ! number of unknowns

   real ( kind = 8 ), dimension(*)  ::  area       ! area of element
   real ( kind = 8 ), dimension(max_nodes, *) :: f_state          ! state solution at given time leve
   real ( kind = 8 ), dimension(max_nodes, *) :: f_target         ! target solution at given time lev
   real ( kind = 8 ), dimension(*)  ::  xc,yc      ! x and y coordinates of nodes
    real ( kind = 8 ), dimension(max_elements, *)  ::  xq, yq         ! x and y quadrature arrays

!
!  LOCAL
!
   integer ( kind = 4 ) ::  i,  id, ip, ipp, iq, iqq, iquad, it, iuk, iukk, iuse,j , nrhs

    real ( kind = 8 ) ::  aij,  ar, x, y,  yy
   real ( kind = 8 ), dimension(4) :: dum    ! dummy array
   real ( kind = 8 ), dimension(4)  ::  l2_norm, l2_error
   real ( kind = 8 ) :: l2_vel,  l2_mag

   real ( kind = 8 ), dimension(4) :: state    ! value of  state unknowns u,v,B1,B2 at quadrature poi
   real ( kind = 8 ), dimension(4) :: statex    ! value of x derivative of state unknowns  at quadrat
   real ( kind = 8 ), dimension(4) :: statey    ! value of y derivative of  state unknowns t quadratu
   real ( kind = 8 ), dimension(4) :: target    ! value of target unknowns u,v,B1,B2 at quadrature po

!
! initialize
!
    l2_norm(1:4) = 0.0
    l2_error(1:4) = 0.0

!
!  matrix assembly triangle by triangle
!
    do   it = 1, n_elements     ! loop over elements
      ar = area(it) /3.0D+00
!
      do  iquad = 1, n_quad    ! loop over quadrature points
        y = yq(it,iquad)
        x = xq(it,iquad)

!
!  Evaluate state solution at quadrature point
!
         call eval_2d  (f_state,  it, max_elements, max_nodes,  &
              n_local, node, x, y, xc, yc,  &
                 state, statex, statey )
!
!  Evaluate target solution at quadrature point
!
         call eval_2d  (f_target,  it, max_elements, max_nodes,  &
              n_local, node, x, y, xc, yc,  &
                 target, dum, dum )
!
          do iuk = 1, 4  ! loop over unknowns
            l2_norm(iuk) = l2_norm(iuk) + ar * ( target(iuk) ) **2
            l2_error(iuk) = l2_error(iuk) +  ar * ( state(iuk) - target(iuk) ) **2
          end do
      end do
    end do

    l2_vel =  sqrt(    l2_error(1)   +  l2_error(2)   ) / sqrt (  l2_norm(1)   +  l2_norm(2)   )
    l2_mag =  sqrt(    l2_error(3)   +  l2_error(4)   ) / sqrt (  l2_norm(3)   +  l2_norm(4)   )

    l2_norm(1:4) = sqrt ( l2_norm(1:4) )
    l2_error(1:4) = sqrt ( l2_error(1:4) )


    write(*,*)
    write(*,*) '  Normalized L^2 errors in state and target solutions  '
    write(*,*)  '     u            v           B1                 B2'
    write(*,'(4e14.5)' ) ( l2_error(iuk)/ l2_norm(iuk) , iuk = 1, 4 )
    write(*,*)   '  L2 errors in velocity         L2 errors in magnetic field  '
    write(*,'(2e20.5)' ) l2_vel,  l2_mag
    write(7,*)
    write(7,*) '  Normalized L^2 errors in state and target solutions  '
    write(7,*)  '     u            v           B1                 B2'
    write(7,'(4e14.5)' ) ( l2_error(iuk)/ l2_norm(iuk) , iuk = 1, 4 )
    write(7,*)   '  L2 errors in velocity         L2 errors in magnetic field  '
    write(7,'(2e20.5)' ) l2_vel,  l2_mag

  return
end
subroutine eval_2d  (f,   it,  max_elements, max_nodes,  &
  n_local, node, x, y, xc, yc, soln, solnx, solny )

!*****************************************************************************80
!
!! EVAL_2D evaluates vector dimensioned (nodes, 4)  at the given point (x,y)
!   Also evaluates derivatives
!
!  eval  evaluates the given 2-d array f(1:n_nodes,1:4) at a point (x,y)
!  for each unknown
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none

  integer ( kind = 4 ) :: max_elements
  integer ( kind = 4 ) :: max_nodes

  real ( kind = 8 ), dimension(max_nodes,*) :: f
  integer ( kind = 4 ) :: it
  integer ( kind = 4 ) :: n_local
  integer ( kind = 4 ), dimension(max_elements,*)  :: node
  real ( kind = 8 ) :: x, y
  real ( kind = 8 ), dimension(*) :: xc, yc
  real ( kind = 8 ), dimension(4) :: soln, solnx, solny
  integer ( kind = 4 ) ::   ip, iq, iuk
  real ( kind = 8 ) :: bb, bx, by

  soln(1:4) = 0.0
  solnx(1:4) = 0.0
  solny(1:4) = 0.0

  do iq = 1, n_local

    call basis_quad ( it, iq, max_elements, node, x, y, xc, yc, bb, bx, by )
    ip = node(it,iq)

    do iuk = 1, 4
      soln(iuk) = soln(iuk) + f(ip, iuk) * bb
      solnx(iuk) = solnx(iuk) + f(ip, iuk) * bx
      solny(iuk) = solny(iuk) + f(ip, iuk) * by
    end do

  end do

  return
end
subroutine eval ( g, id, index, it, max_elements, max_nodes, n_local, &
  node, x, y, xc, yc, ylength, y_inflow, y_outflow, vel, velx, vely, &
  mag, magx, magy )

!*****************************************************************************80
!
!! EVAL evaluates the velocities and the magnetic field
!   and their X and Y derivatives at a
!  given point in a particular triangle.
!
!   EVAL evaluates the velocity and magnetic field at the given point (x,y) if id = 1
!              evaluates solution of adjoint equation if id = 2
!              solution vector is 1-d array
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
!  INPUT
   real ( kind = 8 ), dimension(*) :: g     !  array of unknowns
   integer ( kind = 4 )  ::  id   !  flag to indicate if state (id=1)  or adjoint (id=2) unknowns are to be evalu
  integer ( kind = 4 ) :: it  ! triangle the point is in
  integer ( kind = 4 ) :: max_elements  ! max number of elements
  integer ( kind = 4 ) :: max_nodes  ! max number of nodes
  integer ( kind = 4 ) :: n_local  ! number of local nodes
  integer ( kind = 4 ), dimension(max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
  integer ( kind = 4 ), dimension(max_elements,*)  :: node
  real ( kind = 8 ) :: x, y      !  point to evaluate at
  real ( kind = 8 ), dimension(*) :: xc, yc   !x,  y coordinates
  real ( kind = 8 )  :: ylength,  y_inflow, y_outflow
!
!  OUTPUT
  real ( kind = 8 ), dimension(2) :: vel, velx, vely
  real ( kind = 8 ), dimension(2) :: mag, magx, magy
!
!  Local variables
  integer ( kind = 4 ) :: i, ip, iq, iuk,   iun, kk
  real ( kind = 8 ) ::  yy, ubc
  real ( kind = 8 ) :: bb, bx, by
!
!  FUNCTIONS
  real ( kind = 8 ) :: ubdry

    vel(1:2) = 0. ;   vely(1:2) = 0;   velx(1:2) = 0.
    mag(1:2) = 0. ;   magy(1:2) = 0;   magx(1:2) = 0.
!
  do iq = 1, n_local

    call basis_quad ( it, iq, max_elements, node, x, y, xc, yc,  bb, bx ,by)
    ip = node(it,iq)
    yy = yc(ip)

    if ( id == 1 ) then  ! evaluate state velocities which have some inhomogeneous Dirichlet bc

        do iuk  = 1,2
          iun =  index(ip,iuk )

          if  (iun > 0 ) then
            vel(iuk) = vel(iuk) + bb*g(iun)
            velx(iuk) = velx(iuk) + bx*g(iun)
            vely(iuk) = vely(iuk) + by*g(iun)
          else if (iun < 0) then
            ubc = ubdry ( iun, yy, ylength, y_inflow, y_outflow  )
            vel(iuk) = vel(iuk)+bb*ubc
            velx(iuk) = velx(iuk)+bx*ubc
            vely(iuk) = vely(iuk)+by*ubc
          end if

        end do

      else

        do iuk  = 1,2
          iun =  index(ip,iuk )

          if  (iun > 0 ) then
            vel(iuk) = vel(iuk) + bb*g(iun)
            velx(iuk) = velx(iuk) + bx*g(iun)
            vely(iuk) = vely(iuk) + by*g(iun)
          end if
        end do

      end if    !  end for velocity calculations

    do iuk  = 1,2  ! evaluate magnetic field
      iun =  index(ip, iuk+2 )

      if  (iun > 0 ) then
          mag(iuk) = mag(iuk) + bb*g(iun)
          magx(iuk) = magx(iuk) + bx*g(iun)
          magy(iuk) = magy(iuk) + by*g(iun)
        else if (iun < 0) then   !index for B should never be less than 0
          write(*,*)  '  error in evaluation routine for bc for  B'
          stop
      end if

    end do

  end do

  return
end
subroutine eval_rhs ( g, it, max_elements, n_local, node, x, y, xc, yc, rhs )

!*****************************************************************************80
!
!! EVAL_RHS evaluates the rhs at the given point (x,y)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
!  INPUT
  real ( kind = 8 ), dimension( *) :: g     ! rhs at each node
  integer ( kind = 4 ) :: it  ! triangle the point is in
  integer ( kind = 4 ) :: max_elements  ! max number of elements
  integer ( kind = 4 ) :: n_local  ! number of local nodes
  integer ( kind = 4 ), dimension(max_elements,*)  :: node
  real ( kind = 8 ) :: x, y      !  point to evaluate at
  real ( kind = 8 ), dimension(*) :: xc, yc   !x,  y coordinates

  real ( kind = 8 )  :: rhs
  integer ( kind = 4 ) ::   ip, iq
  real ( kind = 8 ) :: bb, bx, by

  rhs = 0.0

  do iq = 1, n_local

    call basis_quad ( it, iq, max_elements, node, x, y, xc, yc, bb, bx, by )
    ip = node(it,iq)
    rhs = rhs + bb * g(ip)

   end do

  return
end
subroutine file_copy ( old_file_name, new_file_name, ierror )

!*****************************************************************************80
!
!! FILE_COPY makes a copy of a file.
!
!  Discussion:
!
!    The file is assumed to be sequential access, with variable
!    length records.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OLD_FILE_NAME, the name of the file
!    to be copied.
!
!    Input, character ( len = * ) NEW_FILE_NAME, the name of the copy of
!    the file.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, the file names are the same.
!    2, a free unit number could not be found for the old file.
!    3, the routine could not open the old file.
!    4, a free unit number could not be found for the new file.
!    5, the routine could not open the new file.
!
  implicit none

  logical file_exist
  logical file_is_open
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 256 ) line
  character ( len = * ) new_file_name
  integer ( kind = 4 ) new_unit
  character ( len = * ) old_file_name
  integer ( kind = 4 ) old_unit

  ierror = 0
!
!  Does the original file exist?
!
  if ( .not. file_exist ( old_file_name ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old file does not exist.'
    return
  end if
!
!  Is the original file open?
!
  if ( file_is_open ( old_file_name ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old file is open.'
    write ( *, '(a)' ) '  It must be closed before it can be copied.'
    return
  end if
!
!  Make sure the file names aren't the same.
!
  if ( new_file_name == old_file_name ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old and new file names are identical.'
    return
  end if
!
!  Does the new file exist?
!
  if ( file_exist ( new_file_name ) ) then

    if ( file_is_open ( new_file_name ) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
      write ( *, '(a)' ) '  A file is already open with the new name.'
      write ( *, '(a)' ) '  It must be closed before it can be overwritten.'
      return
    end if

    call file_delete ( new_file_name )

  end if
!
!  At this point:
!
!    The old file exists, and is not open.
!    The new file does not exist, and has a different name.
!
!  Open the old file.
!
  call get_unit ( old_unit )

  if ( old_unit == 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not get a unit number for the old file.'
    return
  end if

  open ( unit = old_unit, file = old_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not open the old file:'
    write ( *, '(4x,a)' ) '"' // trim ( old_file_name ) // '".'
    return
  end if
!
!  Open the new file.
!
  call get_unit ( new_unit )

  if ( new_unit == 0 ) then
    ierror = 4
    close ( unit = old_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free unit for the copy file.'
    return
  end if

  open ( unit = new_unit, file = new_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 5
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not open the new file:'
    write ( *, '(4x,a)' ) '"' // trim ( new_file_name ) // '".'
    close ( unit = old_unit )
    return
  end if
!
!  Read an old line, write a new line.
!
  do

    read ( old_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    write ( new_unit, '(a)' ) trim ( line )

  end do

  close ( unit = old_unit )

  endfile ( unit = new_unit )
  close ( unit = new_unit )

  return
end
subroutine file_delete ( file_name )

!*****************************************************************************80
!
!! FILE_DELETE deletes a named file if it exists.
!
!  Discussion:
!
!    You might want to call this routine to get rid of any old copy
!    of a file, before trying to open a new copy with the OPEN argument:
!      status = 'new'.
!
!    It's not always safe to open a file with ' STATUS = 'UNKNOWN' '.
!    For instance, on the SGI, the most recent version of the FORTRAN
!    compiler seems to go crazy when I open an unformatted direct
!    access file this way.  It creates an enormous file (of somewhat
!    random size).  The problem goes away if I delete any old copy
!    using this routine, and then open a fresh copy with
!    'STATUS = 'NEW'.  It's a scary world.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
  implicit none

  logical file_exist
  logical file_is_open
  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical, parameter :: verbose = .false.
!
!  Does the file exist?
!
  if ( .not. file_exist ( file_name ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  There is no file of the given name.'
    return
  end if
!
!  Is the file open?
!
  if ( file_is_open ( file_name ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  The file is currently open.'
    write ( *, '(a)' ) '  It must be closed before it can be deleted.'
    return
  end if
!
!  Get a free unit number.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  A free FORTRAN unit could not be found.'
    return
  end if

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE:'
    write ( *, '(a)' ) '  Deleting "' // trim ( file_name ) // '".'
  end if

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    return
  end if

  close ( unit = iunit, status = 'delete' )

  return
end
function file_exist ( file_name )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  character ( len = * ) file_name
  logical file_exist

  inquire ( file = file_name, exist = file_exist )

  return
end
function file_is_open ( file_name )

!*****************************************************************************80
!
!! FILE_IS_OPEN reports whether a file (specified by filename) is open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_IS_OPEN, is TRUE if the file is open.
!
  implicit none

  character ( len = * ) file_name
  logical file_is_open

  inquire ( file = file_name, opened = file_is_open )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next file name in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
      implicit none
!
      character c
      logical ch_is_digit
      character*(*) file_name
      integer ( kind = 4 ) i
      integer ( kind = 4 ) lens
!
      lens = len ( file_name )

      do i = lens, 1, -1

        c = file_name(i:i)

        if ( ch_is_digit ( c ) ) then

          call digit_inc ( c )

          file_name(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

        end if

      end do

  return
end
subroutine geometry (flag_write, max_elements, max_nodes, n_local,  &
  n_quad, nuk, nx, ny, xlength, ylength, y_in_node, y_out_node, &
  area, index,  n_lband, n_nodes, node, n_elements, n_unknowns,   &
  xc, yc, xq, yq )

!*****************************************************************************80
!
!! GEOMETRY sets up x,y coordinate arrays, node array associating local node number to
!  global node number, index array which associates an unknown to a node,
!  total number of nodes, elements and unknowns,
!  info for quadrature and bandwidth
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
!  INPUT
!
  integer ( kind = 4 ) :: flag_write          !  flag for determining how much is written out
  integer ( kind = 4 )  :: max_elements       ! max number of elements
  integer ( kind = 4 )  :: max_nodes          ! max number of nodes
  integer ( kind = 4 )  ::  n_local          ! number of local  nodes per element
  integer ( kind = 4 )  :: n_quad            ! number of quadrature points for assembly
  integer ( kind = 4 )  :: nuk               ! number of unknowns - (3) u,v,p
  integer ( kind = 4 )  :: nx,  ny           ! number of points in x,y directions (excluding midpoints)
  integer ( kind = 4 ) ::  y_in_node, y_out_node   !  node numbers where in and outflow conditions begin
  real ( kind = 8 ) :: xlength, ylength  ! size of domain
  real ( kind = 8 ) :: y_inflow,  y_outflow   ! y locations of nonzero inflow and outflow
!
!  OUTPUT
!
   real ( kind = 8 ), dimension(*)  ::  area       ! area of element
   integer ( kind = 4 ), dimension(max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
   integer ( kind = 4 ) :: n_lband                             ! lower bandwidth
   integer ( kind = 4 ) :: n_nodes                    ! number of nodes
   integer ( kind = 4 ) :: n_elements                 ! number of triangles
   integer ( kind = 4 ) :: n_unknowns                  ! number of unknowns
   integer ( kind = 4 ), dimension(max_elements,*)  :: node        ! array associating local and global nodes
   real ( kind = 8 ), dimension(*)  ::  xc,yc      ! x and y coordinates of nodes
   real ( kind = 8 ), dimension(max_elements, *)  ::  xq, yq         ! x and y quadrature arrays
!
!  LOCAL
!
   integer ( kind = 4 ) :: icnt, jcnt   ! counters to determine whether at midpoint or vertex
   integer ( kind = 4 ) :: element_cntr    ! counter for elements
   integer ( kind = 4 ) :: node_cntr       ! counter for nodes
   integer ( kind = 4 ) :: unknown_cntr    ! counter for unknowns
   integer ( kind = 4 ) ::  i, ic, ip, ipp, ip1, ip2, ip3, iq, iqq, iuk, iukk,  it, j, jc
   integer ( kind = 4 ) :: n_col, n_row   ! number of nodes in row or column
   real ( kind = 8 ) :: deltax2, deltay2    ! half of delta x and delta y
   real ( kind = 8 ) :: xx, x1, x2, x3, yy, y1, y2, y3
   real ( kind = 8 ) :: xm1, xm2, xm3, ym1, ym2, ym3

  n_row = nx + (nx -1 )
  n_col = ny + (ny - 1)

  deltax2 = 0.5 * (xlength / real (nx-1 ) )  !  half of delta x
  deltay2 = 0.5 * (ylength / real (ny-1 ) )  !  half of delta y

  node_cntr = 0      ! counter for nodes
  element_cntr = 0   ! counter for  elements
  unknown_cntr = 0   ! counter for unknowns
!
!  construct grid numbering vertically
!

  do  ic = 1, n_row
    icnt = mod(ic,2)
    xx = (ic-1) * deltax2
!
    do  jc = 1, n_col
      jcnt = mod(jc,2)
      yy = (jc-1) * deltay2
      node_cntr = node_cntr + 1
      xc( node_cntr) = xx
      yc( node_cntr) = yy
!
      if ( icnt == 1 .and. jcnt == 1 ) then      ! at vertex

          if ( ic /= n_row .and. jc /= n_col ) then  !  not at top or outflow

              element_cntr = element_cntr + 1
              node(element_cntr,1) = node_cntr
              node(element_cntr,2) = node_cntr + 2
              node(element_cntr,3) = node_cntr + 2 * n_col + 2
              node(element_cntr,4) = node_cntr + 1
              node(element_cntr,5) = node_cntr + n_col + 2
              node(element_cntr,6) = node_cntr + n_col + 1

              element_cntr = element_cntr + 1
              node(element_cntr,1) = node_cntr
              node(element_cntr,2) = node_cntr + 2 * n_col + 2
              node(element_cntr,3) = node_cntr + 2 * n_col
              node(element_cntr,4) = node_cntr + n_col + 1
              node(element_cntr,5) = node_cntr + 2 * n_col + 1
              node(element_cntr,6) = node_cntr + n_col

         end if

     end if
!
!  set unknown numbers at node
!  boundary conditions   u,v = 0 on all of boundary
!  B dot n = 0 on boundary; other condition natural
!
        if ( jc == 1 .or. jc == n_col ) then

          index(node_cntr, 1 ) = 0
          index(node_cntr, 2 ) = 0
          if ( ic == 1 .or.  ic == n_row ) then
            index(node_cntr, 3) = 0
          else
            unknown_cntr  = unknown_cntr  + 1
            index(node_cntr, 3) = unknown_cntr
          end if
          index(node_cntr, 4 ) = 0

        else

          if ( ic ==1 ) then
            index(node_cntr, 1 ) = 0
            index(node_cntr, 2 ) = 0
            unknown_cntr  = unknown_cntr  + 1
            index(node_cntr, 3) = 0
            index(node_cntr, 4 ) = unknown_cntr
          else if ( ic == n_row ) then
            index(node_cntr, 1 ) = 0
            index(node_cntr, 2 ) = 0
            unknown_cntr  = unknown_cntr  + 1
            index(node_cntr, 3) = 0
            index(node_cntr, 4 ) = unknown_cntr
          else
            unknown_cntr  = unknown_cntr  + 2
            index(node_cntr, 1 ) = unknown_cntr  - 1
            index(node_cntr, 2 ) = unknown_cntr
            unknown_cntr  = unknown_cntr  + 2
            index(node_cntr, 3) = unknown_cntr  - 1
            index(node_cntr, 4 ) = unknown_cntr
          end if

        end if
!
!  Set unknowns for pressure only at vertices (not midpoints)
!
        if ( icnt == 1 .and.  jcnt == 1 ) then
          unknown_cntr  = unknown_cntr  + 1
          index(node_cntr,5 ) = unknown_cntr
        else
          index(node_cntr,5) = 0
       end if

     end do

   end do
!
!  Record the number of nodes, n_nodes, number of triangles,
!  n_elements, and the number of unknowns,
!
  n_nodes = node_cntr
  n_elements = element_cntr
  n_unknowns = unknown_cntr

!
!  write out arrays if flag_write > 2
!
  if ( flag_write > 2 ) then

      write(*,*)' '
      write(*,*)'    I      XC           YC      INDX U & V, B1, B2, P'
      write(*,*)' '
      do i = 1, n_nodes
        write (*,'(i5,2f12.5,5i5)') &
          i, xc(i), yc(i), index(i,1:5)
       end do
   end if

  if ( flag_write > 2 ) then
      write(*,*)' '
      write(*,*)'    IT    NODE(IT,I),I=1,6)'
      write(*,*)' '
      do it = 1, n_elements
        write (*,'(7i6)') it, node(it,1:6)
       end do

   end if
!
!  Set midpoint quadrature rule information and areas of triangles
!

 do  it = 1, n_elements
    ip1 = node(it,1)
    ip2 = node(it,2)
    ip3 = node(it,3)
    x1 = xc(ip1)
    x2 = xc(ip2)
    x3 = xc(ip3)
    y1 = yc(ip1)
    y2 = yc(ip2)
    y3 = yc(ip3)
    xm1 = .5*(x1+x2)
    xm2 = .5*(x2+x3)
    xm3 = .5*(x3+x1)
    ym1 = .5*(y1+y2)
    ym2 = .5*(y2+y3)
    ym3 = .5*(y3+y1)
    xq(it,1) = xm1
    xq(it,2) = xm2
    xq(it,3) = xm3
    yq(it,1) = ym1
    yq(it,2) = ym2
    yq(it,3) = ym3
    area(it) = 2.*abs(ym1*(xm2-xm3)+ym2*(xm3-xm1)+ym3*(xm1-xm2))
   end do
!
!  Compute the half band width
!
  n_lband = 0

  do  it = 1, n_elements

    do   iq = 1, n_local
      ip = node(it,iq)

      do   iuk = 1, nuk

        i = index(ip,iuk)

        if ( i > 0 ) then

            do iqq = 1, n_local
              ipp = node(it,iqq)

              do  iukk = 1, nuk

        j = index(ipp,iukk)
                if ( j >= i )  n_lband = max( n_lband, j-i)

       end do   ! loop over iukk
      end do    ! loop over iqq

    end if

    end do    !  loop over iuk
    end do    !  loop over iq
   end do    ! loop over it

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A 'free' FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which are
!    often reserved for consol I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine l2norm ( area, f, flag_write, index, max_elements, max_nodes,  &
  n_elements, n_local, n_nodes, node, n_quad, n_unknowns, xc, xq, yc, yq, &
  norm )

!*****************************************************************************80
!
!! L2NORM calculates the L^2 norm of f.
!
!  INPUT:
!    f (# of nodes, 4)  - the nodal values of function
!    geometry variables
!
! OUTPUT:
!   norm(1)  -  Big L^2 norm of first two components of f
!   norm(2)  -  Big L^2 norm of last two components of f
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
!  INPUT
!
   integer ( kind = 4 )  :: flag_write          !  flag for determining how much is written out
   integer ( kind = 4 )  :: max_nodes          ! max number of nodes
   integer ( kind = 4 ), dimension(max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
   integer ( kind = 4 )  :: max_elements       ! max number of elements

   integer ( kind = 4 ) :: n_elements                 ! number of triangles
   integer ( kind = 4 )  ::  n_local          ! number of local  nodes per element

   integer ( kind = 4 )  :: n_nodes           ! number of global nodes
   integer ( kind = 4 ), dimension(max_elements,*)  :: node        ! array associating local and global nodes
    integer ( kind = 4 )  :: n_quad            ! number of quadrature points for assembly
   integer ( kind = 4 ) :: n_unknowns                  ! number of unknowns



   real ( kind = 8 ), dimension(*)  ::  area       ! area of element

   real ( kind = 8 ), dimension(max_nodes, *) :: f          ! nodal values of function


   real ( kind = 8 ), dimension(*)  ::  xc,yc      ! x and y coordinates of nodes
   real ( kind = 8 ), dimension(max_elements, *)  ::  xq, yq         ! x and y quadrature arrays

!  OUTPUT
!
  real ( kind = 8 ), dimension(4)  :: norm        ! L^2 norms of components of f
!
!  LOCAL
!

    integer ( kind = 4 ) ::   iquad, it, iuk

    real ( kind = 8 ) ::   ar, x, y


   real ( kind = 8 ), dimension(4) :: state    ! value of components of f at quadrature point
   real ( kind = 8 ), dimension(4) :: statex    ! value of x derivative of f at quadrature point
   real ( kind = 8 ), dimension(4) :: statey    ! value of y derivative of  f at quadrature point
!
! initialize

    norm(1:4) = 0.0
!
!  matrix assembly triangle by triangle
!
    do it = 1, n_elements

      ar = area(it) / 3.0D+00

      do  iquad = 1, n_quad

        y = yq(it,iquad)
        x = xq(it,iquad)
!
!  Evaluate f at quadrature point
!
        call eval_2d ( f, it, max_elements, max_nodes, n_local, node, &
          x, y, xc, yc, state, statex, statey )

        do iuk = 1, 4
          norm(iuk) = norm(iuk) + ar * ( state(iuk) )**2
        end do

      end do

    end do

    norm(1:4) = sqrt ( norm(1:4) )

  return
end
subroutine mag_nstoke ( area, rhs_cur, coupling_par, dt, flag_write, f_old, &
  f_prev, index, max_elements, max_newton, max_nodes, n_band, n_elements,  &
  n_lband, n_local, n_nodes, node, n_quad, n_simple, nuk, n_unknowns, &
  rey_mag, time_cur, tol_newton, visc, xc, xq, yc, yq, ylength, y_inflow, &
  y_outflow, a, f_new, ipivot, work )

!*****************************************************************************80
!
!! MAG_NSTOKE solves the Navier Stokes equations coupled with magnetic field
!  using Taylor-Hood for velocity/pressure and quadratics on triangles for magnetic field B
!
!  The f_old array contains the old Newton iterate.
!  The f_prev  array contains the solution at previous time step.
!
!  The f_new array contains the right hand side initially and then the
!  current iterate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Reference:
!
!    Max Gunzburger, Catalin Trenchea,
!    Analysis and Discretization of an Optimal Control Problem
!    for the Time-Periodic MHD Equations,
!    Journal of Mathematical Analysis and Applications,
!    Volume 308, Number 2, 2005, pages 440-466.
!
  implicit none
!
!  INPUT
!
   integer ( kind = 4 )  :: flag_write          !  flag for determining how much is written out
   integer ( kind = 4 )  :: max_nodes          ! max number of nodes
   integer ( kind = 4 ), dimension(max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
   integer ( kind = 4 )  :: max_elements       ! max number of elements
   integer ( kind = 4 )  :: max_newton         ! max number of Newton/simple iterations
   integer ( kind = 4 ) :: n_band               ! total bandwidth and leading dimension of a
   integer ( kind = 4 ) :: n_elements                 ! number of triangles
   integer ( kind = 4 )  ::  n_local          ! number of local  nodes per element
   integer ( kind = 4 ) :: n_lband                             ! lower bandwidth
   integer ( kind = 4 )  :: n_nodes           ! number of global nodes
   integer ( kind = 4 ), dimension(max_elements,*)  :: node        ! array associating local and global nodes
     integer ( kind = 4 )  :: n_quad            ! number of quadrature points for assembly
   integer ( kind = 4 )  :: n_simple          ! number of simple iterations
   integer ( kind = 4 ) :: n_unknowns                  ! number of unknowns
   integer ( kind = 4 )  :: nuk               ! number of unknowns - (3) u,v,p



   real ( kind = 8 ), dimension(*)  ::  area       ! area of element

   real ( kind = 8 )  ::   coupling_par               ! coupling parameter
   real ( kind = 8 )  ::   dt                  ! delta t
   real ( kind = 8 ), dimension(*) :: f_old          ! old Newton iterate
   real ( kind = 8 ), dimension(*) :: f_prev         ! solution at previous time step
   real ( kind = 8 )  ::  rey_mag                  ! magnetic Reynolds number
   real ( kind = 8 ), dimension(max_nodes,  *)  ::  rhs_cur   !  rhs control at each node at current
   real ( kind = 8 ) :: time_cur     ! current time
   real ( kind = 8 ) ::  tol_newton     ! tolerance for converence of Newton
   real ( kind = 8 ) ::  visc     ! viscosity = 1/Reynolds number
   real ( kind = 8 ), dimension(*)  ::  xc,yc      ! x and y coordinates of nodes
    real ( kind = 8 ), dimension(max_elements, *)  ::  xq, yq         ! x and y quadrature arrays
   real ( kind = 8 ) :: ylength,  y_inflow,  y_outflow
!
!  OUTPUT
!
  real ( kind = 8 ), dimension(n_unknowns, *) ::  a      ! coefficient matrix
  real ( kind = 8 ), dimension(*) :: f_new           ! current Newton iterate (and rhs)
  integer ( kind = 4 ), dimension(*)  ::  ipivot              ! pivot array for solver
  real ( kind = 8 ), dimension(*) :: work           ! work array for solver
!
!  FUNCTIONS

   real ( kind = 8 ) ::  basis_linear
   real ( kind = 8 ) ::  ubdry
!
!  LOCAL
!
   integer ( kind = 4 ) :: info   ! flag for solver
   integer ( kind = 4 ) :: iter   ! counter for nonlinear iteration
   integer ( kind = 4 ) ::  i,  id, ip, ipp, iq, iqq, iquad, it, iuk, iukk, iuse,j , nrhs
   real ( kind = 8 )  :: par_sim    ! parameter for simple iteration ( = 0 => simple iteration, =0 =>
   real ( kind = 8 ) ::  aij,  ar, x, y,  yy
   real ( kind = 8 ) ::  test, testx, testy, trial, trialx, trialy, test_p, trial_p
   real ( kind = 8 ) ::  diff, normf
   real ( kind = 8 ) ::  rhs
   real ( kind = 8 ), dimension(2) :: vel, velx, vely    ! value of velocity at previous iteration at
   real ( kind = 8 ), dimension(2) :: mag, magx, magy    ! value of magnetic field at previous iterat
   real ( kind = 8 ), dimension(2) :: vel_prev    ! value of velocity at previous time step at quadra
   real ( kind = 8 ), dimension(2) :: mag_prev    ! value of magnetic field at previous time step at
!
! initialize
!
  id = 1   !  tells eval we are evaluating state equations

100    continue

  f_new(1:n_unknowns) =0.0
  par_sim = 0.0
!
!  LOOP over number of Newton iterations
!  first perform n_simple iterations of Simple iteration, then switch to Newton

  do iter = 1, max_newton

    if ( n_simple < iter ) then
      par_sim = 1.0
    end if

    a(1:n_unknowns,1:n_band ) = 0.0

    write (*,*) 'Iteration number ', iter
!
!  matrix assembly triangle by triangle
!
    do it = 1, n_elements

      ar = area(it) / 3.0D+00

      do iquad = 1, n_quad

        x = xq(it,iquad)
        y = yq(it,iquad)
!
!  Evaluate velocities and magnetic field at previous time step
!  at quadrature point
!
        call eval ( f_prev, id,   index, it, max_elements, max_nodes,  &
          n_local, node, x, y, xc, yc, ylength, y_inflow, y_outflow, &
          vel_prev, velx, vely, mag_prev, magx, magy )
!
!  Evaluate velocities and magnetic field at previous Newton at quadrature point
!
        call eval (f_old, id,  index, it, max_elements, max_nodes,  &
          n_local, node, x, y, xc, yc, ylength, y_inflow, y_outflow, &
          vel, velx, vely, mag, magx, magy )

        do iq = 1, n_local

          ip = node(it,iq)
          call basis_quad (it, iq, max_elements, node, x, y, xc, yc, &
            test, testx, testy)

          if ( iq <= 3) then
            test_p = basis_linear( it, iq, 1, max_elements, node, x, y, xc, yc )
          end if

          do iuk = 1, nuk

            i =  index(ip,iuk)

            if ( 0 < i .and. i /= n_unknowns) then

              if ( iuk <= 4 ) then
                call eval_rhs ( rhs_cur(1:max_nodes,iuk),  it,  &
                  max_elements, n_local, node, x, y, xc, yc,  rhs )
              end if

              if ( iuk == 1) then
                f_new(i) = f_new(i) + ar * (   rhs   &
                  +  par_sim* test * (   vel(1)*velx(1)+vel(2)*vely( 1)  )   &
                  + vel_prev(1)*test   / dt   &
                  + par_sim * coupling_par* test *  mag(2) * (  magx(2)-magy(1) )     )
              end if

              if ( iuk == 2) then
                f_new(i) = f_new(i) + ar* (    rhs   &
                  +  par_sim* test * (  vel(1)*velx(2)+vel(2)*vely(2)  )    &
                  + vel_prev(2) * test  / dt   &
                  - par_sim * coupling_par*  mag(1) *test * (  magx(2)-magy(1)  )   )
              end if

              if ( iuk == 3) then
                f_new(i) = f_new(i) +  ar * (    rhs   &
                  +  mag_prev(1)*test    / dt       &
                  - par_sim *   test * ( vely(1)*mag(2) - vely(2)*mag(1)  )  &
                  - par_sim *   test * ( vel(1)*magy(2) - vel(2)*magy(1)  )  )
              end if

              if ( iuk == 4) then
                f_new(i) = f_new(i) +  ar * (      rhs  &
                  +    mag_prev(2) *test / dt   &
                  + par_sim *   test * ( velx(1)*mag(2) - velx(2)*mag(1)  )    &
                  + par_sim *   test * ( vel(1)*magx(2) - vel(2)*magx(1)  )   )
              end if

              do iqq = 1, n_local

                ipp = node(it,iqq)
                yy = yc(ipp)
                call basis_quad ( it, iqq, max_elements, node, x, y, xc, yc, &
                  trial, trialx, trialy)

                if ( iqq <= 3) then
                  trial_p = basis_linear( it, iqq, 1, max_elements, node, &
                    x, y, xc, yc )
                end if

                do iukk = 1, nuk

                  j = index(ipp,iukk)
!
!  0 < J implies there is an unknown at this node,
!  J < 0 implies inhomogeneous Dirichlet.
!
                  if ( j /= 0 ) then

                    aij = 0.0

                    if ( iuk == 1 ) then   ! equation for u

                      if ( iukk == 1 )  then   ! unknown for u

                        aij = test * trial / dt &
                          +  visc * ( testx * trialx + testy * trialy )  &
                          + test * trial * velx(1) * par_sim       &
                          + test * ( trialx * vel(1) +  trialy * vel(2) )

                      else if (iukk == 2) then   ! unknown for v

                        aij = par_sim * test * trial * vely(1)

                      else if (iukk == 3) then   ! unknown for B1

                        aij =  - coupling_par * mag(2) *  trialy * test

                      else if (iukk == 4) then   ! unknown for B2

                        aij = coupling_par * ( test * trial &
                          * ( magx(2)-magy(1) ) *par_sim   &
                          + mag(2) * trialx * test  )

                      else    ! unknown for p

                        aij = -testx * trial_p

                      end if   ! over iukk

                    else if (iuk == 2) then   ! equation for v

                      if (iukk == 1) then

                        aij = par_sim * test * trial * velx(2)

                      else if ( iukk == 2 ) then   ! unknown for v
                        aij =  test * trial / dt &
                          + visc * ( testx * trialx + testy * trialy )   &
                          + test * trial * vely(2) * par_sim    &
                          + test * ( trialy * vel(2) +  trialx * vel(1) )

                      else if (iukk == 3) then   ! unknown for B1
                        aij =  -coupling_par * test * trial &
                          * ( magx(2)-magy(1) ) * par_sim   &
                          +  coupling_par * test * trialy * mag(1)

                      else if (iukk == 4) then   ! unknown for B2
                        aij =  -coupling_par * test * trialx * mag(1)

                      else                ! unknown for p
                        aij = -testy * trial_p
                      end if    ! over iukk

                    else if (iuk == 3) then   ! equation for B1

                      if (iukk == 1) then
                        aij =  - test * trialy * mag(2) &
                          - test * trial * magy(2) ! unknown for u

                      else if (iukk == 2) then   ! unknown for v
                        aij = test * trialy * mag(1) &
                          +   test * trial * magy(1)

                      else if (iukk == 3) then   ! unknown for B1
                        aij =   test * trial / dt  &
                          +  (1.0/rey_mag) * testy * trialy   &
                          +   test * trial * vely(2) * par_sim  &
                          +  test * trialy * vel(2) * par_sim   &
                          +(1.0/rey_mag) * testx * trialx    !  div div term

                      else if (iukk == 4) then   ! unknown for B2
                        aij =  (1.0/rey_mag) * (-testy) *trialx  &
                          -   test * trial * vely(1) * par_sim  &
                          -   test * trialy * vel(1) *par_sim  &
                                    +(1.0/rey_mag) * testx * trialy    !  div div term


                       else                ! unknown for p
                              aij = 0.0
                       end if    ! over iukk

                     else if (iuk == 4) then   ! equation for B2

                       if (iukk == 1) then
                         aij =    test * trial * magx(2)&  ! unknown for u
                           +   test * trialx * mag(2)

                       else if (iukk == 2) then   ! unknown for v
                         aij =  -  test * trialx * mag(1) &
                           -  test * trial * magx(1)

                       else if (iukk == 3) then   ! unknown for B1
                         aij =  (1./ rey_mag) * testx * (-trialy)    &
                           -   test * trial * velx(2) *par_sim &
                           -   test * trialx * vel(2)  *par_sim &
                           +(1./rey_mag) * testy * trialx    !  div div term

                       else if (iukk == 4) then   ! unknown for B2
                         aij =   test * trial / dt &
                           + (1./ rey_mag) * testx * (trialx) &
                           +   test * trial * velx(1) *par_sim   &
                           +   test * trialx * vel(1) *par_sim    &
                           +(1./rey_mag) * testy * trialy    !  div div term

                       else                ! unknown for p
                              aij =  0.0
                       end if    ! over iukk

                     else                      ! equation for pressure
                       if (iukk == 1  ) aij = -trialx * test_p
                       if (iukk == 2  ) aij = - trialy * test_p

                     end if   !  end if over iuk=1,2,3,4,5

                     if ( j > 0 ) then
                       iuse = j - i + n_lband + 1
                       a(i, iuse) = a(i, iuse) + aij*ar    !  add entries into matrix
!
!  Add terms to rhs side due to inhomogeneous bc
!
                     else if  ( j < 0 ) then

                       f_new(i) = f_new(i) &
                         - ar * ubdry ( j, yy, ylength, y_inflow, y_outflow )  &
                         * aij

                     end if

                   end if

                 end do    ! end of loop over unknowns for trial
               end do      ! end of loop over local nodes for trial

            end if          ! end if for i > 0 conditional

          end do        ! end of loop over unknowns for test
        end do          ! end of loop over local nodes for test
      end do            ! end of loop over quadrature points
    end do              ! end of loop over triangles
!
!
!  To avoid singularity of the pressure system, the last pressure
!  is simply assigned a value of 1.
!
    f_new( n_unknowns) = 1.0
    a (  n_unknowns, n_band ) = 1.
!
!  Factor the matrix and solve the linear system.
!

       call banded(  a, f_new, work, ipivot, n_unknowns, n_lband, &
         n_lband, n_unknowns )

!
!  Post process the pressure by computing the mean and subtracting
!  from nodal values of p
!
    call post_process ( area, flag_write, f_new, index, max_elements,    &
             max_nodes,   n_elements, n_nodes, node )

!
!  Check for convergence by looking at normalized difference in
!  successive iterates
!
    diff = 0.
    normf = 0.0

    do i = 1, n_unknowns
      diff = diff + ( f_old(i) - f_new(i) )**2
      normf = normf + ( f_new(i) )**2
    end do

    normf = sqrt(normf)
    diff = sqrt(diff) / normf
    write (*,*) ' normalized difference in iterates  is ', diff

    if (diff < tol_newton) then
        write (*,*) 'Newton iteration for state has converged in  ', iter, ' iterations'
        write(*,*)
        write (7,*) 'Newton iteration for state has converged in  ', iter, ' iterations'
        write(7,*)
        return

      else
        f_old(1:n_unknowns) = f_new(1:n_unknowns)
        f_new(1:n_unknowns) = 0.0

    end if

  end do    ! end loop for Newton iteration

  time_cur = time_cur - dt
  dt = dt / 2.0
  write (*,*) 'Newton iteration did not converge so decreasing delta t to  ', dt

  if ( dt < 0.00005 ) then
    write(*,*) 'reduced delta t too much so stop '
    stop
  end if

  time_cur = time_cur + dt
  f_old(1:n_unknowns) = f_prev(1:n_unknowns)
  go to 100

  return
end
subroutine output ( flag_write, index, mag_file_name, max_nodes,  &
  n_nodes, soln, uv_file_name, xc, yc )

!*****************************************************************************80
!
!! OUTPUT outputs x, y, and velocities at each point
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none

  integer ( kind = 4 )  :: max_nodes

  integer ( kind = 4 )  :: flag_write
  integer ( kind = 4 ), dimension(max_nodes,*)  ::  index
  integer ( kind = 4 )  :: n_nodes
  real ( kind = 8 ), dimension(max_nodes,*) ::   soln
  real ( kind = 8 ), dimension(*)  ::  xc,yc
  character ( len = * ) uv_file_name
  character ( len = * ) mag_file_name
  integer ( kind = 4 ) mag_unit
  integer ( kind = 4 ) :: ip
  integer ( kind = 4 ) uv_unit
  real ( kind = 8 ) ::  x, y

  call get_unit ( uv_unit )
  open ( unit = uv_unit, file = uv_file_name )

  call get_unit ( mag_unit )
  open ( unit = mag_unit, file = mag_file_name )

    do ip = 1, n_nodes

      y = yc(ip)
      x = xc(ip)

      write ( uv_unit, '(4(2x,g14.6))' ) x, y, soln(ip,1), soln(ip,2)
      write ( mag_unit, '(4(2x,g14.6))' ) x, y, soln(ip,3), soln(ip,4)

    end do

  close ( unit = uv_unit )
  close ( unit = mag_unit )

  return
end
subroutine post_process ( area, flag_write,  f,  index, max_elements,    &
  max_nodes,  n_elements, n_nodes, node )

!*****************************************************************************80
!
!! POST_PROCESS resets the pressure, by computing the mean and subtracting
!  that from all the pressures.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none
!
! INPUT
  real ( kind = 8 ), dimension(*) :: area
  real ( kind = 8 ), dimension(*) :: f
  integer ( kind = 4 )  :: flag_write
  integer ( kind = 4 ) :: max_elements, max_nodes
  integer ( kind = 4 ), dimension (max_nodes,*)  ::  index      ! unknown numbering array (ordering: u,v, p)
  integer ( kind = 4 ), dimension (max_elements,*)  ::  node
  integer ( kind = 4 ) :: n_nodes                    ! number of nodes
  integer ( kind = 4 ) :: n_elements                 ! number of triangles c
!
!  Local variables
!
  integer ( kind = 4 ) :: i, ip, iq, it
  real ( kind = 8 ) :: pmean

  pmean = 0.0

  do it = 1, n_elements

    do  iq = 1, 3
      ip = node( it,iq )
      i = index( ip,5 )
      if (i >= 0) pmean = pmean + f(i) * area(it) / 3.0D+00
    end do

  end do

   if (flag_write >= 4 ) write (*,*) 'preset - mean pressure=', pmean

  do ip = 1, n_nodes
    i = index( ip, 5 )
    if ( i > 0 ) f(i) = f(i) - pmean
  end do

  return
end
subroutine r8ud_io ( action, file_name, file_unit, record, n, x )

!*****************************************************************************80
!
!! R8UD_IO reads or writes fixed size vectors, using R8UD protocol.
!
!  Discussion:
!
!    It is assumed that the user wants to store and retrieve a number
!    of vectors.  All the vectors are of the same size, and the user
!    always specifies an index or record number when writing or retrieving
!    a particular vector.  At any time, the user can write a new vector,
!    or retrieve any vector that has been written to the file earlier.
!
!    Data type: R8 ( real ( kind = 8 ) or real ( kind = 8 ) ).
!    Format:    U  ( unformatted or binary )
!    Access:    D  ( direct access )
!
!    The first call to this routine for a given file should be with
!    the 'C' action, specifying FILE_NAME and the vector size N.
!    The output of this call is FILE_UNIT, an integer which must be
!    included on all subsequent calls that manipulate that file.
!
!    To put information into the file, use the 'W' action, specifying
!    the FILE_NAME and FILE_UNIT, as well as the RECORD number,
!    the size N of the vector, and the vector X(1:N) itself.  (Note
!    that every vector written to the file must have the same size.)
!
!    To get information from the file, use the 'R' action, specifying
!    the FILE_NAME and FILE_UNIT, as well as the RECORD number,
!    the size N of the vector.  The routine will return the given vector
!    in X(1:N).
!
!    When you are done with the file, you can close and save it with the
!    'F' action, if you specify the FILE_NAME and FILE_UNIT.
!
!    You can, instead, close and delete the file, using the 'D'
!    action and FILE_NAME and FILE_UNIT.
!
!
!    Things that can go wrong:
!
!    * You should create a file before doing anything else to it.
!
!    * You should only read a record if you have already written it.
!
!    * You should always use exactly the same set of values FILE_NAME,
!      FILE_UNIT and N, when accessing the file.
!
!    * Because I want to allow multiple files, I haven't kept track
!      of properties associated with a single file.  In particular,
!      I don't keep track of the maximum record written, which would
!      make it easy for the user to write the next record, or to
!      determine the total size of the file.
!
!    Features:
!
!    * You can read or write or rewrite any record at any time.
!
!    * You can handle two or more files simultaneously, as long as
!      you use the appropriate set of FILE_NAME, FILE_UNIT and N
!      values for each.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ACTION:
!    'c' to create the file FILE_NAME with record size N;
!    'w' to write vector X(1:N) as record RECORD;
!    'r' to read vector X(1:N) which is record RECORD;
!    'f' to close and save the file FILE_NAME;
!    'd' to close and delete the file FILE_NAME.
!    's' to print some statistics about file FILE_NAME.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input/output, integer FILE_UNIT, the unit number associated with the file.
!    For the 'c' command, this value is output.
!    For the 'w', 'r', 'f' and 'd' commands, this value is input.
!
!    Input, integer RECORD, the index of the record to be written or read.
!    (This value is needed for the 'W' and 'R' actions only.)
!
!    Input, integer N, the size of each vector.
!    (This value is needed for the 'C', 'W' and 'R' actions only.)
!
!    Input/output, real ( kind = 8 ) X(N), the vector to be written or read.
!    (This value is input for the 'W' action, output for the 'R' action,
!    and not needed otherwise.)
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: word_length = 2

  character action
  integer ( kind = 4 ), parameter :: byte_length = 4 * word_length
  integer ( kind = 4 ), save :: call_num = 0
  logical file_exist
  logical file_is_open
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) record
  integer ( kind = 4 ), save :: record_length = -1
  real ( kind = 8 ) x(n)
!
!  Create a file of given name.
!
  if ( action(1:1) == 'c' ) then

    if ( file_exist ( file_name ) ) then
      call file_delete ( file_name )
    end if

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a)' ) '  The value of N is not positive.'
      write ( *, '(a)' ) '  Could not create "' // trim ( file_name ) // '".'
      stop
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Creating file "' // trim ( file_name ) // '".'

    call get_unit ( file_unit )
!
!  NOTE:
!    On some machines, the record length of a binary direct access file
!    is in WORDS, and on some machines, it is in BYTES.
!    So on some systems (Alpha's for one), use the first statement.
!    on others (Apple's, for one), use the second.
!
    record_length = n * word_length
!   record_length = n * byte_length

    open ( file = file_name, unit = file_unit, status = 'new', &
      form = 'unformatted', access = 'direct', recl = record_length, &
      iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a,i6)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not create "' // trim ( file_name ) // '".'
      stop
    end if

    call_num = 1
!
!  Close (if necessary) and delete the file.
!
  else if ( action(1:1) == 'd' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The file does not exist.'
      write ( *, '(a)' ) '  Could not delete "' // trim ( file_name ) // '".'
      return
    end if

    if ( file_is_open ( file_name ) ) then
      close ( unit = file_unit )
    end if

    call file_delete ( file_name )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file "' // trim ( file_name ) // &
      '" has been deleted.'

    file_unit = -1
    record_length = -1
!
!  Close (if necessary) and save the file.
!
  else if ( action(1:1) == 'f' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The file does not exist.'
      write ( *, '(a)' ) '  Could not close and save "' &
        // trim ( file_name ) // '".'
      return
    end if

    if ( file_is_open ( file_name ) ) then
      close ( unit = file_unit )
    end if
!
!  Read a vector X.
!
  else if ( action(1:1) == 'r' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to READ data from a file'
      write ( *, '(a)' ) '  that does not exist.'
      write ( *, '(a)' ) '  The file name is "' // trim ( file_name ) // '".'
      stop
    end if

    if ( record < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to read an illegal record.'
      write ( *, '(a)' ) '  Could not read from "' // trim ( file_name ) // '".'
      stop
    end if

    read ( unit = file_unit, rec = record, iostat = ios ) x(1:n)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a,i6)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not read from "' // trim ( file_name ) // '".'
      stop
    end if
!
!  Print statistics.
!
  else if ( action(1:1) == 's' ) then

    call_num = call_num + 1

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  File statistics:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file name is "' // trim ( file_name ) // '".'
    write ( *, '(a,i12)' ) '  The length of each vector is N =    ', n
    write ( *, '(a,i12)' ) '  The word length of one data item =  ', word_length
    write ( *, '(a,i12)' ) '  The byte length of one data item =  ', byte_length
    write ( *, '(a,i12)' ) '  The length of each record is =      ', &
      record_length
!
!  Write a vector X.
!
  else if ( action(1:1) == 'w' ) then

    call_num = call_num + 1

    if ( .not. file_exist ( file_name ) ) then

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
        write ( *, '(a)' ) '  The value of N is not positive.'
        write ( *, '(a)' ) '  Could not write to "' &
          // trim ( file_name ) // '".'
        stop
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Creating file "' // trim ( file_name ) // '".'

      call get_unit ( file_unit )

      record_length = n * byte_length

      open ( file = file_name, unit = file_unit, status = 'new', &
        form = 'unformatted', access = 'direct', recl = record_length, &
        iostat = ios )

    end if

    if ( record < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a)' ) '  You have asked to write an illegal record.'
      write ( *, '(a)' ) '  Could not write to "' // trim ( file_name ) // '".'
      stop
    end if

    write ( unit = file_unit, rec = record, iostat = ios ) x(1:n)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UD_IO - Fatal error!'
      write ( *, '(a,i6)' ) '  IO error of type ', ios
      write ( *, '(a)' ) '  Could not write to "' // trim ( file_name ) // '".'
      stop
    end if
!
!  Unrecognized ACTION.
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8UD_IO - Warning!'
    write ( *, '(a)' ) '  The requested action "' // trim ( action ) &
      // '" is not recognized!'
    write ( *, '(a)' ) '  The file name was "' // trim ( file_name ) // '".'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Legal actions include:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "c", create a file;'
    write ( *, '(a)' ) '  "d", delete a file;'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "w", write a vector;'
    write ( *, '(a)' ) '  "r", read a vector;'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "s", print statistics;'

  end if

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine set_initial ( index, max_nodes, n_nodes, n_unknowns, xc, yc, &
  f, f_state )

!*****************************************************************************80
!
!! SET_INITIAL sets the initial condition for the state solution.
!
!  Discussion:
!
!    Set the internal parameter "OPTION" to
!
!    1, to initialize velocity to the negative of the target,
!       and magnetic field to 0;
!    2, to initialize velocity to the negative of the target,
!       and read magnetic initial condition from the file "mag_ic.txt".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Parameters:
!
!    Input, integer INDEX(MAX_NODES,5), contains the numbering of
!    unknowns at nodes.
!
!    Input, integer MAX_NODES, the maximum number of nodes.
!
!    Input, integer N_NODES, the number of nodes.
!
!    Input, integer N_UNKNOWNS, the number of unknowns.
!
!    Input, real ( kind = 8 ) XC(N_NODES), YC(N_NODES), the coordinates
!    of the nodes.
!
!    Output, real ( kind = 8 ) F(N_UNKNOWNS), the initial solution vector.
!
!    Output, real ( kind = 8 ) F_STATE(MAX_NODES,4), the initial solution,
!    rearranged by node.
!
  implicit none

  integer ( kind = 4 ) :: max_nodes
  integer ( kind = 4 ) :: n_nodes
  integer ( kind = 4 ) :: n_unknowns

  real ( kind = 8 ) dummy1
  real ( kind = 8 ) dummy2
  real ( kind = 8 ), dimension(n_unknowns) :: f
  real ( kind = 8 ), dimension(max_nodes,4) :: f_state
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension(max_nodes,5) :: index
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iuk
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: option = 2
  real ( kind = 8 ), dimension(4) :: soln
  real ( kind = 8 ) :: x
  real ( kind = 8 ), dimension(n_nodes) :: xc
  real ( kind = 8 ) :: y
  real ( kind = 8 ), dimension(n_nodes) :: yc

  f(1:n_unknowns) = 0.0
!
!  Option 1:
!    The initial velocity is set to the negative of the target,
!    the initial magnetic field is zero.
!
  if ( option == 1 ) then

    do ip = 1, n_nodes

      x = xc(ip)
      y = yc(ip)

      call target ( x, y, soln(1), soln(2), soln(3), soln(4) )

      soln(1) = -soln(1)
      soln(2) = -soln(2)
      soln(3) = 0.0
      soln(4) = 0.0

      do iuk = 1, 4
        j = index(ip,iuk)
        f(j) = soln(iuk)
      end do

    end do
!
!  Option 2:
!    The initial velocity is set to the negative of the target,
!    the initial magnetic field is read from a file.
!
  else

    call get_unit ( iunit )

    open ( unit = iunit, file = 'mag_ic.txt', status = 'old', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SET_INITIAL - Fatal error!'
      write ( *, '(a)' ) '  Initialization option 2 was set.'
      write ( *, '(a)' ) '  This requires access to the file "mag_ic.txt".'
      write ( *, '(a)' ) '  However, this file could not be opened.'
      stop
    end if

    do ip = 1, n_nodes

      x = xc(ip)
      y = yc(ip)

      call target ( x, y, soln(1), soln(2), dummy1, dummy2 )
      soln(1) = -soln(1)
      soln(2) = -soln(2)

      read ( iunit, * ) dummy1, dummy2, soln(3), soln(4)
      soln(3) = 5.0 * soln(3)
      soln(4) = 5.0 * soln(4)

      do iuk = 1, 4
        j = index(ip,iuk)
        f(j) = soln(iuk)
      end do

    end do

    close ( unit = iunit )

  end if
!
!  Make a node-based copy of the data.
!
  do ip = 1, n_nodes

    do iuk = 1, 4
      j = index(ip,iuk)
      f_state(ip,iuk) = f(j)
    end do

  end do

  return
end
subroutine target ( x, y, target_u, target_v, target_b1, target_b2  )

!*****************************************************************************80
!
!! TARGET sets value of target solution when it is given explicitly by a function
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) ::  x, y
  real ( kind = 8 ) :: target_b1, target_b2, target_u, target_v
  real ( kind = 8 ) :: argx, argy, dphidx, dphidy, phix, phiy

  argx = 2.0 * pi * 0.4D+00 * x
  argy = 2.0 * pi * 0.4D+00 * y

  phix =  ( 1.0 - x )**2 * (  1.0 - cos ( argx ) )

  phiy =  ( 1.0 - y )**2 * (  1.0 - cos ( argy ) )

  dphidx = -2.0 * ( 1.0 - x ) * ( 1.0 - cos ( argx ) ) &
                + ( 1.0 - x )**2  * ( 2.0 * pi * 0.4D+00 * sin ( argx )  )

  dphidy = -2.0 * ( 1.0 - y ) * ( 1.0 - cos ( argy ) ) &
                + ( 1.0 - y )**2  * ( 2.0 * pi * 0.4D+00 * sin ( argy )  )

  target_u =   10.0 * phix * dphidy
  target_v = - 10.0 * phiy * dphidx

  target_b1 = sin ( pi * x ) * cos ( pi * y )
  target_b2 = cos ( pi * x ) * sin ( pi * y )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
function ubdry ( id, y, ylength, y_inflow, y_outflow  )

!*****************************************************************************80
!
!! UBDRY evaluates the velocity inlet profile.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none

  real ( kind = 8 ) :: a
  real ( kind = 8 ) :: b
  real ( kind = 8 ) :: c
  integer ( kind = 4 ) :: id
  real ( kind = 8 ) :: ubdry
  real ( kind = 8 ) :: y
  real ( kind = 8 ) :: y_inflow
  real ( kind = 8 ) :: y_outflow
  real ( kind = 8 ) :: ylength

  if ( id == -1 ) then
    a = -4.0 / ( y_inflow * y_inflow )
    b = 4.0 / y_inflow
    c = 0.0
  else
    a = -4.0 / (  (ylength - y_outflow) ** 2 )
    b =  4.0 * ( ylength + y_outflow) / ( (ylength - y_outflow)**2 )
    c = -4.0 * ylength * y_outflow / ( (ylength - y_outflow)**2 )
  end if

  ubdry = a * y * y + b * y + c

  return
end
subroutine xy_write ( xy_file_name, node_num, xc, yc )

!*****************************************************************************80
!
!! XY_WRITE writes the coordinate data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Parameters:
!
!    Input, character ( len = * ) XY_FILE_NAME, the name of the
!    output file containing the node coordinate values.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NODE_NUM), YC(NODE_NUM), the coordinates
!    of nodes.
!
  implicit none

  integer ( kind = 4 ) :: node_num

  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) node
  real ( kind = 8 ), dimension(node_num) :: xc
  character ( len = * ) xy_file_name
  real ( kind = 8 ), dimension(node_num) :: yc

  call get_unit ( iunit )

  open ( unit = iunit, file = xy_file_name, status = 'replace' )

  do node = 1, node_num

    write ( iunit, '(2x,f12.6,2x,f12.6)' ) xc(node), yc(node)

  end do

  close ( unit = iunit )

  return
end
