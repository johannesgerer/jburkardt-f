program main

!*******************************************************************************
!
!! MAIN is the main program for NAST2D_F90.
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM 1998.
!
!  Translated from C into FORTRAN90 by:
!
!     Dr. C. David Pruett
!     Dept. of Mathematics & Statistics
!     MSC 7803
!     James Madison University
!     Harrisonburg, VA 22807 USA
!     dpruett@math.jmu.edu
!     www.math.jmu.edu/~dpruett
!     540-568-6227
!
  use nrtype
 
  implicit none

  real ( rp ) :: beta
  integer :: cycles
  real ( rp ) :: del_inj
  real ( rp ) :: del_streak
  real ( rp ) :: del_trace
  real ( rp ) :: del_vec
  real ( rp ) :: delt
  real ( rp ) :: delx
  real ( rp ) :: dely
  real ( rp ) :: eps
  real ( rp ), dimension(:,:), allocatable :: f
  integer ( i2b ), dimension(:,:), allocatable :: flag
  real ( rp ), dimension(:,:), allocatable :: g
  real ( rp ) :: gamma
  real ( rp ) :: gx
  real ( rp ) :: gy
  real ( rp ), dimension(:,:), allocatable :: heat
  integer :: i
  integer :: ibound
  integer :: ifull
  integer :: imax
  character ( len = 30 ) :: infile
  integer :: init_case
  character ( len = 30 ) :: inputfile
  integer :: isurf
  integer :: itermax
  integer :: itersor
  integer :: iwrite
  integer :: jmax
  integer :: n
  real ( rp ) :: omega
  character ( len = 30 ) :: outfile
  real ( rp ), dimension(:,:), allocatable :: p
  integer :: p_bound
  type ( particleline ), dimension(:), allocatable :: partlines
  real ( rp ) :: pos1x
  real ( rp ) :: pos1y
  real ( rp ) :: pos2x
  real ( rp ) :: pos2y
  integer :: ppc
  real ( rp ) :: pr
  character ( len = 30 ) :: problem
  real ( rp ), dimension(:,:), allocatable :: psi
  real ( rp ) :: re
  real ( rp ) :: res
  real ( rp ), dimension(:,:), allocatable :: rhs
  character ( len = 30 ) :: streakfile
  real ( rp ) :: t
  real ( rp ) :: t_end
  real ( rp ) :: tau
  real ( rp ), dimension(:,:), allocatable :: temp
  real ( rp ) :: ti
  character ( len = 30 ) :: tracefile
  real ( rp ), dimension(:,:), allocatable :: u
  real ( rp ) :: ui
  real ( rp ), dimension(:,:), allocatable :: v
  character ( len = 30 ) :: vecfile
  real ( rp ) :: vi
  integer :: we
  integer :: wn
  real ( rp ), dimension(:,:), allocatable :: work
  real ( rp ), dimension(:,:), allocatable :: work1
  integer :: ws
  integer :: ww
  real ( rp ) :: xlength
  real ( rp ) :: ylength
  real ( rp ), dimension(:,:), allocatable :: zeta

  include 'interfaces.h'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NAST2D_F90'
  write ( *, '(a)' ) '  A FORTRAN90 version of NAST2D'
  write ( *, '(a)' ) '  A program to solve problems involving'
  write ( *, '(a)' ) '  the 2D time-dependent Navier Stokes equations'
  write ( *, '(a)' ) '  with temperature.'
  write ( *, '(a)' ) ' '
!
!  Initialize.
!
  call read_parameters ( problem, inputfile, infile, outfile, xlength, &
    ylength, imax, jmax, t_end, delt, tau, delx, dely, del_vec, &
    del_trace, del_streak, del_inj, vecfile, tracefile, streakfile, &
    n, pos1x, pos1y, pos2x, pos2y, itermax, eps, omega, gamma, p_bound, &
    re, pr, beta, gx, gy, ui, vi, ti, ww, we, wn, ws )

  itersor = 0
  ifull = 0
  isurf = 0
  ibound = 0
!
!  Allocate arrays.
!
  allocate ( f(0:imax+1,0:jmax+1) )
  allocate ( flag(0:imax+1,0:jmax+1) ) 
  allocate ( g(0:imax+1,0:jmax+1) )
  allocate ( heat(0:imax,0:jmax) )
  allocate ( p(0:imax+1,0:jmax+1) )
  allocate ( psi(0:imax,0:jmax) )
  allocate ( rhs(0:imax+1,0:jmax+1) )
  allocate ( temp(0:imax+1,0:jmax+1) )
  allocate ( u(0:imax+1,0:jmax+1) )
  allocate ( v(0:imax+1,0:jmax+1) )
  allocate ( work(0:imax+1,0:jmax+1) )
  allocate ( work1(1:imax,1:jmax) )
  allocate ( zeta(1:imax-1,1:jmax-1) )
!
!  Read initial values from file "infile" if it exists.
!
  call read_bin ( infile, u, v, p, temp, flag, imax, jmax, init_case )

  if ( init_case == 1 ) then

    call init_uvp ( problem, u, v, p, temp, imax, jmax, ui, vi, ti )

    call init_flag ( problem, flag, imax, jmax, delx, dely, ibound )

  else if ( init_case == 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NAST2D_F90 - Fatal error!'
    write ( *, '(a)' ) '  Read error in READ_BIN.'
    stop

  end if
!
!  Initialize particles for streaklines or particle tracing.
!  Make an additional copy to hold injection points.
!
  if ( streakfile /= "none" .or. &
       tracefile /= "none" ) then

    allocate ( partlines(1:n) )

    do i = 1, n
      partlines(i)%length = 0
      nullify ( partlines(i)%particles%next )
    end do

    call set_particles ( n, pos1x, pos1y, pos2x, pos2y, partlines )

    call inject_particles ( n, partlines )

  end if
!
!  Initialize particles for free-boundary problems.
!
  if ( & 
    problem == "drop" .or. &
    problem == "dam" ) then

    if ( allocated ( partlines ) ) then
      deallocate ( partlines )
    end if

    ppc = 4

    if ( problem == "drop" ) then
      n = 2
    else if ( problem == "dam" ) then
      n = 1
    end if

    allocate ( partlines(1:n) )

    do i = 1, n
      partlines(i)%length = 0
      nullify ( partlines(i)%particles%next )
    end do

    call init_particles ( n, imax, jmax, delx, dely, ppc, problem, &
      u, v, partlines )

    call inject_particles ( n, partlines )

  end if
!
!  Initialize the boundary conditions.
!
  call setbcond ( u, v, p, temp, flag, imax, jmax, ww, we, wn, ws )

  call setspecbcond ( problem, u, v, temp, imax, jmax, ui, vi )

  t = 0.0
  cycles = 0
!
!  Time loop.
!
  do while ( t < t_end )
     
    call comp_delt ( delt, t, imax, jmax, delx, dely, u, v, re, pr, &
      tau, iwrite, del_trace, del_inj, del_streak, del_vec )
!
!  Determine fluid cells for free-boundary problems and set
!  boundary conditions at free surface.
!
    if ( &
      problem == "drop" .or. &
      problem == "dam" ) then

      call mark_cells ( flag, imax, jmax, delx, dely, ifull, isurf, &
        n, partlines )

      call set_uvp_surface ( u, v, p, flag, gx, gy, imax, jmax, &
        re, delx, dely, delt )

    else

      ifull = imax * jmax - ibound

    end if
!
!  Compute the new temperature.
!
    call comp_temp ( u, v, temp, flag, work, imax, jmax, delt, &
      delx, dely, gamma, re, pr )
!
!  Compute tentative velocity field (F,G).
!
    call comp_fg ( u, v, temp, f, g, flag, imax, jmax, & 
      delt, delx, dely, gx, gy, gamma, re, beta )
!
!  Compute right-hand side for pressure equation.
!
    call comp_rhs ( f, g, rhs, flag, imax, jmax, delt, delx, dely )
!
!  Solve the pressure equation by successive over relaxation (SOR).
!
    if ( 0 < ifull ) then

      call poisson ( p, rhs, flag, imax, jmax, delx, dely, &
        eps, itersor, itermax, omega, res, ifull, p_bound )

      write (6, '(a,g10.4,a,g10.4,a,i4,a,g12.6,a,i4,i4,i4)' ) &
        ' t=', t+delt, ' delt=', delt, &
        ' its=', itersor, ' res=', res, ' F,S,B-cells=', ifull, &
        isurf, ibound

    end if
!
!  Compute the new velocity field.
!
    call adap_uv ( u, v, f, g, p, flag, imax, jmax, delt, delx, dely )
!
!  Set the boundary conditions.
!
    call setbcond ( u, v, p, temp, flag, imax, jmax, ww, we, wn, ws )

    call setspecbcond ( problem, u, v, temp, imax, jmax, ui, vi )

    if ( problem == "drop"    .or. problem == "dam"  .or. & 
         problem == "molding" .or. problem == "wave" ) then

      call set_uvp_surface ( u, v, p, flag, gx, gy, imax, jmax, &
        re, delx, dely, delt )

    end if

    if ( iand(iwrite,8 ) /= 0 .and. vecfile  /= "none" ) then

      call comp_psi_zeta ( u, v, psi, zeta, flag, imax, jmax, delx, dely )

      call comp_heat ( u, v, temp, heat, flag, re, pr, &
        imax, jmax, delx, dely )

      call outputvec_bin ( u, v, p, temp, psi, zeta, heat, flag, &
        work1, xlength, ylength, imax, jmax, vecfile )

    end if

    if ( iand ( iwrite, 8 ) /= 0 .and. outfile /= "none" ) then
      call write_bin ( outfile, u, v, p, temp, flag, imax, jmax )
    end if

    if ( tracefile /= "none" ) then
      call particle_tracing ( tracefile, t, imax, jmax, delx, dely, &
        delt, u, v, flag, n, partlines, iwrite )
    end if

    if ( streakfile /= "none" ) then
      call streaklines ( streakfile, iwrite, imax, jmax, delx, dely, &
        delt, t, u, v, flag, n, partlines )
    end if
!
!  Advance the time
!
    t = t + delt
    cycles = cycles + 1

  end do

  call comp_psi_zeta ( u, v, psi, zeta, flag, imax, jmax, delx, dely )

  call comp_heat ( u, v, temp, heat, flag, re, pr, &
    imax, jmax, delx, dely )

  if ( vecfile /= "none" ) then

    call outputvec_bin ( u, v, p, temp, psi, zeta, heat, flag, &
      work1, xlength, ylength, imax, jmax, vecfile )

  end if

  if ( outfile /= "none" ) then
    call write_bin ( outfile, u, v, p, temp, flag, imax, jmax )
  end if
!
!  Output the stream function for visualization.
!
  call output_one_real ( psi, 0, imax, 0, jmax, 'psi_data.txt' )
!
!  Free memory.
!
  deallocate ( f )
  deallocate ( flag )
  deallocate ( g )
  deallocate ( heat )
  deallocate ( p )
  deallocate ( psi )
  deallocate ( rhs )
  deallocate ( temp )
  deallocate ( u )
  deallocate ( v )
  deallocate ( work )
  deallocate ( zeta )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NAST2D:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine timestamp ( )

!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    15 March 2003
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
!
  character ( len = 40 ) string
!
  call timestring ( string )

  write ( *, '(a)' ) trim ( string )

  return
end
subroutine timestring ( string )

!*******************************************************************************
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
