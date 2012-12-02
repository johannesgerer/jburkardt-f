!
!  boundary.f90 interfaces
!
interface
  subroutine setbcond ( u, v, p, temp, flag, imax, jmax, ww, we, wn, ws )

  use nrtype

  implicit none

  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( inout ) :: p
  real ( rp ), dimension(0:,0:), intent ( inout ) :: temp
  real ( rp ), dimension(0:,0:), intent ( inout ) :: u
  real ( rp ), dimension(0:,0:), intent ( inout ) :: v
  integer, intent ( in ) :: we
  integer, intent ( in ) :: wn
  integer, intent ( in ) :: ws
  integer, intent ( in ) :: ww

  end subroutine setbcond
end interface

interface
  subroutine setspecbcond ( problem, u, v, temp, imax, jmax, ui, vi )

  use nrtype

  implicit none

  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  character ( len = * ), intent ( in ) :: problem
  real ( rp ), dimension(0:,0:), intent(inout) :: temp
  real ( rp ), dimension(0:,0:), intent(inout) :: u
  real ( rp ), intent ( in ) :: ui
  real ( rp ), dimension(0:,0:), intent(inout) :: v
  real ( rp ), intent ( in ) :: vi

  end subroutine setspecbcond
end interface
!
!  interfaces for extras.f90
!
interface
  subroutine output_one_real ( field, imin, imax, jmin, jmax, filename )

  use nrtype

  implicit none

  real ( rp ), dimension(0:,0:), intent ( in ) :: field
  character ( len = * ), intent ( in ) :: filename
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: imin
  integer, intent ( in ) :: jmax
  integer, intent ( in ) :: jmin

  end subroutine output_one_real
end interface

interface
  subroutine output_one_integer ( field, imin, imax, jmin, jmax, filename )

  use nrtype

  implicit none

  integer, dimension(0:,0:), intent ( in ) :: field
  character ( len = * ), intent ( in ) :: filename
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: imin
  integer, intent ( in ) :: jmax
  integer, intent ( in ) :: jmin

  end subroutine output_one_integer
end interface

interface
  subroutine output_three_real ( u, v, p, imin, imax, jmin, jmax, filename )

  use nrtype

  implicit none

  character ( len = * ), intent ( in ) :: filename
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: imin
  integer, intent ( in ) :: jmax
  integer, intent ( in ) :: jmin
  real ( rp ), dimension (0:,0:), intent ( in ) :: p
  real ( rp ), dimension (0:,0:), intent ( in ) :: u
  real ( rp ), dimension (0:,0:), intent ( in ) :: v

  end subroutine output_three_real
end interface
!
! init.f90 interfaces
!
interface
  subroutine read_parameters ( problem, inputfile, infile, outfile, xlength, &
    ylength, imax, jmax, t_end, delt, tau, delx, dely, del_vec, del_trace, &
    del_streak, del_inj, vecfile, tracefile, streakfile, n, pos1x, pos1y, &
    pos2x, pos2y, itermax, eps, omega, gamma, p_bound, re, pr, beta, gx, &
    gy, ui, vi, ti, ww, we, wn, ws )

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

  end subroutine read_parameters
end interface

interface
  subroutine init_uvp ( problem, u, v, p, temp, imax, jmax, ui, vi, ti )

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

  end subroutine init_uvp
end interface

interface
  subroutine read_bin ( filename, u, v, p, temp, flag, imax, jmax, check )

  use nrtype

  implicit none

  integer, intent ( out ) :: check
  character ( len = * ), intent ( in ) :: filename
  integer ( i2b ), dimension(0:,0:), intent ( out ) :: flag
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( out ) :: p
  real ( rp ), dimension(0:,0:), intent ( out ) :: temp
  real ( rp ), dimension(0:,0:), intent ( out ) :: u
  real ( rp ), dimension(0:,0:), intent ( out ) :: v

  end subroutine read_bin
end interface

interface
  subroutine write_bin ( filename, u, v, p, temp, flag, imax, jmax )

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

  end subroutine write_bin
end interface

interface
  subroutine init_flag ( problem, flag, imax, jmax, delx, dely, ibound )

  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( inout ) :: flag
  integer, intent ( out ) :: ibound
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  character ( len = * ), intent ( in ) :: problem

  end subroutine init_flag
end interface
!
!  Interfaces for surface.f90
!
interface
  subroutine init_particles ( n, imax, jmax, delx, dely, ppc, problem, &
    u, v, partlines )

  use nrtype

  implicit none

  integer, intent ( inout ) :: n

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  type ( particleline ), dimension(n), intent ( inout ) :: partlines
  integer, intent ( in ) :: ppc
  character ( len = * ), intent ( in ) :: problem
  real ( rp ), dimension(0:,0:), intent ( inout ) :: u
  real ( rp ), dimension(0:,0:), intent ( inout ) :: v

  end subroutine init_particles
end interface

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

interface
  subroutine mark_cells ( flag, imax, jmax, delx, dely, ifull, isurf, &
    n, partlines )

  use nrtype

  implicit none

  integer, intent ( in ) :: n

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( inout ) :: flag
  integer, intent ( out ) :: ifull
  integer, intent ( in ) :: imax
  integer, intent ( out ) :: isurf
  integer, intent ( in ) :: jmax
  type ( particleline ), dimension(n), intent ( inout ), target :: partlines

  end subroutine mark_cells
end interface

interface
  subroutine set_uvp_surface ( u, v, p, flag, gx, gy, imax, jmax, &
    re, delx, dely, delt )

  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), intent ( in ) :: gx
  real ( rp ), intent ( in ) :: gy
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( inout ) :: p
  real ( rp ), intent ( in ) :: re
  real ( rp ), dimension(0:,0:), intent ( inout ) :: u
  real ( rp ), dimension(0:,0:), intent ( inout ) :: v

  end subroutine set_uvp_surface
end interface

!
!  uvp.f90 interfaces
!   
interface

subroutine comp_temp ( u, v, temp, flag, work, imax, jmax, delt, &
  delx, dely, gamma, re, pr )

  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), intent ( in ) :: gamma
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), intent ( in ) :: pr
  real ( rp ), intent ( in ) :: re
  real ( rp ), dimension(0:,0:), intent ( inout ) :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v
  real ( rp ), dimension(0:,0:), intent ( out ) :: work

  end subroutine comp_temp
end interface

interface
  subroutine comp_fg ( u, v, temp, f, g, flag, imax, jmax, delt, delx, dely, &
    gx, gy, gamma, re, beta )

  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: beta
  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ), dimension(0:,0:), intent ( out ) :: f
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( out ) :: g
  real ( rp ), intent ( in ) :: gamma
  real ( rp ), intent ( in ) :: gx
  real ( rp ), intent ( in ) :: gy
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), intent ( in ) :: re
  real ( rp ), dimension(0:,0:), intent ( in ) :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v

  end subroutine comp_fg
end interface

interface
  subroutine comp_rhs ( f, g, rhs, flag, imax, jmax, delt, delx, dely )

  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ), dimension(0:,0:), intent ( in ) :: f
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( in ) :: g
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( out ) :: rhs

  end subroutine comp_rhs
end interface

interface
  subroutine poisson ( p, rhs, flag, imax, jmax, delx, dely, &
    eps, iter, itermax, omega, res, ifull, p_bound )

  use nrtype
     
  implicit none

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ), intent ( in ) :: eps
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer, intent ( in ) :: ifull
  integer, intent ( in ) :: imax
  integer, intent ( inout ) :: iter
  integer, intent ( in ) :: itermax
  integer, intent ( in ) :: jmax
  real ( rp ), intent ( in ) :: omega
  real ( rp ), dimension(0:,0:), intent ( inout ) :: p
  integer, intent ( in ) :: p_bound
  real ( rp ), intent ( out ) :: res
  real ( rp ), dimension(0:,0:), intent ( in ) :: rhs

  end subroutine poisson
end interface

interface
  subroutine adap_uv ( u, v, f, g, p, flag, imax, jmax, delt, delx, dely )

  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  real ( rp ), dimension(0:,0:), intent ( in ) :: f
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( in ) :: g
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( in ) :: p
  real ( rp ), dimension(0:,0:), intent ( out ) :: u
  real ( rp ), dimension(0:,0:), intent ( out ) :: v

  end subroutine adap_uv
end interface

interface
  subroutine comp_delt ( delt, t, imax, jmax, delx, dely, u, v, re, pr, &
    tau, iwrite, del_trace, del_inj, del_streak, del_vec )

  use nrtype

  implicit none
 
  real ( rp ), intent ( in ) :: del_inj
  real ( rp ), intent ( in ) :: del_streak
  real ( rp ), intent ( in ) :: del_trace
  real ( rp ), intent ( in ) :: del_vec
  real ( rp ), intent ( out ) :: delt
  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer, intent ( in ) :: imax
  integer, intent ( out ) :: iwrite
  integer, intent ( in ) :: jmax
  real ( rp ), intent ( in ) :: pr
  real ( rp ), intent ( in ) :: re
  real ( rp ), intent ( in ) :: t
  real ( rp ), intent ( in ) :: tau
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v

  end subroutine comp_delt
end interface
!
! interfaces for visual.f90
!
interface
  subroutine outputvec_bin ( u, v, p, temp, psi, zeta, heat, flag, work1, &
    xlength, ylength, imax, jmax, vecfile )

  use nrtype

  implicit none

  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( in ) :: heat
  integer, intent ( in ) :: imax
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

  end subroutine outputvec_bin
end interface

interface
  subroutine comp_psi_zeta ( u, v, psi, zeta, flag, imax, jmax, delx, dely )

  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), dimension(0:,0:), intent ( out ) :: psi
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v
  real ( rp ), dimension(1:,1:), intent ( out ) :: zeta

  end subroutine comp_psi_zeta
end interface

interface
  subroutine comp_heat ( u, v, temp, heat, flag, re, pr, imax, jmax, delx, &
    dely )

  use nrtype

  implicit none

  real ( rp ), intent ( in ) :: delx
  real ( rp ), intent ( in ) :: dely
  integer ( i2b ), dimension(0:,0:), intent ( in ) :: flag
  real ( rp ), dimension(0:,0:), intent ( out ) :: heat
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: jmax
  real ( rp ), intent ( in ) :: pr
  real ( rp ), intent ( in ) :: re
  real ( rp ), dimension(0:,0:), intent ( in ) :: temp
  real ( rp ), dimension(0:,0:), intent ( in ) :: u
  real ( rp ), dimension(0:,0:), intent ( in ) :: v

  end subroutine comp_heat
end interface

interface
  subroutine set_particles ( n, pos1x, pos1y, pos2x, pos2y, partlines )

  use nrtype

  implicit none

  integer, intent ( in ) :: n

  type ( particleline ), dimension(n), intent ( out ) :: partlines
  real ( rp ), intent ( in ) :: pos1x
  real ( rp ), intent ( in ) :: pos1y
  real ( rp ), intent ( in ) :: pos2x
  real ( rp ), intent ( in ) :: pos2y

  end subroutine set_particles
end interface

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

interface
  subroutine particle_tracing ( tracefile, t, imax, jmax, delx, dely, delt, &
    u, v, flag, n, partlines, iwrite )

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

  end subroutine particle_tracing
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
  subroutine streaklines ( streakfile, iwrite, imax, jmax, delx, dely, delt, &
    t, u, v, flag, n, partlines )

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

  end subroutine streaklines
end interface
