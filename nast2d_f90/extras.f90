subroutine output_one_real ( field, imin, imax, jmin, jmax, filename )

!*******************************************************************************
!
!! OUTPUT_ONE_REAL outputs a single real field.
!
!  Modified:
!
!    14 February 2004
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM, 1998.
!
!  Parameters:
!
!    Input, real FIELD(0:IMAX+1,0:JMAX+1), the field to be written out.
!
!    Input, integer IMIN, IMAX, JMIN, JMAX, specify that 
!    FIELD(IMIN:IMAX,JMIN:JMAX) is to be written out.
!
!    Input, character ( len = * ) FILENAME, the name of the file
!    into which the data is to be written.
!
  use nrtype

  implicit none

  real ( rp ), dimension(0:,0:), intent ( in ) :: field
  character ( len = * ), intent ( in ) :: filename
  integer i
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: imin
  integer j
  integer, intent ( in ) :: jmax
  integer, intent ( in ) :: jmin

  open ( unit = 1, file = filename, status = 'replace', form = 'formatted' )

  do j = jmin, jmax
    do i = imin, imax
      write ( 1, '(1x,i5,1x,i5,1x,f12.6)' ) i, j, field(i,j)
    end do
  end do

  close ( unit = 1 )

  return
end
subroutine output_one_integer ( field, imin, imax, jmin, jmax, filename )

!*******************************************************************************
!
!! OUTPUT_ONE_INTEGER outputs a single integer field.
!
!  Modified:
!
!    14 February 2004
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM, 1998.
!
!  Parameters:
!
!    Input, real FIELD(0:IMAX+1,0:JMAX+1), the field to be written out.
!
!    Input, integer IMIN, IMAX, JMIN, JMAX, specify that 
!    FIELD(IMIN:IMAX,JMIN:JMAX) is to be written out.
!
!    Input, character ( len = * ) FILENAME, the name of the file
!    into which the data is to be written.
!
  use nrtype

  implicit none

  integer, dimension(0:,0:), intent ( in ) :: field
  character ( len = * ), intent ( in ) :: filename
  integer i
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: imin
  integer j
  integer, intent ( in ) :: jmax
  integer, intent ( in ) :: jmin

  open ( unit = 1, file = filename, status = 'replace', form = 'formatted' )

  do j = jmin, jmax
    do i = imin, imax
      write ( 1, '(1x,i5,1x,i5,1x,i6)' ) i, j, field(i,j)
    end do
  end do

  close ( unit = 1 )

  return
end
subroutine output_three_real ( u, v, p, imin, imax, jmin, jmax, filename )

!*******************************************************************************
!
!! OUTPUT_THREE_REAL outputs the velocity and pressure fields.
!
!  Modified:
!
!    14 February 2004
!
!  Reference:
!
!    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
!    Numerical Simulation in Fluid Dynamics,
!    SIAM, 1998.
!
!  Parameters:
!
!    Input, real U(0:IMAX+1,0:JMAX+1), V(0:IMAX+1,0:JMAX+1),
!    P(0:IMAX+1,0:JMAX+1), the velocity and pressure fields.
!
!    Input, integer IMIN, IMAX, JMIN, JMAX, specify the ranges of
!    U, V and P to be written out.
!
!    Input, character ( len = * ) FILENAME, the name of the file
!    into which the data is to be written.
!
  use nrtype

  implicit none

  character ( len = * ), intent ( in ) :: filename
  integer i
  integer, intent ( in ) :: imax
  integer, intent ( in ) :: imin
  integer j
  integer, intent ( in ) :: jmax
  integer, intent ( in ) :: jmin
  real ( rp ), dimension (0:,0:), intent ( in ) :: p
  real ( rp ), dimension (0:,0:), intent ( in ) :: u
  real ( rp ), dimension (0:,0:), intent ( in ) :: v

  open ( unit = 1, file = filename, status = 'replace', form = 'formatted' )

  do j = jmin, jmax
    do i = imin, imax
      write ( 1, '(1x,i5,1x,i5,3(1x,g12.6))' ) i, j, u(i,j), v(i,j), p(i,j)
    end do
  end do

  close ( unit = 1 )

  return
end
