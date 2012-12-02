program main

!*****************************************************************************80
!
!! MAIN is the main program for XYZ_DEMO.
!
!  Discussion:
!
!    XYZ_DEMO demonstrates a sample use of RB_3DM_XYZFILE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxi = 121
  integer ( kind = 4 ), parameter :: maxj = 41
  integer ( kind = 4 ), parameter :: maxk = 21
  integer ( kind = 4 ), parameter :: maxgrid = 7

  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid
  real ( kind = 8 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) z(maxi,maxj,maxk,maxgrid)
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XYZ_DEMO'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a (big) (nasty) PLOT3D file.'
  write ( *, '(a)' ) '  This file is a 3D XYZ file with multiple grids.'

  ierror = 0
  iunit = 1

  open ( unit = iunit, file = '3dm_xyz.dat', form = 'unformatted', &
    status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_DEMO - Fatal error!'
    write ( *, '(a)' ) '  Could not open the XYZ direct access file.'
    stop
  end if

  call r8_b_3dm_xyzfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, x, y, z, ierror )

  close ( unit = iunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_DEMO - Warning!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NGRID = ', ngrid
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, IDIM(I), JDIM(I), KDIM(I)'
  write ( *, '(a)' ) ' '

  do i = 1, ngrid

    write ( *, '(4i6)' ) i, idim(i), jdim(i), kdim(i)

  end do

  xmin = x(1,1,1,1)
  xmax = x(1,1,1,1)
  ymin = y(1,1,1,1)
  ymax = y(1,1,1,1)
  zmin = z(1,1,1,1)
  zmax = z(1,1,1,1)

  do igrid = 1, ngrid
    do i = 1, idim(igrid)
      do j = 1, jdim(igrid)
        do k = 1, kdim(igrid)
          xmin = min ( xmin, x(i,j,k,igrid) )
          xmax = max ( xmax, x(i,j,k,igrid) )
          ymin = min ( ymin, y(i,j,k,igrid) )
          ymax = max ( ymax, y(i,j,k,igrid) )
          zmin = min ( zmin, z(i,j,k,igrid) )
          zmax = max ( zmax, z(i,j,k,igrid) )
        end do
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XYZ ranges:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  XMIN = ', xmin
  write ( *, '(a,g14.6)' ) '  XMAX = ', xmax
  write ( *, '(a,g14.6)' ) '  YMIN = ', ymin
  write ( *, '(a,g14.6)' ) '  YMAX = ', ymax
  write ( *, '(a,g14.6)' ) '  ZMIN = ', zmin
  write ( *, '(a,g14.6)' ) '  ZMAX = ', zmax
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XYZ_DEMO'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
