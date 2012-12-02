program main

!*****************************************************************************80
!
!! MAIN is the main program for PLOT3D_TO_AVS.
!
!  Discussion:
!
!    PLOT3D_TO_AVS converts 3D binary multiple grid XYZB and Q format 
!    to another format.
!
!    The program reads in the XYZB data first, since the "B" field contains
!    the blanking information.  This file is '3dm_xyzb.dat'.
!
!    Then it reads in the Q data from '3dm_qb.txt'.
!
!    Then it writes the XYZ data to the file 'avs_xyz.dat', but only that 
!    data corresponding to unblanked points.  All the X's are written first,
!    then Y, then Z.
!
!    Then it writes the unblanked U, V, and W fields ( 2, 3 and 4 ) 
!    to a file 'avs_uvw.dat'.
!
!    Then it writes the unblanked RHO fields to the file 'avs_rho.dat'.
!
!    The MAXI, MAXJ, MAXK and MAXGRID parameter statements below
!    are set for the particular data used for the demonstration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2007
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

  real ( kind = 4 ) alpha(maxgrid)
  integer ( kind = 4 ) b(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) fsmach(maxgrid)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  character ( len = 80 ) input_file1
  character ( len = 80 ) input_file2
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ngrid
  integer ( kind = 4 ) numpts
  character ( len = 80 ) output_file1
  character ( len = 80 ) output_file2
  character ( len = 80 ) output_file3
  real ( kind = 4 ) q(maxi,maxj,maxk,5,maxgrid)
  real ( kind = 4 ) qmax(5)
  real ( kind = 4 ) qmin(5)
  real ( kind = 4 ) re(maxgrid)
  real ( kind = 4 ) time(maxgrid)
  real ( kind = 4 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) ymax
  real ( kind = 4 ) ymin
  real ( kind = 4 ) z(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) zmax
  real ( kind = 4 ) zmin

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PLOT3D_TO_AVS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read PLOT3D information from two binary files:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    3dm_xyzb.dat: a 3D XYZB file with multiple grids,'
  write ( *, '(a)' ) '    3dm_qb.dat    a 3D Q file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Write the unblanked data into files AVS can handle:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    avs_xyz.dat: the XYZ data;'
  write ( *, '(a)' ) '    avs_uvw.dat: the UVW (Q(2),Q(3),Q(4)) data;'
  write ( *, '(a)' ) '    avs_rho.dat: the RHO (Q(1)) data.'

  ierror = 0

  iunit = 1
  input_file1 = '3dm_xyzb.dat'

  open ( unit = iunit, file = input_file1, form = 'unformatted', &
    status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT3D_TO_AVS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    write ( *, '(a)' ) '  Was expecting "' // trim ( input_file1 ) // '".'
    stop
  end if

  call r4_b_3dm_xyzbfile ( iunit, idim, jdim, kdim, maxi, maxj, &
    maxk, maxgrid, ngrid, x, y, z, b, ierror )

  close ( unit = iunit )

  if ( ierror /= 0 ) then
    write ( *, '(a,i8)' ) 'IERROR = ', ierror
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'NGRID = ', ngrid
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, IDIM(I), JDIM(I), KDIM(I)'
  write ( *, '(a)' ) ' '

  do i = 1, ngrid

    write ( *, '(a)' ) i, idim(i), jdim(i), kdim(i)

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
  write ( *, '(a)' ) 'XYZ ranges:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  XMIN = ', xmin
  write ( *, '(a,g14.6)' ) '  XMAX = ', xmax
  write ( *, '(a,g14.6)' ) '  YMIN = ', ymin
  write ( *, '(a,g14.6)' ) '  YMAX = ', ymax
  write ( *, '(a,g14.6)' ) '  ZMIN = ', zmin
  write ( *, '(a,g14.6)' ) '  ZMAX = ', zmax
!
!  Now read the Q data.
!
  ierror = 0
  iunit = 1
  input_file2 = '3dm_qb.txt'

  open ( unit = iunit, file = input_file2, form = 'unformatted', &
    status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT3D_TO_AVS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    write ( *, '(a)' ) '  Was expecting "' // trim ( input_file2 ) // '".'
    stop
  end if

  call r4_b_3dm_qfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, fsmach, alpha, re, time, q, ierror )

  close ( unit = iunit )

  if ( ierror /= 0 ) then
    write ( *, '(a,i8)' ) 'IERROR = ', ierror
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'NGRID = ', ngrid
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, IDIM(I), JDIM(I), KDIM(I)'
  write ( *, '(a)' ) ' '

  do i = 1, ngrid

    write ( *, '(4i8)' ) i, idim(i), jdim(i), kdim(i)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, FSMACH(I), ALPHA(I), RE(I), TIME(I)'
  write ( *, '(a)' ) ' '

  do i = 1, ngrid

    write ( *, '(i8,4g14.6)' ) i, fsmach(i), alpha(i), re(i), time(i)

  end do

  do l = 1, 5
    qmin(l) = q(1,1,1,1,1)
    qmax(l) = q(1,1,1,1,1)
    do igrid = 1, ngrid
      do i = 1, idim(igrid)
        do j = 1, jdim(igrid)
          do k = 1, kdim(igrid)
            qmin(l) = min ( qmin(l), q(i,j,k,l,igrid) )
            qmax(l) = max ( qmax(l), q(i,j,k,l,igrid) )
          end do
        end do
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Q ranges:'
  write ( *, '(a)' ) ' '
  do l = 1, 5
    write ( *, '(i8,2g14.6)' ) l, qmin(l), qmax(l)
  end do
!
!  Now write out to a file just those X, Y, and Z values which
!  are not blanked.
!
  output_file1 = 'avs_xyz.dat'

  open ( unit = iunit, file = output_file1, form = 'unformatted', &
    status = 'replace', access = 'direct', recl = 1, iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT3D_TO_AVS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the XYZ direct access file'
    write ( *, '(a)' ) '  "' // trim ( output_file1 ) // '".'
    stop
  end if

  numpts = 0
  do igrid = 1, ngrid
    do i = 1, idim(igrid)
      do j = 1, jdim(igrid)
        do k = 1, kdim(igrid)
          if ( b(i,j,k,igrid) /= 0 ) then
            write ( iunit ) x(i,j,k,igrid)
            numpts = numpts + 1
          end if
        end do
      end do
    end do
  end do

  write ( *, '(a,i8)' ) 'Number of unblanked points is ', numpts

  do igrid = 1, ngrid
    do i = 1, idim(igrid)
      do j = 1, jdim(igrid)
        do k = 1, kdim(igrid)
          if ( b(i,j,k,igrid) /= 0 ) then
            write ( iunit ) y(i,j,k,igrid)
          end if
        end do
      end do
    end do
  end do

  do igrid = 1, ngrid
    do i = 1, idim(igrid)
      do j = 1, jdim(igrid)
        do k = 1, kdim(igrid)
          if ( b(i,j,k,igrid) /= 0 ) then
            write ( iunit ) z(i,j,k,igrid)
          end if
        end do
      end do
    end do
  end do

  close ( unit = iunit )
!
!  Now write out to a file just those U, V, and W values which
!  are not blanked.
!
  output_file2 = 'avs_uvw.dat'

  open ( unit = iunit, file = output_file2, form = 'unformatted', &
    status = 'replace', access = 'direct', recl = 1, iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT3D_TO_AVS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the UVW direct access file:'
    write ( *, '(a)' ) '  "' // trim ( output_file2 ) // '".'
    stop
  end if

  do igrid = 1, ngrid
    do i = 1, idim(igrid)
      do j = 1, jdim(igrid)
        do k = 1, kdim(igrid)
          if ( b(i,j,k,igrid) /= 0 ) then
            write ( iunit ) q(i,j,k,2,igrid)
          end if
        end do
      end do
    end do
  end do

  do igrid = 1, ngrid
    do i = 1, idim(igrid)
      do j = 1, jdim(igrid)
        do k = 1, kdim(igrid)
          if ( b(i,j,k,igrid) /= 0 ) then
            write ( iunit ) q(i,j,k,3,igrid)
          end if
        end do
      end do
    end do
  end do

  do igrid = 1, ngrid
    do i = 1, idim(igrid)
      do j = 1, jdim(igrid)
        do k = 1, kdim(igrid)
          if ( b(i,j,k,igrid) /= 0 ) then
            write ( iunit ) q(i,j,k,4,igrid)
          end if
        end do
      end do
    end do
  end do

  close ( unit = iunit )
!
!  Now write out to a file just those RHO values which
!  are not blanked.
!
  output_file3 = 'avs_rho.dat'

  open ( unit = iunit, file = output_file3, form = 'unformatted', &
    status = 'replace', access = 'direct', recl = 1, iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT3D_TO_AVS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the RHO direct access file:'
    write ( *, '(a)' ) '  "' // trim ( output_file3 ) // '".'
    stop
  end if

  do igrid = 1, ngrid
    do i = 1, idim(igrid)
      do j = 1, jdim(igrid)
        do k = 1, kdim(igrid)
          if ( b(i,j,k,igrid) /= 0 ) then
            write ( iunit ) q(i,j,k,1,igrid)
          end if
        end do
      end do
    end do
  end do

  close ( unit = iunit )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PLOT3D_TO_AVS'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r4_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_DIM reads a binary 3D multiple grid file for the dimensions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has 
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID), the number 
!    of nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) maxgrid

  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid

  ierror = 0

  read ( iunit, iostat = ios ) ( idim(igrid), jdim(igrid), &
    kdim(igrid), igrid = 1, ngrid )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3DM_DIM - Error!'
    write ( *, '(a)' ) '  A read error occurred.'
    return
  end if

  return
end
subroutine r4_b_3dm_ngrid ( iunit, ngrid, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_NGRID reads a binary 3D multiple grid file for the number of grids.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has 
!    been opened.
!
!    Output, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ngrid

  ierror = 0

  read ( iunit, iostat = ios ) ngrid

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_B_3DM_DIM - Error!'
    write ( *, '(a)' ) '  A read error occurred.'
    return
  end if

  return
end
subroutine r4_b_3dm_q ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, fsmach, alpha, re, time, q, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_Q reads a binary 3D multiple grid Q file for the parameters and Q data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has 
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID), the number of 
!    nodes in the X, Y, and Z directions in each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of Q.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) FSMACH(MAXGRID), ALPHA(MAXGRID), RE(MAXGRID), TIME(MAXGRID),
!    the values of the free stream Mach number, the angle of attack, in degrees,
!    the Reynolds number for the flow, and the time, for each grid.
!
!    Output, real ( kind = 4 ) Q(MAXI,MAXJ,MAXK,5,MAXGRID), 
!    the Q values of the nodes for each grid.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    nonzero, a read error occurred.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  real ( kind = 4 ) alpha(maxgrid)
  real ( kind = 4 ) fsmach(maxgrid)
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
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) q(maxi,maxj,maxk,5,maxgrid)
  real ( kind = 4 ) re(maxgrid)
  real ( kind = 4 ) time(maxgrid)

  ierror = 0

  do igrid = 1, ngrid

    read ( iunit, iostat = ios ) fsmach(igrid), alpha(igrid), re(igrid), &
      time(igrid)

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DM_Q - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

    read ( iunit, iostat = ios ) (((( q(i,j,k,l,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), l = 1, 5 )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DM_Q - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

  end do

  return
end
subroutine r4_b_3dm_qfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, fsmach, alpha, re, time, q, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_QFILE reads a binary 3D multiple grid Q file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has 
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID), the number of 
!    nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of Q.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Output, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) FSMACH(MAXGRID), ALPHA(MAXGRID), RE(MAXGRID), TIME(MAXGRID),
!    the values of the free stream Mach number, the angle of attack, in degrees,
!    the Reynolds number for the flow, and the time, for each grid.
!
!    Output, real ( kind = 4 ) Q(MAXI,MAXJ,MAXK,5,MAXGRID), 
!    the Q values of the nodes for each grid.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading NGRID.
!    2, MAXGRID < NGRID.
!    3, an error occurred while reading the dimensions.
!    4, some IDIM, JDIM or KDIM is greater than MAXI, MAXJ or MAXK.
!    5, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  real ( kind = 4 ) alpha(maxgrid)
  real ( kind = 4 ) fsmach(maxgrid)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) q(maxi,maxj,maxk,5,maxgrid)
  real ( kind = 4 ) re(maxgrid)
  real ( kind = 4 ) time(maxgrid)

  ierror = 0

  call r4_b_3dm_ngrid ( iunit, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  if ( maxgrid < ngrid ) then
    ierror = 2
    return
  end if

  call r4_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  do i = 1, ngrid

    if ( maxi < idim(i) .or. maxj < jdim(i) .or. maxk < kdim(i) ) then
      ierror = 4
      return
    end if

  end do

  call r4_b_3dm_q ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, fsmach, alpha, re, time, q, ierror )

  if ( ierror /= 0 ) then
    ierror = 5
    return
  end if

  return
end
subroutine r4_b_3dm_xyzb ( iunit, idim, jdim, kdim, maxi, maxj, maxk, maxgrid, &
  ngrid, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_XYZB reads a binary 3D multiple grid XYZB file for the XYZB data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Input, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID), the number of
!    nodes in the X, Y, and Z directions in each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of X, Y and Z.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Input, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 4 ) B(MAXI,MAXJ,MAXK,MAXGRID), the blanking information.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  integer ( kind = 4 ) b(maxi,maxj,maxk,maxgrid)
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
  real ( kind = 4 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  do igrid = 1, ngrid

    read ( iunit, iostat = ios ) ((( x(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( y(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( z(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) ), &
      ((( b(i,j,k,igrid), i = 1, idim(igrid) ), &
      j = 1, jdim(igrid) ), k = 1, kdim(igrid) )

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4_B_3DM_XYZB - Fatal error!'
      write ( *, '(a)' ) '  A read error occurred.'
      return
    end if

  end do

  return
end
subroutine r4_b_3dm_xyzbfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
  maxgrid, ngrid, x, y, z, b, ierror )

!*****************************************************************************80
!
!! R4_B_3DM_XYZBFILE reads a binary 3D multiple grid XYZB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit on which the file has
!    been opened.
!
!    Output, integer ( kind = 4 ) IDIM(MAXGRID), JDIM(MAXGRID), KDIM(MAXGRID), the number of
!    nodes in the X, Y, and Z directions for each grid.
!
!    Input, integer ( kind = 4 ) MAXI, MAXJ, MAXK, the maximum dimension for the first three
!    components of X, Y and Z.
!
!    Input, integer ( kind = 4 ) MAXGRID, the maximum value of NGRID.
!
!    Output, integer ( kind = 4 ) NGRID, the number of grids.
!
!    Output, real ( kind = 4 ) X(MAXI,MAXJ,MAXK,MAXGRID),
!    Y(MAXI,MAXJ,MAXK,MAXGRID),
!    Z(MAXI,MAXJ,MAXK,MAXGRID), the X, Y and Z coordinates of the nodes
!    for each grid.
!
!    Output, integer ( kind = 4 ) B(MAXI,MAXJ,MAXK,MAXGRID), the blanking information.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, an error occurred while reading NGRID.
!    2, MAXGRID < NGRID.
!    3, an error occurred while reading the dimensions.
!    4, some IDIM, JDIM or KDIM is greater than MAXI, MAXJ or MAXK.
!    5, an error occurred while reading the data.
!
  implicit none

  integer ( kind = 4 ) maxgrid
  integer ( kind = 4 ) maxi
  integer ( kind = 4 ) maxj
  integer ( kind = 4 ) maxk

  integer ( kind = 4 ) b(maxi,maxj,maxk,maxgrid)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idim(maxgrid)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jdim(maxgrid)
  integer ( kind = 4 ) kdim(maxgrid)
  integer ( kind = 4 ) ngrid
  real ( kind = 4 ) x(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) y(maxi,maxj,maxk,maxgrid)
  real ( kind = 4 ) z(maxi,maxj,maxk,maxgrid)

  ierror = 0

  call r4_b_3dm_ngrid ( iunit, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    return
  end if

  if ( maxgrid < ngrid ) then
    ierror = 2
    return
  end if

  call r4_b_3dm_dim ( iunit, idim, jdim, kdim, maxgrid, ngrid, ierror )

  if ( ierror /= 0 ) then
    ierror = 3
    return
  end if

  do i = 1, ngrid

    if ( maxi < idim(i) .or. maxj < jdim(i) .or. maxk < kdim(i) ) then
      ierror = 4
      return
    end if

  end do

  call r4_b_3dm_xyzb ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, x, y, z, b, ierror )

  if ( ierror /= 0 ) then
    ierror = 5
    return
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
