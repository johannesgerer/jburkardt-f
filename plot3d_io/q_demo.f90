program main

!*****************************************************************************80
!
!! MAIN is the main program for Q_DEMO.
!
!  Discussion:
!
!    Q_DEMO demonstrates a sample use of RB_3DM_QFILE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxgrid = 7
  integer ( kind = 4 ), parameter :: maxi = 121
  integer ( kind = 4 ), parameter :: maxj = 41
  integer ( kind = 4 ), parameter :: maxk = 21

  real ( kind = 8 ) alpha(maxgrid)
  real ( kind = 8 ) fsmach(maxgrid)
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
  real ( kind = 8 ) q(maxi,maxj,maxk,5,maxgrid)
  real ( kind = 8 ) qmax(5)
  real ( kind = 8 ) qmin(5)
  real ( kind = 8 ) re(maxgrid)
  real ( kind = 8 ) time(maxgrid)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Q_DEMO'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a (big) (nasty) PLOT3D file.'
  write ( *, '(a)' ) '  This file is a 3D Q file with multiple grids.'

  ierror = 0
  iunit = 1

  open ( unit = iunit, file = '3dm_q.dat', form = 'unformatted', &
    status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Q_DEMO - Fatal error!'
    write ( *, '(a)' ) '  Could not open the Q file.'
    stop
  end if

  call r8_b_3dm_qfile ( iunit, idim, jdim, kdim, maxi, maxj, maxk, &
    maxgrid, ngrid, fsmach, alpha, re, time, q, ierror )

  close ( unit = iunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Q_DEMO - Fatal error!'
    write ( *, '(a)' ) '  Error return from RB_3DM_QFILE:'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NGRID = ', ngrid
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, IDIM(I), JDIM(I), KDIM(I)'
  write ( *, '(a)' ) ' '

  do i = 1, ngrid

    write ( *, '(4i6)' ) i, idim(i), jdim(i), kdim(i)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, FSMACH(I), ALPHA(I), RE(I), TIME(I)'
  write ( *, '(a)' ) ' '

  do i = 1, ngrid

    write ( *, '(i6,4g14.6)' ) i, fsmach(i), alpha(i), re(i), time(i)

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
  write ( *, '(a)' ) '  Q ranges:'
  write ( *, '(a)' ) ' '
  do l = 1, 5
    write ( *, '(i6,2g14.6)' ) l, qmin(l), qmax(l)
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Q_DEMO'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
