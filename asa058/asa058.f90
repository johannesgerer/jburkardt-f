subroutine clustr ( x, d, dev, b, f, e, i, j, n, nz, k )

!*****************************************************************************80
!
!! CLUSTR uses the K-means algorithm to cluster data.
!
!  Discussion:
!
!    Given a matrix of I observations on J variables, the
!    observations are allocated to N clusters in such a way that the
!    within-cluster sum of squares is minimised.
!
!  Modified:
!
!    23 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Sparks.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Sparks,
!    Algorithm AS 58:
!    Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(I,J), the observed data.
!
!    Input/output, real ( kind = 8 ) D(K,J), the cluster centers.
!    On input, the user has chosen these.  On output, they have been
!    updated.
!
!    Output, real ( kind = 8 ) DEV(K), the sums of squared deviations
!    of observations from their cluster centers.
!
!    Output, integer ( kind = 4 ) B(I), indicates the cluster to which
!    each observation has been assigned.
!
!    Workspace, real ( kind = 8 ) F(I).
!
!    Output, integer ( kind = 4 ) E(K), the number of observations assigned
!    to each cluster.
!
!    Input, integer ( kind = 4 ) I, the number of observations.
!
!    Input, integer ( kind = 4 ) J, the number of variables.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Input, integer ( kind = 4 ) NZ, the minimum number of observations
!    which any cluster is allowed to have.
!
!    Input, integer ( kind = 4 ) K, the maximum number of clusters.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k

  integer ( kind = 4 ) b(i)
  real ( kind = 8 ), parameter :: big = 1.0D+10
  real ( kind = 8 ) d(k,j)
  real ( kind = 8 ) da
  real ( kind = 8 ) db
  real ( kind = 8 ) dc
  real ( kind = 8 ) de
  real ( kind = 8 ) dev(k)
  integer ( kind = 4 ) e(k)
  real ( kind = 8 ) f(i)
  real ( kind = 8 ) fl
  real ( kind = 8 ) fm
  real ( kind = 8 ) fq
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) ig
  integer ( kind = 4 ) ih
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) il
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz
  real ( kind = 8 ) x(i,j)

  e(1:n) = 0
!
!  For each observation, calculate the distance from each cluster
!  center, and assign to the nearest.
!
  do ic = 1, i

    f(ic) = 0.0D+00
    da = big

    do id = 1, n

      db = 0.0D+00
      do ie = 1, j
        dc = x(ic,ie) - d(id,ie)
        db = db + dc * dc
      end do

      if ( db < da ) then
        da = db
        b(ic) = id
      end if

    end do

    ig = b(ic)
    e(ig) = e(ig) + 1

  end do
!
!  Calculate the mean and sum of squares for each cluster.
!
  dev(1:n) = 0.0D+00
  d(1:n,1:j) = 0.0D+00

  do ic = 1, i
    ig = b(ic)
    d(ig,1:j) = d(ig,1:j) + x(ic,1:j)
  end do

  do ij = 1, j
    do ii = 1, n
      d(ii,ij) = d(ii,ij) / real ( e(ii), kind = 8 )
    end do
  end do

  do ij = 1, j
    do ik = 1, i
      il = b(ik)
      da = x(ik,ij) - d(il,ij)
      db = da * da
      f(ik) = f(ik) + db
      dev(il) = dev(il) + db
    end do
  end do

  do ik = 1, i
    il = b(ik)
    fl = e(il)
    if ( 1 < e(il) ) then
      f(ik) = f(ik) * fl / ( fl - 1.0D+00 )
    end if
  end do
!
!  Examine each observation in turn to see if it should be
!  reassigned to a different cluster.
!
  do

    iw = 0

    do ik = 1, i

      il = b(ik)
      ir = il
!
!  If the number of cluster points is less than or equal to the
!  specified minimum, NZ, then bypass this iteration.
!
      if ( nz < e(il) ) then

        fl = e(il)
        dc = f(ik)

        do in = 1, n

          if ( in /= il ) then

            fm = e(in)
            fm = fm / ( fm + 1.0D+00 )

            de = 0.0D+00
            do ip = 1, j
              da = x(ik,ip) - d(in,ip)
              de = de + da * da * fm
            end do

            if ( de < dc ) then
              dc = de
              ir = in
            end if

          end if

        end do
!
!  Reassignment is made here if necessary.
!
        if ( ir /= il ) then

          fq = e(ir)
          dev(il) = dev(il) - f(ik)
          dev(ir) = dev(ir) + dc
          e(ir) = e(ir) + 1
          e(il) = e(il) - 1

          do is = 1, j
            d(il,is) = ( d(il,is) * fl - x(ik,is) ) / ( fl - 1.0D+00 )
            d(ir,is) = ( d(ir,is) * fq + x(ik,is) ) / ( fq + 1.0D+00 )
          end do

          b(ik) = ir

          do it = 1, i

            ij = b(it)

            if ( ij == il .or. ij == ir ) then
              f(it) = 0.0D+00
              do iu = 1, j
                da = x(it,iu) - d(ij,iu)
                f(it) = f(it) + da * da
              end do
              fl = e(ij)
              f(it) = f(it) * fl / ( fl - 1.0D+00 )
            end if

          end do

          iw = iw + 1

        end if

      end if

    end do
!
!  If any reassignments were made on this pass, then do another pass.
!
    if ( iw == 0 ) then
      exit
    end if

  end do

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
