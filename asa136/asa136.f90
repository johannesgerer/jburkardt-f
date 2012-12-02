subroutine kmns ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, &
  iter, wss, ifault )

!*****************************************************************************80
!
!! KMNS carries out the K-means algorithm.
!
!  Discussion:
!
!    This routine attempts to divide M points in N-dimensional space into
!    K clusters so that the within cluster sum of squares is minimized.
!
!  Modified:
!
!    14 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(M,N), the points.
!
!    Input, integer ( kind = 4 ) M, the number of points.
!
!    Input, integer ( kind = 4 ) N, the number of spatial dimensions.
!
!    Input/output, real ( kind = 8 ) C(K,N), the cluster centers.
!
!    Input, integer ( kind = 4 ) K, the number of clusters.
!
!    Output, integer ( kind = 4 ) IC1(M), the cluster to which each point
!    is assigned.
!
!    Workspace, integer ( kind = 4 ) IC2(M), used to store the cluster which
!    each point is most likely to be transferred to at each step.
!
!    Output, integer ( kind = 4 ) NC(K), the number of points in each cluster.
!
!    Workspace, real ( kind = 8 ) AN1(K).
!
!    Workspace, real ( kind = 8 ) AN2(K).
!
!    Workspace, integer ( kind = 4 ) NCP(K).
!
!    Workspace, real ( kind = 8 ) D(M).
!
!    Workspace, integer ( kind = 4 ) ITRAN(K).
!
!    Workspace, integer ( kind = 4 ) LIVE(K).
!
!    Input, integer ( kind = 4 ) ITER, the maximum number of iterations allowed.
!
!    Output, real ( kind = 8 ) WSS(K), the within-cluster sum of squares
!    of each cluster.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no error was detected.
!    1, at least one cluster is empty after the initial assignment.  A better
!       set of initial cluster centers is needed.
!    2, the allowed maximum number off iterations was exceeded.
!    3, K is less than or equal to 1, or greater than or equal to M.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) aa
  real ( kind = 8 ) an1(k)
  real ( kind = 8 ) an2(k)
  real ( kind = 8 ) c(k,n)
  real ( kind = 8 ) d(m)
  real ( kind = 8 ) da
  real ( kind = 8 ) db
  real ( kind = 8 ) dc
  real ( kind = 8 ) dt(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic1(m)
  integer ( kind = 4 ) ic2(m)
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) il
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itran(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) live(k)
  integer ( kind = 4 ) nc(k)
  integer ( kind = 4 ) ncp(k)
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) temp
  real ( kind = 8 ) wss(k)

  ifault = 0

  if ( k <= 1 .or. m <= k ) then
    ifault = 3
    return
  end if
!
!  For each point I, find its two closest centers, IC1(I) and
!  IC2(I).  Assign the point to IC1(I).
!
  do i = 1, m

    ic1(i) = 1
    ic2(i) = 2

    do il = 1, 2
      dt(il) = 0.0D+00
      do j = 1, n
        da = a(i,j) - c(il,j)
        dt(il) = dt(il) + da * da
      end do
    end do

    if ( dt(2) < dt(1) ) then
      ic1(i) = 2
      ic2(i) = 1
      temp = dt(1)
      dt(1) = dt(2)
      dt(2) = temp
    end if

    do l = 3, k

      db = 0.0D+00
      do j = 1, n
        dc = a(i,j) - c(l,j)
        db = db + dc * dc
      end do

      if ( db < dt(2) ) then

        if ( dt(1) <= db ) then
          dt(2) = db
          ic2(i) = l
        else
          dt(2) = dt(1)
          ic2(i) = ic1(i)
          dt(1) = db
          ic1(i) = l
        end if

      end if

    end do

  end do
!
!  Update cluster centers to be the average of points contained within them.
!
  do l = 1, k
    nc(l) = 0
    do j = 1, n
      c(l,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    l = ic1(i)
    nc(l) = nc(l) + 1
    c(l,1:n) = c(l,1:n) + a(i,1:n)
  end do
!
!  Check to see if there is any empty cluster at this stage.
!
  ifault = 1

  do l = 1, k

    if ( nc(l) == 0 ) then
      ifault = 1
      return
    end if

  end do

  ifault = 0

  do l = 1, k

    aa = real ( nc(l), kind = 8 )
    c(l,1:n) = c(l,1:n) / aa
!
!  Initialize AN1, AN2, ITRAN and NCP.
!
!  AN1(L) = NC(L) / (NC(L) - 1)
!  AN2(L) = NC(L) / (NC(L) + 1)
!  ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
!           = 0 otherwise
!
!  In the optimal-transfer stage, NCP(L) stores the step at which
!  cluster L is last updated.
!
!  In the quick-transfer stage, NCP(L) stores the step at which
!  cluster L is last updated plus M.
!
    an2(l) = aa / ( aa + 1.0D+00 )

    if ( 1.0D+00 < aa ) then
      an1(l) = aa / ( aa - 1.0D+00 )
    else
      an1(l) = r8_huge ( )
    end if

    itran(l) = 1
    ncp(l) = -1

  end do

  indx = 0
  ifault = 2

  do ij = 1, iter
!
!  In this stage, there is only one pass through the data.   Each
!  point is re-allocated, if necessary, to the cluster that will
!  induce the maximum reduction in within-cluster sum of squares.
!
    call optra ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &
      live, indx )
!
!  Stop if no transfer took place in the last M optimal transfer steps.
!
    if ( indx == m ) then
      ifault = 0
      exit
    end if
!
!  Each point is tested in turn to see if it should be re-allocated
!  to the cluster to which it is most likely to be transferred,
!  IC2(I), from its present cluster, IC1(I).   Loop through the
!  data until no further change is to take place.
!
    call qtran ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, &
      itran, indx )
!
!  If there are only two clusters, there is no need to re-enter the
!  optimal transfer stage.
!
    if ( k == 2 ) then
      ifault = 0
      exit
    end if
!
!  NCP has to be set to 0 before entering OPTRA.
!
    do l = 1, k
      ncp(l) = 0
    end do

  end do
!
!  If the maximum number of iterations was taken without convergence,
!  IFAULT is 2 now.  This may indicate unforeseen looping.
!
  if ( ifault == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMNS - Warning!'
    write ( *, '(a)' ) '  Maximum number of iterations reached'
    write ( *, '(a)' ) '  without convergence.'
  end if
!
!  Compute the within-cluster sum of squares for each cluster.
!
  wss(1:k) = 0.0D+00
  c(1:k,1:n) = 0.0D+00

  do i = 1, m
    ii = ic1(i)
    c(ii,1:n) = c(ii,1:n) + a(i,1:n)
  end do

  do j = 1, n
    do l = 1, k
      c(l,j) = c(l,j) / real ( nc(l), kind = 8 )
    end do
    do i = 1, m
      ii = ic1(i)
      da = a(i,j) - c(ii,j)
      wss(ii) = wss(ii) + da * da
    end do
  end do

  return
end
subroutine optra ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &
  live, indx )

!*****************************************************************************80
!
!! OPTRA carries out the optimal transfer stage.
!
!  Discussion:
!
!    This is the optimal transfer stage.
!
!    Each point is re-allocated, if necessary, to the cluster that
!    will induce a maximum reduction in the within-cluster sum of
!    squares.
!
!  Modified:
!
!    14 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(M,N), the points.
!
!    Input, integer ( kind = 4 ) M, the number of points.
!
!    Input, integer ( kind = 4 ) N, the number of spatial dimensions.
!
!    Input/output, real ( kind = 8 ) C(K,N), the cluster centers.
!
!    Input, integer ( kind = 4 ) K, the number of clusters.
!
!    Input/output, integer ( kind = 4 ) IC1(M), the cluster to which each
!    point is assigned.
!
!    Input/output, integer ( kind = 4 ) IC2(M), used to store the cluster
!    which each point is most likely to be transferred to at each step.
!
!    Input/output, integer ( kind = 4 ) NC(K), the number of points in
!    each cluster.
!
!    Input/output, real ( kind = 8 ) AN1(K).
!
!    Input/output, real ( kind = 8 ) AN2(K).
!
!    Input/output, integer ( kind = 4 ) NCP(K).
!
!    Input/output, real ( kind = 8 ) D(M).
!
!    Input/output, integer ( kind = 4 ) ITRAN(K).
!
!    Input/output, integer ( kind = 4 ) LIVE(K).
!
!    Input/output, integer ( kind = 4 ) INDX, the number of steps since a
!    transfer took place.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) al1
  real ( kind = 8 ) al2
  real ( kind = 8 ) alt
  real ( kind = 8 ) alw
  real ( kind = 8 ) an1(k)
  real ( kind = 8 ) an2(k)
  real ( kind = 8 ) c(k,n)
  real ( kind = 8 ) d(m)
  real ( kind = 8 ) da
  real ( kind = 8 ) db
  real ( kind = 8 ) dc
  real ( kind = 8 ) dd
  real ( kind = 8 ) de
  real ( kind = 8 ) df
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic1(m)
  integer ( kind = 4 ) ic2(m)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) itran(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) live(k)
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) nc(k)
  integer ( kind = 4 ) ncp(k)
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) rr
!
!  If cluster L is updated in the last quick-transfer stage, it
!  belongs to the live set throughout this stage.   Otherwise, at
!  each step, it is not in the live set if it has not been updated
!  in the last M optimal transfer steps.
!
  do l = 1, k
    if ( itran(l) == 1) then
      live(l) = m + 1
    end if
  end do

  do i = 1, m

    indx = indx + 1
    l1 = ic1(i)
    l2 = ic2(i)
    ll = l2
!
!  If point I is the only member of cluster L1, no transfer.
!
    if ( 1 < nc(l1)  ) then
!
!  If L1 has not yet been updated in this stage, no need to
!  re-compute D(I).
!
      if ( ncp(l1) /= 0 ) then
        de = 0.0D+00
        do j = 1, n
          df = a(i,j) - c(l1,j)
          de = de + df * df
        end do
        d(i) = de * an1(l1)
      end if
!
!  Find the cluster with minimum R2.
!
     da = 0.0D+00
      do j = 1, n
        db = a(i,j) - c(l2,j)
        da = da + db * db
      end do
      r2 = da * an2(l2)

      do l = 1, k
!
!  If LIVE(L1) <= I, then L1 is not in the live set.   If this is
!  true, we only need to consider clusters that are in the live set
!  for possible transfer of point I.   Otherwise, we need to consider
!  all possible clusters.
!
        if ( ( i < live(l1) .or. i < live(l2) ) .and. &
               l /= l1 .and. l /= ll ) then

          rr = r2 / an2(l)

          dc = 0.0D+00
          do j = 1, n
            dd = a(i,j) - c(l,j)
            dc = dc + dd * dd
          end do

          if ( dc < rr ) then
            r2 = dc * an2(l)
            l2 = l
          end if

        end if

      end do
!
!  If no transfer is necessary, L2 is the new IC2(I).
!
      if ( d(i) <= r2 ) then

        ic2(i) = l2
!
!  Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and
!  L2, and update IC1(I) and IC2(I).
!
      else

        indx = 0
        live(l1) = m + i
        live(l2) = m + i
        ncp(l1) = i
        ncp(l2) = i
        al1 = real ( nc(l1), kind = 8 )
        alw = al1 - 1.0D+00
        al2 = real ( nc(l2), kind = 8 )
        alt = al2 + 1.0D+00
        do j = 1, n
          c(l1,j) = ( c(l1,j) * al1 - a(i,j) ) / alw
          c(l2,j) = ( c(l2,j) * al2 + a(i,j) ) / alt
        end do
        nc(l1) = nc(l1) - 1
        nc(l2) = nc(l2) + 1
        an2(l1) = alw / al1
        if ( 1.0D+00 < alw ) then
          an1(l1) = alw / ( alw - 1.0D+00 )
        else
          an1(l1) = r8_huge ( )
        end if
        an1(l2) = alt / al2
        an2(l2) = alt / ( alt + 1.0D+00 )
        ic1(i) = l2
        ic2(i) = l1

      end if

    end if

    if ( indx == m ) then
      return
    end if

  end do
!
!  ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
!  decreased by M before re-entering OPTRA.
!
  do l = 1, k
    itran(l) = 0
    live(l) = live(l) - m
  end do

  return
end
subroutine qtran ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &
  indx )

!*****************************************************************************80
!
!! QTRAN carries out the quick transfer stage.
!
!  Discussion:
!
!    This is the quick transfer stage.
!
!    IC1(I) is the cluster which point I belongs to.
!    IC2(I) is the cluster which point I is most likely to be
!    transferred to.
!
!    For each point I, IC1(I) and IC2(I) are switched, if necessary, to
!    reduce within-cluster sum of squares.  The cluster centers are
!    updated after each step.
!
!  Modified:
!
!    15 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(M,N), the points.
!
!    Input, integer ( kind = 4 ) M, the number of points.
!
!    Input, integer ( kind = 4 ) N, the number of spatial dimensions.
!
!    Input/output, real ( kind = 8 ) C(K,N), the cluster centers.
!
!    Input, integer ( kind = 4 ) K, the number of clusters.
!
!    Input/output, integer ( kind = 4 ) IC1(M), the cluster to which each
!    point is assigned.
!
!    Input/output, integer ( kind = 4 ) IC2(M), used to store the cluster
!    which each point is most likely to be transferred to at each step.
!
!    Input/output, integer ( kind = 4 ) NC(K), the number of points in
!    each cluster.
!
!    Input/output, real ( kind = 8 ) AN1(K).
!
!    Input/output, real ( kind = 8 ) AN2(K).
!
!    Input/output, integer ( kind = 4 ) NCP(K).
!
!    Input/output, real ( kind = 8 ) D(M).
!
!    Input/output, integer ( kind = 4 ) ITRAN(K).
!
!    Input/output, integer ( kind = 4 ) INDX, counts the number of steps
!    since the last transfer.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) al1
  real ( kind = 8 ) al2
  real ( kind = 8 ) alt
  real ( kind = 8 ) alw
  real ( kind = 8 ) an1(k)
  real ( kind = 8 ) an2(k)
  real ( kind = 8 ) c(k,n)
  real ( kind = 8 ) d(m)
  real ( kind = 8 ) da
  real ( kind = 8 ) db
  real ( kind = 8 ) dd
  real ( kind = 8 ) de
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic1(m)
  integer ( kind = 4 ) ic2(m)
  integer ( kind = 4 ) icoun
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) itran(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) nc(k)
  integer ( kind = 4 ) ncp(k)
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_huge
!
!  In the optimal transfer stage, NCP(L) indicates the step at which
!  cluster L is last updated.   In the quick transfer stage, NCP(L)
!  is equal to the step at which cluster L is last updated plus M.
!
  icoun = 0
  istep = 0

  do

    do i = 1, m

      icoun = icoun + 1
      istep = istep + 1
      l1 = ic1(i)
      l2 = ic2(i)
!
!  If point I is the only member of cluster L1, no transfer.
!
      if ( 1 < nc(l1) ) then
!
!  If NCP(L1) < ISTEP, no need to re-compute distance from point I to
!  cluster L1.   Note that if cluster L1 is last updated exactly M
!  steps ago, we still need to compute the distance from point I to
!  cluster L1.
!
        if ( istep <= ncp(l1) ) then

          da = 0.0D+00
          do j = 1, n
            db = a(i,j) - c(l1,j)
            da = da + db * db
          end do

          d(i) = da * an1(l1)

        end if
!
!  If NCP(L1) <= ISTEP and NCP(L2) <= ISTEP, there will be no transfer of
!  point I at this step.
!
        if ( istep < ncp(l1) .or. istep < ncp(l2) ) then

          r2 = d(i) / an2(l2)

          dd = 0.0D+00
          do j = 1, n
            de = a(i,j) - c(l2,j)
            dd = dd + de * de
          end do
!
!  Update cluster centers, NCP, NC, ITRAN, AN1 and AN2 for clusters
!  L1 and L2.   Also update IC1(I) and IC2(I).   Note that if any
!  updating occurs in this stage, INDX is set back to 0.
!
          if ( dd < r2 ) then

            icoun = 0
            indx = 0
            itran(l1) = 1
            itran(l2) = 1
            ncp(l1) = istep + m
            ncp(l2) = istep + m
            al1 = real ( nc(l1), kind = 8 )
            alw = al1 - 1.0D+00
            al2 = real ( nc(l2), kind = 8 )
            alt = al2 + 1.0D+00
            do j = 1, n
              c(l1,j) = ( c(l1,j) * al1 - a(i,j) ) / alw
              c(l2,j) = ( c(l2,j) * al2 + a(i,j) ) / alt
            end do
            nc(l1) = nc(l1) - 1
            nc(l2) = nc(l2) + 1
            an2(l1) = alw / al1
            if ( 1.0D+00 < alw ) then
              an1(l1) = alw / ( alw - 1.0D+00 )
            else
              an1(l1) = r8_huge ( )
            end if
            an1(l2) = alt / al2
            an2(l2) = alt / ( alt + 1.0D+00 )
            ic1(i) = l2
            ic2(i) = l1

          end if

        end if

      end if
!
!  If no re-allocation took place in the last M steps, return.
!
      if ( icoun == m ) then
        return
      end if

    end do

  end do

end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
