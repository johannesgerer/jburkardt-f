subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
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
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine cludia ( m, t, p, nc, e, d )

!*****************************************************************************80
!
!! CLUDIA clusters data for which a distance matrix has been supplied.
!
!  Discussion:
!
!    This routine requires that an initial clustering be provided.
!
!    The distance matrix T(M,M) can contain the distances or the squares
!    of distances between objects, so that an L1 or L2 minimization
!    is performed.
!
!  Modified:
!
!    26 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 117-118,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of objects to cluster.
!
!    Input, real ( kind = 8 ) T(M,M), a distance matrix which records the
!    distances, or the squares of distances, between all pairs of
!    objects.  T(M,M) should be symmetric, and the diagonal entries
!    should be 0.
!
!    Input/output, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters.
!
!    Output, real ( kind = 8 ) E(NC), the sum, for each cluster, of the pairwise
!    distances between points in the cluster.
!
!    Output, real ( kind = 8 ) D, the sum of the entries in E.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nc

  real ( kind = 8 ) a
  real ( kind = 8 ) aa
  real ( kind = 8 ) b
  real ( kind = 8 ) bb
  real ( kind = 8 ) c(nc)
  real ( kind = 8 ) d
  real ( kind = 8 ) e(nc)
  real ( kind = 8 ) f
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(nc)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s(m)
  real ( kind = 8 ) t(m,m)
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
  real ( kind = 8 ) y
  real ( kind = 8 ) yy

  d = 0.0D+00

  e(1:nc) = 0.0D+00

  do j = 1, nc

    r = 0
    do i = 1, m
      if ( p(i) == j ) then
        r = r + 1
        s(r) = i
      end if
    end do

    if ( 1 < r ) then

      f = 0.0D+00

      do k = 1, r - 1
        u = s(k)
        do l = k + 1, r
          v = s(l)
          f = f + t(u,v)
        end do
      end do

      f = f / real ( r, kind = 8 )
      e(j) = f
      d = d + f

    end if

    q(j) = r

  end do

  i = 0
  it = 0

  do

    i = i + 1

    if ( m < i ) then
      i = i - m
    end if

    if ( it == m ) then
      exit
    end if

    r = p(i)
    u = q(r)

    if ( u <= 1 ) then
      cycle
    end if

    f = real ( u, kind = 8 )

    c(1:nc) = 0.0D+00

    do h = 1, m
      if ( h /= i ) then
        v = p(h)
        c(v) = c(v) + t(i,h)
      end if
    end do

    a = ( f * e(r) - c(r) ) / ( f - 1.0D+00 )
    aa = e(r) - a
    bb = huge ( bb )

    do j = 1, nc

      if ( j /= r ) then

        f = real ( q(j), kind = 8 )
        y = ( f * e(j) + c(j) ) / ( f + 1.0D+00 )
        yy = y - e(j)

        if ( yy < bb ) then
          bb = yy
          b = y
          u = j
        end if

      end if

    end do

    if ( bb < aa ) then
      it = 0
      d = d - aa + bb
      e(r) = a
      e(u) = b
      q(r) = q(r) - 1
      q(u) = q(u) + 1
      p(i) = u
    else
      it = it + 1
    end if

  end do

  return
end
subroutine clusta ( m, n, x, w, p, nc, s, e, d, iflag )

!*****************************************************************************80
!
!! CLUSTA solves the multiple location problem in N dimensions.
!
!  Discussion:
!
!    The algorithm attempts to determine the position of a set of
!    points Y so as to minimize the objective function
!
!      F = sum ( 1 <= J <= NC ) sum ( point I in cluster J )
!        W(I) * dist ( X(I), Y(J) )
!
!    where dist ( X, Y ) is the usual Euclidean distance.
!
!
!    This code was not originally designed to handle the case
!    of a single cluster, which seems to me a perfectly reasonable
!    operation.  I have attempted to correct this omission.
!
!  Modified:
!
!    24 November 2003
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 136-138,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, real ( kind = 8 ) W(M), the weights associated with each data point.
!
!    Input/output, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters.
!
!    Output, real ( kind = 8 ) S(NC,N), the location of the cluster centers.
!    These points are NOT required to be data points.
!
!    Output, real ( kind = 8 ) E(NC), the per-cluster objective functions.
!
!    Output, real ( kind = 8 ) D, the total objective function, which is just
!    the sum of the per-cluster objective functions.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no errors were found.
!    nonzero, an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bu
  real ( kind = 8 ) d
  real ( kind = 8 ) e(nc)
  real ( kind = 8 ), parameter :: eps = 0.001D+00
  real ( kind = 8 ) f
  real ( kind = 8 ) g(m)
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: itmax = 100
  logical itrue
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(nc)
  integer ( kind = 4 ) r
  real ( kind = 8 ) ra
  real ( kind = 8 ) s(nc,n)
  real ( kind = 8 ) sa(n)
  real ( kind = 8 ) sb(n)
  real ( kind = 8 ) t(n)
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) w_sum
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: y

  iflag = 0
!
!  Check cluster assignments.
!
  do i = 1, m
    if ( p(i) < 1 .or. nc < p(i) ) then
      iflag = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CLUSTA - Fatal error!'
      write ( *, '(a)' ) '  Illegal cluster assignment.'
      return
    end if
  end do
!
!  Determine the cluster populations.
!
  call cluster_population ( m, p, nc, q )
!
!  Fail if any cluster is empty.
!
  do j = 1, nc
    if ( q(j) == 0 ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CLUSTA - Fatal error!'
      write ( *, '(a)' ) '  Empty clusters were found.'
      return
    end if
  end do
!
!  Handle the special case of one cluster.
!
  if ( nc == 1 ) then

    is = 0

    call standn ( m, n, x, w, s, eps, itmax, is, d )
    e(1) = d

    return

  end if

  itrue = .true.
  i = 0
  it = 0
  d = 0.0D+00

  do

    b = huge ( b )
    u = 0

    do j = 1, nc

      if ( itrue ) then

      else

        if ( j /= r ) then
          p(i) = j
        else
          p(i) = 0
        end if

      end if
!
!  Call STANDN to find the point T which
!  minimizes the objective function for cluster J.
!
      v = 0
      do h = 1, m
        if ( p(h) == j ) then
          v = v + 1
          g(v) = w(h)
        end if
      end do

      allocate ( y(1:v,1:n) )

      v = 0
      do h = 1, m
        if ( p(h) == j ) then
          v = v + 1
          y(v,1:n) = x(h,1:n)
        end if
      end do

      is = 1

      t(1:n) = 0.0D+00
      w_sum = 0.0D+00
      do h = 1, m
        if ( p(h) == j ) then
          t(1:n) = t(1:n) + w(h) * x(h,1:n)
          w_sum = w_sum + w(h)
        end if
      end do
      t(1:n) = t(1:n) / w_sum

      call standn ( v, n, y, g, t, eps, itmax, is, f )

      deallocate ( y )
!
!  On an initial step,
!  set the cluster energy E(J) to F,
!  set the cluster center S(J,*) to T(*),
!  add the cluster energy to the total energy.
!
      if ( itrue ) then

        e(j) = f
        d = d + f
        s(j,1:n) = t(1:n)
!
!  Set the energy and center for the active cluster J.
!
      else if ( j == r ) then

        a = f
        sa(1:n) = t(1:n)
!
!  For a nonactive cluster, save information for the one with
!  the smallest objective function.
!
      else if ( f <= b ) then

        b = f
        u = j
        sb(1:n) = t(1:n)

      end if

    end do

    if ( itrue ) then

      itrue = .false.

    else

      bu = b - e(u)

      ra = e(r) - a

      if ( ra <= bu ) then

        it = it + 1
        p(i) = r

      else

        it = 0
        e(r) = a
        e(u) = b
        d = d - ra + bu
        p(i) = u
        q(r) = q(r) - 1
        q(u) = q(u) + 1
        s(r,1:n) = sa(1:n)
        s(u,1:n) = sb(1:n)

      end if

    end if

    do

      i = i + 1

      if ( m < i ) then
        i = i - m
      end if

      if ( m <= it ) then
        return
      end if

      r = p(i)

      if ( q(r) /= 1 ) then
        exit
      end if

    end do

  end do

  return
end
subroutine cluster_centroids ( m, n, x, nc, p, s )

!*****************************************************************************80
!
!! CLUSTER_CENTROIDS determines the centroids of a clustering.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 65-67,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Input, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Output, real ( kind = 8 ) S(NC,N), the cluster centroids.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(nc)
  integer ( kind = 4 ) r
  real ( kind = 8 ) s(nc,n)
  real ( kind = 8 ) x(m,n)
!
!  Determine the cluster populations.
!
  call cluster_population ( m, p, nc, q )
!
!  Determine the centroid of each cluster.
!
  s(1:nc,1:n) = 0.0D+00
  do i = 1, m
    r = p(i)
    s(r,1:n) = s(r,1:n) + x(i,1:n)
  end do

  do j = 1, nc
    s(j,1:n) = s(j,1:n) / real ( max ( q(j), 1 ), kind = 8 )
  end do

  return
end
subroutine cluster_medians ( m, n, x, nc, p, median )

!*****************************************************************************80
!
!! CLUSTER_MEDIANS determines the medians of a clustering.
!
!  Modified:
!
!    12 May 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 65-67,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Input, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Output, real ( kind = 8 ) MEDIAN(NC,N), the cluster medians.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) median(nc,n)
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) r_pop
  real ( kind = 8 ) temp(m)
  real ( kind = 8 ) x(m,n)

  do r = 1, nc

    r_pop = 0
    do i = 1, m
      if ( p(i) == r ) then
        r_pop = r_pop + 1
      end if
    end do

    if ( r_pop == 0 ) then
      median(r,1:n) = 0.0D+00
      cycle
    end if
!
!  Compute the median of component J in cluster R.
!
    do j = 1, n

      k = 0

      do i = 1, m

        if ( p(i) == r ) then
          k = k + 1
          temp(k) = x(i,j)
        end if

      end do

      call r8vec_sort_bubble_a ( r_pop, temp )

      if ( mod ( r_pop, 2 ) == 1 ) then
        median(r,j) = temp((r_pop+1)/2)
      else
        median(r,j) = ( temp(r_pop) + temp((r_pop/2)+1) ) / 2.0D+00
      end if

    end do

  end do

  return
end
subroutine cluster_median_distance ( m, n, x, nc, p, d )

!*****************************************************************************80
!
!! CLUSTER_MEDIAN_DISTANCE finds the cluster median distance.
!
!  Discussion:
!
!    The cluster median distance is the sum of the in-cluster
!    median distances.
!
!    The in-cluster median distance is the sum of the L1 norms
!    of the difference between each element of the cluster and
!    the median value of the cluster.
!
!  Modified:
!
!    12 May 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 65-67,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Input, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Output, real ( kind = 8 ) D, the cluster median distance.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) r
  real ( kind = 8 ) median(nc,n)
  real ( kind = 8 ) x(m,n)
!
!  Get the cluster medians.
!
  call cluster_medians ( m, n, x, nc, p, median )
!
!  Sum the in-cluster distances.
!
  d = 0.0D+00

  do i = 1, m
    r = p(i)
    d = d + sum ( abs ( median(r,1:n) - x(i,1:n) ) )
  end do

  return
end
subroutine cluster_population ( m, p, nc, q )

!*****************************************************************************80
!
!! CLUSTER_POPULATION sets the cluster populations from the assignment array.
!
!  Modified:
!
!    29 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Output, integer ( kind = 4 ) Q(NC), the cluster populations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nc

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(nc)

  q(1:nc) = 0

  do i = 1, m
    j = p(i)
    q(j) = q(j) + 1
  end do

  return
end
subroutine cluster_variance ( m, n, x, nc, p, e )

!*****************************************************************************80
!
!! CLUSTER_VARIANCE determines the variances associated with a clustering.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 65-67,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Input, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Output, real ( kind = 8 ) E(NC), the cluster variances.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  real ( kind = 8 ) e(nc)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) r
  real ( kind = 8 ) s(nc,n)
  real ( kind = 8 ) x(m,n)
!
!  Get the cluster centroids.
!
  call cluster_centroids ( m, n, x, nc, p, s )

  e(1:nc) = 0.0D+00

  do i = 1, m
    r = p(i)
    e(r) = e(r) + sum ( ( s(r,1:n) - x(i,1:n) )**2 )
  end do

  return
end
subroutine colper ( a, m, n, p, power )

!*****************************************************************************80
!
!! COLPER seeks a column permutation which maximizes the "bond energy".
!
!  Discussion:
!
!    An M by N integer matrix A is given.  The value of A(I,J) is
!    supposed to indicate some degree of relationship between the
!    I-th row and J-th column, with 0 meaning no relationship, and
!    larger values suggesting a stronger relationship.  The matrix
!    can be assumed to contain ordinal data, such as a rating
!    from 0 to 10.
!
!    The objective function to be maximized is:
!
!      sum ( 1 <= I <= M ) sum ( 1 <= K <= N )
!        B(I,P(K)) * ) ( B(I,P(K-1)) + B(I,P(K+1)) )
!
!    where, if K-1 < 1 or N < K+1, we omit that term.
!
!  Modified:
!
!    06 May 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, pages 205-206,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A(M,N), the data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!
!    Output, integer ( kind = 4 ) P(N), a permutation of the column indices
!    which maximizes the objective function.
!
!    Input, real ( kind = 8 ) POWER, the power to which the entries of A should
!    be raised before processing.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  real ( kind = 8 ) b(m,n)
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) p(n)
  real ( kind = 8 ) power
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v

  b(1:m,1:n) = ( real ( a(1:m,1:n), kind = 8 ) )**power

  do k = 1, n - 1

    h = 0.0D+00

    do u = k + 1, n

      q = p(u)

      do v = 1, k + 1

        g = 0.0D+00

        if ( 1 < v ) then
          l = p(v-1)
          g = g + dot_product ( b(1:m,q), b(1:m,l) )
        end if

        if ( v < n ) then
          l = p(v)
          g = g + dot_product ( b(1:m,q), b(1:m,l) )
        end if

        if ( h < g ) then
          h = g
          r = v
          s = u
        end if

      end do

    end do

    if ( r /= s ) then

      t = p(s)

      do u = s, r+1, -1
        p(u) = p(u-1)
      end do

      p(r) = t

    end if

  end do

  return
end
subroutine data_d_read ( file_name, m, n, x )

!*****************************************************************************80
!
!! DATA_D_READ reads a real data set stored in a file.
!
!  Discussion:
!
!    The data set can be thought of as a real M by N array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row at a time.
!
!    Each row begins on a new line, but may extend over more than
!    one line.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Modified:
!
!    11 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Input, integer ( kind = 4 ) M, the number of data items.
!
!    Input, integer ( kind = 4 ) N, the dimension of each data item.
!
!    Output, real ( kind = 8 ) X(M,N), the data values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) line_num
  real ( kind = 8 ) value
  real ( kind = 8 ) x(m,n)

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old' )

  x(1:m,1:n) = huge ( x(1,1) )

  i = 1
  j = 0
  line_num = 0

  do
!
!  Have we read enough data?
!
    if ( i == m .and. j == n ) then
      exit
    end if
!
!  Have we read too much data?
!
    if ( m < i .or. n < j ) then
      exit
    end if
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else
!
!  LAST points to the last character associated with the previous
!  data value read from the line.
!
      last = 0

      do
!
!  Try to read another value from the line.
!
        call s_to_r8 ( line(last+1:), value, ierror, length )
!
!  If we could not read a new value, it's time to read a new line.
!
        if ( ierror /= 0 ) then
          exit
        end if
!
!  Update the pointer.
!
        last = last + length
!
!  If we read a new value, where do we put it?
!
        j = j + 1

        if ( n < j ) then
          j = 1
          i = i + 1
          if ( m < i ) then
            exit
          end if
        end if

        x(i,j) = value
!
!  If you reached the end of the row, it's time to read a new line.
!
        if ( j == n ) then
          exit
        end if

      end do

    end if

  end do

  close ( unit = input )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DATA_D_READ:'
  write ( *, '(a,i6)' ) '  Number of lines read was ', line_num

  return
end
subroutine data_d_print ( m, n, x, title )

!*****************************************************************************80
!
!! DATA_D_PRINT prints a real data set.
!
!  Discussion:
!
!    The data set can be thought of as a real M by N array.
!
!    Each row of the array corresponds to one data "item".
!
!  Modified:
!
!    29 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of data items.
!
!    Input, integer ( kind = 4 ) N, the dimension of each data item.
!
!    Input, real ( kind = 8 ) X(M,N), the data values.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  character ( len = * ) title
  real ( kind = 8 ) x(m,n)

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n

  call r8mat_print ( m, n, x, ' ' )

  return
end
subroutine data_d_show ( m, n, x, j1, j2, rows, columns )

!*****************************************************************************80
!
!! DATA_D_SHOW makes a typewriter plot of a real data set.
!
!  Modified:
!
!    11 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 145,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of data items.
!
!    Input, integer ( kind = 4 ) N, the dimension of the data items.
!
!    Input, real ( kind = 8 ) X(M,N), the data items.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns of data corresponding to
!    X and Y.
!
!    Input, integer ( kind = 4 ) ROWS, COLUMNS, the number of rows and columns
!    of "pixels" to use.
!
  implicit none

  integer ( kind = 4 ) columns
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rows

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  character string(0:rows+1,0:columns+1)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) x1_max
  real ( kind = 8 ) x1_min
  real ( kind = 8 ) x2_max
  real ( kind = 8 ) x2_min

  if ( rows <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  ROWS <= 0.'
    stop
  end if

  if ( columns <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  COLUMNS <= 0.'
    stop
  end if

  if ( j1 < 1 .or. n < j1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of J1.'
    write ( *, '(a)' ) '  J1 = ', j1
    stop
  end if

  if ( j2 < 1 .or. n < j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of J2.'
    write ( *, '(a)' ) '  J2 = ', j2
    stop
  end if

  string(0:rows+1,0:columns+1) = ' '

  do i = 1, rows
    string(i,0) = '.'
    string(i,columns+1) = '.'
  end do

  do j = 1, columns
    string(0,j) = '.'
    string(rows+1,j) = '.'
  end do

  string(0,0) = '.'
  string(0,columns+1) = '.'
  string(rows+1,0) = '.'
  string(rows+1,columns+1) = '.'

  x1_min = minval ( x(1:m,j1) )
  x1_max = maxval ( x(1:m,j1) )
  x2_min = minval ( x(1:m,j1) )
  x2_max = maxval ( x(1:m,j2) )

  do i = 1, m

    i2 = 1 + nint ( &
      real ( columns, kind = 8 ) * ( x(i,j1) - x1_min ) / ( x1_max - x1_min ) )
    i2 = max ( i2, 1 )
    i2 = min ( i2, columns )

    i1 = 1 + nint ( &
      real ( rows, kind = 8    ) * ( x(i,j2) - x2_min ) / ( x2_max - x2_min ) )
    i1 = max ( i1, 1 )
    i1 = min ( i1, rows )

    if ( string(i1,i2) == ' ' ) then
      string(i1,i2) = '*'
    else if ( string(i1,i2) == '*' ) then
      string(i1,i2) = '@'
    end if

  end do

  do i = rows+1, 0, -1
    write ( *, '(2x,78a1)' ) ( string(i,j), j = 0, columns+1 )
  end do

  return
end
subroutine data_size ( file_name, m, n )

!*****************************************************************************80
!
!! DATA_SIZE counts the size of a data set stored in a file.
!
!  Discussion:
!
!    Blank lines and comment lines (which begin with '#') are ignored).
!
!    All other lines are assumed to represent data items.
!
!    This routine assumes that each data line contains exactly the
!    same number of values, which are separated by spaces.
!
!    (This means this routine cannot handle cases where a data item
!    extends over more than one line, or cases where data is squeezed
!    together with no spaces, or where commas are used as separators,
!    but with space on both sides.)
!
!  Modified:
!
!    11 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Output, integer ( kind = 4 ) M, the number of nonblank, noncomment lines.
!
!    Output, integer ( kind = 4 ) N, the number of values per line.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) n_word

  m = 0
  n_max = - huge ( n_max )
  n_min = huge ( n_min )

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old' )

  do

    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      m = m + 1

      call s_word_count ( line, n_word )

      n_max = max ( n_max, n_word )
      n_min = min ( n_min, n_word )

    end if

  end do

  if ( n_max /= n_min ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Number of words per line varies.'
    write ( *, '(a,i6)' ) '  Minimum is ', n_min
    write ( *, '(a,i6)' ) '  Maximum is ', n_max
    n = 0
  else
    n = n_min
  end if

  close ( unit = input )

  return
end
subroutine dif_inverse ( n, a )

!*****************************************************************************80
!
!! DIF_INVERSE returns the inverse of the second difference matrix.
!
!  Formula:
!
!    if ( I <= J )
!      A(I,J) = I * (N-J+1) / (N+1)
!    else
!      A(I,J) = J * (N-I+1) / (N+1)
!
!  Example:
!
!    N = 4
!
!            4 3 2 1
!    (1/5) * 3 6 4 2
!            2 4 6 3
!            1 2 3 4
!
!  Properties:
!
!    A is symmetric.
!
!    Because A is symmetric, it is normal, so diagonalizable.
!
!    A is "persymmetric" ( A(I,J) = A(N+1-J,N+1-I) ).
!
!    The determinant of A is 1 / (N+1).
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( i <= j ) then
        a(i,j) = real ( i * ( n - j + 1 ), kind = 8 ) / real ( n + 1, kind = 8 )
      else
        a(i,j) = real ( j * ( n - i + 1 ), kind = 8 ) / real ( n + 1, kind = 8 )
      end if
    end do
  end do

  return
end
subroutine dismea ( m, n, x, p, nc_old, nc_new, e, d )

!*****************************************************************************80
!
!! DISMEA constructs a set of hierarchical clusters.
!
!  Discussion:
!
!    On input, the data is assumed to have been partitioned into
!    NC_OLD clusters.
!
!    The routine finds the cluster with greatest variance and
!    subdivides it into two smaller clusters.
!
!    This process is repeated until there are NC_NEW clusters.
!
!  Modified:
!
!    30 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 156-157,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input/output, integer ( kind = 4 ) P(M), contains the cluster assignments.
!
!    Input, integer ( kind = 4 ) NC_OLD, NC_NEW, the initial and final number of
!    clusters.  If the user sets NC_OLD to 0 on input, then the
!    code will take care of initializing P and E.  But if NC_OLD is
!    greater than 0, the user is responsible for setting meaningful
!    values for P and correct values of E(1:NC_OLD) and D.
!
!    Input/output, real ( kind = 8 ) E(NC_NEW), the cluster variances.
!
!    Input/output, real ( kind = 8 ) D, the total variance.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc_new

  real ( kind = 8 ) d
  real ( kind = 8 ) dc
  real ( kind = 8 ) e(nc_new)
  real ( kind = 8 ) ec(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) list(1)
  integer ( kind = 4 ) mc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nc_old
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) pc(m)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) xc(m,n)

  if ( m < nc_new ) then
    return
  end if

  do nc = nc_old+1, nc_new

    if ( nc == 1 ) then

      p(1:m) = 1

      call cluster_variance ( m, n, x, nc, p, e )
!
!  Cluster JC is to be split up.
!
!  Count MC, the number of items in cluster JC.
!  Create XC, a copy of the data in cluster JC.
!
    else
!
!  The semantics of the MAXLOC intrinsic get mighty tiresome when
!  you're doing the simplest case...
!
      list = maxloc ( e(1:nc-1) )
      jc = list(1)

      mc = 0

      do i = 1, m

        if ( p(i) == jc ) then
          mc = mc + 1
          xc(mc,1:n) = x(i,1:n)
        end if

      end do
!
!  Set the mini cluster assignment vector to alternate between 1 and 2.
!
      pc(1:mc) = mod ( pc(1:mc), 2 ) + 1
!
!  Have KMEANS split up cluster JC.
!
      call kmeans ( mc, n, xc, 2, pc, ec, dc )
!
!  Items in mini cluster 2 will be moved to the new cluster NC.
!
      mc = 0

      do i = 1, m
        if ( p(i) == jc ) then
          mc = mc + 1
          if ( pc(mc) == 2 ) then
            p(i) = nc
          end if
        end if
      end do

      d = d - e(jc) + dc

      e(jc) = ec(1)
      e(nc) = ec(2)

    end if

  end do

  return
end
subroutine divgow ( m, n, x, p )

!*****************************************************************************80
!
!! DIVGOW constructs a set of hierarchical clusters by doubling.
!
!  Discussion:
!
!    On input, the data is assumed to have been partitioned into
!    NC_OLD clusters.
!
!    The routine finds the cluster with greatest variance and
!    subdivides it into two smaller clusters.
!
!    This process is repeated until there are NC_NEW clusters.
!
!  Modified:
!
!    05 April 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 164-165,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input/output, integer ( kind = 4 ) P(M), contains the cluster assignments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) p1(m)
  integer ( kind = 4 ) p2(m)
  integer ( kind = 4 ) r
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) x1(m,n)
  real ( kind = 8 ) x2(m,n)
!
!  1 -> 2
!
  nd = 1
  call zweigo ( m, n, x, p )

  nc = 2
  jc = 2

  do

    do j = 2, nc, 2

      j1 = j + nc - 3
      j2 = j1 + 1
      m1 = 0
      m2 = 0

      do i = 1, m

        r = p(i)

        if ( r == j1 ) then
          m1 = m1 + 1
          x1(m1,1:n) = x(i,1:n)
        else if ( r == j2 ) then
          m2 = m2 + 1
          x2(m2,1:n) = x(i,1:n)
        end if

      end do

      if ( m1 <= 2 .or. m2 <= 2 ) then
        go to 99
      end if

      nd = nd + 1
      call zweigo ( m, n, x, p )

      nd = nd + 1
      call zweigo ( m, n, x, p )

      l1 = 1
      l2 = 1

      do i = 1, m

        r = p(i)

        if ( r == j1 ) then
          r = jc + p1(l1)
          l1 = l1 + 1
        else
          r = jc + 2 + p2(l2)
          l2 = l2 + 1
        end if

        p(i) = r

      end do

      jc = jc + 4

    end do

    nc = nc + nc

  end do

99 continue

  p(1:m) = p(1:m) - ( nd - 1 )

  return
end
subroutine emeans ( m, n, x, nc, p )

!*****************************************************************************80
!
!! EMEANS clusters data using a variant of the K-Means algorithm for L1 norms.
!
!  Discussion:
!
!    The data must already have been assigned to initial partitions.
!    This could be done randomly, by RANDP, or by JOINER or LEADER
!    or HMEANS any other way.
!
!    The K-Means algorithm tries to improve the initial partition
!    by a series of exchanges.  Every exchange is guaranteed to reduce
!    the total objective function which is based on the L1 norm.
!
!    The individual per-cluster objective functions are:
!
!      E(J) = sum ( P(I) = J ) sum ( 1 <= K <= N ) ( abs ( X(I,K) - S(J,K) ) )
!
!    and the total objective function is simply their sum:
!
!      D = sum ( 1 <= J <= NC ) E(J)
!
!    The algorithm presented here was incomplete until some missing
!    lines were restored on 12 May 2005.
!
!  Modified:
!
!    12 May 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, pages 149-153,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Input/output, integer ( kind = 4 ) P(M), the cluster assignments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  real ( kind = 8 ) e(nc)
  logical done
  logical even(nc)
  logical evenj
  logical evenr
  real ( kind = 8 ) f
  integer ( kind = 4 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) hh
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jr
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(nc)
  integer ( kind = 4 ) q1
  integer ( kind = 4 ) q2
  integer ( kind = 4 ) qj
  integer ( kind = 4 ) qq(nc+1)
  integer ( kind = 4 ) qqj
  integer ( kind = 4 ) qqr
  integer ( kind = 4 ) qr
  integer ( kind = 4 ) qw
  integer ( kind = 4 ) r
  logical revers
  integer ( kind = 4 ) sa(n)
  integer ( kind = 4 ) sb(n)
  integer ( kind = 4 ) sf(n)
  real ( kind = 8 ) sfk
  real ( kind = 8 ) ss
  integer ( kind = 4 ) t
  real ( kind = 8 ) u(nc,n)
  real ( kind = 8 ) ua(n)
  real ( kind = 8 ) ub(n)
  real ( kind = 8 ) uf(n)
  real ( kind = 8 ) uj
  real ( kind = 8 ) ur
  real ( kind = 8 ) v(nc,n)
  real ( kind = 8 ) va(n)
  real ( kind = 8 ) vb(n)
  real ( kind = 8 ) vf(n)
  real ( kind = 8 ) vj
  real ( kind = 8 ) vr
  real ( kind = 8 ) vv
  integer ( kind = 4 ) w
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) xx(m,n)
  real ( kind = 8 ) y(m)
  real ( kind = 8 ) yz
  integer ( kind = 4 ) z
  integer ( kind = 4 ) zq

  do i = 1, m
    if ( p(i) < 1 .or. nc < p(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EMEANS - Fatal error!'
      write ( *, '(a)' ) '  Illegal cluster assignment.'
      stop
    end if
  end do

  e(1:nc) = 0
!
!  Get the cluster populations Q.
!
  call cluster_population ( m, p, nc, q )

  d = 0.0D+00

  qq(1) = 1

  do j = 1, nc

    qj = q(j)

    if ( q(j) == 0 ) then
      return
    end if

    qqj = qq(j)
    qq(j+1) = qq(j) + q(j)
    q1 = qj / 2
    q2 = ( qj + 1 ) / 2
    even(j) = ( 2 * q1 == qj )

    do k = 1, n

      z = 0

      do i = 1, m
        if ( p(i) == j ) then
          z = z + 1
          y(z) = x(i,k)
        end if
      end do

      do w = 2, qj

        done = .true.
        qw = qj - w + 1

        do z = 1, qw

          if ( y(z+1) < y(z) ) then
            h = y(z+1)
            y(z+1) = y(z)
            y(z) = h
            done = .false.
          end if

        end do

        if ( done ) then
          exit
        end if

      end do

      hh = y(q2)

      u(j,k) = y(q2)

      if ( even(j) ) then
        v(j,k) = y(q2+1)
      else
        v(j,k) = y(q2)
      end if

      h = 0.0D+00
      do z = 1, qj
        yz = y(z)
        zq = qqj + z - 1
        xx(zq,k) = yz
        h = h + abs ( yz - hh )
      end do

      e(j) = e(j) + h
      d = d + h

    end do

  end do

  if ( nc == 1 ) then
    return
  end if
!
!  The exchange algorithm.
!
!  Delete object I from cluster R = P(I) if 1 < Q(R).
!
  i = 0
  it = 0

  do

    i = i + 1
    if ( m < i ) then
      i = i - m
    end if

    if ( it == m ) then
      exit
    end if

    r = p(i)
    qr = q(r)

    if ( q(r) <= 1 ) then
      cycle
    end if

    qqr = qq(r)
    a = 0.0D+00
    g = ( qr + 1 ) / 2 + qqr - 1
    evenr = even ( r )

    do k = 1, n

      ur = u(r,k)
      vr = v(r,k)
      h = x(i,k)

      do w = 1, qr

        t = qr - w + qqr
        if ( h == xx(t,k) ) then
          sa(k) = t
          exit
        end if

      end do

      if ( .not. evenr ) then

        ss = xx(g-1,k)
        vv = xx(g+1,k)

        if ( h == ur ) then
          ua(k) = ss
          va(k) = vv
        else if ( h < ur ) then
          ua(k) = ur
          va(k) = vv
          a = a + abs ( h - ur )
        else if ( ur < h ) then
          ua(k) = ss
          va(k) = ur
          a = a + abs ( h - ur )
        end if

      else

        if ( h <= ur ) then
          ua(k) = vr
          va(k) = vr
          a = a + abs ( h - vr )
        else if ( ur < h ) then
          ua(k) = ur
          va(k) = ur
          a = a + abs ( h - ur )
        end if

      end if

    end do
!
!  Try moving object I from cluster R into the clusters J = 1 through N,
!  and determine the cluster Z that gives the least increase.
!
    b = huge ( b )

    do j = 1, n

      if ( j /= r ) then
        f = 0.0D+00
        qj = q(j)
        done = ( 1 < qj )
        qqj = qq(j)
        g = ( qj + 1 ) / 2 + qqj - 1
        zq = qq(j+1)
        evenj = even(j)

        do k = 1, n

          uj = u(j,k)
          vj = v(j,k)
          h = x(i,k)
          sfk = zq

          do w = 1, qj
            t = qqj + w - 1
            if ( h <= xx(t,k) ) then
              sfk = t
              exit
            end if
          end do

          sf(k) = sfk

          if ( .not. evenj ) then

            ss = h

            if ( h <= uj ) then
              uf(k) = h
              if ( done ) then
                ss = xx(g-1,k)
              end if
              if ( h < ss ) then
                uf(k) = ss
              end if
              vf(k) = vj
              f = f + abs ( h - uj )
            else
              if ( done ) then
                ss = xx(g+1,k)
              end if
              uf(k) = uj
              vf(k) = h
              if ( ss < h ) then
                vf(k) = vj
              end if
              f = f + abs ( h - uj )
            end if

          else

            if ( vj <= h ) then
              uf(k) = vj
              vf(k) = vj
              f = f + abs ( h - vj )
            else if ( uj < h ) then
              uf(k) = h
              vf(k) = h
            else
              uf(k) = uj
              vf(k) = uj
              f = f + abs ( h - uj )
            end if

          end if

        end do

        if ( f <= b ) then
          b = f
          z = j
          sb(1:n) = sf(1:n)
          ub(1:n) = uf(1:n)
          vb(1:n) = vf(1:n)
        end if

      end if

    end do
!
!  If the objective function reduction for cluster R is greater than
!  the objective function increase for cluster Z, then carry out the
!  exchange.
!
    if ( a <= b ) then

      it = it + 1

    else

      e(r) = e(r) - a
      e(z) = e(z) + b
      d = d - a + b
      p(i) = z
      q(r) = q(r) - 1
      q(z) = q(z) + 1
      even(r) = .not. even(r)
      even(z) = .not. even(z)
      revers = ( z < r )

      u(r,1:n) = ua(1:n)
      v(r,1:n) = va(1:n)
      u(z,1:n) = ub(1:n)
      v(z,1:n) = vb(1:n)

      do k = 1, n

        jr = sa(k)
        jz = sb(k)

        if ( r < z ) then

          if ( jr <= jz - 2 ) then

            xx(jr:jz-2,k) = xx(jr+1:jz-1,k)
            xx(jz-1,k) = x(i,k)

          end if

        else

          xx(jz+1:jr,k) = xx(jz:jr-1,k)
          xx(jz,k) = x(i,k)

        end if

      end do

      if ( r < z ) then

        do w = r+1, z
          qq(w) = qq(w) - 1
        end do

      else

        do w = z+1, r
          qq(w) = qq(w) + 1
        end do

      end if

    end if

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
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
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
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
subroutine hiercl ( m, d, method, a, b, h )

!*****************************************************************************80
!
!! HIERCL implements seven agglomerative hierarchical clustering methods.
!
!  Modified:
!
!    04 May 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 180-181,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of objects to be clustered.
!
!    Input/output, real ( kind = 8 ) D(M,M), the distance matrix.
!    D(I,J) is the distance from I to J.  D should be symmetric, and
!    have a zero diagonal.  On output, D has been destroyed.
!
!    Input, integer ( kind = 4 ) METHOD, specifies the method to be used.  METHOD
!    should be between 1 and 7.
!
!    Output, integer ( kind = 4 ) A(M-1), B(M-1), contains the two indices of the
!    two clusters being merged at each step.
!
!    Output, real ( kind = 8 ) H(M-1), the distance between the two clusters
!    being merged at each step.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) a(m-1)
  integer ( kind = 4 ) b(m-1)
  real ( kind = 8 ) d(m,m)
  real ( kind = 8 ) dmin
  real ( kind = 8 ) h(m-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) method
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(m)

  if ( method < 1 .or. 7 < method ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HIERCL - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of METHOD.'
    stop
  end if

  k = 0
  p(1:m) = 0
  q(1:m) = 1

  do k = 1, m

    dmin = huge ( dmin )

    do i = 1, m - 1
      if ( p(i) == 0 ) then
        do j = i+1, m
          if ( p(j) == 0 ) then
            if ( d(i,j) <= dmin ) then
              ic = i
              jc = j
              dmin = d(i,j)
            end if
          end if
        end do
      end if
    end do

    p(jc) = 1
    a(k) = ic
    b(k) = jc
    h(k) = dmin

    do i = 1, m

      if ( i /= ic .and. p(i) == 0 ) then

        j = min ( ic, i )
        l = max ( ic, i )
        k1 = min ( jc, i )
        k2 = max ( jc, i )

        if ( method == 1 ) then

          d(j,l) = min ( d(j,l), d(k1,k2) )

        else if ( method == 2 ) then

          d(j,l) = max ( d(j,l), d(k1,k2) )

        else if ( method == 3 ) then

          d(j,l) = 0.5D+00 * ( d(j,l) + d(k1,k2) )

        else if ( method == 4 ) then

          d(j,l) = 0.5D+00 * ( d(j,l) + d(k1,k2) ) - 0.25D+00 * dmin

        else if ( method == 5 ) then

          d(j,l) = ( q(ic) * d(j,l) + q(jc) * d(k1,k2) ) &
            / real ( q(ic) + q(jc), kind = 8 )

        else if ( method == 6 ) then

          d(j,l) = ( q(ic) * d(j,l) + q(jc) * d(k1,k2) &
            - ( q(ic) * q(jc) * dmin / real ( q(ic) + q(jc), kind = 8 ) ) ) &
            / real ( q(ic) + q(jc), kind = 8 )

        else if ( method == 7 ) then

          d(j,l) = ( ( q(ic) + q(i) ) * d(j,l) + ( q(jc) + q(i) ) * d(k1,k2) &
            - q(i) * dmin ) / real ( q(i) + q(ic) + q(jc), kind = 8 )

        end if

      end if

    end do

    q(ic) = q(ic) + q(jc)

  end do

  return
end
subroutine hmeans ( m, n, x, nc, p )

!*****************************************************************************80
!
!! HMEANS clusters data using the H-Means algorithm.
!
!  Discussion:
!
!    The data must already have been assigned to initial partitions.
!    This could be done randomly, by RANDP, or by JOINER or LEADER
!    or any other way.
!
!    The H-Means algorithm tries to improve the initial partition
!    by a series of exchanges.  Every exchange is guaranteed to reduce
!    the total variance or energy of the clustering.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 65-67,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Input/output, integer ( kind = 4 ) P(M), the cluster assignments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  real ( kind = 8 ) d
  real ( kind = 8 ) dmax
  real ( kind = 8 ) e(nc)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(nc)
  integer ( kind = 4 ) r
  real ( kind = 8 ) s(nc,n)
  real ( kind = 8 ) x(m,n)

  id = 0

  dmax = huge ( dmax )

  do i = 1, m
    if ( p(i) < 1 .or. nc < p(i) ) then
      return
    end if
  end do

  do
!
!  Determine the cluster populations.
!
    call cluster_population ( m, p, nc, q )
!
!  Count the number of empty clusters.
!
    ir = 0
    do j = 1, nc
      if ( q(j) == 0 ) then
        ir = ir + 1
      end if
    end do
!
!  Determine the centroid of each cluster.
!
    s(1:nc,1:n) = 0.0D+00
    do i = 1, m
      r = p(i)
      s(r,1:n) = s(r,1:n) + x(i,1:n)
    end do

    do j = 1, nc
      s(j,1:n) = s(j,1:n) / real ( max ( q(j), 1 ), kind = 8 )
    end do
!
!  Determine the cluster energies or variances.
!
    e(1:nc) = 0.0D+00

    do i = 1, m
      r = p(i)
      e(r) = e(r) + sum ( ( s(r,1:n) - x(i,1:n) )**2 )
    end do
    d = sum ( e(1:nc) )
!
!  Return if there are any empty clusters.
!
    if ( ir /= 0 ) then
      exit
    end if
!
!  If the total energy has not decreased on this step, increment ID.
!  And once ID is 3, bail out.
!
    if ( dmax <= d ) then
      id = id + 1
    end if

    if ( 3 <= id ) then
      exit
    end if
!
!  For each point, find its nearest cluster center and move it to that cluster.
!
    dmax = d

    do i = 1, m

      f = huge ( f )
      r = 0

      do j = 1, nc

        g = sum ( ( s(j,1:n) - x(i,1:n) )**2 )

        if ( g < f ) then
          f = g
          r = j
        end if

      end do

      p(i) = r

    end do

  end do

  return
end
function i4_factorial ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL computes the factorial N!
!
!  Discussion:
!
!    FACTORIAL( N ) = PRODUCT ( 1 <= I <= N ) I
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, I4_FACTORIAL is returned as 1.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL, the factorial of N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) n

  i4_factorial = 1

  do i = 1, n
    i4_factorial = i4_factorial * i
  end do

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector A(I)=I.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_perml ( n, x, q, first )

!*****************************************************************************80
!
!! I4VEC_PERML generates permutations of an I4VEC in lexicographic order.
!
!  Modified:
!
!    05 May 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, pages 197-198,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the vector being permuted.
!
!    Input/output, integer ( kind = 4 ) X(N), the vector to be permuted.
!
!    Input/output, integer ( kind = 4 ) Q(N), information about the permutation.
!
!    Input/output, logical FIRST, should be set to TRUE on the first call.
!    Thereafter, the output value will be FALSE until all permutations
!    have been returned, at which point it will be set to TRUE.
!
  implicit none

  integer ( kind = 4 ) n

  logical first
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q(n)
  integer ( kind = 4 ) x(n)

  if ( first ) then
    first = .false.
    q(1:n-1) = n
  end if

  if ( q(n-1) == n ) then
    q(n-1) = n - 1
    call i4_swap ( x(n), x(n-1) )
    return
  end if

  k = 1
  first = .true.

  do kk = n-1, 1, -1
    if ( q(kk) /= kk ) then
      m = q(kk)
      call i4_swap ( x(m), x(kk) )
      q(kk) = m - 1
      k = kk + 1
      first = .false.
      exit
    end if
    q(kk) = n
  end do

  m = n

  do

    call i4_swap ( x(m), x(k) )
    m = m - 1
    k = k + 1

    if ( m <= k ) then
      exit
    end if

  end do


  return
end
subroutine i4vec_perms ( n, x, first )

!*****************************************************************************80
!
!! I4VEC_PERMS generates permutations of an I4VEC in lexicographic order.
!
!  Discussion:
!
!    The routine computes only the first N!/2 permutations, assuming
!    that the other N!/2 permutations can be produced by symmetry.
!
!    The use of the factorial function to produce a counter means that
!    the routine cannot handle values of N much greater than 13 or 14.
!
!  Modified:
!
!    06 May 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 198,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the vector being permuted.
!    N should be at least 2.
!
!    Input/output, integer ( kind = 4 ) X(N), the vector to be permuted.
!
!    Input/output, logical FIRST, should be set to TRUE on the first call.
!    Thereafter, the output value will be FALSE until the first N!/2
!    permutations have been returned, at which point it will be set to TRUE.
!
  implicit none

  integer ( kind = 4 ) n

  logical first
  logical first2
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ), save :: l = 0
  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: q
  integer ( kind = 4 ) x(n)

  if ( first ) then
    first = .false.
    first2 = .true.
    l = i4_factorial ( n ) / 2
    allocate ( q(1:n) )
  end if

  l = l - 1

  if ( l <= 0 ) then
    first = .true.
    deallocate ( q )
  else
    call i4vec_perml ( n, x, q, first2 )
  end if

  return
end
subroutine joiner ( m, n, x, rho, p, nc )

!*****************************************************************************80
!
!! JOINER uses a very simple cluster assignment algorithm.
!
!  Discussion:
!
!    JOINER implements an ad hoc construction of clusters without any
!    special optimal chararacteristics.
!
!  Algorithm:
!
!    1: NC = 0
!
!    2: If all X's have been assigned, return.
!
!    3: NC = NC + 1
!
!    4: Find an element X(I) which has not been assigned to a cluster,
!       and whose distance to all other unassigned elements is maximum.
!       This element becomes the first element of a new cluster,
!       whose centroid is initially equal to X(I).
!
!    5: Any unassigned elements which are closer to the centroid of the
!       NC-th cluster than to any other, and for which this distance
!       is less than RHO, are now assigned to cluster NC.  Each element
!       added to the cluster requires the recalculation of the centroid.
!
!    6: Go to 2.
!
!  Modified:
!
!    24 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 46-48,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of X.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, real ( kind = 8 ) RHO, a clustering tolerance.  To join a cluster,
!    a new point has to be within RHO of the representative.
!    RHO should not be negative.
!
!    Output, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Output, integer ( kind = 4 ) P(M), the cluster assignments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(m)
  integer ( kind = 4 ) r
  real ( kind = 8 ) rho
  real ( kind = 8 ) s(m,n)
  integer ( kind = 4 ) u
  real ( kind = 8 ) x(m,n)

  nc = 0

  if ( m <= 0 ) then
    return
  end if

  p(1:m) = 0

  if ( rho < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JOINER - Fatal error!'
    write ( *, '(a)' ) '  RHO must be nonnegative.'
    stop
  end if

  do

    r = 0
    u = 0
    h = 0.0D+00
!
!  Determine H, the largest distance between two unassigned points.
!
    do i = 1, m-1
      if ( p(i) == 0 ) then
        u = i
        do j = i+1, m
          if ( p(j) == 0 ) then
            f = sum ( ( x(i,1:n) - x(j,1:n) )**2 )
            if ( h < f ) then
              h = f
              r = i
            end if
          end if
        end do
      end if
    end do
!
!  If we didn't find any pair of points to consider, then there
!  are either 1 or 0 points left.  Take care of the case of 1 point
!  left, and then exit.
!
    if ( r == 0 ) then
      if ( u /= 0 ) then
        nc = nc + 1
        p(u) = nc
      end if
      exit
    end if
!
!  The cluster begins with a single point, R.
!
    nc = nc + 1
    p(r) = nc
    q(nc) = 1
    s(nc,1:n) = x(r,1:n)

    do

      h = huge ( h )
!
!  Search for an unassigned point, and record its distance
!  to the centroid of cluster NC.
!
      do i = 1, m
        if ( p(i) == 0 ) then
          f = sum ( ( s(nc,1:n) - x(i,1:n) )**2 )
          if ( f < h ) then
            h = f
            r = i
          end if
        end if
      end do
!
!  If no unassigned points are near enough to the centroid, exit
!  the loop.
!
      if ( rho < sqrt ( h ) ) then
        exit
      end if
!
!  Otherwise, assign point R to cluster NC, adjust the centroid,
!  and go back to search for more points to add.
!
      p(r) = nc
      q(nc) = q(nc) + 1

      s(nc,1:n) = ( real ( q(nc) - 1, kind = 8 ) * s(nc,1:n) + x(r,1:n) ) &
        / real ( q(nc), kind = 8 )

    end do

  end do

  return
end
subroutine kmeans ( m, n, x, nc, p, e, d )

!*****************************************************************************80
!
!! KMEANS clusters data using the K-Means algorithm.
!
!  Discussion:
!
!    The data must already have been assigned to initial partitions.
!    This could be done randomly, by RANDP, or by JOINER or LEADER
!    or HMEANS any other way.
!
!    The K-Means algorithm tries to improve the initial partition
!    by a series of exchanges.  Every exchange is guaranteed to reduce
!    the total variance or energy of the clustering.
!
!  Modified:
!
!    26 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 72-74,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Input/output, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Output, real ( kind = 8 ) E(NC), the cluster variances.
!
!    Output, real ( kind = 8 ) D, the total variance.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  real ( kind = 8 ) e(nc)
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(nc)
  integer ( kind = 4 ) r
  real ( kind = 8 ) s(nc,n)
  integer ( kind = 4 ) v
  real ( kind = 8 ) x(m,n)

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS - Fatal error!'
    write ( *, '(a)' ) '  M <= 0.'
    stop
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS - Fatal error!'
    write ( *, '(a)' ) '  N <= 0.'
    stop
  end if

  if ( nc <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS - Fatal error!'
    write ( *, '(a)' ) '  NC <= 0.'
    stop
  end if

  if ( m < nc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS - Fatal error!'
    write ( *, '(a)' ) '  M < NC.'
    stop
  end if
!
!  Make sure the cluster assignments are legal.
!
  do i = 1, m
    if ( p(i) < 1 .or. nc < p(i) ) then
      return
    end if
  end do
!
!  If there's just one cluster, we're done.
!
  if ( nc == 1 ) then
    return
  end if
!
!  Determine the cluster populations.
!
  call cluster_population ( m, p, nc, q )
!
!  Count the number of empty clusters.
!
  ir = 0
  do j = 1, nc
    if ( q(j) == 0 ) then
      ir = ir + 1
    end if
  end do
!
!  Determine the centroid of each cluster.
!
  s(1:nc,1:n) = 0.0D+00
  do i = 1, m
    r = p(i)
    s(r,1:n) = s(r,1:n) + x(i,1:n)
  end do

  do j = 1, nc
    s(j,1:n) = s(j,1:n) / real ( max ( q(j), 1 ), kind = 8 )
  end do
!
!  Determine the cluster energies or variances.
!
  e(1:nc) = 0.0D+00

  do i = 1, m
    r = p(i)
    e(r) = e(r) + sum ( ( s(r,1:n) - x(i,1:n) )**2 )
  end do
  d = sum ( e(1:nc) )
!
!  Initialize the loop
!
  i = 0
  it = 0

  do

    i = i + 1

    if ( m < i ) then
      i = i - m
    end if
!
!  If we have examined every point without an exchange, there is
!  nothing more we can do.
!
    if ( it == m ) then
      exit
    end if

    r = p(i)

    if ( q(r) <= 1 ) then
      cycle
    end if

    a = real ( q(r), kind = 8 ) &
      * sum ( ( s(r,1:n) - x(i,1:n) )**2 ) / real ( q(r) - 1, kind = 8 )
    b = huge ( b )

    do j = 1, nc

      if ( j /= r ) then

        f = real ( q(j), kind = 8 ) * sum ( ( s(j,1:n) - x(i,1:n) )**2 ) &
          / real ( q(j) + 1, kind = 8 )

        if ( f <= b ) then
          b = f
          v = j
        end if

      end if

    end do

    if ( a <= b ) then

      it = it + 1

    else

      it = 0
      e(r) = e(r) - a
      e(v) = e(v) + b
      d = d - a + b

      s(r,1:n) = ( real ( q(r), kind = 8 ) * s(r,1:n) - x(i,1:n) ) &
        / real ( q(r) - 1, kind = 8 )
      s(v,1:n) = ( real ( q(v), kind = 8 ) * s(v,1:n) + x(i,1:n) ) &
        / real ( q(v) + 1, kind = 8 )

      p(i) = v
      q(r) = q(r) - 1
      q(v) = q(v) + 1

    end if

  end do

  return
end
subroutine leader ( m, n, x, rho, p, nc )

!*****************************************************************************80
!
!! LEADER uses a very simple cluster assignment algorithm.
!
!  Discussion:
!
!    A clustering tolerance RHO is specified.
!
!    Initially, a single cluster is specified, containing the first point.
!
!    At every stage, each cluster is represented by the first point
!    that was assigned to it.
!
!    When it is time to consider a new point, it is added to the first
!    cluster whose representative is closer than RHO.
!
!    If no such point can be found, a new cluster is formed, and the
!    point becomes its representative.
!
!  Modified:
!
!    23 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 38,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of X.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, real ( kind = 8 ) RHO, a clustering tolerance.  To join a cluster,
!    a new point has to be within RHO of the representative.
!    RHO should not be negative.
!
!    Output, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Output, integer ( kind = 4 ) P(M), the cluster assignments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(m)
  real ( kind = 8 ) rho
  real ( kind = 8 ) x(m,n)

  nc = 0
  p(1:m) = 0

  if ( m <= 0 ) then
    return
  end if

  if ( rho < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEADER - Fatal error!'
    write ( *, '(a)' ) '  RHO must be nonnegative.'
    stop
  end if

  nc = 1
  f(1) = 1

  do i = 1, m

    do j = 1, nc

      if ( f(j) < 1 .or. m < f(j) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEADER - Fatal error!'
        write ( *, '(a,i6)' ) '  I = ', i
        write ( *, '(a,i6)' ) '  J = ', j
        write ( *, '(a,g14.6)' ) '  F(J) = ', f(j)
        write ( *, '(a,i6)' ) '  M = ', m
        stop
      end if

      d = sqrt ( sum ( ( x(f(j),1:n) - x(i,1:n) )**2 ) )

      if ( d <= rho ) then
        p(i) = j
        exit
      end if

    end do
!
!  If point I was not within RHO of any representative,
!  put it in its own cluster.
!
    if ( p(i) == 0 ) then
      nc = nc + 1
      p(i) = nc
      f(nc) = i
    end if

  end do

  return
end
subroutine linker ( m, d, q, p, t )

!*****************************************************************************80
!
!! LINKER contructs a minimal tree for a symmetric distance matrix.
!
!  Discussion:
!
!    For each I from 1 to M-1, a partner point P(I) is sought between
!    1 and M, with P(I) distinct from I, so that P(I) is at minimum
!    distance T(I), so that all points are simply connected.  The
!    resulting connected points form a minimal spanning tree.
!
!    The triples ( I, P(I), T(I) ) are ordered according to the
!    magnitude of T(I), and the ordered triples are stored in
!    Q(I), P(I), T(I).
!
!    The two possible hierarchical cluster processes are identical
!    for this case.  That is, once the minimal spanning tree is
!    set up, a hierarchy of clusters can be formed by repeated splitting
!    where T(I) is largest, or merging where T(I) is smallest.
!
!  Modified:
!
!    03 May 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 173-174,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the matrix.
!
!    Input, real ( kind = 8 ) D(M,M), the distance matrix.  D(I,J) is the
!    distance from I to J.  D should be symmetric, and have a zero diagonal.
!
!    Output, integer ( kind = 4 ) Q(M-1), P(M-1), T(M-1), the index of a point,
!    its nearest neighbor, and the distance between them.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) d(m,m)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p(1:m-1)
  integer ( kind = 4 ) q(1:m-1)
  real ( kind = 8 ) t(1:m-1)
  real ( kind = 8 ) u

  if ( m <= 1 ) then
    return
  end if

  q(1:m-1) = 0
  p(1:m-1) = 0
  t(1:m-1) = huge ( t(1:m-1) )

  j = m

  do i = 1, m-1

    u = huge ( u )

    do k = 1, m-1

      if ( q(k) == 0 ) then

        if ( d(j,k) < t(k) ) then

          t(k) = d(j,k)
          p(k) = j

          if ( t(k) < u ) then
            u = t(k)
            n = k
          end if

        end if

      end if

    end do

    j = n
    q(j) = 1

  end do

  call i4vec_indicator ( m-1, q )

  done = .false.

  do i = 2, m-1

    done = .true.

    do j = 1, m - i

      if ( t(j+1) < t(j) ) then
        call r8_swap ( t(j+1), t(j) )
        call i4_swap ( q(j+1), q(j) )
        call i4_swap ( p(j+1), p(j) )
        done = .false.
      end if

    end do

    if ( done ) then
      exit
    end if

  end do

  return
end
subroutine ordered ( m, x, nc, q, s )

!*****************************************************************************80
!
!! ORDERED clusters one-dimensional ordered data into NC clusters.
!
!  Discussion:
!
!    The input data must be sorted in ascending order.
!
!    A dynamic programming algorithm is used.
!
!    At the moment, I don't believe I am correctly describing the
!    output contents of S.
!
!  Modified:
!
!    29 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 63,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, real ( kind = 8 ) X(M), the data to be clustered.  The data must be
!    sorted in ascending order.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters to create.
!
!    Output, integer ( kind = 4 ) Q(NC,NC), describes the clusters of data.
!    The last row of Q stores the first element in each cluster.
!    Thus cluster 1 contains data items Q(NC,1) through Q(NC,2)-1.
!    Other information is contained in previous rows.  In particular,
!    in row J, columns 1 through J, there is similar information
!    about a partition involving J clusters.
!
!    Output, real ( kind = 8 ) S(M,NC), contains pointwise variances for a
!    number of clusters between 1 and NC.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nc

  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) il
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q(nc,nc)
  integer ( kind = 4 ) r(m,nc)
  logical r8vec_ascends
  real ( kind = 8 ) s(m,nc)
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x(m)
!
!  Verify that the data are ordered.
!
  if ( .not. r8vec_ascends ( m, x ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ORDERED - Fatal error!'
    write ( *, '(a)' ) '  Data is not in ascending order.'
    stop
  end if

  r(1,1:nc) = 1
  r(2:m,1:nc) = 0

  s(1,1:nc) = 0.0D+00
  s(2:m,1:nc) = huge ( 1.0D+00 )

  if ( nc <= 1 ) then
    return
  end if

  do i = 2, m

    t = 0.0D+00
    u = 0.0D+00

    do l = i, 1, -1

      f = x(l)
      t = t + f
      u = u + f**2
      v = u - t**2 / real ( i - l + 1, kind = 8 )
      p = l - 1

      if ( p /= 0 ) then

        do j = 2, nc
          f = s(p,j-1) + v
          if ( f <= s(i,j) ) then
            r(i,j) = l
            s(i,j) = f
          end if
        end do

      end if

    end do

    s(i,1) = v
    r(i,1) = 1

  end do

  do k = nc, 1, -1

    il = m + 1

    do ll = k, 1, -1
      iu = il - 1
      il = r(iu,ll)
      q(k,ll) = il
    end do

  end do

  return
end
subroutine profile ( m, n, a, b, trans, z, ind )

!*****************************************************************************80
!
!! PROFILE seeks an optimal variable ordering for a set of data.
!
!  Discussion:
!
!    To understand what is going on here, suppose we have N objects,
!    each of which is an M vector, and that on one sheet of paper,
!    we plot each of the objects as a function of its vector indices.
!    That is, if object 1 is ( 5.3, 19.6, 34.2), then the corresponding
!    broken line graph is (1, 5.3), (2, 19.6), (3, 34.2).
!
!    Now, having plotted all the items, we have N line graphs.  Any pair
!    of line graphs will cross each other from index I to index I+1 if
!    the relative ordering of their I-th and I+1 values reverses.
!
!    We may influence the number of crossings by reordering the indices.
!    This routine seeks to find the optimal ordering of the indices
!    which produces the minimal number of such crossings.
!
!  Modified:
!
!    06 May 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 199-200,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the objects.
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, real ( kind = 8 ) A(M,N), the objects as columns of a matrix.
!
!    Output, integer ( kind = 4 ) B(N,N), a distance matrix for the objects.
!
!    Input, logical TRANS, if TRUE, specifies that before any other
!    operations, each column of A is to be transformed so that
!    the minimum value is 0 and the maximum value is Z.
!
!    Input, real ( kind = 8 ) Z, is used if TRANS is TRUE on input, and
!    represents the upper bound for entries of A.
!
!    Output, integer ( kind = 4 ) IND(N), the optimal ordering of the variables.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) b(n,n)
  logical first
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  logical trans
  integer ( kind = 4 ) u(n)
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) z

  if ( trans ) then

    do j = 1, n

      v = maxval ( a(1:m,j) )
      w = minval ( a(1:m,j) )

      if ( v == w ) then
        a(1:m,j) = 0.0D+00
      else
        a(1:m,j) = z * ( a(1:m,j) - w ) / ( v - w )
      end if

    end do

  end if

  b(n,n) = 0

  do k = 1, n - 1
    b(k,k) = 0
    do i = k+1, n
      h = 0
      do l = 1, m-1
        v = a(l,i)
        w = a(l,k)
        do j = l + 1, m
          if ( ( v - a(j,i) ) * ( w - a(j,k) ) < 0.0D+00 ) then
            h = h + 1
          end if
        end do
      end do
      b(i,k) = h
      b(k,i) = h
    end do
  end do

  first = .true.

  call i4vec_indicator ( n, u )

  g = huge ( g )

  do

    h = 0
    j = u(1)

    do k = 2, n
      i = u(k)
      h = h + b(j,i)
      j = i
    end do

    if ( h <= g ) then
      g = h
      ind(1:n) = u(1:n)
    end if

    call i4vec_perms ( n, u, first )

    if ( first ) then
      exit
    end if

  end do

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8mat_det ( n, a, det )

!*****************************************************************************80
!
!! R8MAT_DET computes the determinant of an R8MAT.
!
!  Modified:
!
!    27 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 125-127,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) piv(1)

  b(1:n,1:n) = a(1:n,1:n)

  det = 1.0D+00

  do k = 1, n

    piv = maxloc ( abs ( b(k:n,k) ) )

    m = piv(1) + k - 1

    if ( m /= k ) then
      det = - det
      call r8_swap ( b(m,k), b(k,k) )
    end if

    det = det * b(k,k)

    if ( b(k,k) /= 0.0D+00 ) then

      b(k+1:n,k) = - b(k+1:n,k) / b(k,k)

      do j = k+1, n
        if ( m /= k ) then
          call r8_swap ( b(m,j), b(k,j) )
        end if
        b(k+1:n,j) = b(k+1:n,j) + b(k+1:n,k) * b(k,j)
      end do

    end if

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Modified:
!
!    23 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 5
    jhi = min ( jlo + 4, n )
    write ( *, '(a)' ) ' '
    write ( *, '(8x,5(i7,7x))' ) (j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(2x,i6,5g14.6)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
function r8vec_ascends ( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS determines if an R8VEC is (weakly) ascending.
!
!  Example:
!
!    X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!
!    R8VEC_ASCENDS = TRUE
!
!    The sequence is not required to be strictly ascending.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS, is TRUE if the entries of X ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical r8vec_ascends
  real ( kind = 8 ) x(n)

  r8vec_ascends = .false.

  do i = 1, n-1
    if ( x(i+1) < x(i) ) then
      return
    end if
  end do

  r8vec_ascends = .true.

  return
end
subroutine r8vec_sort_bubble_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_BUBBLE_A ascending bubble sorts an R8VEC.
!
!  Discussion:
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n-1
    do j = i+1, n
      if ( a(j) < a(i) ) then
        call r8_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine randp ( m, m0, nc, p, seed )

!*****************************************************************************80
!
!! RANDP randomly partitions a set of M items into N clusters.
!
!  Discussion:
!
!    The code has been modified from the printed version so that the
!    user can specify that all clusters must have at least M0 elements.
!
!  Modified:
!
!    04 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 143,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of items to assign.
!
!    Input, integer ( kind = 4 ) M0, the minimum number of items in each cluster.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters.
!
!    Output, integer ( kind = 4 ) P(M), the cluster to which each item is assigned.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed used by the random
!    number generator.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) urand
!
!  Use the first NC * M0 data items to guarantee that each cluster has
!  at least M0 elements.
!
  if ( m < nc * m0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDP - Fatal error!'
    write ( *, '(a)' ) '  Not enough data to satisfy occupancy requirements.'
    write ( *, '(a)' ) '  Require NC * M0 <= M.'
    write ( *, '(a,i6)' ) '  but NC = ', nc
    write ( *, '(a,i6)' ) '  M0 =    ', m0
    write ( *, '(a,i6)' ) '  M =     ', m
    stop
  end if

  k = 0
  do i = 1, nc
    do j = 1, m0
      k = k + 1
      p(k) = i
    end do
  end do
!
!  Now take care of the remaining data items.
!
  do i = k+1, m

    j = int ( real ( nc, kind = 8 ) * urand ( seed ) ) + 1
    j = min ( j, nc )
    j = max ( j, 1 )

    p(i) = j

  end do

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. nchar <= lchar + 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
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
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine standn ( m, n, x, w, s, eps, itmax, is, f )

!*****************************************************************************80
!
!! STANDN solves the single location problem in N dimensions.
!
!  Discussion:
!
!    The algorithm attempts to determine the position of a point X*
!    so as to minimize the objective function
!
!      F = sum ( 1 <= I <= M ) W(I) * dist ( X(I), X* )
!
!    where dist ( X, Y ) is the usual Euclidean distance.
!
!  Modified:
!
!    28 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 134-136,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, real ( kind = 8 ) W(M), the weights associated with each data point.
!
!    Input/output, real ( kind = 8 ) S(N).  On input, if IS is nonzero, then
!    S is the initial estimate of the location of X*.  On output,
!    S contains the program's estimate of the location.
!
!    Input, real ( kind = 8 ) EPS, a tolerance used in an accuracy test.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations allowed.
!
!    Input, integer ( kind = 4 ) IS, is 0 if the initial guess for the location
!    of X* should be made by the program, or nonzero if the user
!    has supplied a guess for the location in S.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) k
  real ( kind = 8 ) p
  logical ptrue
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) v
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  it = 0
  f = 0.0D+00
!
!  If the user did not supply an initial estimate for the solution S,
!  set S to the weighted centroid.
!
  if ( is == 0 ) then

    s(1:n) = 0.0D+00
    do i = 1, m
      s(1:n) = s(1:n) + w(i) * x(i,1:n)
    end do
    s(1:n) = s(1:n) / sum ( w(1:m) )

    if ( m == 1 ) then
      return
    end if

  end if

  do

    it = it + 1

    if ( itmax < it ) then
      return
    end if

    t(1:n) = 0.0D+00
    f = 0.0D+00
    z = 0.0D+00

    do i = 1, m

      v = w(i)
      p = sum ( ( s(1:n) - x(i,1:n) )**2 )
      ptrue = p < 1.0D-10

      if ( ptrue ) then
        cycle
      end if

      p = sqrt ( p )
      f = f + v * p
      p = v / p
      z = z + p
      t(1:n) = t(1:n) + p * x(i,1:n)

    end do

    if ( ptrue ) then
      is = -1
      return
    end if

    p = 0.0D+00
    v = 0.0D+00
    z = 1.0D+00 / z

    do k = 1, n
      y = t(k) * z
      v = v + abs ( y )
      p = p + abs ( y - s(k) )
      s(k) = y
    end do

    if ( p < eps * v ) then
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
!    May 31 2001   9:45:54.872 AM
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
subroutine transf ( m, l, x )

!*****************************************************************************80
!
!! TRANSF transforms a data set to have zero mean and unit variance.
!
!  Discussion:
!
!    Each of the columns of X is transformed to have mean zero and
!    variance 1.
!
!  Modified:
!
!    21 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 21,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) L, the number of columns of X.
!
!    Input/output, real ( kind = 8 ) X(M,L), the data to be transformed.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  integer ( kind = 4 ) k
  real ( kind = 8 ) q
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) x(m,l)
!
  do k = 1, l

    t = sum ( x(1:m,k) )
    u = sum ( x(1:m,k)**2 )

    q = t / real ( m, kind = 8 )

    s = sqrt ( real ( m - 1, kind = 8 ) / ( u - t * q ) )

    x(1:m,k) = s * ( x(1:m,k) - q )

  end do

  return
end
function urand ( seed )

!*****************************************************************************80
!
!! URAND returns a pseudo-random number uniformly distributed in [0,1].
!
!  Modified:
!
!    12 January 2003
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 143,
!    QA278 S68213.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) URAND, the pseudo-random number.
!
  implicit none

  double precision halfm
  integer ( kind = 4 ), save :: ia
  integer ( kind = 4 ), save :: ic
  integer ( kind = 4 ), parameter :: itwo = 2
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: m2 = 0
  integer ( kind = 4 ), save :: mic
  real ( kind = 8 ), save :: s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) urand

  if ( m2 == 0 ) then

    m = 1

    do

      m2 = m
      m = itwo * m

      if ( m <= m2 ) then
        exit
      end if

    end do

    halfm = dble ( m2 )
    ia = 5 + 8 * int ( halfm * atan ( 1.0D+0 ) / 8.0D+0 )
    ic = 1 + 2 * int ( halfm * ( 0.5D0 - sqrt ( 3.0D+00 ) / 6.0D+00 ) )
    mic = ( m2 - ic ) + m2
    s = 0.5D+00 / halfm

  end if

  seed = seed * ia

  if ( mic < seed ) then
    seed = ( seed - m2 ) - m2
  end if

  seed = seed + ic
  if ( m2 < seed / 2 ) then
    seed = ( seed - m2 ) - m2
  end if

  if ( seed < 0 ) then
    seed = ( seed + m2 ) + m2
  end if

  urand = real ( seed, kind = 8 ) * s

  return
end
subroutine wmeans ( m, n, x, nc, p, det )

!*****************************************************************************80
!
!! WMEANS clusters data using the determinant criterion.
!
!  Discussion:
!
!    The data must already have been assigned to initial partitions.
!    This could be done randomly, by RANDP, or by JOINER or LEADER
!    or HMEANS any other way.
!
!    The W-Means algorithm tries to improve the initial partition
!    by a series of exchanges.  Every exchange is guaranteed to reduce
!    the determinant of the sum, formed over the clusters, of the
!    dyadic products of the differences of the cluster elements
!    with their centroid.
!
!  Modified:
!
!    27 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 125-127,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of data.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Input, integer ( kind = 4 ) NC, the number of clusters created.
!
!    Input/output, integer ( kind = 4 ) P(M), the cluster assignments.
!
!    Output, real ( kind = 8 ) DET, the determinant.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc

  real ( kind = 8 ) am(n,n)
  real ( kind = 8 ) bm(n,n)
  real ( kind = 8 ) det
  real ( kind = 8 ) detn
  real ( kind = 8 ) detv
  real ( kind = 8 ) dm(n,n)
  real ( kind = 8 ) dn(n,n)
  real ( kind = 8 ) dv(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(nc)
  integer ( kind = 4 ) r
  real ( kind = 8 ) s(nc,n)
  integer ( kind = 4 ) v
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)
!
!  Make sure the cluster assignments are legal.
!
  do i = 1, m
    if ( p(i) < 1 .or. nc < p(i) ) then
      return
    end if
  end do
!
!  If there's just one cluster, we're done.
!
  if ( nc == 1 ) then
    return
  end if
!
!  Determine the cluster populations.
!
  call cluster_population ( m, p, nc, q )
!
!  Count the number of empty clusters.
!
  ir = 0
  do j = 1, nc
    if ( q(j) == 0 ) then
      ir = ir + 1
    end if
  end do
!
!  Determine the centroid of each cluster.
!
  s(1:nc,1:n) = 0.0D+00
  do i = 1, m
    r = p(i)
    s(r,1:n) = s(r,1:n) + x(i,1:n)
  end do

  do j = 1, nc
    s(j,1:n) = s(j,1:n) / real ( max ( q(j), 1 ), kind = 8 )
  end do

  dm(1:n,1:n) = 0.0D+00

  do i = 1, m

    r = p(i)
    z(1:n) = x(i,1:n) - s(r,1:n)

    do ii = 1, n
      dm(ii,1:n) = z(ii) * z(1:n)
    end do

  end do

  call r8mat_det ( n, dm, det )

  i = 0
  it = 0

  do

    i = i + 1
    if ( m < i ) then
      i = i - m
    end if

    if ( it == m ) then
      return
    end if

    r = p(i)

    if ( q(r) <= 1 ) then
      cycle
    end if

    z(1:n) = x(i,1:n) - s(r,1:n)

    do ii = 1, n
      am(ii,1:n) = z(ii) * z(1:n)
    end do

    detv = huge ( detv )

    do j = 1, nc

      if ( r /= j ) then

        z(1:n) = x(i,1:n) - s(j,1:n)

        do ii = 1, n
          bm(ii,1:n) = z(ii) * z(1:n)
        end do

        dn(1:n,1:n) = dm(1:n,1:n) &
          - real ( q(r), kind = 8 ) * am(1:n,1:n) &
          / real ( q(r) - 1, kind = 8 ) &
          + real ( q(j), kind = 8 ) * bm(1:n,1:n) &
          / real ( q(j) + 1, kind = 8 )

        call r8mat_det ( n, dn, detn )

        if ( detn <= detv ) then
          detv = detn
          dv(1:n,1:n) = dn(1:n,1:n)
          v = j
        end if

      end if

    end do

    if ( det <= detv ) then

      it = it + 1

    else

      it = 0
      det = detv

      dm(1:n,1:n) = dv(1:n,1:n)

      s(r,1:n) = ( real ( q(r), kind = 8 ) * s(r,1:n) - x(i,1:n) ) &
        / real ( q(r) - 1, kind = 8 )
      s(v,1:n) = ( real ( q(v), kind = 8 ) * s(v,1:n) + x(i,1:n) ) &
        / real ( q(v) + 1, kind = 8 )

      p(i) = v
      q(r) = q(r) - 1
      q(v) = q(v) + 1

    end if

  end do

  return
end
subroutine zweigo ( m, n, x, p )

!*****************************************************************************80
!
!! ZWEIGO organizes a set of data into two clusters.
!
!  Algorithm:
!
!    1. Take two point whose distance is maximum, and put one in
!       each cluster.
!
!    2. If there are no unassigned points, stop.
!
!    3. For each unassigned point, compute its distance to the
!       centroids of clusters 1 and 2.  X1 will be the unassigned
!       point that is nearest cluster 1, and X2 the unassigned point
!       that is nearest cluster 2.  (These could actually be the
!       same point.)
!
!    4. Assign X1 to cluster 1, or X2 to cluster 2, depending on
!       which is closest.  Update the cluster population and
!       centroid value.
!
!    5. Go back to step 2.
!
!  Modified:
!
!    24 April 2002
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Edwards and Cavalli-Sforza,
!    A Method for Cluster Analysis,
!    Biometrics, Volume 21, 1965, pages 262-275.
!
!    Helmuth Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 53-55,
!    QA278 S6813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of X.
!
!    Input, integer ( kind = 4 ) N, the number of columns of X.
!
!    Input, real ( kind = 8 ) X(M,N), the data to be clustered.
!
!    Output, integer ( kind = 4 ) P(M), the cluster assignments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dmax
  integer ( kind = 4 ) g1
  integer ( kind = 4 ) g2
  real ( kind = 8 ) h
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) p(m)
  real ( kind = 8 ) s1(n)
  real ( kind = 8 ) s2(n)
  real ( kind = 8 ) x(m,n)

  if ( m < 1 ) then
    return
  end if

  if ( m == 1 ) then
    p(1) = 1
    return
  end if

  if ( m == 2 ) then

    p(1) = 1

    if ( sqrt ( sum ( ( x(1,1:n) - x(2,1:n) )**2 ) ) == 0.0D+00 ) then
      p(2) = 1
    else
      p(2) = 2
    end if

    return

  end if

  p(1:m) = 0
!
!  Find the two furthest apart points and use these to start the
!  two clusters.
!
  dmax = - huge ( dmax )

  do i = 1, m-1
    do j = i+1, m
      h = sum ( ( x(i,1:n) - x(j,1:n) )**2 )
      if ( dmax < h ) then
        dmax = h
        ic = i
        jc = j
      end if
    end do
  end do

  s1(1:n) = x(ic,1:n)
  s2(1:n) = x(jc,1:n)
  p(ic) = 1
  p(jc) = 2
  g1 = 1
  g2 = 2
!
!  Now find the point which is unassigned, and closest to one of the
!  two centroids.  Add it to the corresponding cluster, and update
!  the cluster centroid and population.
!
  ip = 2

  do

    h = sqrt ( sum ( ( s1(1:n) - s2(1:n) )**2 ) )

    if ( ip == m ) then
      exit
    end if

    ip = ip + 1

    d1 = huge ( d1 )
    d2 = huge ( d2 )

    do i = 1, m

      if ( p(i) == 0 ) then

        h1 = sum ( ( s1(1:n) - x(i,1:n) )**2 )
        h2 = sum ( ( s2(1:n) - x(i,1:n) )**2 )

        if ( h1 < d1 ) then
          d1 = h1
          i1 = i
        end if

        if ( h2 < d2 ) then
          d2 = h2
          i2 = i
        end if

      end if

    end do
!
!  Assign point I1 to cluster 1...
!
    if ( d1 < d2 ) then
      p(i1) = 1
      s1(1:n) = ( real ( g1, kind = 8 ) * s1(1:n) + x(i1,1:n) ) &
        / real ( g1 + 1, kind = 8 )
      g1 = g1 + 1
!
!  ...or assign point I2 to cluster 2.
!
    else
      p(i2) = 2
      s2(1:n) = ( real ( g2, kind = 8 ) * s2(1:n) + x(i2,1:n) ) &
        / real ( g2 + 1, kind = 8 )
      g2 = g2 + 1
    end if

  end do

  return
end
