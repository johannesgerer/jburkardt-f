subroutine bvpexm ( x, m, s, z, mj, n, e, d, it, iflag )

!*****************************************************************************80
!
!! BVPEXM implements the exchange algorithm on binary data for the L1 criterion.
!
!  Discussion:
!
!    The routine requires that the input data be binary data.  That is,
!    each data item X(*,1:S) comprises S integers between 0 and 1.
!
!    The L1 criterion measures the L1 distance of each point from its
!    cluster median, and sums this over all clusters.
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
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, pages 129-130,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    is assigned.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, integer ( kind = 4 ) E(N), the per-cluster value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) D, the total value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations that were taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  integer ( kind = 4 ) a(n,s)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e(n)
  integer ( kind = 4 ) ej
  integer ( kind = 4 ) ep
  integer ( kind = 4 ) eq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) y
  integer ( kind = 4 ) x(m,s)
  integer ( kind = 4 ) z(m)

  do i = 1, m

    j = z(i)

    if ( j < 1 .or. n < j ) then
      iflag = 1
      return
    end if

  end do

  it = 0
  iflag = 0
  m0 = 1
!
!  Determine the cluster populations.
!
  mj(1:n) = 0
  do i = 1, m
    j = z(i)
    mj(j) = mj(j) + 1
  end do

  a(1:n,1:s) = 0
  do i = 1, m

    j = z(i)

    do k = 1, s
      if ( x(i,k) == 1 ) then
        a(j,k) = a(j,k) + 1
      end if
    end do

  end do
!
!  Determine the cluster "energies".
!
  e(1:n) = 0

  do j = 1, n

    l = mj(j)

    if ( l < m0 ) then
      iflag = 2
      return
    end if

    e(j) = 0
    do k = 1, s
      e(j) = e(j) + min ( a(j,k), l - a(j,k) )
    end do

  end do

  d = sum ( e(1:n) )

  if ( n <= 1 ) then
    return
  end if

  is = 0
  i = 0

  do

    is = is + 1

    if ( m < is ) then
      exit
    end if

    do

      i = i + 1

      if ( m < i ) then
        it = it + 1
        if ( 15 < it ) then
          return
        end if
        i = 1
      end if

      p = z(i)
      l = mj(p)

      if ( l <= m0 ) then
        exit
      end if

      c(1:n) = 0

      do k = 1, s

        do j = 1, n

          y = a(j,k)
          l = mj(j)

          if ( j /= p ) then
            l = l + 1
            if ( x(i,k) == 1 ) then
              y = y + 1
            end if
          else
            l = l - 1
            if ( x(i,k) == 1 ) then
              y = y - 1
            end if
          end if

          c(j) = c(j) + min ( y, l - y )

        end do

      end do

      eq = - huge ( eq )
      ep = c(p) - e(p)

      do j = 1, n
        if ( j /= p ) then
          ej = e(j) - c(j)
          if ( eq < ej ) then
            eq = ej
            q = j
          end if
        end if
      end do

      if ( eq <= ep ) then
        exit
      end if

      is = 0
      d = d + ep - eq
      e(p) = c(p)
      e(q) = c(q)
      mj(p) = mj(p) - 1
      mj(q) = mj(q) + 1

      do k = 1, s
        if ( x(i,k) == 1 ) then
          a(p,k) = a(p,k) - 1
          a(q,k) = a(q,k) + 1
        end if
      end do

    end do

  end do

  return
end
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
subroutine clrexm ( a, m, s, b, n, z, m0, mj, x, e, d, it, iflag )

!*****************************************************************************80
!
!! CLREXM implements the exchange algorithm for clusterwise linear regression.
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
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, pages 137-139,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(M,S), the independent variables.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input, real ( kind = 8 ) B(M), the dependent variables.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    is assigned.
!
!    Input, integer ( kind = 4 ) M0, the minimum number of items per cluster.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Output, real ( kind = 8 ) X(N,S), the individual regression coefficients.
!
!    Output, real ( kind = 8 ) E(N), the per-cluster value of the L1 criterion.
!
!    Output, real ( kind = 8 ) D, the total value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations that were taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) a(m,s)
  real ( kind = 8 ) ai(s)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) bi
  real ( kind = 8 ) cj
  real ( kind = 8 ) cp
  real ( kind = 8 ) cq
  real ( kind = 8 ) d
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) ej
  real ( kind = 8 ) ep
  real ( kind = 8 ) eq
  real ( kind = 8 ) f(s,n)
  real ( kind = 8 ) fj(s)
  real ( kind = 8 ) fp(s)
  real ( kind = 8 ) fq(s)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnmax
  integer ( kind = 4 ) ns
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ) r(n,(s*(s-1))/2)
  real ( kind = 8 ) rj((s*(s-1))/2)
  real ( kind = 8 ) rp((s*(s-1))/2)
  real ( kind = 8 ) rq((s*(s-1))/2)
  real ( kind = 8 ) t(s,n)
  real ( kind = 8 ) tj(s)
  real ( kind = 8 ) tp(s)
  real ( kind = 8 ) tq(s)
  real ( kind = 8 ) wi
  real ( kind = 8 ) x(n,s)
  integer ( kind = 4 ) z(m)

  iflag = 0
  it = 0
  d = 0.0D+00
  ns = ( s * ( s - 1 ) ) / 2

  if ( m0 < s ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLREXM - Fatal error!'
    write ( *, '(a)' ) '  M0 < S.'
    iflag = 14
    stop
  end if

  do i = 1, m
    j = z(i)
    if ( j < 1 .or. n < j ) then
      iflag = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CLREXM - Fatal error!'
      write ( *, '(a)' ) '  Z(I) out of bounds.'
      stop
    end if
  end do

  e(1:n) = 0.0D+00
  f(1:s,1:n) = 0.0D+00
  t(1:s,1:n) = 0.0D+00
  r(1:n,1:ns) = 0.0D+00
!
!  Determine the cluster populations.
!
  mj(1:n) = 0
  do i = 1, m
    j = z(i)
    mj(j) = mj(j) + 1
  end do

  do i = 1, m

    j = z(i)
    wi = 1.0D+00
    bi = b(i)
    ai(1:s) = a(i,1:s)

    call inexcl ( s, ai, bi, wi, j, n, f, t, r, fj, tj, rj, &
      e(j), nnmax, kmax )

    f(1:kmax,j) = fj(1:kmax)
    t(1:kmax,j) = tj(1:kmax)
    r(j,1:nnmax) = rj(1:nnmax)

  end do
!
!  Make sure every cluster has enough population.
!
  do j = 1, n

    if ( mj(j) < m0 ) then
      iflag = 2
      return
    end if

  end do

  d = sum ( e(1:n) )

  if ( n <= 1 ) then
    go to 21
  end if

  is = 0
  i = 0

11 continue

  is = is + 1

  if ( m < is ) then
    go to 21
  end if

  do

    i = i + 1

    if ( m < i ) then

      it = it + 1

      if ( 15 < it )then
        go to 21
      end if

      i = 1

    end if

    p = z(i)

    if ( mj(p) <= m0 ) then
      go to 11
    end if

    eq = huge ( eq )

    do j = 1, n

      ej = e(j)
      bi = b(i)
      ai(1:s) = a(i,1:s)

      if ( j == p ) then

        ep = ej
        wi = -1.0D+00

        call inexcl ( s, ai, bi, wi, p, n, f, t, r, fp, tp, rp, &
          ep, nnmax, kmax )

        cp = ep

      else

        wi = 1.0D+00

        call inexcl ( s, ai, bi, wi, j, n, f, t, r, fj, tj, rj, &
          ej, nnmax, kmax )

        cj = ej
        ej = cj - e(j)

        if ( ej < eq ) then
          eq = ej
          cq = cj
          q = j
          fq(1:s) = fj(1:s)
          tq(1:s) = tj(1:s)
          rq(1:ns) = rj(1:ns)
        end if

      end if

    end do

    ep = e(p) - cp

    if ( ep <= eq ) then
      exit
    end if

    is = 0
    mj(p) = mj(p) - 1
    mj(q) = mj(q) + 1
    d = d - ep + eq
    e(p) = cp
    e(q) = cq
    f(1:s,p) = fp(1:s)
    f(1:s,q) = fq(1:s)
    t(1:s,p) = tp(1:s)
    t(1:s,q) = tq(1:s)
    r(p,1:ns) = rp(1:ns)
    r(q,1:ns) = rq(1:ns)
    z(i) = q

  end do

  go to 11

21 continue

  do j = 1, n
    do k = s, 1, -1

      tj(k) = t(k,j)

      if ( k /= s ) then
        nn = ( ( k - 1 ) * ( 2 * s - k ) ) / 2 + 1
        do l = k+1, s
          tj(k) = tj(k) - r(j,nn) * tj(l)
          nn = nn + 1
        end do
      end if

      x(j,k) = tj(k)

    end do
  end do

  return
end
subroutine cluster_d_show ( m, n, x, j1, j2, z, rows, columns )

!*****************************************************************************80
!
!! CLUSTER_D_SHOW makes a typewriter plot of points with associated labels.
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
!    Input, real ( kind = 8 ) X(M,N), the data.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns of X to associate with
!    the X and Y directions.
!
!    Input, integer ( kind = 4 ) Z(M), the clusters associated with the data.
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
  character i4_to_a
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
  integer ( kind = 4 ) z(m)

  if ( rows <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLUSTER_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  ROWS <= 0.'
    stop
  end if

  if ( columns <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLUSTER_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  COLUMNS <= 0.'
    stop
  end if

  if ( j1 < 1 .or. n < j1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLUSTER_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of J1.'
    stop
  end if

  if ( j2 < 1 .or. n < j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLUSTER_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of J2.'
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
  x2_min = minval ( x(1:m,j2) )
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
      string(i1,i2) = i4_to_a ( z(i) )
    else
      string(i1,i2) = '*'
    end if

  end do

  do i = rows+1, 0, -1
    write ( *, '(2x,78a1)' ) ( string(i,j), j = 0, columns+1 )
  end do

  return
end
subroutine cluster_d_print ( m, n, x, j1, j2, z )

!*****************************************************************************80
!
!! CLUSTER_D_PRINT prints out the cluster information.
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
!    Input, real ( kind = 8 ) X(M,N), the data.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns of X to associate with
!    the X and Y directions.
!
!    Input, integer ( kind = 4 ) Z(M), the clusters associated with the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) x(m,n)
  integer ( kind = 4 ) z(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Item  X(J1), X(J2), Cluster'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,i4,2g14.6,i4)' ) i, x(i,j1), x(i,j2), z(i)
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
    stop
  end if

  if ( j2 < 1 .or. n < j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_D_SHOW - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of J2.'
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
    write ( *, '(80a1)' ) ( string(i,j), j = 0, columns+1 )
  end do

  return
end
subroutine data_d_write ( file_name, m, n, x )

!*****************************************************************************80
!
!! DATA_D_WRITE writes a real data set into a file.
!
!  Discussion:
!
!    The data set can be thought of as a real M by N array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row (pair of values) at a time.
!
!    Each row begins on a new line.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Modified:
!
!    16 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to write.
!
!    Input, integer ( kind = 4 ) M, the number of data items.
!
!    Input, integer ( kind = 4 ) N, the number of components in a data item.
!
!    Input, real ( kind = 8 ) X(M,N), the data values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  character ( len = * ) file_name
  character ( len = 20 ) form
  integer ( kind = 4 ) i
  integer ( kind = 4 ) output
  real ( kind = 8 ) x(m,n)

  call get_unit ( output )

  open ( unit = output, file = file_name, status = 'replace' )

  write ( output, '(a)' ) '# ' // trim ( file_name )
  write ( output, '(a)' ) '#'
  write ( output, '(a)' ) '#  Created by DATA_D_WRITE.'
  write ( output, '(a,i6)' ) '#  Number of data records is ', m
  write ( output, '(a,i6)' ) '#  Number of data columns is ', n
  write ( output, '(a)' ) '#'

  write ( form, '(a,i6,a)' ) '(', n, 'g14.6)'

  do i = 1, m

    write ( output, form ) x(i,1:n)

  end do

  close ( unit = output )

  return
end
subroutine data_d2_read ( file_name, m, x1, x2 )

!*****************************************************************************80
!
!! DATA_D2_READ reads a data set of pairs of real numbers stored in a file.
!
!  Discussion:
!
!    The data set can be thought of as a real M by 2 array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row (pair of values) at a time.
!
!    Each row (pair of values) begins on a new line.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Modified:
!
!    02 April 2002
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
!    Output, real ( kind = 8 ) X1(M), X2(M), the data values.
!
  implicit none

  integer ( kind = 4 ) m

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) m2
  real ( kind = 8 ) x1(m)
  real ( kind = 8 ) x2(m)

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old' )

  x1(1:m) = huge ( x1(1) )
  x2(1:m) = huge ( x2(1) )

  m2 = 0

  do

    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      m2 = m2 + 1

      last = 0
      call s_to_r8 ( line(last+1:), x1(m2), ierror, length )

      if ( ierror /= 0 ) then
        exit
      end if

      last = last + length

      call s_to_r8 ( line(last+1:), x2(m2), ierror, length )

      if ( ierror /= 0 ) then
        exit
      end if

      if ( m2 == m ) then
        exit
      end if

    end if

  end do

  close ( unit = input )

  return
end
subroutine data_d2_write ( file_name, m, x1, x2 )

!*****************************************************************************80
!
!! DATA_D2_WRITE writes a data set of pairs of real numbers into a file.
!
!  Discussion:
!
!    The data set can be thought of as a real M by 2 array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row (pair of values) at a time.
!
!    Each row (pair of values) begins on a new line.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Modified:
!
!    15 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to write.
!
!    Input, integer ( kind = 4 ) M, the number of data items.
!
!    Input, real ( kind = 8 ) X1(M), X2(M), the data values.
!
  implicit none

  integer ( kind = 4 ) m

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) output
  real ( kind = 8 ) x1(m)
  real ( kind = 8 ) x2(m)

  call get_unit ( output )

  open ( unit = output, file = file_name, status = 'replace' )

  write ( output, '(a)' ) '# ' // trim ( file_name )
  write ( output, '(a)' ) '#'
  write ( output, '(a)' ) '#  Created by DATA_D2_WRITE.'
  write ( output, '(a,i6)' ) '#  Number of data records is ', m
  write ( output, '(a,i6)' ) '#  Number of data columns is ', 2
  write ( output, '(a)' ) '#'

  do i = 1, m

    write ( output, '(2g14.6)' ) x1(i), x2(i)

  end do

  close ( unit = output )

  return
end
subroutine data_i_print ( m, n, a, title )

!*****************************************************************************80
!
!! DATA_I_PRINT prints an integer matrix.
!
!  Modified:
!
!    12 April 2002
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
!    Input, integer ( kind = 4 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 10
    jhi = min ( jlo + 9, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,10i7)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine data_i_read ( file_name, m, n, x )

!*****************************************************************************80
!
!! DATA_I_READ reads an integer data set stored in a file.
!
!  Discussion:
!
!    The data set can be thought of as an integer M by N array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row at a time.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Each row begins on a new line, but may extend over more than
!    one line.
!
!    Individual data values should be separated by spaces or commas.
!
!  Modified:
!
!    12 April 2002
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
!    Output, integer ( kind = 4 ) X(M,N), the data values.
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
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(m,n)

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old' )

  x(1:m,1:n) = huge ( x(1,1) )

  i = 1
  j = 0

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
        call s_to_i4 ( line(last+1:), value, ierror, length )
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

      end do

    end if

  end do

  close ( unit = input )

  return
end
subroutine data_i_show ( m, n, x, j1, j2, rows, columns )

!*****************************************************************************80
!
!! DATA_I_SHOW makes a typewriter plot of an integer data set.
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
  integer ( kind = 4 ) x(m,n)
  integer ( kind = 4 ) x1_max
  integer ( kind = 4 ) x1_min
  integer ( kind = 4 ) x2_max
  integer ( kind = 4 ) x2_min

  if ( rows <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_I_SHOW - Fatal error!'
    write ( *, '(a)' ) '  ROWS <= 0.'
    stop
  end if

  if ( columns <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_I_SHOW - Fatal error!'
    write ( *, '(a)' ) '  COLUMNS <= 0.'
    stop
  end if

  if ( j1 < 1 .or. n < j1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_I_SHOW - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of J1.'
    stop
  end if

  if ( j2 < 1 .or. n < j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_I_SHOW - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of J2.'
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
      real ( columns, kind = 8 ) * real ( x(i,j1) - x1_min, kind = 8 ) &
      / real ( x1_max - x1_min ), kind = 8 )
    i2 = max ( i2, 1 )
    i2 = min ( i2, columns )

    i1 = 1 + nint ( &
      real ( rows, kind = 8    ) * real ( x(i,j2) - x2_min, kind = 8 ) &
      / real ( x2_max - x2_min, kind = 8 ) )
    i1 = max ( i1, 1 )
    i1 = min ( i1, rows )

    if ( string(i1,i2) == ' ' ) then
      string(i1,i2) = '*'
    else if ( string(i1,i2) == '*' ) then
      string(i1,i2) = '@'
    end if

  end do

  do i = rows + 1, 0, -1
    write ( *, '(80a1)' ) ( string(i,j), j = 0, columns+1 )
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
subroutine detexm ( x, m, s, z, m0, mj, xbar, n, detw, it, iflag )

!*****************************************************************************80
!
!! DETEXM implements the exchange algorithm for the determinant criterion.
!
!  Discussion:
!
!    Extraneous arguments were dropped, by the grace of FORTRAN 90.
!    The calling sequence for this routine is now ALMOST the same as for
!    TRWMDM and TRWEXM.
!
!  Modified:
!
!    03 April 2002
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
!    Ellis Horwood, 1985, page 110,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    is assigned.
!
!    Input, integer ( kind = 4 ) M0, the minimum number of data items in each
!    cluster.  This is usually set to 1.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Output, real ( kind = 8 ) XBAR(N,S), the centers of mass, or means, of
!    each of the clusters.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, real ( kind = 8 ) DETW, the value of the determinant.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) alphaj
  real ( kind = 8 ) alphap
  real ( kind = 8 ) alphaq
  real ( kind = 8 ) det
  real ( kind = 8 ) detw
  real ( kind = 8 ) dist
  real ( kind = 8 ) ej
  real ( kind = 8 ) ep
  real ( kind = 8 ) eps1
  real ( kind = 8 ) eps2
  real ( kind = 8 ) eps3
  real ( kind = 8 ) eq
  real ( kind = 8 ) f(s)
  real ( kind = 8 ) fp
  real ( kind = 8 ) fq
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ), parameter :: r = 0.999D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w(s,s)
  real ( kind = 8 ) wbar(s,s)
  real ( kind = 8 ) x(m,s)
  real ( kind = 8 ) xbar(n,s)
  real ( kind = 8 ) xi(s)
  integer ( kind = 4 ) z(m)
!
  iflag = 0
  it = 0
  eps1 = 1.0D-06
  eps2 = 1.0D-06
  eps3 = 1.0D-06
  detw = 0.0D+00

  call means ( x, m, s, z, m0, mj, xbar, n, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DETEXM - Fatal error!'
    write ( *, '(a)' ) '  Error return from MEANS.'
    return
  end if

  do k = 1, s
    w(k,k:s) = 0.0D+00
  end do

  do i = 1, m

    j = z(i)

    f(1:s) = xbar(j,1:s) - x(i,1:s)

    do k = 1, s
      h = f(k)
      do l = k, s
        w(k,l) = w(k,l) + h * f(l)
      end do
    end do

  end do

  call ldlt ( s, w, eps1, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DETEXM - Fatal error!'
    write ( *, '(a)' ) '  Error return from LDLT.'
    return
  end if

  detw = 1.0D+00
  do k = 1, s
    detw = detw * w(k,k)
  end do

  if ( detw <= eps2 ) then
    iflag = 9
    return
  end if

  if ( n <= 1 ) then
    return
  end if

  is = 0
  i = 0

8  continue

  is = is + 1

  if ( m < is ) then
    return
  end if

9 continue

  i = i + 1

  if ( m < i ) then
    it = it + 1
    if ( 15 < it ) then
      return
    end if
    i = 1
  end if

10 continue

  p = z(i)
  l = mj(p)

  if ( mj(p) <= m0 ) then
    go to 8
  end if

  alphap = - real ( mj(p), kind = 8 ) / real ( mj(p) - 1, kind = 8 )

  xi(1:s) = x(i,1:s)
  f(1:s) = xbar(p,1:s) - x(i,1:s)

  do k = 1, s
    wbar(k,1:k) = w(k,1:k)
  end do

  call update ( s, alphap, f, wbar, det, eps3, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DETEXM - Fatal error!'
    write ( *, '(a)' ) '  Error return from UPDATE.'
    return
  end if

  eq = huge ( eq )

  do j = 1, n

    call distw ( xbar, n, xi, wbar, s, j, dist, f )

    if ( j == p ) then

      ep = - alphap * dist

    else

      l = mj(j)
      alphaj = real ( l, kind = 8 ) / real ( l + 1, kind = 8 )
      ej = alphaj * dist

      if ( ej < eq ) then
        eq = ej
        q = j
        alphaq = alphaj
      end if

    end if

  end do

  if ( ep * r <= eq ) then
    go to 8
  end if

  f(1:s) = xbar(q,1:s) - xi(1:s)

  call update ( s, alphaq, f, wbar, det, eps3, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DETEXM - Fatal error!'
    write ( *, '(a)' ) '  Error return from UPDATE.'
    return
  end if

  if ( detw <= det ) then
    go to 8
  end if

  is = 0
  detw = det
  l = mj(p)
  u = l
  k = mj(q)
  v = k
  fp = 1.0D+00 / ( u - 1.0D+00 )
  fq = 1.0D+00 / ( v + 1.0D+00 )
  mj(p) = l - 1
  mj(q) = k + 1

  do k = 1, s
    t = xi(k)
    xbar(p,k) = ( u * xbar(p,k) - t ) * fp
    xbar(q,k) = ( v * xbar(q,k) + t ) * fq
    w(k,1:k) = wbar(k,1:k)
  end do

  z(i) = q
  go to 9

end
subroutine distw ( xbar, n, xi, a, s, j, distj, v )

!*****************************************************************************80
!
!! DISTW computes the squared generalized distance of points from centers.
!
!  Discussion:
!
!    The squared generalized distance of a point Xi from the center
!    XBARj of its cluster is defined as
!
!      DISTJ = transpose ( Xi - XBARj ) * inverse W(C') * ( Xi - XBARj )
!
!    where W(C') is a matrix.
!
!  Modified:
!
!    01 April 2002
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
!    Ellis Horwood, 1985, page 109,
!    QA278 S68213.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) XBAR(N,S), the centers of mass, or means, of
!    each of the clusters.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Input, real ( kind = 8 ) XI(S), the point whose distance from its cluster
!    center is desired.
!
!    Input, real ( kind = 8 ) A(S,S), contains the LDL or Cholesky decomposition
!    of the the matrix W(C').
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input, integer ( kind = 4 ) J, the cluster to which the point XI belongs.
!
!    Output, real ( kind = 8 ) DISTJ, the squared distance of the point XI from
!    the J-th center.
!
!    Output, real ( kind = 8 ) V(S), ?
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) a(s,s)
  real ( kind = 8 ) distj
  real ( kind = 8 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) v(s)
  real ( kind = 8 ) xbar(n,s)
  real ( kind = 8 ) xi(s)

  distj = 0.0D+00

  do k = 1, s

    h = xbar(j,k) - xi(k)
    do l = 1, k - 1
      h = h - a(k,l) * v(l)
    end do

    v(k) = h

    distj = distj + h**2 / a(k,k)

  end do

  return
end
subroutine dwbexm ( x, m, s, z, m0, mj, xbar, n, beta, detwjb, d, it, &
  iflag )

!*****************************************************************************80
!
!! DWBEXM implements the exchange method for the adaptive distance criterion.
!
!  Modified:
!
!    24 June 2002
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
!    Ellis Horwood, 1985, pages 113-114,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    is assigned.
!
!    Input, integer ( kind = 4 ) M0, the minimum number of data items in each
!    cluster.  This is usually set to 1.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Output, real ( kind = 8 ) XBAR(N,S), the centers of mass, or means, of
!    each of the clusters.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Input, real ( kind = 8 ) BETA, the exponent used in the objective
!    function.  BETA must be greater than 0.
!
!    Output, real ( kind = 8 ) DETWJB(N), the value of the per-cluster
!    objective functions.
!
!    Output, real ( kind = 8 ) D, the value of the objective function.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) d
  logical, parameter :: debug = .false.
  real ( kind = 8 ) detbj
  real ( kind = 8 ) detbp
  real ( kind = 8 ) detbq
  real ( kind = 8 ) detj
  real ( kind = 8 ) detwjb(n)
  real ( kind = 8 ) ej
  real ( kind = 8 ) ep
  real ( kind = 8 ) eps1
  real ( kind = 8 ) eps2
  real ( kind = 8 ) eps3
  real ( kind = 8 ) eq
  real ( kind = 8 ) f(s)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ), parameter :: r = 0.999D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w(n,s,s)
  real ( kind = 8 ) wj(s,s)
  real ( kind = 8 ) wp(s,s)
  real ( kind = 8 ) wq(s,s)
  real ( kind = 8 ) x(m,s)
  real ( kind = 8 ) xbar(n,s)
  integer ( kind = 4 ) z(m)

  if ( beta <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DWBEXM - Fatal error!'
    write ( *, '(a)' ) '  Input value of BETA <= 0.'
    stop
  end if

  d = 0.0D+00
  iflag = 0
  it = 0

  if ( m0 <= s + 1 ) then
    m0 = s + 1
  end if

  if ( m < m0 * n ) then
    iflag = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DWBEXM - Fatal error!'
    write ( *, '(a)' ) '  Minimal population times number of clusters'
    write ( *, '(a)' ) '  exceeds total population available.'
    return
  end if

  eps1 = 1.0D-07
  eps2 = 1.0D-07
  eps3 = 1.0D-07

  call means ( x, m, s, z, m0, mj, xbar, n, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DWBEXM - Fatal error!'
    write ( *, '(a)' ) '  Error return from MEANS.'
    return
  end if

  do j = 1, n

    call wjscat ( x, m, s, z, xbar, n, j, wj )

    call ldlt ( s, wj, eps1, iflag )

    if (iflag /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DWBEXM - Fatal error!'
      write ( *, '(a)' ) '  Error return from LDLT.'
      return
    end if

    h = 1.0D+00
    do k = 1, s
      h = h * wj(k,k)
      w(j,k,1:k) = wj(k,1:k)
    end do

    if ( h < eps2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DWBEXM - Fatal error!'
      write ( *, '(a)' ) '  H is too small.'
      iflag = 5
      return
    end if

    h = h**beta
    detwjb(j) = h
    d = d + h

  end do

  if ( n <= 1 ) then
    return
  end if

  i = 0
  is = 1

  if ( m < is ) then
    return
  end if

  do

    i = i + 1

    if ( m < i ) then
      it = it + 1
      if ( 15 < it ) then
        exit
      end if
      i = 1
    end if

    p = z(i)
    l = mj(p)

    if ( l <= m0 ) then

      is = is + 1

      if ( m < is ) then
        exit
      end if

      cycle

    end if

    v = l
    eq = huge ( eq )
    detbp = detwjb(p)

    do j = 1, n

      u = mj(j)

      do k = 1, s
        f(k) = xbar(j,k) - x(i,k)
        wj(k,1:k) = w(j,k,1:k)
      end do

      if ( j == p ) then
        alpha = -real ( mj(j), kind = 8 ) / real ( mj(j) - 1, kind = 8 )
      else
        alpha =  real ( mj(j), kind = 8 ) / real ( mj(j) + 1, kind = 8 )
      end if

      call update ( s, alpha, f, wj, detj, eps3, iflag )

      if ( iflag /= 0 ) then

        if ( debug ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DWBEXM - Fatal error!'
          write ( *, '(a)' ) '  Error return from UPDATE.'
        end if

        return

      end if

      detbj = detj**beta

      if ( j == p ) then

        ep = detbp - detbj
        detbp = detbj
        do k = 1, s
          wp(k,1:k) = wj(k,1:k)
        end do

      else

        ej = detbj - detwjb(j)

        if ( ej < eq ) then

          eq = ej
          detbq = detbj
          do k = 1, s
            wq(k,1:k) = wj(k,1:k)
          end do

          q = j
          h = u

        end if

      end if

    end do

    if ( ep * r <= eq ) then

      is = is + 1

      if ( m < is ) then
        exit
      end if

    else

      is = 0
      detwjb(p) = detbp
      detwjb(q) = detbq
      d = d - ep + eq

      do k = 1, s
        w(p,k,1:k) = wp(k,1:k)
        w(q,k,1:k) = wq(k,1:k)
      end do

      do k = 1, s
        xbar(p,k) = ( v * xbar(p,k) - x(i,k) ) / ( v - 1.0D+00 )
        xbar(q,k) = ( h * xbar(q,k) + x(i,k) ) / ( h + 1.0D+00 )
      end do

      z(i) = q
      mj(p) = mj(p) - 1
      mj(q) = mj(q) + 1

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
function i4_to_a ( i )

!*****************************************************************************80
!
!! I4_TO_A returns the I-th alphabetic character.
!
!  Examples:
!
!    I  I4_TO_A
!
!   -8  ' '
!    0  ' '
!    1  'A'
!    2  'B'
!   ..
!   26  'Z'
!   27  'a'
!   52  'z'
!   53  ' '
!   99  ' '
!
!  Modified:
!
!    23 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the letter to be returned.
!    0 is a space;
!    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
!    27 through 52 requests 'a' through 'z', (ASCII 97:122);
!
!    Output, character I4_TO_A, the requested alphabetic letter.
!
  implicit none

  integer ( kind = 4 ), parameter :: cap_shift = 64
  integer ( kind = 4 ) i
  character i4_to_a
  integer ( kind = 4 ), parameter :: low_shift = 96

  if ( i <= 0 ) then
    i4_to_a = ' '
  else if ( 1 <= i .and. i <= 26 ) then
    i4_to_a = char ( cap_shift + i )
  else if ( 27 <= i .and. i <= 52 ) then
    i4_to_a = char ( low_shift + i - 26 )
  else if ( 53 <= i ) then
    i4_to_a = ' '
  end if

  return
end
subroutine inexcl ( s, ai, bi, wi, j, n, f, t, r, fj, tj, rj, &
  ss, nnmax, kmax )

!*****************************************************************************80
!
!! INEXCL computes auxilliary arrays F, T and R used to control exchanges.
!
!  Modified:
!
!    10 April 2002
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
!    Ellis Horwood, 1985, page 136,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) S, ?
!
!    Input/output, real ( kind = 8 ) AI(S), ?
!
!    Input/output, real ( kind = 8 ) BI, ?
!
!    Input/output, real ( kind = 8 ) WI, ?
!
!    Input, integer ( kind = 4 ) J, ?
!
!    Input, integer ( kind = 4 ) N, ?
!
!    Input, real ( kind = 8 ) F(S,N), ?
!
!    Input, real ( kind = 8 ) T(S,N), ?
!
!    Input, real ( kind = 8 ) R(N,(S*(S-1))/2), ?
!
!    Output, real ( kind = 8 ) FJ(S), ?
!
!    Output, real ( kind = 8 ) TJ(S), ?
!
!    Output, real ( kind = 8 ) RJ((S*(S-1))/2), ?
!
!    Input/output, real ( kind = 8 ) SS, ?
!
!    Output, integer ( kind = 4 ) NNMAX, the number of values of RJ that were set.
!
!    Output, integer ( kind = 4 ) KMAX, the number of times the loop was executed,
!    and the number of values of TJ and FJ that were set.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) ai(s)
  real ( kind = 8 ) ak
  real ( kind = 8 ) al
  real ( kind = 8 ) bi
  real ( kind = 8 ) dp
  real ( kind = 8 ) f(s,n)
  real ( kind = 8 ) fh
  real ( kind = 8 ) fj(s)
  real ( kind = 8 ) fk
  real ( kind = 8 ) hk
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnmax
  real ( kind = 8 ) r(n,(s*(s-1))/2)
  real ( kind = 8 ) rj((s*(s-1))/2)
  real ( kind = 8 ) rl
  real ( kind = 8 ) ss
  real ( kind = 8 ) t(s,n)
  real ( kind = 8 ) tj(s)
  real ( kind = 8 ) tk
  real ( kind = 8 ) wa
  real ( kind = 8 ) wh
  real ( kind = 8 ) wi

  kmax = 0
  nnmax = 0

  do k = 1, s

    if ( wi == 0.0D+00 ) then
      return
    end if

    ak = ai(k)

    if ( ak == 0.0D+00 ) then
      cycle
    end if

    fk = f(k,j)
    wa = wi * ak
    dp = fk + wa * ak
    hk = 1.0D+00 / dp
    fh = fk * hk
    wh = wa * hk
    wi = wi * fh
    fj(k) = dp
    nn = ( ( k - 1 ) * ( s + s - k ) ) / 2 + 1

    do l = k + 1, s
      al = ai(l)
!
!  NN is out of range...
!
      rl = r(j,nn)
      ai(l) = al - ak * rl
      rj(nn) = fh * rl + wh * al
      nnmax = nn
      nn = nn + 1
    end do

    al = bi
    tk = t(k,j)
    bi = al - ak * tk
    tj(k) = fh * tk + wh * al
    kmax = k

  end do

  ss = ss + wi * bi**2

  return
end
subroutine ldlt ( n, a, eps, iflag )

!*****************************************************************************80
!
!! LDLT computes a Cholesky decomposition of the matrix A.
!
!  Discussion:
!
!    The matrix A is assumed to be positive definite symmetric.
!
!    The Cholesky decomposition of A has the form
!
!    A = L * D * L'
!
!    where L is a unit lower triangular matrix and D is a diagonal matrix.
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
!    Ellis Horwood, 1985, page 108,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).  On input, the matrix to be
!    factored.  On output, the diagonal and lower triangle contain information
!    defining the Cholesky factorization of the input matrix.
!
!    Input, real ( kind = 8 ) EPS, a tolerance.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) eps
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) v(n)

  iflag = 0

  do i = 1, n

    do j = i, n

      h = a(i,j)
      do k = 1, i-1
        h = h - v(k) * a(j,k) * a(i,k)
      end do

      if ( i == j ) then

        if ( h < eps )then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LDLT - Fatal error!'
          write ( *, '(a)' ) '  Matrix is not numerically positive definite.'
          iflag = 4
          return
        end if

        v(i) = h
        a(i,i) = h

      else

        a(j,i) = h / v(i)

      end if

    end do

  end do

  return
end
subroutine means ( x, m, s, z, m0, mj, xbar, n, iflag )

!*****************************************************************************80
!
!! MEANS computes the mean vectors for a given partition.
!
!  Discussion:
!
!    MEANS is used whenever the cluster centers are mean vectors.
!
!  Modified:
!
!    30 March 2002
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
!    Ellis Horwood, 1985, page 100,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the number of columns of data.
!
!    Input, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    belongs.
!
!    Input, integer ( kind = 4 ) M0, the minimum number of data items in each
!    cluster.  This is usually set to 1.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Output, real ( kind = 8 ) XBAR(N,S), the centers of mass, or means, of
!    each of the clusters.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) s
  real ( kind = 8 ) x(m,s)
  real ( kind = 8 ) xbar(n,s)
  integer ( kind = 4 ) z(m)

  do i = 1, m

    j = z(i)

    if ( j < 1 .or. n < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MEANS - Fatal error!'
      write ( *, '(a,i6)' ) '  Illegal cluster index for data item I = ', i
      write ( *, '(a,i6)' ) '  Cluster index is ', j
      iflag = 1
      stop
    end if

  end do

  iflag = 0
!
!  Determine the cluster populations.
!
  mj(1:n) = 0
  do i = 1, m
    j = z(i)
    mj(j) = mj(j) + 1
  end do
!
!  Make sure the populations are large enough.
!
  do j = 1, n

    if ( mj(j) < m0 ) then

      if ( debug ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEANS - Fatal error!'
        write ( *, '(a,i6)' ) '  Problem in cluster ', j
        write ( *, '(a,i6)' ) '  Cluster population is ', mj(j)
        write ( *, '(a,i6)' ) '  Minimum acceptable population is ', m0
      end if

      iflag = 2
      return

    end if

  end do
!
!  Sum up the points in each cluster.
!
  xbar(1:n,1:s) = 0.0D+00
  do i = 1, m
    j = z(i)
    xbar(j,1:s) = xbar(j,1:s) + x(i,1:s)
  end do
!
!  Average the points in each cluster.
!
  do j = 1, n
    xbar(j,1:s) = xbar(j,1:s) / real ( mj(j), kind = 8 )
  end do

  return
end
subroutine median ( b, t, value )

!*****************************************************************************80
!
!! MEDIAN computes quantities needed for the OVPEXM objective function.
!
!  Modified:
!
!    10 April 2002
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
!    Ellis Horwood, 1985, page 127,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) T, the number of entries in B.
!
!    Input, integer ( kind = 4 ) B(T), ?
!
!    Output, integer ( kind = 4 ) VALUE, ???
!
  implicit none

  integer ( kind = 4 ) t

  integer ( kind = 4 ) b(t)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
  integer ( kind = 4 ) value
  integer ( kind = 4 ) z

  z = 0

  do i = 2, t
    z = z + b(i) * ( i - 1 )
  end do

  v = b(1)
  u = sum ( b(2:t) )

  do i = 2, t

    if ( u <= v ) then
      exit
    end if

    z = z + v - u
    u = u - b(i)
    v = v + b(i)

  end do

  value = z

  return
end
subroutine ovpexm ( x, m, s, z, mj, n, e, d, it, t, iflag )

!*****************************************************************************80
!
!! OVPEXM implements the exchange algorithm for the L1 criterion.
!
!  Discussion:
!
!    The routine expects data in integer format.
!
!    The routine requires that the input data be ordinal data.  That is,
!    each data item X(*,1:S) comprises S integers between 1 and T.
!
!    The L1 criterion measures the L1 distance of each point from its
!    cluster median, and sums this over all clusters.
!
!  Modified:
!
!    13 April 2002
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
!    Ellis Horwood, 1985, pages 125-127,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    is assigned.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, integer ( kind = 4 ) E(N), the per-cluster value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) D, the total value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations that were taken.
!
!    Input, integer ( kind = 4 ) T, ?
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t

  integer ( kind = 4 ) a(t,s,n)
  integer ( kind = 4 ) b(t)
  integer ( kind = 4 ) bsum
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e(n)
  integer ( kind = 4 ) ej
  integer ( kind = 4 ) ep
  integer ( kind = 4 ) eq
  integer ( kind = 4 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) x(m,s)
  integer ( kind = 4 ) xik
  integer ( kind = 4 ) z(m)

  do i = 1, m

    j = z(i)

    if ( j < 1 .or. n < j ) then
      iflag = 1
      return
    end if

  end do

  it = 0
  d = 0
  iflag = 0

  if ( t < 2 ) then
    iflag = 13
    return
  end if

  m0 = 1
!
!  Determine the cluster populations.
!
  mj(1:n) = 0
  do i = 1, m
    j = z(i)
    mj(j) = mj(j) + 1
  end do
!
!  Make sure the cluster populations are large enough.
!
  do j = 1, n

    if ( mj(j) < m0 ) then
      iflag = 2
      return
    end if

  end do

  a(1:t,1:s,1:n) = 0

  do i = 1, m

    j = z(i)

    do k = 1, s
      xik = x(i,k)
      do l = 1, t
        if ( xik == l ) then
          a(l,k,j) = a(l,k,j) + 1
        end if
      end do
    end do

  end do
!
!  Determine the cluster "energies".
!
  e(1:n) = 0
  do j = 1, n

    f = 0
    do k = 1, s
      do l = 1, t
        b(l) = a(l,k,j)
      end do
      call median ( b, t, bsum )
      f = f + bsum
    end do

    d = d + f
    e(j) = f

  end do

  if ( n <= 1 ) then
    return
  end if

  is = 0
  i = 0

11  continue

  is = is + 1

  if ( m < is ) then
    return
  end if

12  continue

  i = i + 1
  if ( m < i ) then
    it = it + 1
    if ( 15 < it ) then
      return
    end if
    i = 1
  end if

  p = z(i)

  if ( mj(p) <= m0 ) then
    go to 11
  end if

  do j = 1, n
    c(j) = 0
    do k = 1, s
      xik = x(i,k)
      do l = 1, t
        b(l) = a(l,k,j)
        if ( xik == l ) then
          if ( j == p ) then
            b(l) = b(l) - 1
          else
            b(l) = b(l) + 1
          end if
        end if
      end do
      call median ( b, t, bsum )
      c(j) = c(j) + bsum
    end do
  end do

  eq = - huge ( eq )
  ep = c(p) - e(p)

  do j = 1, n
    if ( j /= p ) then
      ej = e(j) - c(j)
      if ( eq < ej ) then
        eq = ej
        q = j
      end if
    end if
  end do

  if ( ep < eq ) then

    is = 0
    d = d + ep - eq
    e(p) = c(p)
    e(q) = c(q)
    mj(p) = mj(p) - 1
    mj(q) = mj(q) + 1
    do k = 1, s
      xik = x(i,k)
      do l = 1, t
        if ( xik == l ) then
          a(l,k,p) = a(l,k,p) - 1
          a(l,k,q) = a(l,k,q) + 1
        end if
      end do
    end do

    z(i) = q

    go to 12

  end if

  go to 11

end
subroutine ovrexm ( x, m, s, z, mj, n, e, d, it, iflag )
!
!*****************************************************************************80
!
!! OVREXM implements the exchange algorithm for the L1 criterion.
!
!  Discussion:
!
!    The routine expects data in integer format.
!
!    OVSEXM and OVREXM are similar routines.  The main difference is
!    that OVSEXM maintains a sorted copy of the data matrix, and hence
!    requires more space, but runs faster.
!
!    The L1 criterion measures the L1 distance of each point from its
!    cluster median, and sums this over all clusters.
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
!    Ellis Horwood, 1985, pages 121-124,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    is assigned.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, integer ( kind = 4 ) E(N), the per-cluster value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) D, the total value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations that were taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  integer ( kind = 4 ) d
  integer ( kind = 4 ) e(n)
  integer ( kind = 4 ) ej
  integer ( kind = 4 ) ep
  integer ( kind = 4 ) eq
  logical even(n)
  logical evenp
  integer ( kind = 4 ) f
  logical final
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) r1
  integer ( kind = 4 ) u(n,s)
  integer ( kind = 4 ) uj(s)
  integer ( kind = 4 ) ujk
  integer ( kind = 4 ) up(s)
  integer ( kind = 4 ) upk
  integer ( kind = 4 ) uq(s)
  integer ( kind = 4 ) v(n,s)
  integer ( kind = 4 ) vj(s)
  integer ( kind = 4 ) vjk
  integer ( kind = 4 ) vp(s)
  integer ( kind = 4 ) vpk
  integer ( kind = 4 ) vq(s)
  integer ( kind = 4 ) x(m,s)
  integer ( kind = 4 ) xik
  integer ( kind = 4 ) xu
  integer ( kind = 4 ) xv
  integer ( kind = 4 ) y(m)
  integer ( kind = 4 ) yr
  integer ( kind = 4 ) yr1
  integer ( kind = 4 ) z(m)

  do i = 1, m

    j = z(i)

    if ( j < 1 .or. n < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OVREXM - Fatal error!'
      iflag = 1
      return
    end if

  end do

  iflag = 0
  it = 0
  d = 0
  m0 = 1
!
!  Determine the cluster populations.
!
  mj(1:n) = 0
  do i = 1, m
    j = z(i)
    mj(j) = mj(j) + 1
  end do
!
!  Make sure the cluster populations are large enough.
!
  do j = 1, n

    if ( mj(j) < m0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OVREXM - Fatal error!'
      iflag = 2
      return
    end if

  end do

  e(1:n) = 0
  do j = 1, n

    l = mj(j)
    l1 = l / 2
    l2 = ( l + 1 ) / 2
    even(j) = ( 2 * l1 ) == l

    do k = 1, s

      r = 0
      do i = 1, m
        if ( z(i) == j ) then
          r = r + 1
          y(r) = x(i,k)
        end if
      end do

      do p = 2, l

        final = .true.
        q = l - p + 1

        do r = 1, q

          r1 = r + 1
          yr = y(r)
          yr1 = y(r+1)

          if ( yr1 < yr ) then
            y(r+1) = yr
            y(r) = yr1
            final = .false.
          end if

        end do

        if ( final ) then
          exit
        end if

      end do

      h = y(l2)
      u(j,k) = h
      if ( even(j) ) then
        h = y(l2+1)
      end if

      v(j,k) = h
      f = 0
      do r = 1, l
        f = f + abs ( y(r) - h )
      end do

      e(j) = e(j) + f
      d = d + f

    end do

  end do

  if ( n <= 1 ) then
    return
  end if

  is = 0
  i = 0

11 continue

  is = is + 1
  if ( m < is ) then
    return
  end if

12 continue

  i = i + 1

  if ( m < i ) then
    it = it + 1
    if ( 15 < it ) then
      return
    end if
    i = 1
  end if

  p = z(i)
  r = mj(p)

  if ( r <= m0 ) then
    go to 11
  end if

  r1 = r / 2
  ep = 0
  evenp = even(p)

  do k = 1, s

    upk = u(p,k)
    vpk = v(p,k)
    xik = x(i,k)

    if ( .not. evenp ) then

      lu = 0
      lv = 0
      xu = - huge ( xu )
      xv = huge ( xv )

      do l = 1, m

        if ( z(l) == p .and. l /= i ) then

          if ( x(l,k) < upk ) then
            xu = max ( xu, x(l,k) )
            lu = lu + 1
          else if ( upk < x(l,k) ) then
            xv = min ( xv, x(l,k) )
            lv = lv + 1
          end if

        end if

      end do

      if ( lu < r1 ) then
        xu = upk
      end if

      if ( lv < r1 ) then
        xv = upk
      end if

      if ( xik == upk ) then
        up(k) = xu
        vp(k) = xv
      else if ( xik < upk ) then
        up(k) = upk
        vp(k) = xv
        ep = ep + abs ( xik - upk )
      else
        up(k) = xu
        vp(k) = upk
        ep = ep + abs ( xik - upk )
      end if

    else

      if ( xik <= upk ) then
        up(k) = vpk
        vp(k) = vpk
        ep = ep + abs ( xik - vpk )
      else
        up(k) = upk
        vp(k) = upk
        ep = ep + abs ( xik - upk )
      end if

    end if

  end do

20 continue

  eq = huge ( eq )

  do j = 1, n

    if ( j == p ) then
      cycle
    end if

    r = mj(j)
    r1 = r / 2
    ej = 0

    do k = 1, s

      ujk = u(j,k)
      vjk = v(j,k)
      xik = x(i,k)

      if ( .not. even(j) ) then

        lr = 0

        if ( xik <= ujk ) then

          uj(k) = xik

          if ( xik /= ujk .and. r /= 1 ) then

            xu = - huge ( xu )
            do l = 1, m
              if ( z(l) == j ) then
                if ( x(l,k) < ujk ) then
                  lr = lr + 1
                  xu = max ( xu, x(l,k) )
                end if
              end if
            end do

            if ( lr < r1 ) then
              xu = ujk
            end if

            if ( xik < xu ) then
              uj(k) = xu
            end if

          end if

          vj(k) = vjk
          ej = ej + abs ( xik - ujk )

        else

          vj(k) = xik

          if ( r /= 1 ) then

            xv = huge ( xv )
            do l = 1, m
              if ( z(l) == j ) then
                if ( vjk < x(l,k) ) then
                  lr = lr + 1
                  xv = max ( xv, x(l,k) )
                end if
              end if
            end do

            if ( lr < r1 ) then
              xv = vjk
            end if

            if ( xv < xik ) then
              vj(k) = xv
            end if

          end if

          uj(k) = ujk
          ej = ej + abs ( xik - ujk )

        end if

      else

        if ( vjk < xik ) then

          uj(k) = vjk
          vj(k) = vjk
          ej = ej + abs ( xik - vjk )

        else if ( xik <= ujk ) then

          uj(k) = ujk
          vj(k) = ujk
          ej = ej + abs ( xik - ujk )

        else

          uj(k) = xik
          vj(k) = xik

        end if

      end if

    end do

    if ( ej < eq ) then
      eq = ej
      q = j
      do k = 1, s
        uq(k) = uj(k)
        vq(k) = vj(k)
      end do
    end if

  end do

  if ( eq < ep ) then

    e(p) = e(p) - ep
    e(q) = e(q) + eq
    d = d - ep + eq
    mj(p) = mj(p) - 1
    mj(q) = mj(q) + 1
    even(p) = .not. even(p)
    even(q) = .not. even(q)

    do k = 1, s
      u(p,k) = up(k)
      v(p,k) = vp(k)
      u(q,k) = uq(k)
      v(q,k) = vq(k)
    end do

    z(i) = q
    is = 0
    go to 12

  end if

  go to 11

! return
end
subroutine ovsexm ( x, m, s, z, mj, n, e, d, it, iflag )

!*****************************************************************************80
!
!! OVSEXM implements the exchange algorithm for the L1 criterion.
!
!  Problems:
!
!    This routine returns NEGATIVE values for E and D in many cases,
!    which should not be possible.
!
!  Discussion:
!
!    The routine expects data in integer format.
!
!    OVSEXM and OVREXM are similar routines.  The main difference is
!    that OVSEXM maintains a sorted copy of the data matrix, and hence
!    requires more space, but runs faster.
!
!    The L1 criterion measures the L1 distance of each point from its
!    cluster median, and sums this over all clusters.
!
!  Modified:
!
!    12 April 2002
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
!    Ellis Horwood, 1985, page 117-120,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    is assigned.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, integer ( kind = 4 ) E(N), the per-cluster value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) D, the total value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations that were taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  integer ( kind = 4 ) d
  integer ( kind = 4 ) e(n)
  integer ( kind = 4 ) ej
  integer ( kind = 4 ) ep
  integer ( kind = 4 ) eq
  logical even(n)
  logical evenp
  integer ( kind = 4 ) f
  logical final
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) mk(n+1)
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) mq
  integer ( kind = 4 ) mr
  integer ( kind = 4 ) ms
  integer ( kind = 4 ) mt
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  logical revers
  integer ( kind = 4 ) sj(s)
  integer ( kind = 4 ) sp(s)
  integer ( kind = 4 ) sq(s)
  integer ( kind = 4 ) u(n,s)
  integer ( kind = 4 ) uj(s)
  integer ( kind = 4 ) ujk
  integer ( kind = 4 ) up(s)
  integer ( kind = 4 ) upk
  integer ( kind = 4 ) uq(s)
  integer ( kind = 4 ) v(n,s)
  integer ( kind = 4 ) vj(s)
  integer ( kind = 4 ) vjk
  integer ( kind = 4 ) vp(s)
  integer ( kind = 4 ) vpk
  integer ( kind = 4 ) vq(s)
  integer ( kind = 4 ) x(m,s)
  integer ( kind = 4 ) xik
  integer ( kind = 4 ) xu
  integer ( kind = 4 ) xv
  integer ( kind = 4 ) xx(m,s)
  integer ( kind = 4 ) y(m)
  integer ( kind = 4 ) yr
  integer ( kind = 4 ) yr1
  integer ( kind = 4 ) z(m)

  iflag = 0
  it = 0
  d = 0
  m0 = 1

  mj(1:n) = 0
  e(1:n) = 0

  do i = 1, m

    j = z(i)

    if ( j < 1 .or. n < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OVSEXM - Fatal error!'
      write ( *, '(a)' ) '  J is out of bounds.'
      iflag = 1
      return
    end if

    mj(j) = mj(j) + 1

  end do

  mk(1) = 1

  do j = 1, n

    l = mj(j)

    if ( l < m0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OVSEXM - Fatal error!'
      write ( *, '(a)' ) '  Cluster size too small.'
      iflag = 2
      return
    end if

    l1 = l / 2
    l2 = ( l + 1 ) / 2
    even(j) = ( 2 * l1 ) == l
    mk(j+1) = mk(j) + l

    do k = 1, s

      r = 0
      do i = 1, m
        if ( z(i) == j ) then
          r = r + 1
          y(r) = x(i,k)
        end if
      end do

      do q = l-1, 1, -1

        final = .true.

        do r = 1, q

          yr = y(r)
          yr1 = y(r+1)

          if ( yr1 < yr ) then
            y(r+1) = yr
            y(r) = yr1
            final = .false.
          end if

        end do

        if ( final ) then
          exit
        end if

      end do

      u(j,k) = y(l2)

      if ( even(j) ) then
        h = y(l2+1)
      else
        h = y(l2)
      end if

      v(j,k) = h
      f = 0
      do r = 1, l
        yr = y(r)
        mr = mk(j) + r - 1
        xx(mr,k) = yr
        f = f + abs ( yr - h )
      end do

      e(j) = e(j) + f
      d = d + f

    end do

  end do

  if ( n <= 1 ) then
    return
  end if

  is = 0
  i = 0

11 continue

  is = is + 1

  if ( m < is ) then
    return
  end if

12 continue

  i = i + 1

  if ( m < i ) then
    it = it + 1
    if ( 15 < it ) then
      return
    end if
    i = 1
  end if

  p = z(i)
  r = mj(p)

  if ( r <= m0 ) then
    go to 11
  end if

  mp = mk(p)
  g = ( r + 1 ) / 2 + mk(p) - 1
  ep = 0
  evenp = even(p)

  do k = 1, s

    upk = u(p,k)
    vpk = v(p,k)
    xik = x(i,k)

    do ms = mk(p) + mj(p) - 1, mk(p), -1

      if ( xik == xx(ms,k) ) then
        sp(k) = ms
      end if

    end do

    if ( .not. evenp ) then

      xu = xx(g-1,k)
      xv = xx(g+1,k)

      if ( xik == upk ) then
        up(k) = xu
        vp(k) = xv
      else if ( xik < upk ) then
        up(k) = upk
        vp(k) = xv
        ep = ep + abs ( xik - upk )
      else
        up(k) = xu
        vp(k) = upk
        ep = ep + abs ( xik - upk )
      end if

    else

      if ( xik <= upk ) then
        up(k) = vpk
        vp(k) = vpk
        ep = ep + abs ( xik - vpk )
      else
        up(k) = upk
        vp(k) = upk
        ep = ep + abs ( xik - upk )
      end if

    end if

  end do

20 continue

  eq = huge ( eq )

  do j = 1, n

    if ( j == p ) then
      cycle
    end if

    r = mj(j)
    ej = 0
    mr = mk(j)
    final = ( 1 < r )
    g = ( r + 1 ) / 2 + mr - 1
    ms = mk(j+1)

    do k = 1, s
      ujk = u(j,k)
      vjk = v(j,k)
      xik = x(i,k)
      mt = ms
      do l = 1, r
        mp = mr + l - 1
        if ( xik <= xx(mp,k) ) then
          mt = mp
          exit
        end if
      end do

      sj(k) = mt

      if ( .not. even(j) ) then

        if ( xik <= ujk ) then

          uj(k) = xik
          xu = xik

          if ( final ) then
            xu = xx(g-1,k)
          end if

          if ( xik < xu ) then
            uj(k) = xu
          end if

          vj(k) = vjk
          ej = ej + abs ( xik - ujk )

        else

          vj(k) = xik
          uj(k) = ujk
          xv = xik

          if ( final ) then
            xv = xx(g+1,k)
          end if

          if ( xv < xik ) then
            vj(k) = xv
          end if

          ej = ej + abs ( xik - ujk )

        end if

      else

        if ( vjk <= xik ) then

          uj(k) = vjk
          vj(k) = vjk
          ej = ej + abs ( xik - vjk )

        else if ( xik <= ujk ) then

          uj(k) = ujk
          vj(k) = ujk
          ej = ej + abs ( xik - ujk )

        else

          uj(k) = xik
          vj(k) = xik

        end if

      end if

    end do

    if ( ej < eq ) then
      eq = ej
      q = j
      do k = 1, s
        uq(k) = uj(k)
        vq(k) = vj(k)
        sq(k) = sj(k)
      end do
    end if

  end do

  if ( ep <= eq ) then
    go to 11
  end if

  e(p) = e(p) - ep
  e(q) = e(q) + eq
  d = d - ep + eq

  mj(p) = mj(p) - 1
  mj(q) = mj(q) + 1
  even(p) = .not. even(p)
  even(q) = .not. even(q)
  revers = ( q < p )

  do k = 1, s

    u(p,k) = up(k)
    v(p,k) = vp(k)
    u(q,k) = uq(k)
    v(q,k) = vq(k)
    xik = x(i,k)
    mp = sp(k)
    mq = sq(k)

    if ( .not. revers ) then
      if ( mp == mq - 1 ) then
        do l = mp, mq - 2
          xx(l,k) = xx(l+1,k)
        end do
        xx(mq-1,k) = xik
      end if
    else
      h = xx(mq,k)
      xx(mq,k) = xik
      do l = mq + 1, mp
        f = xx(l,k)
        xx(l,k) = h
        h = f
      end do
    end if

  end do
!
!  Update the MK array to account for change in MJ.
!
  if ( .not. revers ) then
    do l = p + 1, q
      mk(l) = mk(l) - 1
    end do
  else
    do l = q + 1, p
      mk(l) = mk(l) + 1
    end do
  end if

  z(i) = q
  is = 0
  go to 12

end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Modified:
!
!    12 September 2004
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
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! RMAT_PRINT_SOME prints some of an R8MAT.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  logical d_is_int
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)') j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine randp ( m, m0, n, z, seed )

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
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, integer ( kind = 4 ) Z(M), the cluster to which each item is assigned.
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
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) urand
  integer ( kind = 4 ) z(m)
!
!  Use the first N * M0 data items to guarantee that each cluster has
!  at least M0 elements.
!
  if ( m < n * m0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDP - Fatal error!'
    write ( *, '(a)' ) '  Not enough data to satisfy occupancy requirements.'
    write ( *, '(a)' ) '  Require N * M0 <= M.'
    write ( *, '(a,i6)' ) '  but N = ', n
    write ( *, '(a,i6)' ) '  M0 =    ', m0
    write ( *, '(a,i6)' ) '  M =     ', m
    stop
  end if

  k = 0
  do i = 1, n
    do j = 1, m0
      k = k + 1
      z(k) = i
    end do
  end do
!
!  Now take care of the remaining data items.
!
  do i = k+1, m

    j = int ( real ( n, kind = 8 ) * urand ( seed ) ) + 1
    j = min ( j, n )
    j = max ( j, 1 )

    z(i) = j

  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads a real number from a string.
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
!  Examples:
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
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
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
subroutine tihexm ( m, t, n, mj, method, z, e, d, it, iflag )

!*****************************************************************************80
!
!! TIHEXM implements the exchange algorithm for center-free critera.
!
!  Discussion:
!
!    This routine attempts to minimize an objective function based on
!    the clustering of data points.  The only means of affecting the
!    objective function is to "exchange", that is, to move one point
!    at a time from one cluster to another.
!
!    The objective function is the sum of functions associated with
!    each cluster.  Each cluster objective function is based on the
!    distances between all pairs of points in the cluster.
!
!    There are three variations of the cluster objective function,
!    based on how the sum of distances is to be weighted.
!
!    In the text, the distance matrix is stored as a compressed vector.
!    For convenience and clarity, the full matrix is used here.
!
!    In the original version of the code, the calculation of EP
!    would blow up if a cluster contained only 1, or only 2 points.
!    This problem has been corrected by setting EP to 0 in such cases.
!
!  Modified:
!
!    31 October 2004
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
!    Ellis Horwood, 1985, pages 133-134,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, real ( kind = 8 ) T(M,M), the distance matrix.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Input, integer ( kind = 4 ) METHOD, specifies the cluster objective
!    function used.
!    1, the sum of distances between all pairs of points in the cluster;
!    2, the sum of distances between all pairs of points in the cluster,
!       divided by the number of points in the cluster;
!    3, the sum of distances between all pairs of points in the cluster,
!       divided by (twice) the number of pairs.
!
!    Output, integer ( kind = 4 ) Z(M), the cluster to which each item
!    is assigned.
!
!    Output, real ( kind = 8 ) E(N), the per-cluster value of the L1 criterion.
!
!    Output, real ( kind = 8 ) D, the total value of the L1 criterion.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations that were taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) bj
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) cj
  real ( kind = 8 ) d
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) ej
  real ( kind = 8 ) ep
  real ( kind = 8 ) eq
  real ( kind = 8 ) f
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) method
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ), parameter :: r = 0.999D+00
  real ( kind = 8 ) t(m,m)
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  integer ( kind = 4 ) z(m)

  iflag = 0
  it = 0
  d = 0.0D+00

  if ( method < 1 .or. 3 < method ) then
    iflag = 10
    return
  end if

  if ( method == 1 ) then
    m0 = 1
  else if ( method == 2 ) then
    m0 = 1
  else
    m0 = 2
  end if
!
!  Determine the cluster populations.
!
  mj(1:n) = 0

  do i = 1, m
    j = z(i)
    if ( j < 1 .or. n < j ) then
      iflag = 1
      return
    end if
    mj(j) = mj(j) + 1
  end do

  do j = 1, n
    if ( mj(j) < m0 ) then
      iflag = 2
      return
    end if
  end do
!
!  Determine the "energies".
!
  e(1:n) = 0.0D+00

  do i = 1, m
    do j = 1, m
      if ( i /= j ) then
        if ( z(i) == z(j) ) then
          e(z(i)) = e(z(i)) + t(i,j)
        end if
      end if
    end do
  end do
!
!  Revise the energy for the given method.
!
  do j = 1, n

    if ( method == 1 ) then
      f = e(j)
    else if ( method == 2 ) then
      f = e(j) / real ( mj(j), kind = 8 )
    else if ( method == 3 ) then
      f = e(j) / real ( mj(j) * ( mj(j) - 1 ), kind = 8 )
    end if

    e(j) = f

  end do

  d = sum ( e(1:n) )

  if ( n <= 1 ) then
    return
  end if

  is = 1
  i = 0

  do

    i = i + 1

    if ( m < i ) then
      it = it + 1
      if ( 15 < it ) then
        exit
      end if
      i = 1
    end if

    p = z(i)
    l = mj(p)

    if ( mj(p) <= m0 ) then

      is = is + 1

      if ( m < is ) then
        exit
      end if

    end if

    v = real ( mj(p), kind = 8 )
    c(1:n) = 0.0D+00

    do h = 1, m

      if ( h /= i ) then
        j = z(h)
        c(j) = c(j) + t(h,i)
      end if

    end do

    eq = huge ( eq )

    do j = 1, n

      bj = e(j)
      cj = c(j)
      u = mj(j)

      if ( j == p ) then

        if ( method == 1 ) then
          ep = cj
        else if ( method == 2 ) then
          if ( v == 1.0D+00 ) then
            ep = 0.0D+00
          else
            ep = ( cj - bj ) / ( v - 1.0D+00 )
          end if
        else if ( method == 3 ) then
          if ( v == 1.0D+00 .or. v == 2.0D+00 ) then
            ep = 0.0D+00
          else
            ep = ( cj - 2.0D+00 * ( v - 1.0D+00 ) * bj ) &
              / ( ( v - 1.0D+00 ) * ( v - 2.0D+00 ) )
          end if
        end if

      else

        if ( method == 1 ) then

          ej = cj

        else if ( method == 2 ) then

          ej = ( cj - bj ) / ( u + 1.0D+00 )

        else if ( method == 3 ) then

          if ( u == 0.0D+00 ) then
            ej = 0.0D+00
          else
            ej = ( cj - 2.0D+00 * u * bj ) / ( u * ( u + 1.0D+00 ) )
          end if

        end if

        if ( ej < eq ) then
          eq = ej
          q = j
          w = u
        end if

      end if

    end do

    if ( eq < ep * r ) then

      is = 0
      d = d - ep + eq

      if ( method == 1 ) then

        e(p) = e(p) - c(p)
        e(q) = e(q) + c(q)

      else if ( method == 2 ) then

        if ( v == 1.0D+00 ) then
          e(p) = 0.0D+00
        else
          e(p) = ( v * e(p) - c(p) ) / ( v - 1.0D+00 )
        end if

        e(q) = ( w * e(q) + c(q) ) / ( w + 1.0D+00 )

      else if ( method == 3 ) then

        if ( v == 1.0D+00 .or. v == 2.0D+00 ) then
          e(p) = 0.0D+00
        else
          e(p) = ( v * ( v - 1.0D+00 ) * e(p) - c(p) ) &
            / ( ( v - 1.0D+00 ) * ( v - 2.0D+00 ) )
        end if

        if ( w == 0.0D+00 ) then
          e(q) = 0.0D+00
        else
          e(q) = ( w * ( w - 1.0D+00 ) * e(q) + c(q) ) / ( w * ( w + 1.0D+00 ) )
        end if

      end if

      mj(p) = mj(p) - 1
      mj(q) = mj(q) + 1
      z(i) = q

    else

      is = is + 1

      if ( m < is ) then
        return
      end if

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
subroutine traces ( x, m, s, z, xbar, n, e, d )

!*****************************************************************************80
!
!! TRACES computes the per-cluster and total variances.
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
!    Ellis Horwood, 1985, page 100,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    belongs.
!
!    Input, real ( kind = 8 ) XBAR(N,S), the centers of mass, or means, of the
!    clusters.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, real ( kind = 8 ) E(N), the per-cluster variance or energy.
!
!    Output, real ( kind = 8 ) D, the total variance or energy.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) d
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,s)
  real ( kind = 8 ) xbar(n,s)
  integer ( kind = 4 ) z(m)

  e(1:n) = 0.0D+00

  do i = 1, m
    j = z(i)
    e(j) = e(j) + sum ( ( xbar(j,1:s) - x(i,1:s) )**2 )
  end do

  d = sum ( e(1:n) )

  return
end
subroutine trafor ( x, m, s )

!*****************************************************************************80
!
!! TRAFOR standardizes a data matrix.
!
!  Discussion:
!
!    The transformed matrix has the property that each column has
!    mean 0 and variance 1.
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
!    Ellis Horwood, 1985, page 99,
!    QA278 S68213.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(M,S), the M by S data matrix.  On output,
!    the columns of the matrix have been shifted and scaled to have
!    mean 0 and variance 1.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the number of columns of data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) s

  real ( kind = 8 ) h
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(m,s)
  real ( kind = 8 ) xbar

  do k = 1, s

    xbar = sum ( x(1:m,k) ) / real ( m, kind = 8 )

    h = sum ( ( x(1:m,k) - xbar )**2 )

    if ( 0.0D+00 < h ) then

      x(1:m,k) = ( x(1:m,k) - xbar ) / sqrt ( h )

    else

      x(1:m,k) = 0.0D+00

    end if

  end do

  return
end
subroutine trwexm ( x, m, s, z, m0, mj, xbar, n, e, d, it, iflag )

!*****************************************************************************80
!
!! TRWEXM implements the exchange method for the variance criterion.
!
!  Modified:
!
!    31 March 2002
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
!    Ellis Horwood, 1985, page 105,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of items of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    belongs.
!
!    Input, integer ( kind = 4 ) M0, the minimum number of data items in each
!    cluster.  This is usually set to 1.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Output, real ( kind = 8 ) XBAR(N,S), the centers of mass, or means, of the
!    clusters.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, real ( kind = 8 ) E(N), the per-cluster variance or energy.
!
!    Output, real ( kind = 8 ) D, the total variance or energy.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) d
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) ej
  real ( kind = 8 ) ep
  real ( kind = 8 ) eq
  real ( kind = 8 ) f
  real ( kind = 8 ) fp
  real ( kind = 8 ) fq
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 15
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ), parameter :: r = 0.999D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  integer ( kind = 4 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x(m,s)
  real ( kind = 8 ) xbar(n,s)
  integer ( kind = 4 ) z(m)
!
  iflag = 0
!
!  Find the means of the clusters.
!
  call means ( x, m, s, z, m0, mj, xbar, n, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRWEXM - Fatal error!'
    write ( *, '(a)' ) '  Error return from MEANS.'
    return
  end if
!
!  Find the per-cluster and total variances.
!
  call traces ( x, m, s, z, xbar, n, e, d )

  if ( n <= 1 ) then
    return
  end if

  it = 0
  is = 0
  i = 0

1 continue

  is = is + 1

  if ( m < is ) then
    return
  end if

2 continue

  i = i + 1

  if ( m < i ) then

    it = it + 1

    if ( it_max < it ) then
      return
    end if

    i = 1

  end if

3 continue

  p = z(i)
  l = mj(p)

  if ( l <= m0 ) then
    go to 1
  end if

  v = l
  eq = huge ( eq )

  do j = 1, n

    h = 0.0D+00
    do k = 1, s
      t = xbar(j,k) - x(i,k)
      h = h + t * t
    end do

    if ( j == p ) then
      f = v / ( v - 1.0D+00 )
      ep = h * f
      cycle
    end if

    u = mj(j)
    f = u / ( u + 1.0D+00 )
    ej = h * f

    if ( ej < eq ) then
      eq = ej
      q = j
      w = u
    end if

  end do

  if ( ep * r <= eq ) then
    go to 1
  end if

  is = 0
  e(p) = e(p) - ep
  e(q) = e(q) + eq
  d = d - ep + eq
  fp = 1.0D+00 / ( v - 1.0D+00 )
  fq = 1.0D+00 / ( w + 1.0D+00 )

  do k = 1, s
    h = x(i,k)
    xbar(p,k) = ( v * xbar(p,k) - h ) * fp
    xbar(q,k) = ( w * xbar(q,k) + h ) * fq
  end do

  z(i) = q
  mj(p) = l - 1
  mj(q) = mj(q) + 1
  go to 2

end
subroutine trwmdm ( x, m, s, z, m0, mj, xbar, n, e, d, it, iflag )

!*****************************************************************************80
!
!! TRWMDM implements the minimal distance method for the variance criterion.
!
!  Modified:
!
!    31 March 2002
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
!    Ellis Horwood, 1985, page 102,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input/output, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    belongs.
!
!    Input, integer ( kind = 4 ) M0, the minimum number of data items in each
!    cluster.  This is usually set to 1.
!
!    Output, integer ( kind = 4 ) MJ(N), the number of data items in each cluster.
!
!    Output, real ( kind = 8 ) XBAR(N,S), the centers of mass, or means, of the
!    clusters.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Output, real ( kind = 8 ) E(N), the per-cluster variance or energy.
!
!    Output, real ( kind = 8 ) D, the total variance or energy.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) d
  logical, parameter :: debug = .false.
  real ( kind = 8 ) dmax
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 15
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mj(n)
  real ( kind = 8 ), parameter :: r = 0.999D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) x(m,s)
  real ( kind = 8 ) xbar(n,s)
  integer ( kind = 4 ) z(m)

  iflag = 0
  dmax = huge ( dmax )
  it = 0

  do

    if ( it_max <= it ) then
      exit
    end if

    it = it + 1
!
!  Find the means of the clusters.
!
    call means ( x, m, s, z, m0, mj, xbar, n, iflag )
!
!  MEANS returns an error flag if, in particular, any
!  cluster has become empty.
!
    if ( iflag /= 0 ) then

      if ( debug ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRWMDM - Fatal error!'
        write ( *, '(a)' ) '  Error return from MEANS.'
      end if

      exit

    end if
!
!  Find the per-cluster and total variances.
!
    call traces ( x, m, s, z, xbar, n, e, d )

    if ( n <= 1 ) then
      exit
    end if

    if ( dmax <= d ) then
      exit
    end if

    dmax = d * r
!
!  Assign each data item to the cluster whose center is nearest.
!
    do i = 1, m

      f = huge ( f )

      do j = 1, n

        h = 0.0D+00
        do k = 1, s
          t = xbar(j,k) - x(i,k)
          h = h + t * t
        end do

        if ( h < f ) then
          f = h
          l = j
        end if

      end do

      z(i) = l

    end do

  end do

  return
end
subroutine update ( s, alpha, v, a, det, eps, iflag )

!*****************************************************************************80
!
!! UPDATE updates the Cholesky decomposition of a matrix.
!
!  Discussion:
!
!    It is assumed that the relationship between the new matrix B
!    and the old matrix A is:
!
!      B = A + alpha * V * V'
!
!    where ALPHA is a scalar, and V is a vector.
!
!    Both A and B are assumed to be positive definite symmetric
!    matrices.
!
!    The Cholesky factorization of A is already known, having been
!    computed by routine LDLT.  It is desired to update the Cholesky
!    factorization of A to obtain that of B.
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
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 108,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) S, the order of the matrices.
!
!    Input/output, real ( kind = 8 ) ALPHA, the modification factor.
!
!    Input/output, real ( kind = 8 ) V(S), the modification vector.
!
!    Input/output, real ( kind = 8 ) A(S,S), the Cholesky factorization, which
!    has been updated on output.
!
!    Output, real ( kind = 8 ) DET, the determinant of the modified matrix.
!    The determinant is required to be greater than the user specified
!    tolerance EPS.
!
!    Input, real ( kind = 8 ) EPS, a tolerance for the size of the determinant.
!    For a reasonable test, EPS should be a small positive value.
!
!    Output, integer ( kind = 4 ) IFLAG, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) s

  real ( kind = 8 ) a(s,s)
  real ( kind = 8 ) alj
  real ( kind = 8 ) alpha
  real ( kind = 8 ) betaj
  logical, parameter :: debug = .false.
  real ( kind = 8 ) dj
  real ( kind = 8 ) djbar
  real ( kind = 8 ) det
  real ( kind = 8 ) eps
  real ( kind = 8 ) h
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real ( kind = 8 ) pj
  real ( kind = 8 ) v(s)

  iflag = 0
  det = 1.0D+00

  do j = 1, s

    pj = v(j)
    h = alpha * pj
    dj = a(j,j)
    djbar = dj + h * pj

    if ( djbar < eps ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UPDATE - Fatal error!'
      write ( *, '(a)' ) '  DJBAR < EPS.'
      write ( *, '(a,g14.6)' ) '  DJBAR = ', djbar
      write ( *, '(a,g14.6)' ) '  EPS   = ', eps
      iflag = 6
      return
    end if

    a(j,j) = djbar
    det = det * djbar

    if ( j /= s ) then

      betaj = h / djbar
      alpha = dj * alpha / djbar

      do l = j+1, s
        alj = a(l,j)
        v(l) = v(l) - pj * alj
        a(l,j) = alj + betaj * v(l)
      end do

    end if

  end do

  if ( det < eps ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UPDATE - Fatal error!'
      write ( *, '(a)' ) '  DET < EPS.'
      write ( *, '(a,g14.6)' ) '  DET = ', det
      write ( *, '(a,g14.6)' ) '  EPS = ', eps
    end if

    iflag = 7

  end if

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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
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

  if ( mic < seed) then
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
subroutine wjscat ( x, m, s, z, xbar, n, j, wj )

!*****************************************************************************80
!
!! WJSCAT calculates the scatter matrix for a given cluster.
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
!    Ellis Horwood, 1985, page 112,
!    QA278 S68213.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(M,S), the M by S data matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows of data.
!
!    Input, integer ( kind = 4 ) S, the spatial dimension of the data.
!
!    Input, integer ( kind = 4 ) Z(M), the cluster to which each data item
!    belongs.
!
!    Output, real ( kind = 8 ) XBAR(N,S), the centers of mass, or means, of the
!    clusters.
!
!    Input, integer ( kind = 4 ) N, the number of clusters.
!
!    Input, integer ( kind = 4 ) J, the cluster whose scatter matrix is desired.
!
!    Output, real ( kind = 8 ) WJ(S,S), the scatter matrix for cluster J.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 8 ) f(s)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) wj(s,s)
  real ( kind = 8 ) x(m,s)
  real ( kind = 8 ) xbar(n,s)
  integer ( kind = 4 ) z(m)

  do k = 1, s
    wj(k,k:s) = 0.0D+00
  end do

  do i = 1, m

    l = z(i)

    if ( l == j ) then

      do k = 1, s
        f(k) = xbar(j,k) - x(i,k)
      end do

      do k = 1, s
        do l = k, s
          wj(k,l) = wj(k,l) + f(k) * f(l)
        end do
      end do

    end if

  end do

  return
end
