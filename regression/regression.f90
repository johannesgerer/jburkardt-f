subroutine a328_li ( a, m, n, b, x, res_norm, iflag )

!*****************************************************************************80
!
!! A328_LI minimizes the L-infinity norm of A*x-b.
!
!  Modified:
!
!    12 February 2003
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 88-92,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Output, real ( kind = 4 ) X(N), the computed solution.
!
!    Output, real ( kind = 4 ) RES_NORM, the L-infinity norm of the 
!    residual vector.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    1, the solution may be inaccurate due to rounding errors;
!    2, condition 1 is not fulfilled.
!    3, condition 2 is not fulfilled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) a1
  real ( kind = 4 ) b(m)
  integer ( kind = 4 ) bb
  real ( kind = 4 ) cc
  real ( kind = 4 ) cnorm
  real ( kind = 4 ) dd
  real ( kind = 4 ) eps
  real ( kind = 4 ) eps2
  real ( kind = 4 ) eps3
  logical finish
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: itmax = 10
  integer ( kind = 4 ) ix(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) l
  real ( kind = 4 ) lam(n+1)
  real ( kind = 4 ) lastep
  integer ( kind = 4 ) lst
  real ( kind = 4 ) p(n+1,n+1)
  integer ( kind = 4 ) r(n+1)
  integer ( kind = 4 ) r_temp
  real ( kind = 4 ) ref
  integer ( kind = 4 ) refset(n+1)
  real ( kind = 4 ) res_norm
  integer ( kind = 4 ) ri
  integer ( kind = 4 ) rj
  real ( kind = 4 ) rv(n+1)
  real ( kind = 4 ) s
  real ( kind = 4 ) snorm
  real ( kind = 4 ) ssum
  double precision sum2
  real ( kind = 4 ) sv(n+1)
  real ( kind = 4 ) t
  real ( kind = 4 ) tab(n+1,m)
  real ( kind = 4 ) w(n+1)
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) xr(n+1)
  real ( kind = 4 ) xx(n+1)

  iflag = 0
  lastep = 0.0E+00
  eps = epsilon ( eps )
  eps3 = - 1.0E+00

  call i4vec_indicator ( n+1, r )

  call i4vec_indicator ( m, ix )

  tab(1:n,1:m) = transpose ( a(1:m,1:n) )
  tab(n+1,1:m) = b(1:m)

  do i = 1, n + 1

    t = 0.0E+00

    do j = i, n + 1

      k = r(j)

      do l = i, m

        ref = tab(k,ix(l))

        if ( abs ( ref ) <= t ) then
          exit
        end if

        s = ref
        t = abs ( ref )
        a1 = j
        bb = l

      end do

    end do

    j = 1

    if ( t == 0.0E+00 ) then
      go to 20
    end if

    lst = r(a1)

    r_temp = r(a1)
    r(a1)  = r(i)
    r(i)   = r_temp

    a1 = ix(bb)

    r_temp = ix(bb)
    ix(bb) = ix(i)
    ix(i)  = r_temp

    do j = i + 1, m
      l = ix(j)
      ref = tab(lst,l) / s
      do k = i+1, n+1
        a1 = r(k)
        tab(a1,l) = tab(a1,l) - tab(a1,a1) * ref
      end do
    end do

  end do

  bb = 1
  a1 = 1

10 continue

  do i = bb, n + 1

    l = ix(i)

    if ( i /= bb ) then
      ii = bb
    else
      ii = 1
    end if

    do j = ii, n+1

      if ( i <= j ) then
        kmax = i - 1
      else
        kmax = j - 1
      end if

      rj = r(j)

      if ( rj == n+1 ) then
        dd = b(l)
      else
        dd = a(l,rj)
      end if

      ssum = - dd
      do k = 1, kmax
        ssum = ssum + p(i,r(k)) * p(k,rj)
      end do
      p(i,rj) = - ssum

    end do

    ref = 0.0E+00

    do j = i, n+1
      t = p(i,r(j))
      if ( ref < abs ( t ) ) then
        ref = abs ( t )
        s = t
        k = j
      end if
    end do

    if ( ref == 0.0E+00 ) then
      j = 1
      exit
    end if

    if ( i == n+1 ) then
      if ( a1 == 1 ) then
        go to 30
      else if ( a1 == 2 ) then
        go to 60
      end if
    end if

    r_temp = r(k)
    r(k)   = r(i)
    r(i)   = r_temp

    do j = i + 1, n + 1
      p(i,r(j)) = p(i,r(j)) / s
    end do

  end do

20 continue

  refset(1:n+1) = ix(1:n+1)

  if ( j == 1 ) then
    iflag = 2
    return
  end if

  if ( j == 2 ) then
    iflag = 3
    return
  end if

30 continue

  do j = bb, n + 1

    rj = r(j)

    if ( rj == n+1 ) then
      dd = - 1.0E+00
    else
      dd = 0.0E+00
    end if

    sv(j) = dd - dot_product ( sv(1:j-1), p(1:j-1,rj) )

  end do

  do j = n + 1, 1, -1
    ssum = - sv(j)
    do k = j+1, n+1
      ssum = ssum + lam(k) * p(k,r(j))
    end do
    lam(j) = - ssum / p(j,r(j))
  end do

  t = sum ( abs ( lam(1:n+1) ) )

  eps2 = 1.0E+00 / t

  if ( eps2 <= lastep ) then
    go to 50
  end if

  lastep = eps2

  do i = 1, n + 1
    xr(i) = sign ( 1.0E+00, lam(i) ) * eps2
  end do

  do i = 1, n + 1
    ssum = - xr(i)
    do j = 1, i-1
      ssum = ssum + w(j) * p(i,r(j))
    end do
    w(i) = - ssum / p(i,r(i))
  end do

  do i = n + 1, 1, -1
    ssum = - w(i)
    do j = i+1, n+1
      ssum = ssum + p(i,r(j)) * xx(r(j))
    end do
    xx(r(i)) = - ssum
  end do

  ref = - xx(n+1)
  xx(1:n) = xx(1:n) / ref

  eps2 = eps2 / ref
  ref = - 1.0E+00

  do j = n + 2, m

    i = ix(j)
    ssum = - b(i)
    do k = 1, n
      ssum = ssum + a(i,k) * xx(k)
    end do

    t = ssum

    if ( ref < abs ( t ) ) then
      ref = abs ( t )
      a1 = j
      s = sign ( 1.0E+00, t )
    end if

  end do

  if ( ref <= eps2 ) then
    go to 60
  end if

40 continue

  k = ix(a1)

  do i = 1, n+1

    if ( r(i) == n+1 ) then
      dd = b(k)
    else
      dd = a(k,r(i))
    end if

    w(i) = dd - dot_product ( w(1:i-1), p(1:i-1,r(i)) )

  end do

  do i = n+1, 1, -1

    ssum = - w(i)

    do j = i+1, n+1
      ssum  = ssum + w(j) * p(j,r(i))
    end do

    w(i) = - ssum / p(i,r(i))

  end do

  ref = lam(n+1)
  bb = n+1

  if ( ref == 0.0E+00 ) then
    j = 2
    go to 20
  end if

  ref = ( w(n+1) / ref ) * s

  do jj = 1, n

    t = lam(jj)

    if ( t == 0.0E+00 ) then
      j = 2
      go to 20
    end if

    t = ( w(jj) / t ) * s

    if ( ref < t ) then
      bb = jj
      ref = t
    end if

  end do

  ix(a1) = ix(bb)
  ix(bb) = k
  a1 = 1
  go to 10

50 continue

  eps2 = lastep

  r_temp = ix(a1)
  ix(a1) = ix(bb)
  ix(bb) = r_temp

  ref = - 1.0E+00

  do j = n+2, m
    i = ix(j)
    sum2 = - dble ( b(i) )
    do k = 1, n
      sum2 = sum2 + dble ( xx(k) ) * dble ( a(i,k) )
    end do
    t = sum2
    if ( ref < abs ( t ) ) then
      ref = abs ( t )
    end if
  end do

60 continue

  lastep = 0.0E+00
  it = 0

  do

    if ( itmax <= it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'A328_LI - Fatal error!'
      write ( *, '(a)' ) '  Number of iterations exceeded.'
      iflag = 2
      return
    end if

    it = it + 1

    cnorm = 0.0E+00
    snorm = 0.0E+00

    do i = 1, n+1

      k = ix(i)
      t = abs ( xx(i) )
      snorm = max ( snorm, t )
      sum2 = - dble ( xr(i) )

      do j = 1, n+1
        if ( j == n+1 ) then
          cc = b(k)
        else
          cc = a(k,j)
        end if
        sum2 = sum2 + dble ( xx(j) ) * dble ( cc )
      end do

      rv(i) = - sum2

    end do

    do i = 1, n+1
      ssum = - rv(i)
      do j = 1, i-1
        ssum = ssum + rv(j) * p(i,r(j))
      end do
      rv(i) = - ssum / p(i,r(i))
    end do

    do i = n+1, 1, -1
      ssum = - rv(i)
      do j = i+1, n+1
        ssum = ssum + p(i,r(j)) * w(r(j))
      end do
      w(r(i)) = - ssum
    end do

    xx(1:n+1) = xx(1:n+1) + w(1:n+1)

    do i = 1, n+1
      cnorm = max ( cnorm, abs ( w(i) ) )
    end do

    if ( cnorm / snorm <= eps ) then
      exit
    end if

  end do

  ref = - xx(n+1)
  xx(1:n) = xx(1:n) / ref

  eps2 = eps2 / ref
  ref = - 1.0E+00

  do j = n+2, m

    i = ix(j)

    sum2 = - dble ( b(i) )
    do k = 1, n
      sum2 = sum2 + dble ( xx(k) ) * dble ( a(i,k) )
    end do

    t = sum2

    if ( ref < abs ( t ) ) then
      ref = abs ( t )
      a1 = j
      s = sign ( 1.0E+00, t )
    end if

  end do

  if ( ref <= eps2 ) then
    finish = .true.
  else if ( eps3 < eps2 ) then
    finish = .false.
  else if ( eps2 <= eps3 ) then
    iflag = 1
    res_norm = eps3
    return
  end if

  eps3 = eps2
  refset(n+1) = ix(n+1)
  refset(1:n) = ix(1:n)

  x(1:n) = xx(1:n)

  if ( .not. finish ) then
    go to 40
  end if

  res_norm = eps3

  return
end
subroutine a478_l1 ( a, m, n, b, eps, rank, x, r, iflag )

!*****************************************************************************80
!
!! A478_L1 minimizes the L1 norm of A*x-b using the modified simplex method.
!
!  Modified:
!
!    15 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 62-64,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) EPS, a tolerance used in the rank determination.
!
!    Output, integer ( kind = 4 ) RANK, a numerical determination of the
!    rank of A.
!
!    Output, real ( kind = 4 ) X(N), the solution vector.
!
!    Output, real ( kind = 4 ) R(M), contains the residuals.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no errors occurred.
!    2, the calculations were prematurely stopped because of rounding
!       errors.
!    3, a fatal error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) a2(m+2,n+2)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) b2(m)
  real ( kind = 4 ) d
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) in
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kl
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) kr
  integer ( kind = 4 ) out
  real ( kind = 4 ) pivot
  real ( kind = 4 ) pmax
  real ( kind = 4 ) pmin
  real ( kind = 4 ) r(m)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) s(m)
  logical stage
  double precision sum2
  logical test
  real ( kind = 4 ) x(n)

  if ( m < n ) then
    iflag = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'A478_L1 - Fatal error!'
    write ( *, '(a)' ) '  This routine requires N <= M.'
    return
  end if
!
!  Initialization.
!
  a2(1:m,1:n) = a(1:m,1:n)

  a2(m+1,1:n+2) = 0.0E+00

  do j = 1, n
    a2(m+2,j) = real ( j, kind = 4 )
  end do

  a2(m+2,n+1:n+2) = 0.0E+00

  do i = 1, m
    a2(i,n+2) = real ( n + i, kind = 4 )
  end do

  a2(1:m,n+1) = b(1:m)

  do i = 1, m
    if ( b(i) < 0.0E+00 ) then
      a2(i,1:n+2) = - a2(i,1:n+2)
    end if
  end do

  b2(1:m) = b(1:m)
  x(1:n) = 0.0E+00

  r(1:m) = 0.0E+00
!
!  Compute the marginal costs.
!
  do j = 1, n + 1
    a2(m+1,j) = sum ( a2(1:m,j) )
  end do
!
!  Stage 1.
!  Determine the vector to enter the basis.
!
  stage = .true.
  kount = 0
  kr = 1
  kl = 1

10 continue

  pmax = - 1.0E+00
  do j = kr, n
    if ( abs ( a2(m+2,j) ) <= n ) then
      d = abs ( a2(m+1,j) )
      if ( pmax < d ) then
        pmax = d
        in = j
      end if
    end if
  end do

  if ( a2(m+1,in) < 0.0E+00 ) then
    a2(m+1:m+2,in) = - a2(m+1:m+2,in)
  end if
!
!  Determine the vector to leave the basis.
!
20 continue

  k = 0
  do i = kl, m
    d = a2(i,in)
    if ( eps < d ) then
      k = k + 1
      b2(k) = a2(i,n+1) / d
      s(k) = i
      test = .true.
    end if
  end do

30 continue

  if ( k <= 0 ) then

    test = .false.

  else

    pmin = huge ( pmin )

    do i = 1, k
      if ( b2(i) < pmin ) then
        j = i
        pmin = b2(i)
        out = s(i)
      end if
    end do

    b2(j) = b2(k)
    s(j) = s(k)
    k = k - 1

  end if
!
!  Check for linear dependence in stage 1.
!
  if ( ( .not. test ) .and. stage ) then

    do i = 1, m + 2
      call r4_swap ( a2(i,kr), a2(i,in) )
    end do

    kr = kr + 1
    go to 40

  end if

  if ( .not. test ) then
    a2(m+2,n+1) = 2.0E+00
    go to 70
  end if

  pivot = a2(out,in)

  if ( eps < a2(m+1,in) - pivot - pivot ) then
    do j = kr, n+1
      d = a2(out,j)
      a2(m+1,j) = a2(m+1,j) - d - d
      a2(out,j) = - d
    end do
    a2(out,n+2) = - a2(out,n+2)
    go to 30
  end if
!
!  Pivot on A2(OUT,IN).
!
  do j = kr, n+1
    if ( j /= in ) then
      a2(out,j) = a2(out,j) / pivot
    end if
  end do

  do j = kr, n+1
    if ( j /= in ) then
      call a478_l1_col ( a2(1,j), a2(1,in), a2(out,j), m+1, out )
    end if
  end do

  do i = 1, m+1
    if ( i /= out ) then
      a2(i,in) = - a2(i,in) / pivot
    end if
  end do

  a2(out,in) = 1.0E+00 / pivot

  call r4_swap ( a2(out,n+2), a2(m+2,in) )

  kount = kount + 1

  if ( .not. stage ) then
    go to 50
  end if
!
!  Interchange rows in stage 1.
!
  kl = kl + 1
  do j = kr, n+2
    call r4_swap ( a2(out,j), a2(kount,j) )
  end do

40 continue

  if ( kount + kr /= n+1 ) then
    go to 10
  end if
!
!  Stage 2.
!
  stage = .false.
!
!  Determine the vector to enter the basis.
!
50 continue

  pmax = - huge ( pmax )

  do j = kr, n

    d = a2(m+1,j)

    if ( d < 0.0E+00 ) then
      if ( -2.0E+00 < d ) then
        go to 60
      end if
      d = - d - 2.0E+00
    end if

    if ( pmax < d ) then
      pmax = d
      in = j
    end if

60 continue

  end do

  if ( eps < pmax ) then

    if ( a2(m+1,in) <= 0.0E+00 ) then

      a2(1:m+2,in) = -a2(1:m+2,in)
      a2(m+1,in) = a2(m+1,in) - 2.0E+00

    end if

    go to 20

  end if
!
!  Prepare output.
!
  do i = 1, kl-1
    if ( a2(i,n+1) < 0.0E+00 ) then
      do j = kr, n+2
        a2(i,j) = - a2(i,j)
      end do
    end if
  end do

  a2(m+2,n+1) = 0.0E+00

  if ( kr == 1 ) then

    do j = 1, n
      d = abs ( a2(m+1,j) )
      if ( d <= eps .or. 2.0E+00 - d <= eps ) then
        go to 70
      end if

    end do

    a2(m+2,n+1) = 1.0E+00

  end if

70 continue

  do i = 1, m

    k = a2(i,n+2)
    d = a2(i,n+1)

    if ( k <= 0 ) then
      k = - k
      d = - d
    end if

    if ( i < kl ) then
      x(k) = d
    else
      k = k - n
      r(k) = d
    end if

  end do

  it = kount
  rank = n + 1 - kr

  iflag = int ( a2(m+2,n+1) + eps )
!
!  As originally written, IFLAG = 0 signals the solution was computed
!  and is not unique, 1 that the solutiion is unique.
!
!  This conflicts with the general role of IFLAG: 0 means no error,
!  nonzero means error.  So I've lumped the two cases together here,
!  and if you're interested, you can unlump them.
!
  if ( iflag == 1 ) then
    iflag = 0
  end if

  a2(m+1,n+1) = sum ( a2(kl:m,n+1) )

  return
end
subroutine a478_l1_col ( v1, v2, s, n, iout )

!*****************************************************************************80
!
!! A478_L1_COL is used by A478_L1.
!
!  Modified:
!
!    03 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) V1(N), the vector to be modified.
!
!    Input, real ( kind = 4 ) V2(N), the vector to be added to V1.
!
!    Input, real ( kind = 4 ) S, the multiplier.
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, integer ( kind = 4 ) IOUT, the index of V1 that is not to be changed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iout
  real ( kind = 4 ) s
  real ( kind = 4 ) v1(n)
  real ( kind = 4 ) v2(n)

  do i = 1, n
    if ( i /= iout ) then
      v1(i) = v1(i) - s * v2(i)
    end if
  end do

  return
end
subroutine a495_li ( a, m, n, b, eps1, eps2, x, rank, res_norm, iflag )

!*****************************************************************************80
!
!! A495_LI minimizes the L-Infinity norm of A*x-b using a simplex method.
!
!  Discussion:
!
!    A modified simplex method is applied to a dual linear program.
!
!    It is permissible for the rank of A to be less than N.
!
!  Modified:
!
!    12 February 2003
!
!  Reference:
!
!    I Barrodale, C Phillips,
!    Algorithm 495: Solution of an Overdetermined System of Linear
!      Equations in the Chebyshev Norm,
!    ACM Transactions on Mathematical Software,
!    Volume 1, pages 264-270, 1975.
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 94-98,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) EPS1, a tolerance for an accuracy test.
!
!    Input, real ( kind = 4 ) EPS2, a tolerance for the minimization condition.
!    Setting EPS2 to 0.0 means that a minimum is sought.  Setting
!    0 < EPS2 allows for a relative error between the norm of the residual
!    that is found and the norm of the residual for the minimal
!    solution.
!
!    Output, real ( kind = 4 ) X(N), the computed solution.
!
!    Output, integer ( kind = 4 ) RANK, an estimate of the rank of the matrix.
!
!    Output, real ( kind = 4 ) RES_NORM, the L-infinity norm of the residual.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    2, the calculations were stopped because of rounding errors.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) at(n+3,m+1)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) bt(m+1)
  real ( kind = 4 ) d
  real ( kind = 4 ) dd
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) level
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) np1mr
  integer ( kind = 4 ) pcol
  real ( kind = 4 ) pivot
  integer ( kind = 4 ) prow
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) val
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) xt(n+3)
!
!  Initialization.
!
  rank = n
  np1mr = n + 1 - rank

  at(1:n,1:m) = transpose ( a(1:m,1:n) )

  call r4vec_indicator ( n, at(1:n,m+1) )

  at(n+1,1:m) = 1.0E+00
  at(n+1,m+1) = 0.0E+00

  at(n+2,1:m) = - b(1:m)
  at(n+2,m+1) = 0.0E+00

  do j = 1, m
    at(n+3,j) = real ( n + j, kind = 4 )
  end do
  at(n+3,m+1) = 0.0E+00

  bt(1:m) = b(1:m)
  bt(m+1) = 0.0E+00

  xt(1:n+3) = 0.0E+00

  it = 0
!
!  The return value IFLAG = 1 is being suppressed in favor of IFLAG = 0.
!
! iflag = 1
  iflag = 0
!
!  Level 1.
!
  level = 1
  k = 0

10 continue

  k = k + 1
  mode = 0
  bt(k:m) = 1.0E+00
!
!  Determine the vector to enter the basis.
!
20 continue

  d = - huge ( d )

  do j = k, m

    if ( bt(j) /= 0.0E+00 ) then
      dd = abs ( at(n+2,j) )
      if ( d < dd ) then
        pcol = j
        d = dd
      end if
    end if

  end do
!
!  Test for zero right-hand side.
!
  if ( k <= 1 ) then
    if ( d <= eps1 ) then
      res_norm = 0.0E+00
      mode = 2
      go to 90
    end if
  end if
!
!  Determine the vector to leave the basis.
!
  d = eps1

  do i = 1, n+1-k
    dd = abs ( at(i,pcol) )
    if ( d < dd ) then
      prow = i
      d = dd
    end if
  end do

  if ( eps1 < d ) then
    go to 80
  end if
!
!  Check for linear dependence in level 1.
!
  bt(pcol) = 0.0E+00

  if ( mode == 1 ) then
    go to 20
  end if

  do j = k, m
    if ( bt(j) /= 0.0E+00 ) then
      do i = 1, n+1-k
        if ( eps1 < abs ( at(i,j) ) ) then
          mode = 1
          go to 20
        end if
      end do
    end if
  end do

  rank = k - 1
  np1mr = n + 1 - rank
  iflag = 0

  go to 40
!
!  Interchange columns PCOL and K in level 1.
!
30 continue

  if ( pcol /= k ) then

    do i = 1, n + 3
      call r4_swap ( at(i,pcol), at(i,k) )
    end do

  end if
!
!  Interchange rows PROW and N+1-K in level 1.
!
  if ( prow /= n+1-k ) then

    do j = 1, m + 1
      call r4_swap ( at(prow,j), at(n+1-k,j) )
    end do

  end if

  if ( k < n ) then
    go to 10
  end if

40 continue

  if ( rank == m ) then
    go to 90
  end if
!
!  Level 2.
!
  level = 2
!
!  Determine the vector to enter the basis.
!
  d = eps1
  do j = rank+1, m
    dd = abs ( at(n+2,j) )
    if ( d < dd ) then
      pcol = j
      d = dd
    end if
  end do
!
!  Compare the Chebyshev error with EPS1.
!
  if ( d <= eps1 ) then
    res_norm = 0.0E+00
    mode = 3
    go to 90
  end if

  if ( -eps1 <= at(n+2,pcol) ) then

    at(n+1,pcol) = 2.0E+00 - at(n+1,pcol)

    do i = np1mr, n+3
      if ( i /= n+1 ) then
        at(i,pcol) = - at(i,pcol)
      end if
    end do

  end if
!
!  Arrange for all entries in the pivot column
!  (except the pivot) to be negative.
!
  do i = np1mr, n

    if ( eps1 <= at(i,pcol) ) then

      do j = 1, m
        at(n+1,j) = at(n+1,j) + 2.0E+00 * at(i,j)
        at(i,j) = - at(i,j)
      end do

      at(i,m+1) = - at(i,m+1)

    end if

  end do

  prow = n + 1
  go to 80

60 continue

  if ( rank + 1 == m ) then
    go to 90
  end if
!
!  Interchange columns PCOL and M in level 2.
!
  if ( pcol /= m ) then

    do i = np1mr, n + 3
      call r4_swap ( at(i,pcol), at(i,m) )
    end do

  end if
!
!  Level 3.
!
  level = 3
!
!  Determine the vector to enter the basis.
!
70 continue

  d = - eps1
  val = 2.0E+00 * at(n+2,m)

  do j = rank + 1, m - 1

    if ( at(n+2,j) < d ) then

      pcol = j
      d = at(n+2,j)
      mode = 0

    else

      dd = val - at(n+2,j)

      if ( dd < d ) then
        mode = 1
        pcol = j
        d = dd
      end if

    end if

  end do

  if ( -eps1 <= d ) then
    go to 90
  end if

  dd = - d / at(n+2,m)

  if ( dd < eps2 ) then
    mode = 4
    go to 90
  end if

  if ( mode /= 0 ) then

    at(np1mr:n+1,pcol) = 2.0E+00 * at(np1mr:n+1,m) - at(np1mr:n+1,pcol)
    at(n+2,pcol) = d
    at(n+3,pcol) = - at(n+3,pcol)

  end if
!
!  Determine the vector to leave the basis.
!
  d = huge ( d )

  do i = np1mr, n + 1
    if ( eps1 < at(i,pcol) ) then
      dd = at(i,m) / at(i,pcol)
      if ( dd < d ) then
        prow = i
        d = dd
      end if
    end if
  end do

  if ( d == huge ( d ) ) then
    iflag = 2
    go to 90
  end if
!
!  Pivot on AT(PROW,PCOL).
!
80 continue

  pivot = at(prow,pcol)

  at(prow,1:m) = at(prow,1:m) / pivot

  do j = 1, m
    if ( j /= pcol ) then
      do i = n+1-rank, n+2
        if ( i /= prow ) then
          at(i,j) = at(i,j) - at(prow,j) * at(i,pcol)
        end if
      end do
    end if
  end do

  at(n+1-rank:n+2,pcol) = - at(n+1-rank:n+2,pcol) / pivot

  at(prow,pcol) = 1.0E+00 / pivot

  call r4_swap ( at(prow,m+1), at(n+3,pcol) )

  it = it + 1

  if ( level == 1 ) then
    go to 30
  else if ( level == 2 ) then
    go to 60
  else
    go to 70
  end if
!
!  Prepare output.
!
90 continue

  bt(1:m) = 0.0E+00

  if ( mode == 2 ) then
    go to 100
  end if

  do j = 1, rank
    k = at(n+3,j)
    xt(k) = at(n+2,j)
  end do

  if ( mode == 3 .or. rank == m ) then
    go to 100
  end if

  do i = np1mr, n+1
    k = abs ( at(i,m+1) ) - real ( n, kind = 4 )
    bt(k) = at(n+2,m) * sign ( 1.0E+00, at(i,m+1) )
  end do

  do j = rank+1, m - 1
    k = abs ( at(n+3,j) ) - real ( n, kind = 4 )
    bt(k) = ( at(n+2,m) - at(n+2,j) ) * sign ( 1.0E+00, at(n+3,j) )
  end do
!
!  Test for a non-unique solution.
!
  do i = np1mr, n+1
    if ( abs ( at(i,m) ) <= eps1 ) then
      iflag = 0
      exit
    end if
  end do

100 continue

  if ( mode /= 2 .and. mode /= 3 ) then
    res_norm = at(n+2,m)
  end if

  if ( rank == m ) then
    res_norm = 0.0E+00
  end if

  if ( mode == 4 ) then
    res_norm = res_norm - d
  end if

  x(1:n) = xt(1:n)

  return
end
subroutine abd_li ( a, m, n, b, eps, rank, x, r, res_norm, iflag )

!*****************************************************************************80
!
!! ABD_LI minimizes the L-infinity norm of A*x-b using a simplex method.
!
!  Discussion:
!
!    Among other peculiarities, the original source code reversed
!    the meaning of M and N.
!
!  Modified:
!
!    15 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 102-104,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!    Note that this routine can handle the case where rank ( A ) < N.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix A.
!
!    Input, B(M), the right hand side.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for an accuracy test.
!
!    Output, integer ( kind = 4 ) RANK, an estimate of the rank of the matrix.
!
!    Output, real ( kind = 4 ) X(N), the solution.
!
!    Output, real ( kind = 4 ) R(M), the residuals.
!
!    Output, real ( kind = 4 ) RES_NORM, the value of the objective function.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    2, a feasible solution could not be found.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) at(n+1,m)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) binv(n+1,n+1)
  real ( kind = 4 ) bv(n+1)
  real ( kind = 4 ) d
  real ( kind = 4 ) e
  real ( kind = 4 ) eps
  real ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibound(m)
  integer ( kind = 4 ) ibv
  integer ( kind = 4 ) icbas(n+1)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) ipart
  integer ( kind = 4 ) irbas(n+1)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) ivo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jin
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kl
  integer ( kind = 4 ) l
  real ( kind = 4 ) piv
  real ( kind = 4 ) pivot
  real ( kind = 4 ) r(m)
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res_norm
  double precision s
  real ( kind = 4 ) thmax
  real ( kind = 4 ) tz
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) zc(m)

  at(1:n,1:m) = transpose ( a(1:m,1:n) )
  at(n+1,1:m) = 1.0E+00

  it = 0
  ipart = 1
  kl = 1
  rank = n
  ibv = 0
  iflag = 0
  itest = 1

  call i4vec_indicator ( n+1, irbas )

  do j = 1, n + 1
    do i = 1, n + 1
      binv(i,j) = 0.0E+00
    end do
    binv(j,j) = 1.0E+00
  end do

  ibound(1:m) = 1

  iout = 0

10 continue

  do

    iout = iout + 1

    if ( n+1 < iout ) then
      go to 40
    end if

    piv = 0.0E+00

    do j = 1, m
      d = abs ( at(iout,j) )
      if ( piv < d ) then
        jin = j
        piv = d
      end if
    end do
!
!  Detection of rank deficiency.
!
    if ( eps < piv ) then
      go to 130
    end if

    if ( iout == n + 1 ) then
      exit
    end if
!
!  Interchange rows IOUT and KL.
!
    do j = 1, m
      call r4_swap ( at(iout,j), at(kl,j) )
    end do

    k = irbas(iout)
    irbas(iout) = irbas(kl)
    irbas(kl) = 0
    x(k) = 0.0E+00

    do j = kl, n
      binv(iout,j) = binv(kl,j)
      binv(kl,j) = 0.0E+00
    end do

    icbas(iout) = icbas(kl)
    icbas(kl) = 0

    do i = kl, n + 1
      binv(i,iout) = binv(i,kl)
      binv(i,kl) = 0.0E+00
    end do

    rank = rank - 1
    kl = kl + 1

  end do

  do j = 1, m

    if ( any ( icbas(kl:n) == j ) ) then

    else
      jin = j
      exit
    end if

  end do

  b(jin) = -b(jin)
  ibound(jin) = -ibound(jin)
  at(kl:n,jin) = -at(kl:n,jin)
  at(n+1,jin) = 2.0E+00 - at(n+1,jin)

  go to 130

40 continue

  k1 = kl
!
!  Part 2 of the algorithm.
!
50 continue

  do i = k1, n+1

    if ( bv(i) <= - eps ) then

      ibv = ibv + 1
      iout = i
      jin = icbas(i)
      b(jin) = - b(jin)
      ibound(jin) = - ibound(jin)
      do l = kl, n+1
        d = bv(l)
        at(l,jin) = d + d - at(l,jin)
      end do
      go to 130

    end if

  end do
!
!  Part 3 of the algorithm.
!
60 continue

  ipart = 3

  do j = 1, m

    if ( any ( icbas(kl:n+1) == j ) ) then

      zc(j) = 0.0E+00

    else

      s = - dble ( b(j) )
      do i = kl, n + 1
        k = icbas(i)
        s = s + dble ( at(i,j) ) * dble ( b(k) )
      end do

      zc(j) = real ( s, kind = 4 )

    end if

  end do

  s = 0.0E+00
  do i = kl, n + 1
    k = icbas(i)
    s = s + dble ( bv(i) ) * dble ( b(k) )
  end do

  res_norm = real ( s, kind = 4 )

  if ( res_norm <= - eps ) then

    do j = 1, m
      b(j) = - b(j)
      ibound(j) = - ibound(j)
      zc(j) = - zc(j)
    end do

    res_norm = - res_norm

    do j = kl, n
      do i = kl, n + 1
        binv(i,j) = - binv(i,j)
      end do
    end do

  end if
!
!  Start of an iteration in part 3 of the algorithm.
!  Determine the vector which enters the basis.
!
80 continue

  ivo = 0
  g = huge ( g )
  tz = 2.0E+00 * res_norm + eps

  do j = 1, m

    do i = kl, n+1
      if ( j == icbas(i) ) then
        go to 110
      end if
    end do

    d = zc(j)

    if ( d < - eps ) then
      go to 90
    end if

    if ( d < tz ) then
      cycle
    end if

    e = tz - d

    if ( g <= e ) then
      cycle
    end if

    ivo = - 1
    g = e
    jin = j
    cycle

90  continue

    e = d

    if ( g <= e ) then
      cycle
    end if

    ivo = 1
    g = e
    jin = j

110 continue

  end do

  if ( ivo == 0 ) then
    go to 150
  end if

  if ( ivo == -1 ) then

    do i = kl, n+1
      d = bv(i)
      at(i,jin) = 2.0E+00 * d - at(i,jin)
    end do

    zc(jin) = 2.0E+00 * res_norm - zc(jin)
    b(jin) = - b(jin)
    ibound(jin) = - ibound(jin)

  end if

  itest = 0
!
!  Determine the vector which leaves the basis.
!
  thmax = huge ( thmax )

  do  i = kl, n+1

    d = at(i,jin)

    if ( d < eps ) then
      cycle
    end if

    g = bv(i) / d

    if ( g <= thmax ) then
      thmax = g
      iout = i
      itest = 1
    end if

  end do

  if ( itest /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ABD_LI - Fatal error!'
    write ( *, '(a)' ) '  Could not find a feasible solution.'
    iflag = 2
    return
  end if
!
!  Performing a Gauss Jordan elimination step.
!
130 continue

  pivot = at(iout,jin)

  at(iout,1:m) = at(iout,1:m) / pivot

  k1 = kl + it

  l = min ( k1, n+1 )

  binv(iout,kl:l) = binv(iout,kl:l) / pivot

  do i = kl, n+1

    if ( i /= iout ) then

      d = at(i,jin)

      if ( eps <= abs ( d ) ) then
        at(i,1:m) = at(i,1:m) - d * at(iout,1:m)
        binv(i,kl:l) = binv(i,kl:l) - d * binv(iout,kl:l)
      end if

    end if

  end do

  do i = kl, n+1
    bv(i) = binv(i,n+1)
  end do

  it = it + 1
  icbas(iout) = jin

  if ( ipart == 3 ) then
    go to 140
  end if

  if ( ipart == 1 .and. ibv == 0 ) then
    go to 10
  end if

  if ( iout == n+1 ) then
    go to 60
  end if

  k1 = iout + 1
  go to 50

140 continue

  d = zc(jin)

  if ( abs ( d ) < eps ) then
    go to 80
  end if

  zc(1:m) = zc(1:m) - d * at(iout,1:m)
  res_norm = res_norm - d * bv(iout)
  go to 80
!
!  Calculate the answer.
!
150 continue

  do j = kl, n+1

    s = 0.0E+00
    do i = kl, n+1
      k = icbas(i)
      s = s + dble ( b(k) ) * dble ( binv(i,j) )
    end do

    k = irbas(j)
!
!  Original code allowed K = N+1.
!
    if ( 1 <= k .and. k <= n ) then
      x(k) = real ( s, kind = 4 )
    end if

  end do

  do j = 1, m
    d = zc(j) - res_norm
    if ( ibound(j) == - 1 ) then
      d = - d
    end if
    r(j) = d
  end do
!
!  Suppress IFLAG = 1 in favor of IFLAG = 0
!
  if ( rank < n ) then
!   iflag = 1
    iflag = 0
    return
  end if

  do i = 1, n+1
    if ( bv(i) < eps ) then
      iflag = 0
!     iflag = 1
      return
    end if
  end do

  return
end
subroutine afk_l1 ( a, m, n, b, eps, x, r, iflag )

!*****************************************************************************80
!
!! AFK_L1 minimizes the L1 norm of A*x-b.
!
!  Discussion:
!
!    The routine uses a modification of the method of A478_L1, a revised
!    simplex method with updating of the LU decomposition.
!
!    It is assumed that A has rank N.
!
!  Modified:
!
!    03 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 68-71,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) EPS, a tolerance used for the accuracy test.
!
!    Output, real ( kind = 4 ) X(N), the computed solution.
!
!    Output, real ( kind = 4 ) R(M), the residuals.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    2, rank(A) < N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) delta(m)
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibase(n)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iin
  integer ( kind = 4 ) index(n)
  integer ( kind = 4 ) inext(m)
  logical intl
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) ipt
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) kkk
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) l
  real ( kind = 4 ) lu(n,n)
  real ( kind = 4 ) price(m)
  real ( kind = 4 ) r(m)
  real ( kind = 4 ) ratio
  real ( kind = 4 ) rho
  real ( kind = 4 ) rhs(n)
  real ( kind = 4 ) sigma(m)
  real ( kind = 4 ) subt
  real ( kind = 4 ) sum2
  real ( kind = 4 ) t
  real ( kind = 4 ) test
  real ( kind = 4 ) tot(n)
  real ( kind = 4 ) val
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) yhat
!
!  Initialize settings.
!
  iflag = 0
  it = 0
!
!  Set up the initial LU decomposition.
!
  tot(1:n) = 0.0E+00

  call i4vec_indicator ( n, index )

  intl = .true.

  k = 1

  call afk_l1_update ( k, intl, a, lu, ibase, index, m, n, eps, iflag2 )

  if ( iflag2 /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AFK_L1 - Fatal error!'
    write ( *, '(a)' ) '  AKF_LI_UPDATE returned nonzero IFLAG2 = ', iflag2
    iflag = 2
    return
  end if

  intl = .false.

  do i = 1, n
    rhs(i) = b(ibase(i))
  end do

  kkk = 0
  call afk_l1_calbet ( kkk, x, rhs, lu, index, n )
!
!  Calculate initial R, TOT and SIGMA.
!
  r(1:m) = b(1:m) - matmul ( a(1:m,1:n), x(1:n) )

  sigma(1:m) = 1.0E+00
  sigma(1:m) = sign ( sigma(1:m), r(1:m) )

  rhs(1:n) = 0.0E+00

  do j = 1, n
    sigma(ibase(j)) = 0.0E+00
  end do

  tot(1:n) = matmul ( sigma(1:m), a(1:m,1:n) )
!
!  The main iterative loop begins.
!
  do

    k = index(1)
    x(1) = tot(k)

    do ii = 2, n
      k = index(ii)
      x(ii) = tot(k)
      do i = 1, ii-1
        x(ii) = x(ii) - lu(k,i) * x(i)
      end do
    end do

    x(n) = x(n) / lu(index(n),n)
    do k1 = n-1, 1, -1
      k = index(k1)
      do k2 = n, k1+1, -1
        x(k1) = x(k1) - lu(k,k2) * x(k2)
      end do
      x(k1) = x(k1) / lu(k,k1)
    end do

    t = 1.0E+00 + eps
    k = 0

    do j = 1, n
      if ( t <= abs ( x(j) ) ) then
        k = j
        t = abs ( x(j) )
        rho = sign ( 1.0E+00, x(j) )
      end if
    end do
!
!  Calculate the optimal X and deviations.
!
    if ( k == 0 ) then

      do i = 1, n
        rhs(i) = b(ibase(i))
      end do

      kkk = 0

      call afk_l1_calbet ( kkk, x, rhs, lu, index, n )

      r(1:m) = b(1:m) - matmul ( a(1:m,1:n), x(1:n) )

      exit

    end if

    kkk = k - 1
    rhs(k) = 1.0E+00

    call afk_l1_calbet ( kkk, x, rhs, lu, index, n )

    rhs(k) = 0.0E+00

    do i = 1, m
      delta(i) = 0.0E+00
      if ( sigma(i) /= 0.0E+00 ) then
        delta(i) = rho * dot_product ( a(i,1:n), x(1:n) )
      end if
    end do
!
!  Perform partial sort of ratios.
!
    t = t * 0.5E+00
    kount = 0
    ratio = huge ( ratio )
    sum2 = 0.5E+00 + eps
    subt = 0.0E+00

    do i = 1, m

      if ( ( delta(i) * sigma(i) ) <= eps ) then
        cycle
      end if

      test = r(i) / delta(i)

      if ( ratio <= test ) then
        cycle
      end if

      sum2 = sum2 + abs ( delta(i) )
!
!  Insert I into list.
!
      if ( ( sum2 - subt ) <= t ) then
        kount = kount + 1
        price(kount) = test
        inext(kount) = i
        cycle
      end if
!
!  Update SUM2 and remove IIN from list.
!
20    continue

      sum2 = sum2 - subt
      ratio = test
      ipt = 0
      kkk = 0
!
!  Identify a new ratio.
!
      do

        kkk = kkk + 1

        if ( kount < kkk ) then
          exit
        end if

        if ( ratio <= price(kkk) ) then
          ratio = price(kkk)
          ipt = kkk
        end if

      end do
!
!  Switch values.
!
      if ( ipt /= 0 ) then

        kkk = inext(ipt)
        subt = abs ( delta(kkk) )

        if ( t <= sum2 - subt ) then
          price(ipt) = price(kount)
          inext(ipt) = inext(kount)
          kount = kount - 1
          go to 20
        end if

        iin = inext(ipt)
        inext(ipt) = i
        price(ipt) = test

      else

        iin = i
        subt = abs ( delta(i) )

      end if

    end do
!
!  Update the basic indicators.
!
    do j = 1, kount
      kkk = inext(j)
      tot(1:n) = tot(1:n) - 2.0E+00 * sigma(kkk) * a(kkk,1:n)
      sigma(kkk) = - sigma(kkk)
    end do

    iout = ibase(k)
    delta(iout) = rho
    tot(1:n) = tot(1:n) - rho * a(iout,1:n) - sigma(iin) * a(iin,1:n)

    sigma(iout) = - rho
    ibase(k) = iin

    call afk_l1_update ( k, intl, a, lu, ibase, index, m, n, eps, iflag2 )

    if ( iflag2 /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AFK_L1 - Fatal error!'
      write ( *, '(a)' ) '  AKF_LI_UPDATE returned nonzero IFLAG2 = ', iflag2
      iflag = 2
      return
    end if

    r(1:m) = r(1:m) - ratio * delta(1:m)
    sigma(iin) = 0.0E+00
!
!  Suppress IFLAG = 1 return.
!
    if ( sum2 - t < eps ) then
      iflag = 0
!     iflag = 1
    else
      iflag = 0
    end if

    it = it + 1

  end do

  return
end
subroutine afk_l1_calbet ( kkk, x, rhs, lu, index, n )

!*****************************************************************************80
!
!! AFK_L1_CALBET is used by AFK_L1 to compute the solution.
!
!  Modified:
!
!    15 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 71,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) KKK.  On input, if KKK is nonzero,
!    it indicates the number of variables which can be immediately set to 0.
!
!    Output, real ( kind = 4 ) X(N), the solution.
!
!    Input, real ( kind = 4 ) RHS(N), the right hand side vector.
!
!    Input, real ( kind = 4 ) LU(N,N), ?
!
!    Input, integer ( kind = 4 ) INDEX(N), ?
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) index(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kkk
  real ( kind = 4 ) lu(n,n)
  real ( kind = 4 ) rhs(n)
  real ( kind = 4 ) x(n)

  if ( kkk /= 0 ) then

    do i = 1, kkk
      k = index(i)
      x(k) = 0.0E+00
    end do

    kkk = kkk + 1

  else

    k = index(1)
    x(k) = rhs(1) / lu(k,1)
    kkk = 2

  end if
!
!  KKK now represents the index of the first variable that has
!  not been solved for yet.
!
  do ii = kkk, n
    k = index(ii)
    x(k) = rhs(ii)
    do i = 1, ii-1
      kk = index(i)
      x(k) = x(k) - x(kk) * lu(kk,ii)
    end do
    x(k) = x(k) / lu(k,ii)
  end do

  do k1 = n-1, 1, -1

    k = index(k1)

    do kk = n, k1+1, -1
      k2 = index(kk)
      x(k) = x(k) - x(k2) * lu(k2,k1)
    end do

  end do

  return
end
subroutine afk_l1_update ( kkk, intl, a, lu, ibase, index, m, n, eps, iflag )

!*****************************************************************************80
!
!! AFK_L1_UPDATE updates the LU decomposition of a matrix.
!
!  Modified:
!
!    24 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 71-72,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) KKK.  On input, KKK indicates that pivot
!    rows are needed for entries KKK to N.  On output, KKK is the
!    label of the last row used.
!
!    Input, logical INTL, basis indicator.  If INTL is FALSE, then
!    IBASE contains the row indices to use.  Otherwise, the routine
!    should simply use rows consecutively, and set up IBASE.
!
!    Input, real ( kind = 4 ) A(M,N), the matrix to be factored.
!
!    Input/output, real ( kind = 4 ) LU(N,N), ?
!
!    Input/output, integer ( kind = 4 ) IBASE(N), a list of rows in the order they
!    should be used.  If INTL is TRUE, then IBASE is an output quantity
!    created by this routine.  If INTL is FALSE, then IBASE is an input
!    quantity, and is used by this routine.
!
!    Input, integer ( kind = 4 ) INDEX(N), ?
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for the size of the pivot.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    1, an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibase(n)
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) index(n)
  logical intl
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kkk
  real ( kind = 4 ) lu(n,n)
  real ( kind = 4 ) pivot
  real ( kind = 4 ) subt

  iflag = 0
  irow = 0

  do ii = kkk, n

    if ( intl ) then
      irow = irow + 1
    else
      irow = ibase(ii)
    end if

    do

      if ( m < irow ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'AFK_L1_UPDATE - Fatal error!'
        write ( *, '(a)' ) '  Reached last row.'
        iflag = 1
        return
      end if

      lu(1:n,ii) = a(irow,1:n)
!
!  Set up representation of incoming row.
!
      do icol = 1, ii-1
        k = index(icol)
        subt = lu(k,ii)
        do i = icol+1, n
          k = index(i)
          lu(k,ii) = lu(k,ii) - subt * lu(k,icol)
        end do
      end do
!
!  Find the maximum entry.
!
      pivot = eps
      kk = 0

      do i = ii, n
        k = index(i)
        if ( pivot < abs ( lu(k,ii) ) ) then
          pivot = abs ( lu(k,ii) )
          kk = i
        end if
      end do

      if ( kk /= 0 ) then
        exit
      end if

      irow = irow + 1

    end do
!
!  Swap order.
!
    call i4_swap ( index(kk), index(ii) )
!
!  Put in columns of LU one at a time.
!
    if ( intl ) then
      ibase(ii) = irow
    end if

    do i = ii+1, n
      k = index(i)
      lu(k,ii) = lu(k,ii) / lu(index(ii),ii)
    end do

  end do

  kkk = irow

  return
end
subroutine avllsq ( a, m, n, b, m2, s, w, itmax, eps1, eps2, eps3, x, &
  fx, iflag )

!*****************************************************************************80
!
!! AVLLSQ carries out average linear regression.
!
!  Discussion:
!
!    We assume that we have a set of M observations of N dimensional
!    data, but that these observations have been grouped into S
!    clusters.  We assume cluster J has a total population of M2(J),
!    and has an associated submatrix A2(J) and sub righthand side B2(J).
!    Moreover, each cluster J has been assigned a weight W(J).
!
!    We seek a single N vector X which minimizes the objective function
!
!      F(X) = sum ( 1 <= J <= S )
!        W(J) * || A2(J) * X - B2(J) ||
!
!    where the norm used is the L2 norm.
!
!    Note that the norm is NOT squared in this objective function.
!
!    On input, the rows of the system matrix A and the entries of the
!    right hand side B must be ordered to correspond to the individual
!    clusters.  That is, rows 1 through M2(1) correspond to cluster 1,
!    M2(1)+1 through M2(1)+M2(2) to cluster 2, and so on.
!
!  Modified:
!
!    13 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 189-190,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the total number of observations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, integer ( kind = 4 ) M2(S), the size of each cluster.  The M2's must
!    sum up to M.
!
!    Input, integer ( kind = 4 ) S, the number of clusters.
!
!    Input, real ( kind = 4 ) W(S), the weight coefficients in the objective 
!    function.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations to take.
!
!    Input, real ( kind = 4 ) EPS1, a tolerance to be used by MGS_L2.
!
!    Input, real ( kind = 4 ) EPS2, a convergence tolerance.
!
!    Input, real ( kind = 4 ) EPS3, a singularity tolerance.
!
!    Output, real ( kind = 4 ) X(N), the solution.
!
!    Output, real ( kind = 4 ) FX, the value of the objective function on the
!    penultimate iteration.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no errors detected.
!    2, the individual values M2 sum to more than M.
!    3, the system is near singular.
!    4, failure in MGS_L2.
!    5, some M2 is less than N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) aa(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) bb(m)
  logical bnew
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  real ( kind = 4 ) eps3
  real ( kind = 4 ) f
  real ( kind = 4 ) fd
  real ( kind = 4 ) fx
  real ( kind = 4 ) g(s)
  real ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m2(s)
  integer ( kind = 4 ) ms
  integer ( kind = 4 ) r
  real ( kind = 4 ) w(s)
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)
!
!  Step 1: Initialization.
!
  iflag = 0
  fx = 0.0E+00
  it = 0
  x(1:n) = 0.0E+00

  do j = 1, s

    if ( m2(j) <= n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AVLLSQ - Fatal error!'
      write ( *, '(a)' ) '  M2(J) <= N for some J.'
      iflag = 5
      return
    end if

  end do

  if ( m < sum ( m2(1:s) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVLLSQ - Fatal error!'
    write ( *, '(a)' ) '  The individual values M2 sum to more than M.'
    iflag = 2
    return
  end if

  g(1:s) = sqrt ( w(1:s) )

  do it = 1, itmax
!
!  Step 2: generation and solving of the weighted linear least squares problem.
!
    aa(1:m,1:n) = a(1:m,1:n)
    bb(1:m) = b(1:m)

    ms = 0

    do j = 1, s

      do i = 1, m2(j)
        r = ms + i
        aa(r,1:n) = g(j) * aa(r,1:n)
        bb(r) = g(j) * bb(r)
      end do

      ms = ms + m2(j)

    end do

    bnew = .false.

    call mgs_l2 ( aa, m, n, bb, eps1, bnew, y, iflag2 )

    if ( iflag2 /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AVLLSQ - Fatal error!'
      write ( *, '(a)' ) '  MGS_L2 returned error flag = ', iflag2
      iflag = 4
      return
    end if

    h = sum ( abs ( y(1:n) - x(1:n) ) )

    x(1:n) = y(1:n)
!
!  Check for convergence.
!
    if ( h < eps2 * sum ( abs ( y(1:n) ) ) ) then
      exit
    end if
!
!  Step 3: calculate the objective function and adjust weights G.
!
    fx = 0.0E+00
    ms = 0

    do j = 1, s

      h = 0.0E+00

      do i = 1, m2(j)

        r = ms + i
        f = b(r) - dot_product ( a(r,1:n), x(1:n) )
        h = h + f**2

      end do
!
!  It almost looks like we're being penalized for a small residual here!
!  The text claims we're checking the size of A*X, but we're clearly
!  checking B-A*X.
!
      if ( h < eps3 ) then
        iflag = 3
        exit
      end if

      h = sqrt ( h )
      fx = fx + w(j) * h
      g(j) = sqrt ( w(j) / h )
      ms = ms + m2(j)

    end do

  end do

  return
end
subroutine blod_l1 ( a, m, n, b, x, res_norm )

!*****************************************************************************80
!
!! BLOD_L1 minimizes the L1 norm of A*x-b.
!
!  Discussion:
!
!    The routine performs fast least-absolute-deviations fitting.
!
!    It is assumed that the rank of A = N < M.
!
!  Modified:
!
!    11 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 75-79,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Output, real ( kind = 4 ) X(N), the computed solution.
!
!    Output, real ( kind = 4 ) RES_NORM, the L1 norm of the residual.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) c
  real ( kind = 4 ) cinit
  real ( kind = 4 ) cmax
  real ( kind = 4 ) cn(m+n)
  real ( kind = 4 ) con
  real ( kind = 4 ) cp(m+n)
  real ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) l(m+n)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) list(m+n)
  integer ( kind = 4 ) lrow
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nzr
  real ( kind = 4 ) r(m+n)
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) shi
  real ( kind = 4 ) slo
  real ( kind = 4 ) sn
  real ( kind = 4 ) sp
  real ( kind = 4 ) sz
  real ( kind = 4 ) u
  real ( kind = 4 ) w(m+n)
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) z(m+n,n+1)
!
!  Convenient constants.
!
  nj = 0
  last = 0
!
!  Copy A, B into Z.
!
  cinit = sum ( b(1:m) )
!
!  Z = ( A | B )
!      ( I | 0 )
!
  z(1:m,1:n) = a(1:m,1:n)
  z(1:m,n+1) = b(1:m)
  z(m+1:m+n,1:n+1) = 0.0
  do i = 1, n
    z(m+i,i) = 1.0E+00
  end do

  it = 0
  l(1:n) = 0
  j = 0
!
!  Top of main loop.
!
!  Find pivot column.
!
  do

    it = it + 1

    cmax = 0.0E+00
    shi = 0.0E+00
    slo = 0.0E+00
    j = last
    jcol = 0

    do j1 = 1, n

      if ( j1 /= j ) then

        call blod_l1_crit ( m, z(1,j1), z(1,n+1), sn, sz, sp )

        c = ( abs ( sp - sn ) - sz )
        d = sp + sn + sz
        u = c / d
        cp(j1) = u
        cn(j1) = c

        if ( cmax < u ) then
          cmax = u
          shi = sp + sz
          slo = sn + sz
          jcol = j1
        end if

      end if

    end do
!
!  Test for convergence.
!
    if ( jcol <= 0 ) then
      exit
    end if

    nj = nj + 1
    last = jcol
!
!  Find LROW, the pivot row.
!
    j = jcol

    if ( slo < shi ) then
      shi = 0.0E+00
      call blod_l1_get1 ( z(1,j), z(1,n+1), m, r, w, list, i1 )
    else
      slo = 0.0E+00
      call blod_l1_get2 ( z(1,j), z(1,n+1), m, r, w, list, i1 )
    end if

    call blod_l1_med3 ( r, w, list, i1, slo, shi, lrow )

    l(j) = lrow
!
!  Update Z by pivoting.
!
    z(1:m+n,j) = z(1:m+n,j) / z(lrow,j)
    z(lrow,j) = 1.0E+00

    do j1 = 1, n+1
      if ( j1 /= j ) then
        z(1:m+n,j1) = z(1:m+n,j1) - z(lrow,j1) * z(1:m+n,j)
        z(lrow,j1) = 0.0E+00
      end if
    end do

  end do
!
!  Determine the value of the objective function.
!
  res_norm = sum ( abs ( z(1:m,n+1) ) )
!
!  Extract the solution.
!
  x(1:n) = -z(m+1:m+n,n+1)

  return
end
subroutine blod_l1_crit ( n, u, v, sn, sz, sp )

!*****************************************************************************80
!
!! BLOD_L1_CRIT returns three sums of entries of U.
!
!  Discussion:
!
!    The routine is used by BLOD_L1.
!
!  Modified:
!
!    03 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 77,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the vectors.
!
!    Input, real ( kind = 4 ) U(N), V(N), the two vectors.
!
!    Output, real ( kind = 4 ) SN, SZ, SP, the sum of the absolute values of U(I)
!    for which U(I) * V(I) is negative, zero, or positive respectively.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 4 ) sn
  real ( kind = 4 ) sp
  real ( kind = 4 ) sz
  real ( kind = 4 ) u(n)
  real ( kind = 4 ) v(n)

  sn = 0.0E+00
  sz = 0.0E+00
  sp = 0.0E+00

  do i = 1, n

    if ( u(i) * v(i) < 0.0E+00 ) then
      sn = sn + abs ( u(i) )
    else if ( u(i) * v(i) == 0.0E+00 ) then
      sz = sz + abs ( u(i) )
    else if ( 0.0E+00 < u(i) * v(i) ) then
      sp = sp + abs ( u(i) )
    end if

  end do

  return
end
subroutine blod_l1_get1 ( u, v, n, r, w, list, i1 )

!*****************************************************************************80
!
!! BLOD_L1_GET1 is used by BLOD_L1.
!
!  Modified:
!
!    03 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 76,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) U(N), ?
!
!    Input, real ( kind = 4 ) V(N), ?
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Output, real ( kind = 4 ) R(N), ?
!
!    Output, real ( kind = 4 ) W(N), ?
!
!    Output, integer ( kind = 4 ) LIST(*), ?
!
!    Output, integer ( kind = 4 ) I1, the number of entries inserted into LIST.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) list(*)
  real ( kind = 4 ) r(n)
  real ( kind = 4 ) u(n)
  real ( kind = 4 ) v(n)
  real ( kind = 4 ) w(n)

  i1 = 0

  do i = 1, n

    if ( ( 0.0E+00 < u(i) .and. 0.0E+00 < v(i) ) .or. &
         ( u(i) < 0.0E+00 .and. v(i) < 0.0E+00 ) ) then

      r(i) = v(i) / u(i)
      w(i) = abs ( u(i) )
      i1 = i1 + 1
      list(i1) = i

    end if

  end do

  return
end
subroutine blod_l1_get2 ( u, v, n, r, w, list, i1 )

!*****************************************************************************80
!
!! BLOD_L1_GET2 is used by BLOD_L1.
!
!  Modified:
!
!    03 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 77,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) U(N), ?
!
!    Input, real ( kind = 4 ) V(N), ?
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Output, real ( kind = 4 ) R(N), ?
!
!    Output, real ( kind = 4 ) W(N), ?
!
!    Output, integer ( kind = 4 ) LIST(*), ?
!
!    Output, integer ( kind = 4 ) I1, the number of entries inserted into LIST.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) list(*)
  real ( kind = 4 ) r(n)
  real ( kind = 4 ) u(n)
  real ( kind = 4 ) v(n)
  real ( kind = 4 ) w(n)

  i1 = 0

  do i = 1, n

    if ( ( 0.0E+00 < u(i) .and. v(i) < 0.0E+00 ) .or. &
         ( u(i) < 0.0E+00 .and. 0.0E+00 < v(i) ) ) then

      r(i) = v(i) / u(i)
      w(i) = abs ( u(i) )
      i1 = i1 + 1
      list(i1)  = i

    end if

  end do

  return
end
subroutine blod_l1_med3 ( r, w, list, n, sslo, sshi, med3 )

!*****************************************************************************80
!
!! BLOD_L1_MED3 is used by BLOD_L1.
!
!  Modified:
!
!    30 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 77-78,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) R(*), ?
!
!    Input, real ( kind = 4 ) W(*), ?
!
!    Input, integer ( kind = 4 ) LIST(N), ?
!
!    Input, integer ( kind = 4 ) N, the number of items in LIST.
!
!    Input, real ( kind = 4 ) SSLO, SSHI, ?
!
!    Output, integer ( kind = 4 ) MED3, ?
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) med3
  integer ( kind = 4 ) mid
  real ( kind = 4 ) r(*)
  real ( kind = 4 ) shi
  real ( kind = 4 ) slo
  real ( kind = 4 ) sshi
  real ( kind = 4 ) sslo
  real ( kind = 4 ) test
  real ( kind = 4 ) thi
  real ( kind = 4 ) tlo
  real ( kind = 4 ) w(*)
  real ( kind = 4 ) xt

  lo = 1
  hi = n
  slo = sslo
  shi = sshi

  do

    if ( hi - lo <= 1 ) then

      med3 = list(lo)

      if ( lo == hi ) then
        exit
      end if

      if ( r(list(hi)) < r(list(lo)) ) then
        call i4_swap ( list(lo), list(hi) )
        med3 = list(lo)
      end if

      if ( slo + w(list(lo)) < shi + w(list(hi)) ) then
        med3 = list(hi)
      end if

      exit

    end if

    mid = ( lo + hi ) / 2

    call i4_swap ( list(mid), list(lo+1) )

    if ( r(list(hi)) < r(list(lo+1)) ) then
      call i4_swap ( list(lo+1), list(hi) )
    end if

    if ( r(list(hi)) < r(list(lo)) ) then
      call i4_swap ( list(lo), list(hi) )
    end if

    if ( r(list(lo)) < r(list(lo+1)) ) then
      call i4_swap ( list(lo+1), list(lo) )
    end if

    med3 = list(lo)
    i = lo + 1
    j = hi
    xt = r(med3)
    tlo = slo
    thi = shi

    do

      tlo = tlo + w(list(i))
      i = i + 1

      if ( r(list(i)) < xt ) then
        cycle
      end if

      do

        thi = thi + w(list(j))
        j = j - 1

        if ( r(list(j)) <= xt ) then
          exit
        end if

      end do

      if ( j <= i ) then
        exit
      end if

      call i4_swap ( list(i), list(j) )

    end do

    test = w(med3)

    if ( i == j ) then
      test = test + w(list(i))
      i = i + 1
      j = j - 1
    end if

    if ( abs ( thi - tlo ) <= test ) then
      exit
    end if

    if ( tlo <= thi ) then
      slo = tlo + test
      lo = i
    else
      shi = thi + test
      lo = lo + 1
      hi = j
    end if

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
!  Examples:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Modified:
!
!    14 August 1999
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

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
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

    digit = - 1

  end if

  return
end
subroutine c01m ( n, m, d )

!*****************************************************************************80
!
!! C01M generates a new combination from an old one efficiently.
!
!  Discussion:
!
!    The routine using C01M needs to generate every one of a set
!    of combinations of M things out of N.  If C01M is given one
!    of these combinations as input, it can produce the "next"
!    combination as output.  Thus, the calling routine can rely
!    on C01M to produce all the combinations, one after another.
!
!    Before the first call, set D(1:N) = 0 to signal C01M to begin.
!    The last combination to be returned will consist of M 0's
!    followed by N-M 1's.
!
!    Moreover, C01M tries to order the combinations efficiently,
!    so that the next combination differs from the previous one
!    only by one exchange of items.  However, C01M is not always
!    able to do this.  (There are other routines that CAN generate
!    all the subsets using single exchanges, a method sometimes
!    called the "revolving door" method.)
!
!  Modified:
!
!    25 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 141,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input, integer ( kind = 4 ) M, the size of the subset.
!
!    Input/output, integer ( kind = 4 ) D(N).  On input, describes the
!    previous subset.  On output, describes the next subset.
!    D(I) is 1 if I is NOT in the subset, and 0 if it is.
!    (I have no explanation for this counterintuitive usage!)
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) jn
  integer ( kind = 4 ) m
!
!  Search for the last 1 in D, and store that index in J.
!
  j = n

  do while ( d(j) /= 1 )

    j = j - 1

    if ( j == 0 ) then
      d(m+1:n) = 1
      return
    end if

  end do
!
!  Record in JN the length of the final string 1000...000 in D.
!
  jn = n - j + 1
!
!  Now let J record the index of the last 0 preceding the last 1.
!  Record in J1 the number of "extra" 1's you had to skip over (except for
!  the 1 we know about).
!
  j1 = 0

  do

    j = j - 1

    if ( j == 0 ) then
      d(1:m) = 0
      d(m+1:n) = 1
      return
    end if

    if ( d(j) == 0 ) then
      exit
    end if

   j1 = j1 + 1

  end do
!
!  So now we know our pattern is
!
!    0 : J1+1 1's : JN-1 0's
!
!  Replace this by
!
!    1 : JN 0's : J1 1's
!
  d(j) = 1
  d(j+1:j+jn) = 0
  d(j+jn+1:j+jn+j1) = 1

  return
end
subroutine con_l1 ( m, l, k, n, a, b, c, d, e, f, code, eps, itmax, x, &
  res, res_norm, iflag )

!*****************************************************************************80
!
!! CON_L1 minimizes the L1 norm of A * X - B subject to linear constraints.
!
!  Discussion:
!
!    The constraints have the form:
!
!      C * X  = D
!      E * X >= F
!
!    The entire system is packed into a single matrix Q with an extra
!    column and two extra rows:
!
!         N 1 1
!
!      M  A B *
!      L  C D *
!      K  E F *
!      1  * * *
!      1  * * *
!
!    Note that the original code solves E * X <= F.  For consistency,
!    I have used E * X >= F, to accord with other routines.  It only
!    took me a week to track down this source of confusion...
!
!  Modified:
!
!    16 May 2002
!
!  Reference:
!
!    I Barrodale, F Roberts,
!    Algorithm 552: Solution of the Constrained L1 Linear Approximation Problem,
!    ACM Transactions on Mathematical Software,
!    Volume 6, pages 231-235, 1980.
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 230-234,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) L, the number of rows of C.
!
!    Input, integer ( kind = 4 ) K, the number of rows of E.
!
!    Input, integer ( kind = 4 ) N, the number of unknowns, and the number of
!    columns in A, C and E.
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) C(L,N), the equality constraint matrix.
!
!    Input, real ( kind = 4 ) D(L), the equality right hand side.
!
!    Input, real ( kind = 4 ) E(K,N), the inequality constraint matrix.
!
!    Input, real ( kind = 4 ) F(K), the inequality right hand side.
!
!    Input, integer ( kind = 4 ) CODE, indicates additional sign constraints that
!    may be placed on the solution X, and on the sign of A*X-B.
!    0, no sign constraints.
!    1, the input value of X indicates the sign constraint on X:
!      X(J) == -1.0, then X(J) must be negative.
!      X(J) == +1.0, then X(J) must be positive.
!      X(J) ==  0.0, there is no constraint on X(J).
!      the input value of RES indicates the sign constraint on X:
!      RES(J) == -1.0, then (B-A*X)(J) must be negative.
!      RES(J) == +1.0, then (B-A*X)(J) must be positive.
!      RES(J) ==  0.0, there is no constraint on (B-A*X)(J).
!
!    Input, real ( kind = 4 ) EPS, a tolerance for an accuracy test.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations to take.
!
!    Input/output, real ( kind = 4 ) X(N).  Used for input if CODE = 1.  On output,
!    contains the solution.
!
!    Input/output, real ( kind = 4 ) RES(M+L+K).  Entries 1 through M used for 
!    input if CODE = 1.  On output, contains the residuals.
!
!    Output, real ( kind = 4 ) RES_NORM, the L1 norm of the residual.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, an optimal solution was found.
!    1, there is no feasible solution.
!    2, premature exit because of rounding errors.
!    3, maximal number of iterations taken.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) c(l,n)
  integer ( kind = 4 ) code
  real ( kind = 4 ) cu(2,n+m+l+k)
  real ( kind = 4 ) cuv
  real ( kind = 4 ) d(l)
  real ( kind = 4 ) e(k,n)
  real ( kind = 4 ) eps
  real ( kind = 4 ) f(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) iphase
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) iu(2,n+m+l+k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) js
  integer ( kind = 4 ) kforce
  integer ( kind = 4 ) kk
  real ( kind = 4 ) pivot
  real ( kind = 4 ) q(m+l+k+2,n+2)
  real ( kind = 4 ) res(m+l+k)
  real ( kind = 4 ) res_norm
  integer ( kind = 4 ) s(m+l+k)
  real ( kind = 4 ) sn
  real ( kind = 4 ) temp
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) z
  real ( kind = 4 ) zu
  real ( kind = 4 ) zv
!
!  Initialization.
!
  iflag = 0
  kforce = 1
  it = 0
  js = 1
  ia = 0

  q(1:m,1:n) = a(1:m,1:n)
  q(1:m,n+1) = b(1:m)

  q(m+1:m+l,1:n) = c(1:l,1:n)
  q(m+1:m+l,n+1) = d(1:l)
!
!  In order to solve E * X >= F, we need to put -F into Q.
!  The original code solved E * X <= F...
!
  q(m+l+1:m+l+k,1:n) = e(1:k,1:n)
  q(m+l+1:m+l+k,n+1) = -f(1:k)
!
!  Set up labels in Q.
!
  call r4vec_indicator ( n, q(m+l+k+2,1:n) )

  do i = 1, m+l+k
    q(i,n+2) = real ( n + i, kind = 4 )
    if ( q(i,n+1) < 0.0E+00 ) then
      q(i,1:n+2) = -q(i,1:n+2)
    end if
  end do
!
!  Set up phase 1 costs.
!
  iphase = 2

  cu(1:2,1:n+m+l+k) = 0.0E+00
  iu(1:2,1:n+m+l+k) = 0

  if ( l /= 0 ) then
    cu(1:2,n+m+1:n+m+l) = 1.0E+00
    iu(1:2,n+m+1:n+m+l) = 1
    iphase = 1
  end if

  if ( k /= 0 ) then

    do j = n+m+l+1, n+m+l+k
      cu(2,j) = 1.0E+00
      iu(2,j) = 1
      if ( q(j-n,n+2) < 0.0E+00 ) then
        iphase = 1
      end if
    end do

  end if
!
!  The user may have input some sign conditions on X.
!
  if ( code /= 0 ) then

    do j = 1, n

      if ( x(j) < 0.0E+00 ) then
        cu(1,j) = 1.0E+00
        iu(1,j) = 1
      else if ( 0.0E+00 < x(j) ) then
        cu(2,j) = 1.0E+00
        iu(2,j) = 1
      end if

    end do

    do j = 1, m

      if ( res(j) < 0.0E+00 ) then
        cu(1,n+j) = 1.0E+00
        iu(1,n+j) = 1
        if ( 0.0E+00 < q(j,n+2) ) then
          iphase = 1
        end if
      else if ( 0.0E+00 < res(j) ) then
        cu(2,n+j) = 1.0E+00
        iu(2,n+j) = 1
        if ( q(j,n+2) < 0.0E+00 ) then
          iphase = 1
        end if
      end if

    end do

  end if

  if ( iphase == 2 ) then
    go to 90
  end if
!
!  Compute the marginal costs.
!
10 continue

  do j = js, n+1

    temp = 0.0E+00

    do i = 1, m + l + k

      ii = q(i,n+2)

      if ( 0 <= ii ) then
        z = cu(1,ii)
      else
        z = cu(2,-ii)
      end if

      temp = temp + dble ( q(i,j) ) * dble ( z )

    end do

    q(m+l+k+1,j) = temp

  end do

  do j = js, n

    ii = q(m+l+k+2,j)

    if ( 0 <= ii ) then
      z = cu(1,ii)
    else
      z = cu(2,-ii)
    end if

    q(m+l+k+1,j) = q(m+l+k+1,j) - z

  end do
!
!  Determine the vector to enter the basis.
!
20 continue

  xmax = 0.0E+00

  if ( n < js ) then
    go to 80
  end if

  do j = js, n

    zu = q(m+l+k+1,j)
    ii = q(m+l+k+2,j)

    if ( ii <= 0 ) then

      ii = -ii
      zv = zu
      zu = - zu - cu(1,ii) - cu(2,ii)

    else

      zv = - zu - cu(1,ii) - cu(2,ii)

    end if

    if ( kforce /= 1 .or. ii <= n ) then

      if ( iu(1,ii) /= 1 ) then
        if ( xmax < zu ) then
          xmax = zu
          in = j
        end if
      end if

      if ( iu(2,ii) /= 1 ) then
        if ( xmax < zv ) then
          xmax = zv
          in = j
        end if
      end if

    end if

  end do

  if ( xmax <= eps ) then
    go to 80
  end if

  if ( q(m+l+k+1,in) /= xmax ) then

    q(1:m+l+k+2,in) = -q(1:m+l+k+2,in)
    q(m+l+k+1,in) = xmax

  end if
!
!  Determine the vector to leave the basis.
!
  if ( iphase /= 1 .and. ia /= 0 ) then

    xmax = 0.0E+00

    do i = 1, ia
      z = abs ( q(i,in) )
      if ( xmax < z ) then
        xmax = z
        iout = i
      end if
    end do

    if ( eps < xmax ) then

      do j = 1, n + 2
        call r4_swap ( q(ia,j), q(iout,j) )
      end do

      iout = ia
      ia = ia - 1
      pivot = q(iout,in)
      go to 70

    end if

  end if

  kk = 0
  do i = 1, m + l + k
    z = q(i,in)
    if ( eps < z ) then
      kk = kk + 1
      res(kk) = q(i,n+1) / z
      s(kk) = i
    end if
  end do

40 continue

  if ( kk == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_L1 - Fatal error!'
    write ( *, '(a)' ) '  Premature exit because of rounding errors.'
    iflag = 2
    go to 120
  end if

  xmin = res(1)
  iout = s(1)
  j = 1

  if ( kk /= 1 ) then

    do i = 2, kk
      if ( res(i) < xmin ) then
        j = i
        xmin = res(i)
        iout = s(i)
      end if
    end do

    res(j) = res(kk)
    s(j) = s(kk)

  end if

  kk = kk - 1
  pivot = q(iout,in)
  ii = q(iout,n+2)

  if ( iphase == 1 ) then
    go to 60
  end if

  if ( ii < 0 ) then
    go to 50
  end if

  if ( iu(2,ii) == 1 ) then
    go to 70
  end if

  go to 60

50 continue

  if ( iu(1,-ii) == 1 ) then
    go to 70
  end if

60 continue

  ii = abs ( ii )
  cuv = cu(1,ii) + cu(2,ii)

  if ( q(m+l+k+1,in) - pivot * cuv <= eps ) then
    go to 70
  end if
!
!  Bypass the intermediate vertices.
!
  do j = js, n+1
    z = q(iout,j)
    q(m+l+k+1,j) = q(m+l+k+1,j) - z * cuv
    q(iout,j) = -z
  end do

  q(iout,n+2) = -q(iout,n+2)
  go to 40
!
!  Gauss-Jordan elimination.
!
70 continue

  if ( itmax <= it ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_L1 - Fatal error!'
    write ( *, '(a)' ) '  Number of iterations exceeded.'
    iflag = 3
    go to 120
  end if

  it = it + 1

  do j = js, n+1
    if ( j /= in ) then
      q(iout,j) = q(iout,j) / pivot
    end if
  end do

  do j = js, n+1
    if ( j /= in ) then
      z = -q(iout,j)
      do i = 1, m+l+k+1
        if ( i /= iout ) then
          q(i,j) = q(i,j) + z * q(i,in)
        end if
      end do
    end if
  end do

  do i = 1, m+l+k+1
    if ( i /= iout ) then
      q(i,in) = - q(i,in) / pivot
    end if
  end do

  q(iout,in) = 1.0E+00 / pivot

  call r4_swap ( q(iout,n+2), q(m+l+k+2,in) )

  ii = abs ( q(m+l+k+2,in) )

  if ( iu(1,ii) /= 0 .and. iu(2,ii) /= 0 ) then

    do i = 1, m + l + k + 2
      call r4_swap ( q(i,in), q(i,js) )
    end do

    js = js + 1

  end if

  go to 20
!
!  Test for optimality.
!
80 continue

  if ( kforce == 0 ) then
    go to 110
  end if

  if ( iphase /= 1 .or. eps < q(m+l+k+1,n+1) ) then
    kforce = 0
    go to 20
  end if
!
!  Set up phase 2 costs.
!
90 continue

  iphase = 2

  cu(1:2,1:n) = 0.0E+00
  cu(1:2,n+1:n+m) = 1.0E+00
  cu(1:2,n+m+1:n+m+l+k) = 0.0E+00

  do i = 1, m + l + k

    ii = q(i,n+2)

    if ( ii <= 0 ) then
      ii = -ii
      if ( iu(2,ii) == 0 ) then
        cycle
      end if
      cu(2,ii) = 0.0E+00
    else
      if ( iu(1,ii) == 0 ) then
        cycle
      end if
      cu(1,ii) = 0.0E+00
    end if

    ia = ia + 1

    do j = 1, n + 2
      call r4_swap ( q(ia,j), q(i,j) )
    end do

  end do

  go to 10

100 continue

  if ( q(m+l+k+1,n+1) <= eps ) then
    go to 90
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CON_L1 - Fatal error!'
  write ( *, '(a)' ) '  There is no feasible solution.'
  iflag = 1
  go to 120

110 continue

  if ( iphase == 1 ) then
    go to 100
  end if
!
!  Prepare output.
!
120 continue

  res_norm = 0.0E+00

  x(1:n) = 0.0E+00
  res(1:m+l+k) = 0.0E+00

  do i = 1, m+l+k

    ii = q(i,n+2)
    sn = 1.0E+00

    if ( ii <= 0 ) then
      ii = -ii
      sn = -1.0E+00
    end if

    if ( ii <= n ) then

      x(ii) = sn * q(i,n+1)

    else

      res(ii-n) = sn * q(i,n+1)

      if ( n+1 <= ii .and. ii <= n + m ) then
        res_norm = res_norm + q(i,n+1)
      end if

    end if

  end do

  return
end
subroutine con_l2 ( m, l, k, n, a, b, c, d, e, f, eps, x, rank1, rank2, &
  res_norm, iflag )

!*****************************************************************************80
!
!! CON_L2 minimizes the L2 norm of A * X - B subject to linear constraints.
!
!  Discussion:
!
!    THIS ROUTINE IS NOT WORKING.
!
!
!    The constraints have the form:
!
!      C * X  = D
!      E * X >= F
!
!  Modified:
!
!    17 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 222-225,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) L, the number of rows of C.
!
!    Input, integer ( kind = 4 ) K, the number of rows of E.
!
!    Input, integer ( kind = 4 ) N, the number of unknowns, and the number of
!    columns in A, C and E.
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) C(L,N), the equality constraint matrix.
!
!    Input, real ( kind = 4 ) D(L), the equality right hand side.
!
!    Input, real ( kind = 4 ) E(K,N), the inequality constraint matrix.
!
!    Input, real ( kind = 4 ) F(K), the inequality right hand side.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for an accuracy test.
!
!    Output, real ( kind = 4 ) X(N), the computed solution.
!
!    Output, integer ( kind = 4 ) RANK1, the calculated rank of matrix C.
!
!    Output, integer ( kind = 4 ) RANK2, the calculated rank of matrix A.
!
!    Output, real ( kind = 4 ) RES_NORM, the L2 norm of the residual.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    1, illegal value of K, L, M, or N.
!    2, iteration limit exceeded in LDP_L2.
!    3, no feasible point for E * X >= F.
!    4, no feasible point for C * X = D.
!    5, error return from LDP_L2.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) c(l,n)
  real ( kind = 4 ) d(l)
  real ( kind = 4 ) e(k,n)
  real ( kind = 4 ) eca(k+l+m,n)
  real ( kind = 4 ) f(k)
  real ( kind = 4 ) fdb(k+l+m)
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) ip(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) mdb
  integer ( kind = 4 ) mm2
  real ( kind = 4 ) r
  integer ( kind = 4 ) rank1
  integer ( kind = 4 ) rank2
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) tau
  real ( kind = 4 ) w((n+1)*(k+2)+2*k+n)
  real ( kind = 4 ) x(n)

  res_norm = 0.0E+00
  iflag = 0
  rank1 = 0
  rank2 = 0

  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_L2 - Fatal error!'
    write ( *, '(a)' ) '  K < 0.'
    iflag = 1
    return
  end if

  if ( l < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_L2 - Fatal error!'
    write ( *, '(a)' ) '  L < 0.'
    iflag = 1
    return
  end if

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_L2 - Fatal error!'
    write ( *, '(a)' ) '  M < 0.'
    iflag = 1
    return
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_L2 - Fatal error!'
    write ( *, '(a)' ) '  N < 0.'
    iflag = 1
    return
  end if
!
!            N
!        K ( E )
!  ECA = L ( C )
!        M ( A )
!
  eca(1:k,1:n) = e(1:k,1:n)
  eca(k+1:k+l,1:n) = c(1:l,1:n)
  eca(k+l+1:k+l+m,1:n) = a(1:m,1:n)

  fdb(1:k) = f(1:k)
  fdb(k+1:k+l) = d(1:l)
  fdb(k+l+1:k+l+m) = b(1:m)
!
!  Elimination of the first components of X.
!
  if ( l /= 0 ) then
!
!  Scaling of the rows of (C,D)
!
    do i = k+1, k+l
      r = sqrt ( sum ( eca(i,1:n)**2 ) )
      if ( r /= 0.0E+00 ) then
        fdb(i) = fdb(i) / r
        eca(i,1:n) = eca(i,1:n) / r
      end if
    end do

    res_norm = maxval ( maxval ( abs ( eca(k+1:k+l,1:n) ), dim = 1 ) )
    tau = 10.0E+00 * sqrt ( eps ) * res_norm
!
!  Special case for M = 0.
!
    if ( m == 0 ) then
!
!  Just trying this.
!
      mdb = max ( l, n )
!     mdb = l

      call hfti ( eca(k+1,1), l, n, fdb(k+1), mdb, 1, tau, rank1, &
        res_norm, ip )

      if ( n <= rank1 .and. tau <= res_norm ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CON_L2 - Fatal error!'
        write ( *, '(a)' ) '  No feasible point for C * X = D.'
        iflag = 4
        return
      end if

      x(1:n) = fdb(k+1:k+n)
      return

    end if
!
!  Triangular decomposition of the L by N system C * X = D.
!
!  This RES_NORM = -1.0 is a special input flag to Spaeth's modified HFTI.
!
    res_norm = -1.0E+00
!
!  Just trying this.
!
    mdb = max ( l, n )
!   mdb = l

    call hfti ( eca(k+1,1), l, n, fdb(k+1), mdb, 1, tau, rank1, &
      res_norm, ip )

    res_norm = 0.0E+00

  end if

  mm2 = min ( l, rank1 )

  if ( n <= mm2 ) then
    go to 20
  end if

  if ( 0 < mm2 ) then
!
!  Column permutations of E and A reflect the reordering of variables in IP.
!
    do i = 1, k
      do j = 1, mm2
        call r4_swap ( eca(i,j), eca(i,ip(j)) )
      end do
    end do

    do i = k+l+1, k+l+m
      do j = 1, mm2
        call r4_swap ( eca(i,j), eca(i,ip(j)) )
      end do
    end do
!
!  Reduction of the number of problem variables, corresponding to C.
!
    do i = 1, m + l + k

      if ( k + 1 <= i .and. i <= k + l ) then
        cycle
      end if
!
!  Determination of A1(tilde).
!
      eca(i,1) = eca(i,1) / eca(k+1,1)

      do j = 2, mm2
        r = dot_product ( eca(i,1:j-1), eca(k+1:k+j-1,j) )
        eca(i,j) = ( eca(i,j) - r ) / eca(k+j,j)
      end do
!
!  Determination of A2(tilde).
!
      do j = mm2 + 1, n
        eca(i,j) = eca(i,j) - dot_product ( eca(i,1:mm2), eca(k+1:k+mm2,j) )
      end do
!
!  Determination of B(tilde).
!
      fdb(i) = fdb(i) - dot_product ( eca(i,1:mm2), fdb(k+1:k+mm2) )

    end do

  end if
!
!  Singular Value Decomposition.
!
  if ( n <= mm2 ) then
    go to 10
  end if
!
!  Decomposition of (E,F) (tilde) in U * S * V'.
!
  call svdrs ( eca(k+l+1,mm2+1), k+l+m, m, n-mm2, fdb(k+l+1), m, 1, &
    w(mm2+1) )
!
!  V over A, U'*F(tilde) over B, S over W
!
!  Determination of RANK2.
!
  rank2 = 1

  if ( mm2 + 1 < n ) then

    r = 10.0E+00 * sqrt ( eps ) * w(mm2+1)

    do i = mm2+1, n
      if ( w(i) <= r ) then
        exit
      end if
      rank2 = i - mm2
    end do

  end if
!
!  Special case: no inequalities, K = 0.
!
  if ( k == 0 ) then

    do j = mm2 + 1, mm2 + rank2
      w(j+n) = fdb(k+j+l-mm2) / w(j)
    end do

    do i = k + l + 1, k + l + n - mm2
      r = 0.0E+00
      do j = mm2+1, mm2 + rank2
        r = r + eca(i,j) * w(j+n)
      end do
      x(i-l+mm2-k) = r
    end do

    res_norm = sqrt ( sum ( fdb(k+l+rank2+1:m+l+k)**2 ) )
!
!  Solution of the reduced problem with inequalities.
!
!  Transformation to a least distance problem.
!
  else

    do i = 1, k

      do j = mm2+1, mm2 + rank2
        r = 0.0E+00
        do j2 = mm2+1, n
          r = r + eca(i,j2) * eca(k+j2+l-mm2,j)
        end do
        w(n+j) = r / w(j)
      end do

      do j = mm2+1, mm2 + rank2
        eca(i,j) = w(n+j)
      end do

      r = 0.0E+00
      do j = k + l + 1, k + l + rank2
        r = r + eca(i,j-l+mm2-k) * fdb(j)
      end do

      fdb(i) = fdb(i) - r

    end do
!
!  Solution of the LDP problem.
!
    call ldp_l2 ( eca(1,mm2+1), k, rank2, fdb(1), x(mm2+1), res_norm, iflag2 )

    if ( iflag2 /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CON_L2 - Fatal error!'
      write ( *, '(a,i6)' ) '  LDP_L2 returned nonzero error flag IFLAG2 = ', &
        iflag2
      iflag = 5
      return
    end if
!
!  Back substitution.
!
    do i = mm2+1, n
      r = 0.0E+00
      do j = mm2+1, mm2 + rank2
        r = r + eca(k+i+l-mm2,j) * ( x(j) + fdb(k+j+l-mm2) ) / w(j)
      end do
      w(i+n) = r
    end do

    do i = mm2+1, n
      x(i) = w(i+n)
    end do

    r = 0.0E+00
    do i = k + l + rank2 + 1, k + l + m
      r = r + fdb(i)**2
    end do

    res_norm = sqrt ( res_norm**2 + r )

  end if
!
!  Computation of the first part of the components of X.
!
10 continue

  if ( mm2 < 1 ) then
    return
  end if

  if ( mm2+1 <= n ) then

    do i = k+1, k+mm2
      fdb(i) = fdb(i) - dot_product ( eca(i,mm2+1:n), x(mm2+1:n) )
    end do

  end if
!
!  Solution of the system in echelon form.
!
20 continue

  x(mm2) = fdb(k+mm2) / eca(k+mm2,mm2)

  do i = mm2-1, 1, -1
    r = dot_product ( eca(i+k,i+1:mm2), x(i+1:mm2) )
    x(i) = ( fdb(i+k) - r ) / eca(i+k,i)
  end do
!
!  Backwards permutation of the components of X.
!
  do i = mm2, 1, -1
    call r4_swap ( x(i), x(ip(i)) )
  end do

  return
end
subroutine con_li ( m, l, k, n, a, b, c, d, e, f, g, eps, x, res_norm, &
  iflag )

!*****************************************************************************80
!
!! CON_LI minimizes the L-infinity norm of A * X - B subject to linear constraints.
!
!  Discussion:
!
!    The constraints have the form:
!
!      C * X = D
!      F <= E * X <= G
!
!  Modified:
!
!    28 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 239-245,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) L, the number of rows of C.
!
!    Input, integer ( kind = 4 ) K, the number of rows of E.
!
!    Input, integer ( kind = 4 ) N, the number of unknowns, and the number
!    of columns in A, C and E.
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) C(L,N), the equality constraint matrix.
!
!    Input, real ( kind = 4 ) D(L), the equality right hand side.
!
!    Input, real ( kind = 4 ) E(K,N), the inequality constraint matrix.
!
!    Input, real ( kind = 4 ) F(K), the lower bounds for the linear inequalities.
!
!    Input, real ( kind = 4 ) G(K), the upper bounds for the linear inequalities.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for the accuracy test.
!
!    Output, real ( kind = 4 ) X(N), the solution.
!
!    Output, real ( kind = 4 ) RES_NORM, the value of the objective function.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    1, the computation failed.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) c(l,n)
  real ( kind = 4 ) d(l)
  real ( kind = 4 ) e(k,n)
  real ( kind = 4 ) eps
  real ( kind = 4 ) f(k)
  real ( kind = 4 ) g(k)
  real ( kind = 4 ) gt(k)
  real ( kind = 4 ) h(m+l+k+2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) id
  integer ( kind = 4 ) iff
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iipt
  integer ( kind = 4 ) il
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) ipt
  integer ( kind = 4 ) is
  integer ( kind = 4 ) iss
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jf
  integer ( kind = 4 ) js
  integer ( kind = 4 ) jss
  integer ( kind = 4 ) klm11
  real ( kind = 4 ) pivot
  real ( kind = 4 ) q(n+3,m+l+k+2)
  real ( kind = 4 ) q_save
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) smc
  real ( kind = 4 ) sminus
  real ( kind = 4 ) splus
  integer ( kind = 4 ) stage
  real ( kind = 4 ) tmc
  real ( kind = 4 ) val
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) xt(n+3)
  real ( kind = 4 ) z
  real ( kind = 4 ) zmax
  real ( kind = 4 ) zmin
  real ( kind = 4 ) zz
!
!  Initialization.
!
  iflag = 0
  it = 0
  klm11 = m + l + k
  jf = l
  is = 1
  js = 1

  q(1:n,1:l) = transpose ( c(1:l,1:n) )
  q(1:n,l+1:l+m) = transpose ( a(1:m,1:n) )
  q(1:n,l+m+1:l+m+k) = transpose ( e(1:k,1:n) )

  do i = 1, n
    q(i,m+l+k+1) = i
  end do

  q(n+1,1:l) = 0.0E+00
  q(n+1,l+1:l+m) = 1.0E+00
  q(n+1,l+m+1:l+m+k) = 0.0E+00

  q(n+1,k+l+m+1) = 0.0E+00
!
!  The structure of the Q matrix:
!
!            L      M       K     1   1
!
!  Q = N  (  C'  |  A'   |  E'  | I | - )
!      1  (  0   |  1    |  0   | 0 | - )
!      1  ( -D   | -B    | -F   | - | - )
!      1  (  N+J |  N+J  |  N+J | - | - )
!
  q(n+2,1:l) = -d(1:l)
  q(n+2,l+1:l+m) = -b(1:m)
  q(n+2,l+m+1:l+m+k) = -f(1:k)

  do j = 1, k + l + m
    q(n+3,j) = real ( n + j, kind = 4 )
  end do

  gt(1:k) = g(1:k) - f(1:k)

  h(1:l) = d(1:l)
  h(l+1:l+m) = b(1:m)
  h(l+m+1:l+m+k) = f(1:k)
!
!  Stage 1.
!
  if ( l == 0 ) then
    go to 20
  end if

  stage = 1
  il = 0

10 continue

  do

    il = il + 1

    if ( il == l + 1 ) then
      exit
    end if
!
!  Determine the vector to enter the basis.
!
    zmax = -1.0E+00

    do j = js, jf
      z = abs ( q(n+2,j) )
      if ( zmax < z ) then
        in = j
        zmax = z
      end if
    end do
!
!  Determine the vector to leave the basis.
!
    z = -1.0E+00

    if ( is /= n + 1 ) then

      do i = is, n
        zz = abs ( q(i,in) )
        if ( z < zz ) then
          iout = i
          z = zz
        end if
      end do

      if ( eps < z ) then
        go to 140
      end if

    end if

    if ( eps < zmax ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CON_LI - Fatal error!'
      write ( *, '(a)' ) '  The computation failed.'
      iflag = 1
      return
    end if
!
!  Linear dependence in stage 1.
!
    if ( in /= js ) then

      do i = is, n + 3
        call r4_swap ( q(i,in), q(i,js) )
      end do

    end if

    js = js + 1

  end do
!
!  Stage 2
!
20 continue

  stage = 2
  jss = l + 1
  iff = n

  do j = l+1, m + l + k
    h(j) = 1.0E+00
  end do

  i2 = 0

30 continue

  i2 = i2 + 1

  if ( iff < is ) then
    go to 60
  end if

  if ( i2 == k + l + 1 ) then
    go to 50
  end if
!
!  Determine the vector to enter the basis.
!
  ipt = -1

  do j = jss, m + l + k

    if ( h(j) /= 0.0E+00 ) then

      id = abs ( q(n+3,j) )

      if ( id <= n + l + m ) then
        smc = -abs ( q(n+2,j) )
        tmc = 1.0E+00
      else
        smc = q(n+2,j)
        id = id - ( n + l + m )
        tmc = gt(id) - smc
      end if

      if ( ipt < 0 ) then
        zmin = smc
        in = j
        ipt = 0
      else if ( smc < zmin ) then
        zmin = smc
        in = j
        ipt = 0
      end if

      if ( tmc < zmin ) then
        zmin = tmc
        in = j
        ipt = 1
      end if

    end if

  end do

  h(in) = 0.0E+00
!
!  Interchange S and T double primes.
!
  if ( ipt /= 0 ) then

    do i = is, n+3
      q(i,in) = -q(i,in)
    end do

    q(n+2,in) = zmin

  end if
!
!  Determine the vector to leave the basis.
!
  zmax = -1.0E+00

  do i = is, iff

    z = abs ( q(i,in) )

    if ( zmax < z ) then
      zmax = z
      iout = i
    end if

  end do

  if ( eps < zmax ) then
    go to 140
  end if
!
!  No pivot in stage 2.
!
  go to 30
!
!  Interchange columns IN and JSS in stage 2.
!
40 continue

  if ( in /= jss ) then

    do i = is, n + 3
      call r4_swap ( q(i,in), q(i,jss) )
    end do

    call r4_swap ( h(in), h(jss) )

  end if

  jss = jss + 1
!
!  Interchange rows IOUT and IFF in stage 2.
!
  if ( iout /= iff ) then

    do j = js, m + l + k + 1
      call r4_swap ( q(iout,j), q(iff,j) )
    end do

  end if

  iff = iff - 1
  go to 30
!
!  Stage 3
!
50 continue

  iss = iff + 1
  go to 70

60 continue

  iss = is

70 continue

  stage = 3
!
!  Check for a square system.
!
  if ( jss == m + l + k + 1 ) then
    go to 150
  end if
!
!  Find most negative modified marginal cost.
!
  zmin = 1.0E+00

  do j = jss, m + l + k

    splus = q(n+2,j)
    sminus = -splus

    do i = iss, n

      id = abs ( q(i,m+l+k+1) )

      if ( n + l + m < id ) then

        id = id - ( n + l + m )

        if ( eps <= q(i,j) ) then
          splus = splus + gt(id) * q(i,j)
        else if ( q(i,j) <= -eps ) then
          sminus = sminus - gt(id) * q(i,j)
        end if

      end if

    end do

    id = abs ( q(n+3,j) )
    id = id - ( n + l + m )

    if ( 0 < id ) then
      sminus = sminus + gt(id)
    end if

    if ( splus < zmin ) then
      zmin = splus
      in = j
      ipt = 0
    end if

    if ( sminus < zmin ) then
      zmin = sminus
      in = j
      ipt = 1
    end if

  end do
!
!  Check modified marginal cost greater zero.
!
  stage = 5

  if ( -eps < zmin ) then
    go to 90
  end if

  stage = 3
!
!  Interchange entering column with its image.
!
  if ( ipt /= 0 ) then

    id = abs ( q(n+3,in) ) - ( n + l + m )

    do i = is, n+3
      q(i,in) = -q(i,in)
    end do

    if ( 0 < id ) then
      q(n+2,in) = q(n+2,in) + gt(id)
    else
      q(n+1,in) = q(n+1,in) + 2.0E+00
    end if

  end if
!
!  Arrange for all entries in pivot column (except pivot) to be negative.
!
  do i = iss, n

    if ( eps <= q(i,in) ) then

      do j = js, m + l + k

        id = abs ( q(i,m+l+k+1) )

        if ( id < n + l + m + 1 ) then
          q(n+1,j) = q(n+1,j) + 2.0E+00 * q(i,j)
        else
          id = id - ( n + l + m )
          q(n+2,j) = q(n+2,j) + gt(id) * q(i,j)
        end if

        q(i,j) = -q(i,j)

      end do

      q(i,m+l+k+1) = -q(i,m+l+k+1)

    end if

  end do

  iout = n + 1
!
!  Check pivot greater than zero.
!
  if ( q(n+1,in) < eps ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_LI - Fatal error!'
    write ( *, '(a)' ) '  The pivot value is too small.'
    iflag = 1
    return
  end if

  go to 140

80 continue
!
!  Interchange columns in stage 3.
!
  if ( in /= m + l + k ) then

    do i = is, n + 3
      call r4_swap ( q(i,in), q(i,m+l+k) )
    end do

  end if

  klm11 = m + l + k - 1
!
!  Stage 4
!
  stage = 4
!
!  Check for N+1 by N system.
!
  if ( jss == m + l + k ) then
    go to 150
  end if
!
!  Determine the vector to enter the basis.
!
90 continue

  zmin = 1.0E+00

  if ( stage == 4 ) then
    val = 2.0E+00 * q(n+2,m+l+k)
  else
    val = 0.0E+00
  end if

  do j = jss, klm11

    splus = q(n+2,j)
    sminus = val - splus
    id = abs ( q(n+3,j) )
    id = id - ( n + l + m )

    if ( 0 < id ) then
      sminus = gt(id) - splus
    end if

    iipt = 0

    if ( sminus < splus ) then
      splus = sminus
      iipt = 1
    end if

    if ( splus < zmin ) then
      zmin = splus
      in = j
      ipt = iipt
    end if

  end do

  if ( -eps <= zmin ) then
    go to 150
  end if

  if ( ipt == 0 ) then
    go to 120
  end if
!
!  Interchange entering column with its image.
!
  id = abs ( q(n+3,in) )
  id = id - ( n + l + m )

  if ( id <= 0 ) then
    go to 100
  end if

  do i = iss, n+3
    q(i,in) = -q(i,in)
  end do

  q(n+2,in) = zmin
  go to 120

100 continue

  if ( stage == 5 ) then
    go to 110
  end if

  do i = iss, n + 1
    q(i,in) = 2.0E+00 * q(i,m+l+k) - q(i,in)
  end do

  q(n+2,in) = zmin
  q(n+3,in) = -q(n+3,in)
  go to 120

110 continue

  do i = iss, n+3
    q(i,in) = -q(i,in)
  end do

  q(n+1,in) = 2.0E+00 + q(n+1,in)
!
!  Determine the vector to enter the basis.
!
120 continue

  iout = 0

  if ( stage == 5 ) then
    go to 130
  end if

  do i = iss, n+1
    if ( eps < q(i,in) ) then
      z = q(i,m+l+k) / q(i,in)
      if ( iout == 0 .or. z < zmin ) then
        iout = i
        zmin = z
      end if
    end if
  end do

  if ( iout == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_LI - Fatal error!'
    write ( *, '(a)' ) '  The computation failed.'
    iflag = 1
    return
  end if

  go to 140

130 continue

  zmax = -1.0E+00

  do i = iss, n
    z = q(i,in)
    if ( eps < z .and. zmax < z ) then
      iout = i
      zmax = z
    end if
  end do

  if ( iout /= 0 ) then
    go to 140
  end if

  if ( q(n+1,in) <= eps ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CON_LI - Fatal error!'
    write ( *, '(a)' ) '  Pivot value is too small.'
    iflag = 1
    return
  end if

  iout = n + 1
  stage = 3
!
!  Gauss-Jordan elimination.
!
140 continue

  pivot = q(iout,in)

  do j = js, m + l + k
    q(iout,j) = q(iout,j) / pivot
  end do

  do j = js, m + l + k
    if ( j /= in ) then
      z = -q(iout,j)

      if ( is <= iout .and. iout <= n+2 ) then
        q_save = q(iout,j)
      end if

      q(is:n+2,j) = q(is:n+2,j) + z * q(is:n+2,in)

      if ( is <= iout .and. iout <= n+2 ) then
        q(iout,j) = q_save
      end if

    end if
  end do

  do i = is, n+2
    q(i,in) = -q(i,in) / pivot
  end do

  q(iout,in) = 1.0E+00 / pivot

  call r4_swap ( q(iout,m+l+k+1), q(n+3,in) )

  it = it + 1

  if ( stage == 1 ) then
!
!  Interchange columns IN and JF.
!
    if ( in /= jf ) then

      do i = is, n + 3
        call r4_swap ( q(i,in), q(i,jf) )
      end do

    end if

    jf = jf - 1
!
!  Interchange rows IOUT and IS.
!
    if ( iout /= is ) then

      do j = js, m + l + k + 1
        call r4_swap ( q(iout,j), q(is,j) )
      end do

    end if

    is = is + 1

    go to 10

  else if ( stage == 2 ) then

    go to 40

  else if ( stage == 3 ) then

    go to 80

  else if ( stage == 4 ) then

    go to 90

  else if ( stage == 5 ) then

    go to 90

  end if
!
!  Prepare output
!
150 continue

  if ( q(n+1,m+l+k+1) == 0.0E+00 ) then
    res_norm = 0.0E+00
  else
    res_norm = q(n+2,m+l+k)
  end if

  xt(1:n) = 0.0E+00

  do j = js, jss - 1
    i = q(n+3,j)
    xt(i) = q(n+2,j)
  end do

  h(1:l) = 0.0E+00

  do i = iss, n + 1

    z = res_norm
    id = q(i,m+l+k+1)

    if ( id /= 0 ) then

      j = abs ( id ) - ( n + l + m )

      if ( 0 < j ) then
        if ( id < 0 ) then
          z = gt(j)
        else
          z = 0.0E+00
        end if
      end if

      j = j + l + m
      h(j) = z * sign ( 1.0E+00, q(i,m+l+k+1) )

    end if

  end do

  do j = jss, m + l + k

    z = res_norm
    id = q(n+3,j)

    if ( id /= 0 ) then

      i = abs ( id ) - ( n + l + m )

      if ( 0 < i ) then
        if ( id < 0 ) then
          z = gt(i)
        else
          z = 0.0E+00
        end if
      end if

      i = i + l + m
      h(i) = ( z - q(n+2,j) ) * sign ( 1.0E+00, q(n+3,j) )

    end if

  end do

  x(1:n) = xt(1:n)

  return
end
subroutine cwlr_l1 ( a, m, n, b, s, eps, z, ml, mu, x, e, d, iflag )

!*****************************************************************************80
!
!! CWLR_L1 minimizes the L1 norm of A*x-b using clustering techniques.
!
!  Discussion:
!
!    The routines UPROWB and DOROWB are not described in the text.
!    Hence, this routine cannot be used as it stands.  Moreover,
!    apparently some changes need to be made to AFK_L1 for it to be
!    used by this routine.  All in all, a sorry show.
!
!  Modified:
!
!    28 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 163-165,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the M by N system matrix.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input, integer ( kind = 4 ) S, the number of clusters.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for an accuracy test in A478_L1.
!
!    Input/output, integer ( kind = 4 ) Z(M).  On input, Z contains an initial
!    partition of the data, that is, an assignment of each data item
!    to one of the clusters between 1 and S.  On output, Z contains
!    an improved partition.
!
!    Input, integer ( kind = 4 ) ML, MU, lower and upper bounds for the number
!    of items per cluster.  ML must be at least 1.
!
!    Output, real ( kind = 4 ) X(S,N), contains the computed parameter vectors.
!
!    Output, real ( kind = 4 ) E(S), the sum of squared residuals.
!
!    Output, real ( kind = 4 ) D, the value of the objective function.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error was detected.
!    1, ML <= N.
!    2, MU <= ML.
!    3, some cluster assignment Z(I) is less than ML.
!    4, some cluster assignment Z(I) is greater than MU.
!    5, an error occurred in A478_LI.
!    6, no convergence after 15 iterations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) aa(m+2,n+2)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) bb(m)
  real ( kind = 4 ) beta(n)
  real ( kind = 4 ) c(m+2,n+2,s)
  real ( kind = 4 ) cb(m,s)
  real ( kind = 4 ) cbdo(m)
  real ( kind = 4 ) cbup(m)
  real ( kind = 4 ) cdo(m+2,n+2)
  real ( kind = 4 ) cj
  real ( kind = 4 ) cp
  real ( kind = 4 ) cq
  real ( kind = 4 ) cup(m+2,n+2)
  real ( kind = 4 ) d
  real ( kind = 4 ) e(s)
  real ( kind = 4 ) ep
  real ( kind = 4 ) eps
  real ( kind = 4 ) eq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) index(m)
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jm
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mj(s)
  integer ( kind = 4 ) mjj
  integer ( kind = 4 ) mjp
  integer ( kind = 4 ) mjq
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res(m)
  real ( kind = 4 ) rhs
  real ( kind = 4 ) rowi(n)
  real ( kind = 4 ), parameter :: rr = 0.999E+00
  real ( kind = 4 ) sad
  integer ( kind = 4 ) ss(m)
  logical stage
  real ( kind = 4 ) x(s,n)
  real ( kind = 4 ) xp(n)
  real ( kind = 4 ) xq(n)
  integer ( kind = 4 ) z(m)

  iflag = 0
  it = 0

  if ( ml <= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CWLR_L1 - Fatal error!'
    write ( *, '(a)' ) '  ML <= N.'
    iflag = 1
    return
  end if

  if ( mu <= ml ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CWLR_L1 - Fatal error!'
    write ( *, '(a)' ) '  MU <= ML.'
    iflag = 2
    return
  end if

  if ( .false. ) then
  if ( m-(s-i)*(n+1) < mu ) then
    iflag = 2
    return
  end if
  end if
!
!  Initialization.
!
  d = 0.0E+00

  do j = 1, s

    jm = 0

    do i = 1, m

      if ( z(i) == j ) then

        jm = jm + 1
        if ( mu < jm ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CWLR_L1 - Fatal error!'
          write ( *, '(a)' ) '  A cluster contains too many items.'
          iflag = 3
          return
        end if

        bb(jm) = b(i)
        aa(jm,1:n) = a(i,1:n)

      end if

    end do

    if ( jm < ml ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_L1 - Fatal error!'
      write ( *, '(a)' ) '  A cluster contains too few items.'
      iflag = 4
      return
    end if

    mj(j) = jm
    stage = .true.
!
!  THIS CALL MUST BE REVISED...ESPECIALLY SINCE I NO LONGER
!  PASS BACK INFORMATION IN A.
!
    call a478_l1 ( aa, m, n, bb, eps, rank, beta, res, iflag2 )

    if ( iflag2 /= 0 ) then
      iflag = 5
      return
    end if

    sad = aa(jm+1,n+1)
    aa(jm+2,n+2) = 0.0E+00
    e(j) = sad

    x(j,1:n) = beta(1:n)

    d = d + sad

    if ( s <= 1 ) then
      return
    end if

    cb(1:jm,j) = bb(1:jm)

    c(1:jm+2,1:n+2,j) = aa(1:jm+2,1:n+2)

  end do
!
!  Exchange steps
!
  is = 0
  i = 0

  do j = 1, s
    ic = 0
    do k = 1, m
      if ( z(k) == j ) then
        ic = ic + 1
        index(k) = ic
      end if
    end do
  end do

  is = is + 1

  if ( m < is ) then
    return
  end if

  do

    i = i + 1

    if ( m < i ) then

      if ( itmax <= it ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CWLR_L1 - Fatal error!'
        write ( *, '(a)' ) '  Number of iterations exceeded.'
        iflag = 6
        exit
      end if

      it = it + 1

      i = 1

    end if

    p = z(i)

    if ( mj(p) <= ml ) then

      is = is + 1

      if ( m < is ) then
        exit
      end if

      cycle

    end if

    cq = huge ( cq )

    do j = 1, s

      mjj = mj(j)

      bb(1:mjj+2) = cb(1:mjj+2,j)
      aa(1:mjj+2,1:n+2) = c(1:mjj+2,1:n+2,j)

      if ( j == p ) then

        jm = mjj - 1
!
!  Source for DOROWB was not available...
!
!       call dorowb ( jm, n, m+2, n+2, aa, index(i), eps )

      else

        jm = mjj + 1

        if ( mu < jm ) then
          cycle
        end if

        rowi(1:n) = a(i,1:n)

        rhs = b(i)
!
!  Source for UPROWB was not available...
!
!       call uprowb ( jm, n, m+2, n+2, aa, rowi, rhs, eps, bb )

      end if

      iflag = 0

      call a478_l1 ( aa, m, n, bb, eps, rank, beta, res, iflag2 )

      if ( iflag2 /= 0 ) then
        iflag = 5
        return
      end if

      sad = aa(jm+1,n+1)
      aa(jm+2,n+2) = 0.0E+00

      if ( j == p ) then

        cp = sad
        xp(1:n) = beta(1:n)
        cbdo(1:jm) = bb(1:jm)
        cdo(1:jm+2,1:n+2) = aa(1:jm+2,1:n+2)

      else

        cj = sad - e(j)

        if ( cj < cq ) then

          cq = cj
          eq = sad
          xq(1:n) = beta(1:n)
          cbup(1:jm) = bb(1:jm)
          cup(1:jm+2,1:n+2) = aa(1:jm+2,1:n+2)
          q = j

        end if

      end if

    end do

    ep = e(p) - cp

    if ( ep * rr <= cq ) then

      is = is + 1

      if ( m < is ) then
        return
      end if
!
!  Updates for a successful exchange.
!
    else

      do k = 1, m
        if ( z(k) == p ) then

          if ( index(k) == mj(p) ) then
            index(k) = index(i)
          end if

          if ( index(k) == mj(p) ) then
            exit
          end if

        end if
      end do

      is = 0
      e(p) = cp
      e(q) = eq
      d = d - ep + cq
      mj(p) = mj(p) - 1
      mj(q) = mj(q) + 1
      index(i) = mj(q)
      z(i) = q

      x(p,1:n) = xp(1:n)
      x(q,1:n) = xq(1:n)

      mjp = mj(p)
      cb(1:mjp,p) = cbdo(1:mjp)
      c(1:mjp+2,1:n+2,p) = cdo(1:mjp+2,1:n+2)

      mjq = mj(q)
      cb(1:mjq,q) = cbup(1:mjq)
      c(1:mjq+2,1:n+2,q) = cup(1:mjq+2,1:n+2)

    end if

  end do

  return
end
subroutine cwlr_l2 ( a, m, n, b, s, eps, z, ml, mu, x, e, d, iflag )

!*****************************************************************************80
!
!! CWLR_L2 minimizes the L2 norm of A*x-b using clustering techniques.
!
!  Discussion:
!
!    On input, the user must have set up an initial partition of the data.
!
!  Modified:
!
!    30 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 149-151,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the M by N system matrix.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input, integer ( kind = 4 ) S, the number of clusters.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for an accuracy test in INEXCL.
!
!    Input/output, integer ( kind = 4 ) Z(M).  On input, Z contains an initial
!    partition of the data, that is, an assignment of each data item
!    to one of the clusters between 1 and NC.  On output, Z contains
!    an improved partition.
!
!    Input/output, integer ( kind = 4 ) ML, MU, the minimum and maximum number of
!    observations to be allowed per cluster.  ML must be greater than N.
!
!    Output, real ( kind = 4 ) X(S,N), contains the computed parameter vectors.
!
!    Output, real ( kind = 4 ) E(S), the sum of squared residuals.
!
!    Output, real ( kind = 4 ) D, the value of the objective function.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error.
!    1, ML is not greater than N.
!    2, MU <= ML.
!    3, unfathomable error, currently disabled.
!    4, a cluster assignment was less than 1 or greater than S.
!    5, population of some cluster is greater than MU.
!    6, population of some cluster is less than ML.
!    7, error in INEXCL.
!    8, iteration limit exceeded.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) ai(n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) bi
  real ( kind = 4 ) cj
  real ( kind = 4 ) cp
  real ( kind = 4 ) cq
  real ( kind = 4 ) d
  real ( kind = 4 ) e(s)
  real ( kind = 4 ) ej
  real ( kind = 4 ) ep
  real ( kind = 4 ) eq
  real ( kind = 4 ) eps
  real ( kind = 4 ) f(n,s)
  real ( kind = 4 ) fj(n)
  real ( kind = 4 ) fp(n)
  real ( kind = 4 ) fq(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mj(s)
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 4 ) r(s,((n-1)*n)/2)
  real ( kind = 4 ) rj(((n-1)*n)/2)
  real ( kind = 4 ) rp(((n-1)*n)/2)
  real ( kind = 4 ) rq(((n-1)*n)/2)
  real ( kind = 4 ), parameter :: rr = 0.995E+00
  integer ( kind = 4 ) s2
  real ( kind = 4 ) t(n,s)
  real ( kind = 4 ) tj(n)
  real ( kind = 4 ) tp(n)
  real ( kind = 4 ) tq(n)
  real ( kind = 4 ) wi
  real ( kind = 4 ) x(s,n)
  integer ( kind = 4 ) z(m)

  it = 0
  d = 0.0E+00
  iflag = 0
  s2 = ( ( n - 1 ) * n ) / 2

  if ( ml <= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
    write ( *, '(a)' ) '  Input value of ML is not legal.'
    write ( *, '(a)' ) '  We have ML <= N.'
    iflag = 1
    return
  end if

  if ( mu <= ml ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
    write ( *, '(a)' ) '  Input values of ML and MU are not legal.'
    write ( *, '(a)' ) '  We have MU <= ML.'
    iflag = 2
    return
  end if
!
!  WHAT IS THIS RESTRICTION HERE FOR?
!
  if ( .false. ) then

    if ( m - (s-1) * (n+1) < mu ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
      write ( *, '(a)' ) '  Input value of MU is not legal.'
      iflag = 3
      return
    end if

  end if

  e(1:s) = 0.0E+00
  f(1:n,1:s) = 0.0E+00
  t(1:n,1:s) = 0.0E+00
  r(1:s,1:s2) = 0.0E+00
!
!  Make sure the cluster assignments are legal.
!
  do i = 1, m

    j = z(i)

    if ( j < 1 .or. s < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
      write ( *, '(a)' ) '  A cluster assignment was illegal.'
      write ( *, '(a)' ) '  We require 1 <= Z(I) <= S.'
      write ( *, '(a,i6)' ) '  We have S = ', s
      write ( *, '(a,i6)' ) '  but for index I = ', i
      write ( *, '(a,i6)' ) '  Z(I) = ', z(i)
      iflag = 4
      return
    end if

  end do
!
!  Determine the cluster populations.
!
  mj(1:s) = 0
  do i = 1, m
    j = z(i)
    mj(j) = mj(j) + 1
  end do

  do i = 1, m

    j = z(i)
    wi = 1.0E+00
    bi = b(i)
    ai(1:n) = a(i,1:n)

    call inexcl ( n, ai, bi, wi, j, s, s2, eps, f, t, r, fj, tj, &
      rj, e(j), iflag2 )

    if ( iflag2 /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
      write ( *, '(a,i6)' ) '  INEXCL returned error flag IFLAG2 = ', iflag2
      iflag = 7
      return
    end if

    f(1:n,j) = fj(1:n)
    t(1:n,j) = tj(1:n)
    r(j,1:s2) = rj(1:s2)

  end do

  d = sum ( e(1:s) )

  do j = 1, s

    if ( mu < mj(j) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
      write ( *, '(a)' ) '  A cluster contains too many items.'
      iflag = 5
      return
    end if

    if ( mj(j) < ml ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
      write ( *, '(a)' ) '  A cluster contains too few items.'
      iflag = 6
      return
    end if

  end do

  if ( 1 < s ) then

    i = 0
    is = 1

    do

      i = i + 1

      if ( m < i ) then

        if ( itmax <= it ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
          write ( *, '(a)' ) '  Number of iterations exceeded.'
          iflag = 8
          exit
        end if

        it = it + 1
        i = 1

      end if

      p = z(i)

      if ( mj(p) <= ml ) then
        is = is + 1
        if ( m < is ) then
          exit
        end if
        cycle
      end if

      eq = huge ( eq )

      do j = i, s

        if ( j /= p .and. mu <= mj(j) ) then
          cycle
        end if

        ej = e(j)
        bi = b(i)
        ai(1:n) = a(i,1:n)

        if ( j == p ) then

          ep = ej
          wi = -1.0E+00

          call inexcl ( n, ai, bi, wi, p, s, s2, eps, f, t, r, fp, &
            tp, rp, ep, iflag2 )

          if ( iflag2 /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
            write ( *, '(a,i6)' ) '  INEXCL returned error flag IFLAG2 = ', &
              iflag2
            iflag = 7
            return
          end if

          cp = ep

        else

          wi = 1.0E+00

          call inexcl ( n, ai, bi, wi, j, s, s2, eps, f, t, r, &
            fj, tj, rj, ej, iflag2 )

          if ( iflag2 /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'CWLR_L2 - Fatal error!'
            write ( *, '(a,i6)' ) '  INEXCL returned error flag IFLAG2 = ', &
              iflag2
            iflag = 7
            return
          end if

          cj = ej
          ej = cj - e(j)

          if ( ej < eq ) then

            eq = ej
            cq = cj
            q = i

            fq(1:n) = fj(1:n)
            tq(1:n) = tj(1:n)
            rq(1:s2) = rj(1:s2)

          end if

        end if

      end do

      ep = e(p) - cp

      if ( ep * rr <= eq ) then

        is = is + 1
        if ( m < is ) then
          exit
        end if

      else

        is = 0
        mj(p) = mj(p) - 1
        mj(q) = mj(q) + 1
        d = d - ep + eq
        e(p) = cp
        e(q) = cq

        f(1:n,p) = fp(1:n)
        f(1:n,q) = fq(1:n)
        t(1:n,p) = tp(1:n)
        t(1:n,q) = tq(1:n)

        r(p,1:s2) = rp(1:s2)
        r(q,1:s2) = rq(1:s2)

        z(i) = q

      end if

    end do

  end if

  do j = 1, s
    do k = n, 1, -1
      tj(k) = t(k,j)
      if ( k /= n ) then
        nn = ((k-1)*(n+n-k))/2 + 1
        do l = k+1, n
          tj(k) = tj(k) - r(j,nn) * tj(l)
          nn = nn + 1
        end do
      end if
      x(j,k) = tj(k)
    end do
  end do

  return
end
subroutine cwlr_li ( a, m, n, b, s, eps, z, ml, mu, x, e, d, iflag )

!*****************************************************************************80
!
!! CWLR_LI minimizes the L-infinity norm of A*X-B.
!
!  Modified:
!
!    30 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 176-177,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the M by N system matrix.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input, integer ( kind = 4 ) S, the number of clusters.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for an accuracy test.
!
!    Input/output, integer ( kind = 4 ) Z(M).  On input, Z contains an initial
!    partition of the data, that is, an assignment of each data item
!    to one of the clusters between 1 and NC.  On output, Z contains
!    an improved partition.
!
!    Input/output, integer ( kind = 4 ) ML, MU, the minimum and maximum number of
!    observations to be allowed per cluster.  ML must be greater than N.
!
!    Output, real ( kind = 4 ) X(S,N), contains the computed parameter vectors.
!
!    Output, real ( kind = 4 ) E(S), the sum of squared residuals.
!
!    Output, real ( kind = 4 ) D, the value of the objective function.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    1, ML <= N.
!    2, MU <= ML.
!    3, unfathomable error, currently disabled.
!    4, a cluster assignment was illegal.
!    5, population of some cluster is greater than MU.
!    6, population of some cluster is less than ML.
!    7, error return from A495_LI.
!    8, iteration limit exceeded.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) at(n+3,m+1)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) bb(m+1)
  real ( kind = 4 ) beta(n+3)
  real ( kind = 4 ) cj
  real ( kind = 4 ) cp
  real ( kind = 4 ) cq
  real ( kind = 4 ) d
  real ( kind = 4 ) e(s)
  real ( kind = 4 ) ep
  real ( kind = 4 ) eps
  real ( kind = 4 ) eq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it2
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jm
  integer ( kind = 4 ) mj(s)
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rank
  real ( kind = 4 ) relerr
  real ( kind = 4 ) res_norm
  real ( kind = 4 ), parameter :: rr  =  0.999E+00
  real ( kind = 4 ) x(s,n)
  real ( kind = 4 ) xp(n)
  real ( kind = 4 ) xq(n)
  integer ( kind = 4 ) z(m)

  iflag = 0
  it = 0

  if ( ml <= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CWLR_LI - Fatal error!'
    write ( *, '(a)' ) '  Input value of ML is not legal.'
    write ( *, '(a)' ) '  We have ML <= N.'
    iflag = 1
    return
  end if

  if ( mu <= ml ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CWLR_LI - Fatal error!'
    write ( *, '(a)' ) '  Input values of ML and MU are not legal.'
    write ( *, '(a)' ) '  We have MU <= ML.'
    iflag = 2
    return
  end if

  if ( .false. ) then
  if ( m-(s-1)*(n+1) < mu ) then
    iflag = 3
    return
  end if
  end if
!
!  Make sure the cluster assignments are legal.
!
  do i = 1, m

    j = z(i)

    if ( j < 1 .or. s < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_LI - Fatal error!'
      write ( *, '(a)' ) '  A cluster assignment was illegal.'
      write ( *, '(a)' ) '  We require 1 <= Z(I) <= S.'
      write ( *, '(a,i6)' ) '  We have S = ', s
      write ( *, '(a,i6)' ) '  but for index I = ', i
      write ( *, '(a,i6)' ) '  Z(I) = ', z(i)
      iflag = 4
      return
    end if

  end do
!
!  Initialization.
!
  d = 0.0E+00

  do j = 1, s

    jm = 0

    do i = 1, m

      if ( z(i) == j ) then

        jm = jm + 1

        if ( mu < jm ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CWLR_LI - Fatal error!'
          write ( *, '(a)' ) '  The population of a cluster exceeds MU.'
          write ( *, '(a,i6)' ) '  MU = ', mu
          write ( *, '(a,i6)' ) '  Cluster index J = ', j
          write ( *, '(a,i6)' ) '  Population is at least ', JM
          iflag = 5
          return
        end if

        bb(jm) = b(i)
        at(1:n,jm) = a(i,1:n)

      end if

    end do

    if ( jm < ml ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_LI - Fatal error!'
      write ( *, '(a)' ) '  The population of a cluster less than ML.'
      write ( *, '(a,i6)' ) '  ML = ', ml
      write ( *, '(a,i6)' ) '  Cluster index J = ', j
      write ( *, '(a,i6)' ) '  Population is only ', JM
      iflag = 6
      return
    end if

    mj(j) = jm
    relerr = 0.0E+00

    call a495_li ( at, jm, n, bb, eps, relerr, beta, rank, res_norm, iflag2 )

    if ( iflag2 == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWLR_LI - Fatal error!'
      write ( *, '(a,i6)' ) '  A495_LI returned error flag IFLAG2 = ', iflag2
      iflag = 7
      return
    end if

    e(j) = res_norm
    d = d + res_norm

    x(j,1:n) = beta(1:n)

  end do

  if ( s <= 1 ) then
    return
  end if
!
!  Exchange steps.
!
  i = 0
  is = 1

  if ( m < is ) then
    return
  end if

  do

    i = i + 1

    if ( m < i ) then

      if ( itmax <= it ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CWLR_LI - Fatal error!'
        write ( *, '(a)' ) '  Number of iterations exceeded.'
        iflag = 8
        return
      end if

      it = it + 1
      i = 1

    end if

    p = z(i)

    if ( mj(p) <= ml ) then

      is = is + 1

      if ( m < is ) then
        return
      end if

      cycle

    end if

    cq = huge ( cq )

    do j = 1, s

      jm = 0

      if ( j /= p ) then
        z(i) = j
      else
        z(i) = 0
      end if

      do r = 1, m

        if ( z(r) == j ) then

          jm = jm + 1

          if ( mu < jm ) then
            go to 130
          end if

          bb(jm) = b(r)
          at(1:n,jm) = a(r,1:n)

        end if

      end do

      relerr = 0.0E+00

      call a495_li ( at, jm, n, bb, eps, relerr, beta, rank, res_norm, iflag2 )

      if ( iflag2 == 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CWLR_LI - Fatal error!'
        write ( *, '(a,i6)' ) '  A495_LI returned error flag IFLAG2 = ', iflag2
        iflag = 7
        return
      end if

      if ( j == p ) then

        cp = res_norm

        xp(1:n) = beta(1:n)

      else

        cj = res_norm - e(j)

        if ( cj < cq ) then
          cq = cj
          eq = res_norm
          xq(1:n) = beta(1:n)
          q = j
        end if

      end if

130 continue

    end do

    z(i) = p
    ep = e(p) - cp
!
!  Updates for a successful exchange.
!
    if ( ep * rr <= cq ) then

      is = is + 1

      if ( m < is ) then
        return
      end if

    else

      is = 0
      e(p) = cp
      e(q) = eq
      d = d - ep + cq
      mj(p) = mj(p) - 1
      mj(q) = mj(q) + 1
      z(i) = q

      x(p,1:n) = xp(1:n)
      x(q,1:n) = xq(1:n)

    end if

  end do

  return
end
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
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
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine example_multi_size ( file_name, m, n, s, file_name2 )

!*****************************************************************************80
!
!! EXAMPLE_MULTI_SIZE returns values for a multiple system file.
!
!  Discussion:
!
!    Call this routine to find out the values of M, N, S, and
!    the name of the first file defining a subsystem.
!
!    Then ALLOCATE the sizes of various arrays, and for then
!    call EXAMPLE_READ S times.
!
!  Modified:
!
!    13 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) M, N, S, the number of rows (or observations),
!    the number columns (or variables), and the number of clusters or
!    subsystems.
!
!    Output, character ( len = * ) FILE_NAME2, the name of the first
!    file containing subsystem information.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a_title_read
  integer ( kind = 4 ) b_title_read
  logical done
  character ( len = * ) file_name
  character ( len = * ) file_name2
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  character ( len = 100 ) line
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrec
  integer ( kind = 4 ) s

  m = 0
  n = 0
  s = 0
  file_name2 = ' '

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXAMPLE_MULTI_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the data file "' &
      // trim ( file_name ) // '".'
    return
  end if

  a_title_read = 0
  b_title_read = 0
  nrec = 0

  do

    read ( file_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXAMPLE_SIZE - Fatal error!'
      write ( *, '(a)' ) '  Bad input in data file "' &
      // trim ( file_name ) // '".'
      stop
    end if

    nrec = nrec + 1
!
!  Blank lines are ignored.
!
    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Comment lines begin with a "#"
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  Read the value of N.
!
    if ( n == 0 ) then

      call s_to_i4 ( line, ncol, ierror, last )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Bad data in example file while reading N.'
        stop
      end if

      n = ncol - 2

    else if ( m == 0 ) then

      call s_to_i4 ( line, m, ierror, last )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Bad data in example file while reading M.'
        write ( *, '(a,i6)' ) '  Record number = ', nrec
        write ( *, '(a)' ) '"' // trim ( line ) // '"'
        stop
      end if

    else if ( s == 0 ) then

      call s_to_i4 ( line, s, ierror, last )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Bad data in example file while reading S.'
        write ( *, '(a,i6)' ) '  Record number = ', nrec
        write ( *, '(a)' ) '"' // trim ( line ) // '"'
        stop
      end if

    else

      file_name2 = line

      exit

    end if

  end do

  close ( unit = file_unit )

  return
end
subroutine example_print ( m, n, a, b, a_title, b_title )

!*****************************************************************************80
!
!! EXAMPLE_PRINT prints data from an example file.
!
!  Discussion:
!
!    The data is printed out in groups of 5 columns at a time.
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows (or observations) and
!    columns (or variables).
!
!    Input, real ( kind = 4 ) A(M,N), B(M), the coefficient and right hand side data.
!
!    Input, character ( len = * ) A_TITLE(N), B_TITLE, the names of
!    the variables and the right hand side.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  character ( len = * ) a_title(n)
  character ( len = 10 ) a_title2(n)
  real ( kind = 4 ) b(m)
  character ( len = * ) b_title
  character ( len = 10 ) b_title2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
!  Make 10-character right-justified copies of the titles.
!
  do i = 1, n
    a_title2(i) = a_title(i)
    a_title2(i) = adjustr ( a_title2(i) )
  end do

  b_title2 = b_title
  b_title2 = adjustr ( b_title2 )

  jlo = 1
  jhi = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXAMPLE_PRINT:'
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n
  write ( *, '(a)' ) ' '

  do jlo = 1, n+1, 5

    jhi = jlo + 4
    jhi = min ( jhi, n+1 )

    if ( jhi <= n ) then
      write ( *, '(a3,2x,5a10)' ) &
        '  #', ( a_title2(j), j = jlo, jhi )
    else
      write ( *, '(a3,2x,5a10)' ) &
        '  #', ( a_title2(j), j = jlo, jhi-1 ), b_title2
    end if

    write ( *, '(a)' ) ' '

    do i = 1, m

      if ( jhi <= n ) then
        write ( *, '(i3,2x,5f10.3)' ) i, ( a(i,j), j = jlo, jhi )
      else
        write ( *, '(i3,2x,5f10.3)' ) i, ( a(i,j), j = jlo, jhi-1 ), b(i)
      end if

    end do

  end do

  return
end
subroutine example_read ( file_name, m, n, a, b, a_title, b_title )

!*****************************************************************************80
!
!! EXAMPLE_READ reads data from an example file.
!
!  Modified:
!
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows (or observations) and
!    columns (or variables).
!
!    Output, real ( kind = 4 ) A(M,N), B(M), the coefficient and righthand 
!    side data.
!
!    Output, character ( len = * ) A_TITLE(N), B_TITLE, the names of
!    the variables and the right hand side.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  character ( len = * ) a_title(n)
  integer ( kind = 4 ) a_title_read
  real ( kind = 4 ) b(m)
  character ( len = * ) b_title
  integer ( kind = 4 ) b_title_read
  logical, parameter :: debug = .false.
  logical done
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) last
  character ( len = 100 ) line
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrec
  real ( kind = 4 ) r

  m2 = 0
  n2 = 0

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXAMPLE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the data file.'
    return
  end if

  i = 0
  a_title_read = -1
  b_title_read = 0
  nrec = 0

  do

    read ( file_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    nrec = nrec + 1
!
!  Blank lines are ignored.
!
    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Comment lines begin with a "#"
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  Read the value of N.
!
    if ( n2 == 0 ) then

      call s_to_i4 ( line, ncol, ierror, last )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_READ - Fatal error!'
        write ( *, '(a)' ) '  Bad data in example file while reading N.'
        stop
      end if

      n2 = ncol - 2

      if ( n2 /= n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_READ - Fatal error!'
        write ( *, '(a,i6)' ) '  Input value of N =     ', n
        write ( *, '(a,i6)' ) '  Data file value of N = ', n2
        stop
      end if

    else if ( m2 == 0 ) then

      call s_to_i4 ( line, m2, ierror, last )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_READ - Fatal error!'
        write ( *, '(a)' ) '  Bad data in example file while reading M.'
        stop
      end if

      if ( m2 /= m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_READ - Fatal error!'
        write ( *, '(a,i6)' ) '  Input value of M =     ', m
        write ( *, '(a,i6)' ) '  Data file value of M = ', m2
        stop
      end if
!
!  Read the variable labels, one per line.
!
    else if ( a_title_read < n ) then

      a_title_read = a_title_read + 1
!
!  Skip "Index" title.
!
      if ( a_title_read == 0 ) then

      else
        a_title(a_title_read) = line
      end if

    else if ( b_title_read < 1 ) then

      b_title_read = b_title_read + 1
      b_title = line

    else

      done = .true.

      i = i + 1
!
!  Skip the index record.
!
      call r4_next ( line, r, done )

      do j = 1, n
        call r4_next ( line, r, done )
        a(i,j) = r
      end do

      call r4_next ( line, r, done )
      b(i) = r

    end if

    if ( 0 < m2 .and. m2 <= i ) then
      exit
    end if

  end do

  close ( unit = file_unit )

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXAMPLE_READ - Note:'
    write ( *, '(a,i6)' ) '  Number of records read was ', nrec
  end if

  return
end
subroutine example_size ( file_name, m, n )

!*****************************************************************************80
!
!! EXAMPLE_SIZE returns the values of M and N in an example file.
!
!  Discussion:
!
!    Call this routine to find out the values of M and N, then
!    ALLOCATE the sizes of various arrays, and call EXAMPLE_READ.
!
!  Modified:
!
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) M, N, the number of rows (or observations) 
!    and columns (or variables).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a_title_read
  integer ( kind = 4 ) b_title_read
  logical done
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  character ( len = 100 ) line
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrec

  m = 0
  n = 0

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXAMPLE_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the data file.'
    return
  end if

  a_title_read = 0
  b_title_read = 0
  nrec = 0

  do

    read ( file_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXAMPLE_SIZE - Fatal error!'
      write ( *, '(a)' ) '  Bad data in example file.'
      stop
    end if

    nrec = nrec + 1
!
!  Blank lines are ignored.
!
    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Comment lines begin with a "#"
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  Read the value of N.
!
    if ( n == 0 ) then

      call s_to_i4 ( line, ncol, ierror, last )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Bad data in example file while reading N.'
        stop
      end if

      n = ncol - 2

    else if ( m == 0 ) then

      call s_to_i4 ( line, m, ierror, last )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXAMPLE_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Bad data in example file while reading M.'
        write ( *, '(a,i6)' ) '  Record number = ', nrec
        write ( *, '(a)' ) '"' // trim ( line ) // '"'
        stop
      end if

      exit

    end if

  end do

  close ( unit = file_unit )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Examples:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  return
end
subroutine g1 ( a, b, c, s, sig )

!*****************************************************************************80
!
!! G1 computes an orthogonal rotation matrix.
!
!  Discussion:
!
!    The routine is given a vector
!
!      ( A )
!      ( B )
!
!    and computes quantities C and S so that
!
!      (  C  S )   ( A )   ( SQRT ( A**2 + B**2 ) )
!      ( -S  C ) * ( B ) = ( 0                    )
!
!  Modified:
!
!    30 November 2000
!
!  Reference:
!
!    Charles Lawson and Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974,
!    Revised edition, SIAM, 1995.
!    QA275.L38
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the entries of the vector.
!
!    Output, real ( kind = 4 ) C, S, the entries of the rotation matrix.
!
!    Output, real ( kind = 4 ) SIG, the value of SQRT ( A**2 + B**2 ).
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  real ( kind = 4 ) s
  real ( kind = 4 ) sig
  real ( kind = 4 ) xr
  real ( kind = 4 ) yr

  if ( abs ( b ) < abs ( a ) ) then
    xr = b / a
    yr = sqrt ( 1.0E+00 + xr**2 )
    c = sign ( 1.0E+00 / yr, a )
    s = c * xr
    sig = abs ( a ) * yr
  else if ( b /= 0.0E+00 ) then
    xr = a / b
    yr = sqrt ( 1.0E+00 + xr**2 )
    s = sign ( 1.0E+00 / yr, b )
    c = s * xr
    sig = abs ( b ) * yr
  else
    c = 0.0E+00
    s = 1.0E+00
    sig = 0.0E+00
  end if

  return
end
subroutine gen ( m, n, a, b, iseed )

!*****************************************************************************80
!
!! GEN generates a random matrix A and right hand side B.
!
!  Modified:
!
!    03 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 45,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of observations, equations,
!    or rows.
!
!    Input, integer ( kind = 4 ) N, the number of variables, or columns.
!
!    Output, real ( kind = 4 ) A(M,N), the M by N matrix.
!
!    Output, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input/output, integer ( kind = 4 ) ISEED, a seed value for the random
!    number generator.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) g
  real ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iseed
  integer ( kind = 4 ) j
  real ( kind = 4 ) s
  real ( kind = 4 ) uniform_01_sample
  real ( kind = 4 ) v

  do i = 1, m

    s = 0.0E+00

    do j = 1, n
      h = uniform_01_sample ( iseed )
      s = s + h
      a(i,j) = h
    end do

    g = uniform_01_sample ( iseed )

    if ( 0.5E+00 <= g ) then
      v = - 1.0E+00
    else
      v = + 1.0E+00
    end if

    b(i) = s + n * uniform_01_sample ( iseed ) * 0.02E+00 * v

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
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
subroutine givr_l2 ( a, m, n, b, eps, x, iflag )

!*****************************************************************************80
!
!! GIVR_L2 minimizes the L2 norm of A*x-b using fast Givens rotations.
!
!  Modified:
!
!    22 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 35,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N) contains the M by N system matrix.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for the accuracy test.
!
!    Output, real ( kind = 4 ) X(N), the calculated solution.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    1, the accuracy test failed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) al
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) bi
  real ( kind = 4 ) dp
  real ( kind = 4 ) eps
  real ( kind = 4 ) f(n)
  real ( kind = 4 ) fh
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nn
  real ( kind = 4 ) r(((n-1)*n)/2)
  real ( kind = 4 ) rl
  real ( kind = 4 ) t(n)
  real ( kind = 4 ) wh
  real ( kind = 4 ) wi
  real ( kind = 4 ) x(n)

  iflag = 0

  n2 = ( ( n - 1 ) * n ) / 2
  f(1:n) = 0.0E+00
  t(1:n) = 0.0E+00
  r(1:n2) = 0.0E+00

  do i = 1, m

    bi = b(i)
    wi = 1.0E+00

    do k = 1, n

      if ( wi == 0.0E+00 ) then
        exit
      end if

      if ( a(i,k) /= 0.0E+00 ) then

        dp = f(k) + wi * a(i,k)**2

        if ( abs ( dp ) <= eps ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GIVR_L2 - Fatal error!'
          write ( *, '(a)' ) '  The pivot value is too small.'
          iflag = 1
          return
        end if

        fh = f(k) / dp
        wh = wi * a(i,k) / dp
        wi = wi * f(k) / dp

        f(k) = dp

        nn = ( ( k - 1 ) * ( n + n - k ) ) / 2 + 1

        do l = k+1, n
          al = a(i,l)
          rl = r(nn)
          a(i,l) = al - a(i,k) * rl
          r(nn) = fh * rl + wh * al
          nn = nn + 1
        end do

        al = bi
        bi = al - a(i,k) * t(k)
        t(k) = fh * t(k) + wh * al

      end if

    end do

  end do

  do k = n, 1, -1

    x(k) = t(k)

    nn = ( ( k - 1 ) * ( n + n - k ) ) / 2 + 1
    do l = k+1, n
      x(k) = x(k) - r(nn) * x(l)
      nn = nn + 1
    end do

  end do

  return
end
subroutine hfti ( a, m, n, b, mdb, nb, eps, rank, res_norm, ip )

!*****************************************************************************80
!
!! HFTI minimizes the L2 norm of A*x-b using Householder transformations.
!
!  Discussion:
!
!    The routine uses Householder forward triangulation with column
!    interchanges.
!
!    The array B is overburdened with the task of being both an
!    input and output quantity, of different dimensions, or none at all
!    (if NB = 0).
!
!  Modified:
!
!    17 May 2002
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974,
!    QA275.L38
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 38,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(M,N).  On input, the M by N system 
!    matrix.  On output, the information in A has been overwritten by data
!    used to solve the system.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input/output, real ( kind = 4 ) B(MDB,NB), on input, the NB right hand 
!    side vectors of length M, and on output, the NB least squares solution 
!    vectors of length N.
!
!    Input, integer ( kind = 4 ) MDB, the leading dimension of B, which must
!    satisfy max ( M, N ) <= MDB.
!
!    Input, integer ( kind = 4 ) NB, the number of systems to be solved.
!
!    Input, real ( kind = 4 ) EPS, a tolerance used in the pseudorank test.
!
!    Output, integer ( kind = 4 ) RANK, the numerically determined rank of the
!    matrix A.
!
!    Output, real ( kind = 4 ) RES_NORM(NB), the Euclidean norms of the 
!    residual vectors for the least squares problems.
!
!    Output, integer ( kind = 4 ) IP(N), the pivot vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mdb
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(mdb,nb)
  real ( kind = 4 ) eps
  real ( kind = 4 ), parameter :: factor = 0.001E+00
  real ( kind = 4 ) g(n)
  real ( kind = 4 ) h(n)
  real ( kind = 4 ) hmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ldiag
  integer ( kind = 4 ) lmax
  real ( kind = 4 ) r4_diff
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res_norm(nb)
  real ( kind = 4 ) sm
  real ( kind = 4 ) tmp

  k = 0
  hmax = 0.0E+00
  ldiag = min ( m, n )

  if ( ldiag <= 0 ) then
    rank = 0
    return
  end if

  do j = 1, ldiag
!
!  Update squared column lengths and find LMAX.
!
    if ( 1 < j ) then

      lmax = j
      do l = j, n
        h(l) = h(l) - a(j-1,l)**2
        if ( h(lmax) < h(l) ) then
          lmax = l
        end if
      end do

      if ( 0.0E+00 < r4_diff ( hmax + factor * h(lmax), hmax ) ) then
        go to 10
      end if

    end if
!
!  Set LMAX to the index of the column with largest L2 norm.
!
    lmax = j

    do l = j, n

      h(l) = sum ( a(j:m,l)**2 )

      if ( h(lmax) < h(l) ) then
        lmax = l
      end if

    end do

    hmax = h(lmax)
!
!  Do column interchanges if needed.
!
10  continue

    ip(j) = lmax

    if ( ip(j) /= j ) then

      do i = 1, m
        call r4_swap ( a(i,j), a(i,lmax) )
      end do

      h(lmax) = h(j)

    end if
!
!  Compute the J-th transformation and apply it to A and B.
!
!  Note that A(1,J+1) can be an illegal memory reference if J = N,
!  but in that case, A(1,J+1) can be replaced by anything, since
!  it won't be used.
!
    if ( j < n ) then
      call h12 ( 1, j, j+1, m, a(1,j), 1, h(j), a(1,j+1), 1, m, n-j )
    else
      call h12 ( 1, j, j+1, m, a(1,j), 1, h(j), a(1,1), 1, m, n-j )
    end if

    call h12 ( 2, j, j+1, m, a(1,j), 1, h(j), b, 1, mdb, nb )

  end do
!
!  Determine the pseudorank K using the tolerance EPS.
!
  k = ldiag

  do j = 1, ldiag
    if ( abs ( a(j,j) ) <= eps ) then
      k = j - 1
      exit
    end if
  end do
!
!  Compute the norms of the residual vectors.
!
  do jb = 1, nb
    res_norm(jb) = sqrt ( sum ( b(k+1:m,jb)**2 ) )
  end do
!
!  Special for pseudorank = 0.
!
  if ( k <= 0 ) then

    b(1:n,1:nb) = 0.0E+00

    rank = k
    return

  end if
!
!  If the pseudorank is less than N, compute the Householder
!  decomposition of the first K rows.
!
  if ( k /= n ) then
    do i = k, 1, -1
      call h12 ( 1, i, k+1, n, a(i,1), m, g(i), a, m, 1, i-1 )
    end do
  end if

  do jb = 1, nb
!
!  Solve the K by K triangular system.
!
    do i = k, 1, -1
      sm = dot_product ( a(i,i+1:k), b(i+1:k,jb) )
      b(i,jb) = ( b(i,jb) - sm ) / a(i,i)
    end do
!
!  Complete the computation of the solution vector.
!
    if ( k < n ) then

!
!  DEBUG: DIMENSION 1 OUT OF RANGE.
!
      if ( k+1 < 1 .or. n < k+1 .or. n < 1 .or. mdb < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HFTI - DEBUG:'
        write ( *, '(a,i6)' ) 'MDB = ', mdb
        write ( *, '(a,i6)' ) 'K+1 = ', k
        write ( *, '(a,i6)' ) 'N   = ', n
        write ( *, '(a,i6)' ) 'JB  = ', jb
        stop
      end if

      b(k+1:n,jb) = 0.0E+00

      do i = 1, k
        call h12 ( 2, i, k+1, n, a(i,1), m, g(i), b(1,jb), 1, mdb, 1 )
      end do

    end if
!
!  Reorder the solution vector to compensate for the column interchanges.
!
    do j = ldiag, 1, -1

      if ( ip(j) /= j ) then
        l = ip(j)
        call r4_swap ( b(l,jb), b(j,jb) )
      end if

    end do

  end do
!
!  The solution vectors X are now in the first N rows of the array B.
!
  rank = k

  return
end
subroutine hfti_l2 ( a, m, n, b, eps, x )

!*****************************************************************************80
!
!! HFTI_L2 minimizes the L2 norm of A*x-b using Householder transformations.
!
!  Discussion:
!
!    The routine uses Householder forward triangulation with column
!    interchanges.
!
!    I split off this version, which is called by the user directly,
!    from another version which is the workhorse of several routines,
!    because of the awkward, brittle and confusing interface.
!
!  Modified:
!
!    17 May 2002
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974,
!    QA275.L38
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 38,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(M,N).  On input, the M by N system 
!    matrix.  On output, the information in A has been overwritten by data
!    used to solve the system.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input/output, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input, real ( kind = 4 ) EPS, a tolerance used in the pseudorank test.
!
!    Output, real ( kind = 4 ) X(N), the solution.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) eps
  real ( kind = 4 ), parameter :: factor = 0.001E+00
  real ( kind = 4 ) g(n)
  real ( kind = 4 ) h(n)
  real ( kind = 4 ) hmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ldiag
  integer ( kind = 4 ) lmax
  real ( kind = 4 ) r4_diff
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) sm
  real ( kind = 4 ) tmp
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) xt(n+m)

  k = 0
  hmax = 0.0E+00
  ldiag = min ( m, n )
  xt(1:m) = b(1:m)

  if ( ldiag <= 0 ) then
    return
  end if

  do j = 1, ldiag
!
!  Update squared column lengths and find LMAX.
!
    if ( 1 < j ) then

      lmax = j
      do l = j, n
        h(l) = h(l) - a(j-1,l)**2
        if ( h(lmax) < h(l) ) then
          lmax = l
        end if
      end do

      if ( 0.0E+00 < r4_diff ( hmax + factor * h(lmax), hmax ) ) then
        go to 10
      end if

    end if
!
!  Set LMAX to the index of the column with largest L2 norm.
!
    lmax = j

    do l = j, n

      h(l) = sum ( a(j:m,l)**2 )

      if ( h(lmax) < h(l) ) then
        lmax = l
      end if

    end do

    hmax = h(lmax)
!
!  Do column interchanges if needed.
!
10  continue

    ip(j) = lmax

    if ( ip(j) /= j ) then

      do i = 1, m
        call r4_swap ( a(i,j), a(i,lmax) )
      end do

      h(lmax) = h(j)

    end if
!
!  Compute the J-th transformation and apply it to A and B.
!
!  Note that A(1,J+1) can be an illegal memory reference if J = N,
!  but in that case, A(1,J+1) can be replaced by anything, since
!  it won't be used.
!
    if ( j < n ) then
      call h12 ( 1, j, j+1, m, a(1,j), 1, h(j), a(1,j+1), 1, m, n-j )
    else
      call h12 ( 1, j, j+1, m, a(1,j), 1, h(j), a(1,1), 1, m, n-j )
    end if

    call h12 ( 2, j, j+1, m, a(1,j), 1, h(j), xt, 1, 1, 1 )

  end do
!
!  Determine the pseudorank K using the tolerance EPS.
!
  k = ldiag

  do j = 1, ldiag
    if ( abs ( a(j,j) ) <= eps ) then
      k = j - 1
      exit
    end if
  end do
!
!  Compute the norms of the residual vectors.
!
  res_norm = sqrt ( sum ( xt(k+1:m)**2 ) )
!
!  Special for pseudorank = 0.
!
  if ( k <= 0 ) then
    x(1:n) = 0.0E+00
    return
  end if
!
!  If the pseudorank is less than N, compute the Householder
!  decomposition of the first K rows.
!
  if ( k /= n ) then
    do i = k, 1, -1
      call h12 ( 1, i, k+1, n, a(i,1), m, g(i), a, m, 1, i-1 )
    end do
  end if
!
!  Solve the K by K triangular system.
!
  do i = k, 1, -1
    sm = dot_product ( a(i,i+1:k), xt(i+1:k) )
    xt(i) = ( xt(i) - sm ) / a(i,i)
  end do
!
!  Complete the computation of the solution vector.
!
  if ( k < n ) then

    xt(k+1:n) = 0.0E+00

    do i = 1, k
      call h12 ( 2, i, k+1, n, a(i,1), m, g(i), xt, 1, 1, 1 )
    end do

  end if
!
!  Reorder the solution vector to compensate for the column interchanges.
!
  do j = ldiag, 1, -1

    if ( ip(j) /= j ) then
      l = ip(j)
      call r4_swap ( xt(l), xt(j) )
    end if

  end do

  x(1:n) = xt(1:n)

  return
end
subroutine h12 ( mode, lpivot, l1, m, u, iue, up, c, ice, icv, ncv )

!*****************************************************************************80
!
!! H12 constructs or applies a Householder transformation.
!
!  Discussion:
!
!    A Householder transformation H is specified by a pivot vector U,
!    and has the form:
!
!      H = I + U * U' / S
!
!    where S is a scalar.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974,
!    QA275.L38
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MODE:
!    1, constructs a Householder transformation and applies it to
!      vectors in C.
!    2, applies a previously constructed transformation to the vectors
!      in C.
!
!    Input, integer ( kind = 4 ) LPIVOT, the index of the pivot element.
!
!    Input, integer ( kind = 4 ) L1, the transformation will be constructed to
!    zero out elements indexed from L1 through M.  L1 should satisfy
!    1 < L1 <= M.
!
!    Input, integer ( kind = 4 ) M, the column dimension of U.
!
!    Input/output, real ( kind = 4 ) U(*).  If MODE = 1, then on input U 
!    contains the pivot vector, and on output will contain information defining 
!    the Householder transformation.  If MODE = 2, then U should contain
!    the data computed on a previous call with MODE = 1.
!
!    Input, integer ( kind = 4 ) IUE, the storage increment in U between
!    successive elements.
!
!    Input/output, real ( kind = 4 ) UP.  On output from a call with MODE = 1, 
!    UP contains data used for the Householder transformation.  This same
!    value must be the input value of UP for a call with MODE = 2.
!
!    Input/output, real ( kind = 4 ) C(*), contains NCV vectors, each of 
!    length M, to which the Householder transformation is to be applied.
!    The exact arrangement of data in C is specified by ICE and ICV.
!
!    Input, integer ( kind = 4 ) ICE, the storage increment between
!    succesive elements of vectors in C.
!
!    Input, integer ( kind = 4 ) ICV, the storage increment between the first
!    elements of succesive vectors in C.
!
!    Input, integer ( kind = 4 ) NCV, the number of vectors in C to be 
!    transformed.  If NCV <= 0, no operations will be done on C.
!
  implicit none

  integer ( kind = 4 ) iue
  integer ( kind = 4 ) m

  real ( kind = 4 ) b
  real ( kind = 4 ) c(*)
  real ( kind = 4 ) cl
  real ( kind = 4 ) clinv
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ice
  integer ( kind = 4 ) icv
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lpivot
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) ncv
  real ( kind = 4 ) sm
  real ( kind = 4 ) u(iue,m)
  real ( kind = 4 ) up

  if ( lpivot <= 0 ) then
    return
  else if ( l1 <= lpivot ) then
    return
  else if ( m < l1 ) then
    return
  end if

  cl = abs ( u(1,lpivot) )

  if ( mode == 2 ) then

    if ( cl <= 0.0E+00 ) then
      return
    end if

  else
!
!  Construct the transformation.
!
    do j = l1, m
      cl = max ( cl, abs ( u(1,j) ) )
    end do

    if ( cl <= 0.0E+00 ) then
      return
    end if

    clinv = 1.0E+00 / cl
    sm = ( u(1,lpivot) * clinv )**2

    do j = l1, m
      sm = sm + ( u(1,j) * clinv )**2
    end do

    cl = cl * sqrt ( sm )

    if ( 0.0E+00 < u(1,lpivot) ) then
      cl = - cl
    end if

    up = u(1,lpivot) - cl
    u(1,lpivot) = cl

  end if
!
!  Apply the transformation I+U*U'/B to C.
!
  if ( ncv <= 0 ) then
    return
  end if

  b = up * u(1,lpivot)
!
!  B must be nonpositive here.
!
  if ( 0.0E+00 <= b ) then
    return
  end if

  b = 1.0E+00 / b
  i2 = 1 - icv + ice * ( lpivot - 1 )
  incr = ice * ( l1 - lpivot )

  do j = 1, ncv

    i2 = i2 + icv
    i3 = i2 + incr
    i4 = i3

    sm = c(i2) * up
    do i = l1, m
      sm = sm + c(i3) * u(1,i)
      i3 = i3 + ice
    end do

    if ( sm /= 0.0E+00 ) then

      sm = sm * b
      c(i2) = c(i2) + sm * up
      do i = l1, m
        c(i4) = c(i4) + sm * u(1,i)
        i4 = i4 + ice
      end do

    end if

  end do

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
subroutine i4_random ( ilo, ihi, i )

!*****************************************************************************80
!
!! I4_RANDOM returns a random integer in a given range.
!
!  Modified:
!
!    23 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ILO, IHI, the minimum and maximum
!    acceptable values.
!
!    Output, integer ( kind = 4 ) I, the randomly chosen integer.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  real ( kind = 4 ) r
  real ( kind = 4 ), parameter :: rhi = 1.0E+00
  real ( kind = 4 ), parameter :: rlo = 0.0E+00

  call r4_random ( rlo, rhi, r )

  i = ilo + int ( r * real ( ihi + 1 - ilo, kind = 4 ) )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two integer values.
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
subroutine i4_to_s_zero ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_ZERO converts an integer to a string, with zero padding.
!
!  Examples:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  000001
!        -1  -00001
!         0  000000
!      1952  001952
!    123456  123456
!   1234567  ******  <-- Not enough room!
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
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right justified, and zero padded.
!    If there is not enough space, the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  Working from right to left, strip off the digits of the integer
!  and place them into S(ILO:IHI).
!
  ipos = ihi

  do while ( ival /= 0 .or. ipos == ihi )

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

  end do
!
!  Fill the empties with zeroes.
!
  do i = ilo, ipos
    s(i:i) = '0'
  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
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
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine icmgs_l2 ( a, m, n, b, eps, x, iflag )

!*****************************************************************************80
!
!! ICMGS_L2 uses modifed Gram Schmidt on a problem with nonzero intercept.
!
!  Discussion:
!
!    ICMGS_L2 minimizes the L2 norm || (E,A)*(x0,x)' - b ||.
!
!    Here E is the vector of all 1's, and x0 is the scalar intercept variable.
!
!    This routine can only be used for models with an intercept.
!    The matrix A and vector B will be destroyed.
!
!  Modified:
!
!    30 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 31,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(M,N).  On input, A contains the system 
!    matrix.  On output, A contains information about the factorization of the 
!    matrix, which can be used to efficiently solve more systems.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input/output, real ( kind = 4 ) B(M), on input, the right hand side vector.
!    On output, B has been overwritten by other information.
!
!    Input, real ( kind = 4 ) EPS, accuracy tolerance for the orthogonalization.
!    If the squared euclidean length of the vector to be treated in the
!    K-th step is smaller than EPS, then the routine is halted with
!    IFLAG = 1.  Recommended value: EPS =  10**(-2*T+4).
!
!    Output, real ( kind = 4 ) X(N+1), the calculated solution.
!
!    Output, integer ( kind = 4 ) IFLAG.
!    0, no error occurred, the calculation was completed.
!    1, the calculation has halted because the system seemed singular.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) as(n)
  real ( kind = 4 ) b(m)
  logical, parameter :: bnew = .false.
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) k
  real ( kind = 4 ) s
  real ( kind = 4 ) x(n+1)

  iflag = 0

  do k = 1, n

    s = sum ( a(1:m,k) ) / real ( m, kind = 4 )

    as(k) = s
    a(1:m,k) = a(1:m,k) - s

  end do

  call mgs_l2 ( a, m, n, b, eps, bnew, x, iflag2 )

  if ( iflag2 /= 0 ) then
    iflag = iflag2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ICMGS_L2 - Fatal error!'
    write ( *, '(a,i6)' ) '  MGS_L2 returned nonzero IFLAG = ', iflag2
    return
  end if

  s = sum ( b(1:m) ) / real ( m, kind = 4 ) - dot_product ( as(1:n), x(1:n) )
!
!  Shift the X values.
!
  do k = n+1, 2, -1
    x(k) = x(k-1)
  end do

  x(1) = s

  return
end
subroutine inexcl ( n, ai, bi, wi, j, s, s2, eps, f, t, r, fj, tj, rj, ss, &
  iflag )

!*****************************************************************************80
!
!! INEXCL computes auxilliary arrays F, T and R used to control exchanges.
!
!  Modified:
!
!    29 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 151,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, ?
!
!    AI ?
!
!    BI ?
!
!    Input/output, real ( kind = 4 ) WI, ?
!
!    J, ?
!
!    S, ?
!
!    S2, ?
!
!    EPS, ?
!
!    F, ?
!
!    T, ?
!
!    R, ?
!
!    FJ, ?
!
!    TJ, ?
!
!    RJ, ?
!
!    SS, ?
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    ?
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s2

  real ( kind = 4 ) ai(n)
  real ( kind = 4 ) ak
  real ( kind = 4 ) al
  real ( kind = 4 ) bi
  real ( kind = 4 ) dp
  real ( kind = 4 ) eps
  real ( kind = 4 ) f(n,s)
  real ( kind = 4 ) fh
  real ( kind = 4 ) fj(n)
  real ( kind = 4 ) fk
  real ( kind = 4 ) hk
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nn
  real ( kind = 4 ) r(s,s2)
  real ( kind = 4 ) rj(s2)
  real ( kind = 4 ) rl
  real ( kind = 4 ) ss
  real ( kind = 4 ) t(n,s)
  real ( kind = 4 ) tj(n)
  real ( kind = 4 ) tk
  real ( kind = 4 ) wa
  real ( kind = 4 ) wh
  real ( kind = 4 ) wi

  iflag = 0

  do k = 1, n

    if ( wi == 0.0E+00 ) then

      if ( k /= n ) then

        nn = ( (k-1) * (n+n-k) ) / 2 + 1

        do l = nn, s2
          rj(l) = r(j,l)
        end do

      end if

      fj(k:n) = f(k:n,j)
      tj(k:n) = t(k:n,j)

      return

    end if

    ak = ai(k)

    if ( ak == 0.0E+00 ) then

      nn = ( (k-1) * (n+n-k) ) / 2 + 1

      do l = k+1, n
        rj(nn) = r(j,nn)
        nn = nn + 1
      end do

      fj(k) = f(k,j)
      tj(k) = t(k,j)

    else

      fk = f(k,j)
      wa = wi * ak
      dp = fk + wa * ak

      if ( abs ( dp ) <= eps ) then
        iflag = 5
        return
      end if

      hk = 1.0E+00 / dp
      fh = fk * hk
      wh = wa * hk
      wi = wi * fh
      fj(k) = dp

      nn = ( (k-1) * (n+n-k) ) / 2 + 1

      do l = k+1, n
        al = ai(l)
        rl = r(j,nn)
        ai(l) = al - ak * rl
        rj(nn) = fh * rl + wh * al
        nn = nn + 1
      end do

      al = bi
      tk = t(k,j)
      bi = al - ak * tk
      tj(k) = fh * tk + wh * al

    end if

  end do

  ss = ss + wi * bi * bi

  return
end
subroutine ldp_l2 ( a, m, n, b, x, xnorm, iflag )

!*****************************************************************************80
!
!! LDP_L2 implements least distance programming algorithm.
!
!  Discussion:
!
!    The problem is to determine an N vector X with minimum L2 norm,
!    subject to
!
!      A * X >= B
!
!    where A is an M by N matrix, and B is an M vector.
!
!  Modified:
!
!    14 May 2002
!
!  Reference:
!
!    Charles Lawson and Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974,
!    Revised edition, SIAM, 1995.
!    QA275.L38
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the M by N matrix.
!    There is no restriction on the rank of the matrix.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in G.
!    M may be greater or less than N.  N must be greater than 0.
!
!    Input, real ( kind = 4 ) B(M), the right hand side of the constraint.
!
!    Output, real ( kind = 4 ) X(N), the computed solution.
!
!    Output, real ( kind = 4 ) XNORM, the value of the L2 norm of X.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error.
!    1, M <= 0, meaning no constraints, so X = 0 is solution.
!    2, N <= 0.
!    3, maximum number of iterations taken without convergence.
!    4, the inequality constraints A * X >= B are incompatible.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) f(n+1)
  real ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jf
  real ( kind = 4 ) r4_diff
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) w(n+1,m)
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) xnorm
  real ( kind = 4 ) y(m)

  iflag = 0

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LDP_L2 - Fatal error!'
    write ( *, '(a)' ) '  N <= 0.'
    iflag = 2
    return
  end if

  x(1:n) = 0.0E+00
  xnorm = 0.0E+00

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LDP_L2 - Warning!'
    write ( *, '(a)' ) '  M <= 0.'
    iflag = 1
    return
  end if
!
!  Find Y that minimizes the L2 norm of
!
!    ( A' ) * Y - ( 0 )
!    ( B' )       ( 1 )
!
!  with 0 <= Y.
!
  w(1:n,1:m) = transpose ( a(1:m,1:n) )
  w(n+1,1:m) = b(1:m)

  f(1:n) = 0.0E+00
  f(n+1) = 1.0E+00

  call nn_l2 ( w, n+1, m, f, y, res_norm, iflag2 )

  if ( iflag2 /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LDP_L2 - Fatal error!'
    write ( *, '(a)' ) '  NN_L2 returned error flag IFLAG2 = ', iflag2
    iflag = 4
    return
  end if

  if ( res_norm <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LDP_L2 - Fatal error!'
    write ( *, '(a)' ) '  NN_L2 reports incompatibility in constraints.'
    iflag = 4
    return
  end if

  fac = 1.0E+00 - dot_product ( b(1:m), y(1:m) )

  if ( r4_diff ( 1.0E+00 + fac, 1.0E+00 ) == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LDP_L2 - Fatal error!'
    write ( *, '(a)' ) '  Nearly singular system.'
    iflag = 4
    return
  end if

  x(1:n) = matmul ( y(1:m), a(1:m,1:n) ) / fac

  xnorm = sqrt ( sum ( x(1:n)**2 ) )

  return
end
subroutine mgs_l2 ( a, m, n, b, eps, bnew, x, iflag )

!*****************************************************************************80
!
!! MGS_L2 minimizes the L2 norm of A*x-b using the modified Gram-Schmidt method.
!
!  Discussion:
!
!    Several right-hand sides B may successively be processed.
!
!    ...Actually, that's not true unless we save R, which we don't anymore!
!
!  Modified:
!
!    26 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 28,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(M,N).  On input, A contains the system 
!    matrix.  On output, A contains information about the factorization of the 
!    matrix, which can be used to efficiently solve more systems.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input/output, real ( kind = 4 ) B(M), on input, the right hand side vector.
!    On output, B has been overwritten by other information.
!
!    Input, real ( kind = 4 ) EPS, accuracy tolerance for the orthogonalization.
!    If the squared euclidean length of the vector to be treated in the
!    K-th step is smaller than EPS, then the routine is halted with
!    IFLAG = 1.  Recommended value: EPS =  10**(-2*T+4).
!
!    Input, logical BNEW.  For BNEW = .false. the method is applied to
!    (A, B), and the least squares solution X is calculated.  After such
!    a call, using the factored A and R from the previous call, and
!    setting BNEW  = .true., X can easily be determined for another
!    right-hand side B without having to recompute the factorization of A.
!
!    Output, real ( kind = 4 ) X(N), the calculated solution.
!
!    Output, integer ( kind = 4 ) IFLAG.
!    0, no error occurred, the calculation was completed.
!    1, the calculation has halted because the system seemed singular.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  logical bnew
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) r(n,n)
  real ( kind = 4 ) x(n)
!
!  Carry out the Gram Schmidt orthogonalization.
!
!  R(K,K) will contain the reciprocal values of the squared lengths of
!  the orthogonal vectors.
!
!  Factoring out D = 1 / DIAG(R) we then have
!    A' * A = R' * D * R
!  where now the R matrix has unit diagonal.
!
  iflag = 0

  do k = 1, n

    if ( .not. bnew ) then
      r(k,k) = 1.0E+00 / sum ( a(1:m,k)**2 )
    end if
!
!  Process the right hand side.
!
    x(k) = r(k,k) * dot_product ( b(1:m), a(1:m,k) )

    if ( k /= n ) then

      b(1:m) = b(1:m) - a(1:m,k) * x(k)

      if ( .not. bnew ) then

        do j = k+1, n

          r(k,j) = dot_product ( a(1:m,k), a(1:m,j) ) * r(k,k)

          a(1:m,j) = a(1:m,j) - a(1:m,k) * r(k,j)

        end do

      end if

    end if

  end do
!
!  Back substitution.
!
  do i = n-1, 1, -1
    x(i) = x(i) - dot_product ( r(i,i+1:n), x(i+1:n) )
  end do

  return
end
subroutine nn_l1 ( m, l, k, n, a, b, c, d, e, h, eps, itmax, x, res, &
  res_norm, iflag )

!*****************************************************************************80
!
!! NN_L1 minimizes the L1 norm of A * X - B with linear constraints and X >=0.
!
!  Discussion:
!
!    The constraints have the form:
!
!      C * X  = D
!      E * X >= F
!          X >= 0
!
!  Modified:
!
!    16 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 252,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) L, the number of rows of C.
!
!    Input, integer ( kind = 4 ) K, the number of rows of E.
!
!    Input, integer ( kind = 4 ) N, the number of unknowns, and the number of
!    columns in A, C and E.
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) C(L,N), the equality constraint matrix.
!
!    Input, real ( kind = 4 ) D(L), the equality right hand side.
!
!    Input, real ( kind = 4 ) E(K,N), the inequality constraint matrix.
!
!    Input, real ( kind = 4 ) H(K), the inequality right hand side.
!
!    Input, real ( kind = 4 ) EPS, a tolerance for an accuracy test.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations to take.
!
!    Output, real ( kind = 4 ) X(N+2), contains the solution in entries 1 
!    through N.
!
!    Output, real ( kind = 4 ) RES(M+L+K), contains the residuals.
!
!    Output, real ( kind = 4 ) RES_NORM, the L1 norm of the residual.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, an optimal solution was found.
!    1, there is no feasible solution.
!    2, premature exit because of rounding errors.
!    3, maximal number of iterations taken.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) c(l,n)
  integer ( kind = 4 ) code
  real ( kind = 4 ) d(l)
  real ( kind = 4 ) e(k,n)
  real ( kind = 4 ) eps
  real ( kind = 4 ) h(k)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) itmax
  real ( kind = 4 ) res(m+l+k)
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) x(n)
!
!  By setting CODE to 1, we can place sign constraints on the solution.
!
  code = 1
  x(1:n) = 1.0E+00
  res(1:m+l+k) = 0.0E+00

  call con_l1 ( m, l, k, n, a, b, c, d, e, h, code, eps, itmax, x, &
    res, res_norm, iflag )

  return
end
subroutine nn_l2 ( a, m, n, b, x, res_norm, iflag )

!*****************************************************************************80
!
!! NN_L2 minimizes the L2 norm of A * X - B with X >=0.
!
!  Discussion:
!
!    Given an M by N matrix A and an M vector B, compute an
!    N vector X that solves the least squares problem (that is,
!    minimizes the L2 norm of the residual)
!
!      A * X = B
!
!    subject to
!
!      X >= 0.
!
!    Householder transformations are applied, which essentially
!    multiply the system by an orthogonal matrix Q:
!
!      Q * A * X = Q * B
!
!  Modified:
!
!    03 May 2002
!
!  Reference:
!
!    Charles Lawson and Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974,
!    Revised edition, SIAM, 1995.
!    QA275.L38
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(M,N).  On input, the M by N matrix.
!    On output, the product matrix Q * A, where Q is an M by M orthogonal
!    matrix generated implicitly by this routine.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input/output, real ( kind = 4 ) B(M).  On input, contains the right
!    hand side vector.  On output, contains Q * B.
!
!    Output, real ( kind = 4 ) X(N), the computed solution.
!
!    Output, real ( kind = 4 ) RES_NORM, the L2 or Euclidean norm of 
!    the residual.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, the solution has been computed successfully.
!    1, the dimensions of the problem are bad, M <= 0 or N <= 0.
!    2, iteration count exceeded.  More than 3*N iterations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) alpha
  real ( kind = 4 ) asave
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) cc
  real ( kind = 4 ) dummy(1)
  real ( kind = 4 ), parameter :: factor = 0.01E+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) index(n)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) iz1
  integer ( kind = 4 ) iz2
  integer ( kind = 4 ) izmax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) l
  integer ( kind = 4 ) npp1
  integer ( kind = 4 ) nsetp
  real ( kind = 4 ) r4_diff
  real ( kind = 4 ) res_norm
  integer ( kind = 4 ) rtnkey
  real ( kind = 4 ) sm
  real ( kind = 4 ) ss
  real ( kind = 4 ) t
  real ( kind = 4 ) temp
  real ( kind = 4 ) unorm
  real ( kind = 4 ) up
  real ( kind = 4 ) w(n)
  real ( kind = 4 ) wmax
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) ztest
  real ( kind = 4 ) zz(m)

  iflag = 0

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NN_L2 - Fatal error!'
    write ( *, '(a)' ) '  M <= 0.'
    iflag = 1
    return
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NN_L2 - Fatal error!'
    write ( *, '(a)' ) '  N <= 0.'
    iflag = 1
    return
  end if

  it = 0
  itmax = 3 * n
!
!  Initialize INDEX and X.
!
  x(1:n) = 0.0E+00

  call i4vec_indicator ( n, index )

  iz2 = n
  iz1 = 1
  nsetp = 0
  npp1 = 1
!
!  The main loop begins here.
!
10 continue
!
!  Quit if all coefficients are already in the solution.
!  or if M columns of A have been triangularized.
!
  if ( iz2 < iz1 .or. m <= nsetp ) then
    go to 70
  end if
!
!  Compute components of the dual (negative gradient) vector W.
!
  do iz = iz1, iz2
    j = index(iz)
    w(j) = dot_product ( b(npp1:m), a(npp1:m,j) )
  end do
!
!  Find the largest positive W.
!
  do

    wmax = 0.0E+00
    do iz = iz1, iz2
      j = index(iz)
      if ( wmax < w(j) ) then
        wmax = w(j)
        izmax = iz
      end if
    end do
!
!  If WMAX <= 0, go to termination.
!  This indicates satisfaction of the Kuhn-Tucker conditions.
!
    if ( wmax <= 0.0E+00 ) then
      go to 70
    end if

    iz = izmax
    j = index(iz)
!
!  The sign of W(J) is OK for J to be moved to set P.
!
!  Begin the transformation and check new diagonal element to avoid
!  near linear dependence.
!
    asave = a(npp1,j)

    call h12 ( 1, npp1, npp1+1, m, a(1,j), 1, up, dummy, 1, 1, 0 )

    unorm = sqrt ( sum ( a(1:nsetp,j)**2 ) )

    if ( 0.0E+00 < r4_diff ( unorm + abs ( a(npp1,j) ) * factor, unorm ) ) then
!
!  Column J is sufficiently independent.
!
!  Copy B into ZZ, update ZZ, and solve for ZTEST, the proposed
!  new value for X(N).
!
      zz(1:m) = b(1:m)

      call h12 ( 2, npp1, npp1+1, m, a(1,j), 1, up, zz, 1, 1, 1 )

      ztest = zz(npp1) / a(npp1,j)
!
!  See if ZTEST is positive.
!
      if ( 0.0E+00 < ztest ) then
        exit
      end if

    end if
!
!  Reject J as a candidate to be moved from set Z to set P.
!  Restore A(NPP1,J), set W(J) = 0, and loop back to test dual
!  coefficients again.
!
    a(npp1,j) = asave
    w(j) = 0.0E+00

  end do
!
!  The index J = INDEX(IZ)  has been selected to be moved from
!  set Z to set P.
!
!  Update B, update indices,  Apply Householder transformations to columns
!  in new set Z, zero subdiagonal elements in column J, set W(J) = 0.
!
  b(1:m) = zz(1:m)
  index(iz) = index(iz1)
  index(iz1) = j
  iz1 = iz1+1
  nsetp = npp1
  npp1 = npp1+1

  if ( iz1 <= iz2) then
    do jz = iz1,iz2
      jj = index(jz)
      call h12 ( 2, nsetp, npp1, m, a(1,j), 1, up, a(1,jj), 1, m, 1 )
    end do
  end if

  if ( nsetp /= m) then
    a(npp1:m,j) = 0.0E+00
  end if

  w(j) = 0.0E+00
!
!  Solve the triangular system.
!  Store the solution temporarily in ZZ.
!
  rtnkey = 1
  go to 80

40 continue
!
!  Secondary loop begins here.
!
!  Iteration counter.
!
50 continue

  if ( itmax <= it ) then
    iflag = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NN_L2 - Warning!'
    write ( *, '(a)' ) '  Number of iterations exceeded.'
    go to 70
  end if

  it = it + 1
!
!  See if all new constrained coefficients are feasible.
!  If not, compute ALPHA.
!
  alpha = 2.0E+00

  do ip = 1, nsetp
    l = index(ip)
    if ( zz(ip) <= 0.0E+00 ) then
      t = -x(l) / ( zz(ip) - x(l) )
      if ( t < alpha ) then
        alpha = t
        jj = ip
      end if
    end if
  end do
!
!  If all new constrained coefficients are feasible then ALPHA will
!  still be 2.  If so exit from secondary loop to main loop.
!
  if ( alpha == 2.0E+00 ) then

    do ip = 1, nsetp
      i = index(ip)
      x(i) = zz(ip)
    end do

    go to 10

  end if
!
!  Otherwise use ALPHA which will be between 0.0 and 1.0 to
!  interpolate between the old X and the new ZZ.
!
  do ip = 1, nsetp
    l = index(ip)
    x(l) = x(l) + alpha * ( zz(ip) - x(l) )
  end do
!
!  Modify A and B and the index arrays to move coefficient I
!  from set P to set Z.
!
  i = index(jj)

60 continue

  x(i) = 0.0E+00

  if ( jj /= nsetp ) then

    jj = jj + 1

    do j = jj, nsetp

      ii = index(j)
      index(j-1) = ii
      call g1 ( a(j-1,ii), a(j,ii), cc, ss, a(j-1,ii) )
!
!  Apply procedure g2 (cc,ss,a(j-1,l),a(j,l))
!
      a(j,ii) = 0.0E+00

      do l = 1, n
        if ( l /= ii ) then
          temp = a(j-1,l)
          a(j-1,l) =   cc * temp + ss * a(j,l)
          a(j,l)    = -ss * temp + cc * a(j,l)
        end if
      end do
!
!  Apply procedure g2 (cc,ss,b(j-1),b(j))
!
      temp = b(j-1)
      b(j-1) =   cc * temp + ss * b(j)
      b(j)    = -ss * temp + cc * b(j)

    end do

  end if

  npp1 = nsetp
  nsetp = nsetp - 1
  iz1 = iz1 - 1
  index(iz1) = i
!
!  See if the remaining coefficients in set P are feasible.  They should
!  be because of the way ALPHA was determined.
!  If any are infeasible it is due to round-off error. Any
!  that are nonpositive will be set to zero
!  and moved from set P to set Z.
!
  do jj = 1, nsetp
    i = index(jj)
    if ( x(i) <= 0.0E+00 ) then
      go to 60
    end if
  end do
!
!  Copy B into ZZ, then solve again and loop back.
!
  zz(1:m) = b(1:m)

  rtnkey = 2
  go to 80
!
!  End of secondary loop.
!
!
!  End of the main loop
!
!  Come here for termination.
!  Compute the norm of the final residual vector.
!
70 continue

  sm = sum ( b(npp1:m)**2 )

  if ( m < npp1 ) then
    w(1:n) = 0.0E+00
  end if

  res_norm = sqrt ( sm )
  return
!
!  The following block of code is used as an internal subroutine
!  to solve the triangular system, putting the solution in ZZ.
!
80 continue

  do ip = nsetp, 1, -1

    if ( ip /= nsetp ) then
      do ii = 1, ip
        zz(ii) = zz(ii) - a(ii,jj) * zz(ip+1)
      end do
    end if

    jj = index(ip)
    zz(ip) = zz(ip) / a(ip,jj)

  end do
!
!  Something wrong here...Maybe 50 should be 330???
!
  if ( rtnkey == 1 ) then
    go to 40
  else
    go to 50
  end if

end
subroutine nn_li ( k, l, m, n, q, h, g, eps, x, res_norm, iflag )

!*****************************************************************************80
!
!! NN_LI minimizes the L-infinity norm of A*X-B, linear constraints and 0<=X.
!
!  Modified:
!
!    07 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 254,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, ?
!
!    Input, integer ( kind = 4 ) L, ?
!
!    Input, integer ( kind = 4 ) M, ?
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    ?, real Q(N+3,M+L+K+2), ?
!
!    Input/?, real H(M+L+K+2), ?
!
!    Output, real ( kind = 4 ) G(M), ?
!
!    Input, real ( kind = 4 ) EPS, ?
!
!    Output, real ( kind = 4 ) X(N+3), ?
!
!    Output, real ( kind = 4 ) RES_NORM, ?
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error detected.
!    1, the computation failed in CON_LI.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) big
  real ( kind = 4 ) eps
  real ( kind = 4 ) g(m)
  real ( kind = 4 ) h(m+l+k+2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 4 ) q(n+3,m+l+k+2)
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) x(n+3)

  big = maxval ( h(1:m) )

  big = big * 100.0E+00

  q(1:m,k+1:k+m) = 0.0E+00
  do i = 1, m
    q(i,k+i) = 1.0E+00
  end do

  h(k+1:k+m) = 0.0E+00
  g(1:m) = big
!
!  This call to CON_LI must be corrected.
!  The number of arguments isn't even right!.
!
! call con_li ( k, l, m, n, q, h, g, eps, x, res_norm, iflag )

  return
end
subroutine normal_l2 ( a, m, n, b, lambda, x )

!*****************************************************************************80
!
!! NORMAL_L2 minimizes the L2 norm of A*x-b using the normal equations.
!
!  Discussion:
!
!    The norm used is the L2 or Euclidean vector norm.
!
!    The coefficients of the normal equations are computed with double
!    precision, but are rounded to single precision before continuing
!    the calculation.
!
!    Given an overdetermined M by N linear system A * x = b, where
!    N < M, and A has maximal rank N, it is possible to define the
!    normal equations by multiplying through by the transpose:
!
!      A' * A * x = A' * b
!
!    This N by N system may then be solved via Gauss elimination (or,
!    since the system matrix is now symmetric, Cholesky factorization).
!    The resulting solution vector x is the minimizer of the L2 norm
!    of the residual r = A * x - b.
!
!    However, the method of normal equations suffers from several problems:
!
!    * the method fails, or becomes very inaccurate, if the rank of A
!      is less than N or "nearly" so;
!
!    * the condition number of the system matrix (A'*A) is the square
!      of the condition number of A.  Hence, the normal equations are
!      solved with much less accuracy than might be expected.
!
!  Modified:
!
!    22 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 23,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N) contains the M by N system matrix.  Note 
!    that this routine does not alter the matrix A in any way.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input, real ( kind = 4 ) LAMBDA, a quantity to be added to the diagonal 
!    of the normal equation matrix.  You may set this to zero.
!
!    Output, real ( kind = 4 ) X(N), the calculated solution.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) c(n,n)
  real ( kind = 4 ) d(n)
  real ( kind = 4 ) lambda
  real ( kind = 4 ) x(n)

  x(1:n) = 0.0E+00
!
!  Compute C = A'*A.
!
  c(1:n,1:n) = matmul ( transpose ( a(1:m,1:n) ), a(1:m,1:n) )
!
!  Add LAMBDA to the diagonal, so now C = A'*A + LAMBDA * I.
!
  call r4mat_diag_add_scalar ( n, c, lambda )
!
!  Compute D = A' * B.
!
  d(1:n) = matmul ( transpose ( a(1:m,1:n) ), b(1:m) )
!
!  Compute the Cholesky decomposition of C = L*L'.
!
  call r4mat_cholesky_factor ( n, c, c )
!
!  Solve C*X = L*L'*X = D.
!
  call r4mat_cholesky_solve ( n, c, d, x )

  return
end
subroutine npart_enum ( n, npart, npartitions )

!*****************************************************************************80
!
!! NPART_ENUM enumerates the number of partitions of N with NPART parts.
!
!  Modified:
!
!    23 January 1999
!
!  Reference:
!
!    Algorithm 3.6,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 74.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    Normally N must be positive, but for this routine any
!    N is allowed.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    Normally, 1 <= NPART <= N is required,
!    but for this routine any value of NPART is allowed.
!
!    Output, integer ( kind = 4 ) NPARTITIONS is the number of partitions of N
!    with NPART parts.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) npart
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) p(0:n,0:n)

  if ( n <= 0 ) then

    npartitions = 0

  else if ( npart <= 0 .or. n < npart ) then

    npartitions = 0

  else

    call npart_table ( n, npart, n, p )

    npartitions = p(n,npart)

  end if

  return
end
subroutine npart_rsf_lex_random ( n, npart, a )

!*****************************************************************************80
!
!! NPART_RSF_LEX_RANDOM returns a random RSF NPART partition.
!
!  Modified:
!
!    12 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Output, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) rank

  call npart_enum ( n, npart, npartitions )

  call i4_random ( 1, npartitions, rank )

  call npart_rsf_lex_unrank ( rank, n, npart, a )

  return
end
subroutine npart_rsf_lex_unrank ( rank, n, npart, a )

!*****************************************************************************80
!
!! NPART_RSF_LEX_UNRANK unranks an RSF NPART partition in the lex ordering.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Algorithm 3.9,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 78.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the partition.
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Output, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) npartcopy
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) p(0:n,0:npart)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rankcopy
!
!  Check.
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input N is illegal.'
    stop
  end if

  if ( npart < 1 .or. n < npart ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input NPART is illegal.'
    stop
  end if

  call npart_enum ( n, npart, npartitions )

  if ( rank < 0 .or. npartitions < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if
!
!  Get the table of partitions of N with NPART parts.
!
  call npart_table ( n, npart, n, p )

  a(1:npart) = 0

  rankcopy = rank
  ncopy = n
  npartcopy = npart

  do while ( 0 < ncopy )

    if ( rankcopy < p(ncopy-1,npartcopy-1) ) then
      a(npart+1-npartcopy) = a(npart+1-npartcopy) + 1
      ncopy = ncopy - 1
      npartcopy = npartcopy - 1
    else
      do i = 1, npartcopy
        a(npart+1-i) = a(npart+1-i) + 1
      end do
      rankcopy = rankcopy - p(ncopy-1,npartcopy-1)
      ncopy = ncopy - npartcopy
    end if

  end do

  return
end
subroutine npart_table ( n, npart, nmax, p )

!*****************************************************************************80
!
!! NPART_TABLE tabulates the number of partitions of N having NPART parts.
!
!  Modified:
!
!    17 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Algorithm 3.5,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 72.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input, integer ( kind = 4 ) NMAX, the leading dimension of P.
!
!    Output, integer ( kind = 4 ) P(0:NMAX,0:NPART), P(I,J) is the number of
!    partitions of I having J parts.
!
  implicit none

  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p(0:nmax,0:npart)

  p(0,0) = 1
  p(1:n,0) = 0

  do i = 1, n
    do j = 1, npart
      if ( i < j ) then
        p(i,j) = 0
      else if ( i < 2 * j ) then
        p(i,j) = p(i-1,j-1)
      else
        p(i,j) = p(i-1,j-1) + p(i-j,j)
      end if
    end do
  end do

  return
end
subroutine orth_l1 ( a, m, n, b, eps1, eps2, itmax, x, iflag )

!*****************************************************************************80
!
!! ORTH_L1 carries out orthogonal regression in the L1 norm.
!
!  Discussion:
!
!    For the matrix Z = ( A | -b ), this routine tries to find an
!    N+1 vector X with unit L2 norm which minimizes the L1 norm of
!
!      Z * X = A * X(1:N) - B * X(N+1).
!
!  Modified:
!
!    17 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 282,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) EPS1, a tolerance used by CON_L1.
!
!    Input, real ( kind = 4 ) EPS2, a convergence tolerance.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations.
!
!    Output, real ( kind = 4 ) X(N+1), the computed solution.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no errors detected.
!    1, the iteration limit was exceeded.
!    2, there was an error in CON_L1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) a2(m,n+1)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) b2(m)
  real ( kind = 4 ) c(1,n+1)
  integer ( kind = 4 ) code
  real ( kind = 4 ) d(1)
  real ( kind = 4 ) e(1,1)
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  real ( kind = 4 ) h(1)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) itmax2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) k
  real ( kind = 4 ) q(m+1,n+2)
  real ( kind = 4 ) res(m+1)
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) u
  real ( kind = 4 ) v
  real ( kind = 4 ) x(n+1)
  real ( kind = 4 ) x_old(n+1)

  iflag = 0

  it = 0
  itmax2 = 10 * ( m + 1 )

  x(1:n) = 0.0E+00
  x(n+1) = 1.0E+00

  do

    if ( itmax <= it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ORTH_L1 - Fatal error!'
      write ( *, '(a)' ) '  Number of iterations exceeded.'
      iflag = 1
      return
    end if

    it = it + 1

    x_old(1:n+1) = x(1:n+1)

    a2(1:m,1:n) = a(1:m,1:n)
    a2(1:m,n+1) = -b(1:m)

    b2(1:m) = 0.0E+00

    c(1,1:n+1) = x(1:n+1)
    d(1) = 1.0E+00

    code = 0
    l = 1
    k = 0
!
!  I believe that the number of variables here is N+1!
!
!    A -B  X  0
!    X  Y  Y  1
!
    call con_l1 ( m, l, k, n+1, a2, b2, c, d, e, h, code, eps1, itmax2, &
      x, res, res_norm, iflag2 )

    if ( iflag2 /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ORTH_L1 - Fatal error!'
      write ( *, '(a,i6)' ) '  CON_L1 returned nonzero error flag IFLAG2 = ', &
        iflag2
      iflag = 2
      return
    end if

    x(1:n+1) = x(1:n+1) / sqrt ( sum ( x(1:n+1)**2 ) )

    u = sum ( abs ( x(1:n+1) ) )
    v = sum ( abs ( x(1:n+1) - x_old(1:n+1) ) )

    if ( v <= eps2 * u ) then
      exit
    end if

  end do

  return
end
subroutine orth_l2 ( a, m, n, b, eps1, eps2, itmax, x, start, x_start, iflag )

!*****************************************************************************80
!
!! ORTH_L2 carries out orthogonal regression in the L2 norm.
!
!  Discussion:
!
!    For the matrix Z = ( A | -b ), this routine tries to find an N+1
!    vector X with unit L2 norm which minimizes the L2 norm of Z*X.
!
!    The power method for inverse ( Z * Z' ) is used.
!
!  Modified:
!
!    15 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 266-167,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) EPS1, a tolerance for an accuracy test in the 
!    modified Gram Schmidt algorithm.
!
!    Input, real ( kind = 4 ) EPS2, a tolerance for convergence of the iteration.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations.
!
!    Input, logical START, specifies if a starting point is given.
!    TRUE, then the user has supplied a starting vector in X_START;
!    FALSE, no starting vector has been supplied.
!
!    Output, real ( kind = 4 ) X(N+1), the solution.
!
!    Input, real ( kind = 4 ) X_START(N+1), should be set to a starting point 
!    if the input value of START is TRUE.  Otherwise, the program will set
!    X_START(1:N+1) to 1 for a starting point.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no errors detected.
!    1, the EPS1 accuracy test failed.
!    2, there was no convergence within ITMAX iterations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  real ( kind = 4 ) h
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) r(n+1,n+1)
  real ( kind = 4 ) s
  logical start
  real ( kind = 4 ) u
  real ( kind = 4 ) x(n+1)
  real ( kind = 4 ) x_start(n+1)
  real ( kind = 4 ) z(m,n+1)

  z(1:m,1:n) = a(1:m,1:n)
  z(1:m,n+1) = - b(1:m)

  it = 0
  iflag = 0

  do k = 1, n + 1

    s = sum ( ( z(1:m,k) )**2 )

    if ( s < eps1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ORTH_L2 - Fatal error!'
      write ( *, '(a)' ) '  A column has an L2 norm that is too small.'
      iflag = 1
      return
    end if

    h = 1.0E+00 / s

    r(k,k) = h

    do j = k+1, n + 1

      r(k,j) = h * dot_product ( z(1:m,k), z(1:m,j) )

      z(1:m,j) = z(1:m,j) - z(1:m,k) * r(k,j)

    end do

  end do

  if ( .not. start ) then
    x_start(1:n+1) = 1.0E+00
  end if

  do

    if ( itmax <= it .and. 1 < itmax ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ORTH_L2 - Fatal error!'
      write ( *, '(a)' ) '  Number of iterations exceeded.'
      iflag = 2
      exit
    end if

    it = it + 1

    do k = 1, n + 1
      x(k) = x_start(k) - dot_product ( x(1:k-1), r(1:k-1,k) )
    end do

    do k = n + 1, 1, -1
      x(k) = r(k,k) * x(k) - dot_product ( r(k,k+1:n+1), x(k+1:n+1) )
    end do

    s = maxval ( abs ( x(1:k) ) )

    u = sum ( abs ( x(1:n+1) - maxval ( abs ( x(1:n+1) ) ) * x_start(1:n+1) ) )

    x_start(1:n+1) = x(1:n+1) / maxval ( abs ( x(1:n+1) ) )

    if ( u <= eps2 * sum ( abs ( x(1:n+1) ) ) ) then
      exit
    end if

  end do

  x(1:n+1) = x_start(1:n+1) / sqrt ( sum ( x_start(1:n+1)**2 ) )

  return
end
subroutine orth_li ( a, m, n, b, eps1, eps2, itmax, x, iflag )

!*****************************************************************************80
!
!! ORTH_LI carries out orthogonal regression in the L-infinity norm.
!
!  Discussion:
!
!    THIS ROUTINE IS NOT PRODUCING SATISFACTORY RESULTS.
!
!    For the matrix Z = ( A | -b ), this routine tries to find a vector X
!    with unit L-infinity norm which minimizes the L-infinity norm of Z*X.
!
!  Modified:
!
!    30 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 287,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) EPS1, a tolerance used in CON_LI.
!
!    Input, real ( kind = 4 ) EPS2, a tolerance used for the convergence of
!    the iterations.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations.
!
!    Output, real ( kind = 4 ) X(N+1), the computed solution.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no errors were detected.
!    1, an error was detected in CON_LI.
!    2, the iterations did not converge.
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 0
  integer ( kind = 4 ), parameter :: l = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) a2(m,n+1)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) b2(m)
  real ( kind = 4 ) c(l,n+1)
  real ( kind = 4 ) d(l)
  real ( kind = 4 ) e(1,n+1)
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  real ( kind = 4 ) f(1)
  real ( kind = 4 ) g(1)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) v
  real ( kind = 4 ) x(n+1)
  real ( kind = 4 ) x_old(n+1)

  iflag = 0
  it = 0

  x(1:n) = 0.0E+00
  x(n+1) = 1.0E+00

  do

    x_old(1:n+1) = x(1:n+1)

    if ( itmax <= it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ORTH_LI - Fatal error!'
      write ( *, '(a)' ) '  Number of iterations exceeded.'
      iflag = 2
      return
    end if

    it = it + 1

    a2(1:m,1:n) = a(1:m,1:n)
    a2(1:m,n+1) = -b(1:m)
    b2(1:m) = 0.0E+00

    c(1,1:n+1) = x(1:n+1)
    d(1) = 1.0E+00

    call con_li ( m, l, k, n+1, a2, b2, c, d, e, f, g, eps1, x, res_norm, &
     iflag2 )

    if ( iflag2 /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ORTH_LI - Fatal error!'
      write ( *, '(a,i6)' ) '  CON_LI returned nonzero IFLAG2 = ', iflag2
      iflag = 1
      exit
    end if

    x(1:n+1) = x(1:n+1) / sqrt ( sum ( x(1:n+1)**2 ) )

    v = sum ( abs ( x(1:n+1) - x_old(1:n+1) ) )

    if ( v <= eps2 * sum ( abs ( x(1:n+1) ) ) ) then
      exit
    end if

  end do

  return
end
subroutine orth_lm ( z, m, n, ic, q0, q, sigma, iflag )

!*****************************************************************************80
!
!! ORTH_LM is a least squares solver for linear manifolds.
!
!  Discussion:
!
!    A set Z of M points in N+1 dimensions is given.
!
!    Vectors Q(0) through Q(N+1) in N+1 dimensions are to be determined
!    such that, for the linear manifolds
!
!      X = Q0 + T(1) * Q(1) + ... + T(N+1) * Q(N+1)
!
!    the sum of squared distances
!
!      sum ( 1 <= I <= M )
!        || Q(0) + sum ( 1 <= J <= N+1 ) T(J) * Q(J) - Z(I) ||**2
!
!    is minimized.  Here, the norm is the usual Euclidean or L2 norm.
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 296,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) Z(M,N+1), contains the values of the points
!    in its rows.
!
!    Input, integer ( kind = 4 ) M, the number of points.
!
!    Input, integer ( kind = 4 ) N, the number of vectors to compute is N+1.
!
!    Input, logical IC, is TRUE if Q0 is to be included, and FALSE
!    otherwise.
!
!    Output, real ( kind = 4 ) Q0(N+1), will contain Q0 if IC is TRUE.
!
!    Output, real ( kind = 4 ) Q(M,N+1), will contain the Q vectors as columns.
!
!    Output, real ( kind = 4 ) SIGMA(N+1), the singular values of Z.
!
!    Output, ingeger IFLAG, error flag.
!    0, no error detected.
!    K, the K-th singular value could not be detected after 30 iterations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical ic
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) k
  real ( kind = 4 ) q(m,n+1)
  real ( kind = 4 ) q0(n+1)
  real ( kind = 4 ) s
  real ( kind = 4 ) sigma(n+1)
  real ( kind = 4 ) z(m,n+1)

  if ( ic ) then

    do k = 1, n+1
      q0(k) = sum ( z(1:m,k) ) / real ( m, kind = 4 )
      z(1:m,k) = z(1:m,k) - q0(k)
    end do

  end if

  call svd ( m, n+1, z, sigma, .false., z, .true., q, iflag )

  return
end
subroutine orth_lp ( a, m, n, b, p, eps1, eps2, eps3, eps4, eps5, &
  itmax1, itmax2, x, iflag )

!*****************************************************************************80
!
!! ORTH_LP carries out orthogonal regression in the LP norm.
!
!  Discussion:
!
!    For the matrix Z = ( A | -b ), this routine tries to find a vector X
!    with unit L2 norm which minimizes the LP norm of Z*X.
!
!  Modified:
!
!    06 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 274-275,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) P, indicates the norm to be used.  P should be 
!    positive.
!
!    Input, real ( kind = 4 ) EPS1, a tolerance for ORTH_L2.
!
!    Input, real ( kind = 4 ) EPS2, a tolerance for ORTH_L2.
!
!    Input, real ( kind = 4 ) EPS3, a tolerance for the convergence of the 
!    iteration.
!
!    Input, real ( kind = 4 ) EPS4, ?
!
!    Input, real ( kind = 4 ) EPS5, ?
!
!    Input, integer ( kind = 4 ) ITMAX1, the maximum number of inner iterations.
!
!    Input, integer ( kind = 4 ) ITMAX2, the maximum number of outer iterations.
!
!    Output, real ( kind = 4 ) X(N+1), the computed solution.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error occurred.
!    1, the accuracy test failed in ORTH_L2.
!    2, there was no convergence in the inner iterations in ITMAX1 iterations.
!    3, P <= 1 is not allowed.
!    4, there was no convergence in the outer iterations in ITMAX2 iterations.
!    5, ?
!    6, the condition defined by EPS5 was violated.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) aw(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) bw(m)
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  real ( kind = 4 ) eps3
  real ( kind = 4 ) eps4
  real ( kind = 4 ) eps5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax1
  integer ( kind = 4 ) itmax2
  integer ( kind = 4 ) itmax3
  real ( kind = 4 ) p
  real ( kind = 4 ) s
  logical start
  real ( kind = 4 ) t
  real ( kind = 4 ) w(m)
  real ( kind = 4 ) x(n+1)
  real ( kind = 4 ) x_old(n+1)
  real ( kind = 4 ) x_start(n+1)
!
  it = 0
  iflag = 0

  if ( p <= 1.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ORTH_LP - Fatal error!'
    write ( *, '(a)' ) '  P <= 1 is not allowed.'
    iflag = 3
    return
  end if

  w(1:m) = 1.0E+00

  if ( p == 2.0E+00 ) then

    start = .false.

  else

    start = .true.

    x_start(1:n) = 0.0E+00
    x_start(n+1) = 1.0E+00

    x_old(1:n+1) = x_start(1:n+1)

  end if

  do

    if ( itmax2 <= it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ORTH_LP - Fatal error!'
      write ( *, '(a)' ) '  Number of iterations exceeded.'
      iflag = 4
      exit
    end if

    it = it + 1

    do i = 1, m
      aw(i,1:n) = w(i) * a(i,1:n)
    end do

    do i = 1, m
      bw(i) = w(i) * b(i)
    end do

    if ( 2.0E+00 < p ) then
      itmax3 = 1
    else
      itmax3 = itmax1
    end if

    call orth_l2 ( aw, m, n, bw, eps1, eps2, itmax3, x, start, x_start, iflag2 )

    if ( iflag2 /= 0 .and. iflag2 /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ORTH_LP - Fatal error!'
      write ( *, '(a)' ) '   ORTH_L2 returned nonzero IFLAG2 = ', iflag2
      iflag = 1
      exit
    end if

    if ( p == 2.0E+00 ) then
      exit
    end if

    if ( 2.0E+00 < p ) then

      t = 1.0E+00 / dot_product ( x(1:n+1), x_old(1:n+1) )

      x(1:n+1) = ( ( p - 2.0E+00 ) * x_old(1:n+1) + t * x(1:n+1) ) &
        / ( p - 1.0E+00 )

      x(1:n+1) = x(1:n+1) / sqrt ( sum ( x(1:n+1)**2 ) )

    end if

    s = sum ( abs ( x(1:n+1) ) )
    t = sum ( abs ( x(1:n+1) - x_old(1:n+1) ) )

    x_old(1:n+1) = x(1:n+1)
    x_start(1:n+1) = x(1:n+1)

    if ( t <= eps3 * s ) then
      exit
    end if

    do i = 1, m

      t = abs ( dot_product ( a(i,1:n), x(1:n) ) - b(i) * x(n+1) )

      if ( t < eps4 ) then
        iflag = 5
        return
      end if

      w(i) = t**( ( p - 2.0E+00 ) / 2.0E+00 )

      if ( w(i) <= eps5 ) then
        iflag = 6
        return
      end if

    end do

  end do

  return
end
subroutine perm_random2 ( n, p )

!*****************************************************************************80
!
!! PERM_RANDOM2 selects a random permutation of N objects.
!
!  Discussion:
!
!    The input values of P are used as labels; that is, the I-th object
!    is labeled P(I).
!
!  Modified:
!
!    12 May 2002
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), on input, a list of labels.
!    On output, the list has been permuted randomly.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) t

  do i = 1, n
    call i4_random ( i, n, j )
    t    = p(i)
    p(i) = p(j)
    p(j) = t
  end do

  return
end
function pythag ( a, b )

!*****************************************************************************80
!
!! PYTHAG computes SQRT ( A**2 + B**2 ) carefully.
!
!  Discussion:
!
!    The formula
!
!      PYTHAG = sqrt ( A**2 + B**2 )
!
!    is reasonably accurate, but can fail if, for example, A**2 is larger
!    than the machine overflow.  The formula can lose most of its accuracy
!    if the sum of the squares is very large or very small.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Modified:
!
!    02 March 2000
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the two legs of a right triangle.
!
!    Output, real ( kind = 4 ) PYTHAG, the length of the hypotenuse.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) p
  real ( kind = 4 ) pythag
  real ( kind = 4 ) r
  real ( kind = 4 ) s
  real ( kind = 4 ) t
  real ( kind = 4 ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0E+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

    do

      t = 4.0E+00 + r

      if ( t == 4.0E+00 ) then
        exit
      end if

      s = r / t
      u = 1.0E+00 + 2.0E+00 * s
      p = u * p
      r = ( s / u )**2 * r

    end do

  end if

  pythag = p

  return
end
subroutine qrbd ( q, e, nn, v, mdv, nrv, c, mdc, ncc, iflag )

!*****************************************************************************80
!
!! QRBD uses the QR algorithm for the singular values of a bidiagonal matrix.
!
!  Discussion:
!
!    The bidiagonal matrix D, with diagonal entries Q(1:N) and
!    superdiagonal entries E(2:N) is pre- and post-multiplied by
!    elementary rotation matrices R and P so that
!
!      Rk*...*R2*R1 * D * P1'*P2'*...*Pk' = diag(S1,...,Sn)
!
!  Reference:
!
!    Reinsch and Golub,
!    Singular Value Decomposition and Least Squares Solutions,
!    Numerische Mathematik,
!    Volume 14, 1970.
!
!    Charles Lawson and Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974,
!    Revised edition, SIAM, 1995.
!    QA275.L38
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) Q(NN).  On input, the diagonal
!    entries of the matrix.  On output, the singular values, which
!    are nonnegative, listed in nonincreasing order.
!
!    Input/output, real ( kind = 4 ) E(NN).  On input, E(2:NN) contains
!    the superdiagonal entries of the matrix.
!
!    Input, integer ( kind = 4 ) NN, the order of the matrix.
!
!    Input/output, real ( kind = 4 ) V(MDV,NN), an NRV by NN matrix
!    whose input value is to be postmultiplied by the transformations
!    P1' * P2' * ... * PK.
!
!    Input, integer ( kind = 4 ) MDV, the leading dimension of V.
!
!    Input, integer ( kind = 4 ) NRV, the number of rows in V.
!
!    Input/output, real ( kind = 4 ) C(MDC,NCC), an NN by NCC matrix
!    whose input value is to be premultiplied by the transformations
!    Rm * ... * R2 * R1.
!
!    Input, integer ( kind = 4 ) MDC, the leading dimension of C, which must be
!    at least NN.
!
!    Input, integer ( kind = 4 ) NCC, the number of columns in C.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error.
!    1, there was an error.
!
  implicit none

  integer ( kind = 4 ) mdc
  integer ( kind = 4 ) mdv
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nn

  real ( kind = 4 ) c(mdc,ncc)
  real ( kind = 4 ) cs
  real ( kind = 4 ) dnorm
  real ( kind = 4 ) e(nn)
  real ( kind = 4 ) f
  logical fail
  real ( kind = 4 ) g
  real ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nqrs
  integer ( kind = 4 ) nrv
  real ( kind = 4 ) q(nn)
  real ( kind = 4 ) r4_diff
  real ( kind = 4 ) small
  real ( kind = 4 ) sn
  real ( kind = 4 ) t
  real ( kind = 4 ) temp
  real ( kind = 4 ) v(mdv,nn)
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z
!
  n = nn
  iflag = 0

  if ( n <= 0 ) then
    return
  end if

  fail = .false.
  nqrs = 0
  e(1) = 0.0E+00

  dnorm = 0.0E+00
  do j = 1, n
    dnorm = max ( dnorm, abs ( q(j) ) + abs ( e(j) ) )
  end do

  do k = n, 1, -1
!
!  Test for splitting or rank deficiencies..
!  First test for last diagonal term, Q(K), being small.
!
20  continue

    if ( k == 1 ) then
      go to 50
    end if

    if ( r4_diff ( dnorm + q(k), dnorm ) /= 0.0E+00 ) then
      go to 50
    end if
!
!  Since Q(K) is small we make a special pass to transform E(K) to zero.
!
!  Transformation constructed to zero out position (I,K).
!
    cs = 0.0E+00
    sn = -1.0E+00

    do ii = 2, k

      i = k + 1 - ii
      f =     -sn * e(i+1)
      e(i+1) = cs * e(i+1)
      call g1 ( q(i), f, cs, sn, q(i) )
!
!  Accumulate the right hand transformations in V.
!
      do j = 1, nrv
        temp = v(j,i)
        v(j,i) =   cs * temp + sn * v(j,k)
        v(j,k)  = -sn * temp + cs * v(j,k)
      end do

    end do
!
!  The matrix is now bidiagonal, and of lower order
!  since E(K)  ==  zero.
!
50  continue

    do l = k, 1, -1

      if ( r4_diff ( dnorm+e(l), dnorm ) == 0.0E+00 ) then
        go to 100
      end if

      if ( r4_diff ( dnorm+q(l-1), dnorm ) == 0.0E+00 ) then
        go to 70
      end if

    end do
!
!  This loop can't complete since E(1) = zero.
!
    go to 100
!
!  Cancellation of E(L), 1 < L.
!
70  continue

    cs = 0.0E+00
    sn = -1.0E+00

    do i = l, k

      f =   -sn * e(i)
      e(i) = cs * e(i)

      if ( r4_diff ( dnorm+f, dnorm ) == 0.0E+00 ) then
        exit
      end if

      call g1 ( q(i), f, cs, sn, q(i) )
!
!  Accumulate the left hand transformations in C.
!
      do j = 1, ncc
        temp = c(i,j)
        c(i,j)   =   cs * temp + sn * c(l-1,j)
        c(l-1,j)  = -sn * temp + cs * c(l-1,j)
      end do

    end do
!
!  Test for convergence.
!
100 continue

    z = q(k)

    if ( l == k ) then
      go to 170
    end if
!
!  Shift from bottom 2 by 2 minor of B'*B.
!
    x = q(l)
    y = q(k-1)
    g = e(k-1)
    h = e(k)
    f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0E+00 * h * y )
    g = sqrt ( 1.0E+00 + f**2 )

    if ( 0.0E+00 <= f ) then
      t = f + g
    else
      t = f - g
    end if

    f = ( ( x - z ) * ( x + z ) + h * ( y / t - h ) ) / x
!
!  Next QR sweep.
!
    cs = 1.0E+00
    sn = 1.0E+00

    do i = l+1, k

      g = e(i)
      y = q(i)
      h = sn * g
      g = cs * g
      call g1 ( f, h, cs, sn, e(i-1) )
      f =  x * cs + g * sn
      g = -x * sn + g * cs
      h = y * sn
      y = y * cs
!
!  Accumulate rotations from the right in V.
!
      do j = 1, nrv
        temp      = v(j,i-1)
        v(j,i-1)  =  cs * temp + sn * v(j,i)
        v(j,i)    = -sn * temp + cs * v(j,i)
      end do

      call g1 ( f, h, cs, sn, q(i-1) )

      f =  cs * g + sn * y
      x = -sn * g + cs * y
!
!  Accumulate rotations from the left in C.
!
      do j = 1, ncc
        temp      = c(i-1,j)
        c(i-1,j)  =  cs * temp + sn * c(i,j)
        c(i,j)    = -sn * temp + cs * c(i,j)
      end do

    end do

    e(l) = 0.0E+00
    e(k) = f
    q(k) = x
    nqrs = nqrs + 1

    if ( nqrs <= 10 * n ) then
      go to 20
    end if
!
!  Return to 'test for splitting'.
!
    small = abs ( e(k) )
    i = k
!
!  If failure to converge set smallest magnitude
!  term in off-diagonal to zero.  Continue on.
!
    do j = l, k
      temp = abs ( e(j) )
      if ( temp /= 0.0E+00 ) then
        if ( temp < small ) then
          small = temp
          i = j
        end if
      end if
    end do

    e(i) = 0.0E+00
    nqrs = 0
    fail = .true.
    go to 20
!
!  Cutoff for convergence failure.  NQRS will be 2*N usually.
!
170 continue

    if ( z < 0.0E+00 ) then
      q(k) = -z
      v(1:nrv,k) = - v(1:nrv,k)
    end if

190 continue

  end do
!
!  Convergence. Q(K) is made nonnegative.
!
  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    if ( q(i-1) < q(i) ) then
      go to 220
    end if
  end do

  if ( fail ) then
    iflag = 1
  end if

  return
!
!  Every singular value is in order.
!
220 continue

  do i = 2, n

    t = q(i-1)
    k = i - 1

    do j = i, n
      if ( t < q(j) ) then
        t = q(j)
        k = j
      end if
    end do

    if ( k == i - 1 ) then
      cycle
    end if

    q(k) = q(i-1)
    q(i-1) = t

    do j = 1, ncc
      call r4_swap ( c(i-1,j), c(k,j) )
    end do

    do j = 1, nrv
      call r4_swap ( v(j,i-1), v(j,k) )
    end do

  end do
!
!  End of ordering algorithm.
!
  if ( fail ) then
    iflag = 1
  end if

  return
end
function r4_diff ( x, y )

!*****************************************************************************80
!
!! R4_DIFF computes the difference ( X - Y ) of two real numbers.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, Y, the numbers whose difference is desired.
!
!    Output, real ( kind = 4 ) R4_DIFF, the difference of the two numbers.
!
  implicit none

  real ( kind = 4 ) r4_diff
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  r4_diff = x - y

  return
end
subroutine r4_next ( s, r, done )

!*****************************************************************************80
!
!! R4_NEXT "reads" real numbers from a string, one at a time.
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
!    Input, character ( len = * ) S, a string, presumably containing real
!    numbers.  These may be separated by spaces or commas.
!
!    Output, real ( kind = 4 ) R.  If DONE is FALSE, then R contains the
!    "next" real value read from the string.  If DONE is TRUE, then
!    R is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh string, the user should set DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another real
!    value was read, or TRUE if no more reals could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ), save :: next = 1
  real ( kind = 4 ) r
  character ( len = * ) s

  r = 0.0E+00

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( len ( s ) < next ) then
    done = .true.
    return
  end if

  call s_to_r4 ( s(next:), r, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine r4_random ( rlo, rhi, r )

!*****************************************************************************80
!
!! R4_RANDOM returns a random real in a given range.
!
!  Modified:
!
!    06 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) RLO, RHI, the minimum and maximum values.
!
!    Output, real ( kind = 4 ) R, the randomly chosen value.
!
  implicit none

  real ( kind = 4 ) r
  real ( kind = 4 ) rhi
  real ( kind = 4 ) rlo
  real ( kind = 4 ) t
!
!  Pick T, a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R in ( RLO, RHI ).
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r4_swap ( x, y )

!*****************************************************************************80
!
!! R4_SWAP swaps two R4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real    ( kind = 4 ) x
  real    ( kind = 4 ) y
  real    ( kind = 4 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r4mat_cholesky_factor ( n, a, c )

!*****************************************************************************80
!
!! R4MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is a lower triangular matrix L such that:
!
!      A = L * L'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix A.
!
!    Input, real ( kind = 4 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 4 ) C(N,N), the N by N lower triangular
!    Cholesky factor.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n,n)
  real    ( kind = 4 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 4 ) sum2

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0E+00

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then
        if ( sum2 <= 0.0E+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R4MAT_CHOLESKY_FACTOR - Fatal error!'
          write ( *, '(a)' ) '  Matrix is not positive definite.'
          stop
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0E+00 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0E+00
        end if
      end if

    end do

  end do

  return
end
subroutine r4mat_cholesky_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R4MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix A.
!
!    Input, real ( kind = 4 ) A(N,N), the N by N Cholesky factor of the
!    system matrix.
!
!    Input, real ( kind = 4 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 4 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n,n)
  real    ( kind = 4 ) b(n)
  real    ( kind = 4 ) x(n)
!
!  Solve L * y = b.
!
  call r4mat_l_solve ( n, a, b, x )
!
!  Solve L' * x = y.
!
  call r4mat_lt_solve ( n, a, x, x )

  return
end
subroutine r4mat_diag_add_scalar ( n, a, s )

!*****************************************************************************80
!
!! R4MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns.
!
!    Input/output, real ( kind = 4 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 4 ) S, the value to be added to the diagonal
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 4 ) s

  do i = 1, n
    a(i,i) = a(i,i) + s
  end do

  return
end
subroutine r4mat_indicator ( m, n, table )

!*****************************************************************************80
!
!! R4MAT_INDICATOR sets up an "indicator" R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!    The value of each entry suggests its location, as in:
!
!      11  12  13  14
!      21  22  23  24
!      31  32  33  34
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, real ( kind = 4 ) TABLE(M,N), the table.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  real    ( kind = 4 ) table(m,n)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do j = 1, n
      table(i,j) = real ( fac * i + j, kind = 4 )
    end do
  end do

  return
end
subroutine r4mat_l_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R4MAT_L_SOLVE solves a lower triangular linear system.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix A.
!
!    Input, real ( kind = 4 ) A(N,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 4 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 4 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n,n)
  real    ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  real    ( kind = 4 ) x(n)
!
!  Solve L * x = b.
!
  do i = 1, n
    x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
  end do

  return
end
subroutine r4mat_lt_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R4MAT_LT_SOLVE solves a transposed lower triangular linear system.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!    Given the lower triangular matrix A, the linear system to be solved is:
!
!      A' * x = b
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 4 ) A(N,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 4 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 4 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n,n)
  real    ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  real    ( kind = 4 ) x(n)
!
!  Solve L'*x = b.
!
  do i = n, 1, -1
    x(i) = ( b(i) - dot_product ( x(i+1:n), a(i+1:n,i) ) ) / a(i,i)
  end do

  return
end
subroutine r4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R4MAT_PRINT prints an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2008
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
!    Input, real ( kind = 4 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 4 ) a(m,n)
  character ( len = * )  title

  call r4mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R4MAT_PRINT_SOME prints some of an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 4 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 4 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r4vec_bin ( n, a, nbin, bin_min, bin_max, bin, bin_limit )

!*****************************************************************************80
!
!! R4VEC_BIN bins a real vector, returning the population of each bin.
!
!  Discussion:
!
!    The user specifies minimum and maximum bin values, BIN_MIN and
!    BIN_MAX, and the number of bins, NBIN.  This determines a
!    "bin width":
!
!      H = ( BIN_MAX - BIN_MIN ) / NBIN
!
!    so that bin I will count all entries X(J) such that
!
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!
!    The array does NOT have to be sorted.
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input, real ( kind = 4 ) A(N), an (unsorted) array to be binned.
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.  Two extra bins, #0
!    and #NBIN+1, count extreme values.
!
!    Input, real ( kind = 4 ) BIN_MIN, BIN_MAX, define the range and size of 
!    the bins.  BIN_MIN and BIN_MAX must be distinct.
!    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
!    this, but proper results will be computed if BIN_MIN > BIN_MAX.
!
!    Output, integer ( kind = 4 ) BIN(0:NBIN+1).
!    BIN(0) counts entries of A less than BIN_MIN.
!    BIN(NBIN+1) counts entries greater than or equal to BIN_MAX.
!    For 1 <= I <= NBIN, BIN(I) counts the entries X(J) such that
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!    where H is the bin spacing.
!
!    Output, real ( kind = 4 ) BIN_LIMIT(0:NBIN), the "limits" of the bins.
!    BIN(I) counts the number of entries X(J) such that
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin

  real ( kind = 4 ) a(n)
  integer ( kind = 4 ) bin(0:nbin+1)
  real ( kind = 4 ) bin_limit(0:nbin)
  real ( kind = 4 ) bin_max
  real ( kind = 4 ) bin_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) t

  if ( bin_max == bin_min ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4VEC_BIN - Fatal error!'
    write ( *, '(a)' ) '  BIN_MIN = BIN_MAX.'
    stop
  end if

  bin(0:nbin+1) = 0

  do i = 1, n

    t = ( a(i) - bin_min ) / ( bin_max - bin_min )

    if ( t < 0.0E+00 ) then
      j = 0
    else if ( 1.0E+00 <= t ) then
      j = nbin + 1
    else
      j = 1 + int ( real ( nbin, kind = 4 ) * t )
    end if

    bin(j) = bin(j) + 1

  end do
!
!  Compute the bin limits.
!
  do i = 0, nbin
    bin_limit(i) = ( real ( nbin - i, kind = 4 ) * bin_min &
      + real ( i, kind = 4 ) * bin_max ) / real ( nbin, kind = 4 )
  end do

  return
end
subroutine r4vec_indicator ( n, a )

!*****************************************************************************80
!
!! R4VEC_INDICATOR sets an R4VEC to the indicator vector.
!
!  Discussion:
!
!    An R4VEC is a vector of R4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 4 )
  end do

  return
end
function r4vec_norm_l1 ( n, a )

!*****************************************************************************80
!
!! R4VEC_NORM_L1 returns the L1 norm of an R4VEC.
!
!  Discussion:
!
!    An R4VEC is a vector of R4's.
!
!    The vector L1 norm is defined as:
!
!      R4VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 4 ) A(N), the vector whose L1 norm is desired.
!
!    Output, real ( kind = 4 ) R4VEC_NORM_L1, the L1 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n)
  real    ( kind = 4 ) r4vec_norm_l1

  r4vec_norm_l1 = sum ( abs ( a(1:n) ) )

  return
end
function r4vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! R4VEC_NORM_L2 returns the L2 norm of an R4VEC.
!
!  Discussion:
!
!    An R4VEC is a vector of R4's.
!
!    The vector L2 norm is defined as:
!
!      R4VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 4 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 4 ) R4VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n)
  real    ( kind = 4 ) r4vec_norm_l2

  r4vec_norm_l2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
function r4vec_norm_li ( n, a )

!*****************************************************************************80
!
!! R4VEC_NORM_LI returns the L-oo norm of an R4VEC.
!
!  Discussion:
!
!    An R4VEC is a vector of R4's.
!
!    The vector L-oo norm is defined as:
!
!      R4VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 4 ) A(N), the vector whose L-oo norm is desired.
!
!    Output, real ( kind = 4 ) R4VEC_NORM_LI, the L-oo norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n)
  real    ( kind = 4 ) r4vec_norm_li

  r4vec_norm_li = maxval ( abs ( a(1:n) ) )

  return
end
function r4vec_norm_lp ( n, a, p )

!*****************************************************************************80
!
!! R4VEC_NORM_LP returns the LP norm of an R4VEC.
!
!  Discussion:
!
!    An R4VEC is a vector of R4's.
!
!    The vector LP norm is defined as:
!
!      R4VEC_NORM_LP = ( sum ( 1 <= I <= N ) ( abs ( A(I) ) )**P )**(1/P).
!
!    Usually, the LP norms with
!      1 <= P <= oo
!    are of interest.  This routine allows
!      0 < P <= Huge ( P ).
!    If P = Huge ( P ), then the L-oo norm is returned, which is
!    simply the maximum of the absolute values of the vector components.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 4 ) A(N), the vector whose LP norm is desired.
!
!    Input, real ( kind = 4 ) P, the index of the norm.
!
!    Output, real ( kind = 4 ) R4VEC_NORM_LP, the LP norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) a(n)
  real    ( kind = 4 ) p
  real    ( kind = 4 ) r4vec_norm_lp

  if ( p <= 0.0E+00 ) then
    r4vec_norm_lp = -1.0E+00
  else if ( p == huge ( p ) ) then
    r4vec_norm_lp = maxval ( abs ( a(1:n) ) )
  else if ( p == 1.0E+00 ) then
    r4vec_norm_lp = sum ( abs ( a(1:n) ) )
  else if ( p == 2.0E+00 ) then
    r4vec_norm_lp = sqrt ( sum ( a(1:n)**2 ) )
  else
    r4vec_norm_lp = ( sum ( ( abs ( a(1:n) ) )**p ) )**( 1.0E+00 / p )
  end if

  return
end
subroutine r4vec_print ( n, a, title )

!*****************************************************************************80
!
!! R4VEC_PRINT prints an R4VEC.
!
!  Discussion:
!
!    An R4VEC is a vector of R4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  real      ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine random_partition ( m, n, z )

!*****************************************************************************80
!
!! RANDOM_PARTITION generates a random partition.
!
!  Discussion:
!
!    The partition assigns each of M objects to one of N sets.
!    Every object belongs to exactly one set, but some sets
!    may be empty, or have several elements.
!
!  Modified:
!
!    26 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 153,
!    ISBN 0-12-656460-4.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of objects.
!
!    Input, integer ( kind = 4 ) N, the number of sets.
!
!    Output, integer ( kind = 4 ) Z(M), contains, for each object, the index of
!    the set, between 1 and N, containing the object.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 4 ) urand
  integer ( kind = 4 ) z(m)

  do i = 1, m

    z(i) = int ( real ( n, kind = 4 ) * urand ( ) ) + 1

  end do

  return
end
subroutine random_partition2 ( m, n, ml, mu, z )

!*****************************************************************************80
!
!! RANDOM_PARTITION2 generates a random partition with occupancy constraints.
!
!  Discussion:
!
!    The partition assigns each of M objects to one of N sets.
!    Every object belongs to exactly one set.  Every set will contain
!    at least ML objects, and no more than MU objects.
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 153,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of objects.
!
!    Input, integer ( kind = 4 ) N, the number of sets.
!
!    Output, integer ( kind = 4 ) Z(M), contains, for each object, the index of
!    the set, between 1 and N, containing the object.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  real ( kind = 4 ) urand
  integer ( kind = 4 ) z(m)

  if ( ml < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_PARTITION2 - Fatal error!'
    write ( *, '(a)' ) '  ML < 0.'
    write ( *, '(a,i6)' ) '  ML = ', ml
    stop
  end if
!
!  Take care of the minimum occupancy requirement.
!
  if ( m < n * ml ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_PARTITION2 - Fatal error!'
    write ( *, '(a)' ) '  M < ML * N.'
    write ( *, '(a,i6)' ) '  N  = ', n
    write ( *, '(a,i6)' ) '  ML = ', ml
    write ( *, '(a,i6)' ) '  M  = ', m
    stop
  end if

  k = 0

  do i = 1, n
    do j = 1, ml
      k = k + 1
      z(k) = i
    end do
  end do

  if ( mu * n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_PARTITION2 - Fatal error!'
    write ( *, '(a)' ) '  MU is too small.'
    write ( *, '(a)' ) '  MU * N < M.'
    write ( *, '(a,i6)' ) '  N  = ', n
    write ( *, '(a,i6)' ) '  MU = ', mu
    write ( *, '(a,i6)' ) '  M  = ', m
    stop
  end if

  if ( ml < mu ) then
!
!  We need to find a random partition of M - ML * N into N parts,
!  with the constraint that the largest part is no greater than MU - ML.
!
    do

      call npart_rsf_lex_random ( m-ml*n, n, a )

      if ( maxval ( a(1:n) ) <= mu - ml ) then
        exit
      end if

    end do
!
!  Assign the remaining entries of Z according to the partition.
!
    do i = 1, n
      do j = 1, a(i)
        k = k + 1
        z(k) = i
      end do
    end do

  end if
!
!  Now randomly permute the vector.
!
  call perm_random2 ( m, z )

  return
end
subroutine regr_lp ( a, m, n, b, p, eps1, eps2, eps3, itmax, x, iflag )

!*****************************************************************************80
!
!! REGR_LP minimizes the LP norm of A*x-b for P > 1.
!
!  Discussion:
!
!    Least squares minimization corresponds to P = 2.
!
!  Modified:
!
!    25 April 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 51,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(M,N).  On input, A contains the system
!    matrix.  On output, A contains information about the factorization of 
!    the matrix, which can be used to efficiently solve more systems.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input/output, real ( kind = 4 ) B(M), on input, the right hand side vector.
!    On output, B has been overwritten by other information.
!
!    Input, real ( kind = 4 ) P, determines the LP norm to be used.  That is,
!    || X || = ( sum ( 1 <= I <= N ) X(I)**P )**(1/P)
!    For this routine, P must be greater than 1.
!
!    Input, real ( kind = 4 ) EPS1, the tolerance for the accuracy test in MGS_L2.
!
!    Input, real ( kind = 4 ) EPS2, the tolerance for convergence of the iteration.
!
!    Input, real ( kind = 4 ) EPS3.  If during the iterations, one component of
!    A * X - B becomes less than EPS3 in modulus, then the iteration
!    is stopped with return IFLAG = 4.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations.
!
!    Output, real ( kind = 4 ) X(N), the approximate solution.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error was observed.
!    1, there was an error in MGS_L2.
!    2, the input value of P was less than or equal to 1.
!    3, ITMAX iterations were taken without convergence.
!    4, see remarks under EPS3.  X contains an approximate solution.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  logical bnew
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  real ( kind = 4 ) eps3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) k
  real ( kind = 4 ) p
  real ( kind = 4 ) r(n,n)
  real ( kind = 4 ) s
  real ( kind = 4 ) t
  real ( kind = 4 ) w(m)
  real ( kind = 4 ) wa(m,n)
  real ( kind = 4 ) wb(m)
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  iflag = 0

  if ( p <= 1.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REGR_LP - Fatal error!'
    write ( *, '(a)' ) '  Input value of P less than or equal to 1.'
    iflag = 2
    return
  end if

  it = 0

  w(1:m) = 1.0E+00
  y(1:n) = 0.0E+00
!
!  Next iteration.
!
  do

    if ( itmax <= it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REGR_LP - Warning!'
      write ( *, '(a)' ) '  No convergence after ITMAX iterations.'
      iflag = 3
      return
    end if

    it = it + 1

    wb(1:m) = w(1:m) * b(1:m)

    do i = 1, m
      wa(i,1:n) = w(i) * a(i,1:n)
    end do

    bnew = .false.

    call mgs_l2 ( wa, m, n, wb, eps1, bnew, x, iflag2 )

    if ( iflag2 /= 0 ) then
      iflag = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REGR_LP - Warning!'
      write ( *, '(a)' ) '  MGS_L2 returned an error flag.'
      exit
    end if
!
!  If the user requested P = 2, then that's exactly what MGS_L2 returned,
!  so we're done.
!
    if ( p == 2.0E+00 ) then
      iflag = 0
      exit
    end if

    if ( 2.0E+00 <= p ) then
      x(1:n) = ( ( p - 2.0E+00 ) * y(1:n) + x(1:n) ) / ( p - 1.0E+00 )
    end if
!
!  Compute the solution norm and stepsize.
!
    s = sum ( abs ( x(1:n) ) )

    t = sum ( abs ( x(1:n) - y(1:n) ) )

    y(1:n) = x(1:n)
!
!  If the stepsize was relatively large, prepare for another iteration.
!
    if ( t < eps2 * s ) then
      iflag = 0
      exit
    end if

    do i = 1, m

      s = abs ( dot_product ( a(i,1:n), x(1:n) ) - b(i) )

      if ( s < eps3 ) then
        iflag = 4
      else
        w(i) = s**( ( p - 2.0E+00 ) / 2.0E+00 )
      end if

    end do

    if ( iflag /= 0 ) then
      exit
    end if

  end do

  return
end
subroutine residual ( a, m, n, b, x, nbin, r, rmin, rmean, rmax, r_bin_count, &
  v, vmin, vmean, vmax, v_bin_count )

!*****************************************************************************80
!
!! RESIDUAL calculates the residual vector A*X-B and related information.
!
!  Discussion:
!
!    The routine calculates the residual vector:
!
!      R(1:M)  =  A(1:M,1:N) * X(1:N) - B(1:M)
!
!    and its minimum, maximum and mean.  It also divides the interval
!    ( -|RMIN|, RMAX) into NBIN intervals, and returns the number of
!    entries of R that lie in each interval.
!
!    Similarly, the relative residual vector is computed:
!
!      V(I) = 100 * R(I) / ( A(I,1:N) * X(1:N) )
!
!    and corresponding data determined.
!
!  Modified:
!
!    12 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 109-110,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N) contains the M by N system matrix.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Output, real ( kind = 4 ) X(N), the calculated solution.
!
!    Input, integer ( kind = 4 ) NBIN, the number of intervals to use.
!
!    Output, real ( kind = 4 ) R(M), the residuals.
!
!    Output, real ( kind = 4 ) RMIN, the minimum residual value.
!
!    Output, real ( kind = 4 ) RMEAN, the mean residual value.
!
!    Output, real ( kind = 4 ) RMAX, the maximum residual value.
!
!    Output, integer ( kind = 4 ) R_BIN_COUNT(0:NBIN+1), the number of residual
!    vector entries in each of the NBIN intervals from -|RMIN| to RMAX.
!    The 0-th and NBIN+1 entries count extreme values.
!
!    Output, real ( kind = 4 ) V(M), the relative residuals.
!
!    Output, real ( kind = 4 ) VMIN, the minimum relative residual.
!
!    Output, real ( kind = 4 ) VMEAN, the mean relative residual.
!
!    Output, real ( kind = 4 ) VMAX, the maximum relative residual.
!
!    Output, integer ( kind = 4 ) V_BIN_COUNT(0:NBIN+1), the number of residual
!    vector entries in each of the NBIN intervals from -|VMIN| to VMAX.
!    The 0-th and NBIN+1 entries count extreme values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) ax
  real ( kind = 4 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mv
  real ( kind = 4 ) r(m)
  integer ( kind = 4 ) r_bin_count(0:nbin+1)
  real ( kind = 4 ) r_bin_limit(0:nbin)
  real ( kind = 4 ) rmax
  real ( kind = 4 ) rmean
  real ( kind = 4 ) rmin
  real ( kind = 4 ) v(m)
  integer ( kind = 4 ) v_bin_count(0:nbin+1)
  real ( kind = 4 ) v_bin_limit(0:nbin)
  real ( kind = 4 ) vmax
  real ( kind = 4 ) vmean
  real ( kind = 4 ) vmin
  real ( kind = 4 ) x(n)

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
  rmin = minval ( r(1:m) )
  rmax = maxval ( r(1:m) )
  rmean = sum ( r(1:m) ) / real ( m, kind = 4 )

  call r4vec_bin ( m, r, nbin, rmin, rmax, r_bin_count, r_bin_limit )

  mv = 0
  vmean = 0.0E+00
  vmin =   huge ( vmin )
  vmax = - huge ( vmax )

  do i = 1, m

    ax = dot_product ( a(i,1:n), x(1:n) )

    if ( ax /= 0.0E+00 ) then

      v(i) = ( r(i) / ax ) * 100.0E+00
      vmin = min ( vmin, v(i) )
      vmax = max ( vmax, v(i) )
      vmean = vmean + v(i)
      mv = mv + 1

    else if ( r(i) == 0.0E+00 ) then

      v(i) = 0.0E+00
      vmin = min ( vmin, v(i) )
      vmax = max ( vmax, v(i) )
      vmean = vmean + v(i)
      mv = mv + 1

    else

      v(i) = huge ( v(i) )

    end if

  end do

  vmean = vmean / real ( mv, kind = 4 )

  call r4vec_bin ( m, v, nbin, vmin, vmax, v_bin_count, v_bin_limit )

  return
end
subroutine robust ( a, m, n, b, method, eps1, eps2, eps3, itmax, &
  x, sr, iflag )

!*****************************************************************************80
!
!! ROBUST carries out robust regression, with eight choices for the method.
!
!  Modified:
!
!    02 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, pages 196-197,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N and greater than 3.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input, integer ( kind = 4 ) METHOD, the method to be used.
!    1, Andrews
!    2, Biweight
!    3, Cauchy
!    4, Fair
!    5, Huber
!    6, Logistic
!    7, Talwar
!    8, Welsch
!
!    Input, real ( kind = 4 ) EPS1, the tolerance for the accuracy test 
!    in MGS_L2.
!
!    Input, real ( kind = 4 ) EPS2, the tolerance for convergence of
!    the iteration.
!
!    Input, real ( kind = 4 ) EPS3, tolerance for an accuracy test.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations.
!
!    Output, real ( kind = 4 ) X(N), the approximate solution.
!
!    Output, real ( kind = 4 ) SR, the sum of
!    W(1:M) * ( 0.6745 * R(1:M) / ALPHA )**2.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error was observed.
!    1, there was an error in MGS_L2.
!    2, ITMAX iterations were taken without convergence.
!    3, Illegal or inconsistent values of M, N or METHOD.
!    4, see remarks under EPS3.  X contains an approximate solution.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) ar(m)
  real ( kind = 4 ) armed
  real ( kind = 4 ) b(m)
  logical bnew
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  real ( kind = 4 ) eps3
  logical final
  logical finit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) method
  integer ( kind = 4 ) mr
  real ( kind = 4 ) res(m)
  real ( kind = 4 ) s
  real ( kind = 4 ) sr
  real ( kind = 4 ) t
  real ( kind = 4 ) u
  real ( kind = 4 ) ua
  real ( kind = 4 ) w(m)
  real ( kind = 4 ) wa(m,n)
  real ( kind = 4 ) wb(m)
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)
!
  iflag = 0

  if ( m <= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROBUST - Fatal error!'
    write ( *, '(a)' ) '  M <= N.'
    iflag = 3
    return
  end if

  if ( m <= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROBUST - Fatal error!'
    write ( *, '(a)' ) '  M <= 3.'
    iflag = 3
    return
  end if

  if ( method < 1 .or. 8 < method ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROBUST - Fatal error!'
    write ( *, '(a)' ) '  METHOD < 1 or 8 < METHOD.'
    iflag = 3
    return
  end if

  it = 0
  sr = 0.0E+00
  finit = .false.

  w(1:m) = 1.0E+00
  y(1:n) = 0.0E+00

  do

    if ( itmax <= it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROBUST - Fatal error!'
      write ( *, '(a)' ) '  Number of iterations exceeded.'
      iflag = 2
      return
    end if

    it = it + 1

    wb(1:m) = w(1:m) * b(1:m)
    do i = 1, m
      wa(i,1:n) = w(i) * a(i,1:n)
    end do

    bnew = .false.
    call mgs_l2 ( wa, m, n, wb, eps1, bnew, x, iflag2 )

    if ( iflag2 /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROBUST - Fatal error!'
      write ( *, '(a)' ) '  MGS_L2 returned nonzero IFLAG = ', iflag2
      iflag = 1
      return
    end if

    t = sum ( abs ( x(1:n) - y(1:n) ) )
    s = sum ( abs ( x(1:n) ) )

    y(1:n) = x(1:n)

    if ( t < eps2 * s ) then
      finit = .true.
    end if

    mr = 0

    res(1:m) = b(1:m) - matmul ( a(1:m,1:n), x(1:n) )

    if ( it == 1 ) then

      ar(1:m) = res(1:m)

      do j = 2, m

        final = .true.

        do k = 1, m - j + 1

          if ( ar(k+1) < ar(k) ) then
            call r4_swap ( ar(k), ar(k+1) )
            final = .false.
          end if

        end do

        if ( final ) then
          exit
        end if

      end do

      m1 = m / 2
      m2 = (m+1) / 2

      if ( mod ( m, 2 ) == 0 ) then
        armed = 0.5E+00 * ( ar(m2) + ar(m2+1) )
      else
        armed = ar(m2)
      end if

      ar(1:m) = abs ( res(1:m) - armed )

      do j = 2, m

        final = .true.

        do k = 1, m - j + 1

          if ( ar(k+1) < ar(k) ) then
            call r4_swap ( ar(k), ar(k+1) )
            final = .false.
          end if

        end do

        if ( final ) then
          exit
        end if

      end do

      if ( mod ( m, 2 ) == 0 ) then
        armed = 0.5E+00 * ( ar(m2) + ar(m2+1) )
      else
        armed = ar(m2)
      end if

      if ( armed <= eps3 ) then
        iflag = 4
        return
      end if

      armed = 0.6745E+00 / armed

    end if

    do i = 1, m

      u = res(i) * armed
      w(i) = 0.0E+00

      if ( method == 1 ) then

        if ( 4.207E+00 < abs ( u ) ) then

        else if ( eps3 <= abs ( u ) ) then

          ua = u / 1.339E+00
          w(i) = sin ( ua ) / ua

        else if ( abs ( u ) < eps3 ) then

          w(i) = 1.0E+00

        end if

      else if ( method == 2 ) then

        if ( abs ( u ) <= 4.685E+00 ) then
          w(i) = ( 1.0E+00 - ( u / 4.685E+00 )**2 )**2
        end if

      else if ( method == 3 ) then

        w(i) = 1.0E+00 / ( 1.0E+00 + ( u / 2.385E+00 )**2 )

      else if ( method == 4 ) then

        w(i) = 1.0E+00 / ( 1.0E+00 + abs ( u ) / 1.4E+00 )

      else if ( method == 5 ) then

        if ( 1.345E+00 < abs ( u ) ) then
          w(i) = 1.345E+00 / abs ( u )
        else
          w(i) = 1.0E+00
        end if

      else if ( method == 6 ) then

        ua = u / 1.205E+00

        if ( ua < eps3 ) then
          w(i) = 1.0E+00
        else
          w(i) = tanh ( ua ) / ua
        end if

      else if ( method == 7 ) then

        if ( abs ( u ) <= 2.795E+00 ) then
          w(i) = 1.0E+00
        end if

      else if ( method == 8 ) then

        ua = ( u / 2.985E+00 )**2

        if ( 30.0E+00 < ua ) then
          w(i) = 0.0E+00
        else
          w(i) = exp ( -ua )
        end if

      end if

      w(i) = sqrt ( w(i) )

    end do

    if ( finit ) then
      exit
    end if

  end do

  sr = armed**2 * dot_product ( w(1:m)**2, res(1:m)**2 )

  return
end
subroutine rr_l1 ( a, m, n, b, eps, x, lambda, iflag )

!*****************************************************************************80
!
!! RR_L1 carries out ridge regression in the L1 norm.
!
!  Discussion:
!
!    The routine calls A478_L1.
!
!    The routine minimizes the L1 norm of the objective function:
!
!    ( A          )  * X - ( B )
!    ( Lambda * I )        ( 0 )
!
!    for user specified LAMBDA.
!
!  Modified:
!
!    15 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 211,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/, real A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, real ( kind = 4 ) EPS, a tolerance used in the rank determination.
!
!    Output, real ( kind = 4 ) X(N), the solution vector.
!
!    Input, real ( kind = 4 ) LAMBDA, the value of LAMBDA in the objective function.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, an optimal solution was found.  It is probably not unique.
!    1, an optimal solution was found, and it is unique.
!    2, the calculations were prematurely stopped because of rounding
!       errors.
!    3, a fatal error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) a2(m+n,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) b2(m+n)
  real ( kind = 4 ) e2(m+n)
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  real ( kind = 4 ) lambda
  integer ( kind = 4 ) rank
  real ( kind = 4 ) x(n)

  a2(1:m,1:n) = a(1:m,1:n)

  a2(m+1:m+n,1:n) = 0.0E+00
  do i = 1, n
    a2(m+i,i) = lambda
  end do

  b2(1:m) = b(1:m)
  b2(m+1:m+n) = 0.0E+00

  call a478_l1 ( a2, m+n, n, b2, eps, rank, x, e2, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RR_L1 - Fatal error!'
    write ( *, '(a)' ) '  A478_L1 returned nonzero IFLAG = ', iflag
  end if

  return
end
subroutine rr_l2 ( a, m, n, b, eps, x, lambda, iflag )

!*****************************************************************************80
!
!! RR_L2 carries out ridge regression in the L2 norm.
!
!  Discussion:
!
!    The routine calls MGS_L2.
!
!    The routine minimizes the L2 norm of the objective function:
!
!    ( A          )  * X - ( B )
!    ( Lambda * I )        ( 0 )
!
!    for user specified LAMBDA.
!
!  Modified:
!
!    03 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 208,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(M,N).  On input, A contains the system matrix.
!    On output, A contains information about the factorization of the matrix,
!    which can be used to efficiently solve more systems.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input/output, real ( kind = 4 ) B(M), on input, the right hand side vector.
!    On output, B has been overwritten by other information.
!
!    Input, real ( kind = 4 ) EPS, accuracy tolerance for the orthogonalization.  If
!    the squared euclidean length of the vector to be treated in the
!    K-th step is smaller than EPS, then the routine is halted with
!    IFLAG = 1.  Recommended value: EPS =  10**(-2*T+4).
!
!    Output, real ( kind = 4 ) X(N), the calculated solution.
!
!    Input, real ( kind = 4 ) LAMBDA, the value of LAMBDA in the objective function.
!
!    Output, integer ( kind = 4 ) IFLAG.
!    0, no error occurred, the calculation was completed.
!    1, the calculation has halted because the system seemed singular.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) a2(m+n,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) b2(m+n)
  logical bnew
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 4 ) lambda
  real ( kind = 4 ) x(n)
!
  a2(1:m,1:n) = a(1:m,1:n)
  b2(1:m) = b(1:m)

  a2(m+1:m+n,1:n) = 0.0E+00
  do i = 1, n
    a2(m+i,i) = lambda
  end do

  b2(m+1:m+n) = 0.0E+00

  bnew = .false.

  call mgs_l2 ( a2, m+n, n, b2, eps, bnew, x, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RR_L2 - Fatal error!'
    write ( *, '(a)' ) '  MGS_L2 returned nonzero IFLAG = ', iflag
  end if

  return
end
subroutine rr_li ( a, m, n, b, eps1, eps2, x, lambda, iflag )

!*****************************************************************************80
!
!! RR_LI carries out ridge regression in the L-infinity norm.
!
!  Discussion:
!
!    The routine calls A495_LI.
!
!    The routine minimizes the L-infinity norm of the objective function:
!
!    ( A          )  * X - ( B )
!    ( Lambda * I )        ( 0 )
!
!    for user specified LAMBDA.
!
!  Modified:
!
!    15 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 214,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N).  The system matrix.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input/output, real ( kind = 4 ) B(M), on input, the right hand side vector.
!    On output, B has been overwritten by other information.
!
!    Input, real ( kind = 4 ) EPS1, a tolerance for an accuracy test.
!
!    Input, real ( kind = 4 ) EPS2, a tolerance for the minimization condition.
!    Setting EPS2 to 0.0 means that a minimum is sought.  Setting
!    0 < EPS2 allows for a relative error between the norm of the residual
!    that is found and the norm of the residual for the minimal
!    solution.
!
!    Output, real ( kind = 4 ) X(N), the calculated solution.
!
!    Input, real ( kind = 4 ) LAMBDA, the value of LAMBDA in the objective function.
!
!    Output, integer ( kind = 4 ) IFLAG.
!    0, no error occurred, the calculation was completed.
!    1, the calculation has halted because the system seemed singular.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) a2(m+n,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) b2(m+n)
  real ( kind = 4 ) eps1
  real ( kind = 4 ) eps2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 4 ) lambda
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) x(n)

  a2(1:m,1:n) = a(1:m,1:n)

  a2(m+1:m+n,1:n) = 0.0E+00

  do i = 1, n
    a2(m+i,i) = lambda
  end do

  b2(1:m) = b(1:m)
  b2(m+1:m+n) = 0.0E+00

  call a495_li ( a2, m+n, n, b2, eps1, eps2, x, rank, res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RR_LI - Fatal error!'
    write ( *, '(a)' ) '  A495_LI returned nonzero IFLAG = ', iflag
  end if

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
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
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
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
subroutine s_to_r4 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R4 reads a real number from a string.
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
!    Output, real ( kind = 4 ) R, the real value that was read from the string.
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
  real ( kind = 4 ) r
  real ( kind = 4 ) rbot
  real ( kind = 4 ) rexp
  real ( kind = 4 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0E+00
  lchar = - 1
  isgn = 1
  rtop = 0.0E+00
  rbot = 1.0E+00
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
        rtop = 10.0E+00 * rtop + real ( ndig, kind = 4 )
      else if ( ihave == 5 ) then
        rtop = 10.0E+00 * rtop + real ( ndig, kind = 4 )
        rbot = 10.0E+00 * rbot
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
    rexp = 1.0E+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0E+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0E+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine scr ( a, m, n, nl, nu, first, na, aa, kbit )

!*****************************************************************************80
!
!! SCR selects M by NA submatrices from an M by N matrix.
!
!  Bogosity:
!
!    The user input NL was an INPUT/OUTPUT quantity, which was incremented
!    up to NU by the program as it uses larger and larger subsets.  This
!    is an unacceptable violation of the intuitive interface.  Even after
!    I identified this poor programming, I was trapped by it, passing in
!    a value of NL = 2 as a parameter, which got incremented to 4, and
!    then polluted other values which had been set to 2!
!
!    DON'T TRY TO BE "EFFICIENT" BY USING A NATURAL INPUT QUANTITY AS
!    A SHORTCUT TO PASSING BACK OUTPUT INFORMATION, OR TO BE USED AS
!    SAVED MEMORY!
!
!  Discussion:
!
!    SCR is able to list, one at a time, all the M by NA submatrices of an
!    M by N matrix A, where NA is allowed to range between NL and NU.
!
!    Each time SCR is called, a new possiblity is returned.
!
!    The user signals the first call by setting FIRST to TRUE, and
!    the KBIT array to 0.
!
!    The routine warns the user that there are no more submatrices
!    by setting the return value of FIRST to FALSE.
!
!  Modified:
!
!    13 May 2002
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 126,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the matrix from which submatrices are selected.
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, integer ( kind = 4 ) NL, NU, the limits on NA.  It is required that
!    1 <= NL <= NA <= NU <= N.
!
!    Input/output, logical FIRST.  On the first call to SCR, set FIRST
!    to TRUE.  The output value of FIRST will be set to FALSE until the last
!    possible submatrix has been determined.  On the next call, the routine
!    will return with FIRST set to TRUE.
!
!    Output, integer ( kind = 4 ) NA, the number of columns that were selected to
!    construct the current submatrix.
!
!    Output, real ( kind = 4 ) AA(M,NU), the selected M by NA submatrix.
!
!    Input/output, KBIT(N), indicates which columns of A are used in
!    the submatrix.  KBIT(I) = 1 means column I is in the submatrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nu

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) aa(m,nu)
  logical first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kbit(n)
  integer ( kind = 4 ) ksum
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nl
  integer ( kind = 4 ), save :: nl_save = 0
!
!  On the first call for a given problem, make a copy of NL that you
!  can modify.
!
  if ( first ) then
    nl_save = nl
    first = .false.
  end if

  do

    k = 1
!
!  Modify the current KBIT configuration.
!
!  Starting with the first index, set the initial string of
!  contiguous 1's to 0.
!
!  Exit when you reach the end of the vector, or you see
!  the first 0 entry, which should be reset to 1.
!
    do

      if ( kbit(k) == 0 ) then
        kbit(k) = 1
        exit
      end if

      kbit(k) = 0
      k = k + 1

      if ( n < k ) then
        exit
      end if

    end do
!
!  Add up the bits.
!
    ksum = sum ( kbit(1:n) )

    if ( ksum == 0 ) then

      nl_save = nl_save + 1

      if ( nu < nl_save ) then
        first = .true.
        exit
      end if

    else if ( ksum /= nl_save ) then

      cycle

    else if ( ksum == nl_save ) then

      na = 0

      do k = 1, n

        if ( kbit(k) /= 0 ) then
          na = na + 1
          aa(1:m,na) = a(1:m,k)
        end if

      end do

      exit

    end if

  end do

  return
end
subroutine scrf_l1 ( a, m, n, b, na, eps, x, kbit, res_norm, iflag )

!*****************************************************************************80
!
!! SCRF_L1 minimizes the L1 norm of A*X-B using NA variables out of N.
!
!  Discussion:
!
!    For a problem with N variables, the routine checks all regression
!    possibilities with a smaller, fixed number of variables NA,
!    and selects the first combination with minimum sum of absolute
!    deviations.
!
!  Modified:
!
!    26 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 140,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(M,N), the system matrix.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix A.
!
!    Input, real ( kind = 4 ) B(M), the right hand side.
!
!    Input, integer ( kind = 4 ) NA, the number of variables to be selected.
!    1 <= NA <= N.
!
!    Input, real ( kind = 4 ) EPS, a tolerance used in the rank determination.
!
!    Output, real ( kind = 4 ) X(NA), the computed solution using only the
!    selected variables.
!
!    Output, integer ( kind = 4 ) KBIT(N), contains a 1 for those variables which
!    are included, and a 0 for those which are excluded.
!
!    Output, real ( kind = 4 ) RES_NORM, the residual.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error occurred.
!    1, an error occurred in A478_L1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) aa(m,na)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) bb(m)
  real ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kbit(n)
  integer ( kind = 4 ) kbit_best(n)
  integer ( kind = 4 ) nb
  real ( kind = 4 ) r(m)
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) res_norm_best
  real ( kind = 4 ) x(na)
  real ( kind = 4 ) x_best(na)

  if ( na < 1 .or. n < na ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCRF_L1 - Fatal error!'
    write ( *, '(a)' ) '  NA is out of bounds.'
    write ( *, '(a)' ) '  We require 1 <= NA <= N.'
    write ( *, '(a,i6)' ) '  NA = ', na
    write ( *, '(a,i6)' ) '  N  = ' , n
    stop
  end if

  res_norm_best = huge ( res_norm_best )

  kbit(1:n) = 0

  do
!
!  Select the variables.
!
    call c01m ( n, n-na, kbit )
!
!  Assemble the submatrix.
!
    nb = 0

    do k = 1, n

      if ( kbit(k) == 1 ) then
        nb = nb + 1
        aa(1:m,nb) = a(1:m,k)
      end if

    end do

    bb(1:m) = b(1:m)
!
!  Solve the subproblem.
!
    call a478_l1 ( aa, m, na, bb, eps, rank, x, r, iflag2 )

    if ( iflag2 /= 0 ) then
      iflag = 1
      return
    end if

    res_norm = sum ( abs ( r(1:m) ) )

    if ( res_norm < res_norm_best ) then

      res_norm_best = res_norm
      x_best(1:na) = x(1:na)
      kbit_best(1:n) = kbit(1:n)

    end if

    if ( sum ( kbit(1:na) ) == na ) then
      exit
    end if

  end do

  res_norm = res_norm_best
  x(1:na) = x_best(1:na)
  kbit(1:n) = kbit_best(1:n)

  return
end
subroutine svd ( m, n, a, s, matu, u, matv, v, iflag )

!*****************************************************************************80
!
!! SVD computes the singular value decomposition for a real matrix.
!
!  Discussion:
!
!    This subroutine determines the singular value decomposition
!
!      A = U * S * V'
!
!    of a real M by N rectangular matrix.
!
!    U will be an M by N matrix with orthogonal columns;
!    S will be an N by N diagonal matrix with nonnegative entries;
!    V will be an N by N orthogonal matrix.
!
!    Householder bidiagonalization and a variant of the QR algorithm are used.
!
!  Modified:
!
!    22 April 2002
!
!  Reference:
!
!    Golub and Reinsch,
!    Numerische Mathematik,
!    Volume 14, 1970, pages 403-420.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A and U.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A and U,
!    and the order of V.
!
!    Input, real ( kind = 4 ) A(M,N), the M by N matrix to be decomposed.
!
!    Output, real ( kind = 4 ) S(N), the singular values of A.  These are the
!    diagonal elements of S.  They are unordered.  If an error exit is
!    made, the singular values should be correct for indices
!    IERR+1, IERR+2,..., N.
!
!    Input, logical MATU, should be set to TRUE if the U matrix in the
!    decomposition is desired, and to FALSE otherwise.
!
!    Output, real ( kind = 4 ) U(M,N), contains the matrix U, with orthogonal columns,
!    of the decomposition, if MATU has been set to TRUE.  Otherwise
!    U is used as a temporary array.  U may coincide with A.
!    If an error exit is made, the columns of U corresponding
!    to indices of correct singular values should be correct.
!
!    Input, logical MATV, should be set to TRUE if the V matrix in the
!    decomposition is desired, and to FALSE otherwise.
!
!    Output, real ( kind = 4 ) V(N,N), the orthogonal matrix V of the decomposition if
!    MATV has been set to TRUE.  Otherwise V is not referenced.
!    V may also coincide with A if U is not needed.  If an error
!    exit is made, the columns of V corresponding to indices of
!    correct singular values should be correct.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, for normal return,
!    K, if the K-th singular value has not been determined after 30 iterations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) c
  real ( kind = 4 ) f
  real ( kind = 4 ) g
  real ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: itmax = 30
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) l1
  logical matu
  logical matv
  integer ( kind = 4 ) mn
  real ( kind = 4 ) pythag
  real ( kind = 4 ) rv1(n)
  real ( kind = 4 ) s(n)
  real ( kind = 4 ) scale
  real ( kind = 4 ) ss
  real ( kind = 4 ) tst1
  real ( kind = 4 ) tst2
  real ( kind = 4 ) u(m,n)
  real ( kind = 4 ) v(n,n)
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  iflag = 0
  u(1:m,1:n) = a(1:m,1:n)
!
!  Householder reduction to bidiagonal form.
!
  g = 0.0E+00
  scale = 0.0E+00
  x = 0.0E+00

  do i = 1, n

    l = i + 1
    rv1(i) = scale * g
    g = 0.0E+00
    ss = 0.0E+00
    scale = 0.0E+00

    if ( i <= m ) then

      scale = sum ( abs ( u(i:m,i) ) )

      if ( scale /= 0.0E+00 ) then

        u(i:m,i) = u(i:m,i) / scale

        ss = sum ( u(i:m,i)**2 )

        f = u(i,i)
        g = - sign ( sqrt ( ss ), f )
        h = f * g - ss
        u(i,i) = f - g

        if ( i /= n ) then

          do j = l, n
            ss = dot_product ( u(i:m,i), u(i:m,j) )
            u(i:m,j) = u(i:m,j) + ss * u(i:m,i) / h
          end do

        end if

        u(i:m,i) = scale * u(i:m,i)

      end if

    end if

    s(i) = scale * g
    g = 0.0E+00
    ss = 0.0E+00
    scale = 0.0E+00

    if ( i <= m .and. i /= n ) then

      scale = sum ( abs ( u(i,l:n) ) )

      if ( scale /= 0.0E+00 ) then

        u(i,l:n) = u(i,l:n) / scale
        ss = sum ( u(i,l:n)**2 )
        f = u(i,l)
        g = - sign ( sqrt ( ss ), f )
        h = f * g - ss
        u(i,l) = f - g
        rv1(l:n) = u(i,l:n) / h

        if ( i /= m ) then

          do j = l, m

            ss = dot_product ( u(j,l:n), u(i,l:n) )

            u(j,l:n) = u(j,l:n) + ss * rv1(l:n)

          end do

        end if

        u(i,l:n) = scale * u(i,l:n)

      end if

    end if

    x = max ( x, abs ( s(i) ) + abs ( rv1(i) ) )

  end do
!
!  Accumulation of right-hand transformations.
!
  if ( matv ) then

    do i = n, 1, -1

      if ( i /= n ) then

         if ( g /= 0.0E+00 ) then

          v(l:n,i) = ( u(i,l:n) / u(i,l) ) / g

          do j = l, n

            ss = dot_product ( u(i,l:n), v(l:n,j) )

            v(l:n,j) = v(l:n,j) + ss * v(l:n,i)

          end do

        end if

        v(i,l:n) = 0.0E+00
        v(l:n,i) = 0.0E+00

      end if

      v(i,i) = 1.0E+00
      g = rv1(i)
      l = i

    end do

  end if
!
!  Accumulation of left-hand transformations.
!
  if ( matu ) then

    mn = min ( m, n )

    do i = min ( m, n ), 1, -1

      l = i + 1
      g = s(i)

      if ( i /= n ) then
        u(i,l:n) = 0.0E+00
      end if

      if ( g /= 0.0E+00 ) then

        if ( i /= mn ) then

          do j = l, n
            ss = dot_product ( u(l:m,i), u(l:m,j) )
            f = ( ss / u(i,i) ) / g
            u(i:m,j) = u(i:m,j) + f * u(i:m,i)
          end do

        end if

        u(i:m,i) = u(i:m,i) / g

      else

        u(i:m,i) = 0.0E+00

      end if

      u(i,i) = u(i,i) + 1.0E+00

    end do

  end if
!
!  Diagonalization of the bidiagonal form.
!
  tst1 = x

  do kk = 1, n

     k1 = n - kk
     k = k1 + 1
     it = 0
!
!  Test for splitting.
!
10   continue

     do ll = 1, k

       l1 = k - ll
       l = l1 + 1
       tst2 = tst1 + abs ( rv1(l) )

       if ( tst2 == tst1 ) then
         go to 20
       end if

       tst2 = tst1 + abs ( s(l1) )

       if ( tst2 == tst1 ) then
         exit
       end if

     end do
!
!  Cancellation of rv1(l) if L greater than 1.
!
     c = 0.0E+00
     ss = 1.0E+00

     do i = l, k

       f = ss * rv1(i)
       rv1(i) = c * rv1(i)
       tst2 = tst1 + abs ( f )

       if ( tst2 == tst1 ) then
         go to 20
       end if

       g = s(i)
       h = pythag ( f, g )
       s(i) = h
       c = g / h
       ss = -f / h

       if ( matu ) then

         do j = 1, m
           y = u(j,l1)
           z = u(j,i)
           u(j,l1) = y * c + z * ss
           u(j,i) = -y * ss + z * c
         end do

       end if

    end do
!
!  Test for convergence.
!
20 continue

    z = s(k)

    if ( l == k ) then
      go to 30
    end if
!
!  Shift from bottom 2 by 2 minor.
!
    if ( itmax <= it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SVD - Fatal error!'
      write ( *, '(a)' ) '  Number of iterations exceeded.'
      iflag = k
      return
    end if

    it = it + 1
    x = s(l)
    y = s(k1)
    g = rv1(k1)
    h = rv1(k)
    f = 0.5E+00 * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y )
    g = pythag ( f, 1.0E+00 )
    f = x - ( z / x ) * z + ( h / x ) * ( y / ( f + sign ( g, f ) ) - h)
!
!  Next QR transformation.
!
    c = 1.0E+00
    ss = 1.0E+00

    do i1 = l, k1

      i = i1 + 1
      g = rv1(i)
      y = s(i)
      h = ss * g
      g = c * g
      z = pythag ( f, h )
      rv1(i1) = z
      c = f / z
      ss = h / z
      f = x * c + g * ss
      g = -x * ss + g * c
      h = y * ss
      y = y * c

      if ( matv ) then

        do j = 1, n
          x = v(j,i1)
          z = v(j,i)
          v(j,i1) = x * c + z * ss
          v(j,i) = -x * ss + z * c
        end do

      end if

      z = pythag ( f, h )
      s(i1) = z
!
!  Rotation can be arbitrary if Z is zero.
!
      if ( z /= 0.0E+00 ) then
        c = f / z
        ss = h / z
      end if

      f = c * g + ss * y
      x = -ss * g + c * y

      if ( matu ) then

        do j = 1, m
          y = u(j,i1)
          z = u(j,i)
          u(j,i1) = y * c + z * ss
          u(j,i) = -y * ss + z * c
        end do

      end if

    end do

    rv1(l) = 0.0E+00
    rv1(k) = f
    s(k) = x
    go to 10
!
!  Convergence.
!
30 continue

    if ( z <= 0.0E+00 ) then

      s(k) = - z

      if ( matv ) then
        v(1:n,k) = - v(1:n,k)
      end if

    end if

  end do

  return
end
subroutine svdr_l2 ( a, m, n, b, eps, x, cond, rank, iflag )

!*****************************************************************************80
!
!! SVDR_L2 minimizes the L2 norm of A*x-b using the SVD.
!
!  Discussion:
!
!    The singular value decomposition (SVD) of A is used to determine
!    a numerical value for the L2 condition number COND of the matrix A,
!    the estimated rank of A and the minimal norm solution X of the
!    least squares problem A * X = B.
!
!  Modified:
!
!    24 April 2002
!
!  Reference:
!
!    Golub and Van Loan,
!    Matrix Computations,
!    The Johns Hopkins University Press, 1983 page 140.
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 40,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(M,N).  On input, A contains the system matrix.
!    On output, A contains information about the factorization of the matrix,
!    which can be used to efficiently solve more systems.
!
!    Input, integer ( kind = 4 ) M, the first dimension of A, the number of
!    observations or equations.  M must be greater than N.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 4 ) B(M), the right hand side vector.
!
!    Input, real ( kind = 4 ) EPS, used for an accuracy test.  Singular values less
!    than or equal to EPS will be regarded as zero.
!
!    Output, real ( kind = 4 ) X(N), the least squares minimizing solution.
!
!    Output, real ( kind = 4 ) COND, the ratio of the largest to smallest singular
!    value if the matrix A has maximal rank N, or HUGE(COND) if one
!    or more singular values are zero.
!
!    Output, integer ( kind = 4 ) RANK, the numerically determined rank of A.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    0, no error occurred.
!    nonzero, the IFLAG-th singular value could not be determined after
!      30 iterations.  The computation was not completed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) b(m)
  real ( kind = 4 ) cond
  real ( kind = 4 ) eps
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) k
  logical, parameter :: matu = .true.
  logical, parameter :: matv = .true.
  integer ( kind = 4 ) rank
  real ( kind = 4 ) sigma(n)
  real ( kind = 4 ) smax
  real ( kind = 4 ) smin
  real ( kind = 4 ) u(m,n)
  real ( kind = 4 ) v(n,n)
  real ( kind = 4 ) x(n)
!
!  Get the singular value decomposition.
!
  call svd ( m, n, a, sigma, matu, u, matv, v, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDR_L2 - Warning!'
    write ( *, '(a,i6)' ) '  SVD returned with nonzero IFLAG = ', iflag
    return
  end if
!
!  Determine the largest singular value.
!
  smax = maxval ( sigma(1:n) )
  smin = minval ( sigma(1:n) )

  if ( smin < eps ) then
    smin = 0.0E+00
  end if

  x(1:n) = 0.0E+00

  rank = 0

  do k = 1, n

    if ( eps * smax < sigma(k) ) then

      rank = rank + 1

      x(1:n) = x(1:n) + v(1:n,k) * dot_product ( b(1:m), u(1:m,k) ) / sigma(k)

    end if

  end do

  if ( rank < n ) then
    cond = huge ( cond )
  else
    cond = smax / smin
  end if

  return
end
subroutine svdrs ( a, mda, m, n, b, mdb, nb, s )

!*****************************************************************************80
!
!! SVDRS: singular value decomposition with a right side vector.
!
!  Discussion:
!
!    This subroutine computes the singular value decomposition of the
!    M by N matrix A and optionally applies the transformations
!    from the left to the NB column vectors of the M by NB matrix B.
!
!    Either N <= M or M < N is permitted.
!
!    The singular value decomposition of A is of the form
!
!      A  =  U * S * V'
!
!    where U is M by M orthogonal, S is M by N diagonal with the
!    diagonal terms nonnegative and ordered from large to small, and
!    V is N by N orthogonal.  Note that these matrices also satisfy
!
!      S = U' * A * V
!
!    The matrix V is returned in the leading N rows and
!    columns of the array A.
!
!    The singular values, i.e. the diagonal terms of the matrix S,
!    are returned in the array S.  If M < N, positions M+1
!    through N of S will be set to zero.
!
!    The product matrix G = U' * B overwrites the input matrix B.
!
!    If the user wishes to obtain a minimum length least squares
!    solution of the linear system
!
!      A * x ~ = ~ b
!
!    the solution X can be constructed, following use of this subroutine,
!    by computing the sum for I = 1, ..., R of the outer products
!
!      (column I of V) * 1/S(I) * (row I of G)
!
!    Here R denotes the pseudorank of A which the user may choose
!    in the range 0 through min ( M, N ) based on the sizes of the
!    singular values.
!
!    This code gives special treatment to rows and columns that are
!    entirely zero.  This causes certain zero singular values to appear as
!    exact zeros rather than as about macheps times the largest singular value.
!    It similarly cleans up the associated columns of U and V.
!
!  Algorithm:
!
!    1. Exchange columns of A to pack nonzero columns to the left.
!       Set N_COPY = number of nonzero columns.
!       Use locations A(1,N), A(1,N-1),..., A(1,N_COPY+1) to record the
!       column permutations.
!
!    2. Exchange rows of A to pack nonzero rows to the top.
!       Quit packing if find N_COPY nonzero rows.  Make same row exchanges
!       in B.  Set M_COPY so that all nonzero rows of the permuted A
!       are in first M_COPY rows.  If M_COPY <= N_COPY then all M_COPY rows are
!       nonzero.  If M_COPY > N_COPY then the first N_COPY rows are known
!       to be nonzero, and rows NCOPY+1 through M_COPY may be zero or nonzero.
!
!    3. Apply the original algorithm to the M_COPY by N_COPY problem.
!
!    4. Move permutation record from A(,) to S(I), I = N_COPY+1,...,N.
!
!    5. Build V up from N_COPY by N_COPY to N by N by placing ones on
!       the diagonal and zeros elsewhere.  This is only partly done
!       explicitly.  It is completed during step 6.
!
!    6. Exchange rows of V to compensate for column exchanges of step 2.
!
!    7. Place zeros in S(N_COPY+1:N)  to represent zero singular values.
!
!  Reference:
!
!    Charles Lawson and Richard Hanson,
!    Solving Least Squares Problems,
!    Prentice-Hall, 1974,
!    Revised edition, SIAM, 1995.
!    QA275.L38
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(MDA,N).
!    On input, the M by N matrix A.
!    On output contains the N by N matrix V.
!
!    Input, integer ( kind = 4 ) MDA, the leading dimension of A.
!    MDA must be at least the maximum of M and N.
!
!    Input, integer ( kind = 4 ) M, the number of rows of A and B.
!    M must be at least 1.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A,
!    the number of rows and columns of the matrix V.
!    N must be greater than 0.
!
!    Input/output, real ( kind = 4 ) B(MDB,NB).  If 0 < NB, this array
!    must contain an M by NB matrix on input and will contain the
!    M by NB product matrix, G = U' * B on output.
!
!    Input, integer ( kind = 4 ) MDB, the leading dimension of B, which must be
!    at least M.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides; that is, the
!    number of columns of data in B.  If no right hand sides are being
!    supplied, set NB = 0.
!
!    Output, real ( kind = 4 ) S(N), the singular values of A, with the ordering
!    S(1) >= S(2) >= ... >= S(N) >= 0.
!    If M < N the singular values indexed from M+1 through N will be zero.
!
!  Local parameters:
!
!    real WORK(N,2).
!    locations 1 thru N will hold the off-diagonal terms of
!    the bidiagonal matrix for subroutine QRBD.  Locations N+1
!    thru 2*N will save info from one call to the next of H12.
!
  implicit none

  integer ( kind = 4 ) mda
  integer ( kind = 4 ) mdb
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 4 ) a(mda,n)
  real ( kind = 4 ) b(mdb,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_copy
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) n_copy
  integer ( kind = 4 ) ns
  integer ( kind = 4 ) nsp1
  real ( kind = 4 ) s(n)
  real ( kind = 4 ) t
  real ( kind = 4 ) work(n,2)

  n_copy = n

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDRS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of M = ', m
    stop
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDRS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of N = ', n
    stop
  end if

  if ( mda < m .or. mda < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDRS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of MDA = ', mda
    stop
  end if

  if ( mdb < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDRS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of MDB = ', mdb
    stop
  end if
!
!  If column J is entirely zero, exchange it with column N.
!
  j = n_copy

  do while ( 1 <= j )

    if ( all ( a(1:m,j) == 0.0E+00 ) ) then

      if ( j /= n_copy ) then
        a(1:m,j) = a(1:m,n_copy)
      end if

      a(1,n_copy) = j
      n_copy = n_copy - 1

    end if

    j = j - 1

  end do

  ns = 0
!
!  If N_COPY = 0 then A is entirely zero and the SVD computation can be skipped.
!
  if ( n_copy /= 0 ) then
!
!  Pack nonzero rows to the top.
!  Quit packing if find N_COPY nonzero rows.
!
    i = 1
    m_copy = m

    do while ( i <= n_copy .and. i < m_copy )

      if ( a(i,i) == 0.0E+00 ) then

        if ( all ( a(i,1:n_copy) == 0.0E+00 ) ) then

          do j = 1, nb
            call r4_swap ( b(i,j), b(m_copy,j) )
          end do

          a(i,1:n_copy) = a(m_copy,1:n_copy)

          if ( m_copy <= n_copy ) then
            a(m_copy,1:n_copy) = 0.0E+00
          end if

           m_copy = m_copy - 1

        else
          i = i + 1
        end if

      else
        i = i + 1
      end if

    end do
!
!  Begin the SVD algorithm:
!
!  1: Reduce the matrix to upper bidiagonal form with Householder
!     transformations
!       H(N_COPY)...H(1) * A * Q(1)...Q(N_COPY-2) = (D',0)'
!     where D is upper bidiagonal.
!
!  2: Apply H(N_COPY)...H(1) to B.  Here H(N_COPY)...H(1)*B replaces
!     B in storage.
!
!  3: The matrix product W = Q(1)...Q(N_COPY-2) overwrites the first
!     N_COPY rows of A in storage.
!
!  4: An SVD for D is computed.  Here K rotations Ri and Pi are
!     computed so that
!       Rk...R1*D*P1'...Pk' = diag(S1,...,Sm_copy)
!     to working accuracy.  The S(I) are nonnegative and nonincreasing.
!     Here Rk...R1*B overwrites B in storage while
!     A*P1'...Pk' overwrites A in storage.
!
!  5: It follows that, with the proper definitions, U'*B overwrites B,
!     while V overwrites the first N_COPY row and columns of A.
!
    l = min ( m_copy, n_copy )
!
!  The following loop reduces A to upper bidiagonal and
!  also applies the premultiplying transformations to B.
!
    do j = 1, l
!
!  Certain HARMLESS illegal memory references will occur in the
!  following lines.
!
      if ( j < m_copy ) then

        mode = 1
        if ( 0 < n_copy-j ) then
          call h12 ( mode, j, j+1, m_copy, a(1,j), 1, t, a(1,j+1), 1, mda, &
            n_copy-j )
        else
          call h12 ( mode, j, j+1, m_copy, a(1,j), 1, t, a(1,1), 1, mda, &
            n_copy-j )
        end if

        mode = 2
        call h12 ( mode, j, j+1, m_copy, a(1,j), 1, t, b, 1, mdb, nb )

      end if

      if ( j < n_copy-1 ) then

        mode = 1
        call h12 ( mode, j+1, j+2, n_copy, a(j,1), mda, work(j,2), a(j+1,1), &
          mda, 1, m_copy-j )

      end if

    end do
!
!  Copy the bidiagonal matrix into S and WORK for QRBD.
!
    do j = 2, l
      s(j) = a(j,j)
      work(j,1) = a(j-1,j)
    end do
    s(1) = a(1,1)

    ns = n_copy

    if ( m_copy < n_copy ) then
      ns = m_copy + 1
      s(ns) = 0.0E+00
      work(ns,1) = a(m_copy,m_copy+1)
    end if
!
!  Construct the explicit N_COPY by N_COPY product matrix,
!  W = Q1*Q2*...*QL*I in A.
!
    do i = n_copy, 1, -1

      if ( i <= min ( m_copy, n_copy-2 ) ) then
        mode = 2
        call h12 ( mode, i+1, i+2, n_copy, a(i,1), mda, work(i,2), a(1,i+1), &
          1, mda, n_copy-i )
      end if

      a(i,1:n_copy) = 0.0E+00
      a(i,i) = 1.0E+00

    end do
!
!  Compute the singular value decomposition of the bidiagonal matrix.
!
    call qrbd ( s(1), work(1,1), ns, a, mda, n_copy, b, mdb, nb, iflag2 )

    if ( iflag2 == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SVDRS - Warning:'
      write ( *, '(a)' ) '  Full accuracy was not attained in the singular'
      write ( *, '(a)' ) '  value decomposition of the bidiagonal matrix.'
    end if

  end if

  s(ns+1:n_copy) = 0.0E+00
!
!  Extract the record of permutations and store zeros.
!
  s(n_copy+1:n) = a(1,n_copy+1:n)
  a(1:n_copy,n_copy+1:n) = 0.0E+00
!
!  Permute rows and set zero singular values.
!
  do k = n_copy+1, n

    i = int ( s(k) )

    s(k) = 0.0E+00

    a(k,1:n) = a(i,1:n)

    a(i,1:n) = 0.0E+00
    a(i,k) = 1.0E+00

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
function uniform_01_sample ( iseed )

!*****************************************************************************80
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!  Formula:
!
!    ISEED = ISEED * (7**5) mod (2**31 - 1)
!    RANDOM = ISEED * / ( 2**31 - 1 )
!
!  Modified:
!
!    01 March 1999
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ISEED, the integer "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  ISEED should not be zero.
!
!    Output, real ( kind = 4 ) UNIFORM_01_SAMPLE, a random value between 0 and 1.
!
!  Local Parameters:
!
!    IA = 7**5
!    IB = 2**15
!    IB16 = 2**16
!    IP = 2**31-1
!
  implicit none

  integer ( kind = 4 ), parameter :: ia = 16807
  integer ( kind = 4 ), parameter :: ib15 = 32768
  integer ( kind = 4 ), parameter :: ib16 = 65536
  integer ( kind = 4 ), parameter :: ip = 2147483647
  integer ( kind = 4 ) iprhi
  integer ( kind = 4 ) iseed
  integer ( kind = 4 ) ixhi
  integer ( kind = 4 ) k
  integer ( kind = 4 ) leftlo
  integer ( kind = 4 ) loxa
  real ( kind = 4 ) uniform_01_sample
!
!  Don't let ISEED be 0.
!
  if ( iseed == 0 ) then
    iseed = ip
  end if
!
!  Get the 15 high order bits of ISEED.
!
  ixhi = iseed / ib16
!
!  Get the 16 low bits of ISEED and form the low product.
!
  loxa = ( iseed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are
!  essential.
!
  iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( iseed < 0 ) then
    iseed = iseed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = real ( iseed, kind = 4 ) * 4.656612875E-10

  return
end
function urand ( )

!*****************************************************************************80
!
!! URAND returns a uniformly distributed pseudo random number.
!
!  Modified:
!
!    28 March 2000
!
!  Reference:
!
!    Helmuth Spaeth,
!    Mathematical Algorithms for Linear Regression,
!    Academic Press, 1991, page 153,
!    ISBN 0-12-656460-4.
!
!  Parameters:
!
!    Output, real ( kind = 4 ) URAND, a uniformly distributed pseudo random number.
!
  implicit none

  integer ( kind = 4 ), save :: ix = 15021
  integer ( kind = 4 ), save :: iy = 23456
  integer ( kind = 4 ), save :: iz = 30156
  real ( kind = 4 ) urand

  ix = 171 * mod ( ix, 177 ) -  2 * ( ix / 177 )
  iy = 172 * mod ( iy, 176 ) - 35 * ( iy / 176 )
  iz = 170 * mod ( iz, 178 ) - 63 * ( iz / 178 )

  if ( ix < 0 ) then
    ix = ix + 30269
  end if

  if ( iy < 0 ) then
    iy = iy + 30307
  end  if

  if ( iz < 0 ) then
    iz = iz + 30323
  end if

  urand = mod ( &
      real ( ix, kind = 4 ) / 30269.0E+00 &
    + real ( iy, kind = 4 ) / 30307.0E+00 &
    + real ( iz, kind = 4 ) / 30323.0E+00, 1.0E+00 )

  return
end
