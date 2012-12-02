subroutine cl1 ( k, l, m, n, klmd, klm2d, nklmd, n2d, q, kode, toler, iter, &
x, res, error, cu, iu, s )

!*****************************************************************************80
!
!! CL1 solves a linear programming problem.
!
!  Discussion:
!
!    This routine uses a modification of the simplex
!    method of linear programming to calculate an l1 solution
!    to a k by n system of linear equations
!      ax=b
!    subject to l linear equality constraints
!      cx=d
!    and m linear inequality constraints
!      ex.le.f.
!
!    if your fortran compiler permits a single column of a two
!    dimensional array to be passed to a one dimensional array
!    through a subroutine call, considerable savings in
!    execution time may be achieved through the use of 
!    subroutine col which operates on column vectors.
!
!    see comments following statement labelled 440 for
!    instructions on the implementation of this modification.
!
!  Modified:
!
!    27 December 2007
!
!  Author:
!
!    Ian Barrodale,
!    Frank Roberts.
!
!  Reference:
!
!    Ian Barrodale, Frank Roberts,
!    Algorithm 552:
!    Solution of the Constrained L1 Approximation Problem,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 2, March 1980, pages 231-235.
!
!  Parameters:
!
!    Input, integer K, the number of rows of the matrix A.
!    1 <= K.
!
!    Input, integer L, the number of rows of the matrix C.
!    0 <= L.
!
!    Input, integer M, the number of rows of the matrix E.
!    0 <= M.
!
!    n      number of columns of the matrices a,c,e (n.ge.1).
!
!    klmd   set to at least k+l+m for adjustable dimensions.
!
!    klm2d  set to at least k+l+m+2 for adjustable dimensions.
!
!    nklmd  set to at least n+k+l+m for adjustable dimensions.
!
!   n2d    set to at least n+2 for adjustable dimensions
!
!    q      two dimensional real array with klm2d rows and
!           at least n2d columns.
!           on entry the matrices a,c and e, and the vectors
!           b,d and f must be stored in the first k+l+m rows
!           and n+1 columns of q as follows
!                a b
!            q = c d
!                e f
!           these values are destroyed by the subroutine.
!
!    kode   a code used on entry to, and exit
!           from, the subroutine.
!           on entry, this should normally be set to 0.
!           however, if certain nonnegativity constraints
!           are to be included implicitly, rather than
!           explicitly in the constraints ex.le.f, then kode
!           should be set to 1, and the nonnegativity
!           constraints included in the arrays x and
!           res (see below).
!           on exit, kode has one of the
!           following values
!             0- optimal solution found,
!             1- no feasible solution to the
!                constraints,
!             2- calculations terminated
!                prematurely due to rounding errors,
!             3- maximum number of iterations reached.
!
!    toler  a small positive tolerance. empirical
!           evidence suggests toler = 10**(-d*2/3),
!           where d represents the number of decimal
!           digits of accuracy available. essentially,
!           the subroutine cannot distinguish between zero
!           and any quantity whose magnitude does not exceed
!           toler. in particular, it will not pivot on any
!           number whose magnitude does not exceed toler.
!
!    iter   on entry iter must contain an upper bound on
!           the maximum number of iterations allowed.
!           a suggested value is 10*(k+l+m). on exit iter
!           gives the number of simplex iterations.
!
!    x      one dimensional real array of size at least n2d.
!           on exit this array contains a
!           solution to the l1 problem. if kode=1
!           on entry, this array is also used to include
!           simple nonnegativity constraints on the
!          variables. the values -1, 0, or 1
!           for x(j) indicate that the j-th variable
!           is restricted to be .le.0, unrestricted,
!           or .ge.0 respectively.
!
!    res    one dimensional real array of size at least klmd.
!           on exit this contains the residuals b-ax
!           in the first k components, d-cx in the
!           next l components (these will be =0),and
!           f-ex in the next m components. if kode=1 on
!           entry, this array is also used to include simple
!           nonnegativity constraints on the residuals
!           b-ax. the values -1, 0, or 1 for res(i)
!           indicate that the i-th residual (1.le.i.le.k) is
!           restricted to be .le.0, unrestricted, or .ge.0
!           respectively.
!
!    error  on exit, this gives the minimum sum of
!           absolute values of the residuals.
!
!    cu     a two dimensional real array with two rows and
!           at least nklmd columns used for workspace.
!
!    iu     a two dimensional integer array with two rows and
!           at least nklmd columns used for workspace.
!
!    s      integer array of size at least klmd, used for
!           workspace.
!
  implicit none

  integer klm2d
  integer klmd
  integer n2d

  real cu
  real cuv
  real error
  integer i
  integer ia
  integer ii
  integer iimn
  integer in
  integer iu
  integer j
  integer jmn
  integer jpn
  integer js
  integer k
  integer kk
  integer klm
  integer l
  integer m
  integer n
  integer n1
  integer n2
  integer nk
  integer nk1
  integer nkl
  real pivot
  real q
  real res
  integer s
  real sn
  double precision sum
  real toler
  real tpivot
  real x
  real xmax
  real xmin
  real z
  real zu
  real zv

      integer iout, iter, klm1, klm2, kode, nklm, nkl1
      integer maxit, nklmd, iphase, kforce, iineg
      dimension q(klm2d,n2d), x(n2d), res(klmd)
      dimension cu(2,nklmd), iu(2,nklmd), s(klmd)
!
!  Initialization.
!
  maxit = iter
  n1 = n + 1
  n2 = n + 2
  nk = n + k
  nk1 = nk + 1
  nkl = nk + l
  nkl1 = nkl + 1
  klm = k + l + m
  klm1 = klm + 1
  klm2 = klm + 2
  nklm = n + klm
  kforce = 1
  iter = 0
  js = 1
  ia = 0
!
!  Set up labels in Q.
!
  do j = 1, n
    q(klm2,j) = j
  end do

  do i = 1, klm
    q(i,n2) = n + i
    if ( q(i,n1) .lt. 0.0E+00 ) then
      q(i,1:n2) = - q(i,1:n2)
    end if
  end do
!
!  Set up phase 1 costs.
!
  iphase = 2

  cu(1:2,1:nklm) = 0.0E+00
  iu(1:2,1:nklm) = 0

  if ( l .ne. 0 ) then

    cu(1:2,nk1:nkl) = 1.0E+00
    iu(1:2,nk1:nkl) = 1
    iphase = 1

  end if

  if ( m .ne. 0 ) then

    cu(2,nkl1:nklm) = 1.0E+00
    iu(2,nkl1:nklm) = 1

    do j = nkl1, nklm
      if (q(j-n,n2) .lt. 0.0E+00 ) then
        iphase = 1
      end if
    end do

  end if

  if ( kode .ne. 0 ) then

      do j=1,n
         if (x(j)) 90, 110, 100
   90    cu(1,j) = 1.0
         iu(1,j) = 1
         go to 110
  100    cu(2,j) = 1.0
         iu(2,j) = 1
  110    continue
      end do

      do j=1,k
         jpn = j + n
         if (res(j)) 120, 140, 130
  120    cu(1,jpn) = 1.0
         iu(1,jpn) = 1
         if (q(j,n2).gt.0.0) iphase = 1
         go to 140
  130    cu(2,jpn) = 1.0
         iu(2,jpn) = 1
         if (q(j,n2).lt.0.0) iphase = 1
  140    continue
      end do

  end if

  if (iphase.eq.2) go to 500
!
!  Compute the marginal costs.
!
  160 continue

      do j=js,n1
         sum = 0.0D+00
         do i=1,klm
            ii = q(i,n2)
            if ( 0 .le. ii ) then
              z = cu(1,ii)
            else
              iineg = -ii
              z = cu(2,iineg)
            end if
            sum = sum + dble(q(i,j))*dble(z)
         end do
         q(klm1,j) = sum
      end do

      do j=js,n
         ii = q(klm2,j)
         if (ii.lt.0) go to 210
         z = cu(1,ii)
         go to 220
  210    iineg = -ii
         z = cu(2,iineg)
  220    q(klm1,j) = q(klm1,j) - z
      end do
!
!  Determine the vector to enter the basis.
!
  240 continue

      xmax = 0.0E+00
      if (js.gt.n) go to 490
      do j=js,n
         zu = q(klm1,j)
         ii = q(klm2,j)
         if (ii.gt.0) go to 250
         ii = -ii
         zv = zu
         zu = -zu - cu(1,ii) - cu(2,ii)
         go to 260
  250    zv = -zu - cu(1,ii) - cu(2,ii)
  260    if (kforce.eq.1 .and. ii.gt.n) go to 280
         if (iu(1,ii).eq.1) go to 270
         if (zu.le.xmax) go to 270
         xmax = zu
         in = j
  270    if (iu(2,ii).eq.1) go to 280
         if (zv.le.xmax) go to 280
         xmax = zv
         in = j
  280    continue
      end do

      if (xmax.le.toler) go to 490

      if ( q(klm1,in) .ne. xmax) then
        q(1:klm2,in) = -q(1:klm2,in)
        q(klm1,in) = xmax
      end if
!
!  Determine the vector to leave the basis.
!
! 300 continue

      if ( iphase .ne. 1 .and. ia .ne. 0 ) then

        xmax = 0.0E+00
        do i=1,ia
          z = abs(q(i,in))
          if ( xmax .lt. z ) then
            xmax = z
            iout = i
          end if
        end do

        if ( toler .lt. xmax ) then

          do j = 1, n2
            z = q(ia,j)
            q(ia,j) = q(iout,j)
            q(iout,j) = z
          end do

          iout = ia
          ia = ia - 1
          pivot = q(iout,in)
          go to 420

        end if

      end if

      kk = 0

      do i=1,klm
        z = q(i,in)
        if ( toler .lt. z ) then
          kk = kk + 1
          res(kk) = q(i,n1)/z
          s(kk) = i
        end if
      end do

  350 continue

      if (kk.le.0) then
        kode = 2
        go to 590
      end if

      xmin = res(1)
      iout = s(1)
      j = 1

      if ( 1 < kk ) then

        do i=2,kk
          if (res(i).lt.xmin) then
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
      ii = q(iout,n2)
      if (iphase.eq.1) go to 400
      if (ii.lt.0) go to 390
      if (iu(2,ii).eq.1) go to 420
      go to 400
  390 iineg = -ii
      if (iu(1,iineg).eq.1) go to 420
  400 ii = iabs(ii)
      cuv = cu(1,ii) + cu(2,ii)
      if (q(klm1,in)-pivot*cuv.le.toler) go to 420
!
!  Bypass intermediate vertices.
!
      do j=js,n1
         z = q(iout,j)
         q(klm1,j) = q(klm1,j) - z*cuv
         q(iout,j) = -z
      end do

      q(iout,n2) = -q(iout,n2)
      go to 350
!
!  Gauss-Jordan elimination.
!
  420 continue

      if ( maxit .le. iter ) then
        kode = 3
        go to 590
      end if

      iter = iter + 1

      do j = js, n1
        if ( j .ne. in ) then
          q(iout,j) = q(iout,j) / pivot
        end if
      end do
!
!  If permitted, use subroutine COL of the description
!  section and replace the following seven statements down
!  to and including statement number 460 by..
!
      if ( .false. ) then

         do j=js,n1
           if ( j .ne. in ) then
             z = -q(iout,j)
             call col(q(1,j), q(1,in), z, iout, klm1)
           end if
         end do

      else

        do j=js,n1
          if ( j.ne.in) then
            z = -q(iout,j)
            do i=1,klm1
              if (i.ne.iout) then
                q(i,j) = q(i,j) + z*q(i,in)
              end if
            end do
          end if
        end do

      end if

      tpivot = -pivot

      do i=1,klm1
        if (i.ne.iout) then
          q(i,in) = q(i,in)/tpivot
        end if
      end do

      q(iout,in) = 1.0E+00 / pivot
      z = q(iout,n2)
      q(iout,n2) = q(klm2,in)
      q(klm2,in) = z
      ii = abs(z)

      if (iu(1,ii).eq.0 .or. iu(2,ii).eq.0) then
        go to 240
      end if

      do i=1,klm2
         z = q(i,in)
         q(i,in) = q(i,js)
         q(i,js) = z
      end do
      js = js + 1
      go to 240
!
!  Test for optimality.
!
  490 continue

      if ( kforce .eq. 0 ) go to 580
      if (iphase.eq.1 .and. q(klm1,n1).le.toler) go to 500
      kforce = 0
      go to 240
!
!  Set up phase 2 costs.
!
  500 continue

      iphase = 2
      cu(1:2,1:nklm) = 0.0E+00
      cu(1:2,n1:nk) = 1.0E+00

      do i=1,klm
         ii = q(i,n2)
         if (ii.gt.0) go to 530
         ii = -ii
         if (iu(2,ii).eq.0) go to 560
         cu(2,ii) = 0.
         go to 540
  530    if (iu(1,ii).eq.0) go to 560
         cu(1,ii) = 0.
  540    ia = ia + 1
         do j=1,n2
            z = q(ia,j)
            q(ia,j) = q(i,j)
            q(i,j) = z
         end do
  560    continue
      end do

      go to 160
  570 if (q(klm1,n1).le.toler) go to 500
      kode = 1
      go to 590
  580 if (iphase.eq.1) go to 570
!
!  Prepare output.
!
  kode = 0

  590 continue

  sum = 0.0D+00
  x(1:n) = 0.0E+00
  res(1:klm) = 0.0E+00

  do i = 1, klm

    ii = q(i,n2)
    sn = 1.0E+00

    if ( ii .le. 0 ) then
      ii = -ii
      sn = -1.0E+00
    end if

    if ( ii .le. n ) then
      x(ii) = sn * q(i,n1)
    else
      res(ii-n) = sn * q(i,n1)
      if ( n1 .le. ii .and. ii .le. nk ) then
        sum = sum + dble ( q(i,n1) )
      end if
    end if

  end do

  error = sum

  return
end
subroutine col ( v1, v2, xmlt, notrow, k )

!*****************************************************************************80
!
!! COL adds a multiple of vector V2 to vector V1.
!
!  Discussion:
!
!    Elements 1 through K of V2 are added to V1, excluding
!    element NOTROW.
!
!  Modified:
!
!    27 December 2007
!
!  Parameters:
!
!    Input/output, real V1(K), ?
!
!    Input, real V2(K), ?
!
!    Input, real XMLT, the multiplier.
!
!    Input, integer NOTROW, the entry of V1 that is not to be changed.
!
!    Input, integer K, the dimension of V1 and V2.
!
  implicit none

  integer k

  integer notrow
  real v1(k)
  real v1_save
  real v2(k)
  real xmlt

  v1_save = v1(notrow)
  v1(1:k) = v1(1:k) + xmlt * v2(1:k)
  v1(notrow) = v1_save

  return
end
subroutine l1orth ( z, mdim, m, mp3dim, mp1, np1, np3, mp1np1, eps1, eps2, &
  itmax, iflag, y, q, r, c, ic, ir, e )

!*****************************************************************************80
!
!! L1ORTH carries out orthogonal regression for the L1 norm.
!
!  Discussion:
!
!    We assume we have an overdetermined M by N linear system
!
!      A * X = B.
!
!    We seek a solution X which minimizes the sum of orthogonal
!    distances:
!
!      F1(X) = ( sum ( 1 <= I <= M ) abs ( A(I,1:N) * X(1:N) - B(I) ) )
!            / sqrt ( 1 + X(1:N) * X(1:N) )
!
!    which can be reformulated as seek a vector Y of unit L2 norm
!    which minimizes the L1 norm:
!
!      || ( A | -b ) * Y ||.
!
!    This routine repeatedly calls CL1, which implements an
!    algorithm of Barrodale and Roberts for constrained L1 minimization.
!
!  Modified:
!
!    09 February 2003
!
!  Reference:
!
!    I Barrodale and F Roberts,
!    Algorithm 552: 
!    Solution of the Constrained L1 Approximation Problem,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 2, pages 231-235, 1980.
!
!    Helmut Spaeth,
!    Mathematical Algorithms for Linear Regresssion,
!    Academic Press, 1993, pages 280-285.
!
!  Parameters:
!
!    Input, real Z(MDIM,NP1), contains the system matrix A in the
!    first M rows and N columns, and -B (the right hand side
!    of the linear system) in column N+1.
!
!    Input, integer MDIM, the first dimension of Z, which must
!    be at least M.
!
!    Input, integer M, the number of rows of the matrix A.
!
!    Input, integer MP3DIM, >= m + 3.
!
!    Input, integer MP1, = m + 1.
!
!    Input, integer NP1, = n + 1.
!
!    Input, integer NP3, = n + 3.
!
!    Input,  integer MP1NP1, = ( m + 1 ) + ( n + 1 ).
!
!    Input, real EPS1, an error tolerance for CL1.
!
!    Input, real EPS2, an error tolerance for the outer iteration.
!
!    Input/output, integer ITMAX.  On input, the maximum number
!    of outer iterations that may be taken.  On output, the number
!    of outer iterations that were actually taken.
!
!    Output, integer IFLAG, error flag.
!    0, normal exit.
!    1, 2 or 3: an error return from Cl1.
!    4, no convergence with tolerance EPS2 after ITMAX iterations.
!
!    Workspace, real Y(NP3).
!
!    Workspace, real Q(MP3DIM,NP3).
!
!    Workspace, real R(MP1).
!
!    Workspace, real C(2,MP1NP1).
!
!    Workspace, real IC(2,MP1NP1).
!
!    Workspace, real IR(MP1).
!
!    Workspace, real E(NP1).
!
  implicit none

  integer mdim
  integer mp1
  integer mp1np1
  integer mp3dim
  integer np1
  integer np3

  real c(2,mp1np1)
  real e(np1)
  real eps1
  real eps2
  real error
  integer i
  integer ic(2,mp1np1)
  integer iflag
  integer ir(mp1)
  integer it
  integer itmax
  integer itn
  integer itnorm
  integer k
  integer m
  integer n
  integer np2
  real q(mp3dim,np3)
  real r(mp1)
  real s
  real u
  real v
  real y(np3)
  real yk
  real z(mdim,np1)

  it = 0
  itnorm = 10 * mp1
  n = np1 - 1
  np2 = np1 + 1

  y(1:n) = 0.0E+00
  y(n+1) = 1.0E+00

  e(1:n) = 0.0E+00
  e(n+1) = 1.0E+00

  do

    it = it + 1

    if ( itmax < it ) then
      iflag = 4
      return
    end if

    q(1:m,1:np1) = z(1:m,1:np1)
    q(1:m,np2) = 0.0E+00
    q(mp1,1:np1) = y(1:np1)
    q(mp1,np2) = 1.0E+00

    itn = itnorm
    iflag = 0

    call cl1 ( m, 1, 0, np1, mp1, mp3dim, mp1np1, np3, q, iflag, &
      eps1, itn, y, r, error, c, ic, ir )

    if ( iflag /= 0 ) then
      return
    end if

    s = 0.0E+00
    do k = 1, np1
      yk = y(k)
      s = s + yk * yk
    end do

    s = 1.0E+00 / sqrt ( s )
    u = 0.0E+00
    v = 0.0E+00

    do k = 1, np1
      yk = s * y(k)
      y(k) = yk
      u = u + abs ( yk )
      v = v + abs ( yk - e(k) )
      e(k) = yk
    end do

    if ( v <= eps2 * u ) then
      exit
    end if

  end do

  itmax = it

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
