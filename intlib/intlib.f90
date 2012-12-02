subroutine avint ( ntab, xtab, ytab, a, b, result )

!*****************************************************************************80
!
!! AVINT estimates the integral of unevenly spaced data.
!
!  Discussion:
!
!    The data is given as NTAB pairs of values 
!    ( XTAB(1:NTAB), YTAB(1:NTAB) ).
!
!    The quadrature method uses overlapping parabolas and smoothing.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Paul Hennion,
!    Algorithm 77:
!    Interpolation, Differentiation and Integration,
!    Communications of the ACM,
!    Volume 5, page 96, 1962.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of entries in XTAB and
!    YTAB.  NTAB must be at least 2.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, real ( kind = 8 ) YTAB(NTAB), the function values,
!    YTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ba
  real ( kind = 8 ) bb
  real ( kind = 8 ) bc
  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  real ( kind = 8 ) cc
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inlft
  integer ( kind = 4 ) inrt
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) istop
  real ( kind = 8 ) result
  real ( kind = 8 ) slope
  real ( kind = 8 ) syl
  real ( kind = 8 ) syl2
  real ( kind = 8 ) syl3
  real ( kind = 8 ) syu
  real ( kind = 8 ) syu2
  real ( kind = 8 ) syu3
  real ( kind = 8 ) term1
  real ( kind = 8 ) term2
  real ( kind = 8 ) term3
  real ( kind = 8 ) total
  real ( kind = 8 ) x1
  real ( kind = 8 ) x12
  real ( kind = 8 ) x13
  real ( kind = 8 ) x2
  real ( kind = 8 ) x23
  real ( kind = 8 ) x3
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) ytab(ntab)

  result = 0.0D+00

  if ( a == b ) then
    return
  end if

  if ( b < a ) then
  end if

  if ( ntab < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Fatal error!'
    write ( *, '(a,i8)' ) '  NTAB is less than 3.  NTAB = ', ntab
    stop
  end if

  do i = 2, ntab
 
    if ( xtab(i) <= xtab(i-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AVINT - Fatal error!'
      write ( *, '(a)' ) '  XTAB(I) is not greater than XTAB(I-1).'
      write ( *, '(a,i8)' ) '  Here, I = ', I
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
      write ( *, '(a,g14.6)' ) '  XTAB(I) =   ', xtab(i)
      stop
    end if
 
  end do
!
!  Special case for NTAB = 2.
!
  if ( ntab == 2 ) then
    slope = ( ytab(2) - ytab(1) ) / ( xtab(2) - xtab(1) )
    fa = ytab(1) + slope * ( a - xtab(1) )
    fb = ytab(2) + slope * ( b - xtab(2) )
    result = 0.5D+00 * ( fa + fb ) * ( b - a )
    return
  end if

  if ( xtab(ntab-2) < a .or. b < xtab(3) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Fatal error!'
    write ( *, '(a)' ) '  There were less than 3 function values'
    write ( *, '(a)' ) '  between the limits of integration.'
    stop
  end if

  i = 1
  do

    if ( a <= xtab(i) ) then
      exit
    end if

    i = i + 1

  end do

  inlft = i

  i = ntab

  do

    if ( xtab(i) <= b ) then
      exit
    end if

    i = i - 1

  end do

  inrt = i

  if ( inrt - inlft < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Fatal error!'
    write ( *, '(a)' ) '  There were less than 3 function values'
    write ( *, '(a)' ) '  between the limits of integration.'
    stop
  end if

  if ( inlft == 1 ) then
    istart = 2
  else
    istart = inlft
  end if

  if ( inrt == ntab ) then
    istop = ntab - 1
  else
    istop = inrt
  end if

  total = 0.0D+00

  syl = a
  syl2 = syl * syl
  syl3 = syl2 * syl

  do i = istart, istop

    x1 = xtab(i-1)
    x2 = xtab(i)
    x3 = xtab(i+1)

    x12 = x1 - x2
    x13 = x1 - x3
    x23 = x2 - x3

    term1 =   ( ytab(i-1) ) / ( x12 * x13 )
    term2 = - ( ytab(i)   ) / ( x12 * x23 )
    term3 =   ( ytab(i+1) ) / ( x13 * x23 )

    ba = term1 + term2 + term3
    bb = - ( x2 + x3 ) * term1 - ( x1 + x3 ) * term2 - ( x1 + x2 ) * term3
    bc = x2 * x3 * term1 + x1 * x3 * term2 + x1 * x2 * term3

    if ( i == istart ) then
      ca = ba
      cb = bb
      cc = bc
    else
      ca = 0.5D+00 * ( ba + ca )
      cb = 0.5D+00 * ( bb + cb )
      cc = 0.5D+00 * ( bc + cc )
    end if

    syu = x2
    syu2 = syu * syu
    syu3 = syu2 * syu

    total = total + ca * ( syu3 - syl3 ) / 3.0D+00 &
                  + cb * ( syu2 - syl2 ) / 2.0D+00 &
                  + cc * ( syu  - syl )
    ca = ba
    cb = bb
    cc = bc

    syl  = syu
    syl2 = syu2
    syl3 = syu3

  end do

  syu = b
  syu2 = syu * syu
  syu3 = syu2 * syu

  result = total + ca * ( syu3 - syl3 ) / 3.0D+00 &
                 + cb * ( syu2 - syl2 ) / 2.0D+00 &
                 + cc * ( syu  - syl  )
  
  return
end
subroutine cadre ( func, a, b, abserr, relerr, error, result, ind )

!*****************************************************************************80
!
!! CADRE estimates the integral of F(X) from A to B.
!
!  Discussion:
!
!    CADRE is the Cautious Adaptive Romberg Extrapolator.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Carl DeBoor, John Rice,
!    CADRE: An algorithm for numerical quadrature,
!    in Mathematical Software,
!    edited by John Rice,
!    Academic Press, 1971.
!    ISBN: 012587250X,
!    LC: QA1.M766,
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external parameter in the
!    calling program, write a function routine of the form 
!      FUNCTION FUNC ( X ) 
!    which evaluates the function at X, and pass the name of the function
!    in FUNC.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, real ( kind = 8 ) ABSERR, the absolute error tolerance.
!
!    Input, real ( kind = 8 ) RELERR, the relative error tolerance.
!
!    Output, real ( kind = 8 ) ERROR, an estimate of the absolute error.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!    Output, integer ( kind = 4 ) IND, reliability indicator.
!    If IND <= 2, RESULT is very reliable.  Higher values of
!    IND indicate less reliable values of RESULT.
!
  implicit none

  integer ( kind = 4 ), parameter :: mxstge = 30
  integer ( kind = 4 ), parameter :: maxtbl = 10
  integer ( kind = 4 ), parameter :: maxts = 2049

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) ait(maxtbl)
  logical              aitken
  real ( kind = 8 ), parameter :: aitlow = 1.1D+00
  real ( kind = 8 ), parameter :: aittol = 0.1D+00
  real ( kind = 8 ) astep
  real ( kind = 8 ) b
  real ( kind = 8 ) beg
  real ( kind = 8 ) begin(mxstge)
  real ( kind = 8 ) bma
  real ( kind = 8 ) curest
  real ( kind = 8 ) dif(maxtbl)
  real ( kind = 8 ) diff
  real ( kind = 8 ) end
  real ( kind = 8 ) ergoal
  real ( kind = 8 ) erra
  real ( kind = 8 ) errer
  real ( kind = 8 ) error
  real ( kind = 8 ) errr
  real ( kind = 8 ) est(mxstge)
  real ( kind = 8 ) fbeg
  real ( kind = 8 ) fbeg2
  real ( kind = 8 ) fend
  real ( kind = 8 ) fextm1
  real ( kind = 8 ) fextrp
  real ( kind = 8 ) finis(mxstge)
  real ( kind = 8 ) fn
  real ( kind = 8 ) fnsize
  real ( kind = 8 ), external :: func
  logical h2conv
  real ( kind = 8 ) h2next
  real ( kind = 8 ) h2tfex
  real ( kind = 8 ), parameter :: h2tol = 0.15D+00
  real ( kind = 8 ) hovn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ibegs(mxstge)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) istage
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) istep2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lm1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nnleft
  real ( kind = 8 ) prever
  real ( kind = 8 ) r(maxtbl)
  logical reglar
  logical reglsv(mxstge)
  real ( kind = 8 ) relerr
  real ( kind = 8 ) result
  logical right
  real ( kind = 8 ) rn(4)
  real ( kind = 8 ) rnderr
  real ( kind = 8 ) sing
  real ( kind = 8 ) singnx
  real ( kind = 8 ) slope
  real ( kind = 8 ) stage
  real ( kind = 8 ) step
  real ( kind = 8 ) stepmn
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sumabs
  real ( kind = 8 ) t(maxtbl,maxtbl)
  real ( kind = 8 ) tabs
  real ( kind = 8 ) tabtlm
  real ( kind = 8 ), parameter :: tljump = 0.01D+00
  real ( kind = 8 ) ts(2049)
  real ( kind = 8 ) vint

  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  begin(1:mxstge) = 0.0D+00
  est(1:mxstge) = 0.0D+00
  finis(1:mxstge) = 0.0D+00
  ibegs(1:mxstge) = 0
  reglsv(1:mxstge) = .false.
 
  vint = 0.0D+00
 
  rn(1:4) = (/ 0.7142005D+00, 0.3466282D+00, 0.8437510D+00, 0.1263305D+00 /)
 
  rnderr = epsilon ( rnderr )
  result = 0.0D+00
  error = 0.0D+00
  ind = 1
  bma = abs ( b - a )
  errr = min ( 0.1D+00, max ( abs ( relerr ), 10.0D+00*rnderr) )
  erra = abs ( abserr )
  stepmn = max ( bma / 2**mxstge, max ( bma, abs ( a ), abs ( b ) ) * rnderr )
  stage = 0.5D+00
  istage = 1
  curest = 0.0D+00
  fnsize = 0.0D+00
  prever = 0.0D+00
  reglar = .false.
  beg = a
  fbeg = func ( beg ) / 2.0D+00
  ts(1) = fbeg
  ibeg = 1
  end = b
  fend = func ( end ) / 2.0D+00
  ts(2) = fend
  iend = 2
 
10 continue
 
  right = .false.
 
20 continue

  step = end - beg
  astep = abs ( step )
 
  if ( astep < stepmn ) then
    ind = 5
    result = curest + vint
    return
  end if
 
  t(1,1) = fbeg + fend
  tabs = abs ( fbeg ) + abs ( fend )
  l = 1
  n = 1
  h2conv = .false.
  aitken = .false.
  go to 40
 
30 continue
 
40 continue
 
  lm1 = l
  l = l + 1
  n2 = n * 2
  fn = n2
  istep = ( iend - ibeg ) / n

  if ( 1 < istep ) then
    go to 60
  end if

  ii = iend
  iend = iend + n

  if ( maxts < iend ) then
    go to 440
  end if

  hovn = step / fn
 
  iii = iend
  do i = 1, n2, 2
    ts(iii) = ts(ii)
    ts(iii-1) = func ( end - i * hovn )
    iii = iii-2
    ii = ii-1
  end do
 
  istep = 2
 
60 continue
 
  istep2 = ibeg + istep / 2
 
  sum1 = 0.0D+00
  sumabs = 0.0D+00
  do i = istep2, iend, istep
    sum1 = sum1 + ts(i)
    sumabs = sumabs + abs ( ts(i) )
  end do
 
  t(l,1) = t(l-1,1) / 2.0D+00 + sum1 / fn
  tabs = tabs / 2.0D+00 + sumabs / fn
 
  n = n2
  it = 1
  vint = step * t(l,1)
  tabtlm = tabs * rnderr
  fnsize = max ( fnsize, abs ( t(l,1) ) )
  ergoal = max ( astep * rnderr * fnsize, &
    stage * max ( erra , errr * abs ( curest + vint ) ) )
  fextrp = 1.0D+00
  do i = 1, lm1
    fextrp = fextrp * 4.0D+00
    t(i,l) = t(l,i) - t(l-1,i)
    t(l,i+1) = t(l,i) + t(i,l) / ( fextrp - 1.0D+00 )
  end do
 
  errer = astep * abs ( t(1,l) )

  if ( 2 < l ) then
    go to 90
  end if

  if ( abs ( t(1,2) ) <= tabtlm ) then
    go to 290
  end if

  go to 40
 
90 continue
 
  do i = 2, lm1

    if ( tabtlm < abs ( t(i-1,l) ) ) then
      diff = t(i-1,lm1) / t(i-1,l)
    else
      diff = 0.0D+00
    end if

    t(i-1,lm1) = diff

  end do
 
  if ( abs ( 4.0D+00 - t(1,lm1) ) <= h2tol ) then
    go to 130
  end if

  if ( t(1,lm1) == 0.0D+00 ) then
    go to 120
  end if

  if ( abs ( 2.0D+00 - abs ( t(1,lm1) ) ) < tljump ) then
    go to 280
  end if

  if ( l == 3 ) then
    go to 30
  end if

  h2conv = .false.

  if ( abs ( ( t(1,lm1) - t(1,l-2) ) / t(1,lm1) ) <= aittol ) then
    go to 160
  end if
  
  if ( .not. reglar .and. l == 4 ) then
    go to 30
  end if
 
120 continue
 
  if ( errer <= ergoal ) then
    go to 310
  end if

  go to 380

130 continue

  if ( .not. h2conv ) then
    aitken = .false.
    h2conv = .true.
  end if

140 continue

  fextrp = 4.0D+00

150 continue

  it = it + 1
  vint = step * t(l,it)
  errer = abs ( step / ( fextrp - 1.0D+00 ) * t(it-1,l))

  if ( errer <= ergoal ) then
    go to 340
  end if

  if ( it == lm1 ) then
    go to 270
  end if

  if ( t(it,lm1) == 0.0D+00 ) then
    go to 150
  end if

  if ( t(it,lm1) <= fextrp ) then
    go to 270
  end if

  if ( abs ( t(it,lm1) / 4.0D+00 - fextrp ) / fextrp < aittol ) then
    fextrp = fextrp * 4.0D+00
  end if

  go to 150
 
160 continue

  if ( t(1,lm1) < aitlow ) then
    go to 380
  end if
 
  if ( .not. aitken ) then
    h2conv = .false.
    aitken = .true.
  end if
 
170 continue

  fextrp = t(l-2,lm1)

  if ( 4.5D+00 < fextrp ) then
    go to 140
  end if

  if ( fextrp < aitlow ) then
    go to 380
  end if

  if ( h2tol < abs ( fextrp - t(l-3,lm1) ) / t(1,lm1) ) then
    go to 380
  end if

  sing = fextrp
  fextm1 = fextrp - 1.0D+00

  ait(1) = 0.0D+00
  do i = 2, l
    ait(i) = t(i,1) + (t(i,1)-t(i-1,1)) / fextm1
    r(i) = t(1,i-1)
    dif(i) = ait(i) - ait(i-1)
  end do

  it = 2

190 continue

  vint = step * ait(l)

200 continue

  errer = errer / fextm1
 
  if ( errer <= ergoal ) then
    ind = max ( ind, 2 )
    go to 340
  end if
 
210 continue

  it = it + 1

  if ( it == lm1 ) then
    go to 270
  end if

  if ( it <= 3 ) then
    h2next = 4.0D+00
    singnx = 2.0D+00 * sing
  end if

  if ( h2next < singnx ) then
    go to 230
  end if

  fextrp = singnx
  singnx = 2.0D+00 * singnx
  go to 240

230 continue

  fextrp = h2next
  h2next = 4.0D+00 * h2next

240 continue
 
  do i = it, lm1
    if ( tabtlm < abs ( dif(i+1) ) ) then
      r(i+1) = dif(i) / dif(i+1)
    else
      r(i+1) = 0.0D+00
    end if
  end do
 
  h2tfex = -h2tol * fextrp

  if ( r(l) - fextrp < h2tfex ) then
    go to 270
  end if

  if ( r(l-1) - fextrp < h2tfex ) then
    go to 270
  end if

  errer = astep * abs ( dif(l) )
  fextm1 = fextrp - 1.0D+00
  do i = it, l
    ait(i) = ait(i)+dif(i) / fextm1
    dif(i) = ait(i)-ait(i-1)
  end do
 
  go to 190
 
270 continue

  fextrp = max ( prever / errer, aitlow )
  prever = errer
  if ( l < 5 ) then
    go to 40
  end if

  if ( 2.0D+00 < l - it .and. istage < mxstge ) then
    go to 370
  end if

  if ( errer / fextrp**( maxtbl - l ) < ergoal ) then
    go to 40
  end if

  go to 370
 
280 continue

  if ( ergoal < errer ) then
    go to 370
  end if

  diff = abs ( t(1,l) ) * 2.0D+00 * fn
  go to 340
 
290 continue

  slope = ( fend - fbeg ) * 2.0D+00
  fbeg2 = fbeg * 2.0D+00
 
  do i = 1, 4
    diff = abs ( func ( beg + rn(i) * step ) - fbeg2 - rn(i) * slope )
    if ( tabtlm < diff ) then
      go to 330
    end if
  end do
 
  go to 340
 
310 continue

  slope = ( fend - fbeg ) * 2.0D+00
  fbeg2 = fbeg * 2.0D+00
  i = 1
 
320 continue

  diff = abs ( func ( beg + rn(i) * step ) - fbeg2 - rn(i) * slope )
 
330 continue

  errer = max ( errer, astep * diff )

  if ( ergoal < errer ) then
    go to 380
  end if

  i = i+1

  if ( i <= 4 ) then
    go to 320
  end if

  ind = 3
 
340 continue

  result = result + vint
  error = error + errer
 
350 continue

  if ( right ) then
    go to 360
  end if

  istage = istage - 1

  if ( istage == 0 ) then
    return
  end if

  reglar = reglsv(istage)
  beg = begin(istage)
  end = finis(istage)
  curest = curest - est(istage+1) + vint
  iend = ibeg - 1
  fend = ts(iend)
  ibeg = ibegs(istage)
  go to 400
 
360 continue

  curest = curest + vint
  stage = stage * 2.0D+00
  iend = ibeg
  ibeg = ibegs(istage)
  end = beg
  beg = begin(istage)
  fend = fbeg
  fbeg = ts(ibeg)
  go to 10
 
370 continue

  reglar = .true.
 
380 continue
 
  if ( istage == mxstge ) then
    ind = 5
    result = curest + vint
    return
  end if
 
390 continue

  if ( right ) then
    go to 410
  end if

  reglsv(istage+1) = reglar
  begin(istage) = beg
  ibegs(istage) = ibeg
  stage = stage / 2.0D+00

400 continue

  right = .true.
  beg = ( beg + end ) / 2.0D+00
  ibeg = ( ibeg + iend ) / 2
  ts(ibeg) = ts(ibeg) / 2.0D+00
  fbeg = ts(ibeg)
  go to 20

410 continue

  nnleft = ibeg - ibegs(istage)
  if ( maxts <= end + nnleft ) then
    go to 440
  end if

  iii = ibegs(istage)
  ii = iend
  do i = iii, ibeg
    ii = ii + 1
    ts(ii) = ts(i)
  end do
 
  do i = ibeg, ii
    ts(iii) = ts(i)
    iii = iii + 1
  end do
 
  iend = iend + 1
  ibeg = iend - nnleft
  fend = fbeg
  fbeg = ts(ibeg)
  finis(istage) = end
  end = beg
  beg = begin(istage)
  begin(istage) = end
  reglsv(istage) = reglar
  istage = istage + 1
  reglar = reglsv(istage)
  est(istage) = vint
  curest = curest + est(istage)
  go to 10

440 continue

  ind = 4

460 continue

  result = curest + vint

  return
end
subroutine chinsp ( func, a, b, epsin, epsout, result )

!*****************************************************************************80
!
!! CHINSP estimates an integral using a modified Clenshaw-Curtis scheme.
!
!  Discussion:
!
!    The integral is approximated by Chebyshev polyonomials over each
!    subinterval.  These are integrated to give the approximate integral.
!    If the error estimate is unsatisfactory, the integration is repeated
!    with smaller intervals.
!
!    The internal parameter NUPPER is currently set to 9,
!    corresponding to 1024 subintervals for the unfolded integral,
!    and 1025 function evaluations.  This parameter may be changed
!    if necessary.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Tore Havie,
!    On a Modification of the Clenshaw Curtis Quadrature Rule,
!    BIT,
!    Volume 9, Number 4, December 1969, pages 338-350.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      FUNCTION FUNC ( X )
!    which evaluates the function at the point X.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, real ( kind = 8 ) EPSIN, the relative error tolerance.
!
!    Output, real ( kind = 8 ) EPSOUT, estimated integration error.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: nupper = 9

  real ( kind = 8 ) a
  real ( kind = 8 ) a0
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) acof(257)
  real ( kind = 8 ) alf
  real ( kind = 8 ) alfnj
  real ( kind = 8 ) alfno
  real ( kind = 8 ) b
  real ( kind = 8 ) bcof(257)
  real ( kind = 8 ) bet
  real ( kind = 8 ) betnj
  real ( kind = 8 ) betno
  real ( kind = 8 ) bounds
  real ( kind = 8 ) ccof(513)
  real ( kind = 8 ) cof
  real ( kind = 8 ) cofmax
  real ( kind = 8 ) const1
  real ( kind = 8 ) const2
  real ( kind = 8 ) deln
  real ( kind = 8 ) deltan
  real ( kind = 8 ) epsin
  real ( kind = 8 ) epsout
  real ( kind = 8 ) error
  real ( kind = 8 ) etank
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) gamman
  real ( kind = 8 ) hnstep
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ksign
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncof
  integer ( kind = 4 ) nhalf
  integer ( kind = 4 ) nn
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) rk
  real ( kind = 8 ) rn
  real ( kind = 8 ) rnderr
  real ( kind = 8 ) rounde
  real ( kind = 8 ) tend
  real ( kind = 8 ) tnew
  real ( kind = 8 ) triarg
  real ( kind = 8 ) umid
  real ( kind = 8 ) wmean
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xplus
  real ( kind = 8 ) xsink

  if ( a == b ) then
    result = 0.0D+00
    return
  end if
!
!  ROUNDE = RNDERR*(R1+R2*N), where R1, R2 are two empirical constants.
!
!  Set coefficients in formula for accumulated roundoff error.
!  N is the current number of function values used.
!
  rnderr = epsilon ( 1.0D+00 )
 
  r1 = 1.0D+00
  r2 = 2.0D+00
  error = epsin
!
!  Integration interval parameters.
!
  alf = 0.5D+00 * ( b - a )
  bet = 0.5D+00 * ( b + a )
!
!  Parameters for trigonometric recurrence relations.
!
  triarg = atan ( 1.0D+00 )
  alfno = -1.0D+00
!
!  Parameters for integration stepsize and loops.
!
  rn = 2.0D+00
  n = 2
  nhalf = 1
  hnstep = 1.0D+00
!
!  Initial calculation for the end-point approximation.
!
  const1 = 0.5D+00 * ( func ( a ) + func ( b ) )
  const2 = func ( bet )
  acof(1) = 0.5D+00 * ( const1 + const2 )
  acof(2) = 0.5D+00 * ( const1 - const2 )
  bcof(2) = acof(2)
  tend = 2.0D+00 * ( acof(1) - acof(2) / 3.0D+00 )
!
!  Start actual calculations.
!
  do i = 1, nupper
!
!  Compute function values.
!
    const1 = -sin ( triarg )
    const2 = 0.5D+00 * alfno / const1
    alfno = const1
    betno = const2
    gamman = 1.0D+00 - 2.0D+00 * alfno**2
    deltan = -2.0D+00 * alfno * betno
    bcof(1) = 0.0D+00
 
    do j = 1, nhalf
      alfnj = gamman * const1 + deltan * const2
      betnj = gamman * const2 - deltan * const1
      xplus = alf * alfnj + bet
      xmin = -alf * alfnj + bet
      ccof(j) = func ( xplus ) + func ( xmin )
      bcof(1) = bcof(1) + ccof(j)
      const1 = alfnj
      const2 = betnj
    end do
 
    bcof(1) = 0.5D+00 * hnstep * bcof(1)
!
!  Calculation of first B-coefficient finished compute the higher
!  coefficients if NHALF greater than one.
!
    if ( nhalf <= 1 ) then
      go to 60
    end if

    const1 = one
    const2 = 0.0D+00
    ncof = nhalf - 1
    ksign = -1
 
    do k = 1, ncof
!
!  Compute trigonometric sum for B-coefficient.
!
      etank = gamman * const1 - deltan * const2
      xsink = gamman * const2 + deltan * const1
      cof = 2.0D+00 * ( 2.0D+00 * etank**2 - 1.0D+00 )
      a2 = 0.0D+00
      a1 = 0.0D+00
      a0 = ccof(nhalf)
 
      do j = 1, ncof
        a2 = a1
        a1 = a0
        index = nhalf-j
        a0 = ccof(index) + cof * a1 - a2
      end do
 
      bcof(k+1) = hnstep * ( a0 - a1 ) * etank
      bcof(k+1) = ksign * bcof(k+1)
      ksign = -ksign
      const1 = etank
      const2 = xsink
 
    end do
!
!  Calculation of B-coefficients finished.
!
!  Compute new modified mid-point approximation when the interval
!  of integration is divided in N equal sub intervals.
!
60  continue
 
    umid = 0.0D+00
    rk = rn
    nn = nhalf + 1
    do k = 1, nn
      index = nn + 1 - k
      umid = umid + bcof(index) / ( rk**2 - one )
      rk = rk - 2.0D+00
    end do
 
    umid = -2.0D+00 * umid
!
!  Compute new C-coefficients for end-point approximation and largest
!  absolute value of coefficients.
!
    nn = n+2
    cofmax = 0.0D+00
 
    do j = 1, nhalf
      index = nn - j
      ccof(j) = 0.5D+00 * ( acof(j) + bcof(j) )
      ccof(index) = 0.5D+00 * ( acof(j) - bcof(j) )
      const1 = abs ( ccof(j) )
      cofmax = max ( cofmax, const1 )
      const1 = abs ( ccof(index) )
      cofmax = max ( cofmax, const1 )
    end do
 
    ccof(nhalf+1) = acof(nhalf+1)
!
!  Calculation of new coefficients finished.
!
!  Compute new end-point approximation when the interval of
!  integration is divided in 2N equal sub intervals.
!
    wmean = 0.5D+00 * ( tend + umid )
    bounds = 0.5D+00 * ( tend - umid )
    deln = 0.0D+00
    rk = 2.0D+00 * rn
    do j = 1, nhalf
      index = n + 2 - j
      deln = deln + ccof(index) / ( rk**2 - one )
      rk = rk - 2.0D+00
    end do
 
    deln = -2.0D+00 * deln
    tnew = wmean + deln
    epsout = abs ( bounds / tnew )

    if ( cofmax < rnderr ) then
      go to 160
    end if

    rounde = rnderr * ( r1 + r2 * rn )
    epsout = max ( epsout, rounde )
    error = max ( error, rounde )

    if ( error < epsout ) then
      go to 160
    end if
!
!  Required accuracy obtained or the maximum number of function
!  values used without obtaining the required accuracy.
!
120 continue
 
    n = 2 * n + 1
    tend = alf * ( tend + deln )
    umid = alf * ( umid + deln )
    deln = alf * deln
    result = alf * tnew
    return
!
!  If I = NUPPER then the required accuracy is not obtained.
!
160 continue
 
    if ( i == nupper ) then
      go to 120
    end if
 
    acof(1:n) = ccof(1:n)
    acof(n+1) = ccof(n+1)
    bcof(n+1) = ccof(n+1)
    tend = tnew
    nhalf = n
    n = 2 * n
    rn = 2.0D+00 * rn
    hnstep = 0.5D+00 * hnstep
    triarg = 0.5D+00 * triarg
 
  end do
 
  return
end
subroutine class ( kind, n, alpha, beta, b, a, muzero )

!*****************************************************************************80
!
!! CLASS sets recurrence coeeficients for various orthogonal polynomials.
!
!  Discussion:
!
!    CLASS supplies the coefficients A(J), B(J) of the recurrence relation
!
!      B(J)*P(J) (X) = (X-A(J))*P(J-1)(X) - B(J-1)*P(J-2)(X)
!
!    for the various classical (normalized) orthogonal polynomials,
!    and the zero-th moment
!
!      MUZERO = Integral W(X) DX
!
!    of the given polynomial's weight function W(X).  Since the
!    polynomials are orthonormalized, the tridiagonal matrix is
!    guaranteed to be symmetric.
!
!  Modified:
!
!    18 December 2002
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA, parameters needed for Laguerre 
!    and Jacobi polynomials.
!
!    Input, integer ( kind = 4 ) KIND, specifies which polynomial is to be handled:
!    1: Legendre polynomials P(X) on (-1, +1),
!    W(X) = 1.
!    2: Chebyshev polynomials of the first kind T(X) on (-1, +1),
!    W(X) = 1 / SQRT(1 - X*X)
!    3: Chebyshev polynomials of the second kind U(X) on (-1, +1),
!    W(X) = SQRT(1 - X*X)
!    4: Hermite polynomials H(X) on (-infinity,+infinity),
!    W(X) = EXP(-X**2)
!    5: Jacobi polynomials P(ALPHA,BETA)(X) on (-1, +1),
!    W(X) = (1-X)**ALPHA + (1+X)**BETA,
!    ALPHA and BETA greater than -1.
!    6: Laguerre polynomials, L(ALPHA)(X) on (0, +infinity),
!    W(X) = EXP(-X) * X**ALPHA,
!    ALPHA greater than -1.
!
!    Input, integer ( kind = 4 ) N, specifies the number of coefficients to
!    calculate.
!
!    Input, real ( kind = 8 ) ALPHA, the value of the ALPHA parameter,
!    required only for Jacobi or Laguerre polynomials.
!
!    Input, real ( kind = 8 ) BETA, the value of the BETA parameter,
!    required only for Jacobi polynomials.
!
!    Output, real ( kind = 8 ) B(N-1), the offdiagonal coefficients.
!
!    Output, real ( kind = 8 ) A(N), the diagonal coefficients.
!
!    Output, real ( kind = 8 ) MUZERO, the zero-th moment, Integral W(X) DX,
!    of the polynomial's weight function over its interval of
!    definition.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) abi
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n-1)
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ) muzero
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
!
!  KIND = 1:
!
!  Legendre polynomials P(X) on (-1, +1),
!  W(X) = 1.
!
  if ( kind == 1 ) then
 
    muzero = 2.0D+00
 
    a(1:n) = 0.0D+00
 
    do i = 1, n-1
      b(i) = real ( i, kind = 8 ) &
        / sqrt ( 4.0D+00 * real ( i * i, kind = 8 ) - 1.0D+00 )
    end do
!
!  KIND = 2:
!
!  Chebyshev polynomials of the first kind T(X) on (-1, +1),
!  W(X) = 1 / SQRT(1 - X*X)
!
  else if ( kind == 2 ) then
 
    muzero = pi
    a(1:n) = 0.0D+00
    b(1) = sqrt ( 0.5D+00 )
    b(1:n-1) = 0.5D+00
!
!  KIND = 3:
!
!  Chebyshev polynomials of the second kind U(X) on (-1, +1),
!  W(X) = SQRT(1 - X*X)
!
  else if ( kind == 3 ) then
 
    muzero = pi / 2.0D+00
    a(1:n) = 0.0D+00
    b(1:n-1) = 0.5D+00
!
!  KIND = 4:
!
!  Hermite polynomials H(X) on (-infinity,+infinity),
!  W(X) = EXP(-X**2)
!
  else if ( kind == 4 ) then
 
    muzero = sqrt ( pi )
    a(1:n) = 0.0D+00
    do i = 1, n-1
      b(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
    end do
!
!  KIND = 5:
!
!  Jacobi polynomials P(ALPHA,BETA)(X) on (-1, +1),
!  W(X) = (1-X)**ALPHA + (1+X)**BETA,
!  ALPHA and BETA greater than -1
!
  else if ( kind == 5 ) then
 
    muzero = 2.0D+00**( alpha + beta + 1.0D+00 ) * gamma ( alpha + 1.0D+00 ) &
      * gamma ( beta + 1.0D+00 ) / gamma ( 2.0D+00 + alpha + beta )
 
    do i = 1, n
      a(i) = ( beta**2 - alpha**2 ) / &
        ( ( 2.0D+00 * real ( i - 1, kind = 8 ) + alpha + beta ) &
        * ( 2.0D+00 * real ( i, kind = 8 ) + alpha + beta ) )
    end do
 
    abi = 2.0D+00 + alpha + beta
    b(1) = sqrt ( 4.0D+00 * ( 1.0D+00 + alpha ) * ( 1.0D+00 + beta ) &
      / ( ( abi + 1.0D+00 ) * abi * abi ) )
 
    do i = 2, n-1
      abi = real ( 2 * i ) + alpha + beta
      b(i) = sqrt ( 4.0D+00 * real ( i, kind = 8 ) &
        * ( real ( i, kind = 8 ) + alpha ) &
        * ( real ( i, kind = 8 ) + beta ) &
        * ( real ( i, kind = 8 ) + alpha + beta ) / &
        ( ( abi * abi - 1.0D+00 ) * abi * abi ) )
    end do
!
!  KIND = 6:
!
!  Laguerre polynomials
!
!  L(ALPHA)(X) on (0, +infinity),
!  W(X) = EXP(-X) * X**ALPHA,
!  ALPHA greater than -1.
!
  else if ( kind == 6 ) then
 
    muzero = gamma ( alpha + 1.0D+00 )
 
    do i = 1, n
      a(i) = 2.0D+00 * real ( i, kind = 8 ) - 1.0D+00 + alpha
    end do
 
    do i = 1, n-1
      b(i) = sqrt ( real ( i, kind = 8 ) * ( real ( i, kind = 8 ) + alpha ) )
    end do
 
  end if
 
  return
end
subroutine cspint ( ntab, xtab, ftab, a, b, y, e, work, result )

!*****************************************************************************80
!
!! CSPINT estimates the integral of a tabulated function.
!
!  Discussion:
!
!    The routine is given the value of a function F(X) at a set of 
!    nodes XTAB, and estimates
!
!      Integral ( A <= X <= B ) F(X) DX
!
!    by computing the cubic natural spline S(X) that interpolates
!    F(X) at the nodes, and then computing
!
!      Integral ( A <= X <= B ) S(X) DX
!
!    exactly.
!
!    Other output from the program includes the definite integral
!    from X(1) to X(I) of S(X), and the coefficients necessary for
!    the user to evaluate the spline S(X) at any point.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), contains the points at which the
!    function was evaluated.  The XTAB's must be distinct and
!    in ascending order.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated values of
!    the function, FTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) A, lower limit of integration.
!
!    Input, real ( kind = 8 ) B, upper limit of integration.
!
!    Output, real ( kind = 8 ) Y(3,NTAB), will contain the coefficients
!    of the interpolating natural spline over each subinterval.
!    For XTAB(I) <= X <= XTAB(I+1),
!      S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I))
!                   + Y(2,I)*(X-XTAB(I))**2
!                   + Y(3,I)*(X-XTAB(I))**3
!
!    Output, real ( kind = 8 ) E(NTAB), E(I) = the definite integral from
!    XTAB(1) to XTAB(I) of S(X).
!
!    Workspace, real ( kind = 8 ) WORK(NTAB).
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e(ntab)
  real ( kind = 8 ) ftab(ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) s
  real ( kind = 8 ) term
  real ( kind = 8 ) u
  real ( kind = 8 ) work(ntab)
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) y(3,ntab)

  if ( ntab < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSPINT - Fatal error!'
    write ( *, '(a,i8)' ) '  NTAB must be at least 3, but input NTAB = ', ntab
    stop
  end if
 
  do i = 1, ntab-1
 
    if ( xtab(i+1) <= xtab(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CSPINT - Fatal error!'
      write ( *, '(a)' ) '  Nodes not in strict increasing order.'
      write ( *, '(a,i8)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
      write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
      stop
    end if
 
  end do
 
  s = 0.0D+00
  do i = 1, ntab-1
    r = ( ftab(i+1) - ftab(i) ) / ( xtab(i+1) - xtab(i) )
    y(2,i) = r - s
    s = r
  end do
 
  result = 0.0D+00
  s = 0.0D+00
  r = 0.0D+00
  y(2,1) = 0.0D+00
  y(2,ntab) = 0.0D+00
 
  do i = 2, ntab-1
    y(2,i) = y(2,i) + r * y(2,i-1)
    work(i) = 2.0D+00 * ( xtab(i-1) - xtab(i+1) ) - r * s
    s = xtab(i+1) - xtab(i)
    r = s / work(i)
  end do
 
  do j = 2, ntab-1
    i = ntab+1-j
    y(2,i) = ( ( xtab(i+1) - xtab(i) ) * y(2,i+1) - y(2,i) ) / work(i)
  end do
 
  do i = 1, ntab-1
    s = xtab(i+1) - xtab(i)
    r = y(2,i+1) - y(2,i)
    y(3,i) = r / s
    y(2,i) = 3.0D+00 * y(2,i)
    y(1,i) = ( ftab(i+1) - ftab(i) ) / s - ( y(2,i) + r ) * s
  end do
 
  e(1) = 0.0D+00
  do i = 1, ntab-1
    s = xtab(i+1)-xtab(i)
    term = ((( y(3,i) * 0.25D+00 * s + y(2,i) / 3.0D+00 ) * s &
      + y(1,i) * 0.5D+00 ) * s + ftab(i) ) * s
    e(i+1) = e(i) + term
  end do
!
!  Determine where the endpoints A and B lie in the mesh of XTAB's.
!
  r = a
  u = 1.0D+00
 
  do j = 1, 2
!
!  The endpoint is less than or equal to XTAB(1).
!
    if ( r <= xtab(1) ) then
      result = result - u * ( ( r - xtab(1) ) * y(1,1) * 0.5D+00 &
        + ftab(1) ) * ( r - xtab(1) )
!
!  The endpoint is greater than or equal to XTAB(NTAB).
!
    else if ( xtab(ntab) <= r ) then

      result = result -u * ( e(ntab) + ( r - xtab(ntab) ) &
        * ( ftab(ntab) + 0.5D+00 * ( ftab(ntab-1) &
        + ( xtab(ntab) - xtab(ntab-1) ) * y(1,ntab-1) ) &
        * ( r - xtab(ntab) )))
!
!  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
!
    else

      do i = 1, ntab-1
 
        if ( r <= xtab(i+1) ) then
          r = r - xtab(i)
          result = result - u * ( e(i) + ( ( ( &
              y(3,i) * 0.25D+00  * r &
            + y(2,i) / 3.0D+00 ) * r &
            + y(1,i) * 0.5D+00 ) * r + ftab(i) ) * r )
          go to 120
        end if
 
      end do
 
    end if
 
  120   continue
 
    u = -1.0D+00
    r = b
 
  end do
 
  return
end
subroutine cubint ( ntab, xtab, ftab, ia, ib, result, error )

!*****************************************************************************80
!
!! CUBINT approximates an integral using cubic interpolation of data.
!
!  Discussion:
!
!    The integral to be approximated is
! 
!      Integral ( XTAB(IB) <= X <= XTAB(IA) ) F(X) DX
!
!    The routine estimates the error in integration.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Philip Gill, GF Miller,
!    An algorithm for the integration of unequally spaced data,
!    The Computer Journal, 
!    Number 15, Number 1, 1972, pages 80-83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, integer ( kind = 4 ) IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer ( kind = 4 ) IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real ( kind = 8 ) ERROR, an estimate of the error in
!    integration.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) c
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) error
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) h3
  real ( kind = 8 ) h4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r4
  real ( kind = 8 ) result
  real ( kind = 8 ) s
  real ( kind = 8 ) term
  real ( kind = 8 ) xtab(ntab)

  result = 0.0D+00
  error = 0.0D+00
 
  if ( ia == ib ) then
    return
  end if
 
  if ( ntab < 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  NTAB must be at least 4, but input NTAB = ', ntab
    stop
  end if
 
  if ( ia < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IA must be at least 1, but input IA = ', ia
    stop
  end if
 
  if ( ntab < ia ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IA must be <= NTAB, but input IA = ', ia
    stop
  end if
 
  if ( ib < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IB must be at least 1, but input IB = ', ib
    stop
  end if
 
  if ( ntab < ib ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IB must be <= NTAB, but input IB = ', ib
    stop
  end if
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
  if ( ib < ia ) then
    ind = -1
    it = ib
    ib = ia
    ia = it
  else
    ind = 1
  end if
 
  s = 0.0D+00
  c = 0.0D+00
  r4 = 0.0D+00
  j = ntab-2
  if ( ia < ntab-1 .or. ntab == 4 ) then
    j = max ( 3, ia )
  end if

  k = 4
  if ( 2 < ib .or. ntab == 4 ) then
    k = min ( ntab, ib + 2 ) - 1
  end if
 
  do i = j, k
 
    if ( i <= j ) then
 
      h2 = xtab(j-1) - xtab(j-2)
      d3 = ( ftab(j-1) - ftab(j-2) ) / h2
      h3 = xtab(j) - xtab(j-1)
      d1 = ( ftab(j) - ftab(j-1) ) / h3
      h1 = h2 + h3
      d2 = ( d1 - d3 ) / h1
      h4 = xtab(j+1) - xtab(j)
      r1 = ( ftab(j+1) - ftab(j) ) / h4
      r2 = ( r1 - d1 ) / ( h4 + h3 )
      h1 = h1 + h4
      r3 = (r2-d2) / h1
 
      if ( ia <= 1 ) then
        result = h2 * ( ftab(1) + h2 * ( 0.5D+00 * d3 - h2 &
          * ( d2 / 6.0D+00 -(h2+h3+h3)*r3/12.0D+00)))
        s = -h2**3 * (h2*(3.0D+00*h2+5.0D+00*h4)+10.0D+00*h3*h1) / 60.0D+00
      end if
 
    else
 
      h4 = xtab(i+1) - xtab(i)
      r1 = ( ftab(i+1) - ftab(i) ) / h4
      r4 = h4 + h3
      r2 = ( r1 - d1 ) / r4
      r4 = r4 + h2
      r3 = ( r2 - d2 ) / r4
      r4 = ( r3 - d3 ) / ( r4 + h1 )
 
    end if
 
    if ( ia < i .and. i <= ib ) then
 
      term = h3 * ( ( ftab(i) + ftab(i-1) ) * 0.5D+00 &
        -h3 * h3 * ( d2 + r2 + ( h2 - h4 ) * r3 ) / 12.0D+00 )
      result = result + term
      c = h3**3 * ( 2.0D+00 * h3 * h3 &
        + 5.0D+00 * ( h3 * ( h4 + h2 ) + 2.0D+00 * h2 * h4 ) ) / 120.0D+00
      error = error + (c+s)*r4
 
      if ( i /= j ) then
        s = c
      else
        s = s + c + c
      end if
 
    else
 
      error = error + r4 * s
 
    end if
 
    if ( k <= i ) then
 
      if ( ntab <= ib ) then
        term = h4 * ( ftab(ntab) - h4 * ( 0.5 * r1 &
          + h4 * ( r2 / 6.0D+00 + ( h3 + h3 + h4 ) * r3 / 12.0D+00 )))
        result = result + term
        error = error - h4**3 * r4 * &
          ( h4 * ( 3.0D+00 * h4 + 5.0D+00 * h2 ) &
          + 10.0D+00 * h3 * ( h2 + h3 + h4 ) ) / 60.0D+00
      end if
 
      if ( ntab-1 <= ib ) then
        error = error + s * r4
      end if

    else

      h1 = h2
      h2 = h3
      h3 = h4
      d1 = r1
      d2 = r2
      d3 = r3
    end if
 
  end do
!
!  Restore original values of IA and IB, reverse signs
!  of RESULT and ERROR, to account for integration
!  that proceeded from high X to low X.
!
  if ( ind /= 1 ) then
    it = ib
    ib = ia
    ia = it
    result = -result
    error = -error
  end if
 
  return
end
subroutine filon_cos ( ntab, ftab, a, b, t, result )

!*****************************************************************************80
!
!! FILON_COS uses Filon's method on integrals with a cosine factor.
!
!  Discussion:
!
!    The integral to be approximated has the form:
!
!      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
!
!    where T is user specified.
!
!    The function is interpolated over each subinterval by
!    a parabolic arc.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Chase, Lloyd Fosdick,
!    An Algorithm for Filon Quadrature,
!    Communications of the Association for Computing Machinery,
!    Volume 12, Number 8, August 1969, pages 453-457.
!
!    Stephen Chase, Lloyd Fosdick,
!    Algorithm 353:
!    Filon Quadrature,
!    Communications of the Association for Computing Machinery,
!    Volume 12, Number 8, August 1969, pages 457-458.
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of data points.
!    NTAB must be odd, and greater than 1.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the value of the function
!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) T, the multiplier of the X argument of the cosine.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c2n
  real ( kind = 8 ) c2nm1
  real ( kind = 8 ) cost
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ) gamma
  real ( kind = 8 ) h
  real ( kind = 8 ) result
  real ( kind = 8 ) sint
  real ( kind = 8 ) t
  real ( kind = 8 ) theta
  real ( kind = 8 ) xtab(ntab)

  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  if ( ntab <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILON_COS - Fatal error!'
    write ( *, '(a)' ) '  NTAB < 2'
    write ( *, '(a,i8)' ) '  NTAB = ', ntab
    stop
  end if
 
  if ( mod ( ntab, 2 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILON_COS - Fatal error!'
    write ( *, '(a)' ) '  NTAB must be odd.'
    write ( *, '(a,i8)' ) '  NTAB = ', ntab
    stop
  end if
!
!  Set up a vector of the NTAB X values.
! 
  call r8vec_even ( ntab, a, b, xtab )

  h = ( b - a ) / real ( ntab - 1, kind = 8 )

  theta = t * h
  sint = sin ( theta )
  cost = cos ( theta )

  if ( 6.0D+00 * abs ( theta ) <= 1.0D+00 ) then

    alpha = 2.0D+00 * theta**3 /   45.0D+00 &
          - 2.0D+00 * theta**5 /  315.0D+00 &
          + 2.0D+00 * theta**7 / 4725.0D+00
  
    beta =  2.0D+00            /     3.0D+00 &
          + 2.0D+00 * theta**2 /    15.0D+00 &
          - 4.0D+00 * theta**4 /   105.0D+00 &
          + 2.0D+00 * theta**6 /   567.0D+00 &
          - 4.0D+00 * theta**8 / 22275.0D+00

    gamma = 4.0D+00            /      3.0D+00 &
          - 2.0D+00 * theta**2 /     15.0D+00 &
          +           theta**4 /    210.0D+00 &
          -           theta**6 /  11340.0D+00

  else

    alpha = ( theta**2 + theta * sint * cost &
      - 2.0D+00 * sint**2 ) / theta**3

    beta = ( 2.0D+00 * theta + 2.0D+00 * theta * cost**2 &
      - 4.0D+00 * sint * cost ) / theta**3

    gamma = 4.0D+00 * ( sint - theta * cost ) / theta**3
  
  end if

  c2n = sum ( ftab(1:ntab:2) * cos ( t * xtab(1:ntab:2) ) ) &
    - 0.5D+00 * ( ftab(ntab) * cos ( t * xtab(ntab) ) &
                + ftab(1) * cos ( t * xtab(1) ) )

  c2nm1 = sum ( ftab(2:ntab-1:2) * cos ( t * xtab(2:ntab-1:2) ) )
 
  result = h * ( &
      alpha * ( ftab(ntab) * sin ( t * xtab(ntab) ) & 
              - ftab(1)    * sin ( t * xtab(1) ) ) &
    + beta * c2n &
    + gamma * c2nm1 )

  return
end
subroutine filon_sin ( ntab, ftab, a, b, t, result )

!*****************************************************************************80
!
!! FILON_SIN uses Filon's method on integrals with a sine factor.
!
!  Discussion:
!
!    The integral to be approximated has the form
!
!      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
!
!    where T is user specified.
!
!    The function is interpolated over each subinterval by
!    a parabolic arc.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Chase, Lloyd Fosdick,
!    An Algorithm for Filon Quadrature,
!    Communications of the Association for Computing Machinery,
!    Volume 12, Number 8, August 1969, pages 453-457.
!
!    Stephen Chase, Lloyd Fosdick,
!    Algorithm 353:
!    Filon Quadrature,
!    Communications of the Association for Computing Machinery,
!    Volume 12, Number 8, August 1969, pages 457-458.
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of data points, 
!    including the endpoints.  NTAB must be odd, and greater than 1.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the value of the function
!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) T, multiplier of the X argument of the sine.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) cost
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ) gamma
  real ( kind = 8 ) h
  real ( kind = 8 ) result
  real ( kind = 8 ) s2n
  real ( kind = 8 ) s2nm1
  real ( kind = 8 ) sint
  real ( kind = 8 ) t
  real ( kind = 8 ) theta
  real ( kind = 8 ) xtab(ntab)

  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  if ( ntab <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
    write ( *, '(a)' ) '  NTAB < 2'
    write ( *, '(a,i8)' ) '  NTAB = ',ntab
    stop
  end if
 
  if ( mod ( ntab, 2 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
    write ( *, '(a)' ) '  NTAB must be odd.'
    write ( *, '(a,i8)' ) '  NTAB = ',ntab
    stop
  end if
!
!  Set up a vector of the NTAB X values.
! 
  call r8vec_even ( ntab, a, b, xtab )

  h = ( b - a ) / real ( ntab - 1, kind = 8 )
  theta = t * h

  sint = sin ( theta )
  cost = cos ( theta )

  if ( 6.0D+00 * abs ( theta ) <= 1.0D+00 ) then

    alpha = 2.0D+00 * theta**3 /   45.0D+00 &
          - 2.0D+00 * theta**5 /  315.0D+00 &
          + 2.0D+00 * theta**7 / 4725.0D+00
  
    beta =  2.0D+00            /     3.0D+00 &
          + 2.0D+00 * theta**2 /    15.0D+00 &
          - 4.0D+00 * theta**4 /   105.0D+00 &
          + 2.0D+00 * theta**6 /   567.0D+00 &
          - 4.0D+00 * theta**8 / 22275.0D+00

    gamma = 4.0D+00            /      3.0D+00 &
          - 2.0D+00 * theta**2 /     15.0D+00 &
          +           theta**4 /    210.0D+00 &
          -           theta**6 /  11340.0D+00

  else
 
    alpha = ( theta**2 + theta * sint * cost &
      - 2.0D+00 * sint**2 ) / theta**3

    beta = ( 2.0D+00 * theta + 2.0D+00 * theta * cost**2 &
      - 4.0D+00 * sint * cost ) / theta**3

    gamma = 4.0D+00 * ( sint - theta * cost ) / theta**3
 
  end if
  
  s2n = sum ( ftab(1:ntab:2) * sin ( t * xtab(1:ntab:2) ) ) &
    - 0.5D+00 * ( ftab(ntab) * sin ( t * xtab(ntab) ) &
                + ftab(1) * sin ( t * xtab(1) ) )

  s2nm1 = sum ( ftab(2:ntab-1:2) * sin ( t * xtab(2:ntab-1:2) ) )

  result = h * ( &
      alpha * ( ftab(1) * cos ( t * xtab(1) ) &
              - ftab(ntab) * cos ( t * xtab(ntab) ) ) &
    + beta * s2n &
    + gamma * s2nm1 )
 
  return
end
function gamma ( x )

!*****************************************************************************80
!
!! GAMMA calculates the Gamma function for a real argument X.
!
!  Definition:
!
!    GAMMA(X) = Integral ( 0 <= T <= Infinity ) T**(X-1) EXP(-T) DT
!
!  Recursion:
!
!    GAMMA(X+1) = X * GAMMA(X)
!
!  Special values:
!
!    GAMMA(0.5) = SQRT(PI)
!    If N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the GAMMA
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for X .GE. 12 are from reference 2.
!    The accuracy achieved depends on the arithmetic system, the
!    compiler, the intrinsic functions, and proper selection of the
!    machine-dependent constants.
!
!  Machine-dependent constants:
!
!    BETA: radix for the floating-point representation.
!    MAXEXP: the smallest positive power of BETA that overflows.
!    XBIG: the largest argument for which GAMMA(X) is representable
!      in the machine, i.e., the solution to the equation
!      GAMMA(XBIG) = BETA**MAXEXP.
!    XMININ: the smallest positive floating-point number such that
!      1/XMININ is machine representable.
!
!    Approximate values for some important machines are:
!
!                               BETA       MAXEXP        XBIG
!
!    CRAY-1         (S.P.)        2         8191        966.961
!    Cyber 180/855
!      under NOS    (S.P.)        2         1070        177.803
!    IEEE (IBM/XT,
!      SUN, etc.)   (S.P.)        2          128        35.040
!    IEEE (IBM/XT,
!      SUN, etc.)   (D.P.)        2         1024        171.624
!    IBM 3033       (D.P.)       16           63        57.574
!    VAX D-Format   (D.P.)        2          127        34.844
!    VAX G-Format   (D.P.)        2         1023        171.489
!
!                               XMININ
!
!    CRAY-1         (S.P.)   1.84D-2466
!    Cyber 180/855
!      under NOS    (S.P.)   3.14D-294
!    IEEE (IBM/XT,
!      SUN, etc.)   (S.P.)   1.18D-38
!    IEEE (IBM/XT,
!      SUN, etc.)   (D.P.)   2.23D-308
!    IBM 3033       (D.P.)   1.39D-76
!    VAX D-Format   (D.P.)   5.88D-39
!    VAX G-Format   (D.P.)   1.12D-308
!
!  Author:
!
!    William Cody and L. Stoltz,
!    Applied Mathematics Division,
!    Argonne National Laboratory,
!    Argonne, Illinois, 60439.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics, 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) GAMMA, the value of the function.  The program
!    returns the a huge value for singularities or when overflow would occur.
!    The computation is believed to be free of underflow and overflow.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) fact
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: PI = &
    3.14159265358979323846264338327950288419716939937510D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353e+01, &
     3.15350626979604161529144e+02, &
    -1.01515636749021914166146e+03, &
    -3.10777167157231109440444e+03, &
     2.25381184209801510330112e+04, &
     4.75584627752788110767815e+03, &
    -1.34659959864969306392456e+05, &
    -1.15132259675553483497211e+05 /)
  real ( kind = 8 ), parameter :: SQRTPI = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum1
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: XBIG = 35.040D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: XMININ = 1.18D-38
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    gamma = y - y1

    if ( gamma /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - PI / sin ( PI * gamma )
      y = y + 1.0D+00

    else

      gamma = huge ( gamma )
      return

    end if

  end if
!
!  Argument < EPS
!
  if ( y < epsilon ( y ) ) then

    if ( XMININ <= y ) then
      gamma = 1.0D+00 / y
    else
      gamma = huge ( gamma )
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0D+00 < argument < 1.0D+00
!
    if ( y < 1.0D+00 ) then
      z = y
      y = y + 1.0D+00
!
!  1.0D+00 < argument < 12.0, reduce argument if necessary.
!
    else
      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00
    end if
!
!  Evaluate approximation for 1.0D+00 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    gamma = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0D+00 < argument < 1.0.
!
    if ( y1 < y ) then
      gamma = gamma / y1
!
!  Adjust result for case  2.0D+00 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        gamma = gamma * y
        y = y + 1.0D+00
      end do

    end if
!
!  Evaluate for 12 <= argument.
!
  else

    if ( y <= XBIG ) then

      ysq = y**2
      sum1 = c(7)
      do i = 1, 6
        sum1 = sum1 / ysq + c(i)
      end do
      sum1 = sum1 / y - y + SQRTPI
      sum1 = sum1 + ( y - 0.5D+00 ) * log ( y )
      gamma = exp ( sum1 )

    else

      gamma = huge ( gamma )
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    gamma = - gamma
  end if

  if ( fact /= 1.0D+00 ) then
    gamma = fact / gamma
  end if

  return
end
subroutine gaus8 ( func, a, b, err, result, ierr )

!*****************************************************************************80
!
!! GAUS8 estimates the integral of a function.
!
!  Discussion:
!
!    GAUS8 integrates real functions of one variable over finite
!    intervals using an adaptive 8-point Legendre-Gauss
!    algorithm.
!
!    GAUS8 is intended primarily for high accuracy integration or
!    integration of smooth functions.
!
!  Modified:
!
!    30 October 2000
!
!  Author:
!
!    Ron Jones,
!    Sandia National Laboratory,
!    Los Alamos, New Mexico
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, name of external function to
!    be integrated.  This name must be in an external statement in the
!    calling program.  FUNC must be a function of one real argument.  The value
!    of the argument to FUNC is the variable of integration
!    which ranges from A to B.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input/output, real ( kind = 8 ) ERR.
!    On input, ERR is a requested pseudorelative error
!    tolerance.  Normally pick a value of ABS ( ERR ) so that
!    STOL < ABS ( ERR ) <= 1.0D-3 where STOL is the single precision
!    unit roundoff.
!    RESULT will normally have no more error than
!    ABS ( ERR ) times the integral of the absolute value of
!    FUN(X).  Usually, smaller values for ERR yield more
!    accuracy and require more function evaluations.
!    A negative value for ERR causes an estimate of the
!    absolute error in RESULT to be returned in ERR.  Note that
!    ERR must be a variable (not a constant) in this case.
!    Note also that the user must reset the value of ERR
!    before making any more calls that use the variable ERR.
!    On output, ERR will be an estimate of the absolute error
!    in RESULT if the input value of ERR was negative.  ERR is
!    unchanged if the input value of ERR was non-negative.
!    The estimated error is solely for information to the user
!    and should not be used as a correction to the computed integral.
!
!    Output, real ( kind = 8 ) RESULT, the computed value of the integral.
!
!    Output, integer ( kind = 4 ) IERR, a status code.
!    Normal Codes:
!     1 RESULT most likely meets requested error tolerance, or A = B.
!    -1 A and B are too nearly equal to allow normal
!        integration.  RESULT is set to zero.
!     Abnormal Code:
!     2 RESULT probably does not meet requested error tolerance.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) aa(30)
  real ( kind = 8 ) ae
  real ( kind = 8 ) anib
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) ce
  real ( kind = 8 ) ee
  real ( kind = 8 ) ef
  real ( kind = 8 ) eps
  real ( kind = 8 ) err
  real ( kind = 8 ) est
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) g8
  real ( kind = 8 ) gl
  real ( kind = 8 ) glr
  real ( kind = 8 ) gr(30)
  real ( kind = 8 ) h
  real ( kind = 8 ) hh(30)
  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: kml = 6
  integer ( kind = 4 ), save :: kmx = 5000
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmn
  integer ( kind = 4 ) lmx
  integer ( kind = 4 ) lr(30)
  integer ( kind = 4 ) mxl
  integer ( kind = 4 ) nbits
  integer ( kind = 4 ) nib
  integer ( kind = 4 ), save :: nlmn = 1
  integer ( kind = 4 ) nlmx
  real ( kind = 8 ) result
  real ( kind = 8 ) tol
  real ( kind = 8 ) vl(30)
  real ( kind = 8 ) vr
  real ( kind = 8 ), save :: w1 = 3.62683783378361983D-01
  real ( kind = 8 ), save :: w2 = 3.13706645877887287D-01
  real ( kind = 8 ), save :: w3 = 2.22381034453374471D-01
  real ( kind = 8 ), save :: w4 = 1.01228536290376259D-01
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: x1 = 1.83434642495649805D-01
  real ( kind = 8 ), save :: x2 = 5.25532409916328986D-01
  real ( kind = 8 ), save :: x3 = 7.96666477413626740D-01
  real ( kind = 8 ), save :: x4 = 9.60289856497536232D-01
!
!  Warning!  Statement function!
!
  g8(x,h) = h * ( ( & 
                w1 * ( func ( x - x1 * h ) + func ( x + x1 * h ) )   &
              + w2 * ( func ( x - x2 * h ) + func ( x + x2 * h ) ) ) &
            + ( w3 * ( func ( x - x3 * h ) + func ( x + x3 * h ) )   &
              + w4 * ( func ( x - x4 * h ) + func ( x + x4 * h ) ) ) )

  if ( a == b ) then
    err = 0.0D+00
    result = 0.0D+00
    return
  end if
 
  if ( icall /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAUS8 - Fatal error!'
    write ( *, '(a)' ) '  GAUS8 was called recursively.'
    stop
  end if

  icall = 1
!
!  DIGITS ( X ) = number of base 2 digits in representation of X.
!
  k = digits ( result )

  anib = log10 ( 2.0D+00 ) * real ( k, kind = 8 ) / 0.30102000D+00
  nbits = int ( anib )
  nlmx = min ( 30, ( nbits * 5 ) / 8 )
  result = 0.0D+00
  ierr = 1
  ce = 0.0D+00
  result = 0.0D+00
  lmx = nlmx
  lmn = nlmn
 
  if ( b /= 0.0D+00 ) then

    if ( sign ( 1.0D+00, b ) * a <= 0.0D+00 ) then
      go to 10
    end if

    c = abs ( 1.0D+00 - a / b )
    if ( 0.1D+00 < c ) then
      go to 10
    end if

    if ( c <= 0.0D+00 ) then
      go to 140
    end if

    anib = 0.5D+00 - log ( c ) / 0.69314718D+00
    nib = int ( anib )
    lmx = min ( nlmx, nbits-nib-7 )

    if ( lmx < 1 ) then
      go to 130
    end if

    lmn = min ( lmn, lmx )

  end if
 
10    continue
 
  tol = max ( abs ( err ), 2.0D+00**(5-nbits)) / 2.0D+00
  if ( err == 0.0D+00 ) then
    tol = sqrt ( epsilon ( 1.0D+00 ) )
  end if

  eps = tol
  hh(1) = ( b - a ) / 4.0D+00
  aa(1) = a
  lr(1) = 1
  l = 1
  est = g8 ( aa(l) + 2.0D+00 * hh(l), 2.0D+00 * hh(l) )
  k = 8
  area = abs ( est )
  ef = 0.5D+00
  mxl = 0
!
!  Compute refined estimates, estimate the error, etc.
!
20 continue
 
  gl = g8 ( aa(l) + hh(l), hh(l) )
  gr(l) = g8 ( aa(l) + 3.0D+00 * hh(l), hh(l) )
  k = k + 16
  area = area + ( abs ( gl ) + abs ( gr(l) ) - abs ( est ) )
 
  glr = gl + gr(l)
  ee = abs ( est - glr ) * ef
  ae = max ( eps * area, tol * abs ( glr ) )

  if ( ee - ae <= 0.0D+00 ) then
    go to 40
  else
    go to 50
  end if
 
30 continue
 
  mxl = 1
 
40 continue
 
  ce = ce + ( est - glr )
 
  if ( lr(l) <= 0 ) then
    go to 60
  else
    go to 80
  end if
!
!  Consider the left half of this level
!
50 continue

  if ( kmx < k ) then
    lmx = kml
  end if

  if ( lmx <= l ) then
    go to 30
  end if

  l = l + 1
  eps = eps * 0.5D+00
  ef = ef / sqrt ( 2.0D+00 )
  hh(l) = hh(l-1) * 0.5D+00
  lr(l) = -1
  aa(l) = aa(l-1)
  est = gl
  go to 20
!
!  Proceed to right half at this level
!
60 continue

  vl(l) = glr

70 continue

  est = gr(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 4.0D+00 * hh(l)
  go to 20
!
!  Return one level
!
80 continue

  vr = glr

90 continue

  if ( l <= 1 ) then
    go to 120
  end if

  l = l - 1
  eps = eps * 2.0D+00
  ef = ef * sqrt ( 2.0D+00 )
 
  if ( lr(l) <= 0 ) then
    vl(l) = vl(l+1) + vr
    go to 70
  else
    vr = vl(l+1) + vr
    go to 90
  end if
!
!  Exit
!
120   continue
 
  result = vr

  if ( mxl == 0 .or. abs ( ce ) <= 2.0D+00 * tol * area ) then
    go to 140
  end if

  ierr = 2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAUS8 - Warning!'
  write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
  icall = 0

  if ( err < 0.0D+00 ) then
    err = ce
  end if

  return
 
130   continue
 
  ierr = -1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAUS8 - Warning!'
  write ( *, '(a)' ) '  A and B are too close to carry out integration.'
 
140   continue
 
  icall = 0
 
  if ( err < 0.0D+00 ) then
    err = ce
  end if
 
  return
end
subroutine gausq2 ( n, d, e, z, ierr )

!*****************************************************************************80
!
!! GAUSQ2 finds the eigenvalues of a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    GAUSQ2 finds the eigenvalues and first components of the
!    eigenvectors of a symmetric tridiagonal matrix by the implicit QL
!    method.
!
!    GAUSQ2 is a translation of an ALGOL procedure as modified by
!    Dubrulle.
!
!    GAUSQ2 is a modified version of the EISPACK routine IMTQL2.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!    Augustin Dubrulle,
!    A short note on the implicit QL algorithm for symmetric
!    tridiagonal matrices,
!    Numerische Mathematik,
!    Volume 15, Number 5, September 1970, page 450.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N).
!    On input, D contains the diagonal elements of the matrix.
!    On output, D contains the eigenvalues in ascending order.
!    If an error exit is made, the eigenvalues are correct but
!    unordered for indices 1, 2, ..., IERR-1;
!
!    Input/output, real ( kind = 8 ) E(N).
!    On input, E contains the subdiagonal elements of the input matrix
!    in its first N-1 positions.  E(N) is arbitrary.
!    On output, E has been destroyed.
!
!    Input/output, real ( kind = 8 ) Z(N).
!    On input, Z contains the first row of the identity matrix.
!    On output, Z contains the first components of the orthonormal
!    eigenvectors of the symmetric tridiagonal matrix.  If an error exit is
!    made, Z contains the eigenvectors associated with the stored
!    eigenvalues.
!
!    Output, integer ( kind = 4 ) IERR.
!    0, for normal return,
!    J, if the j-th eigenvalue has not been determined after 30 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) epmach
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) temp
  real ( kind = 8 ) z(n)

  epmach = epsilon ( epmach )
 
  ierr = 0

  if ( n == 1 ) then
    return
  end if
 
  e(n) = 0.0D+00
 
  do l = 1, n
 
    j = 0
!
!  Look for a small sub-diagonal element
!
    do m = l, n

      if ( m == n ) then
        exit
      end if

      if ( abs ( e(m) ) <= epmach * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
        exit
      end if

    end do
 
10  continue
 
    p = d(l)

    if ( m == l ) then
      go to 20
    end if
 
    if ( j == 30 ) then
      ierr = l
      return
    end if
 
    j = j + 1
!
!  Form shift
!
    g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
    r = sqrt ( g * g + 1.0D+00 )
    g = d(m) - p + e(l) / ( g + sign ( r, g ) )
    s = 1.0D+00
    c = 1.0D+00
    p = 0.0D+00
    mml = m - l
 
    do ii = 1, mml
 
      i = m - ii
      f = s * e(i)
      b = c * e(i)
 
      if ( abs ( g ) <= abs ( f ) ) then
 
        c = g / f
        r = sqrt ( c * c + 1.0D+00 )
        e(i+1) = f * r
        s = 1.0D+00 / r
        c = c * s
 
      else
 
        s = f / g
        r = sqrt ( s * s + 1.0D+00 )
        e(i+1) = g * r
        c = 1.0D+00 / r
        s = s * c
 
      end if
 
      g = d(i+1) - p
      r = ( d(i) - g ) * s + 2.0D+00 * c * b
      p = s * r
      d(i+1) = g + p
      g = c * r - b
!
!  Form the first component of the vector.
!
      f = z(i+1)
      z(i+1) = s * z(i) + c * f
      z(i) = f * z(i) - s * f
    end do
 
    d(l) = d(l) - p
    e(l) = g
    e(m) = 0.0D+00
    go to 10
 
20  continue
 
  end do
!
!  Order the eigenvalues and eigenvectors.
!
  do ii = 2, n
 
    i = ii - 1
    k = i
    p = d(i)
 
    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do
 
    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p

      temp = z(i)
      z(i) = z(k)
      z(k) = temp

    end if
 
  end do
 
  return
end
subroutine gaussq ( kind, norder, alpha, beta, kpts, endpts, b, xtab, &
  weight )

!*****************************************************************************80
!
!! GAUSSQ computes a Gauss quadrature rule.
!
!  Discussion:
!
!    GAUSSQ computes the nodes and weights for Gaussian-type quadrature
!    rules with pre-assigned nodes.
!
!    These are used when one wishes to approximate
!
!      Integral ( A <= X <= B )  F(X) W(X) DX
!
!    by
!
!      Sum ( 1 <= J <= NORDER ) WEIGHT(I) * F(XTAB(I))
!
!
!    GAUSSQ includes six integration rules that are applicable
!    to this problem, for particular weight functions and particular
!    intervals, including infinite and semi-infinite intervals.
!
!    Associated with each weight function W(X) is a set of
!    orthogonal polynomials.  The nodes XTAB are just the zeroes
!    of the proper NORDER-th degree polynomial.
!
!    GAUSSQ allows the user to modify the rule to require that
!    one or both of the endpoints of the interval are to be
!    included as quadrature nodes.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Gene Golub, John Welsch,
!    Calculation of Gaussian Quadrature Rules,
!    Mathematics of Computation,
!    Volume 23, Number 106, April 1969, pages 221-230.
!
!    Gene Golub,
!    Some Modified Matrix Eigenvalue Problems,
!    SIAM Review,
!    Volume 15, Number 2, Part 1, April 1973, pages 318-334.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice-Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, chooses the rule to be calculated.
!    1:  Legendre quadrature,
!         W(X) = 1
!         on (-1, 1)
!    2:  Chebyshev quadrature of the first kind
!         W(X) = 1/SQRT(1 - X*X)
!         on (-1, +1)
!    3:  Chebyshev quadrature of the second kind
!         W(X) = SQRT(1 - X*X)
!         on (-1, 1)
!    4:  Hermite quadrature,
!         W(X) = EXP(-X*X)
!         on (-infinity, +infinity)
!    5:  Jacobi quadrature,
!         W(X) = (1-X)**ALPHA * (1+X)**BETA
!         on (-1, 1),
!         -1 < ALPHA, -1 < BETA.
!         Note that KIND = 2 and 3 are a special case of this.
!    6:  Generalized Laguerre quadrature,
!         W(X) = EXP(-X)*X**ALPHA
!         on (0, +infinity),
!         -1 < ALPHA
!
!    Input, integer ( kind = 4 ) NORDER, the number of points used for the quadrature rule.
!
!    Input, real ( kind = 8 ) ALPHA, is only required for Gauss-Jacobi and
!    Gauss-Laguerre quadrature.  Its value is ignored in other cases.
!
!    Input, real ( kind = 8 ) BETA, is only required for Gauss-Jacobi
!    quadrature.  Its value is ignored in other cases.
!
!    Input, integer ( kind = 4 ) KPTS, is normally zero.
!    If KPTS is nonzero, it signals that one or both of the
!    endpoints of the interval is required to be a node.
!    This is called Gauss-Radau or Gauss-Lobatto quadrature.
!    Then KPTS is the number of endpoints that must be
!    included, either 1 or 2.
!
!    Input, real ( kind = 8 ) ENDPTS(2).
!    If KPTS is 1 or 2, ENDPTS contains the locations of the
!    endpoints to be fixed.
!
!    Workspace, real ( kind = 8 ) B(NORDER).
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the nodes for the quadrature rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights for the 
!    quadrature rule.
!
  implicit none

  integer ( kind = 4 ) norder

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(norder)
  real ( kind = 8 ) beta
  real ( kind = 8 ) endpts(2)
  real ( kind = 8 ) gam
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) kpts
  real ( kind = 8 ) muzero
  real ( kind = 8 ) solve
  real ( kind = 8 ) t1
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
!  Get the diagonal coefficients XTAB(1:NORDER) and off-diagonal
!  coefficients B(1:NORDER-1) and MUZERO.
!
  call class ( kind, norder, alpha, beta, b, xtab, muzero )
!
!  The matrix of coefficients is assumed to be symmetric.
!  The array XTAB contains the diagonal elements, the array
!  B the off-diagonal elements.
!  Make appropriate changes in the lower right 2 by 2 submatrix.
!
!  If KPTS = 1, only XTAB(NORDER) must be changed.
!
  if ( kpts == 1 ) then
 
    xtab(norder) = endpts(1) &
      + solve ( endpts(1), norder, xtab, b ) * b(norder-1)**2
!
!  If KPTS = 2, XTAB(NORDER) and B(NORDER-1) must be recomputed.
!
  else if ( kpts == 2 ) then
 
    gam = solve ( endpts(1), norder, xtab, b )
    t1 = ( ( endpts(1) - endpts(2) ) &
      / ( solve ( endpts(2), norder, xtab, b ) - gam ) )
    b(norder-1) = sqrt ( t1 )
    xtab(norder) = endpts(1) + gam * t1
 
  end if
!
!  The indices of the elements of B run from 1 to NORDER-1.
!  The value of B(NORDER) is of no importance.
!
!  Now compute the eigenvalues of the symmetric tridiagonal
!  matrix, which has been modified as necessary.
!
!  The method used is a QL-type method with origin shifting.
!
  weight(1) = 1.0D+00
  weight(2:norder) = 0.0D+00
 
  call gausq2 ( norder, xtab, b, weight, ierr )
 
  do i = 1, norder
    weight(i) = muzero * weight(i)**2
  end do
 
  return
end
subroutine hiordq ( ntab, delt, y, work, result )

!*****************************************************************************80
!
!! HIORDQ approximates the integral of a function using equally spaced data.
!
!  Discussion:
!
!    The method applies the trapezoidal rule to various subsets of the
!    data, and then applies Richardson extrapolation.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    Alan Kaylor Cline,
!    Department of Computer Science,
!    University of Texas at Austin.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, number of data points.
!
!    Input, real ( kind = 8 ) DELT, the spacing between the X values of the
!    data.  The actual X values are not needed!
!
!    Input, real ( kind = 8 ) Y(NTAB), the Y values of the data.
!
!    Work array, real ( kind = 8 ) WORK(2*(NTAB-1)).  The actual minimum amount
!    of workspace required is two times the number of integer
!    divisors of NTAB-1.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) delt
  real ( kind = 8 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jbak
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  real ( kind = 8 ) result
  real ( kind = 8 ) sum2
  real ( kind = 8 ) sum1
  real ( kind = 8 ) work(2*(ntab-1))
  real ( kind = 8 ) y(ntab)
!
!  Determine initial trapezoidal rule
!
  sum1 = ( y(1) + y(ntab) ) / 2.0D+00
  j = -1
 
  do k = 1, ntab-1
!
!  Check if K divides NTAB-1
!
    if ( ( ( ntab - 1 ) / k ) * k == ntab - 1 ) then
!
!  Determine the K-point trapezoidal rule.
!
      sum2 = -sum1
      do i = 1, ntab, (ntab-1)/k
        sum2 = sum2 + y(i)
      end do
 
      j = j + 2
      work(j) = delt * sum2 * real ( ( ntab - 1 ) / k, kind = 8 )
      work(j+1) = real ( ( ( ntab - 1 ) / k )**2, kind = 8 )
!
!  Apply Richardson extrapolation.
!
      if ( k /= 1 ) then
 
        do jj = 3, j, 2
          jbak = j+1-jj
          fac = work(j+1) / ( work(j+1) - work(jbak+1) )
          work(jbak) = work(jbak+2) + fac * ( work(jbak) - work(jbak+2) )
        end do
 
      end if
 
    end if
 
  end do
 
  result = work(1)
 
  return
end
subroutine iratex ( func, a, b, epsin, epsout, result, ind )

!*****************************************************************************80
!
!! IRATEX estimates the integral of a function.
!
!  Discussion:
!
!    IRATEX estimates the integral from A to B of F(X), using the
!    trapezoidal rule for several stepsizes H.
!
!    Then a rational function of H*H is fitted to these results, and 
!    an estimate of the integral is computed by extrapolating the 
!    results to a stepsize of zero.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Roland Bulirsch, Josef Stoer,
!    Fehlerabschaetzungen und Extrapolation mit rationaled Funktionen
!    bei Verfahren vom Richardson-Typus,
!    (Error estimates and extrapolation with rational functions
!    in processes of Richardson type),
!    Numerische Mathematik,
!    Volume 6, Number 1, December 1964, pages 413-427.
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC.
!    FUNC is the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      FUNCTION FUNC ( X ) 
!    which evaluates the function at the point X.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, real ( kind = 8 ) EPSIN, the requested relative error tolerance.
!
!    Output, real ( kind = 8 ) EPSOUT, an estimate of the error in the
!    integration.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!    Output, integer ( kind = 4 ) IND, error return flag.
!    IND = 0 if the requested accuracy was not achieved,
!    IND = 1 if the accuracy was achieved.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) ba
  logical bo
  logical bu
  real ( kind = 8 ) c
  real ( kind = 8 ) d(6)
  real ( kind = 8 ) d1
  real ( kind = 8 ) ddt
  real ( kind = 8 ) den
  real ( kind = 8 ) dt(7)
  real ( kind = 8 ) e
  real ( kind = 8 ) ent
  real ( kind = 8 ) epsin
  real ( kind = 8 ) epsout
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) gr
  real ( kind = 8 ) hm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  logical odd
  real ( kind = 8 ) rnderr
  real ( kind = 8 ) result
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t2a
  real ( kind = 8 ) ta
  real ( kind = 8 ) tab
  real ( kind = 8 ) tb
  real ( kind = 8 ) tnt
  real ( kind = 8 ) v
  real ( kind = 8 ) w

  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  rnderr = epsilon ( 1.0D+00 )
  epsin = max ( epsin, 8.0D+00 * rnderr )
  ind = 0
  n = 2
  np1 = 3
  ba = b - a
  t1 = 0.0D+00
  gr = 0.0D+00
  sm = 0.0D+00
  t2a = 0.5D+00 * ( func ( a ) + func ( b ) )
  t2 = t2a
  tb = abs ( t2a )
  c = t2 * ba

  dt(1) = c
  dt(2:7) = 0.0D+00
 
  odd = .true.
  bu = .false.
 
  do m = 1, 15
 
    bo = ( 7 <= m )
    hm = ba / real ( n, kind = 8 )
!
!  N+1 is odd.
!
    if ( odd ) then
 
      do i = 1, n, 2
        arg = a + real ( i, kind = 8 ) * hm
        w = func ( arg )
        t2 = t2 + w
        tb = abs ( w ) + tb
      end do
 
      ent = t2
      tab = tb * abs ( hm )
      d(1) = 16.0D+00 / 9.0D+00
      d(3) = 4.0D+00 * d(1)
      d(5) = 4.0D+00 * d(3)
!
!  N+1 is even.
!
    else
 
      do i = 1, n, 6
        w = real ( i, kind = 8 ) * hm
        t1 = t1 + func ( a + w ) + func ( b - w )
      end do
 
      ent = t1 + t2a
      t2a = t2
      d(1) = 2.25D+00
      d(3) = 9.0D+00
      d(5) = 36.0D+00
 
    end if
 
    ddt = dt(1)
    t = ent*hm
    dt(1) = t
    ent = t
 
    if ( m < 7 ) then
      mr = m
      w = real ( n * n, kind = 8 )
      d(m) = w
    else
      mr = 6
      d(6) = 64.0D+00
      w = 144.0D+00
    end if
 
    do i = 1, mr
 
      d1 = d(i) * ddt
      den = d1 - ent
      e = ent - ddt
      tnt = ent
      v = 0.0D+00
      ent = 0.0D+00
 
      if ( epsin <= abs ( den ) ) then
        e = e / den
        v = tnt * e
        ent = d1 * e
        t = v + t
        ddt = dt(i+1)
      end if
 
      dt(i+1) = v
 
    end do
 
    ta = c
    c = t
    result = c
    if ( .not. bo ) then
      t = t - v
    end if

    v = t - ta
    t = v + t
    epsout = abs ( v )
 
    if ( ta < t ) then
      d1 = ta
      ta = t
      t = d1
    end if
 
    bo = bo .or. ( ta < gr .and. sm < t )
 
    if ( bu .and. bo .and. epsout < epsin * w * tab ) then
      go to 140
    end if
 
    gr = ta + epsin
    sm = t - epsin
    odd = .not. odd
    n = np1
    np1 = n + 1
    bu = bo
    d(2) = 4.0D+00
    d(4) = 16.0D+00
 
  end do
 
  bo = .false.
 
  140 continue
 
  v = rnderr * tab
 
  epsout = max ( epsout, v )

  if ( bo ) then
    ind = 1
  end if

  return
end
subroutine monte_carlo ( func, a, b, n, result )

!*****************************************************************************80
!
!! MONTE_CARLO estimates the integral of a function by Monte Carlo.
!
!  Modified:
!
!    18 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, name of external function to be
!    integrated.  This name must be in an external statement in the calling
!    program.  FUNC must be a function of one real argument.  The value
!    of the argument to FUNC is the variable of integration
!    which ranges from A to B.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, integer ( kind = 4 ) N, the number of points to use.
!
!    Output, real ( kind = 8 ) RESULT, the computed value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) x(n)

  call random_number ( harvest = x(1:n) )

  x(1:n) = a + ( b - a ) * x(1:n)

  result = 0.0D+00
  do i = 1, n
    result = result + func ( x(i) )
  end do

  result = result / real ( n, kind = 8 )

  return
end
subroutine plint ( ntab, xtab, ftab, a, b, result )

!*****************************************************************************80
!
!! PLINT approximates the integral of unequally spaced data.
!
!  Discussion:
!
!    The method uses piecewise linear interpolation.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 2.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), the function values, 
!    FTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) ftab(ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ind
  real ( kind = 8 ) result
  real ( kind = 8 ) slope
  real ( kind = 8 ) syl
  real ( kind = 8 ) xtab(ntab)

  result = 0.0D+00
!
!  Check the parameters:
!
  if ( ntab < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLINT - Fatal error!'
    write ( *, '(a,i8)' ) '  NTAB < 2, NTAB = ', ntab
    stop
  end if
 
  do i = 2, ntab
    if ( xtab(i) <= xtab(i-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PLINT - Fatal error!'
      write ( *, '(a)' ) '  Nodes not in strict increasing order.'
      write ( *, '(a,i8)' ) '  XTAB(I) <= XTAB(I-1) for I = ', i
      write ( *, '(a,g14.6)' ) '  XTAB(I)   = ', xtab(i)
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
      stop
    end if
  end do

  if ( a == b ) then
    return
  end if
!
!  If B < A, temporarily switch A and B, and store sign.
!
  if ( b < a ) then
    syl = b
    b = a
    a = syl
    ind = -1
  else
    syl = a
    ind = 1
  end if
!
!  Find ILO and IHI so that A <= XTAB(ILO) <= XTAB(IHI) <= B
!  with the possible exception that A and B may be in the same
!  interval, or completely to the right or left of the XTAB's.
!
  ilo = ntab + 1

  do i = 1, ntab
    if ( a <= xtab(i) ) then
      ilo = i
      exit
    end if
  end do

  ihi = 0

  do i = ntab, 1, -1
    if ( xtab(i) <= b ) then
      ihi = i
      exit
    end if
  end do
!
!  Treat special cases where A, B lie both to left or both to right
!  of XTAB interval, or in between same pair of XTAB's.
!
  if ( ihi == 0 ) then

    slope = ( ftab(2) - ftab(1) ) / ( xtab(2) - xtab(1) )
    fa = ftab(1) + slope * ( a - xtab(1) )
    fb = ftab(1) + slope * ( b - xtab(1) )
    result = 0.5D+00 * ( b - a ) * ( fa + fb )

  else if ( ilo == ntab + 1 ) then

    slope = ( ftab(ntab) - ftab(ntab-1) ) / ( xtab(ntab) - xtab(ntab-1) )
    fa = ftab(ntab-1) + slope * ( a - xtab(ntab-1) )
    fb = ftab(ntab-1) + slope * ( b - xtab(ntab-1) )
    result = 0.5D+00 * ( b - a ) * ( fa + fb )

  else if ( ihi + 1 == ilo ) then

    slope = ( ftab(ilo) - ftab(ihi) ) / ( xtab(ilo) - xtab(ihi) )
    fa = ftab(ihi) + slope * ( a - xtab(ihi) )
    fb = ftab(ihi) + slope * ( b - xtab(ihi) )
    result = 0.5D+00 * ( b - a ) * ( fa + fb )

  else
!
!  Carry out approximate integration.  We know that ILO is no greater
!  than IHI-1, but equality is possible; A and B may be on either side
!  of a single XTAB(I).  That's OK, then the loop below won't be executed
!  at all.
!
    result = 0.0D+00
    do i = ilo, ihi-1
      result = result + 0.5D+00 * ( xtab(i+1) - xtab(i) ) &
        * ( ftab(i) + ftab(i+1) )
    end do
!
!  Add contribution from A-ILO and IHI-B.
!  Still have to watch out if ILO = 1 or IHI=NTAB...
!
    if ( ilo == 1 ) then
      slope = ( ftab(2) - ftab(1) ) / ( xtab(2) - xtab(1) )
      fa = ftab(1) + slope * ( a - xtab(1) )
      result = result + 0.5D+00 * ( xtab(ilo) - a ) * ( fa + ftab(ilo) )
    else
      slope = ( ftab(ilo) - ftab(ilo-1) ) / ( xtab(ilo) - xtab(ilo-1) )
      fa = ftab(ilo-1) + slope * ( a - xtab(ilo-1) )
      result = result + 0.5D+00 * ( xtab(ilo) - a ) * ( fa + ftab(ilo) )
    end if
 
    if ( ihi == ntab ) then
      slope = ( ftab(ntab) - ftab(ntab-1) ) / ( xtab(ntab) - xtab(ntab-1) )
      fb = ftab(ntab-1) + slope * ( b - xtab(ntab-1) )
      result = result + 0.5D+00 * ( b - xtab(ntab) ) * ( fb + ftab(ntab) )
    else
      slope = ( ftab(ihi+1) - ftab(ihi) ) / ( xtab(ihi+1) - xtab(ihi) )
      fb = ftab(ihi) + slope * ( b - xtab(ihi) )
      result = result + 0.5D+00 * ( b - xtab(ihi) ) * ( fb + ftab(ihi) )
    end if

  end if
!
!  Restore original values of A and B, reverse sign of integral
!  because of earlier switch.
! 
  if ( ind /= 1 ) then
    ind = 1
    syl = b
    b = a
    a = syl
    result = -result
  end if
 
  return
end
subroutine qnc79 ( func, a, b, err, result, ierr, k )

!*****************************************************************************80
!
!! QNC79 approximates the integral of F(X) using Newton-Cotes quadrature.
!
!  Discussion:
!
!    QNC79 is a general purpose program for evaluation of one
!    dimensional integrals  of user defined functions.  QNC79 will
!    pick its own points for evaluation of the integrand and these
!    will vary from problem to problem.
!
!    Thus QNC79 is not designed to integrate over data sets.
!
!    Moderately smooth integrands will be integrated efficiently
!    and reliably.  For problems with strong singularities,
!    oscillations etc., the user may wish to use more sophisticated
!    routines such as those in QUADPACK.
!
!    One measure of the reliability of QNC79 is the output parameter
!    K, giving the number of integrand evaluations that were needed.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      FUNCTION FUNC ( X )
!    which evaluates the function at the point X.
!
!    Input, real ( kind = 8 ) A, lower limit of integral.
!
!    Input, real ( kind = 8 ) B, upper limit of integral.
!
!    Input, real ( kind = 8 ) ERR, is a requested error tolerance.
!    Normally pick a value, 0 < ERR < 1.0D-03.
!
!    Output, real ( kind = 8 ) RESULT, computed value of the integral.
!    Hopefully RESULT is accurate to within ERR times the
!    integral of ABS ( FUNC ( X ) ).
!
!    Output, integer ( kind = 4 ) IERR, a status code
!     1 RESULT most likely meets requested error tolerance.
!    -1 A and B are too nearly equal to allow normal integration.
!     2 RESULT probably does not meet requested error tolerance.
!
!    Output, integer ( kind = 4 ) K, the number of function evaluations
!    actually used to do the integration.
!    A value of K .GT. 1000 indicates a difficult problem.
!    Other programs may be more efficient.
!    QNC79 will gracefully give up if K exceeds 2000.
!
  implicit none

  integer ( kind = 4 ), parameter :: kmx = 2000

  real ( kind = 8 ) a
  real ( kind = 8 ) aa(40)
  real ( kind = 8 ) ae
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) bank
  real ( kind = 8 ) blocal
  real ( kind = 8 ) c
  real ( kind = 8 ) ce
  real ( kind = 8 ) ee
  real ( kind = 8 ) ef
  real ( kind = 8 ) eps
  real ( kind = 8 ) err
  real ( kind = 8 ) f(13)
  real ( kind = 8 ) f1(40)
  real ( kind = 8 ) f2(40)
  real ( kind = 8 ) f3(40)
  real ( kind = 8 ) f4(40)
  real ( kind = 8 ) f5(40)
  real ( kind = 8 ) f6(40)
  real ( kind = 8 ) f7(40)
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) hh(40)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: kml = 7
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmn
  integer ( kind = 4 ) lmx
  integer ( kind = 4 ) lr(40)
  integer ( kind = 4 ) nbits
  integer ( kind = 4 ) nib
  integer ( kind = 4 ), save :: nlmn = 2
  integer ( kind = 4 ) nlmx
  real ( kind = 8 ) q13
  real ( kind = 8 ) q7
  real ( kind = 8 ) q7l
  real ( kind = 8 ) q7r(40)
  real ( kind = 8 ) result
  real ( kind = 8 ) test
  real ( kind = 8 ) tol
  real ( kind = 8 ) vl(40)
  real ( kind = 8 ) vr
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) w4

  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  if ( icall /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Fatal error!'
    write ( *, '(a)' ) '  QNC79 was called recursively!'
    stop
  end if
 
  icall = 1
  w1 = 41.0D+00 / 140.0D+00
  w2  = 216.0D+00 / 140.0D+00
  w3 = 27.0D+00 / 140.0D+00
  w4  = 272.0D+00 / 140.0D+00
!
!  DIGITS ( X ) = number of base 2 digits in representation of X.
!
  nbits = int ( log10 ( 2.0D+00 ) &
    * real ( digits ( result ), kind = 8 ) / 0.30102000 )

  nlmx = min ( 40, ( nbits * 4 ) / 5 )
  result = 0.0D+00
  ierr = 1
  ce = 0.0D+00
  lmx = nlmx
  lmn = nlmn

  if ( b == 0.0D+00 ) then
    go to 3
  end if

  if ( sign ( 1.0D+00, b ) * a <= 0.0D+00 ) then
    go to 3
  end if

  c = abs ( 1.0D+00 - a / b )

  if ( 0.1D+00 < c ) then
    go to 3
  end if
 
  if ( c <= 0.0D+00 ) then
    ierr  = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Fatal error!'
    write ( *, '(a)' ) '  A and B are too close.'
    stop
  end if
 
  nib = int ( 0.5D+00 - log ( c ) / log(2.0D+00) )
  lmx = min ( nlmx, nbits-nib-4 )

  if ( lmx < 2 ) then
    go to 32
  end if

  lmn = min ( lmn, lmx )
 
3 continue
 
  tol = max ( abs ( err ), 2.0D+00**(5-nbits) )
  if ( err == 0.0D+00 ) then
    tol = sqrt ( epsilon ( tol ) )
  end if

  eps = tol
  hh(1) = ( b - a ) / 12.0D+00
  aa(1) = a
  lr(1) = 1
 
  do i = 1, 11, 2
    f(i) = func ( a + real ( i - 1, kind = 8 ) * hh(1) )
  end do
 
  blocal = b
  f(13) = func ( blocal )
  k = 7
  l = 1
  area = 0.0D+00
  q7 = 0.0D+00
  ef = 256.0D+00 / 255.0D+00
  bank = 0.0D+00
!
!  Compute refined estimates, estimate the error, etc.
!
5 continue
 
  do i = 2, 12, 2
    f(i) = func ( aa(l) + real ( i - 1, kind = 8 ) * hh(l) )
  end do
 
  k = k + 6
!
!  Compute left and right half estimates.
!
  q7l = hh(l) * ( ( w1 * ( f(1) + f(7) )   &
                  + w2 * ( f(2) + f(6) ) ) &
                + ( w3 * ( f(3) + f(5) ) + w4 * f(4) ) )

  q7r(l) = hh(l) * ( ( w1 * ( f(7) + f(13) ) + w2 * ( f(8) + f(12) ) ) &
                + ( w3 * ( f(9) + f(11) ) + w4 * f(10) ) )
!
!  Update estimate of integral of absolute value.
!
  area = area + ( abs ( q7l ) + abs ( q7r(l) ) - abs ( q7 ) )
!
!  Do not bother to test convergence before minimum refinement level.
!
  if ( l < lmn ) then
    go to 11
  end if
!
!  Estimate the error in new value for whole interval, Q13.
!
  q13 = q7l + q7r(l)
  ee = abs ( q7 - q13 ) * ef
!
!  Compute nominal allowed error.
!
  ae = eps * area
!
!  Borrow from bank account, but not too much.
!
  test = min ( ae + 0.8D+00 * bank, 10.0D+00 * ae )
!
!  Don't ask for excessive accuracy.
!
  test = max ( test, tol * abs ( q13 ), 0.00003D+00 * tol * area )
!
!  Now, did this interval pass or not?
!
  if ( ee <= test ) then
    go to 8
  end if

  go to 10
!
!  Have hit max refinement level - penalize the cumulative error.
!
7 continue
 
  ce = ce + ( q7 - q13 )
  go to 9
!
!  On good intervals accumulate the theoretical estimate.
!
8 continue

  ce = ce + ( q7 - q13 ) / 255.0D+00
!
!  Update the bank account.  Don't go into debt.
!
9 continue

  bank = bank + ( ae - ee )
  if ( bank < 0.0D+00 ) then
    bank = 0.0D+00
  end if
!
!  Did we just finish a left half or a right half?
!
  if ( lr(l) <= 0 ) then
    go to 15
  end if

  go to 20
!
!  Consider the left half of the next deeper level.
!
10 continue
 
  if ( kmx < k ) then
    lmx = min ( kml, lmx )
  end if

  if ( lmx <= l ) then
    go to 7
  end if

11 continue

  l = l + 1
  eps = eps * 0.5D+00

  if ( l <= 17 ) then
    ef = ef / sqrt ( 2.0D+00 )
  end if

  hh(l) = hh(l-1) * 0.5D+00
  lr(l) = -1
  aa(l) = aa(l-1)
  q7 = q7l
  f1(l) = f(7)
  f2(l) = f(8)
  f3(l) = f(9)
  f4(l) = f(10)
  f5(l) = f(11)
  f6(l) = f(12)
  f7(l) = f(13)
  f(13) = f(7)
  f(11) = f(6)
  f(9)  = f(5)
  f(7)  = f(4)
  f(5)  = f(3)
  f(3)  = f(2)
  go to 5
!
!  Proceed to right half at this level
!
15 continue

  vl(l) = q13

16 continue

  q7 = q7r(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 12.0D+00 * hh(l)
  f(1)  = f1(l)
  f(3)  = f2(l)
  f(5)  = f3(l)
  f(7)  = f4(l)
  f(9)  = f5(l)
  f(11) = f6(l)
  f(13) = f7(l)
  go to 5
!
!  Left and right halves are done, so go back up a level
!
20 continue

  vr = q13

22 continue

  if ( l <= 1 ) then
    go to 30
  end if

  if ( l <= 17 ) then
    ef = ef * sqrt ( 2.0D+00 )
  end if

  eps = eps * 2.0D+00
  l = l-1
 
  if ( lr(l) <= 0 ) then
    vl(l) = vl(l+1) + vr
    go to 16
  else
    vr = vl(l+1) + vr
    go to 22
  end if
 
   30 continue
 
  if ( 2.0D+00 * tol * area < abs ( ce ) ) then
    ierr = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Warning!'
    write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
  end if
 
32    continue
 
  result = vr
 
  icall = 0
 
  return
end
subroutine quad ( func, a, b, abserr, relerr, nleast, nmost, work, &
  result )

!*****************************************************************************80
!
!! QUAD approximates the integral of F(X) by Romberg integration.
!
!  Discussion:
!
!    The integration is repeated until convergence of the results, 
!    or the maximum number of steps is taken.  The Romberg method 
!    uses the trapezoidal rule, subinterval halving, and Richardson 
!    extrapolation.
!
!    Convergence is declared if either of the following occurs:
!
!      ABS ( RESULT - INTEGRAL ) < ABSERR
!
!    or
!
!      RESULT = integral ( A <= X <= B ) (1+Y(X)) * FUNC ( X ) DX 
!
!    for some function Y with ABS ( Y(X) ) < RELERR+RNDERR  with RNDERR the
!    machine rounding error.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    CF Dunkl,
!    Romberg quadrature to prescribed accuracy,
!    SHARE file number 7090-1481 TYQUAD
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to 
!    be integrated.  The user must declare the name an external parameter in the
!    calling program, pass the name of the function in FUNC,
!    and write a function of the form FUNCTION FUNC ( X ) which
!    evaluates the function at the point X.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, real ( kind = 8 ) ABSERR, the absolute error tolerance.
!
!    Input, real ( kind = 8 ) RELERR, the relative error tolerance.
!
!    Input, integer ( kind = 4 ) NLEAST, the least number of times the integration
!    is to be carried out before the convergence test is made.
!    A value 3 <= NLEAST <= 15 is suggested.
!
!    Input, integer ( kind = 4 ) NMOST, the most number of times the
!    integration is to be carried out.
!
!    Workspace, real ( kind = 8 ) WORK(NMOST+1).
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) nmost

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) f
  real ( kind = 8 ) fcna
  real ( kind = 8 ) fcnb
  real ( kind = 8 ) fcnxi
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nleast
  integer ( kind = 4 ) nx
  real ( kind = 8 ) qx1
  real ( kind = 8 ) qx2
  real ( kind = 8 ) relerr
  real ( kind = 8 ) result
  real ( kind = 8 ) rnderr
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sumabs
  real ( kind = 8 ) t
  real ( kind = 8 ) tabs
  real ( kind = 8 ) work(nmost+1)
  real ( kind = 8 ) x

  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  rnderr = epsilon ( 1.0D+00 )
  qx1 = 0.0
  qx2 = 0.0
  h = b - a
  fcna = func ( a )
  fcnb = func ( b )
  tabs = abs ( h ) * ( abs ( fcna ) + abs ( fcnb ) ) / 2.0D+00
  t = h * ( fcna + fcnb ) / 2.0D+00
  nx = 1
 
  do i = 1, nmost
 
    h = 0.5 * h
 
    sum1 = 0.0D+00
    sumabs = 0.0D+00
    do j = 1, nx
      fcnxi = func ( a + h * real ( 2 * j - 1, kind = 8 ) )
      sumabs = sumabs + abs ( fcnxi )
      sum1 = sum1 + fcnxi
    end do
 
    t = 0.5D+00 * t + h * sum1
    tabs = tabs / 2.0D+00 + abs ( h ) * sumabs
    work(i) = 2.0D+00 * ( t + h * sum1 ) / 3.0D+00

    if ( 1 < i ) then
!
!  Construct difference table for Richardson extrapolation.
!
      f = 4.0D+00
      do j = 2, i
        k = i + 1 - j
        f = f * 4.0D+00
        work(k) = work(k+1) + ( work(k+1) - work(k) ) / ( f - 1.0D+00 )
      end do
!
!  Perform acceptance check.
!
      if ( nleast <= i ) then

        x = abs ( work(1) - qx2 ) + abs ( qx2 - qx1 )

        if ( x <= 3.0D+00 * tabs * ( abs ( relerr ) + rnderr ) .or. &
             x <= 3.0D+00 * abs ( abserr ) ) then
          result = work(1)
          return
        end if

      end if
!
!  Save old result, perform bisection, repeat.
!
      qx1 = qx2

    end if
 
    qx2 = work(1)
    nx = nx * 2
 
  end do

  result = work(1)

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns N values, evenly spaced between ALO and AHI.
!
!  Modified:
!
!    17 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine rminsp ( func, a, b, epsin, epsout, iop, result )

!*****************************************************************************80
!
!! RMINSP approximates the integral of a function using Romberg integration.
!
!  Discussion:
!
!    Both the midpoint and trapezoidal rule are used,
!    the intervals are repeatedly bisected, and Richardson
!    extrapolation is carried out to achieve a high accuracy.
!
!    RMINSP can carry out a cosine transform of the integral.  The
!    only effect this has is to handle a function F(X) which has
!    singularities near the endpoints.  This transform is based on
!    the fact that
!
!      Integral ( -1 <= X <= 1 ) F(X) DX
!
!    equals
!
!      Integral ( 0 <= Z <= PI ) F(COS(Z)) * SIN(Z)  DZ
!
!    If suitable accuracy is not achieved, the internal variable
!    NUPPER might be increased.  Its current value of 9 corresponds
!    to a maximum of 1024 subintervals and 1025 function evaluations.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Robert Kubik,
!    Algorithm 257:
!    Havie Integrator,
!    Communications of the ACM,
!    Volume 8, Number 6, June 1965, page 381.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!       FUNCTION FUNC ( X )
!    which evaluates the function at the point X.
!
!    Input, real ( kind = 8 ) A, lower limit of integration.
!
!    Input, real ( kind = 8 ) B, upper limit of integration.
!
!    Input, real ( kind = 8 ) EPSIN, requested relative error tolerance.
!
!    Output, real ( kind = 8 ) EPSOUT, estimated achieved relative error.
!
!    Input, integer ( kind = 4 ) IOP, method switch:
!    1, Use ordinary algorithm.
!    2, Use cosine transformation to decrease effect of
!       singularities near the endpoints.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: nupper = 9

  real ( kind = 8 ) a
  real ( kind = 8 ) acof(11)
  real ( kind = 8 ) alf
  real ( kind = 8 ) alfnj
  real ( kind = 8 ) alfno
  real ( kind = 8 ) ar
  real ( kind = 8 ) b
  real ( kind = 8 ) bcof(nupper+1)
  real ( kind = 8 ) bet
  real ( kind = 8 ) betnj
  real ( kind = 8 ) betno
  real ( kind = 8 ) const1
  real ( kind = 8 ) const2
  real ( kind = 8 ) deltan
  real ( kind = 8 ) endpts
  real ( kind = 8 ) epsin
  real ( kind = 8 ) epsout
  real ( kind = 8 ) error
  real ( kind = 8 ), parameter :: fac1 = 0.411233516712057D+00
  real ( kind = 8 ), parameter :: fac2 = 0.822467033441132D+00
  real ( kind = 8 ) factor
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) gamman
  real ( kind = 8 ) hnstep
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) iop
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhalf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rn
  real ( kind = 8 ) rnderr
  real ( kind = 8 ) result
  real ( kind = 8 ) rounde
  real ( kind = 8 ) tend
  real ( kind = 8 ) triarg
  real ( kind = 8 ) umid
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xplus

  result = 0.0D+00
 
  if ( a == b ) then
    return
  end if
!
!  Set coefficients in formula for accumulated roundoff error,
!  rounde = rnderr*(r1+r2*n), where r1, r2 are two empirical constants
!  and n is the current number of function values used.
!
  rnderr = epsilon ( 1.0D+00 )
 
  if ( iop == 2 ) then
    r1 = 50.0D+00
  else
    r1 = 1.0D+00
  end if

  if ( iop == 1 ) then
    r2 = 0.02D+00
  else
    r2 = 2.0D+00
  end if

  error = epsin
!
!  Initial calculations.
!
  alf = 0.5D+00 * ( b - a )
  bet = 0.5D+00 * ( b + a )
  acof(1) = func ( a ) + func ( b )
  bcof(1) = func ( bet )
!
!  Modified Romberg algorithm, ordinary case.
!
  if ( iop /= 2 ) then

    hnstep = 2.0D+00
    bcof(1) = hnstep * bcof(1)
    factor = 1.0D+00
!
!  Modified Romberg, cosine transformed case.
!
  else

    hnstep = pi
    ar = fac1
    endpts = acof(1)
    acof(1) = fac2 * acof(1)
    bcof(1) = hnstep * bcof(1) - ar * endpts
    factor = 4.0D+00
    ar = ar / 4.0D+00
    triarg = pi / 4.0D+00
    alfno = -1.0D+00
  end if
 
  hnstep = 0.5D+00 * hnstep
  nhalf = 1
  n = 2
  rn = 2.0D+00
  acof(1) = 0.5D+00 * ( acof(1) + bcof(1) )
  acof(2) = acof(1) - ( acof(1) - bcof(1) ) / ( 4.0D+00 * factor - 1.0D+00 )
!
!  End of initial calculation.
!
!  Start actual calculations.
!
  do i = 1, nupper
 
    umid = 0.0D+00
!
!  Modified Romberg algorithm, ordinary case.
!  compute first element in mid-point formula for ordinary case
!
    if ( iop == 1 ) then
 
      alfnj = 0.5D+00 * hnstep
 
      do j = 1, nhalf
        xplus = alf * alfnj + bet
        xmin = -alf * alfnj + bet
        umid = umid + func ( xplus ) + func ( xmin )
        alfnj = alfnj + hnstep
      end do
 
      umid = hnstep * umid
!
!  Modified Romberg algorithm, cosine transformed case
!  compute first element in mid-point formula for cosine transformed
!  Romberg algorithm.
!
    else if ( iop == 2 ) then
 
      const1 = -sin(triarg)
      const2 = 0.5D+00 * alfno / const1
 
      alfno = const1
      betno = const2
      gamman = 1.0D+00 - 2.0D+00 * alfno**2
      deltan = -2.0D+00 * alfno * betno
 
      do j = 1, nhalf
        alfnj = gamman * const1 + deltan * const2
        betnj = gamman * const2 - deltan * const1
        xplus = alf * alfnj + bet
        xmin = -alf * alfnj + bet
        umid = umid + betnj * ( func ( xplus ) + func ( xmin ) )
        const1 = alfnj
        const2 = betnj
      end do
 
      umid = hnstep * umid - ar * endpts
      ar = ar / 4.0D+00
 
    end if
!
!  Modified Romberg algorithm, calculate (i+1)-th row in the U table.
!
    const1 = 4.0D+00 * factor
    index = i+1
 
    do j = 2, i+1
      tend = umid + ( umid - bcof(j-1) ) / ( const1 - 1.0D+00 )
      bcof(j-1) = umid
      umid = tend
      const1 = 4.0D+00 * const1
    end do
 
    bcof(i+1) = tend
    xplus = const1
!
!  Calculation of (i+1)-th row in the U table is finished
!
!  Test to see if the required accuracy has been obtained.
!
    epsout = 1.0D+00
    iout = 1
 
    do j = 1, index
 
      const1 = 0.5D+00 * ( acof(j) + bcof(j) )
      const2 = 0.5D+00 * abs ( ( acof(j) - bcof(j) ) / const1 )
 
      if ( const2 <= epsout ) then
        epsout = const2
        iout = j
      end if
 
      acof(j) = const1
 
    end do
!
!  Testing on accuracy finished
!
    if ( iout == index ) then
      iout = iout + 1
    end if

    acof(index+1) = acof(index) - ( acof(index) - bcof(index) ) &
      / ( xplus - 1.0D+00 )

    rounde = rnderr * ( r1 + r2 * rn)

    epsout = max ( epsout, rounde )
    error = max ( error, rounde )
 
    if ( epsout <= error ) then
      go to 10
    end if
 
    nhalf = n
    n = 2 * n
    rn = 2.0D+00 * rn
    hnstep = 0.5D+00 * hnstep

    if ( 1 < iop ) then
      triarg = 0.5D+00 * triarg
    end if
 
  end do
!
!  Accuracy not reached with maximum number of subdivisions.
!
  n = nhalf
!
!  Calculation for modified Romberg algorithm finished.
!
  10  continue
 
  n = 2 * n
  index = index - 1
  n = n + 1
  j = iout

  if ( index <= j - 1 ) then
    j = index
  end if

  tend = alf * ( 2.0D+00 * acof(j) - bcof(j) )
  umid = alf * bcof(j)
  result = alf * acof(iout)

  return
end
subroutine simp ( func, a, b, eps, result )

!*****************************************************************************80
!
!! SIMP approximates the integral of a function by an adaptive Simpson's rule.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    James Lyness,
!    Algorithm 379:
!    SQUANK - Simpson quadrature used adaptively, noise killed,
!    Communications of the ACM,
!    Volume 13, Number 4, April 1970, pages 260-263.
!
!    William McKeeman, Lawrence Tesler,
!    Algorithm 182:
!    Nonrecursive adaptive integration,
!    Communications of the ACM,
!    Volume 6, 1963, page 315.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      FUNCTION FUNC ( X )
!    which evaluates the function at the point X.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, real ( kind = 8 ) EPS, the requested error tolerance.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxlev = 30

  real ( kind = 8 ) a
  real ( kind = 8 ) a1
  real ( kind = 8 ) absar
  real ( kind = 8 ) b
  real ( kind = 8 ) da
  real ( kind = 8 ) dx(maxlev)
  real ( kind = 8 ) ep
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsp(maxlev)
  real ( kind = 8 ) est
  real ( kind = 8 ) est1
  real ( kind = 8 ) est2(maxlev)
  real ( kind = 8 ) est3(maxlev)
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2(maxlev)
  real ( kind = 8 ) f3(maxlev)
  real ( kind = 8 ) f4(maxlev)
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fbp(maxlev)
  real ( kind = 8 ) fm
  real ( kind = 8 ) fmp(maxlev)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lvl
  integer ( kind = 4 ) nrtr(maxlev)
  real ( kind = 8 ) pval(maxlev,3)
  real ( kind = 8 ) result
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sx
  real ( kind = 8 ) x2(maxlev)
  real ( kind = 8 ) x3(maxlev)

  result = 0.0D+00
 
  if ( a == b ) then
    return
  end if
 
  ep = eps
  a1 = a
  nrtr(1:maxlev) = 0
  pval(1:maxlev,1:3) = 0.0D+00
 
  lvl = 0
  absar = 0.0D+00
  est = 0.0D+00
  da = b - a1

  fa = func ( a1 )
  fm = 4.0D+00 * func ( ( a1 + b ) * 0.5D+00 )
  fb = func ( b )
!
!  1 = RECUR
!
   30 continue
 
  lvl = lvl + 1
  dx(lvl) = da / 3.0D+00
  sx = dx(lvl) / 6.0D+00
  f1 = 4.0D+00 * func ( 0.5D+00 * dx(lvl) + a1 )
  x2(lvl) = a1 + dx(lvl)
  f2(lvl) = func ( x2(lvl) )
  x3(lvl) = x2(lvl) + dx(lvl)
  f3(lvl) = func ( x3(lvl) )
  epsp(lvl) = ep
  f4(lvl) = 4.0D+00 * func ( dx(lvl) * 0.5D+00 + x3(lvl) )
  fmp(lvl) = fm
  est1 = sx * (fa+f1+f2(lvl))
  fbp(lvl) = fb
  est2(lvl) = sx * ( f2(lvl) + f3(lvl) + fm )
  est3(lvl) = sx * ( f3(lvl) + f4(lvl) + fb )
  sum1 = est1 + est2(lvl) + est3(lvl)
  absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) &
    + abs ( est3(lvl) )

  if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) then
    go to 40
  end if

  if ( maxlev <= lvl ) then
    go to 50
  end if
!
!  2 = UP
!
40 continue
 
  if ( 1 < lvl ) then
    lvl = lvl-1
  end if

  l = nrtr(lvl)

  if ( l == 0 ) then
    go to 50
  end if

  pval(lvl,l) = sum1

  if ( l == 1 ) then
    go to 60
  end if

  if ( l == 2 ) then
    go to 70
  end if

  if ( l == 3 ) then
    go to 80
  end if
 
50 continue
 
  nrtr(lvl) = 1
  est = est1
  fm = f1
  fb = f2(lvl)
  ep = epsp(lvl) / 1.7D+00
  da = dx(lvl)
  go to 30
 
60 continue
 
  nrtr(lvl) = 2
  fa = f2(lvl)
  fm = fmp(lvl)
  fb = f3(lvl)
  est = est2(lvl)
  a1 = x2(lvl)
  ep = epsp(lvl) / 1.7D+00
  da = dx(lvl)
  go to 30
 
70 continue
 
  nrtr(lvl) = 3
  fa = f3(lvl)
  fm = f4(lvl)
  fb = fbp(lvl)
  est = est3(lvl)
  a1 = x3(lvl)
  ep = epsp(lvl) / 1.7D+00
  da = dx(lvl)
  go to 30
 
80 continue
 
  sum1 = pval(lvl,1) + pval(lvl,2) + pval(lvl,3)

  if ( 1 < lvl ) then
    go to 40
  end if
 
90 continue
 
  result = sum1
 
  return
end
subroutine simpne ( ntab, x, y, result )

!*****************************************************************************80
!
!! SIMPNE approximates the integral of unevenly spaced data.
!
!  Discussion:
!
!    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
!    to the data and integrates that exactly.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, number of data points.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) X(NTAB), contains the X values of the data,
!    in order.
!
!    Input, real ( kind = 8 ) Y(NTAB), contains the Y values of the data.
!
!    Output, real ( kind = 8 ) RESULT.
!    RESULT is the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) del(3)
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) feints
  real ( kind = 8 ) g(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) pi(3)
  real ( kind = 8 ) result
  real ( kind = 8 ) sum1
  real ( kind = 8 ) x(ntab)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y(ntab)

  result = 0.0D+00
 
  if ( ntab <= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIMPNE - Fatal error!'
    write ( *, '(a)' ) '  NTAB <= 2.'
    stop
  end if
 
  n = 1
 
  do
 
    x1 = x(n)
    x2 = x(n+1)
    x3 = x(n+2)
    e = x3 * x3- x1 * x1
    f = x3 * x3 * x3 - x1 * x1 * x1
    feints = x3 - x1

    del(1) = x3 - x2
    del(2) = x1 - x3
    del(3) = x2 - x1

    g(1) = x2 + x3
    g(2) = x1 + x3
    g(3) = x1 + x2

    pi(1) = x2 * x3
    pi(2) = x1 * x3
    pi(3) = x1 * x2
 
    sum1 = 0.0D+00
    do i = 1, 3
      sum1 = sum1 + y(n-1+i) * del(i) &
        * ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
    end do
    result = result - sum1 / ( del(1) * del(2) * del(3) )
 
    n = n + 2

    if ( ntab <= n + 1 ) then
      exit
    end if

  end do
 
  if ( mod ( ntab, 2 ) /= 0 ) then
    return
  end if

  n = ntab - 2
  x3 = x(ntab)
  x2 = x(ntab-1)
  x1 = x(ntab-2)
  e = x3 * x3 - x2 * x2
  f = x3 * x3 * x3 - x2 * x2 * x2
  feints = x3 - x2

  del(1) = x3 - x2
  del(2) = x1 - x3
  del(3) = x2 - x1

  g(1) = x2 + x3
  g(2) = x1 + x3
  g(3) = x1 + x2

  pi(1) = x2 * x3
  pi(2) = x1 * x3
  pi(3) = x1 * x2
 
  sum1 = 0.0D+00
  do i = 1, 3
    sum1 = sum1 + y(n-1+i) * del(i) * &
      ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
  end do
 
  result = result - sum1 / ( del(1) * del(2) * del(3) )
 
  return
end
subroutine simpsn ( ntab, h, y, result )

!*****************************************************************************80
!
!! SIMPSN approximates the integral of evenly spaced data.
!
!  Discussion:
!
!    Simpson's rule is used.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of data points.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) H, specifies the increment between the
!    X values.  Note that the actual X values are not needed,
!    just the constant spacing!
!
!    Input, real ( kind = 8 ) Y(NTAB), the data.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral
!    from the first to the last point.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) del(3)
  real ( kind = 8 ) f
  real ( kind = 8 ) g(3)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) pii(3)
  real ( kind = 8 ) result
  real ( kind = 8 ) sum1
  real ( kind = 8 ) y(ntab)

  result = 0.0D+00
 
  if ( ntab <= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIMPSN - Fatal error!'
    write ( *, '(a,i8)' ) '  NTAB < 2, NTAb = ', ntab
    stop
  end if
 
  if ( mod ( ntab, 2 ) == 0 ) then
    n = ntab - 1
  else
    n = ntab
  end if
 
  result = y(1) + y(n) + 4.0D+00 * y(n-1)
  do i = 2, n-2, 2
    result = result + 4.0D+00 * y(i) + 2.0D+00 * y(i+1)
  end do
  result = h * result / 3.0D+00
 
  if ( mod ( ntab, 2 ) == 1 ) then
    return
  end if
 
  f = h**3
  del(1) = h
  del(2) = -2.0D+00 * h
  del(3) = h
  g(1) = h
  g(2) = 0.0D+00
  g(3) = -h
  pii(1) = 0.0D+00
  pii(2) = -h**2
  pii(3) = 0.0D+00
  n = n-1
 
  sum1 = 0.0D+00
  do i = 1, 3
    sum1 = sum1 + y(n-1+i) * del(i) * &
      ( f / 3.0D+00 - g(i) * 0.5D+00 * h**2 + pii(i) * h )
  end do
 
  result = result + 0.5D+00 * sum1 / h**3
 
  return
end
function solve ( shift, n, a, b )

!*****************************************************************************80
!
!! SOLVE solves a special linear system.
!
!  Discussion:
!
!    SOLVE solves for the N-th component of the solution DELTA to the equation
!
!      (Jn - shift*Identity) * DELTA  = En,
!
!    En is the vector of all zeroes except for 1 in the N-th position.
!
!    The matrix Jn is symmetric tridiagonal, with diagonal
!    elements A(I), off-diagonal elements B(I).  This equation
!    must be solved to obtain the appropriate changes in the lower
!    2 by 2 submatrix of coefficients for orthogonal polynomials.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SHIFT, the value of the factor that multiplies
!    the identity matrix, in the definition of the system matrix.
!
!    Input, integer ( kind = 4 ) N, the index of the desired component.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) shift
  real ( kind = 8 ) solve

  alpha = a(1) - shift
  do i = 2, n-1
    alpha = a(i) - shift - b(i-1)**2 / alpha
  end do
 
  solve = 1.0D+00 / alpha
 
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
subroutine wedint ( ntab, h, ftab, result )

!*****************************************************************************80
!
!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, is the number of data points.  
!    (NTAB-1) must be divisible by 6.
!
!    Input, real ( kind = 8 ) H, is the spacing between the points at which
!    the data was evaluated.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated data values.
!
!    Output, real ( kind = 8 ) RESULT, is the approximation to the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  real ( kind = 8 ) result

  result = 0.0D+00
 
  if ( ntab <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDINT - Fatal error!'
    write ( *, '(a)' ) '  NTAB < 2'
    write ( *, '(a,i8)' ) '  NTAB = ', ntab
    stop
  end if
 
  if ( mod ( ntab, 6 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDINT - Fatal error!'
    write ( *, '(a)' ) '  NTAB must equal 6*N+1 for some N!'
    stop
  end if
 
  do i = 1, ntab-6, 6
    result = result & 
      +           ftab(i)   &
      + 5.0D+00 * ftab(i+1) &
      +           ftab(i+2) &
      + 6.0D+00 * ftab(i+3) &
      +           ftab(i+4) &
      + 5.0D+00 * ftab(i+5) &
      +           ftab(i+6)
  end do
 
  result = 3.0D+00 * h * result / 10.0D+00
 
  return
end
