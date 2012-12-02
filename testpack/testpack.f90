program main

!*****************************************************************************80
!
!! MAIN is the main program for testing ADAPT and TESTPACK.
!
!  Discussion:
!
!    ADAPT is a multidimensional quadrature program by Alan Genz.
!
!    TESTPACK is a collection of six test integrand functions.
!
!    The routine MULTST tests ADAPT using the test integrands.
!
!  Modified:
!
!    21 March 2007
!
!  Author:
!
!    Original FORTRAN77 version by Alan Genz.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndiml = 5
  integer ( kind = 4 ), parameter :: tstlim = 6
  integer ( kind = 4 ), parameter :: tstmax = 6

  external adapt
  real ( kind = 8 ), dimension ( tstmax ) :: difclt = (/ &
    110.0D+00, 600.0D+00, 600.0D+00, 100.0D+00, 150.0D+00, &
    100.0D+00 /)
  real ( kind = 8 ), dimension ( tstmax ) :: expnts = (/ &
    1.5D+00, 2.0D+00, 2.0D+00, 1.0D+00, 2.0D+00, 2.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: maxpts = 10000
  integer ( kind = 4 ), dimension ( ndiml ) :: ndims = (/ &
    2, 3, 4, 6, 8 /)
  integer ( kind = 4 ) :: nsamp = 20
  real ( kind = 8 ) :: rel_tol = 1.0E-06
  character ( len = 6 ) :: sbname = 'ADAPT '
  integer ( kind = 4 ), dimension ( tstlim ) :: tstfns = (/ &
    1, 2, 3, 4, 5, 6 /)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TESTPACK'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call MULTST, which can test a routine that'
  write ( *, '(a)' ) '  is designed to estimate multidimensional'
  write ( *, '(a)' ) '  integrals, by numerical quadrature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The routine to be tested here is called ADAPT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The test integrands are Genz''s standard set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MULTST, ADAPT and the test integrands were'
  write ( *, '(a)' ) '  written in FORTRAN77 by Alan Genz.'

  call multst ( nsamp, tstlim, tstfns, tstmax, difclt, &
    expnts, ndiml, ndims, sbname, adapt, rel_tol, maxpts )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TESTPACK'
  write ( *, '(a)' ) '  Normal end of execution'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine adapt ( ndim, a, b, minpts, maxpts, functn, rel_tol, &
  itest, alpha, beta, lenwrk, wrkstr, relerr, finest, ifail )

!*****************************************************************************80
!
!! ADAPT carries out adaptive multidimensional quadrature.
!
!  Modified:
!
!    21 March 2007
!
!  Author:
!
!    Original FORTRAN77 version by Alan Genz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDIM, the number of variables.
!    2 <= NDIM.
!
!    Input, real ( kind = 8 ) A(NDIM), the lower limits of integration.
!
!    Input, real ( kind = 8 ) B(NDIM), the upper limits of integration.
!
!    Input/output, integer ( kind = 4 ) MINPTS, the minimum number of function
!    evaluations to be allowed,  MINPTS must not exceed MAXPTS.  If MINPTS < 0
!    then the routine assumes a previous call has been made with the same 
!    integrand and continues that calculation.
!
!    Input, integer ( kind = 4 ) MAXPTS, the maximum number of function
!    evaluations allowed, which must be at least RULCLS, where
!    RULCLS = 2**NDIM + 2 * NDIM**2 + 2 * NDIM + 1, when NDIM <= 15 and
!    RULCLS = ( NDIM * ( 14 - NDIM * ( 6 - 4 * NDIM ) ) ) / 3 + 1,
!    when 15 < NDIM.
!    for NDIM  =  2   3   4   5   6   7   8   9   10   11   12
!    RULCLS   =  17  33  57  93 149 241 401 693 1245 2313 4409
!    A suggested starting value for MAXPTS is 100*RULCLS.  If
!    this is not large enough for the required accuracy, then
!    MAXPTS and LENWRK should be increased accordingly.
!
!    Input, external FUNCTN, the user-defined function
!    to be integrated.  It must have the form
!      subroutine functn ( indx, ndim, z, alpha, beta, f )
!    where
!      INDX is the index of the test function,
!      NDIM is the spatial dimension,
!      Z(1:NDIM) is the evaluation point,
!      ALPHA(1:NDIM) is a set of parameters,
!      BETA(1:NDIM) is a set of parameters.
!      F is the function value.
!
!    Input, real ( kind = 8 ) REL_TOL, the user's requested relative accuracy.
!
!    Input, integer ( kind = 4 ) ITEST, the index of the test.
!
!    Input, real ( kind = 8 ) ALPHA(NDIM), BETA(NDIM), parameters
!    associated with the integrand function.
!
!    Input, integer ( kind = 4 ) LENWRK, the length of the array WRKSTR.
!    The routine needs (2*NDIM+3)*(1+MAXPTS/RULCLS)/2 for LENWRK if
!    MAXPTS function calls are used.
!
!    Input/output, real ( kind = 8 ) WRKSTR(LENWRK).  This array does not
!    need to be set or inspected by the user.  However, the output value of
!    WKRSTR from one call may be needed by the program on a followup call
!    if the input value of MINPTS < 0, which signals that another calculation
!    is requested for the same integrand.
!
!    Output, real ( kind = 8 ) RELERR, the estimated relative accuracy
!    of the integral estimate.
!
!    Output, real ( kind = 8 ) FINEST, the estimated value of integral.
!
!    Output, integer ( kind = 4 ) IFAIL
!    * 0, for normal exit, when estimated relative error RELERR is less
!    than REL_TOL, and with MAXPTS or less function calls made.
!    * 1, if MAXPTS was too small for ADAPT to obtain the required relative
!    error REL_TOL.  In this case ADAPT returns a value of FINEST with
!    estimated relative error RELERR.
!    * 2, if LENWRK was too small for MAXPTS function calls.  In
!    this case ADAPT returns a value of FINEST with estimated error
!    RELERR using the working storage available, but RELERR is likely to
!    be greater than REL_TOL.
!    * 3, if NDIM < 2 or MAXPTS < MINPTS or MAXPTS < RULCLS.
!
  implicit none

  integer ( kind = 4 ) lenwrk
  integer ( kind = 4 ) ndim

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) alpha(ndim)
  real ( kind = 8 ) b(ndim)
  real ( kind = 8 ) beta(ndim)
  real ( kind = 8 ) center(ndim)
  real ( kind = 8 ) df1
  real ( kind = 8 ) df2
  real ( kind = 8 ) dif
  real ( kind = 8 ) difmax
  integer ( kind = 4 ) divaxn
  integer ( kind = 4 ) divaxo
  integer ( kind = 4 ) divflg
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ) finest
  integer ( kind = 4 ) funcls
  external functn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifail
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) index2
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) lambda2
  real ( kind = 8 ) lambda4
  real ( kind = 8 ) lambda5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) maxpts
  integer ( kind = 4 ) minpts
  integer ( kind = 4 ) n
  real ( kind = 8 ) ratio
  real ( kind = 8 ) rel_tol
  real ( kind = 8 ) relerr
  real ( kind = 8 ) rgncmp
  real ( kind = 8 ) rgnerr
  integer ( kind = 4 ), save :: rgnstr = 0
  real ( kind = 8 ) rgnval
  real ( kind = 8 ) rgnvol
  integer ( kind = 4 ) rulcls
  integer ( kind = 4 ), save :: sbrgns = 0
  integer ( kind = 4 ) sbtmpp
  integer ( kind = 4 ) subrgn
  integer ( kind = 4 ) subtmp
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) sum3
  real ( kind = 8 ) sum4
  real ( kind = 8 ) sum5
  real ( kind = 8 ) weit1
  real ( kind = 8 ) weit2
  real ( kind = 8 ) weit3
  real ( kind = 8 ) weit4
  real ( kind = 8 ) weit5
  real ( kind = 8 ) weitp1
  real ( kind = 8 ) weitp2
  real ( kind = 8 ) weitp3
  real ( kind = 8 ) weitp4
  real ( kind = 8 ) width(ndim)
  real ( kind = 8 ) widthl(ndim)
  real ( kind = 8 ) wrkstr(lenwrk)
  real ( kind = 8 ) z(ndim)

  ifail = 3
  relerr = 1.0D+00
  funcls = 0

  if ( ndim < 2 ) then

    minpts = 0
    wrkstr(lenwrk-1) = sbrgns

    relerr = 1.0D+00
    finest = 0.0D+00
    ifail = 3

    return
  end if

  if ( maxpts < minpts ) then

    minpts = 0
    wrkstr(lenwrk-1) = sbrgns

    relerr = 1.0D+00
    finest = 0.0D+00
    ifail = 3

    return
  end if

  if ( ndim <= 15 ) then
    rulcls = 2**ndim + 2 * ndim * ndim + 2 * ndim + 1
  else if ( 15 < ndim ) then
    rulcls = 1 + ( ndim * ( 12 + ( ndim - 1 ) &
      * ( 6 + ( ndim - 2 ) * 4 ) ) ) / 3
  end if

  if ( maxpts < rulcls ) then

    relerr = 1.0D+00
    finest = 0.0D+00
    ifail = 3

    return
  end if
!
!  Initialization.
!
  rgnstr = 2 * ndim + 3
  divaxo = 0
!
!  Basic rule initialization.
!
  lambda5 = 9.0D+00 / 19.0D+00

  if ( ndim <= 15 ) then

    lambda4 = 9.0D+00 / 10.0D+00
    lambda2 = 9.0D+00 / 70.0D+00
    weit5 = 1.0D+00 / ( 3.0D+00 * lambda5 )**3 / 2.0D+00**ndim

  else

    ratio = real ( ndim - 2, kind = 8 ) / 9.0D+00

    lambda4 = ( 1.0D+00 / 5.0D+00 - ratio ) &
      / ( 1.0D+00 / 3.0D+00 - ratio / lambda5 )

    ratio = ( 1.0D+00 - lambda4 / lambda5 ) &
      * real ( ndim - 1, kind = 8 ) * ratio / 6.0D+00

    lambda2 = ( 1.0D+00 / 7.0D+00 - lambda4 / 5.0D+00 - ratio ) &
      / ( 1.0D+00 / 5.0D+00 - lambda4 / 3.0D+00 - ratio / lambda5 )

    weit5 = 1.0D+00 / ( 6.0D+00 * lambda5 )**3

  end if

  weit4 = ( 1.0D+00 / 15.0D+00 - lambda5 / 9.0D+00 ) &
    / ( 4.0D+00 * ( lambda4 - lambda5 ) * lambda4**2 )

  weit3 = ( 1.0D+00 / 7.0D+00 - ( lambda5 + lambda2 ) / 5.0D+00 &
    + lambda5 * lambda2 / 3.0D+00 ) / ( 2.0D+00 * lambda4 &
    * ( lambda4 - lambda5 ) * ( lambda4 - lambda2 ) ) &
    - 2.0D+00 * real ( ndim - 1, kind = 8 ) * weit4

  weit2 = ( 1.0D+00 / 7.0D+00 - ( lambda5 + lambda4 ) / 5.0D+00 &
    + lambda5 * lambda4 / 3.0D+00 ) / ( 2.0D+00 * lambda2 &
    * ( lambda2 - lambda5 ) * ( lambda2 - lambda4 ) )

  if ( ndim <= 15 ) then
    weit1 = 1.0D+00 - 2.0D+00 * real ( ndim, kind = 8 ) &
      * ( weit2 + weit3 + real ( ndim - 1, kind = 8 ) * weit4 ) &
      - 2.0D+00**ndim * weit5
  else
    weit1 = 1.0D+00 - 2.0D+00 * real ( ndim, kind = 8 ) &
      * ( weit2 + weit3 + real ( ndim - 1, kind = 8 ) * &
      ( weit4 + 2.0D+00 * real ( ndim - 2, kind = 8 ) * weit5 / 3.0D+00 ) )
  end if

  weitp4 = 1.0D+00 / ( 6.0D+00 * lambda4 )**2

  weitp3 = ( 1.0D+00 / 5.0D+00 - lambda2 / 3.0D+00 ) / &
    ( 2.0D+00 * lambda4 * ( lambda4 - lambda2 ) ) &
    - 2.0D+00 * real ( ndim - 1, kind = 8 ) * weitp4

  weitp2 = ( 1.0D+00 / 5.0D+00 - lambda4 / 3.0D+00 ) &
    / ( 2.0D+00 * lambda2 * ( lambda2 - lambda4 ) )

  weitp1 = 1.0D+00 - 2.0D+00 * real ( ndim, kind = 8 ) * &
    ( weitp2 + weitp3 + real ( ndim - 1, kind = 8 ) * weitp4 )

  ratio = lambda2 / lambda4

  lambda5 = sqrt ( lambda5 )
  lambda4 = sqrt ( lambda4 )
  lambda2 = sqrt ( lambda2 )
!
!  End basic rule initialization.
!
  if ( minpts < 0 ) then

    sbrgns = int ( wrkstr(lenwrk-1) )
    divflg = 0
    subrgn = rgnstr
    wrkstr(lenwrk) = wrkstr(lenwrk) - wrkstr(subrgn)
    finest = finest - wrkstr(subrgn-1)
    divaxo = int ( wrkstr(subrgn-2) )

    do j = 1, ndim
      subtmp = subrgn - 2 * ( j + 1 )
      center(j) = wrkstr(subtmp+1)
      width(j) = wrkstr(subtmp)
    end do

    width(divaxo) = width(divaxo) / 2.0D+00
    center(divaxo) = center(divaxo) - width(divaxo)

  else

    width(1:ndim) = ( b(1:ndim) - a(1:ndim) ) / 2.0D+00
    center(1:ndim) = a(1:ndim) + width(1:ndim)

    finest = 0.0D+00
    wrkstr(lenwrk) = 0.0D+00
    divflg = 1
    subrgn = rgnstr
    sbrgns = rgnstr

  end if
!
!  Begin basic rule.
!
  do

    rgnvol = 2.0D+00**ndim * product ( width(1:ndim) )

    z(1:ndim) = center(1:ndim)

    call functn ( itest, ndim, z, alpha, beta, sum1 )
!
!  Compute symmetric sums of functn(lambda2,0,0,...,0) and
!  functn(lambda4,0,0,...,0), and maximum fourth difference.
!
    difmax = -1.0D+00
    sum2 = 0.0D+00
    sum3 = 0.0D+00

    do j = 1, ndim

      z(j) = center(j) - lambda2 * width(j)
      call functn ( itest, ndim, z, alpha, beta, f1 )
      z(j) = center(j) + lambda2 * width(j)
      call functn ( itest, ndim, z, alpha, beta, f2 )
      widthl(j) = lambda4 * width(j)
      z(j) = center(j) - widthl(j)
      call functn ( itest, ndim, z, alpha, beta, f3 )
      z(j) = center(j) + widthl(j)
      call functn ( itest, ndim, z, alpha, beta, f4 )
      sum2 = sum2 + f1 + f2
      sum3 = sum3 + f3 + f4
      df1 = f1 + f2 - 2.0D+00 * sum1
      df2 = f3 + f4 - 2.0D+00 * sum1
      dif = abs ( df1 - ratio * df2 )

      if ( difmax < dif ) then
        difmax = dif
        divaxn = j
      end if

      z(j) = center(j)

    end do

    if ( sum1 == sum1 + difmax / 8.0D+00 ) then
      divaxn = mod ( divaxo, ndim ) + 1
    end if
!
!  Compute symmetric sum of functn(lambda4,lambda4,0,0,...,0).
!
    sum4 = 0.0D+00

    do j = 2, ndim

      do k = j, ndim

        do l = 1, 2

          widthl(j-1) = -widthl(j-1)
          z(j-1) = center(j-1) + widthl(j-1)

          do m = 1, 2

            widthl(k) = -widthl(k)
            z(k) = center(k) + widthl(k)

            call functn ( itest, ndim, z, alpha, beta, f1 )

            sum4 = sum4 + f1

          end do

        end do

        z(k) = center(k)

      end do

      z(j-1) = center(j-1)

    end do
!
!  If NDIM < 16 compute symmetric sum of functn(lambda5,lambda5,...,lambda5).
!
    if ( ndim <= 15 ) then

      sum5 = 0.0D+00

      widthl(1:ndim) = -lambda5 * width(1:ndim)
      z(1:ndim) = center(1:ndim) + widthl(1:ndim)

      do

        call functn ( itest, ndim, z, alpha, beta, f1 )
        sum5 = sum5 + f1

        j = ndim

        do

          widthl(j) = - widthl(j)
          z(j) = center(j) + widthl(j)

          if ( 0.0D+00 <= widthl(j) ) then
            exit
          end if

          j = j - 1

          if ( j < 1 ) then
            exit
          end if

        end do

        if ( j < 1 ) then
          exit
        end if

      end do
!
!  If 15 < NDIM, compute symmetric sum of
!  FUNCTN(lambda5,lambda5,lambda5,0,0,...,0).
!
    else

      sum5 = 0.0D+00

      widthl(1:ndim) = lambda5 * width(1:ndim)

      do i = 3, ndim
        do j = i, ndim
          do k = j, ndim

            do l = 1, 2
              widthl(i-2) = -widthl(i-2)
              z(i-2) = center(i-2) + widthl(i-2)
              do m = 1, 2
                widthl(j-1) = -widthl(j-1)
                z(j-1) = center(j-1) + widthl(j-1)
                do n = 1, 2
                  widthl(k) = -widthl(k)
                  z(k) = center(k) + widthl(k)
                  call functn ( itest, ndim, z, alpha, beta, f1 )
                  sum5 = sum5 + f1
                end do
              end do
            end do

            z(k) = center(k)

          end do

          z(j-1) = center(j-1)

        end do

        z(i-2) = center(i-2)

      end do

    end if
!
!  Compute fifth and seventh degree rules and error.
!
    rgncmp = rgnvol * ( weitp1 * sum1 &
                      + weitp2 * sum2 &
                      + weitp3 * sum3 &
                      + weitp4 * sum4 )

    rgnval = rgnvol * ( weit1 * sum1 &
                      + weit2 * sum2 &
                      + weit3 * sum3 &
                      + weit4 * sum4 &
                      + weit5 * sum5 )

    rgnerr = abs ( rgnval - rgncmp )
!
!  End basic rule.
!
    finest = finest + rgnval
    wrkstr(lenwrk) = wrkstr(lenwrk) + rgnerr
    funcls = funcls + rulcls
!
!  Place results of basic rule into partially ordered list
!  according to subregion error.
!
!  When DIVFLG = 0, start at the top of the list and move down the
!  list tree to find the correct position for the results from the
!  first half of the recently divided subregion.
!
    if ( divflg /= 1 ) then

      do

        subtmp = 2 * subrgn

        if ( sbrgns < subtmp ) then
          exit
        end if

        if ( subtmp /= sbrgns ) then
          sbtmpp = subtmp + rgnstr
          if ( wrkstr(subtmp) < wrkstr(sbtmpp) ) then
            subtmp = sbtmpp
          end if
        end if

        if ( wrkstr(subtmp) <= rgnerr ) then
          exit
        end if

        do k = 1, rgnstr
          wrkstr(subrgn-k+1) = wrkstr(subtmp-k+1)
        end do

        subrgn = subtmp

      end do
!
!  When DIVFLG = 1 start at bottom right branch and move up list
!  tree to find correct position for results from second half of
!  recently divided subregion.
!
    else

      do

        subtmp = ( subrgn / ( 2 * rgnstr ) ) * rgnstr

        if ( subtmp < rgnstr ) then
          exit
        end if

        if ( rgnerr <= wrkstr(subtmp) ) then
          exit
        end if

        do k = 1, rgnstr
          index1 = subrgn - k + 1
          index2 = subtmp - k + 1
          wrkstr(index1) = wrkstr(index2)
        end do

        subrgn = subtmp

      end do

    end if
!
!  Store results of basic rule in correct position in list.
!
    wrkstr(subrgn) = rgnerr
    wrkstr(subrgn-1) = rgnval
    wrkstr(subrgn-2) = divaxn

    do j = 1, ndim
      subtmp = subrgn - 2 * ( j + 1 )
      wrkstr(subtmp+1) = center(j)
      wrkstr(subtmp) = width(j)
    end do
!
!  When DIVFLG = 0 prepare for second application of basic rule.
!
    if ( divflg /= 1 ) then
      center(divaxo) = center(divaxo) + 2.0D+00 * width(divaxo)
      sbrgns = sbrgns + rgnstr
      subrgn = sbrgns
      divflg = 1
      cycle
    end if
!
!  End ordering and storage of basic rule results.
!  Make checks for possible termination of routine.
!
    relerr = 1.0D+00

    if ( wrkstr(lenwrk) <= 0.0D+00 ) then
      wrkstr(lenwrk) = 0.0D+00
    end if

    if ( abs ( finest ) /= 0.0D+00 ) then
      relerr = wrkstr(lenwrk) / abs ( finest )
    end if

    if ( 1.0D+00 < relerr ) then
      relerr = 1.0D+00
    end if

    if ( lenwrk < sbrgns + rgnstr + 2 ) then
      ifail = 2
    end if

    if ( maxpts < funcls + 2 * rulcls ) then
      ifail = 1
    end if

    if ( relerr < rel_tol .and. minpts <= funcls ) then
      ifail = 0
    end if

    if ( ifail < 3 ) then
      minpts = funcls
      wrkstr(lenwrk-1) = sbrgns
      exit
    end if
!
!  Prepare to use basic rule on each half of subregion with largest
!  error.
!
    divflg = 0
    subrgn = rgnstr
    wrkstr(lenwrk) = wrkstr(lenwrk) - wrkstr(subrgn)
    finest = finest - wrkstr(subrgn-1)
    divaxo = int ( wrkstr(subrgn-2) )

    do j = 1, ndim
      subtmp = subrgn - 2 * ( j + 1 )
      center(j) = wrkstr(subtmp+1)
      width(j) = wrkstr(subtmp)
    end do

    width(divaxo) = width(divaxo) / 2.0D+00
    center(divaxo) = center(divaxo) - width(divaxo)

  end do

  return
end
subroutine genz_function ( indx, ndim, z, alpha, beta, f )

!*****************************************************************************80
!
!! GENZ_FUNCTION evaluates one of the test integrand functions.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    Original FORTRAN77 version by Alan Genz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration:
!    Recent Developments, Software and Applications,
!    edited by Patrick Keast, Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX, the index of the test function.
!
!    Input, integer ( kind = 4 ) NDIM, the spatial dimension.
!
!    Input, real ( kind = 8 ) Z(NDIM), the point at which the integrand
!    is to be evaluated.
!
!    Input, real ( kind = 8 ) ALPHA(NDIM), BETA(NDIM), parameters
!    associated with the integrand function.
!
!    Output, real ( kind = 8 ) F, the value of the test function.
!
  implicit none

  integer ( kind = 4 ) ndim

  real ( kind = 8 ) alpha(ndim)
  real ( kind = 8 ) beta(ndim)
  real ( kind = 8 ) f
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323844D+00
  real ( kind = 8 ) total
  real ( kind = 8 ) value
  real ( kind = 8 ) z(ndim)

  value = 0.0D+00
!
!  Oscillatory.
!
  if ( indx == 1 ) then

    total = 2.0D+00 * pi * beta(1) + sum ( z(1:ndim) )
    value = cos ( total )
!
!  Product Peak.
!
  else if ( indx == 2 ) then

    total = product ( &
      1.0D+00 / alpha(1:ndim)**2 + ( z(1:ndim) - beta(1:ndim) )**2 )

    value = 1.0D+00 / total
!
!  Corner Peak.
!
  else if ( indx == 3 ) then
!
!  For this case, the BETA's are used to randomly select
!  a corner for the peak.
!
    total = 1.0D+00
    do j = 1, ndim
      if ( beta(j) < 0.5D+00 ) then
        total = total + z(j)
      else
        total = total + alpha(j) - z(j)
      end if
    end do
    value = 1.0D+00 / total**( ndim + 1 )
!
!  Gaussian.
!
!  The G95 compiler fails with a floating exception if the argument
!  to EXP is too small, hence we have to do some tedious handholding.
!
  else if ( indx == 4 ) then

    total = sum ( ( alpha(1:ndim) * ( z(1:ndim) - beta(1:ndim) ) )**2 )

    if ( log ( epsilon ( total ) ) < -total ) then
      value = exp ( - total )
    else
      value = 0.0D+00
    end if
!
!  C0 Function.
!
  else if ( indx == 5 ) then

    total = dot_product ( alpha(1:ndim), abs ( z(1:ndim) - beta(1:ndim) ) )
    value = exp ( - total )
!
!  Discontinuous.
!
  else if ( indx == 6 ) then

    if ( any ( beta(1:ndim) < z(1:ndim) ) ) then

      value = 0.0D+00

    else

      total = dot_product ( alpha(1:ndim), z(1:ndim) )

      value = exp ( total )

    end if

  end if

  f = value

  return
end
function genz_integral ( indx, ndim, a, b, alpha, beta )

!*****************************************************************************80
!
!! GENZ_INTEGRAL computes the exact integrals of the test functions.
!
!  Modified:
!
!    26 May 2007
!
!  Author:
!
!    Original FORTRAN77 version by Alan Genz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration:
!    Recent Developments, Software and Applications,
!    edited by Patrick Keast, Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX, the index of the test.
!
!    Input, integer ( kind = 4 ) NDIM, the spatial dimension.
!
!    Input, real ( kind = 8 ) A(NDIM), B(NDIM), the lower and upper limits
!    of integration.
!
!    Input, real ( kind = 8 ) ALPHA(NDIM), BETA(NDIM), parameters
!    associated with the integrand function.
!
!    Output, real ( kind = 8 ) GENZ_INTEGRAL, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) ndim

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) ab
  real ( kind = 8 ) alpha(ndim)
  real ( kind = 8 ) b(ndim)
  real ( kind = 8 ) beta(ndim)
  real ( kind = 8 ) genz_integral
  real ( kind = 8 ) genz_phi
  integer ( kind = 4 ) ic(ndim)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isum
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323844D+00
  integer ( kind = 4 ) rank
  real ( kind = 8 ) s
  real ( kind = 8 ) sgndm
  real ( kind = 8 ) total
  real ( kind = 8 ) value
!
!  Oscillatory
!
  if ( indx == 1 ) then

    value = 0.0D+00
!
!  Generate all sequences of NDIM 0's and 1's.
!
    rank = 0

    do

      call tuple_next ( 0, 1, ndim, rank, ic )

      if ( rank == 0 ) then
        exit
      end if

      total = 2.0D+00 * pi * beta(1)
      do j = 1, ndim
        if ( ic(j) /= 1 ) then
          total = total + alpha(j)
        end if
      end do

      isum = sum ( ic(1:ndim) )

      s = 1 + 2 * ( ( isum / 2 ) * 2 - isum )

      if ( mod ( ndim, 2 ) == 0 ) then
        value = value + s * cos ( total )
      else
        value = value + s * sin ( total )
      end if

    end do

    if ( 1 < mod ( ndim, 4 ) ) then
      value = - value
    end if
!
!  Product Peak.
!
  else if ( indx == 2 ) then

    value = 1.0D+00

    do j = 1, ndim
      value = value * alpha(j) * ( &
          atan ( ( 1.0D+00 - beta(j) ) * alpha(j) ) &
        + atan (           + beta(j)   * alpha(j) ) )
    end do
!
!  Corner Peak.
!
  else if ( indx == 3 ) then

    value = 0.0D+00

    sgndm = 1.0D+00
    do j = 1, ndim
      sgndm = - sgndm / real ( j, kind = 8 )
    end do

    rank = 0

    do

      call tuple_next ( 0, 1, ndim, rank, ic )

      if ( rank == 0 ) then
        exit
      end if

      total = 1.0D+00

      do j = 1, ndim
        if ( ic(j) /= 1 ) then
          total = total + alpha(j)
        end if
      end do

      isum = sum ( ic(1:ndim) )

      s = 1 + 2 * ( ( isum / 2 ) * 2 - isum )
      value = value + s / total

    end do

    value = value * sgndm
!
!  Gaussian.
!
  else if ( indx == 4 ) then

    value = 1.0D+00

    ab = sqrt ( 2.0D+00 )
    do j = 1, ndim
      value = value * ( sqrt ( pi ) / alpha(j) ) * &
        (   genz_phi ( ( 1.0D+00 - beta(j) ) * ab * alpha(j) ) &
          - genz_phi (           - beta(j)   * ab * alpha(j) ) )
    end do
!
!  C0 Function.
!
  else if ( indx == 5 ) then

    value = 1.0D+00
    do j = 1, ndim
      ab = alpha(j) * beta(j)
      value = value * &
        ( 2.0D+00 - exp ( - ab ) - exp ( ab - alpha(j) ) ) &
        / alpha(j)
    end do
!
!  Discontinuous.
!
  else if ( indx == 6 ) then

    value = 1.0D+00
    do j = 1, ndim
      value = value * ( exp ( alpha(j) * beta(j) ) - 1.0D+00 ) &
        / alpha(j)
    end do

  end if

  genz_integral = value

  return
end
function genz_name ( indx )

!*****************************************************************************80
!
!! GENZ_NAME returns the name of a Genz test integrand.
!
!  Modified:
!
!    26 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX, the index of the test integrand.
!
!    Output, character ( len = 13 ) GENZ_NAME, the name of the test integrand.
!
  implicit none

  integer ( kind = 4 ) indx
  character ( len = 13 ) genz_name

  if ( indx == 1 ) then
    genz_name = 'Oscillatory  '
  else if ( indx == 2 ) then
    genz_name = 'Product Peak '
  else if ( indx == 3 ) then
    genz_name = 'Corner Peak  '
  else if ( indx == 4 ) then
    genz_name = 'Gaussian     '
  else if ( indx == 5 ) then
    genz_name = 'C0 Function  '
  else if ( indx == 6 ) then
    genz_name = 'Discontinuous'
  else
    genz_name = '?????????????'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  GENZ_NAME - Fatal error!'
    write ( *, '(a)' ) '  1 <= INDX <= 6 is required.'
    stop
  end if

  return
end
function genz_phi ( z )

!*****************************************************************************80
!
!! GENZ_PHI estimates the normal cumulative density function.
!
!  Discussion:
!
!    The approximation is accurate to 1.0E-07.
!
!    This routine is based upon algorithm 5666 for the error function,
!    from Hart et al.
!
!  Modified:
!
!    13 March 2007
!
!  Author:
!
!    Original FORTRAN77 version by Alan Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Z, a value which can be regarded as the distance,
!    in standard deviations, from the mean.
!
!    Output, real ( kind = 8 ) GENZ_PHI, the integral of the normal PDF
!    from negative infinity to Z.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) ROOTPI, despite the name, is actually the
!    square root of TWO * pi.
!
  implicit none

  real ( kind = 8 ) expntl
  real ( kind = 8 ) genz_phi
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: p0 = 220.2068679123761D+00
  real ( kind = 8 ), parameter :: p1 = 221.2135961699311D+00
  real ( kind = 8 ), parameter :: p2 = 112.0792914978709D+00
  real ( kind = 8 ), parameter :: p3 = 33.91286607838300D+00
  real ( kind = 8 ), parameter :: p4 = 6.373962203531650D+00
  real ( kind = 8 ), parameter :: p5 = 0.7003830644436881D+00
  real ( kind = 8 ), parameter :: p6 = 0.03526249659989109D+00
  real ( kind = 8 ), parameter :: q0 = 440.4137358247522D+00
  real ( kind = 8 ), parameter :: q1 = 793.8265125199484D+00
  real ( kind = 8 ), parameter :: q2 = 637.3336333788311D+00
  real ( kind = 8 ), parameter :: q3 = 296.5642487796737D+00
  real ( kind = 8 ), parameter :: q4 = 86.78073220294608D+00
  real ( kind = 8 ), parameter :: q5 = 16.06417757920695D+00
  real ( kind = 8 ), parameter :: q6 = 1.755667163182642D+00
  real ( kind = 8 ), parameter :: q7 = 0.08838834764831844D+00
  real ( kind = 8 ), parameter :: rootpi = 2.506628274631001D+00
  real ( kind = 8 ) z
  real ( kind = 8 ) zabs

  zabs = abs ( z )
!
!  12 < |Z|.
!
  if ( 12.0D+00 < zabs ) then

    p = 0.0D+00

  else
!
!  |Z| <= 12
!
    expntl = exp ( - zabs * zabs / 2.0D+00 )
!
!  |Z| < 7
!
    if ( zabs < 7.0D+00 ) then

      p = expntl * (((((( &
                  p6 &
         * zabs + p5 ) &
         * zabs + p4 ) &
         * zabs + p3 ) &
         * zabs + p2 ) &
         * zabs + p1 ) &
         * zabs + p0 ) / ((((((( &
                  q7 &
         * zabs + q6 ) &
         * zabs + q5 ) &
         * zabs + q4 ) &
         * zabs + q3 ) &
         * zabs + q2 ) &
         * zabs + q1 ) &
         * zabs + q0 )
!
!  CUTOFF <= |Z|
!
    else

      p = expntl / ( &
        zabs + 1.0D+00 / ( &
        zabs + 2.0D+00 / ( &
        zabs + 3.0D+00 / ( &
        zabs + 4.0D+00 / ( &
        zabs + 0.65D+00 ))))) / rootpi

    end if

  end if

  if ( 0.0D+00 < z ) then
    p = 1.0D+00 - p
  end if

  genz_phi = p

  return
end
function genz_random ( seed )

!*****************************************************************************80
!
!! GENZ_RANDOM is a portable random number generator
!
!  Modified:
!
!    21 March 2007
!
!  Author:
!
!    Original FORTRAN77 version by Linus Schrage.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Linus Schrage,
!    A More Portable Fortran Random Number Generator,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 2, June 1979, pages 132-138.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) GENZ_RANDOM is a pseudorandom value.
!
  implicit none

  integer ( kind = 4 ), parameter :: a = 16807
  integer ( kind = 4 ), parameter :: b15 = 32768
  integer ( kind = 4 ), parameter :: b16 = 65536
  integer ( kind = 4 ) fhi
  real ( kind = 8 ) genz_random
  integer ( kind = 4 ) k
  integer ( kind = 4 ) leftlo
  integer ( kind = 4 ), parameter :: p = 2147483647
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) xalo
  integer ( kind = 4 ) xhi

  xhi = seed / b16
  xalo = ( seed - xhi * b16 ) * a
  leftlo = xalo / b16
  fhi = xhi * a + leftlo
  k = fhi / b15

  seed = ( &
           ( &
             ( xalo - leftlo * b16 ) - p &
           ) &
         + ( fhi - k * b15 ) * b16 &
         ) + k

  if ( seed < 0 ) then
    seed = seed + p
  end if

  genz_random = real ( seed, kind = 8 ) / real ( p, kind = 8 )

  return
end
subroutine multst ( nsamp, tstlim, tstfns, tstmax, difclt, expnts, ndiml, &
  ndims, sbname, subrtn, rel_tol, maxpts )

!*****************************************************************************80
!
!! MULTST tests a multidimensional integration routine.
!
!  Discussion:
!
!    The routine uses the Genz test integrand functions, with
!    the user selecting the particular subset of test integrands,
!    the set of difficulty factors, and the spatial dimensions.
!
!  Modified:
!
!    21 March 2007
!
!  Author:
!
!    Original FORTRAN77 version by Alan Genz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration:
!    Recent Developments, Software and Applications,
!    edited by Patrick Keast, Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NSAMP, the number of samples.
!    1 <= NSAMP.
!
!    Input, integer ( kind = 4 ) TSTLIM, the number of test integrands.
!
!    Input, integer ( kind = 4 ) TSTFNS(TSTLIM), the indices of the test 
!    integrands.  Each index is between 1 and 6.
!
!    Input, integer ( kind = 4 ) TSTMAX, the number of difficulty levels 
!    to be tried.
!
!    Input, real ( kind = 8 ) DIFCLT(TSTMAX), difficulty levels.
!
!    Input, real ( kind = 8 ) EXPNTS(TSTMAX), the difficulty exponents.
!
!    Input, integer ( kind = 4 ) NDIML, the number of sets of variable sizes.
!
!    Input, integer ( kind = 4 ) NDIMS(NDIML), the number of variables 
!    for the integrals in each test.
!
!    Input, character ( len = 6 ) SBNAME, the name of the integration
!    subroutine to be tested.
!
!    Input, external SUBRTN, the integration subroutine to be tested.
!
!    Input, real ( kind = 8 ) REL_TOL, the relative error tolerance.
!
!    Input, integer ( kind = 4 ) MAXPTS, the maximum number of integrand calls
!    for all tests.
!
  implicit none

  integer ( kind = 4 ), parameter :: mxtsfn = 6
  integer ( kind = 4 ) nsamp
  integer ( kind = 4 ) tstlim
  integer ( kind = 4 ) tstmax

  real ( kind = 8 ), allocatable, dimension ( : ) :: a
  real ( kind = 8 ), allocatable, dimension ( : ) :: alpha
  real ( kind = 8 ), allocatable, dimension ( : ) :: b
  real ( kind = 8 ), allocatable, dimension ( : ) :: beta
  real ( kind = 8 ) callsa(mxtsfn,mxtsfn)
  real ( kind = 8 ) callsb(mxtsfn,mxtsfn)
  real ( kind = 8 ) concof
  real ( kind = 8 ) dfact
  real ( kind = 8 ) dfclt
  real ( kind = 8 ) difclt(tstmax)
  integer ( kind = 4 ) digits
  real ( kind = 8 ) errest
  real ( kind = 8 ) errlog
  real ( kind = 8 ) ersacb(mxtsfn,mxtsfn)
  real ( kind = 8 ) ersact(mxtsfn,mxtsfn)
  real ( kind = 8 ) ersdsb(mxtsfn,mxtsfn)
  real ( kind = 8 ) ersdsc(mxtsfn,mxtsfn)
  real ( kind = 8 ) ersesb(mxtsfn,mxtsfn)
  real ( kind = 8 ) ersest(mxtsfn,mxtsfn)
  real ( kind = 8 ) ersrel(mxtsfn,mxtsfn)
  real ( kind = 8 ) estlog
  real ( kind = 8 ) exn
  real ( kind = 8 ) expnts(tstmax)
  real ( kind = 8 ) expons(mxtsfn)
  real ( kind = 8 ) finest
  external genz_function
  real ( kind = 8 ) genz_integral
  character ( len = 13 ) genz_name
  real ( kind = 8 ) genz_random
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idfclt(mxtsfn)
  integer ( kind = 4 ) ifail
  integer ( kind = 4 ) ifails
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lenwrk
  integer ( kind = 4 ) maxpts
  real ( kind = 8 ) medacb(mxtsfn)
  real ( kind = 8 ) medacb_med(3)
  real ( kind = 8 ) medact(nsamp)
  real ( kind = 8 ) medact_med(3)
  real ( kind = 8 ) medcla(mxtsfn)
  real ( kind = 8 ) medcla_med(3)
  real ( kind = 8 ) medclb(mxtsfn)
  real ( kind = 8 ) medclb_med(3)
  real ( kind = 8 ) medcls(nsamp)
  real ( kind = 8 ) medcls_med(3)
  real ( kind = 8 ) meddsb(mxtsfn)
  real ( kind = 8 ) meddsb_med(3)
  real ( kind = 8 ) meddsc(nsamp)
  real ( kind = 8 ) meddsc_med(3)
  real ( kind = 8 ) medesb(mxtsfn)
  real ( kind = 8 ) medesb_med(3)
  real ( kind = 8 ) medest(nsamp)
  real ( kind = 8 ) medest_med(3)
  real ( kind = 8 ) medrel
  real ( kind = 8 ) medrll(nsamp)
  real ( kind = 8 ) medrll_med(3)
  integer ( kind = 4 ) minpts
  integer ( kind = 4 ) n
  character ( len = 13 ) name
  integer ( kind = 4 ) nconf
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) ndiml
  integer ( kind = 4 ) ndims(ndiml)
  integer ( kind = 4 ) ndimv
  real ( kind = 8 ) qality
  real ( kind = 8 ) qallty(nsamp)
  real ( kind = 8 ) qallty_med(3)
  real ( kind = 8 ) qualty(mxtsfn,mxtsfn)
  integer ( kind = 4 ) rcalsa
  integer ( kind = 4 ) rcalsb
  real ( kind = 8 ) rel_tol
  real ( kind = 8 ) relerr
  integer ( kind = 4 ) rulcls
  character ( len = 6 )  sbname
  integer ( kind = 4 ) seed
  real ( kind = 8 ) small
  external subrtn
  real ( kind = 8 ) tactrb(mxtsfn)
  real ( kind = 8 ) tactrb_med(3)
  real ( kind = 8 ) tactrs(mxtsfn)
  real ( kind = 8 ) tactrs_med(3)
  real ( kind = 8 ) tcalsa(mxtsfn)
  real ( kind = 8 ) tcalsa_med(3)
  real ( kind = 8 ) tcalsb(mxtsfn)
  real ( kind = 8 ) tcalsb_med(3)
  real ( kind = 8 ) terdsb(mxtsfn)
  real ( kind = 8 ) terdsb_med(3)
  real ( kind = 8 ) terdsc(mxtsfn)
  real ( kind = 8 ) terdsc_med(3)
  real ( kind = 8 ) testrb(mxtsfn)
  real ( kind = 8 ) testrb_med(3)
  real ( kind = 8 ) testrs(mxtsfn)
  real ( kind = 8 ) testrs_med(3)
  real ( kind = 8 ) tqualt(mxtsfn)
  real ( kind = 8 ) tqualt_med(3)
  real ( kind = 8 ) total
  real ( kind = 8 ) trelib(mxtsfn)
  real ( kind = 8 ) trelib_med(3)
  integer ( kind = 4 ) tstfns(tstlim)
  real ( kind = 8 ) value
  real ( kind = 8 ), allocatable, dimension ( : ) :: wrkstr
!
!  Initialize and compute confidence coefficient.
!
  concof = 0.0D+00
  nconf = max ( 1, ( 2 * nsamp ) / 5 - 2 )

  do i = 1, nconf
    concof = 1.0D+00 + real ( nsamp - nconf + i, kind = 8 ) * concof &
      / real ( nconf - i + 1, kind = 8 )
  end do

  concof = 1.0D+00 - concof / real ( 2**( nsamp - 1 ), kind = 8 )

  seed = 123456

  small = epsilon ( small )

  idfclt(1:tstlim) = difclt(tstfns(1:tstlim))
  expons(1:tstlim) = expnts(tstfns(1:tstlim))
!
!  Begin main loop for different numbers of variables.
!
  do ndimv = 1, ndiml

    ndim = ndims(ndimv)

    allocate ( a(1:ndim) )
    allocate ( alpha(1:ndim) )
    allocate ( b(1:ndim) )
    allocate ( beta(1:ndim) )

    if ( ndim <= 15 ) then
      rulcls = 2**ndim + 2 * ndim**2 + 2 * ndim + 1
    else
      rulcls = ( ndim * ( 14 - ndim * ( 6 - 4 * ndim ) ) ) / 3 + 1
    end if

    lenwrk = ( 2 * ndim + 3 ) * ( 1 + maxpts / rulcls ) / 2
    allocate ( wrkstr(1:lenwrk) )

    if ( mod ( ndimv - 1, 6 ) == 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a)' ) '  Test results with', &
        nsamp, ' samples per test'
      write ( *, '(a)' ) ' '

      write ( *, '(a,10i6)' ) '  Difficulty levels', idfclt(1:tstlim)
      write ( *, '(a,10f6.1)' ) '      Exponents    ', expons(1:tstlim)

      digits = int ( - log10 ( rel_tol ) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i3,a,i8)') '  Requested digits = ', digits, &
        ', Maximum values = ', maxpts
      write ( *, '(a)' ) ' '
      write ( *, '(2x,a,a,f5.2)' ) sbname, &
        ' tests, variable results with confidence', concof
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) &
        ' Vari-  Integrand     Correct digits   Relia-  Wrong', &
        '   Integrand   Quality Total'
      write ( *, '(a,a)' ) &
        ' ables              Estimated   Actual bility Digits', &
        '    Values             Fails'
      write ( *, '(a)' ) ' '

    end if
!
!  Begin loop for different test integrands.
!
    do it = 1, tstlim

      itest = tstfns(it)
      exn = expnts(itest)
      dfclt = difclt(itest)

      a(1:ndim) = 0.0D+00
      b(1:ndim) = 1.0D+00

      ifails = 0
      medrel = 0
!
!  Begin loop for different samples.
!
      do k = 1, nsamp

        ifail = 1
!
!  Choose the integrand function parameters at random.
!
        do n = 1, ndim
          alpha(n) = genz_random ( seed )
          beta(n) = genz_random ( seed )
        end do
!
!  Modify ALPHA to account for difficulty parameter.
!
        total = sum ( alpha(1:ndim) )
        dfact = total * real ( ndim, kind = 8 )**exn / dfclt
        alpha(1:ndim) = alpha(1:ndim) / dfact
!
!  For tests 1 and 3, we modify the value of B.
!
        if ( itest == 1 .or. itest == 3 ) then
          b(1:ndim) = alpha(1:ndim)
        end if
!
!  For test 6, we modify the value of BETA.
!
        if ( itest == 6 ) then
          beta(3:ndim) = 1.0D+00
        end if
!
!  Get the exact value of the integral.
!
        value = genz_integral ( itest, ndim, a, b, alpha, beta )
!
!  Call the integration subroutine.
!
        minpts = 4 * 2**ndim

        call subrtn ( ndim, a, b, minpts, maxpts, genz_function, rel_tol, &
          itest, alpha, beta, lenwrk, wrkstr, errest, finest, ifail )

        relerr = abs ( ( finest - value ) / value )
        ifails = ifails + min ( ifail, 1 )
        relerr = max ( min ( 1.0D+00, relerr ), small )
        errlog = max ( 0.0D+00, -log10 ( relerr ) )
        errest = max ( min ( 1.0D+00, errest ), small )
        estlog = max ( 0.0D+00, -log10 ( errest ) )
        meddsc(k) = max ( 0.0D+00, estlog - errlog )
        medest(k) = estlog
        medact(k) = errlog
        medcls(k) = minpts

        if ( relerr <= errest ) then
          medrel = medrel + 1
        end if

      end do
!
!  End loop for different samples and compute medians.
!
      call r8vec_median ( nsamp, medest, medest_med )
      call r8vec_median ( nsamp, medact, medact_med )
      call r8vec_median ( nsamp, medcls, medcls_med )
      call r8vec_median ( nsamp, meddsc, meddsc_med )

      medrel = medrel / real ( nsamp, kind = 8 )

      trelib(it) = medrel

      tactrs(it) = medact_med(2)
      testrs(it) = medest_med(2)
      terdsc(it) = meddsc_med(2)
      tcalsa(it) = medcls_med(2)

      tcalsb(it) = medcls_med(3)
      tactrb(it) = medact_med(3)
      testrb(it) = medest_med(3)
      terdsb(it) = meddsc_med(3)

      ersrel(itest,ndimv) = medrel

      ersest(itest,ndimv) = medest_med(2)
      ersact(itest,ndimv) = medact_med(2)
      ersdsc(itest,ndimv) = meddsc_med(2)

      ersesb(itest,ndimv) = medest_med(3)
      ersacb(itest,ndimv) = medact_med(3)
      ersdsb(itest,ndimv) = meddsc_med(3)

      callsa(itest,ndimv) = medcls_med(2)

      callsb(itest,ndimv) = medcls_med(3)

      qality = 0.0D+00

      if ( medcls_med(1) /= 0.0D+00 ) then
        qality = ( medact_med(1) + 1.0D+00 ) * &
          ( medest_med(1) + 1.0D+00 - meddsc_med(1) ) / log ( medcls_med(1) )
      end if

      tqualt(it) = qality
      qualty(itest,ndimv) = qality
      rcalsa = int ( medcls_med(2) )
      rcalsb = int ( medcls_med(3) )
      name = genz_name ( itest )

      write ( *, &
        '(i4,2x,a14,f4.1,f5.1,f5.1,f5.1,f5.2,f4.1,f4.1,i7,i8,f6.2,i5)' ) &
        ndim, name, medest_med(2), medest_med(3), medact_med(2), &
        medact_med(3), medrel, meddsc_med(2), meddsc_med(3), rcalsa, rcalsb, &
        qality, ifails

    end do
!
!  End loop for different test integrands.
!
    call r8vec_median ( tstlim, tactrs, tactrs_med )
    call r8vec_median ( tstlim, trelib, trelib_med )
    call r8vec_median ( tstlim, testrs, testrs_med )
    call r8vec_median ( tstlim, terdsc, terdsc_med )
    call r8vec_median ( tstlim, tactrb, tactrb_med )
    call r8vec_median ( tstlim, testrb, testrb_med )
    call r8vec_median ( tstlim, terdsb, terdsb_med )
    call r8vec_median ( tstlim, tqualt, tqualt_med )
    call r8vec_median ( tstlim, tcalsa, tcalsa_med )
    call r8vec_median ( tstlim, tcalsb, tcalsb_med )

    rcalsa = int ( tcalsa_med(1) )
    rcalsb = int ( tcalsb_med(1) )

    write ( *, &
      '(i4,2x,a14,f4.1,f5.1,f5.1,f5.1,f5.2,f4.1,f4.1,i7,i8,f6.2)' ) &
      ndim, 'Medians       ', testrs_med(1), testrb_med(1), tactrs_med(1), &
      tactrb_med(1), trelib_med(1), terdsc_med(1), terdsb_med(1), rcalsa, &
      rcalsb, tqualt_med(1)

    write ( *, '(a)' ) ' '

    deallocate ( a )
    deallocate ( alpha )
    deallocate ( b )
    deallocate ( beta )
    deallocate ( wrkstr )

  end do
!
!  End loop for different numbers of variables.
!
  if ( 1 < ndiml ) then

    write ( *, '(a)' ) ' '
    write ( *, '(6x,a,a,12i3)' ) sbname, &
      ' Test integrand medians for variables', ndims(1:ndiml)

    write ( *, '(a)' ) ' '
    write ( *, '(a,a)' ) &
      '        Integrand     Correct digits   Relia-  Wrong', &
      '   Integrand   Quality'
    write ( *, '(a,a)' ) &
      '                    Estimated   Actual bility digits', &
      '     Values'
    write ( *, '(a)' ) ' '

    do it = 1, tstlim

      itest = tstfns(it)

      medact(1:ndiml) = ersact(itest,1:ndiml)
      medest(1:ndiml) = ersest(itest,1:ndiml)
      meddsc(1:ndiml) = ersdsc(itest,1:ndiml)
      medacb(1:ndiml) = ersacb(itest,1:ndiml)
      medesb(1:ndiml) = ersesb(itest,1:ndiml)
      meddsb(1:ndiml) = ersdsb(itest,1:ndiml)
      medrll(1:ndiml) = ersrel(itest,1:ndiml)
      qallty(1:ndiml) = qualty(itest,1:ndiml)
      medcla(1:ndiml) = callsa(itest,1:ndiml)
      medclb(1:ndiml) = callsb(itest,1:ndiml)

      call r8vec_median ( ndiml, medrll, medrll_med )
      call r8vec_median ( ndiml, medact, medact_med )
      call r8vec_median ( ndiml, medest, medest_med )
      call r8vec_median ( ndiml, meddsc, meddsc_med )
      call r8vec_median ( ndiml, medacb, medacb_med )
      call r8vec_median ( ndiml, medesb, medesb_med )
      call r8vec_median ( ndiml, meddsb, meddsb_med )
      call r8vec_median ( ndiml, qallty, qallty_med )
      call r8vec_median ( ndiml, medcla, medcla_med )
      call r8vec_median ( ndiml, medclb, medclb_med )

      rcalsa = int ( medcla_med(1) )
      rcalsb = int ( medclb_med(1) )
      name = genz_name ( itest )

      write ( *, '(6x,a14,f4.1,f5.1,f5.1,f5.1,f5.2,f4.1,f4.1,i7,i8,f6.2)' ) &
        name, medest_med(1), medesb_med(1), medact_med(1), &
        medacb_med(1), medrll_med(1), meddsc_med(1), meddsb_med(1), &
        rcalsa, rcalsb, qallty_med(1)

      tactrs(it) = medact_med(1)
      testrs(it) = medest_med(1)
      terdsc(it) = meddsc_med(1)
      tactrb(it) = medacb_med(1)
      testrb(it) = medesb_med(1)
      terdsb(it) = meddsb_med(1)
      tcalsa(it) = medcla_med(1)
      tcalsb(it) = medclb_med(1)
      trelib(it) = medrll_med(1)
      tqualt(it) = qallty_med(1)

    end do

    call r8vec_median ( tstlim, tactrs, tactrs_med )
    call r8vec_median ( tstlim, testrs, testrs_med )
    call r8vec_median ( tstlim, terdsc, terdsc_med )
    call r8vec_median ( tstlim, tactrb, tactrb_med )
    call r8vec_median ( tstlim, testrb, testrb_med )
    call r8vec_median ( tstlim, terdsb, terdsb_med )
    call r8vec_median ( tstlim, trelib, trelib_med )
    call r8vec_median ( tstlim, tqualt, tqualt_med )
    call r8vec_median ( tstlim, tcalsa, tcalsa_med )
    call r8vec_median ( tstlim, tcalsb, tcalsb_med )

    rcalsa = int ( tcalsa_med(1) )
    rcalsb = int ( tcalsb_med(1) )

    write ( *, '(a,f4.1,f5.1,f5.1,f5.1,f5.2,f4.1,f4.1,i7,i8,f6.2)' ) &
      '      Global medians', testrs_med(1), testrb_med(1), tactrs_med(1), &
      tactrb_med(1), trelib_med(1), terdsc_med(1), terdsb_med(1), rcalsa, &
      rcalsb, tqualt_med(1)

    write ( *, '(a)' ) ' '

  end if

  return
end
subroutine r8vec_median ( n, r, rmed )

!*****************************************************************************80
!
!! R8VEC_MEDIAN estimates the median of an R8VEC.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    Original FORTRAN77 version by Alan Genz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the array.
!
!    Input, real ( kind = 8 ) R(N), the array to be examined.
!
!    Output, real ( kind = 8 ) RMED(3).  RMED(1) contains the median,
!    RMED(2) and RMED(3) specify the confidence interval.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) nconf
  integer ( kind = 4 ) nd
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmed(3)

  do j = 1, n

    kmax = j

    do k = j + 1, n
      if ( r(kmax) < r(k) ) then
        kmax = k
      end if
    end do

    rmax = r(kmax)
    r(kmax) = r(j)
    r(j) = rmax

  end do

  nd = n / 2

  if ( mod ( n, 2 ) == 0 ) then
    rmed(1) = ( r(nd) + r(nd+1) ) / 2.0D+00
  else
    rmed(1) = r(nd+1)
  end if

  nconf = max ( 1, ( 2 * n ) / 5 - 2 )

  rmed(2) = r(n-nconf+1)
  rmed(3) = r(nconf)

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
subroutine tuple_next ( m1, m2, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT computes the next element of a tuple space.
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between M1 and M2.  The elements are produced one at a time.
!    The first element is
!      (M1,M1,...,M1),
!    the second element is
!      (M1,M1,...,M1+1),
!    and the last element is
!      (M2,M2,...,M2)
!    Intermediate elements are produced in lexicographic order.
!
!  Example:
!
!    N = 2, M1 = 1, M2 = 3
!
!    INPUT        OUTPUT
!    -------      -------
!    Rank  X      Rank   X
!    ----  ---    -----  ---
!    0     * *    1      1 1
!    1     1 1    2      1 2
!    2     1 2    3      1 3
!    3     1 3    4      2 1
!    4     2 1    5      2 2
!    5     2 2    6      2 3
!    6     2 3    7      3 1
!    7     3 1    8      3 2
!    8     3 2    9      3 3
!    9     3 3    0      0 0
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M1, M2, the minimum and maximum entries.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input/output, integer ( kind = 4 ) RANK, counts the elements.
!    On first call, set RANK to 0.  Thereafter, the output value of RANK
!    will indicate the order of the element returned.  When there are no
!    more elements, RANK will be returned as 0.
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( m2 < m1 ) then
    rank = 0
    return
  end if

  if ( rank <= 0 ) then

    x(1:n) = m1
    rank = 1

  else

    rank = rank + 1
    i = n

    do

      if ( x(i) < m2 ) then
        x(i) = x(i) + 1
        exit
      end if

      x(i) = m1

      if ( i == 1 ) then
        rank = 0
        x(1:n) = m1
        exit
      end if

      i = i - 1

    end do

  end if

  return
end

