program main

!*****************************************************************************80
!
!! MAIN is the main program for SPECFUN_PRB1.
!
!  Discussion:
!
!    SPECFUN_PRB1 calls sample problems for the SPECFUN library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPECFUN_PRB1'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPECFUN library.'

  call alg_test ( )
  call daw_test ( )
  call ei_test ( )
  call r8_erf_test ( )
  call r8_gamma_test ( )
  call i0_test ( )
  call i1_test ( )
  call j0_test ( )
  call j1_test ( )
  call k0_test ( )
  call k1_test ( )
  call r8_psi_test ( )
  call ri_test ( )
  call rj_test ( )
  call rk_test ( )
  call ry_test ( )
  call y0_test ( )
  call y1_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPECFUN_PRB1'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine alg_test ( )

!*****************************************************************************80
!
!! ALG_TEST tests DLGAMA.
!
!  Discussion:
!
!    Accuracy tests compare function values against values
!    generated with the duplication formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Performance evaluation of programs related to the real gamma function,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 46-54.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) all9
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) cl
  real ( kind = 8 ) del
  real ( kind = 8 ) dlgama
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) half
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) p875
  real ( kind = 8 ) p3125
  real ( kind = 8 ) p625
  real ( kind = 8 ) p6875
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sixten
  real ( kind = 8 ) ten
  real ( kind = 8 ) two
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xc
  real ( kind = 8 ) xl
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) x99
  real ( kind = 8 ) xp99
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz
!
!  C1 = 0.5 - LN(SQRT(PI))
!  C2 = LN(2)
!  C3 = LN(2) - 11/16
!
  data c1 /-7.2364942924700087072D-02 /
  data c2 /6.9314718055994530942D-01/
  data c3 /5.6471805599453094172D-03/
  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two /2.0D+00/
  data ten /10.0D+00/
  data sixten /16.0D+00/
  data p6875 /0.6875D+00/
  data p875 /0.875D+00 /
  data p3125 /1.3125D+00/
  data p625 / 1.625D+00/
  data x99 / -999.0D+00 /
  data xp99 /0.99D+00/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  a = zero
  b = p875
  n = 2000
  xn = real ( n, kind = 8 )
  jt = 0
!
!  Determine largest argument for DLGAMA by iteration.
!
  cl = xp99 * xmax
  z = -cl / x99

   80 continue

  zz = cl / ( log ( z ) - one )

  if ( abs ( zz / z - one ) > two * beta * eps ) then
    z = zz
    go to 80
  end if

  cl = zz
!
!  Random argument accuracy tests.
!
  do j = 1, 3

    k1 = 0
    k3 = 0
    xc = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Use the duplication formula.
!
      if ( j /= 3 ) then

        if ( j == 1 ) then
          z = x + half
          x = z - half
          y = x + x
        else
          x = x + x
          x = x * half
          y = ( x + x ) - one
          z = x - half
        end if

        u = dlgama ( x )
        w = ( y - half ) - half
        zz = anint ( w * sixten ) / sixten
        w = w - zz
        v = ((( half - zz * p6875 ) - c1 ) - w * p6875 ) &
          - c3 * ( w + zz )
        v = ( ( v + dlgama ( y ) ) - dlgama ( z ) )

      else

        z = x * half + half
        y = z - half
        x = y + y
        u = dlgama ( x )
        v = ( c1 + ( ( x - half ) - half ) * c2 ) &
          + dlgama ( y ) + dlgama ( z ) - half

      end if
!
!  Accumulate the results.
!
      w = ( u - v ) / u

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        xc = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    k2 = n - k3 - k1
    r7 = sqrt ( r7 / xn )

    if ( j == 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) &
        'Test of LGAMA(X) vs LN(2*SQRT(PI))-2X*LN(2)', &
        '+LGAMA(2X)-LGAMA(X+1/2)'
    else if ( j == 2 ) then
      write ( *,1001)
    else
      write ( *,1002)
    end if

    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    else
      w = all9
    end if

    write ( *,1021) r6,ibeta,w,xc
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    else
      w = all9
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for the next test.
!
    a = p3125
    b = p625

    if ( j == 2 ) then
      a = two + two
      b = ten + ten
    end if

  end do
!
!  Special tests.
!  First test with special arguments.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of special arguments'
  write ( *, '(a)' ) ' '
  z = eps
  zz = dlgama ( z )
  write ( *,1041) z,zz
  z = half
  zz = dlgama ( z )
  write ( *,1041) z,zz
  z = one
  zz = dlgama ( z )
  write ( *,1041) z,zz
  z = two
  zz = dlgama ( z )
  write ( *,1041) z,zz
!
!  Test of error returns.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of Error Returns:'
  write ( *, '(a)' ) ' '
  z = xmin
  write ( *,1053) z
  zz = dlgama ( z )
  write ( *,1061) zz
  z = cl
  write ( *,1053) z
  zz = dlgama ( z )
  write ( *,1061) zz
  z = -one
  write ( *,1052) z
  zz = dlgama ( z )
  write ( *,1061) zz
  z = zero
  write ( *,1052) z
  zz = dlgama ( z )
  write ( *,1061) zz
  z = xp99 * xmax
  write ( *,1052) z
  zz = dlgama ( z )
  write ( *,1061) zz
  write ( *, '(a)' ) '  This concludes the tests.'
  return
 1001 format('1Test of LGAMA(X) vs LN(2*SQRT(PI))-(2X-1)*LN(2)+', &
      'LGAMA(X-1/2)-LGAMA(2X-1)'//)
 1002 format('1Test of LGAMA(X) vs -LN(2*SQRT(PI))+X*LN(2)+', &
      'LGAMA(X/2)+LGAMA(X/2+1/2)'//)
 1010 format(I7,' Random arguments were tested from the interval (', &
      F5.1,',',F5.1,')'//)
 1011 format('  LGAMA(X) was larger',I6,' times,'/ &
      14X,' agreed',I6,' times, and'/ &
      10X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1041 format(' LGAMA (',E13.6,') = ',E13.6//)
 1052 format(' LGAMA will be called with the argument',E13.6,/ &
      ' This should trigger an error message'//)
 1053 format(' LGAMA will be called with the argument',E13.6,/ &
      ' This should not trigger an error message'//)
 1061 format(' LGAMA returned the value',E13.6///)
end
subroutine daw_test ( )

!*****************************************************************************80
!
!! DAW_TEST tests DAW.
!
!  Discussion:
!
!    Accuracy test compare function values against a local
!    Taylor's series expansion.  Derivatives are generated
!    from the recurrence relation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) all9
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) daw
  real ( kind = 8 ) del
  real ( kind = 8 ) delta
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) forten
  real ( kind = 8 ) half
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) p(0:14)
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sixten
  real ( kind = 8 ) two
  real ( kind = 8 ) t1
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xkay
  real ( kind = 8 ) xl
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) x99
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two / 2.0d0 /
  data forten / 14.0d0 /
  data sixten / 1.6d1 /
  data x99 /-999.0d0 /
  data delta /6.25d-2 /
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  ait = real ( it, kind = 8 )
  albeta = log ( beta )
  a = delta
  b = one
  jt = 0
!
!  Random argument accuracy tests based on local Taylor expansion.
!
  do j = 1, 4

    n = 2000
    xn = real ( n, kind = 8 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Purify arguments.
!
      y = x - delta
      w = sixten * y
      t1 = w + y
      y = t1 - w
      x = y + delta
!
!  Use Taylor's Series Expansion.
!
      p(0) = daw ( y )
      z = y + y
      p(1) = one - z * p(0)
      xkay = two

      do ii = 2, 14
        p(ii) = - ( z * p(ii-1) + xkay * p(ii-2) )
        xkay = xkay + two
      end do

      zz = p(14)
      xkay = forten

      do ii = 1, 14
        zz = zz * delta / xkay + p(14-ii)
        xkay = xkay - one
      end do

      z = daw ( x )
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    k2 = n - k1 - k3
    r7 = sqrt ( r7 / xn )
    write ( *,1000)
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta
    w = x99

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    end if

    write ( *,1021) r6,ibeta,w,x1
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
    w = x99

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b
    b = b + b

    if ( j == 1 ) then
      b = b + half
    end if

  end do
!
!  Special tests.  First check values for negative arguments.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '

  write ( *,1030) ibeta

  do i = 1, 10
    x = ren ( j ) * ( two + two )
    b = daw ( x )
    a = b + daw ( -x )
    if ( a * b /= zero ) then
      a = ait + log ( abs ( a / b ) ) / albeta
    end if
    write ( *,1031) x,a
    x = x + del
  end do
!
!  Next, test with special arguments.
!
  write ( *,1040)
  z = xmin
  zz = daw ( z )
  write ( *,1041) zz
!
!  Test of error return for arguments > xmax.  First, determine xmax.
!
  if ( half < xmin * xmax ) then
    xbig = half / xmin
  else
    xbig = xmax
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of Error Returns:'
  write ( *, '(a)' ) ' '

  z = xbig * ( one - delta * delta )
  write ( *,1052) z
  zz = daw ( z )
  write ( *,1062) zz
  z = xbig
  write ( *,1053) z
  zz = daw ( z )
  write ( *,1062) zz
  w = one + delta * delta

  if ( w < xmax / xbig ) then
    z = xbig * w
    write ( *,1053) z
    zz = daw ( z )
    write ( *,1062) zz
  end if

  write ( *, '(a)' ) '  This concludes the tests.'
  return
 1000 format('1Test of Dawson''s Integral vs Taylor expansion'//)
 1010 format(I7,' Random arguments were tested from the interval ', &
      '(',F5.2,',',F5.2,')'//)
 1011 format('  F(X) was larger',I6,' times,'/ &
      10X,' agreed',I6,' times, and'/ &
      6X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1030 format(7X,'Estimated loss of base',i3,' significant digits in'// &
         8X,'X',10X,'F(x)+F(-x)'/)
 1031 format(3X,F7.3,F16.2)
 1040 format(//' Test of special arguments'//)
 1041 format('  F(XMIN) = ',E24.17/)
 1052 format(' DAW will be called with the argument',E13.6,/ &
      ' This should not underflow'//)
 1053 format(' DAW will be called with the argument',E13.6,/ &
      ' This may underflow'//)
 1062 format(' DAW returned the value',E13.6///)

end
subroutine ei_test ( )

!*****************************************************************************80
!
!! EI_TEST tests EI, EONE and EXPEI.
!
!  Discussion:
!
!    Accuracy test compare function values against local Taylor's
!    series expansions.  Derivatives for Ei(x) are generated from
!    the recurrence relation using a technique due to Gautschi.
!
!    Special argument tests are run with the
!    related functions E1(x) and exp(-x)Ei(x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
!    Walter Gautschi, BJ Klein,
!    Recursive computation of certain derivatives - A study of error
!    propagation,
!    Communications of the ACM,
!    Volume 13, 1970, pages 7-9.
!
!    Walter Gautschi, BJ Klein,
!    Remark on Algorithm 282,
!    Communications of the ACM,
!    Volume 13, 1970, pages 53-54.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c1
  real ( kind = 8 ) d(0:25)
  real ( kind = 8 ) del
  real ( kind = 8 ) dx
  real ( kind = 8 ) en
  real ( kind = 8 ) ei
  real ( kind = 8 ) eone
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) expei
  real ( kind = 8 ) fiv12
  real ( kind = 8 ) four
  real ( kind = 8 ) fourth
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  integer ( kind = 4 ) n1
  real ( kind = 8 ) one
  real ( kind = 8 ) p0625
  real ( kind = 8 ) rem
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) six
  real ( kind = 8 ) sum
  real ( kind = 8 ) ten
  real ( kind = 8 ) two
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xc
  real ( kind = 8 ) xden
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlge
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnp1
  real ( kind = 8 ) xnum
  real ( kind = 8 ) x0
  real ( kind = 8 ) x99
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero

  data zero / 0.0D+00 /
  data one / 1.0D+00 /
  data fourth / 0.25d0 /
  data four / 4.0d0/
  data six / 6.0d0/
  data ten / 10.0d0 /
  data x0 / 0.3725d0/
  data x99 / -999.0d0 /
  data p0625 / 0.0625d0/
  data fiv12 / 512.0d0/
  data rem / -7.424779065800051695596d-5/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  dx = -p0625
  a = fourth + dx
  b = x0 + dx
  n = 2000
  n1 = 25
  xn = real ( n, kind = 8 )
  jt = 0
!
!  Random argument accuracy tests.
!
  do j = 1, 8

    k1 = 0
    k3 = 0
    xc = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      y = del * ren ( jt ) + xl
      x = y - dx
      y = x + dx
!
!  Test Ei against series expansion.
!
      v = ei ( x )
      z = ei ( y )
      sum = zero
      u = x
      call dsubn ( u, n1, xmax, d )
      en = real ( n1, kind = 8 ) + one
      sum = d(n1) * dx / en

      do k = n1, 1, -1
        en = en - one
        sum = ( sum + d(k-1) ) * dx / en
      end do

      u = v + sum
!
!  Accumulate results.
!
      w = z - u
      w = w / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        xc = y
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    k2 = n - k3 - k1
    r7 = sqrt ( r7 / xn )
    write ( *,1000)
    write ( *,1010) n,a,b
    write ( *,1011) k1
    write ( *,1015) k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    else
      w = x99
    end if

    write ( *,1021) r6,ibeta,w,xc
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    else
      w = x99
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    if ( j == 1 ) then
      dx = -dx
      a = x0 + dx
      b = six
    else if ( j <= 4 ) then
      a = b
      b = b + b
    else if ( j == 5 ) then
      a = -fourth
      b = -one
    else if ( j == 6 ) then
      a = b
      b = -four
    else
      a = b
      b = -ten
    end if

  end do
!
!  Special tests.  First, check accuracy near the zero of Ei(x).
!
  write ( *,1040)
  x = ( four - one ) / ( four + four )
  y = ei ( x )
  write ( *,1041) x,y
  z = ( ( y - ( four + one ) / fiv12 ) - rem ) / y

  if ( z /= zero ) then
    w = log ( abs ( z ) ) / albeta
  else
    w = x99
  end if

  write ( *,1042) z,ibeta,w
  w = max ( ait + w, zero )
  write ( *,1022) ibeta,w
!
!  Check near XBIG, the largest argument acceptable to EONE, i.e.,
!  the negative of the smallest argument acceptable to EI.
!  Determine XBIG with Newton iteration on the equation
!  EONE(x) = XMIN.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of Error Returns:'
  write ( *, '(a)' ) ' '
  two = one + one
  v = sqrt ( eps )
  c1 = real ( minexp, kind = 8 ) * log ( beta )
  xn = -c1

  320 continue

  xnum = -xn - log ( xn ) + log ( one + one / xn ) - c1
  xden = -( xn * xn + xn + xn + two ) / ( xn * ( xn + one ) )
  xnp1 = xn - xnum / xden
  w = ( xn - xnp1 ) / xnp1

  if ( abs ( w ) > v ) then
    xn = xnp1
    go to 320
  end if

  xbig = xnp1
  x = aint ( ten * xbig ) / ten
  write ( *,1052) x
  y = eone ( x )
  write ( *,1062) y
  x = xbig * ( one + v )
  write ( *,1053) x
  y = eone ( x )
  write ( *,1062) y
!
!  Check near XMAX, the largest argument acceptable to EI.  Determine
!  XLGE with Newton iteration on the equation
!  EI(x) = XMAX.
!
  c1 = real ( maxexp, kind = 8 ) * log ( beta )
  xn = c1

  330 continue

  xnum = xn - log ( xn ) + log ( one + one / xn ) - c1
  xden = ( xn * xn - two ) / ( xn * ( xn + one ) )
  xnp1 = xn - xnum / xden
  w = ( xn - xnp1 ) / xnp1

  if ( abs ( w ) > v ) then
    xn = xnp1
    go to 330
  end if

  xlge = xnp1
  x = aint ( ten * xlge ) / ten
  write ( *,1054) x
  y = ei ( x )
  write ( *,1064) y
  x = xlge * ( one + v )
  write ( *,1055) x
  y = ei ( x )
  write ( *,1064) y
!
!  Check with XHUGE, the largest acceptable argument for EXPEI.
!
  if ( xmin * xmax <= one ) then
    x = xmax
  else
    x = one / xmin
  end if

  write ( *,1056) x
  y = expei ( x )
  write ( *,1065) y
  x = zero
  write ( *,1055) x
  y = ei ( x )
  write ( *,1064) y
  write ( *, '(a)' ) '  This concludes the tests.'
  return
 1000 format('1Test of Ei(x) vs series expansion'//)
 1010 format(I7,' Random arguments were tested from the interval (', &
      F7.3,',',F7.3,')'//)
 1011 format('     EI(X) was larger',I6,' times,')
 1015 format(14X,' agreed',I6,' times, and'/ &
      10X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1040 format(//' Test of special arguments'//)
 1041 format('   EI (',E13.6,') = ',E13.6//)
 1042 format(' The relative error is',E15.4,' = ',I4,' **',F7.2/)
 1052 format(' EONE will be called with the argument',E13.6,/ &
      ' This should not underflow'//)
 1053 format(' EONE will be called with the argument',E13.6,/ &
      ' This should underflow'//)
 1054 format(' EI will be called with the argument',E13.6,/ &
      ' This should not overflow'//)
 1055 format(' EI will be called with the argument',E13.6,/ &
      ' This should overflow'//)
 1056 format(' EXPEI will be called with the argument',E13.6,/ &
      ' This should not underflow'//)
 1062 format(' EONE returned the value',E13.6///)
 1064 format(' EI returned the value',E13.6///)
 1065 format(' EXPEI returned the value',E13.6///)
end
subroutine r8_erf_test ( )

!*****************************************************************************80
!
!! R8_ERF_TEST tests R8_ERF and related functions.
!
!  Discussion:
!
!    Accuracy test compare function values against local Taylor's
!    series expansions.  Derivatives for erfc(x) are expressed as
!    repeated integrals of erfc(x).  These are generated from the
!    recurrence relation using a technique due to Gautschi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    William Cody,
!    Performance evaluation of programs related to the error and
!    complementary error functions,
!    ACM Transactions on Mathematical Software,
!    Volume 16, Number 1, March 1990, pages 29-37.
!
!    Walter Gautschi,
!    Evaluation of the repeated integrals of the coerror function,
!    ACM Transactions on Mathematical Software,
!    Volume 3, 1977, pages 240-252.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  real ( kind = 8 ) epscon
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) ff
  real ( kind = 8 ) f0
  real ( kind = 8 ) half
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  real ( kind = 8 ) one
  real ( kind = 8 ) r
  real ( kind = 8 ) ren
  real ( kind = 8 ) r1(500)
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) r8_erf
  real ( kind = 8 ) r8_erfc
  real ( kind = 8 ) r8_erfcx
  real ( kind = 8 ) sc
  real ( kind = 8 ) sixten
  real ( kind = 8 ) thresh
  real ( kind = 8 ) ten
  real ( kind = 8 ) two
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xc
  real ( kind = 8 ) xl
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) xn1
  real ( kind = 8 ) x99
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz
!
!  C1 = 1/sqrt(pi)
!
  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two / 2.0d0 /
  data ten /10.0d0/
  data sixten /16.0d0/
  data thresh /0.46875d0/
  data x99 / -999.0d0/
  data c1 /5.6418958354775628695d-1/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  c = abs ( ait * albeta ) + ten
  a = zero
  b = thresh
  n = 2000
  xn = real ( n, kind = 8 )
  jt = 0
  n0 = ( ibeta / 2 ) * ( it + 5 ) / 6 + 4
!
!  Determine largest argument for ERFC test by Newton iteration.
!
  c2 = log ( xmin ) + log ( one / c1 )
  xbig = sqrt ( -c2 )

   50 continue

  x = xbig
  f0 = x * x
  ff = f0 + half / f0 + log ( x ) + c2
  f0 = x + x + one / x - one / ( x * f0 )
  xbig = x - ff / f0

  if ( abs ( x - xbig ) / x > ten * eps ) then
    go to 50
  end if
!
!  Random argument accuracy tests.
!
  do j = 1, 5

    k1 = 0
    k3 = 0
    xc = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Test erf against double series expansion.
!
      if ( j == 1 ) then

        f0 = real ( n0, kind = 8 )
        ff = f0 + f0 + one
        z = x * x
        w = z + z
        u = zero
        v = zero

        do k = 1, n0
          u = -z / f0 * ( one + u )
          v = w / ff * ( one + v )
          f0 = f0 - one
          ff = ff - two
        end do

        v = c1 * ( x + x ) * ( ( ( u * v + ( u + v ) ) &
          + half ) + half )
        u = r8_erf ( x )
!
!  Test erfc or scaled erfc against expansion in repeated
!  integrals of the coerror function.
!
      else

        z = x + half
        x = z - half
        r = zero

        if ( x <= one ) then
          n0 = 499
        else
          n0 = min ( 499, int ( c / ( abs ( log ( z ) ) ) ) )
        end if

        n1 = n0
        xn1 = real ( n1 + 1, kind = 8 )

        do k = 1, n0
          r = half / ( z + xn1 * r )
          r1(n1) = r
          n1 = n1 - 1
          xn1 = xn1 - one
        end do

        ff = c1 / ( z + r1(1) )

        if ( ( j / 2 ) * 2 == j ) then
          f0 = r8_erfc ( z ) * exp ( x + half * half )
          u = r8_erfc ( x )
        else
          f0 = r8_erfcx ( z )
          u = r8_erfcx ( x )
        end if

        sc = f0 / ff
!
!  Scale things to avoid premature underflow.
!
        epscon = f0
        ff = sixten * ff / eps

        do n1 = 1, n0
          ff = r1(n1) * ff
          r1(n1) = ff * sc
          if ( r1(n1) < epscon ) then
            k = n1
            go to 111
          end if
        end do

  111       continue

        v = r1(k)

        do n1 = 1, k-1
          v = v + r1(k-n1)
        end do
!
!  Remove scaling here.
!
        v = v * eps / sixten + f0

      end if
!
!  Accumulate results.
!
      w = ( u - v ) / u

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        xc = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    k2 = n - k3 - k1
    r7 = sqrt ( r7 / xn )

    if ( j == 1 ) then
      write ( *,1000)
      write ( *,1010) n,a,b
      write ( *,1011) k1
    else if ( (j/2)*2 == j ) then
      write ( *,1001)
      write ( *,1010) n,a,b
      write ( *,1012) k1
    else
      write ( *,1002)
      write ( *,1010) n,a,b
      write ( *,1013) k1
    end if

    write ( *,1015) k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    else
      w = x99
    end if

    write ( *,1021) r6,ibeta,w,xc
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    else
      w = x99
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    if ( j == 1 ) then
      a = b
      b = two
    else if ( j == 3 ) then
      a = b
      b = aint(xbig*sixten) / sixten - half
    else if ( j == 4 ) then
      b = ten + ten
    end if

  end do
!
!  Special tests.
!  First check values for negative arguments.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1030) ibeta
  x = zero
  del = -half

  do i = 1, 10

    u = r8_erf ( x )
    a = u + r8_erf ( -x )

    if ( a * u /= zero ) then
      a = ait + log ( abs ( a / u ) ) / albeta
    end if

    v = r8_erfc ( x )
    b = u + v - one

    if ( b /= zero ) then
      b = ait + log ( abs ( b ) ) / albeta
    end if

    w = r8_erfcx ( x )
    c = aint ( x * sixten ) / sixten
    r = ( x - c ) * ( x + c )
    c = ( exp ( c * c ) * exp ( r ) * v - w ) / w

    if ( c /= zero ) then
      c = max ( zero, ait + log ( abs ( c ) ) / albeta)
    end if

    write ( *,1031) x,a,b,c
    x = x + del

  end do
!
!  Next, test with special arguments.
!
  write ( *,1040)
  z = xmax
  zz = r8_erf ( z )
  write ( *,1041) z,zz
  z = zero
  zz = r8_erf ( z )
  write ( *,1041) z,zz
  zz = r8_erfc ( z )
  write ( *,1042) z,zz
  z = -xmax
  zz = r8_erfc ( z )
  write ( *,1042) z,zz
!
!  Test of error returns.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of Error Returns:'
  write ( *, '(a)' ) ' '
  w = xbig
  z = w * ( one - half * half )
  write ( *,1052) z
  zz = r8_erfc ( z )
  write ( *,1062) zz
  z = w * ( one + ten * eps )
  write ( *,1053) z
  zz = r8_erfc ( z )
  write ( *,1062) zz
  w = xmax

  if ( c1 < xmax * xmin ) then
    w = c1/xmin
  end if

  z = w * ( one - one / sixten )
  write ( *,1054) z
  zz = r8_erfcx ( z )
  write ( *,1064) zz
  w = -sqrt ( log ( xmax / two ) )
  z = w * ( one - one / ten )
  write ( *,1055) z
  zz = r8_erfcx ( z )
  write ( *,1064) zz
  z = w * ( one + ten * eps )
  write ( *,1056) z
  zz = r8_erfcx ( z )
  write ( *,1064) zz
  write ( *, '(a)' ) '  This concludes the tests.'
  return

 1000 format('1Test of erf(x) vs double series expansion'//)
 1001 format(///' Test of erfc(x) vs exp(x+1/4) SUM i^n erfc(x+1/2)'//)
 1002 format('1Test of exp(x*x) erfc(x) vs SUM i^n erfc(x+1/2)'//)
 1010 format(I7,' Random arguments were tested from the interval (', &
      F7.3,',',F7.3,')'//)
 1011 format('    ERF(X) was larger',I6,' times,')
 1012 format('   ERFC(X) was larger',I6,' times,')
 1013 format('  ERFCX(X) was larger',I6,' times,')
 1015 format(14X,' agreed',I6,' times, and'/ &
      10X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1030 format(7X,'Estimated loss of base',i3,'significant digits in'// &
         3X,'X',5X,'Erf(x)+Erf(-x)',3X,'Erf(x)+Erfc(x)-1', &
         3X,'Erfcx(x)-exp(x*x)*erfc(x)'/)
 1031 format(F7.3,3F16.2)
 1040 format(//' Test of special arguments'//)
 1041 format('   ERF (',E13.6,') = ',E13.6//)
 1042 format('  ERFC (',E13.6,') = ',E13.6//)
 1052 format(' ERFC will be called with the argument',E13.6,/ &
      ' This should not underflow'//)
 1053 format(' ERFC will be called with the argument',E13.6,/ &
      ' This may underflow'//)
 1054 format(' ERFCX will be called with the argument',E13.6,/ &
      ' This should not underflow'//)
 1055 format(' ERFCX will be called with the argument',E13.6,/ &
      ' This should not overflow'//)
 1056 format(' ERFCX will be called with the argument',E13.6,/ &
      ' This may overflow'//)
 1062 format(' ERFC returned the value',E13.6///)
 1064 format(' ERFCX returned the value',E13.6///)
end
subroutine r8_gamma_test ( )

!*****************************************************************************80
!
!! R8_GAMMA_TEST tests R8_GAMMA.
!
!  Discussion:
!
!    Accuracy tests compare function values against values
!    generated with the duplication formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Performance evaluation of programs related to the real gamma function,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 46-54.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) alnx
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) cl
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) half
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  integer ( kind = 4 ) nx
  real ( kind = 8 ) one
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) ren
  real ( kind = 8 ) ten
  real ( kind = 8 ) two
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xc
  real ( kind = 8 ) xl
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xminv
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xxn
  real ( kind = 8 ) xp
  real ( kind = 8 ) xph
  real ( kind = 8 ) x99
  real ( kind = 8 ) xp99
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data c1 /2.8209479177387814347d-1/
  data c2 /9.1893853320467274178d-1/
  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two / 2.0D+00 /
  data ten /10.0d0/
  data x99 / -999.0d0 /
  data xp99 / 0.99d0 /
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  a = zero
  b = two
  n = 2000
  xn = real ( n, kind = 8 )
  jt = 0
!
!  Determine smallest argument for GAMMA.
!
  if ( xmin * xmax < one ) then
    xminv = one / xmax
  else
    xminv = xmin
  end if
!
!  Determine largest argument for GAMMA by Newton iteration.
!
  cl = log ( xmax )
  xp = half * cl
  cl = c2 - cl

   50 continue

  x = xp
  alnx = log ( x )
  xnum = ( x - half ) * alnx - x + cl
  xp = x - xnum / ( alnx - half / x )

  if ( abs ( xp - x ) / x >= ten * eps ) then
    go to 50
  end if

  cl = xp
!
!  Random argument accuracy tests.
!
  do j = 1, 4

    k1 = 0
    k3 = 0
    xc = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Use duplication formula for X not close to the zero.
!
      xph = x * half + half
      xp = xph - half
      x = xp + xp
      nx = int ( x )
      xxn = real ( nx, kind = 8 )
      c = ( two**nx ) * ( two**( x - xxn ) )
      z = r8_gamma ( x )
      zz = ( ( c * c1 ) * r8_gamma ( xp ) ) * r8_gamma ( xph )
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        xc = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    k2 = n - k3 - k1
    r7 = sqrt ( r7 / xn )
    write ( *,1000)
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    else
      w = x99
    end if

    write ( *,1021) r6,ibeta,w,xc
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    else
      w = x99
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b

    if ( j == 1 ) then
      b = ten
    else if ( j == 2 ) then
      b = cl - half
    else
      a = -( ten - half ) * half
      b = a + half
    end if

  end do
!
!  Special tests.
!
!  First test with special arguments.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1040)
  x = -half
  y = r8_gamma ( x )
  write ( *,1041) x, y
  x = xminv / xp99
  y = r8_gamma ( x )
  write ( *,1041) x, y
  x = one
  y = r8_gamma ( x )
  write ( *,1041) x, y
  x = two
  y = r8_gamma ( x )
  write ( *,1041) x, y
  x = cl * xp99
  y = r8_gamma ( x )
  write ( *,1041) x, y
!
!  Test of error returns.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of Error Returns:'
  write ( *, '(a)' ) ' '
  x = -one
  write ( *,1052) x
  y = r8_gamma ( x )
  write ( *,1061) y
  x = zero
  write ( *,1052) x
  y = r8_gamma ( x )
  write ( *,1061) y
  x = xminv * ( one - eps )
  write ( *,1052) x
  y = r8_gamma ( x )
  write ( *,1061) y
  x = cl * ( one + eps )
  write ( *,1052) x
  y = r8_gamma ( x )
  write ( *,1061) y
  write ( *, '(a)' ) '  This concludes the tests.'
  return

 1000 format('1Test of GAMMA(X) vs Duplication Formula'//)
 1010 format(I7,' Random arguments were tested from the interval (', &
      F7.3,',',F7.3,')'//)
 1011 format(' GAMMA(X) was larger',I6,' times,'/ &
      13X,' agreed',I6,' times, and'/ &
      9X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1040 format(//' Test of special arguments'//)
 1041 format(' GAMMA (',E13.6,') = ',E13.6//)
 1052 format(' GAMMA will be called with the argument',E13.6,/ &
      ' This should trigger an error message'//)
 1061 format(' GAMMA returned the value',E13.6///)
end
subroutine i0_test ( )

!*****************************************************************************80
!
!! I0_TEST tests BESI0 and BESEI0.
!
!  Discussion:
!
!    Accuracy tests compare function values against values
!    generated with the multiplication formula for small
!    arguments and values generated from a Taylor's Series
!    Expansion using Amos' Ratio Scheme for initial values
!    for large arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Donald Amos,
!    Computation of Modified Bessel Functions and Their Ratios,
!    Mathematics of Computation,
!    Volume 28, Number 24, January 1974.
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) ak
  real ( kind = 8 ) akk
  real ( kind = 8 ) albeta
  real ( kind = 8 ) arr(8,6)
  real ( kind = 8 ) ateten
  real ( kind = 8 ) b
  real ( kind = 8 ) besei0
  real ( kind = 8 ) besi0
  real ( kind = 8 ) beta
  real ( kind = 8 ) bot
  real ( kind = 8 ) c
  real ( kind = 8 ) const
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) delta
  real ( kind = 8 ) e
  real ( kind = 8 ) eight
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) f
  real ( kind = 8 ) half
  real ( kind = 8 ) hund
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) mborg
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) mb2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) ovrchk
  real ( kind = 8 ) one28
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) temp
  real ( kind = 8 ) top
  real ( kind = 8 ) two
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u(560)
  real ( kind = 8 ) u2(560)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xbad
  real ( kind = 8 ) xj1
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xlarge
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnine
  real ( kind = 8 ) x1
  real ( kind = 8 ) x99
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two / 2.0d0/
  data eight /8.0d0/
  data xnine /9.0d0 /
  data sixten /1.6d1/
  data ateten /1.8d1/
  data hund /1.0d2/
  data one28 /1.28d2/
  data x99 /-999.0D+00/
  data xlam /1.03125d0/
  data xlarge /1.0d4/
  data c /0.9189385332d0/

  data arr /0.0d0,1.0d0,-1.0d0,1.0d0,-2.0d0,1.0d0,-3.0d0,1.0d0, &
            -999.0d0,-999.0d0,-999.0d0,3.0d0,-12.0d0,9.0d0,-51.0d0, &
            18.0d0,-5040.0d0,720.0d0,0.0d0,-999.0d0,-999.0d0, &
            60.0d0,-360.0d0,345.0d0,-1320.0d0,192.0d0,-120.0d0, &
            24.0d0,0.0d0,-999.0d0,-999.0d0,2520.0d0,-96.0d0,15.0d0, &
            -33.0d0,7.0d0,-6.0d0,2.0d0,0.0d0,-999.0d0,-4.0d0,1.0d0, &
            -3.0d0,1.0d0,-2.0d0,1.0d0,-1.0d0,1.0d0/
!
!  Statement functions.
!
  top(x) = x - half*log(x) + log(one+(one/eight-xnine/ &
     one28/x)/x)
  bot(x) = -(sixten*x+ateten) / (((one28*x+sixten)*x+ &
     xnine)*x) + one - half/x
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  ait = real ( it, kind = 8 )
  albeta = log ( beta )
  a = zero
  b = two
  jt = 0
  const = c + log ( xmax )
  delta = xlam - one
  f = ( xlam - one ) * ( xlam + one ) * half
!
!  Random argument accuracy tests.
!
  do j = 1, 4
!
!  Calculate the number of terms needed for convergence of the
!  series by using Newton's iteration on the asymptotic form of
!  the multiplication theorem.
!
    xbad = b
    d = ait * albeta - c + one
    e = log ( xbad * f ) + one
    akk = one

  100   continue

    ak = akk

    z = d + e * ak - ( ak + half ) * log ( ak + one )
    zz = e - ( ak + half ) / ( ak + one ) - log ( ak + one )
    akk = ak - z / zz

    if ( abs ( ak - akk ) > hund * eps * ak ) then
      go to 100
    end if

    mborg = int ( akk ) + 1
    n = 2000
    xn = real ( n, kind = 8 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Carefully purify arguments.
!
      if ( j == 1 ) then
        y = x / xlam
      else
        y = x - delta
      end if

      temp = sixten * y
      t1 = temp + y
      t1 = temp + t1
      y = t1 - temp
      y = y - temp

      if ( j == 1 ) then
        x = y * xlam
      else
        x = y + delta
      end if
!
!  Use Amos's Ratio Scheme.
!
      d = f * y
      mb = mborg + mborg
      mb2 = mb - 1
      xmb = real ( mb2, kind = 8 )
      temp = ( xmb + one + half ) * ( xmb + one + half )
      u2(mb) = y / ( xmb + half + sqrt ( temp + y * y ) )
!
!  Generate ratios using recurrence.
!
      do ii = 2, mb
        ovrchk = xmb / ( y * half )
        u2(mb2) = one / ( ovrchk + u2(mb2+1) )
        xmb = xmb - one
        mb2 = mb2 - 1
      end do

      u(1) = besi0 ( y )

      if ( j == 1 ) then
!
!  Accuracy test is based on the multiplication theorem.
!
        mb = mb - mborg

        do ii = 2, mb
          u(ii) = u(ii-1) * u2(ii-1)
        end do
!
!  Accurate summation.
!
        mb = mb - 1
        xmb = real ( mb, kind = 8 )
        sum = u(mb+1)
        ind = mb

        do ii = 2, mb
          sum = sum * d / xmb + u(ind)
          ind = ind - 1
          xmb = xmb - one
        end do

        zz = sum * d + u(ind)

      else
!
!  Accuracy test is based on Taylor's Series Expansion.
!
        u(2) = u(1) * u2(1)
        mb = 8
        j1 = mb
        xj1 = real ( j1 + 1, kind = 8 )
        iexp = 0
!
!  Accurate summation.
!
        do ii = 1, mb

          j2 = 1

  160         continue

          j2 = j2 + 1

          if ( arr(j1,j2) /= x99 ) then
            go to 160
          end if

          j2 = j2 - 1
          t1 = arr(j1,j2)
          j2 = j2 - 1
!
!  Group I0 terms in the derivative.
!
          if ( j2 == 0 ) then
            go to 168
          end if

  165         continue

          t1 = t1 / ( y * y ) + arr(j1,j2)
          j2 = j2 - 1

          if ( j2 >= 1 ) then
            go to 165
          end if

  168         continue

          if ( iexp == 1 ) then
            t1 = t1 / y
          end if

          j2 = 6

  170         continue

          j2 = j2 - 1
          if ( arr(ii,j2) /= x99 ) then
            go to 170
          end if

          j2 = j2 + 1
          t2 = arr(ii,j2)
          j2 = j2 + 1

          if ( iexp == 0 ) then
            iexp = 1
          else
            iexp = 0
          end if
!
!  Group I1 terms in the derivative.
!
          if ( j2 == 7 ) then
            go to 177
          end if

  175         continue

          t2 = t2 / ( y * y ) + arr(ii,j2)
          j2 = j2 + 1
          if ( j2 <= 6 ) then
            go to 175
          end if

  177         continue

          if ( iexp == 1 ) then
            t2 = t2 / y
          end if

          if ( j1 == 8 ) then
            sum = u(1) * t1 + u(2) * t2
          else
            sum = sum * ( delta / xj1 ) + ( u(1) * t1 + u(2) * t2 )
          end if

          j1 = j1 - 1
          xj1 = real ( j1 + 1, kind = 8 )

        end do

        zz = sum * delta + u(1)

      end if

      z = besi0 ( x )
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    n = k1 + k2 + k3
    r7 = sqrt ( r7 / xn )

    if ( j == 1 ) then
      write ( *,1000)
    else
      write ( *,1001)
    end if

    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    else
      w = x99
    end if

    write ( *,1021) r6,ibeta,w,x1
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    else
      w = x99
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b
    b = b + b

    if ( j == 1 ) then
      b = b + b - half
    end if

  end do
!
!  Test of error returns.
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1031)
  y = besi0 ( xmin )
  write ( *,1032) y
  y = besi0 ( zero )
  write ( *,1033) 0,y
  x = -one * ren ( jt )
  y = besi0 ( x )
  write ( *,1034) x,y
  x = -x
  y = besi0(x)
  write ( *,1034) x,y
  y = besei0 ( xmax )
  write ( *,1035) y
!
!  Determine largest safe argument for unscaled functions.
!
  write ( *, 1036)
  xa = log ( xmax )

  330 continue

  xb = xa - ( top ( xa ) - const ) / bot ( xa )

  if ( eps < abs ( xb - xa ) / xb ) then
    xa = xb
    go to 330
  end if

  xlarge = xb / xlam
  y = besi0 ( xlarge )
  write ( *,1034) xlarge,y
  xlarge = xb * xlam
  y = besi0 ( xlarge )
  write ( *,1034) xlarge,y
  write ( *, 1037)

  return
 1000 format('1Test of I0(X) vs Multiplication Theorem'//)
 1001 format('1Test of I0(X) vs Taylor series'//)
 1010 format(I7,' Random arguments were tested from the interval ', &
      '(',F5.2,',',F5.2,')'//)
 1011 format(' I0(X) was larger',I6,' times,'/ &
      10X,' agreed',I6,' times, and'/ &
      6X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1031 format(' Test with extreme arguments'/)
 1032 format(' I0(XMIN) = ',E24.17/)
 1033 format(' I0(',I1,') = ',E24.17/)
 1034 format(' I0(',E24.17,' ) = ',E24.17/)
 1035 format(' E**-X * I0(XMAX) = ',E24.17/)
 1036 format(' Tests near the largest argument for unscaled functions'/)
 1037 format(' This concludes the tests.')
end
subroutine i1_test ( )

!*****************************************************************************80
!
!! I1_TEST tests BESI1 and BESEI1.
!
!  Discussion:
!
!    Accuracy tests compare function values against values
!    generated with the multiplication formula for small
!    arguments and values generated from a Taylor's Series
!    Expansion using Amos's Ratio Scheme for initial values
!    for large arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Donald Amos,
!    Computation of Modified Bessel Functions and Their Ratios,
!    Mathematics of Computation,
!    Volume 28, Number 24, January 1974.
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) ak
  real ( kind = 8 ) akk
  real ( kind = 8 ) albeta
  real ( kind = 8 ) arr(8,7)
  real ( kind = 8 ) b
  real ( kind = 8 ) besei1
  real ( kind = 8 ) besi1
  real ( kind = 8 ) beta
  real ( kind = 8 ) bot
  real ( kind = 8 ) c
  real ( kind = 8 ) const
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) delta
  real ( kind = 8 ) e
  real ( kind = 8 ) eight
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) f
  real ( kind = 8 ) half
  real ( kind = 8 ) hund
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) mborg
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) mb2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) ovrchk
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) temp
  real ( kind = 8 ) three
  real ( kind = 8 ) top
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u(560)
  real ( kind = 8 ) u2(560)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xbad
  real ( kind = 8 ) xj1
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xlarge
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) x99
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data three / 3.0d0/
  data eight /8.0d0/
  data sixten /1.6d1/
  data hund /1.0d2/
  data x99 /-999.0d0 /
  data xlam /1.03125d0/
  data xlarge /1.0d4/
  data c / 0.9189385332d0/

  data  arr /1.0d0,-1.0d0,1.0d0,-2.0d0,1.0d0,-3.0d0,1.0d0,-4.0d0, &
            -999.0d0,-999.0d0,3.0d0,-12.0d0,9.0d0,-51.0d0,18.0d0, &
            -132.0d0,40320.0d0,0.0d0,-999.0d0,-999.0d0,60.0d0, &
            -360.0d0,345.0d0,-2700.0d0,10440.0d0,-5040.0d0,720.0d0, &
            0.0d0,-999.0d0,-999.0d0,2520.0d0,-20160.0d0,729.0d0, &
            -1320.0d0,192.0d0,-120.0d0,24.0d0,0.0d0,-999.0d0, &
            -999.0d0,26.0d0,-96.0d0,15.0d0,-33.0d0,7.0d0,-6.0d0, &
            2.0d0,0.0d0,1.0d0,-4.0d0,1.0d0,-3.0d0,1.0d0,-2.0d0, &
            1.0d0,-1.0d0 /
!
!  Statement functions.
!
  top(x) = x - half*log(x) + log(one-three/(eight*x))

  bot(x) = three / ((eight*x-three)*x) + one - half/x
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  ait = real ( it, kind = 8 )
  albeta = log ( beta )
  a = zero
  b = one
  jt = 0
  const = c + log ( xmax )
  delta = xlam - one
  f = ( xlam - one ) * ( xlam + one ) * half
!
!  Random argument accuracy tests.
!
  do j = 1, 4
!
!  Calculate the number of terms needed for convergence of the
!  series by using Newton's iteration on the asymptotic form of
!  the multiplication theorem.
!
    xbad = b
    d = ait * albeta - c + one
    e = log ( xbad * f ) + one
    akk = one

  100   continue

    ak = akk
    z = d + e * ak - ( ak + half ) * log ( ak + one )
    zz = e - ( ak + half ) / ( ak + one ) - log ( ak + one )
    akk = ak - z / zz

    if ( abs ( ak - akk ) > hund * eps * ak ) then
      go to 100
    end if

    mborg = int ( akk ) + 1
    n = 2000
    xn = real ( n, kind = 8 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Carefully purify arguments.
!
      if ( j == 1 ) then
        y = x / xlam
      else
        y = x - delta
      end if

      w = sixten * y
      t1 = w + y
      t1 = w + t1
      y = t1 - w
      y = y - w

      if ( j == 1 ) then
        x = y * xlam
      else
        x = y + delta
      end if
!
!  Use Amos's Ratio Scheme.
!
      d = f * y
      mb = mborg + mborg
      mb2 = mb - 1
      xmb = real ( mb2, kind = 8 )
      temp = ( xmb + one + half ) * ( xmb + one + half )
      u2(mb) = y / ( xmb + half + sqrt ( temp + y * y ) )
!
!  Generate ratios using recurrence.
!
      do ii = 2, mb
        ovrchk = xmb / ( y * half )
        u2(mb2) = one / ( ovrchk + u2(mb2+1) )
        xmb = xmb - one
        mb2 = mb2 - 1
      end do

      u(2) = besi1 ( y )
      u(1) = u(2) / u2(1)
!
!  Accuracy test is based on the multiplication theorem.
!
      if ( j == 1 ) then

        mb = mb - mborg

        do ii = 3, mb
          u(ii) = u(ii-1) * u2(ii-1)
        end do
!
!  Accurate summation.
!
        mb = mb - 1
        xmb = real ( mb - 1, kind = 8 )
        sum = u(mb+1)
        ind = mb

        do ii = 2, mb
          sum = sum * d / xmb + u(ind)
          ind = ind - 1
          xmb = xmb - one
        end do

        zz = xlam * sum

      else
!
!  Accuracy test is based on Taylor's Series Expansion.
!
        mb = 8
        j1 = mb
        xj1 = real ( j1 + 1, kind = 8 )
        iexp = 1
!
!  Accurate summation.
!
        do ii = 1, mb

          j2 = 1

  160         continue

          j2 = j2 + 1

          if ( arr(j1,j2) /= x99 ) then
            go to 160
          end if

          j2 = j2 - 1
          t1 = arr(j1,j2)
          j2 = j2 - 1
!
!  Group I0 terms in the derivative.
!
          if ( j2 == 0 ) then
            go to 168
          end if

  165         continue

          t1 = t1 / ( y * y ) + arr(j1,j2)
          j2 = j2 - 1

          if ( j2 >= 1 ) then
            go to 165
          end if

  168         continue

          if ( iexp == 1 ) then
            t1 = t1 / y
          end if

          j2 = 7

  170         continue

          j2 = j2 - 1

          if ( arr(ii,j2) /= x99 ) then
            go to 170
          end if

          j2 = j2 + 1
          t2 = arr(ii,j2)
          j2 = j2 + 1

          if ( iexp == 0 ) then
            iexp = 1
          else
            iexp = 0
          end if
!
!  Group I1 terms in the derivative.
!
          if ( j2 == 8 ) then
            go to 177
          end if

  175         continue

          t2 = t2 / ( y * y ) + arr(ii,j2)
          j2 = j2 + 1

          if ( j2 <= 7 ) then
            go to 175
          end if

  177         continue

          if ( iexp == 1 ) then
            t2 = t2 / y
          end if

          if ( j1 == 8 ) then
            sum = u(1) * t1 + u(2) * t2
          else
            sum = sum * ( delta / xj1 ) + ( u(1) * t1 + u(2) * t2 )
          end if

          j1 = j1 - 1
          xj1 = real ( j1 + 1, kind = 8 )

        end do

        zz = sum * delta + u(2)

      end if

      z = besi1 ( x )
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    n = k1 + k2 + k3
    r7 = sqrt ( r7 / xn )

    if ( j == 1 ) then
      write ( *,1000)
    else
      write ( *,1001)
    end if

    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    else
      w = x99
    end if

    write ( *,1021) r6,ibeta,w,x1
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    else
      w = x99
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b
    b = b + b

    if ( j == 1 ) then
      b = b + b + three + half
    end if

  end do
!
!  Test of error returns.
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1031)
  y = besi1 ( xmin )
  write ( *,1032) y
  y = besi1 ( zero )
  write ( *,1033) 0,y
  x = -one * ren ( jt )
  y = besi1 ( x )
  write ( *,1034) x,y
  x = -x
  y = besi1 ( x )
  write ( *,1034) x,y
  y = besei1 ( xmax )
  write ( *,1035) y
!
!  Determine largest safe argument for unscaled functions.
!
  write ( *, 1036)
  xa = log ( xmax )

  330 continue

  xb = xa - ( top ( xa ) - const ) / bot ( xa )

  if ( eps < abs ( xb - xa ) / xb ) then
    xa = xb
    go to 330
  end if

  xlarge = xb / xlam
  y = besi1 ( xlarge )
  write ( *,1034) xlarge,y
  xlarge = xb * xlam
  y = besi1 ( xlarge )
  write ( *,1034) xlarge,y
  write ( *, 1037)
  return

 1000 format('1Test of I1(X) vs Multiplication Theorem'//)
 1001 format('1Test of I1(X) vs Taylor series'//)
 1010 format(I7,' Random arguments were tested from the interval ', &
      '(',F5.2,',',F5.2,')'//)
 1011 format(' I1(X) was larger',I6,' times,'/ &
      10X,' agreed',I6,' times, and'/ &
      6X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1031 format(' Test with extreme arguments'/)
 1032 format(' I1(XMIN) = ',E24.17/)
 1033 format(' I1(',I1,') = ',E24.17/)
 1034 format(' I1(',E24.17,' ) = ',E24.17/)
 1035 format(' E**-X * I1(XMAX) = ',E24.17/)
 1036 format(' Tests near the largest argument for unscaled functions'/)
 1037 format(' This concludes the tests.')
end
subroutine j0_test ( )

!*****************************************************************************80
!
!! J0_TEST tests BESJ0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) all9
  real ( kind = 8 ) b
  real ( kind = 8 ) besj0
  real ( kind = 8 ) besj1
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj0
  real ( kind = 8 ) bj0p(6,10)
  real ( kind = 8 ) bj1
  real ( kind = 8 ) bj1p(6,10)
  real ( kind = 8 ) cj0
  real ( kind = 8 ) cj1
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) delta
  real ( kind = 8 ) eight
  real ( kind = 8 ) elev
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) four
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) term
  real ( kind = 8 ) twenty
  real ( kind = 8 ) two56
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xi(2)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xm
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) yinv
  real ( kind = 8 ) ysq
  real ( kind = 8 ) yx(2)
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz
!
!  Mathematical constants
!
  data zero / 0.0D+00 /
  data one / 1.0D+00 /
  data four / 4.0D+00 /
  data delta /0.0625d0/
  data eight /8.0d0 /
  data twenty /20.0d0 /
  data all9 /-999.0d0 /
  data two56 /256.0d0/
  data sixten /16.0d0 /
  data elev /11.0d0/
!
!  Coefficients for Taylor expansion
!
  data bj0p / 0.0d0,1.8144d6,-2.3688d5,1.0845d4,-2.7d2,5.0d0, &
            -1.8144d5,2.394d4,-1.125d3,30.0d0,-1.0d0,0.0d0, &
             0.0d0,2.016d4,-2.7d3,1.32d2,-4.0d0,0.0d0, &
            -2.52d3,3.45d2,-18.0d0,1.0d0,0.0d0,0.0d0, &
             0.0d0,3.6d2,-51.0d0,3.0d0,0.0d0,0.0d0, &
            -60.0d0,9.0d0,-1.0d0,0.0d0,0.0d0,0.0d0, &
             0.0d0,12.0d0,-2.0d0,0.0d0,0.0d0,0.0d0, &
            -3.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
             0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            -1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/

  data bj1p/-3.6288d6,9.2736d5,-6.201d4,1.965d3,-40.0d0,1.0d0, &
             3.6288d5,-9.324d4,6.345d3,-2.1d2,5.0d0,0.0d0, &
            -4.032d4,1.044d4,-7.29d2,26.0d0,-1.0d0,0.0d0, &
             5.04d3,-1.32d3,96.0d0,-4.0d0,0.0d0,0.0d0, &
            -7.2d2,1.92d2,-15.0d0,1.0d0,0.0d0,0.0d0, &
             1.2d2,-33.0d0,3.0d0,0.0d0,0.0d0,0.0d0, &
            -24.0d0,7.0d0,-1.0d0,0.0d0,0.0d0,0.0d0, &
             6.0d0,-2.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            -2.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
             1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
!
!  Zeroes of J0
!
  data xi /616.0d0,1413.0d0/
  data yx(1)/-7.3927648221700192757d-4/, &
       yx(2)/-1.8608651797573879013d-4/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  jt = 0
  b = zero
!
!  Random argument accuracy tests based on a Taylor expansion.
!
  do j = 1, 3

    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    n = 2000
    a = b

    if ( j == 1 ) then
      b = four
    else if ( j == 2 ) then
      b = eight
    else
      b = twenty
      n = 2000
    end if

    xn = real ( n, kind = 8 )
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Carefully purify arguments and evaluate identities.
!
      y = x - delta
      w = sixten * y
      y = ( w + y ) - w
      x = y + delta
      sum = zero
      term = zero
      bj0 = besj0 ( y )
      z = besj0 ( x )
      d = delta

      if ( abs ( z ) < abs ( bj0 ) ) then
        cj0 = x
        x = y
        y = cj0
        cj0 = bj0
        bj0 = z
        z = cj0
        d = -d
      end if

      bj1 = besj1 ( y )
      yinv = one / y
      ysq = one / ( y * y )
      xm = elev
!
!  Evaluate (12-II)th derivative at Y.
!
      do ii = 1, 10

        cj0 = bj0p(1,ii)
        cj1 = bj1p(1,ii)
        jj = ( 11 - ii ) / 2 + 1

        do iii = 2, jj
          cj0 = cj0 * ysq + bj0p(iii,ii)
          cj1 = cj1 * ysq + bj1p(iii,ii)
        end do

        if ( ( ii / 2 ) * 2 /= ii ) then
          cj0 = cj0 * yinv
        else
          cj1 = cj1 * yinv
        end if

        term = cj0 * bj0 + cj1 * bj1
        sum = ( sum + term ) * d / xm
        xm = xm - one

      end do

      sum = ( sum - bj1 ) * d + bj0
      zz = sum
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Process and output statistics.
!
    n = k1 + k2 + k3
    r7 = sqrt ( r7 / xn )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Test of J0(X) vs Taylor expansion'
    write ( *, '(a)' ) ' '
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    else
      w = all9
    end if

    write ( *,1021) r6,ibeta,w,x1
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    else
      w = all9
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

  end do
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1031) ibeta

  do i = 1, 2

    x = xi(i) / two56
    y = besj0 ( x )
    t = ( y - yx(i) ) / yx(i)

    if ( t /= zero ) then
      w = log ( abs ( t ) ) / albeta
    else
      w = all9
    end if

    w = max ( ait + w, zero )
    write ( *,1032) x,y,w

  end do
!
!  Test of error returns.
!
  write ( *,1033)
  x = xmax
  write ( *,1034) x
  y = besj0 ( x )
  write ( *,1036) y
  write ( *, '(a)' ) '  This concludes the tests.'
  return

 1010 format(I7,' random arguments were tested from the interval ', &
   1H(,F5.1,1H,,F5.1,1H)//)
 1011 format(' ABS(J0(X)) was larger',I6,' times', / &
       15X,' agreed',I6,' times, and'/ &
     11X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number.' //)
 1021 format(' The maximum relative error of',E15.4,3H = ,I4,3H **, &
    F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
    ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      3H = ,I4,3H **,F7.2)
 1031 format(' Accuracy near zeros'//10X,'X',15X,'BESJ0(X)', &
      13X,'Loss of base',I3,' digits'/)
 1032 format(E20.10,E25.15,8X,F7.2/)
 1033 format(//' Test with extreme arguments'///)
 1034 format(' J0 will be called with the argument ',E17.10/ &
       ' This may stop execution.'//)
 1036 format(' J0 returned the value',E25.17/)
end
subroutine j1_test ( )

!*****************************************************************************80
!
!! J1_TEST tests BESJ1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) all9
  real ( kind = 8 ) b
  real ( kind = 8 ) besj0
  real ( kind = 8 ) besj1
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj0
  real ( kind = 8 ) bj0p(6,10)
  real ( kind = 8 ) bj1
  real ( kind = 8 ) bj1p(6,10)
  real ( kind = 8 ) cj0
  real ( kind = 8 ) cj1
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) delta
  real ( kind = 8 ) eight
  real ( kind = 8 ) elev
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) four
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) term
  real ( kind = 8 ) thrten
  real ( kind = 8 ) twenty
  real ( kind = 8 ) two
  real ( kind = 8 ) two56
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xi(2)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xm
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) yinv
  real ( kind = 8 ) ysq
  real ( kind = 8 ) yx(2)
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz
!
!  Mathematical constants
!
  data zero   / 0.0D+00 /
  data one    / 1.0D+00 /
  data four   / 4.0D+00 /
  data delta  / 0.0625d0/
  data eight  / 8.0d0 /
  data twenty / 20.0d0 /
  data all9   /-999.0d0/
  data two56  / 256.0d0/
  data two    / 2.0d0 /
  data thrten / 13.0d0/
  data sixten / 16.0d0/
  data elev   / 11.0d0/
!
!  Coefficients for Taylor expansion
!
  data bj0p/1.99584d7,-2.58552d6,1.16235d5,-2.775d3,4.5d1,-1.0d0, &
             0.0d0,-1.8144d6,2.3688d5,-1.0845d4,2.7d2,-5.0d0, &
             1.8144d5,-2.394d4,1.125d3,-30.0d0,1.0d0,0.0d0, &
             0.0d0,-2.016d4,2.7d3,-1.32d2,4.0d0,0.0d0, &
             2.52d3,-3.45d2,18.0d0,-1.0d0,0.0d0,0.0d0, &
             0.0d0,-3.6d2,51.0d0,-3.0d0,0.0d0,0.0d0, &
             60.0d0,-9.0d0,1.0d0,0.0d0,0.0d0,0.0d0, &
             0.0d0,-12.0d0,2.0d0,0.0d0,0.0d0,0.0d0, &
             3.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
             0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
  data bj1p/-3.99168d7,1.016064d7,-6.7095d5,2.067d4,-3.9d2,6.0d0, &
             3.6288d6,-9.2736d5,6.201d4,-1.965d3,40.0d0,-1.0d0, &
            -3.6288d5,9.324d4,-6.345d3,2.1d2,-5.0d0,0.0d0, &
             4.032d4,-1.044d4,7.29d2,-26.0d0,1.0d0,0.0d0, &
            -5.04d3,1.32d3,-96.0d0,4.0d0,0.0d0,0.0d0, &
             7.2d2,-1.92d2,15.0d0,-1.0d0,0.0d0,0.0d0, &
            -1.2d2,33.0d0,-3.0d0,0.0d0,0.0d0,0.0d0, &
             24.0d0,-7.0d0,1.0d0,0.0d0,0.0d0,0.0d0, &
            -6.0d0,2.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
             2.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
!
!  Zeroes of J1
!
  data xi/981.0d0,1796.0d0/

  data yx(1)/-1.3100393001327972376d-4/, &
       yx(2)/ 1.1503460702301698285d-5/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  jt = 0
  b = zero
!
!  Random argument accuracy tests based on a Taylor expansion.
!
  do j = 1, 4

    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    n = 2000
    a = b

    if ( j == 1 ) then
      b = one
    else if ( j == 2 ) then
      b = four
    else if ( j == 3 ) then
      b = eight
    else
      b = twenty
    end if

    xn = real ( n, kind = 8 )
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Use traditional Maclaurin series for small arguments.
!
      if ( j == 1 ) then

        y = x / two
        sum = y
        xm = thrten

        do ii = 1, 12
          sum = sum * y / xm
          xm = xm - one
          sum = ( one - sum / xm ) * y
        end do

        zz = sum
        z = besj1 ( x )

      else
!
!  Use local Taylor series elsewhere.  First, purify arguments.
!
        y = x - delta
        w = sixten * y
        y = ( w + y ) - w
        x = y + delta
        sum = zero
        term = zero
        bj1 = besj1 ( y )
        z = besj1 ( x )
        d = delta

        if ( abs ( z ) < abs ( bj1 ) ) then
          cj1 = x
          x = y
          y = cj1
          cj1 = bj1
          bj1 = z
          z = cj1
          d = -d
        end if

        bj0 = besj0 ( y )
        yinv = one / y
        ysq = one / ( y * y )
        xm = elev
!
!  Evaluate (12-II)th derivative at Y.
!
        do ii = 1, 10

          cj0 = bj0p(1,ii)
          cj1 = bj1p(1,ii)
          jj = ( 12 - ii ) / 2 + 1

          do iii = 2, jj
            cj0 = cj0 * ysq + bj0p(iii,ii)
            cj1 = cj1 * ysq + bj1p(iii,ii)
          end do

          if ( ( ii / 2 ) * 2 == ii ) then
            cj0 = cj0 * yinv
          else
            cj1 = cj1 * yinv
          end if

          term = cj0 * bj0 + cj1 * bj1
          sum = ( sum + term ) * d / xm
          xm = xm - one

        end do

        sum = ( sum + bj0 - bj1 * yinv ) * d + bj1
        zz = sum

      end if
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Process and output statistics.
!
    n = k1 + k2 + k3
    r7 = sqrt ( r7 / xn )

    if ( j == 1 ) then
      write ( *,1000)
    else
      write ( *,1001)
    end if

    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    else
      w = all9
    end if

    write ( *,1021) r6,ibeta,w,x1
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    else
      w = all9
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

  end do
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1031) ibeta

  do i = 1, 2
    x = xi(i) / two56
    y = besj1 ( x )
    t = ( y - yx(i) ) / yx(i)

    if ( t /= zero ) then
      w = log ( abs ( t ) ) / albeta
    else
      w = all9
    end if

    w = max ( ait + w, zero )
    write ( *,1032) x,y,w

  end do
!
!  Test of error returns.
!
  write ( *,1033)
  x = xmax
  write ( *,1034) x
  y = besj1 ( x )
  write ( *,1036) y
  write ( *, '(a)' ) '  This concludes the tests.'
  return

 1000 format('1Test of J1(X) VS Maclaurin expansion'  //)
 1001 format('1Test of J1(X) VS local Taylor expansion'  //)
 1010 format(I7,' random arguments were tested from the interval ', &
   1H(,F5.1,1H,,F5.1,1H)//)
 1011 format(' ABS(J1(X)) was larger',I6,' times', / &
       15X,' agreed',I6,' times, and'/ &
     11X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number.' //)
 1021 format(' The maximum relative error of',E15.4,3H = ,I4,3H **, &
    F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
    ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      3H = ,I4,3H **,F7.2)
 1031 format(' Accuracy near zeros'//10X,'X',15X,'BESJ1(X)', &
      13X,'Loss of base',I3,' digits'/)
 1032 format(E20.10,E25.15,8X,F7.2/)
 1033 format(//' Test with extreme arguments'///)
 1034 format(' J1 will be called with the argument ',E17.10/ &
       ' This may stop execution.'//)
 1036 format(' J1 returned the value',E25.17/)
end
subroutine k0_test ( )

!*****************************************************************************80
!
!! K0_TEST tests BESK0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  User defined functions:
!
!    BOT, TOP
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) all9
  real ( kind = 8 ) amaxexp
  real ( kind = 8 ) ateten
  real ( kind = 8 ) b
  real ( kind = 8 ) besek0
  real ( kind = 8 ) besk0
  real ( kind = 8 ) besk1
  real ( kind = 8 ) beta
  real ( kind = 8 ) bot
  real ( kind = 8 ) c
  real ( kind = 8 ) const
  real ( kind = 8 ) del
  real ( kind = 8 ) eight
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) five
  real ( kind = 8 ) half
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) oneneg
  real ( kind = 8 ) one28
  real ( kind = 8 ) pi
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  logical sflag
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  logical tflag
  real ( kind = 8 ) top
  real ( kind = 8 ) twenty
  real ( kind = 8 ) two
  real ( kind = 8 ) u(0:559)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xden
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xlarge
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnine
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two / 2.0d0 /
  data eight /8.0d0/
  data xnine /9.0d0 /
  data ateten /18.0d0 /
  data twenty /20.0d0/
  data one28 /128.0d0/
  data five / 5.0d0 /
  data oneneg / -1.0d0 /
  data xden / 16.0d0 /
  data all9 /-999.0d0/
  data pi /3.141592653589793d0/

  top(x) = -x - half*log(two*x) + log(one-(one/eight-xnine/ &
     one28/x)/x)

  bot(x) = (xden*x-ateten) / (((one28*x-xden)*x+xnine)*x) &
     - one - half/x
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  amaxexp = real ( maxexp, kind = 8 )
  jt = 0
  b = eps
  xlam = ( xden - one ) / xden
  const = half * log ( pi ) - log ( xmin )
!
!  Random argument accuracy tests.
!
  do j = 1, 3

    sflag = ( j == 1 .and. amaxexp / ait <= five )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    n = 2000

    if ( sflag ) then
      b = sqrt ( eps )
    end if

    a = b

    if ( j == 1 ) then
      b = one
    else if ( j == 2 ) then
      b = eight
    else
      b = twenty
    end if

    xn = real ( n, kind = 8 )
    del = (b - a) / xn
    xl = a
!
!  Accuracy test is based on the multiplication theorem.
!
    do i = 1, n

      x = del * ren ( jt ) + xl
      y = x / xlam
      w = xden * y
      y = ( w + y ) - w
      x = y * xlam
      u(0) = besk0 ( y )
      u(1) = besk1 ( y )
      tflag = sflag .and. ( y < half )

      if ( tflag ) then
        u(0) = u(0) * eps
        u(1) = u(1) * eps
      end if

      mb = 1
      xmb = one
      y = y * half
      t = u(0) * eps
      w = ( one - xlam ) * ( one + xlam )
      c = w * y

      do ii = 2, 60

        t = xmb * t / c
        z = u(ii-1)

        if ( z < t ) then
          go to 120
        else if ( u(ii-1) > one ) then
          if ( xmb / y > ( xmax / u(ii-1) ) ) then
            xl = xl + del
            a = xl
            go to  200
          end if
        end if

        u(ii) = xmb / y * u(ii-1) + u(ii-2)
        xmb = xmb + one
        mb = mb + 1

      end do

  120     continue

      sum = u(mb)
      ind = mb

      do ii = 1, mb
        ind = ind - 1
        sum = sum * w * y / xmb + u(ind)
        xmb = xmb - one
      end do

      zz = sum

      if ( tflag ) then
        zz = zz / eps
      end if

      z = besk0 ( x )
      y = z

      if ( u(0) > y ) then
        y = u(0)
      end if

      w = ( z - zz ) / y

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

  200     continue

    end do

    n = k1 + k2 + k3
    xn = real ( n, kind = 8 )
    r7 = sqrt ( r7 / xn )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Test of K0(X) vs Multiplication Theorem'
    write ( *, '(a)' ) ' '
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta
    w = all9

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    end if

    if ( j == 3 ) then
      write ( *,1024) r6,ibeta,w,x1
    else
      write ( *,1021) r6,ibeta,w,x1
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
    w = all9

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    end if

    if ( j == 3 ) then
      write ( *,1025) r7,ibeta,w
    else
      write ( *,1023) r7,ibeta,w
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

  end do
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1031)
  y = besk0 ( xmin )
  write ( *,1032) y
  y = besk0 ( zero )
  write ( *,1033) 0,y
  x = ren ( jt ) * oneneg
  y = besk0 ( x )
  write ( *,1034) x,y
  y = besek0 ( xmax )
  write ( *,1035) y
  xa = log ( xmax )

  330 continue

  xb = xa - ( top ( xa ) + const ) / bot ( xa )

  if ( eps < abs ( xb - xa ) / xb ) then
    xa = xb
    go to 330
  end if

  xlarge = xb * xlam
  y = besk0 ( xlarge )
  write ( *,1034) xlarge,y
  xlarge = xb * ( xnine / eight )
  y = besk0 ( xlarge )
  write ( *,1034) xlarge,y
  return


 1010 format(I7,' random arguments were tested from the interval (', &
      F5.1,',',F5.1,')'//)
 1011 format(' ABS(K0(X)) was larger',I6,' times,'/ &
      20X,' agreed',I6,' times, and'/ &
      16X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number.'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1024 format(' The maximum absolute error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1025 format(' The root mean square absolute error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1031 format(//' Test with extreme arguments'/)
 1032 format(' K0(XMIN) = ',E24.17/)
 1033 format(' K0(',I1,') = ',E24.17/)
 1034 format(' K0(',E24.17,' ) = ',E24.17/)
 1035 format(' E**X * K0(XMAX) = ',E24.17/)
end
subroutine k1_test ( )

!*****************************************************************************80
!
!! K1_TEST tests BESK1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) all9
  real ( kind = 8 ) amaxexp
  real ( kind = 8 ) ateten
  real ( kind = 8 ) b
  real ( kind = 8 ) besek1
  real ( kind = 8 ) besk0
  real ( kind = 8 ) besk1
  real ( kind = 8 ) beta
  real ( kind = 8 ) bot
  real ( kind = 8 ) c
  real ( kind = 8 ) const
  real ( kind = 8 ) del
  real ( kind = 8 ) eight
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) fiften
  real ( kind = 8 ) five
  real ( kind = 8 ) four8
  real ( kind = 8 ) half
  real ( kind = 8 ) hund
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) oneneg
  real ( kind = 8 ) one28
  real ( kind = 8 ) pi
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  logical sflag
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  logical tflag
  real ( kind = 8 ) thirty
  real ( kind = 8 ) three
  real ( kind = 8 ) top
  real ( kind = 8 ) twenty
  real ( kind = 8 ) two
  real ( kind = 8 ) u(0:559)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xden
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xlarge
  real ( kind = 8 ) xleast
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnine
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz
!
!  Mathematical constants
!
  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two / 2.0d0 /
  data eight /8.0d0/
  data xnine / 9.0d0 /
  data ateten /18.0d0 /
  data twenty /20.0d0 /
  data one28 /128.0d0/
  data five /5.0d0 /
  data oneneg /-1.0d0 /
  data xden / 16.0d0 /
  data all9 /-999.0d0/
  data three /3.0d0 /
  data fiften / 15.0d0 /
  data thirty /30.0d0 /
  data four8 /48.0d0/
  data pi /3.141592653589793d0/
  data hund /100.0d0/
!
!  Machine-dependent constant
!
  data xleast /2.23d-308/

  top(x) = -x - half*log(two*x) + log(one+(three/eight-fiften/ &
     one28/x)/x)

  bot(x) = - one - half/x + ((- four8*x + thirty) / &
     (((one28*x+four8)*x-fiften)*x))
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  amaxexp = real ( maxexp, kind = 8 )
  jt = 0
  b = eps
  xlam = ( xden - one ) / xden
  const = half * log (  pi ) - log ( xmin )
!
!  Random argument accuracy tests.
!
  do j = 1, 3

    sflag = ( j == 1 .and. amaxexp / ait <= five )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    n = 2000
    a = b

    if ( j == 1 ) then
      b = one
    else if ( j == 2 ) then
      b = eight
    else
      b = twenty
    end if

    xn = real ( n, kind = 8 )
    del = ( b - a ) / xn
    xl = a
!
!  Accuracy test is based on the multiplication theorem.
!
    do i = 1, n

      x = del * ren ( jt ) + xl
      y = x / xlam
      w = xden * y
      y = ( w + y ) - w
      x = y * xlam
      u(0) = besk0 ( y )
      u(1) = besk1 ( y )
      tflag = sflag .and. ( y < half )

      if ( tflag ) then
        u(0) = u(0) * eps
        u(1) = u(1) * eps
      end if

      mb = 1
      xmb = one
      y = y * half
      w = ( one - xlam ) * ( one + xlam )
      c = w * y
      t = u(0) + c * u(1)
      t1 = eps / hund

      do ii = 2, 60

        z = u(ii-1)

        if ( z / t1 < t ) then
          go to 120
        else if ( u(ii-1) > one ) then
          if ( ( xmb / y ) > ( xmax / u(ii-1) ) ) then
            xl = xl + del
            a = xl
            go to  200
          end if
        end if

        u(ii) = xmb / y * u(ii-1) + u(ii-2)

        if ( t1 > one / eps ) then
          t = t * t1
          t1 = one
        end if

        t1 = xmb * t1 / c
        xmb = xmb + one
        mb = mb + 1

      end do

  120     continue

      sum = u(mb)
      ind = mb
      mb = mb - 1

      do ii = 1, mb
        xmb = xmb - one
        ind = ind - 1
        sum = sum * w * y / xmb + u(ind)
      end do

      zz = sum

      if ( tflag ) then
        zz = zz / eps
      end if

      zz = zz * xlam
      z = besk1 ( x )
      y = z

      if ( u(0) > y ) then
        y = u(0)
      end if

      w = ( z - zz ) / y

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

  200     continue

    end do

    n = k1 + k2 + k3
    xn = real ( n, kind = 8 )
    r7 = sqrt ( r7 / xn )
    write ( *,1000)
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta
    w = all9

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    end if

    if ( j == 3 ) then
      write ( *,1024) r6,ibeta,w,x1
    else
      write ( *,1021) r6,ibeta,w,x1
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
    w = all9

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    end if

    if ( j == 3 ) then
      write ( *,1025) r7,ibeta,w
    else
      write ( *,1023) r7,ibeta,w
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

  end do
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1031)
  y = besk1 ( xleast )
  write ( *,1032) y
  y = besk1 ( xmin )
  write ( *,1036) y
  y = besk1 ( zero )
  write ( *,1033) 0,y
  x = ren ( jt ) * oneneg
  y = besk1 ( x )
  write ( *,1034) x,y
  y = besek1 ( xmax )
  write ( *,1035) y
  xa = log ( xmax )

  330 continue

  xb = xa - ( top ( xa ) + const ) / bot ( xa )

  if ( eps < abs ( xb - xa ) / xb ) then
    xa = xb
    go to 330
  end if

  xlarge = xb * xlam
  y = besk1 ( xlarge )
  write ( *,1034) xlarge,y
  xlarge = xb * ( xnine / eight )
  y = besk1 ( xlarge )
  write ( *,1034) xlarge,y

  return
 1000 format('1Test of K1(X) vs Multiplication Theorem'//)
 1010 format(I7,' random arguments were tested from the interval (', &
      F5.1,',',F5.1,')'//)
 1011 format(' ABS(K1(X)) was larger',I6,' times,'/ &
      20X,' agreed',I6,' times, and'/ &
      16X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number.'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1024 format(' The maximum absolute error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6)
 1025 format(' The root mean square absolute error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1031 format(//' Test with extreme arguments'/)
 1032 format(' K1(XLEAST) = ',E24.17/)
 1033 format(' K1(',I1,') = ',E24.17/)
 1034 format(' K1(',E24.17,' ) = ',E24.17/)
 1035 format(' E**X * K1(XMAX) = ',E24.17/)
 1036 format(' K1(XMIN) = ',E24.17/)
end
subroutine r8_psi_test ( )

!*****************************************************************************80
!
!! R8_PSI_TEST tests R8_PSI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Performance evaluation of programs related to the real gamma function,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 46-54.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) all9
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) del
  real ( kind = 8 ) eight
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) half
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) one7
  real ( kind = 8 ) one6
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) three
  real ( kind = 8 ) twenty
  real ( kind = 8 ) y
  real ( kind = 8 ) v0
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xh
  real ( kind = 8 ) xl
  real ( kind = 8 ) xl2
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) xx
  real ( kind = 8 ) x0
  real ( kind = 8 ) x01
  real ( kind = 8 ) x1
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zh
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data three / 3.0D+00/
  data eight / 8.0d0 /
  data twenty / 20.0d0 /
  data all9 /-999.0d0/
  data xl2 /6.9314718055994530942d-1/
  data one7 / -17.625D+00 /
  data one6 / -16.875d0/
  data x0 /374.0d0/
  data x01 /256.0d0/
  data v0 /-6.7240239024288040437d-04/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  ait = real ( it, kind = 8 )
  jt = 0
!
!  Random argument accuracy tests.
!
  do j = 1, 4

    k1 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    n = 2000

    if ( j == 1 ) then
      a = zero
      b = one
    else if ( j == 2 ) then
      a = b + b
      b = eight
    else if ( j == 3 ) then
      a = b
      b = twenty
    else
      a = one7
      b = one6
      n = 500
    end if

    xn = real ( n, kind = 8 )
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Carefully purify arguments and evaluate identity.
!
      xx = x * half
      xh = xx + half
      xx = xh - half
      x = xx + xx
      z = r8_psi ( x )
      zh = r8_psi ( xh )
      zz = r8_psi ( xx )
      zz = ( zz + zh ) * half + xl2
!
!  Accumulate results.
!
      w = ( zz - z ) / zz

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Process and output statistics.
!
    k2 = n - k3 - k1
    r7 = sqrt ( r7 / xn )
    if ( 2 * ( j / 2 ) /= j ) then
      write ( *,1000)
    end if
    write ( *,1001)
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    else
      w = all9
    end if

    write ( *,1021) r6,ibeta,w,x1
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    else
      w = all9
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

  end do
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  x = x0 / x01
  y = r8_psi ( x )
  z = ( y - v0 ) / v0

  if ( z /= zero ) then
    w = log ( abs ( z ) ) / albeta
  else
    w = all9
  end if

  w = max ( ait + w, zero )
  write ( *,1031) x,y,ibeta,w
  write ( *,1033)

  if ( xmax * xmin >= one ) then
    x = xmin
  else
    x = one / xmax
  end if

  write ( *,1035) x
  y = r8_psi ( x )
  write ( *,1036) y
  x = xmax
  write ( *,1035) x
  y = r8_psi ( x )
  write ( *,1036) y
!
!  Test of error returns
!
  write ( *,1037)
  x = zero
  write ( *,1034) x
  y = r8_psi ( x )
  write ( *,1036) y
  x = -three / eps
  write ( *,1034) x
  y = r8_psi ( x )
  write ( *,1036) y
  write ( *, '(a)' ) '  This concludes the tests.'
  return

 1000 format('1')
 1001 format(' Test of PSI(X) vs (PSI(X/2)+PSI(X/2+1/2))/2 + ln(2)' &
   //)
 1010 format(I7,' random arguments were tested from the interval ', &
   1H(,F5.1,1H,,F5.1,1H)//)
 1011 format(' ABS(PSI(X)) was larger',I6,' times', / &
       21X,' agreed',I6,' times, and'/ &
     17X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number.' //)
 1021 format(' The maximum relative error of',E15.4,3H = ,I4,3H **, &
    F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
    ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      3H = ,I4,3H **,F7.2)
 1031 format(' Accuracy near positive zero'//' PSI(',E14.7,') = ', &
      E24.17/13X,'Loss of base',I3,' digits = ',F7.2/)
 1033 format(//' Test with extreme arguments'/)
 1034 format(' PSI will be called with the argument ',E17.10/ &
       ' This may stop execution.'/)
 1035 format(' PSI will be called with the argument ',E17.10/ &
       ' This should not stop execution.'/)
 1036 format(' PSI returned the value',E25.17//)
 1037 format(//' Test of error returns'//)
end
subroutine ri_test ( )

!*****************************************************************************80
!
!! RI_TEST tests RIBESL.
!
!  Discussion:
!
!    Two different accuracy tests are used.  In the first interval,
!    function values are compared against values generated with the
!    multiplication formula, where the Bessel values used in the
!    multiplication formula are obtained from the function program.
!    In the remaining intervals, function values are compared
!    against values generated with a local Taylor series expansion.
!    Derivatives in the expansion are expressed in terms of the
!    first two Bessel functions, which are in turn obtained from
!    the function program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) ak
  real ( kind = 8 ) akk
  real ( kind = 8 ) albeta
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alphsq
  real ( kind = 8 ) a1
  real ( kind = 8 ) ar1(11,6)
  real ( kind = 8 ) ar2(13,9)
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) delta
  real ( kind = 8 ) deriv
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) f
  real ( kind = 8 ) g(5)
  real ( kind = 8 ) half
  real ( kind = 8 ) hund
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) ize
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) last
  integer ( kind = 4 ) m
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) mborg
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) mvr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncalc
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) ndx(24)
  integer ( kind = 4 ) ndx2(8)
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  integer ( kind = 4 ) nk
  integer ( kind = 4 ) no1
  integer ( kind = 4 ) num
  real ( kind = 8 ) one
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) six
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) ten
  real ( kind = 8 ) two
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u(560)
  real ( kind = 8 ) u2(560)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xbad
  real ( kind = 8 ) xj1
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xlarge
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) x99
  real ( kind = 8 ) y
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two / 2.0d0 /
  data six /6.0d0/
  data ten / 10.0d0 /
  data sixten /1.6d1 /
  data hund /1.0d2 /
  data x99 /-999.0d0/
  data xlam /1.03125d0/
  data xlarge /1.0d4/
  data c /0.9189385332d0/
!
!  Arrays related to expansion of the derivatives in terms
!  of the first two Bessel functions.
!
  data  ndx/9,7,5,3,1,8,6,4,2,7,5,3,1,6,4,2,5,3,1,4,2,3,1,2/

  data  ndx2/5,9,13,16,19,21,23,24/

  data  ar1/0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,1.0d0,3.0d0,0.0d0,-2.0d0, &
            -1.2d1,0.0d0,1.0d0,0.0d0,-1.0d0,1.0d0,2.0d0,0.0d0, &
            -2.0d0,-6.0d0,1.0d0,7.0d0,2.4d1,0.0d0,0.0d0,1.0d0,0.0d0, &
            -3.0d0,0.0d0,2.0d0,1.1d1,0.0d0,-1.2d1,-5.0d1,0.0d0, &
            0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,-6.0d0,0.0d0,2.0d0, &
            3.5d1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0, &
            0.0d0,0.0d0,-1.0d1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            0.0d0,0.0d0,0.0d0,0.0d0,1.0d0/

  data  ar2/1.0d0,9.0d0,6.0d1,0.0d0,-3.0d0,-5.1d1,-3.6d2,0.0d0, &
            1.0d0,1.8d1,3.45d2,2.52d3,0.0d0,0.0d0,-3.0d0,-3.3d1, &
            -1.2d2,1.0d0,1.5d1,1.92d2,7.2d2,0.0d0,-4.0d0,-9.6d1, &
            -1.32d3,-5.04d3,0.0d0,3.0d0,7.8d1,2.74d2,0.0d0,-2.7d1, &
            -5.7d2,-1.764d3,0.0d0,4.0d0,2.46d2,4.666d3,1.3068d4, &
            0.0d0,0.0d0,-1.8d1,-2.25d2,0.0d0,3.0d0,1.5d2,1.624d3, &
            0.0d0,0.0d0,-3.6d1,-1.32d3,-1.3132d4,0.0d0,0.0d0,3.0d0, &
            8.5d1,0.0d0,0.0d0,-4.5d1,-7.35d2,0.0d0,0.0d0,6.0d0, &
            5.5d2,6.769d3,0.0d0,0.0d0,0.0d0,-1.5d1,0.0d0,0.0d0, &
            3.0d0,1.75d2,0.0d0,0.0d0,0.0d0,-6.0d1,-1.96d3,0.0d0, &
            0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,-2.1d1,0.0d0,0.0d0, &
            0.0d0,4.0d0,3.22d2,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-2.8d1,0.0d0,0.0d0, &
            0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            0.0d0,1.0d0/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  ait = real ( it, kind = 8 )
  albeta = log ( beta )
  a = zero
  b = two
  jt = 0
  delta = xlam - one
  f = delta * ( xlam + one ) * half
!
!  Random argument accuracy tests.
!
  do j = 1, 4
!
!  Determine the number of terms needed for convergence of the series
!  used in the multiplication theorem.  Use Newton iteration on the
!  asymptotic form of the convergence check for I0(X).
!
    xbad = b
    d = ait * albeta - c + one
    e = log ( xbad * f ) + one
    akk = one

  100   continue

    ak = akk
    z = d + e * ak - ( ak + half ) * log ( ak + one )
    zz = e - ( ak + half ) / ( ak + one ) - log ( ak + one )
    akk = ak - z/zz

    if ( abs ( ak - akk ) > hund * eps * ak ) then
      go to 100
    end if

    mborg = int ( akk ) + 1
    n = 2000
    xn = real ( n, kind = 8 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    a1 = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      mb = mborg
      x = del * ren ( jt ) + xl
      alpha = ren ( jt )
      ize = 1
!
!  Carefully purify arguments.
!
      if ( j == 1 ) then
        y = x / xlam
      else
        y = x - delta
      end if

      w = sixten * y
      t1 = w + y
      t1 = w + t1
      y = t1 - w
      y = y - w

      if ( j == 1 ) then
        x = y * xlam
      else
        x = y + delta
      end if

      call ribesl ( y, alpha, mb, ize, u2, ncalc )

      if ( j == 1 ) then
!
!  Accuracy test is based on the multiplication theorem.
!
        d = f * y
        mb = ncalc - 2
        xmb = real ( mb, kind = 8 )
        sum = u2(mb+1)
        ind = mb

        do ii = 2, mb
          sum = sum * d / xmb + u2(ind)
          ind = ind - 1
          xmb = xmb - one
        end do

        zz = sum * d + u2(ind)
        zz = zz * xlam**alpha

      else
!
!  Accuracy test is based on local Taylor's series expansion.
!
        ysq = y * y
        alphsq = alpha * alpha
        mb = 8
        j1 = mb
        xj1 = real ( j1 + 1, kind = 8 )
        iexp = 0
        nk = 13
        num = 2

        do ii = 1, mb

          if ( nk == 0 ) then
            nk = 11
            num = 1
          end if

          k = 9 - j1

          if ( k > 1 ) then
            no1 = ndx2(k-1) + 1
          else
            no1 = 1
          end if

          mvr = no1
          last = ndx2(k)
          k = last - no1 + 1
!
!  Group I(ALPHA) terms in the derivative.
!
          do iii = 1, k

            j2 = ndx(mvr)

            if ( num == 1 ) then
              g(iii) = ar1(nk,j2)
            else
              g(iii) = ar2(nk,j2)
            end if

            if ( j2 > 1 ) then

  157             continue

              j2 = j2 - 1

              if ( num == 1 ) then
                g(iii) = g(iii) * alpha + ar1(nk,j2)
              else
                g(iii) = g(iii) * alpha + ar2(nk,j2)
              end if

              if ( j2 > 1 ) then
                go to 157
              end if

            end if

            mvr = mvr + 1
            nk = nk - 1

          end do

          t1 = g(1)
          do iii = 2, k
            t1 = t1 / ysq + g(iii)
          end do

          if ( iexp == 1 ) then
            t1 = t1 / y
          end if
!
!  Group I(ALPHA+1) terms in the derivative.
!
          iexp = 1 - iexp
          nk = nk + k
          mvr = no1
          kk = k

          do iii = 1, k

            j2 = ndx(mvr)
            m = mod ( j2, 2 )

            if ( m == 1 ) then
              j2 = j2 - 1
            end if

            if ( j2 >= 2 ) then

              if ( num == 1 ) then
                g(iii) = ar1(nk,j2)
              else
                g(iii) = ar2(nk,j2)
              end if

  163             continue

              j2 = j2 - 2

              if ( j2 >= 2 ) then

                if ( num == 1 ) then
                  g(iii) = g(iii) * alphsq + ar1(nk,j2)
                else
                  g(iii) = g(iii) * alphsq + ar2(nk,j2)
                end if

                go to 163
              end if

            else

              kk = iii - 1

            end if

            mvr = mvr + 1
            nk = nk - 1

          end do

          t2 = g(1)
          do iii = 2, kk
            t2 = t2 / ysq + g(iii)
          end do

          if ( iexp == 1 ) then
            t2 = t2 / y
          end if

          deriv = u2(1) * t1 + u2(2) * t2

          if ( j1 == 8 ) then
            sum = deriv
          else
            sum = sum * delta / xj1 + deriv
          end if

          j1 = j1 - 1
          xj1 = xj1 - one

        end do

        zz = sum * delta + u2(1)

      end if

      mb = 2
      call ribesl ( x, alpha, mb, ize, u, ncalc )
      z = u(1)
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
       r6 = w
       x1 = x
       a1 = alpha
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    k2 = n - k1 - k3
    r7 = sqrt ( r7 / xn )

    if ( j == 1 ) then
      write ( *,1000)
    else
      write ( *,1001)
    end if

    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    else
      w = x99
    end if

    write ( *,1021) r6,ibeta,w,x1,a1
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    else
      w = x99
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b
    b = b + b

    if ( j == 2 ) then
      b = ten
    end if

  end do
!
!  Test of error returns.
!
!  First, test with bad parameters.
!
  write ( *, 2006)
  x = one
  alpha = one + half
  mb = 5
  ize = 2
  u(1) = zero
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  alpha = half
  mb = -mb
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  mb = -mb
  ize = 5
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
!
!  Last tests are with extreme parameters.
!
  x = zero
  alpha = ren ( jt )
  mb = 2
  ize = 1
  u(1) = zero
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  alpha = zero
  mb = 2
  u(1) = zero
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  alpha = one
  mb = 2
  u(1) = zero
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  x = -one
  alpha = half
  mb = 5
  ize = 2
  u(1) = zero
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2012) x
  write ( *, 2013)
  write ( *, 2014) ncalc,u(1)
!
!  Determine largest safe argument for scaled functions.
!
  write ( *, 2015)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This RIBESL test will be skipped.'
  write ( *, '(a)' ) 'It causes a floating exception.'

  if ( .false. ) then
    x = xlarge * ( one - sqrt ( sqrt ( eps ) ) )
    ize = 2
    mb = 2
    u(1) = zero
    call ribesl ( x, alpha, mb, ize, u, ncalc )
  end if

  write ( *, 2012) x
  write ( *, 2014) ncalc,u(1)
  x = xlarge * ( one + sqrt ( sqrt ( eps ) ) )
  mb = 2
  u(1) = zero
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2012) x
  write ( *, 2013)
  write ( *, 2014) ncalc,u(1)
!
!  Determine largest safe argument for unscaled functions.
!
  write ( *, 2016)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This RIBESL test will be skipped.'
  write ( *, '(a)' ) 'It causes a floating exception.'

  if ( .false. ) then
    n = int ( log ( xmax ) )
    z = real ( n, kind = 8 )
    x = z * ( one - sqrt ( sqrt ( eps ) ) )
    ize = 1
    mb = 2
    u(1) = zero
    call ribesl ( x, alpha, mb, ize, u, ncalc )
  end if

  write ( *, 2012) x
  write ( *, 2014) ncalc,u(1)
  x = z * ( one + sqrt ( sqrt ( eps ) ) )
  mb = 2
  u(1) = zero
  call ribesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2012) x
  write ( *, 2013)
  write ( *, 2014) ncalc,u(1)
  write ( *, 2020)
  return
 1000 format('1Test of I(X,ALPHA) vs Multiplication Theorem'//)
 1001 format('1Test of I(X,ALPHA) vs Taylor series'//)
 1010 format(I7,' Random arguments were tested from the interval ', &
      '(',F5.2,',',F5.2,')'//)
 1011 format(' I(X,ALPHA) was larger',I6,' times,'/ &
      15X,' agreed',I6,' times, and'/ &
      11X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6,' and NU =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 2006 format('1Check of Error Returns'/// &
      ' The following summarizes calls with indicated parameters'// &
      ' NCALC different from MB indicates some form of error'// &
      ' See documentation for RIBESL for details'// &
      7X,'ARG',12X,'ALPHA',6X,'MB',3X,'IZ',7X,'RES',6X,'NCALC'//)
 2011 format(2E15.7,2I5,E15.7,I5//)
 2012 format(' RIBESL will be called with the argument',E13.6)
 2013 format(' This should trigger an error message.')
 2014 format(' NCALC returned the value',I5/ &
      ' and RIBESL returned the value',E13.6/)
 2015 format(' Tests near the largest argument for scaled functions'/)
 2016 format(' Tests near the largest argument for unscaled functions'/)
 2020 format(' This concludes the tests.')
end
subroutine rj_test ( )

!*****************************************************************************80
!
!! RJ_TEST tests RJBESL.
!
!  Discussion:
!
!    Two different accuracy tests are used.  In the first interval,
!    function values are compared against values generated with the
!    multiplication formula, where the Bessel values used in the
!    multiplication formula are obtained from the function program.
!    In the remaining intervals, function values are compared
!    against values generated with a local Taylor series expansion.
!    Derivatives in the expansion are expressed in terms of the
!    first two Bessel functions, which are in turn obtained from
!    the function program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
!  Acknowledgement:
!
!    This program is a minor modification of the test driver for RIBESL
!    whose primary author was Laura Stoltz.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alphsq
  real ( kind = 8 ) a1
  real ( kind = 8 ) ar1(11,6)
  real ( kind = 8 ) ar2(13,9)
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) delta
  real ( kind = 8 ) deriv
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) f
  real ( kind = 8 ) fxmx
  real ( kind = 8 ) g(5)
  real ( kind = 8 ) half
  real ( kind = 8 ) hund
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) last
  integer ( kind = 4 ) m
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) mvr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncalc
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) ndx(24)
  integer ( kind = 4 ) ndx2(8)
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  integer ( kind = 4 ) nk
  integer ( kind = 4 ) no1
  integer ( kind = 4 ) num
  real ( kind = 8 ) one
  real ( kind = 8 ) ren
  real ( kind = 8 ) rtpi
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) six
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) ten
  real ( kind = 8 ) two
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u(560)
  real ( kind = 8 ) u2(560)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xj1
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xlarge
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) x99
  real ( kind = 8 ) y
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data two / 2.0d0/
  data six /6.0d0/
  data ten /10.0d0 /
  data sixten /1.6d1 /
  data hund /1.0d2/
  data x99 /-999.0d0/
  data xlam / 1.03125d0/
  data xlarge / 1.0d4/
  data c /0.9189385332d0 /
  data rtpi /0.6366d0/
!
!  Arrays related to expansion of the derivatives in terms
!  of the first two Bessel functions.
!
  data  ndx/9,7,5,3,1,8,6,4,2,7,5,3,1,6,4,2,5,3,1,4,2,3,1,2/

  data  ndx2/5,9,13,16,19,21,23,24/

  data ar1/0.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,1.0d0,-3.0d0,0.0d0,-2.0d0, &
           1.2d1,0.0d0,1.0d0,0.0d0,-1.0d0,-1.0d0,2.0d0,0.0d0, &
           2.0d0,-6.0d0,1.0d0,-7.0d0,2.4d1,0.0d0,0.0d0,1.0d0,0.0d0, &
           -3.0d0,0.0d0,-2.0d0,1.1d1,0.0d0,1.2d1,-5.0d1,0.0d0, &
           0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,-6.0d0,0.0d0,-2.0d0, &
           3.5d1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0, &
           0.0d0,0.0d0,-1.0d1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
           0.0d0,0.0d0,0.0d0,0.0d0,1.0d0/

  data ar2/-1.0d0,9.0d0,-6.0d1,0.0d0,3.0d0,-5.1d1,3.6d2,0.0d0, &
            1.0d0,-1.8d1,3.45d2,-2.52d3,0.0d0,0.0d0,-3.0d0,3.3d1, &
            -1.2d2,-1.0d0,1.5d1,-1.92d2,7.2d2,0.0d0,4.0d0,-9.6d1, &
            1.32d3,-5.04d3,0.0d0,3.0d0,-7.8d1,2.74d2,0.0d0,-2.7d1, &
            5.7d2,-1.764d3,0.0d0,-4.0d0,2.46d2,-4.666d3,1.3068d4, &
            0.0d0,0.0d0,1.8d1,-2.25d2,0.0d0,3.0d0,-1.5d2,1.624d3, &
            0.0d0,0.0d0,-3.6d1,1.32d3,-1.3132d4,0.0d0,0.0d0,-3.0d0, &
            8.5d1,0.0d0,0.0d0,4.5d1,-7.35d2,0.0d0,0.0d0,6.0d0, &
            -5.5d2,6.769d3,0.0d0,0.0d0,0.0d0,-1.5d1,0.0d0,0.0d0, &
            -3.0d0,1.75d2,0.0d0,0.0d0,0.0d0,6.0d1,-1.96d3,0.0d0, &
            0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,-2.1d1,0.0d0,0.0d0, &
            0.0d0,-4.0d0,3.22d2,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-2.8d1,0.0d0,0.0d0, &
            0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            0.0d0,1.0d0/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  ait = real ( it, kind = 8 )
  albeta = log ( beta )
  a = zero
  b = two
  jt = 0
  delta = xlam - one
  f = delta * ( xlam + one ) * half
!
!  Random argument accuracy tests.
!
  do j = 1, 4

    n = 2000
    xn = real ( n, kind = 8 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    a1 = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Carefully purify arguments.
!
  110     continue

      alpha = ren ( jt )

      if ( j == 1 ) then
        y = x / xlam
      else
        y = x - delta
      end if

      w = sixten * y
      t1 = w + y
      t1 = w + t1
      y = t1 - w
      y = y - w

      if ( j == 1 ) then
        x = y * xlam
        mb = 11
      else
        x = y + delta
        mb = 2
      end if

      call rjbesl ( y, alpha, mb, u2, ncalc )
      call rjbesl ( x, alpha, mb, u, ncalc )
!
!  Accuracy test is based on the multiplication theorem.
!
      if ( j == 1 ) then

        d = -f * y
        mb = ncalc - 2
        xmb = real ( mb, kind = 8 )
        sum = u2(mb+1)
        ind = mb

        do ii = 2, mb
          sum = sum * d / xmb + u2(ind)
          ind = ind - 1
          xmb = xmb - one
        end do

        zz = sum * d + u2(ind)
        zz = zz * xlam**alpha
!
!  Accuracy test is based on local Taylor's series expansion.
!
      else

        if ( abs ( u(1) ) < abs ( u2(1) ) ) then

          z = x
          x = y
          y = z
          delta = x - y

          do ii = 1, mb
            z = u(ii)
            u(ii) = u2(ii)
            u2(ii) = z
          end do

        end if
!
!  Filter out cases where function values or derivatives are small.
!
        w = sqrt ( rtpi / x ) / sixten
        z = min ( abs ( u2(1) ), abs ( u2(2) ) )

        if ( z < w) then
          go to 110
        end if

        ysq = y * y
        alphsq = alpha * alpha
        mb = 8
        j1 = mb
        xj1 = real ( j1 + 1, kind = 8 )
        iexp = 0
        nk = 13
        num = 2

        do ii = 1, mb

          if ( nk == 0 ) then
            nk = 11
            num = 1
          end if

          k = 9 - j1

          if ( k > 1 ) then
            no1 = ndx2(k-1) + 1
          else
            no1 = 1
          end if

          mvr = no1
          last = ndx2(k)
          k = last - no1 + 1
!
!  Group J(ALPHA) terms in the derivative.
!
          do iii = 1, k

            j2 = ndx(mvr)

            if ( num == 1 ) then
              g(iii) = ar1(nk,j2)
            else
              g(iii) = ar2(nk,j2)
            end if

            if ( j2 > 1 ) then

  157             continue

              j2 = j2 - 1

              if ( num == 1 ) then
                g(iii) = g(iii) * alpha + ar1(nk,j2)
              else
                g(iii) = g(iii) * alpha + ar2(nk,j2)
              end if

              if ( j2 > 1) then
                go to 157
              end if

            end if

            mvr = mvr + 1
            nk = nk - 1

          end do

          t1 = g(1)

          do iii = 2, k
            t1 = t1 / ysq + g(iii)
          end do

          if ( iexp == 1 ) then
            t1 = t1 / y
          end if
!
!  Group J(ALPHA+1) terms in the derivative.
!
          iexp = 1 - iexp
          nk = nk + k
          mvr = no1
          kk = k

          do iii = 1, k

            j2 = ndx(mvr)
            m = mod ( j2, 2 )

            if ( m == 1 ) then
              j2 = j2 - 1
            end if

            if ( j2 >= 2 ) then

              if ( num == 1 ) then
                g(iii) = ar1(nk,j2)
              else
                g(iii) = ar2(nk,j2)
              end if

  163             continue

              j2 = j2 - 2

              if ( j2 >= 2 ) then

                if ( num == 1 ) then
                  g(iii) = g(iii) * alphsq + ar1(nk,j2)
                else
                  g(iii) = g(iii) * alphsq + ar2(nk,j2)
                end if

                go to 163
              end if

            else

              kk = iii - 1

            end if

            mvr = mvr + 1
            nk = nk - 1

          end do

          t2 = g(1)
          do iii = 2, kk
            t2 = t2 / ysq + g(iii)
          end do

          if ( iexp == 1 ) then
            t2 = t2 / y
          end if

          deriv = u2(1) * t1 - u2(2) * t2

          if ( j1 == 8 ) then
            sum = deriv
          else
            sum = sum * delta / xj1 + deriv
          end if

          j1 = j1 - 1
          xj1 = xj1 - one

        end do

        zz = sum * delta + u2(1)

      end if

      mb = 2
      z = u(1)
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
        a1 = alpha
        fxmx = z
      end if

      r7 = r7 + w * w
      xl = xl + del

    end do
!
!  Gather and print statistics for test.
!
    k2 = n - k1 - k3
    r7 = sqrt ( r7 / xn )

    if ( j == 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) &
        '  Test of J(X,ALPHA) vs Multiplication Theorem'
      write ( *, '(a)' ) ' '

    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) &
        '  Test of J(X,ALPHA) vs Taylor series'
      write ( *, '(a)' ) ' '
    end if

    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    else
      w = x99
    end if

    write ( *,1021) r6,ibeta,w,x1,a1
    write ( *,1024) fxmx
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    else
      w = x99
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b

    if ( j == 1 ) then
      b = ten
    else if ( j == 2 ) then
      b = b + b
    else if ( j == 3 ) then
      a = a + ten
      b = a + ten
    end if

  end do
!
!  Test of error returns.
!
!  First, test with bad parameters.
!
  write ( *, 2006)
  x = one
  alpha = one + half
  mb = 5
  call rjbesl ( x, alpha, mb, u, ncalc )
  write ( *, 2011) x,alpha,mb,u(1),ncalc
  alpha = half
  mb = -mb
  call rjbesl ( x, alpha, mb, u, ncalc )
!
!  Last tests are with extreme parameters.
!
  x = zero
  alpha = one
  mb = 2
  call rjbesl ( x, alpha, mb, u, ncalc )
  write ( *, 2011) x,alpha,mb,u(1),ncalc
  x = -one
  alpha = half
  mb = 5
  call rjbesl ( x, alpha, mb, u, ncalc )
  write ( *, 2011) x,alpha,mb,u(1),ncalc
!
!  Determine largest safe argument.
!
  write ( *, 2015)
  x = xlarge * ( one - sqrt ( sqrt ( eps ) ) )
  mb = 2
  call rjbesl ( x, alpha, mb, u, ncalc )
  write ( *, 2012) x
  write ( *, 2014) ncalc,u(1)
  x = xlarge * ( one + sqrt ( sqrt ( eps ) ) )
  mb = 2
  call rjbesl ( x, alpha, mb, u, ncalc )
  write ( *, 2012) x
  write ( *, 2013)
  write ( *, 2014) ncalc,u(1)
  write ( *, 2020)
  return

 1010 format(I7,' Random arguments were tested from the interval ', &
      '(',F5.2,',',F5.2,')'//)
 1011 format(' J(X,ALPHA) was larger',I6,' times,'/ &
      15X,' agreed',I6,' times, and'/ &
      11X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6,' and NU =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1024 format(4x,'with J(X,ALPHA) = ',E13.6)
 2006 format('1Check of Error Returns'/// &
      ' The following summarizes calls with indicated parameters'// &
      ' NCALC different from MB indicates some form of error'// &
      ' See documentation for RJBESL for details'// &
      7X,'ARG',12X,'ALPHA',6X,'MB',6X,'B(1)',6X,'NCALC'//)
 2011 format(2E15.7,I5,E15.7,I5//)
 2012 format(' RJBESL will be called with the argument',E13.6)
 2013 format(' This should trigger an error message.')
 2014 format(' NCALC returned the value',I5/ &
      ' and RJBESL returned U(1) = ',E13.6/)
 2015 format(' Tests near the largest acceptable argument for RJBESL'/)
 2020 format(' This concludes the tests.')
end
subroutine rk_test ( )

!*****************************************************************************80
!
!! RK_TEST tests RKBESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) alpha
  real ( kind = 8 ) amaxexp
  real ( kind = 8 ) a1
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) eight
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) five
  real ( kind = 8 ) half
  real ( kind = 8 ) hund
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) ize
  integer ( kind = 4 ) iz1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncalc
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  logical sflag
  logical tflag
  real ( kind = 8 ) one
  real ( kind = 8 ) ovrchk
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) ten
  real ( kind = 8 ) tinyx
  real ( kind = 8 ) t1
  real ( kind = 8 ) u(560)
  real ( kind = 8 ) u2(560)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnine
  real ( kind = 8 ) x1
  real ( kind = 8 ) x99
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data five / 5.0d0 /
  data eight / 8.0d0/
  data xnine /9.0d0 /
  data ten /10.0d0 /
  data hund / 1.0d2 /
  data x99 /-999.0d0 /
  data sixten /1.6d1 /
  data xlam /0.9375d0/
  data c / 2.2579d-1/
  data tinyx /1.0d-10/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  ait = real ( it, kind = 8 )
  albeta = log ( beta )
  amaxexp = real ( maxexp, kind = 8 )
  a = eps
  b = one
  jt = 0
!
!  Random argument accuracy tests.
!
  do j = 1, 3

    sflag = ( j == 1 .and. amaxexp / ait <= five )
    n = 2000
    xn = real ( n, kind = 8 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    a1 = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
      alpha = ren ( jt )
!
!  Accuracy test is based on the multiplication theorem.
!
      ize = 1
      mb = 3
!
!  Carefully purify arguments.
!
      y = x / xlam
      w = sixten * y
      t1 = w + y
      y = t1 - w
      x = y * xlam
      call rkbesl ( y, alpha, mb, ize, u2, ncalc )
      tflag = sflag .and. ( y < half )

      if ( tflag ) then
        u2(1) = u2(1) * eps
        u2(2) = u2(2) * eps
      end if

      mb = 1
      xmb = zero
      w = ( one - xlam ) * ( one + xlam ) * half
      d = w * y
      t = u2(1) + d * u2(2)
      t1 = eps / hund
!
!  Generate terms using recurrence.
!
      do ii = 3, 35

        xmb = xmb + one
        t1 = xmb * t1 / d
        z = u2(ii-1)
        ovrchk = ( xmb + alpha ) / ( y * half )

        if ( z / t1 < t ) then
          go to 120
        else if ( u2(ii-1) > one ) then
          if ( ovrchk > ( xmax / u2(ii-1) ) ) then
            xl = xl + del
            a = xl
            go to 200
          end if
        end if

        u2(ii) = ovrchk * u2(ii-1) + u2(ii-2)

        if ( t1 > one / eps ) then
          t = t * t1
          t1 = one
        end if

        mb = mb + 1

      end do
!
!  Accurate summation.
!
      xmb = xmb + one

  120     continue

      sum = u2(mb+1)
      ind = mb

      do ii = 1, mb
        sum = sum * d / xmb + u2(ind)
        ind = ind - 1
        xmb = xmb - one
      end do

      zz = sum
      if ( tflag ) then
        zz = zz / eps
      end if
      zz = zz * xlam ** alpha
      mb = 2
      call rkbesl ( x, alpha, mb, ize, u, ncalc )
      z = u(1)

      y = z
      if ( u2(1) > y ) then
        y = u2(1)
      end if
!
!  Accumulate results.
!
      w = ( z - zz ) / y

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
        a1 = alpha
        iz1 = ize
      end if

      r7 = r7 + w * w
      xl = xl + del
  200     continue

    end do
!
!  Gather and print statistics for test.
!
    n = k1 + k2 + k3
    r7 = sqrt ( r7 / xn )
    write ( *,1000)
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta
    w = x99

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    end if

    if ( j == 3 ) then
      write ( *,1024) r6,ibeta,w,x1,a1,iz1
    else
      write ( *,1021) r6,ibeta,w,x1,a1,iz1
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
    w = x99

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    end if

    if ( j == 3 ) then
      write ( *,1025) r7,ibeta,w
    else
      write ( *,1023) r7,ibeta,w
    end if
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b
    b = b + b
    if ( j == 1 ) then
      b = ten
    end if

  end do
!
!  Test of error returns.
!
!  First, test with bad parameters.
!
  write ( *, 2006)
  x = -one
  alpha = half
  mb = 5
  ize = 2
  u(1) = zero
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  x = -x
  alpha = one + half
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  alpha = half
  mb = -mb
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  mb = -mb
  ize = 5
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
!
!  Last tests are with extreme parameters.
!
  x = xmin
  alpha = zero
  mb = 2
  ize = 2
  u(1) = zero
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  x = tinyx * ( one - sqrt ( eps ) )
  mb = 20
  u(1) = zero
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  x = tinyx * ( one + sqrt ( eps ) )
  mb = 20
  u(1) = zero
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
!
!  Determine largest safe argument for unscaled functions.
!
  z = log ( xmin )
  w = z - c
  zz = - z - ten

  350 continue

  z = zz
  zz = one / ( eight * z )
  a = z + log ( z ) * half + zz * ( one - xnine * half * zz ) + w
  b = one + ( half - zz * ( one - xnine * zz ) ) / z
  zz = z - a / b

  if ( abs ( z - zz ) > hund * eps * z ) then
    go to 350
  end if

  x = zz * xlam
  ize = 1
  mb = 2
  u(1) = zero
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  x = zz
  mb = 2
  u(1) = zero
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  x = ten / eps
  ize = 2
  u(1) = zero
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  x = xmax
  u(1) = zero
  call rkbesl ( x, alpha, mb, ize, u, ncalc )
  write ( *, 2011) x,alpha,mb,ize,u(1),ncalc
  return

 1000 format('1Test of K(X,ALPHA) vs Multiplication Theorem'//)
 1010 format(I7,' Random arguments were tested from the interval ', &
      '(',F5.2,',',F5.2,')'//)
 1011 format(' K(X,ALPHA) was larger',I6,' times,'/ &
      15X,' agreed',I6,' times, and'/ &
      11X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6,', NU =',E13.6, &
      ' and IZE =',I2)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1024 format(' The maximum absolute error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6,', NU =',E13.6, &
      ' and IZE =',I2)
 1025 format(' The root mean square absolute error was',E15.4, &
      ' = ',I4,' **',F7.2)
 2006 format('1Check of Error Returns'/// &
      ' The following summarizes calls with indicated parameters'// &
      ' NCALC different from MB indicates some form of error'// &
      ' See documentation for RKBESL for details'// &
      7X,'ARG',12X,'ALPHA',6X,'MB',3X,'IZ',7X,'RES',6X,'NCALC'//)
 2011 format(2E15.7,2I5,E15.7,I5//)
end
subroutine ry_test ( )

!*****************************************************************************80
!
!! RY_TEST tests RYBESL.
!
!  Discussion:
!
!    Two different accuracy tests are used.  In the first interval,
!    function values are compared against values generated with the
!    multiplication formula, where the Bessel values used in the
!    multiplication formula are obtained from the function program.
!    In the remaining intervals, function values are compared
!    against values generated with a local Taylor series expansion.
!    Derivatives in the expansion are expressed in terms of the
!    first two Bessel functions, which are in turn obtained from
!    the function program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
!    William Cody, Laura Stoltz,
!    The Use of Taylor series to test accuracy of function programs,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 55-63.
!
!  Acknowledgement:
!
!    This program is a minor modification of the test
!    driver for RIBESL whose primary author was Laura Stoltz.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alphsq
  real ( kind = 8 ) a1
  real ( kind = 8 ) ar1(11,6)
  real ( kind = 8 ) ar2(13,9)
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) delta
  real ( kind = 8 ) deriv
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) f
  real ( kind = 8 ) fxmx
  real ( kind = 8 ) g(5)
  real ( kind = 8 ) half
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) last
  integer ( kind = 4 ) m
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) mvr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncalc
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) ndx(24)
  integer ( kind = 4 ) ndx2(8)
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  integer ( kind = 4 ) nk
  integer ( kind = 4 ) no1
  integer ( kind = 4 ) num
  real ( kind = 8 ) one
  real ( kind = 8 ) onep25
  real ( kind = 8 ) p875
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) ten
  real ( kind = 8 ) two
  real ( kind = 8 ) twobpi
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u(20)
  real ( kind = 8 ) u2(20)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xj1
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) x99
  real ( kind = 8 ) y
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero   / 0.0D+00 /
  data half   / 0.5D+00 /
  data one    / 1.0D+00 /
  data onep25 / 1.25d0 /
  data two    / 2.0d0/
  data p875   / 0.875d0 /
  data ten    / 10.0d0 /
  data sixten / 1.6d1/
  data x99    /-999.0d0/
  data xlam   / 1.03125d0/
  data twobpi / 0.6366d0/
!
!  Arrays related to expansion of the derivatives in terms
!  of the first two Bessel functions.
!
  data  ndx/9,7,5,3,1,8,6,4,2,7,5,3,1,6,4,2,5,3,1,4,2,3,1,2/
  data  ndx2/5,9,13,16,19,21,23,24/
  data ar1/0.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,1.0d0,-3.0d0,0.0d0,-2.0d0, &
           1.2d1,0.0d0,1.0d0,0.0d0,-1.0d0,-1.0d0,2.0d0,0.0d0, &
           2.0d0,-6.0d0,1.0d0,-7.0d0,2.4d1,0.0d0,0.0d0,1.0d0,0.0d0, &
           -3.0d0,0.0d0,-2.0d0,1.1d1,0.0d0,1.2d1,-5.0d1,0.0d0, &
           0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,-6.0d0,0.0d0,-2.0d0, &
           3.5d1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0, &
           0.0d0,0.0d0,-1.0d1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
           0.0d0,0.0d0,0.0d0,0.0d0,1.0d0/
  data ar2/-1.0d0,9.0d0,-6.0d1,0.0d0,3.0d0,-5.1d1,3.6d2,0.0d0, &
            1.0d0,-1.8d1,3.45d2,-2.52d3,0.0d0,0.0d0,-3.0d0,3.3d1, &
            -1.2d2,-1.0d0,1.5d1,-1.92d2,7.2d2,0.0d0,4.0d0,-9.6d1, &
            1.32d3,-5.04d3,0.0d0,3.0d0,-7.8d1,2.74d2,0.0d0,-2.7d1, &
            5.7d2,-1.764d3,0.0d0,-4.0d0,2.46d2,-4.666d3,1.3068d4, &
            0.0d0,0.0d0,1.8d1,-2.25d2,0.0d0,3.0d0,-1.5d2,1.624d3, &
            0.0d0,0.0d0,-3.6d1,1.32d3,-1.3132d4,0.0d0,0.0d0,-3.0d0, &
            8.5d1,0.0d0,0.0d0,4.5d1,-7.35d2,0.0d0,0.0d0,6.0d0, &
            -5.5d2,6.769d3,0.0d0,0.0d0,0.0d0,-1.5d1,0.0d0,0.0d0, &
            -3.0d0,1.75d2,0.0d0,0.0d0,0.0d0,6.0d1,-1.96d3,0.0d0, &
            0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,-2.1d1,0.0d0,0.0d0, &
            0.0d0,-4.0d0,3.22d2,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-2.8d1,0.0d0,0.0d0, &
            0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
            0.0d0,1.0d0/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  ait = real ( it, kind = 8 )
  albeta = log ( beta )
  a = zero
  b = two
  jt = 0
  delta = xlam - one
  f = delta * ( xlam + one ) * half
!
!  Random argument accuracy tests.
!
  do j = 1, 4

    n = 2000
    xn = real ( n, kind = 8 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    a1 = zero
    r6 = zero
    r7 = zero
    del = ( b - a ) / xn
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl

  110     continue

      alpha = ren ( jt )
!
!  Carefully purify arguments.
!
      if ( j <= 1 ) then
        y = x / xlam
      else
        y = x - delta
      end if

      w = sixten * y
      t1 = w + y
      t1 = w + t1
      y = t1 - w
      y = y - w

      if ( j <= 1 ) then
        x = y * xlam
        mb = 15
      else
        x = y + delta
        mb = 2
      end if

      call rybesl ( y, alpha, mb, u2, ncalc )
      call rybesl ( x, alpha, mb, u, ncalc )
!
!  Accuracy test is based on the multiplication theorem.
!  First, filter out cases where function values are small.
!
      if ( j <= 1 ) then

        if ( sign ( one, u(1) ) * sign ( one, u2(1) ) &
          < zero ) then
          go to 110
        end if

        d = -f * y
        mb = ncalc - 1
        xmb = real ( mb, kind = 8 )
        sum = u2(mb+1)
        z = sum
        ind = mb

        do ii = 2, mb
          sum = sum * d / xmb + u2(ind)
          z = z * d / xmb
          ind = ind - 1
          xmb = xmb - one
        end do

        z = z * d
!
!  Check for convergence.
!
        if ( abs ( z / u2(ind) ) > eps ) then
          xl = xl + del
          go to 200
        end if

        zz = sum * d + u2(ind)
!
!  Check for numerical stability.
!
        d = abs ( zz / u2(ind) )

        if ( d > onep25 .or. d < p875 ) then
          go to 110
        end if

        zz = zz * xlam**alpha
!
!  Accuracy test is based on local Taylor's series expansion.
!  First, filter out cases where function values or derivatives
!  are small.
!
      else

        w = min ( abs ( u(1) ), abs ( u2(1) ), abs ( u2(2) ) )

        if ( w < sqrt ( twobpi / x ) / sixten ) then
          go to 110
        end if

        if ( abs ( u(1) ) < abs ( u2(1) ) ) then

          z = x
          x = y
          y = z
          delta = x - y

          do ii = 1, 9
            z = u(ii)
            u(ii) = u2(ii)
            u2(ii) = z
          end do

        end if

        ysq = y * y
        alphsq = alpha * alpha
        mb = 8
        j1 = mb
        xj1 = real ( j1 + 1, kind = 8 )
        iexp = 0
        nk = 13
        num = 2

        do ii = 1, mb

          if ( nk == 0 ) then
            nk = 11
            num = 1
          end if

          k = 9 - j1

          if ( k > 1 ) then
            no1 = ndx2(k-1) + 1
          else
            no1 = 1
          end if

          mvr = no1
          last = ndx2(k)
          k = last - no1 + 1
!
!  Group I(ALPHA) terms in the derivative.
!
          do iii = 1, k

            j2 = ndx(mvr)

            if ( num == 1 ) then
              g(iii) = ar1(nk,j2)
            else
              g(iii) = ar2(nk,j2)
            end if

            if ( j2 > 1 ) then

157               continue

              j2 = j2 - 1

              if ( num == 1 ) then
                g(iii) = g(iii) * alpha + ar1(nk,j2)
              else
                g(iii) = g(iii) * alpha + ar2(nk,j2)
              end if

              if ( j2 > 1 ) then
                go to 157
              end if

            end if

            mvr = mvr + 1
            nk = nk - 1

          end do

          t1 = g(1)

          do iii = 2, k
            t1 = t1 / ysq + g(iii)
          end do

          if ( iexp == 1 ) then
            t1 = t1 / y
          end if
!
!  Group I(ALPHA+1) terms in the derivative.
!
          iexp = 1 - iexp
          nk = nk + k
          mvr = no1
          kk = k

          do iii = 1, k

            j2 = ndx(mvr)
            m = mod ( j2, 2 )

            if ( m == 1 ) then
              j2 = j2 - 1
            end if

            if ( j2 >= 2 ) then

              if ( num == 1 ) then
                g(iii) = ar1(nk,j2)
              else
                g(iii) = ar2(nk,j2)
              end if

  163             continue

              j2 = j2 - 2

              if ( j2 >= 2 ) then
                if ( num == 1 ) then
                  g(iii) = g(iii) * alphsq + ar1(nk,j2)
                else
                  g(iii) = g(iii) * alphsq + ar2(nk,j2)
                end if
                go to 163
              end if

            else

              kk = iii - 1

            end if

            mvr = mvr + 1
            nk = nk - 1

          end do

          t2 = g(1)
          do iii = 2, kk
            t2 = t2 / ysq + g(iii)
          end do

          if ( iexp == 1 ) then
            t2 = t2 / y
          end if

          deriv = u2(1) * t1 - u2(2) * t2

          if ( j1 == 8 ) then
            sum = deriv
          else
            sum = sum * delta / xj1 + deriv
          end if

          j1 = j1 - 1
          xj1 = xj1 - one

        end do

        zz = sum * delta + u2(1)

      end if

      mb = 2
      z = u(1)
!
!  Accumulate results.
!
      w = ( z - zz ) / z

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
        a1 = alpha
        fxmx = z
      end if

      r7 = r7 + w * w
      xl = xl + del
  200     continue

    end do
!
!  Gather and print statistics for test.
!
    n = k1 + k2 + k3
    r7 = sqrt ( r7 / xn )

    if ( j <= 1 ) then
      write ( *,1000)
    else
      write ( *,1001)
    end if

    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3

    if ( n /= 2000 ) then
      write ( *,1012) 2000 - n
    end if

    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( r6 ) / albeta
    else
      w = x99
    end if

    write ( *,1021) r6,ibeta,w,x1,a1
    write ( *,1024) fxmx
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( r7 ) / albeta
    else
      w = x99
    end if

    write ( *,1023) r7,ibeta,w
    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b

    if ( j == 1 ) then
      b = ten
    else if ( j == 2 ) then
      b = b + b
    else if ( j == 3 ) then
      a = b + ten
      b = a + ten
    end if

  end do
!
!  Test of error returns.
!
!  First, test with bad parameters.
!
  write ( *, 2006)
  x = one
  alpha = one + half
  mb = 5
  call rybesl ( x, alpha, mb, u, ncalc )
  write ( *, 2011) x,alpha,mb,u(1),ncalc
  alpha = half
  mb = -mb
  call rybesl ( x, alpha, mb, u, ncalc )
  write ( *, 2011) x,alpha,mb,u(1),ncalc
!
!  Tests with small parameters.
!
  if ( xmin * xmax > one ) then
    x = xmin
  else
    x = one / xmax
  end if

  alpha = zero
  mb = 2
  call rybesl ( x, alpha, mb, u, ncalc )
  write ( *, 2011) x,alpha,mb,u(1),ncalc
  x = x + x + x
  mb = 2
  call rybesl ( x, alpha, mb, u, ncalc )
  write ( *, 2011) x,alpha,mb,u(1),ncalc
  alpha = one - eps
  mb = 2
  call rybesl ( x, alpha, mb, u, ncalc )
  write ( *, 2011) x,alpha,mb,u(1),ncalc
!
!  Last tests are with large parameters.
!
  write ( *, 2015)
  x = half / sqrt ( eps )
  alpha = half
  mb = 2
  call rybesl ( x, alpha, mb, u, ncalc )
  write ( *, 2012) x,alpha
  write ( *, 2014) ncalc,u(1)
  x = x * sixten
  mb = 2
  call rybesl ( x, alpha, mb, u, ncalc )
  write ( *, 2012) x,alpha
  write ( *, 2013)
  write ( *, 2014) ncalc,u(1)
  write ( *, 2020)
  return
 1000 format('1Test of Y(X,ALPHA) vs Multiplication Theorem'//)
 1001 format('1Test of Y(X,ALPHA) vs Taylor series'//)
 1010 format(I7,' Random arguments were tested from the interval ', &
      '(',F5.2,',',F5.2,')'//)
 1011 format(' Y(X,ALPHA) was larger',I6,' times,'/ &
      15X,' agreed',I6,' times, and'/ &
      11X,'was smaller',I6,' times.'//)
 1012 format(' NOTE: first ',I3,' arguments in test interval skipped'/ &
      7x,'because multiplication theorem did not converge'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number'//)
 1021 format(' The maximum relative error of',E15.4,' = ',I4,' **', &
      F7.2/4X,'occurred for X =',E13.6,' and NU =',E13.6)
 1022 format(' The estimated loss of base',I4, &
      ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      ' = ',I4,' **',F7.2)
 1024 format(4x,'with Y(X,ALPHA) = ',E13.6)
 2006 format('1Check of Error Returns'/// &
      ' The following summarizes calls with indicated parameters'// &
      ' NCALC different from MB indicates some form of error'// &
      ' See documentation for RYBESL for details'// &
      7X,'ARG',12X,'ALPHA',6X,'MB',6X,'B(1)',6X,'NCALC'//)
 2011 format(2E15.7,I5,E15.7,I5//)
 2012 format(' RYBESL will be called with the arguments',2E13.6)
 2013 format(' This should trigger an error message.')
 2014 format(' NCALC returned the value',I5/ &
      ' and RYBESL returned U(1) = ',E13.6/)
 2015 format(' Tests near the largest acceptable argument for RYBESL'/)
 2020 format(' This concludes the tests.')
end
subroutine y0_test ( )

!*****************************************************************************80
!
!! Y0_TEST tests BESY0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) aleps
  real ( kind = 8 ) all9
  real ( kind = 8 ) b
  real ( kind = 8 ) besy0
  real ( kind = 8 ) besy1
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) del
  real ( kind = 8 ) eight
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) five5
  real ( kind = 8 ) half
  real ( kind = 8 ) hund
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) one
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  logical sflag
  real ( kind = 8 ) sgn
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  logical tflag
  real ( kind = 8 ) three
  real ( kind = 8 ) twenty
  real ( kind = 8 ) two56
  real ( kind = 8 ) u(560)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xi(3)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) yx(3)
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz

  data zero / 0.0D+00 /
  data one / 1.0D+00 /
  data three / 3.0d0 /
  data five5 /5.5d0/
  data eight / 8.0d0 /
  data twenty / 20.0d0 /
  data all9 / -999.0d0 /
  data two56 / 256.0d0/
  data half /0.5d0 /
  data sixten /16.0d0 /
  data xlam /0.9375d0 /
  data hund /100.0d0/

  data xi / &
    228.0d0, &
    1013.0d0, &
    1814.0d0 /

  data yx / &
    -2.6003142722933487915d-3, &
     2.6053454911456774983d-4, &
    -3.4079448714795552077d-5 /
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  aleps = log ( eps )
  ait = real ( it, kind = 8 )
  jt = 0
  a = eps
  b = three
  n = 2000
!
!  Random argument accuracy tests based on the multiplication theorem.
!
  do j = 1, 4

    sflag = ( j == 1 ) .and. ( maxexp / it <= 5 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    sgn = one
    del = ( b - a ) / real ( n, kind = 8 )
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Carefully purify arguments.
!
      y = x / xlam
      w = sixten * y
      y = ( w + y ) - w
      x = y * xlam
!
!  Generate Bessel functions with forward recurrence.
!
      u(1) = besy0 ( y )
      u(2) = besy1 ( y )
      tflag = sflag .and. ( y < half )

      if ( tflag ) then
        u(1) = u(1) * eps
        u(2) = u(2) * eps
      end if

      mb = 1
      xmb = one
      y = y * half
      w = ( one - xlam ) * ( one + xlam )
      c = w * y
      t = abs ( u(1) + c * u(2) )
      t1 = eps / hund

      do ii = 3, 60

        z = abs ( u(ii-1) )

        if ( z / t1 < t ) then
          go to 120
        end if

        if ( y < xmb ) then
          if ( z > xmax * ( y / xmb ) ) then
            a = x
            xl = xl + del
            go to  200
          end if
        end if

        u(ii) = xmb / y * u(ii-1) - u(ii-2)

        if ( t1 > one / eps ) then
          t = t * t1
          t1 = one
        end if

        t1 = xmb * t1 / c
        xmb = xmb + one
        mb = mb + 1

      end do
!
!  Evaluate Bessel series expansion.
!
120       continue

      sum = u(mb)
      ind = mb

      do ii = 2, mb
        ind = ind - 1
        xmb = xmb - one
        sum = sum * w * y / xmb + u(ind)
      end do

      zz = sum

      if ( tflag ) then
        zz = zz / eps
      end if

      z = besy0 ( x )
      y = z

      if ( abs ( u(1) ) > abs ( y ) ) then
        y = u(1)
      end if
!
!  Accumulate results.
!
      w = ( z - zz ) / y

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

  200     continue

    end do
!
!  Gather and print statistics.
!
    n = k1 + k2 + k3
    r7 = sqrt ( r7 / real ( n, kind = 8 ) )
    write ( *,1000)
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta
    w = all9

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    end if

    if ( j == 4 ) then
      write ( *,1024) r6,ibeta,w,x1
    else
      write ( *,1021) r6,ibeta,w,x1
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
    w = all9

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    end if

    if ( j == 4 ) then
      write ( *,1025) r7,ibeta,w
    else
      write ( *,1023) r7,ibeta,w
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w
!
!  Initialize for next test.
!
    a = b

    if ( j == 1 ) then
      b = five5
      n = 2000
    else if ( j == 2 ) then
      b = eight
      n = 2000
    else
      b = twenty
      n = 500
    end if

  end do
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1031) ibeta

  do i = 1, 3
    x = xi(i) / two56
    y = besy0 ( x )
    w = all9
    t = ( y - yx(i) ) / yx(i)
    if ( t /= zero ) then
      w = log ( abs ( t ) ) / albeta
    end if

    w = max ( ait + w, zero )
    write ( *,1032) x,y,w
  end do
!
!  Test of error returns.
!
  write ( *,1033)
  x = xmin
  write ( *,1035) x
  y = besy0 ( x )
  write ( *,1036) y
  x = zero
  write ( *,1034) x
  y = besy0 ( x )
  write ( *,1036) y
  x = xmax
  write ( *,1034) x
  y = besy0 ( x )
  write ( *,1036) y
  write ( *, '(a)' ) '  This concludes the tests.'
  return
 1000 format('1Test of Y0(X) VS Multiplication Theorem'  //)
 1010 format(I7,' random arguments were tested from the interval ', &
   1H(,F5.1,1H,,F5.1,1H)//)
 1011 format(' ABS(Y0(X)) was larger',I6,' times', / &
     15X,' agreed',I6,' times, and'/ &
     11X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number.' //)
 1021 format(' The maximum relative error of',E15.4,3H = ,I4,3H **, &
    F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
    ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      3H = ,I4,3H **,F7.2)
 1024 format(' The maximum absolute error of',E15.4,3H = ,I4,3H **, &
    F7.2/4X,'occurred for X =',E13.6)
 1025 format(' The root mean square absolute error was',E15.4, &
      3H = ,I4,3H **,F7.2)
 1031 format(' Accuracy near zeros'//10X,'X',15X,'BESY0(X)', &
      13X,'Loss of base',I3,' digits'/)
 1032 format(E20.10,E25.15,8X,F7.2/)
 1033 format(//' Test with extreme arguments'/)
 1034 format(' Y0 will be called with the argument ',E17.10/ &
       ' This may stop execution.'//)
 1035 format(' Y0 will be called with the argument ',E17.10/ &
       ' This should not stop execution.'//)
 1036 format(' Y0 returned the value',E25.17/)
end
subroutine y1_test ( )

!*****************************************************************************80
!
!! Y1_TEST tests BESY1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, George Zazi.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Laura Stoltz,
!    Performance evaluation of programs for certain Bessel functions,
!    ACM Transactions on Mathematical Software,
!    Volume 15, Number 1, March 1989, pages 41-48.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ait
  real ( kind = 8 ) albeta
  real ( kind = 8 ) aleps
  real ( kind = 8 ) all9
  real ( kind = 8 ) b
  real ( kind = 8 ) besy0
  real ( kind = 8 ) besy1
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) del
  real ( kind = 8 ) eight
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  real ( kind = 8 ) four
  real ( kind = 8 ) half
  real ( kind = 8 ) hund
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) mb
  integer ( kind = 4 ) mbm1
  integer ( kind = 4 ) mbp1
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndum
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  logical sflag
  logical tflag
  real ( kind = 8 ) one
  real ( kind = 8 ) ren
  real ( kind = 8 ) r6
  real ( kind = 8 ) r7
  real ( kind = 8 ) sgn
  real ( kind = 8 ) sixten
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) twenty
  real ( kind = 8 ) two56
  real ( kind = 8 ) t1
  real ( kind = 8 ) u(560)
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xi(2)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xlam
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmb
  real ( kind = 8 ) xmin
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) yx(2)
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  real ( kind = 8 ) zz
!
!  Mathematical constants
!
  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data four / 4.0D+00/
  data eight  /8.0d0 /
  data twenty /20.0d0 /
  data all9 /-999.0d0 /
  data two56 /256.0d0/
  data sixten / 16.0d0 /
  data xlam / 0.9375d0/
  data hund /100.0d0/
!
!  Zeroes of Y1
!
  data xi/562.0d0,1390.0d0/
  data yx(1)/-9.5282393097722703132d-4/, &
       yx(2)/-2.1981830080598002741d-6/
!
!  Determine machine parameters and set constants.
!
  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, &
     maxexp, eps, epsneg, xmin, xmax )

  beta = real ( ibeta, kind = 8 )
  albeta = log ( beta )
  aleps = log ( eps )
  ait = real ( it, kind = 8 )
  jt = 0
  b = eps
!
!  Random argument accuracy tests based on the multiplication theorem.
!
  do j = 1, 3

    sflag = ( j == 1 ) .and. ( maxexp / it <= 5 )
    k1 = 0
    k2 = 0
    k3 = 0
    x1 = zero
    r6 = zero
    r7 = zero
    n = 2000
    a = b

    if ( j == 1 ) then
      b = four
    else if ( j == 2 ) then
      b = eight
    else
      b = twenty
      n = 500
    end if

    sgn = one
    del = ( b - a ) / real ( n, kind = 8 )
    xl = a

    do i = 1, n

      x = del * ren ( jt ) + xl
!
!  Carefully purify arguments and evaluate identities.
!
      y = x / xlam
      w = sixten * y
      y = ( w + y ) - w
      x = y * xlam
      u(1) = besy0 ( y )
      u(2) = besy1 ( y )
      tflag = sflag .and. ( y < half )

      if ( tflag ) then
        u(1) = u(1) * eps
        u(2) = u(2) * eps
      end if

      mb = 1
      xmb = one
      y = y * half
      w = ( one - xlam ) * ( one + xlam )
      c = w * y
      t = abs ( xlam * u(2) )
      t1 = eps / hund

      do ii = 3, 60

        u(ii) = xmb / y * u(ii-1) - u(ii-2)
        t1 = xmb * t1 / c
        xmb = xmb + one
        mb = mb + 1
        mbp1 = mb + 1
        z = abs ( u(ii) )

        if ( z / t1 < t ) then
          go to 120
        else if ( y < xmb ) then
          if ( z > xmax * ( y / xmb ) ) then
            a = x
            xl =xl + del
            go to  200
          end if
        end if

        if ( t1 > one / eps ) then
          t = t * t1
          t1 = one
        end if

      end do

  120     continue

      sum = u(mbp1)
      ind = mbp1
      mbm1 = mb - 1

      do ii = 1, mbm1
        ind = ind - 1
        xmb = xmb - one
        sum = sum * w * y / xmb + u(ind)
      end do

      zz = sum * xlam
      if ( tflag ) then
        zz = zz / eps
      end if

      z = besy1 ( x )
      y = z

      if ( abs ( u(2) ) > abs ( y ) ) then
        y =  u(2)
      end if
!
!  Accumulate results.
!
      w = ( z - zz ) / y

      if ( w > zero ) then
        k1 = k1 + 1
      else if ( w < zero ) then
        k3 = k3 + 1
      else
        k2 = k2 + 1
      end if

      w = abs ( w )

      if ( w > r6 ) then
        r6 = w
        x1 = x
      end if

      r7 = r7 + w * w
      xl = xl + del

  200     continue

    end do
!
!  Process and output statistics.
!
    n = k1 + k2 + k3
    r7 = sqrt ( r7 / real ( n, kind = 8 ) )
    write ( *,1000)
    write ( *,1010) n,a,b
    write ( *,1011) k1,k2,k3
    write ( *,1020) it,ibeta

    if ( r6 /= zero ) then
      w = log ( abs ( r6 ) ) / albeta
    else
      w = all9
    end if

    if ( j == 4 ) then
      write ( *,1024) r6,ibeta,w,x1
    else
      write ( *,1021) r6,ibeta,w,x1
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

    if ( r7 /= zero ) then
      w = log ( abs ( r7 ) ) / albeta
    else
      w = all9
    end if

    if ( j == 4 ) then
      write ( *,1025) r7,ibeta,w
    else
      write ( *,1023) r7,ibeta,w
    end if

    w = max ( ait + w, zero )
    write ( *,1022) ibeta,w

  end do
!
!  Special tests.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special Tests:'
  write ( *, '(a)' ) ' '
  write ( *,1031) ibeta

  do i = 1, 2

    x = xi(i) / two56
    y = besy1 ( x )
    t = ( y - yx(i) ) / yx(i)

    if ( t /= zero ) then
      w = log ( abs ( t ) ) / albeta
    else
      w = all9
    end if

    w = max ( ait + w, zero )
    write ( *,1032) x,y,w

  end do
!
!  Test of error returns.
!
  write ( *,1033)
  x = xmin

  if ( xmin * xmax < one ) then
    write ( *,1034) x
  else
    write ( *,1035) x
  end if

  y = besy1 ( x )
  write ( *,1036) y
  x = -one
  write ( *,1034) x
  y = besy1 ( x )
  write ( *,1036) y
  x = xmax
  write ( *,1034) x

  if ( .false. ) then
    y = besy1 ( x )
    write ( *,1036) y
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YUP, it stopped execution all right.'
    write ( *, '(a)' ) '  So let''s skip this test now!'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This concludes the tests.'

  return

 1000 format('1Test of Y1(X) VS Multiplication Theorem'  //)
 1010 format(I7,' random arguments were tested from the interval ', &
   1H(,F5.1,1H,,F5.1,1H)//)
 1011 format(' ABS(Y1(X)) was larger',I6,' times', / &
       15X,' agreed',I6,' times, and'/ &
     11X,'was smaller',I6,' times.'//)
 1020 format(' There are',I4,' base',I4, &
      ' significant digits in a floating-point number.' //)
 1021 format(' The maximum relative error of',E15.4,3H = ,I4,3H **, &
    F7.2/4X,'occurred for X =',E13.6)
 1022 format(' The estimated loss of base',I4, &
    ' significant digits is',F7.2//)
 1023 format(' The root mean square relative error was',E15.4, &
      3H = ,I4,3H **,F7.2)
 1024 format(' The maximum absolute error of',E15.4,3H = ,I4,3H **, &
    F7.2/4X,'occurred for X =',E13.6)
 1025 format(' The root mean square absolute error was',E15.4, &
      3H = ,I4,3H **,F7.2)
 1031 format(' Accuracy near zeros'//10X,'X',15X,'BESY1(X)', &
      13X,'Loss of base',I3,' digits'/)
 1032 format(E20.10,E25.15,8X,F7.2/)
 1033 format(//' Test with extreme arguments'///)
 1034 format(' Y1 will be called with the argument ',E17.10/ &
       ' This may stop execution.'//)
 1035 format(' Y1 will be called with the argument ',E17.10/ &
       ' This should not stop execution.'//)
 1036 format(' Y1 returned the value',E25.17/)
end
