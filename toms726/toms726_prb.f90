program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS726_PRB.
!
!  Modified:
!
!    11 January 2008
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS726_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS726 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
  call test11 ( )

  call test99 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS726_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01
!
!  Discussion:
!
!    This test generates the first N beta-coefficients in the recurrence
!    relation for the orthogonal polynomials relative to the weight
!    function
!
!      ((1-om2*x**2)*(1-x**2))**(-1/2)  on (-1,1)
!
!    for om2 = .1(.2).9,.99,.999, both in single and real ( kind = 8 ),
!    using modified moments if MODMOM = true and ordinary moments
!    otherwise.  In the former case, N = 80, in the latter, N = 20.  Printed
!    are the double-precision values of the coefficients along with the
!    relative errors of the single-precision values.
!
!    This test relates to Example 3.1 of the companion paper, where
!    orthogonal polynomials are generated relative to a weight
!    function on (-1,1) having square root singularities at
!    1, -1, 1/omega, -1/omega, with  omega  between 0 and 1.
!
!  Modified:
!
!    11 January 2008
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
  implicit none

  real a(159)
  real alpha(80)
  real b(159)
  real beta(80)
  real ( kind = 8 ) d(80)
  real ( kind = 8 ) d0(80)
  real ( kind = 8 ) da(159)
  real ( kind = 8 ) dalpha(80)
  real ( kind = 8 ) db(159)
  real ( kind = 8 ) dbeta(80)
  real ( kind = 8 ) deps
  real ( kind = 8 ) dnu(160)
  real ( kind = 8 ) dom2
  real ( kind = 8 ), dimension ( 7 ) :: doom2 = (/ &
    0.1D+00, 0.3D+00, 0.5D+00, 0.7D+00, 0.9D+00, &
    0.99D+00, 0.999D+00  /)
  real ( kind = 8 ) drr(80)
  real ( kind = 8 ) ds(80)
  real eps
  real errb
  real f(80)
  real f0(80)
  real fnu(160)
  integer iderr
  integer ierr
  integer iom
  integer k
  logical modmom
  real om2
  integer n
  real rr(80)
  real s(80)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'

  modmom = .true.
  eps = epsilon ( eps )
  deps = epsilon ( deps )

  if ( modmom ) then
    n = 80
  else
    n = 20
  end if

  do iom = 1, 7

    dom2 = doom2(iom)
    om2 = sngl(dom2)
!
!  Compute the modified or ordinary moments using Eqs. (3.7) and (3.9)
!  of the companion paper.  On machines with limited exponent range, some
!  of the high-order modified moments may underflow, without this having
!  any deteriorating effect on the accuracy.
!
    call mm_01_r4 ( n, eps, modmom, om2, fnu, ierr, f, f0, rr )

    call mm_01_r8 ( n, deps, modmom, dom2, dnu, iderr, d, d0, drr )

    if ( ierr /= 0 .or. iderr /= 0 ) then
      write(*,2) ierr,iderr,om2
    2     format(/5x,'ierr in fmm01 = ',i1,'  iderr in dmm01 = ',i1, &
    '  for om2 = ',f8.4/)
      cycle
    end if
!
!  Generate the recursion coefficients for the polynomials defining the
!  modified or ordinary moments.
!
    if ( modmom ) then
      call recur_r4 ( 2 * n - 1, 3, 0.0E+00, 0.0E+00, a, b, ierr )
      call recur_r8 ( 2 * n - 1, 3, 0.0D+00, 0.0D+00, da, db, iderr )
    else
      a(1:2*n-1) = 0.0E+00
      b(1:2*n-1) = 0.0E+00
      da(1:2*n-1) = 0.0D+00
      db(1:2*n-1) = 0.0D+00
    end if
!
!  Compute the desired recursion coefficients by means of the modified
!  Chebyshev algorithm.
!
    call cheb_r4 ( n, a, b, fnu, alpha, beta, s, ierr )
!
!  On machines with limited single-precision exponent range, the routine
!  CHEB_R4  may generate an underflow exception, which however is harmless
!  and can be ignored.
!
    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, iderr )

    write(*,3) ierr,iderr
    3   format(/5x,'ierr in CHEB_R4 = ',i3,'  iderr in CHEB_R8 = ',i3/)
    write(*,4)
    4   format(/5x,'k',14x,'dbeta(k)'/)

    do k = 1, n
      if ( iderr == 0 .or. k - 1 < abs ( iderr ) ) then
        if ( ierr == 0 .or. k - 1 < abs ( ierr ) ) then
          errb = sngl ( abs ( real ( beta(k), kind = 8 ) - dbeta(k) ) / dbeta(k) )
          if ( k == 1 ) then
            write(*,5) k - 1, dbeta(k),errb,om2
    5       format(1x,i5,d36.28,e12.4,'   om2  = ',f6.3)
          else
            write(*,6) k - 1, dbeta(k), errb
    6       format(1x,i5,d36.28,e12.4)
          end if
        else
          write(*,7) k - 1, dbeta(k)
    7     format(1x,i5,d36.28)
        end if
      end if
    end do

    write ( *, '(a)' ) ' '

  end do

  return
end
subroutine mm_01_r4 ( n, eps, modmom, om2, fnu, ierr, f, f0, rr )

!*****************************************************************************80
!
!! MM_01_R4 generates the modified or ordinary weights of the weight function.
!
!  Discussion:
!
!    The weight function is
!
!      ((1-om2*x**2) * (1-x**2))**(-1/2)  on (-1,1)
!
!    Equations (3.7) and (3.9) of the companion paper are used.
!
  integer n

  logical done
  real fnu(2*n)
  real f(n)
  real f0(n)
  integer ierr
  integer k
  logical modmom
  integer nd
  integer nu
  real pi
  real q
  real q1
  real rr(n)

  ierr = 0
  nd = 2 * n
  pi = 4.0E+00 * atan ( 1.0E+00 )
!
!  Compute the Fourier coefficients of ((1-om2*sin(theta)**2))**(-1/2)
!  as minimal solution of a three-term recurrence relation as described
!  on pp.310-311 of W. Gautschi,On generating orthogonal polynomials'',
!  SIAM J. Sci. Statist. Comput. 3, 1982, 289-317.
!
  q = om2 / ( 2.0E+00 - om2 + 2.0E+00 * sqrt ( 1.0E+00 - om2 ) )
  q1 = ( 1.0E+00 + q * q ) / q
  f(1:n) = 0.0E+00
  nu = 2 * n

  do

    nu = nu + 10
    f0(1:n) = f(1:n)

    if ( 500 < nu ) then
      ierr = 1
      return
    end if

    r = 0.0E+00
    s = 0.0E+00

    do k = 1, nu

      n1 = nu - k + 1
      fn1 = real ( n1, kind = 4 )
      r = - ( fn1 - 0.5E+00 ) / ( fn1 * q1 + ( fn1 + 0.5 ) * r )
      s = r * ( 2.0E+00 + s )

      if ( n1 <= n ) then
        rr(n1) = r
      end if

    end do

    c0 = 1.0E+00 / ( 1.0E+00 + s )

    f(1) = rr(1) * c0
    do k = 2, n
      f(k) = rr(k) * f(k-1)
    end do

    done = .true.

    do k = 1, n
      if ( abs ( f(k) - f0(k) ) > eps * abs ( f(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do
!
!  Compute the desired modified or ordinary moments in terms of
!  the above Fourier coefficients.
!
  fnu(1) = pi * c0
  if ( n == 1 ) then
    return
  end if

  fnu(2) = 0.0E+00
  if ( n == 2 ) then
    return
  end if

  if ( modmom ) then

    c = 2.0E+00 * pi

    do k = 3, 2 * n - 1, 2
      k1 = ( k - 1 ) / 2
      c = - 0.25E+00 * c
      fnu(k) = c * f(k1)
      fnu(k+1) = 0.0E+00
    end do

  else

    c = 0.5E+00 * pi
    fnu(3) = c * ( c0 - f(1) )
    fnu(4) = 0.0E+00
    c = - c

    do k = 5, 2 * n - 1, 2

      k1 = ( k - 1 ) / 2
      k1m1 = k1 - 1
      c = - 0.25E+00 * c
      c1 = 1.0E+00
      sum = f(k1)

      do i = 1, k1m1
        c1 = - c1 * real ( 2 * k1 - i + 1, kind = 4 ) / real ( i, kind = 4 )
        sum = sum + c1 * f(k1-i)
      end do

      c1 = - c1 * real ( k1 + 1, kind = 4 ) / real ( 2 * k1, kind = 4 )
      sum = sum + c1 * c0
      fnu(k) = c * sum
      fnu(k+1) = 0.0E+00

    end do

  end if

  return
end
subroutine mm_01_r8 ( n, deps, modmom, dom2, dnu, ierrd, d, d0, drr )

!*****************************************************************************80
!
!! MM_01_R8 generates the modified or ordinary weights of the weight function.
!
!  Discussion:
!
!    The weight function is
!
!      ((1-om2*x**2) * (1-x**2))**(-1/2)  on (-1,1)
!
!    Equations (3.7) and (3.9) of the companion paper are used.
!
  integer n

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) d0(n)
  real ( kind = 8 ) dc
  real ( kind = 8 ) dc0
  real ( kind = 8 ) dc1
  real ( kind = 8 ) deps
  real ( kind = 8 ) dn1
  real ( kind = 8 ) dnu(2*n)
  real ( kind = 8 ) dom2
  logical done
  real ( kind = 8 ) dpi
  real ( kind = 8 ) dq
  real ( kind = 8 ) dq1
  real ( kind = 8 ) dr
  real ( kind = 8 ) drr(n)
  real ( kind = 8 ) ds
  real ( kind = 8 ) dsum
  integer i
  integer ierrd
  integer k
  integer k1
  logical modmom
  integer n1
  integer nd
  integer nud

  ierrd = 0
  nd = 2 * n
  dpi = 4.0D+00 * atan ( 1.0D+00 )
  dq = dom2 / ( 2.0D+00 - dom2 + 2.0D+00 * sqrt ( 1.0D+00 - dom2 ) )
  dq1 = ( 1.0D+00 + dq * dq ) / dq
  d(1:n) = 0.0D+00
  nud = 2 * n

  do

    nud = nud + 10

    d0(1:n) = d(1:n)

    if ( 1000 < nud ) then
      ierrd = 1
      return
    end if

    dr = 0.0D+00
    ds = 0.0D+00

    do k = 1, nud
      n1 = nud - k + 1
      dn1 = real ( n1, kind = 8 )
      dr = - ( dn1 - 0.5D+00 ) / ( dn1 * dq1 + ( dn1 + 0.5D+00 ) * dr )
      ds = dr * ( 2.0D+00 + ds )
      if ( n1 <= n ) then
        drr(n1) = dr
      end if
    end do

    dc0 = 1.0D+00 / ( 1.0D+00 + ds )
    d(1) = drr(1) * dc0

    do k = 2, n
      d(k) = drr(k) * d(k-1)
    end do

    done = .true.

    do k = 1, n
      if ( abs ( d(k) - d0(k) ) > deps * abs ( d(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do

  dnu(1) = dpi * dc0

  if ( n == 1 ) then
    return
  end if

  dnu(2) = 0.0D+00

  if ( n == 2 ) then
    return
  end if

  if ( modmom ) then

    dc = 2.0D+00 * dpi

    do k = 3, 2 * n - 1, 2
      k1 = ( k - 1 ) / 2
      dc = - 0.25D+00 * dc
      dnu(k) = dc * d(k1)
      dnu(k+1) = 0.0D+00
    end do

  else

    dc = 0.5D+00 * dpi
    dnu(3) = dc * ( dc0 - d(1) )
    dnu(4) = 0.0D+00
    dc = - dc
    do k = 5, 2 * n - 1, 2
      k1 = ( k - 1 ) / 2
      k1m1 = k1 - 1
      dc = - 0.25D+00 * dc
      dc1 = 1.0D+00
      dsum = d(k1)
      do i = 1, k1m1
        dc1 = - dc1 * real ( 2 * k1 - i + 1, kind = 8 ) / real ( i, kind = 8 )
        dsum = dsum + dc1 * d(k1-i)
      end do
      dc1 = - dc1 * real ( k1 + 1, kind = 8 ) / real ( 2 * k1, kind = 8 )
      dsum = dsum + dc1 * dc0
      dnu(k) = dc * dsum
      dnu(k+1) = 0.0D+00
    end do

  end if

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 works with a weight function X^SIGMA * Log ( X ).
!
!  Discussion:
!
!    This test relates to Example 3.2, where orthogonal polynomials are
!    generated relative to a weight function on (0,1) having
!    a logarithmic singularity at the origin as well as an
!    algebraic singularity with exponent  sigma  greater than -1.
!
!    This test generates the first n recursion coefficients for the
!    orthogonal polynomials relative to the weight function
!
!      x**sigma * ln(1/x)  on (0,1],  sigma = -.5, 0, .5,
!
!    where n = 100 when using modified (Legendre) moments, and n = 12 when
!    using ordinary moments. It prints the double-precision values of the
!    coefficients as well as the relative errors of the respective single-
!    precision values and the maximum relative errors.
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
  real a(199)
  real alpha(100)
  real b(199)
  real beta(100)
  real ( kind = 8 ) da(199)
  real ( kind = 8 ) dalpha(100)
  real ( kind = 8 ) db(199)
  real ( kind = 8 ) dbeta(100)
  real ( kind = 8 ) dnu(200)
  real ( kind = 8 ) ds(100)
  real ( kind = 8 ) dsigma
  real fnu(200)
  integer ierr
  logical intexp
  logical modmom
  integer n
  real s(100)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Weight function = x**sigma * ln(1/x)'

  modmom = .true.
!
!  Generate the recursion coefficients for the polynomials defining the
!  modified or ordinary moments.
!
  if ( modmom ) then
    n = 100
    call recur_r4 ( 2 * n - 1, 2, 0.0E+00, 0.0E+00, a, b, ierr )
    call recur_r8 ( 2 * n - 1, 2, 0.0D+00, 0.0D+00, da, db, iderr )
  else
    n = 12
    a(1:2*n-1) = 0.0E+00
    b(1:2*n-1) = 0.0E+00
    da(1:2*n-1) = 0.0D+00
    db(1:2*n-1) = 0.0D+00
  end if

  do is = 1, 3

    dsigma = - 0.5D+00 + 0.5D+00 * real ( is - 1, kind = 8 )
    sigma = real ( dsigma, kind = 4 )

    if ( is == 2 ) then
      intexp = .true.
    else
      intexp = .false.
    end if
!
!  Compute the modified or ordinary moments using Eqs. (3.12) and
!  (3.11) of the companion paper.  On machines with limited exponent
!  range, some of the high-order modified moments may underflow, without
!  this having any deteriorating effect on the accuracy.
!
    call mm_02_r4 ( n, modmom, intexp, sigma, fnu )
    call mm_02_r8 ( n, modmom, intexp, dsigma, dnu )
!
!  Compute the desired recursion coefficients by means of the modified
!  Chebyshev algorithm.
!
    call cheb_r4 ( n, a, b, fnu, alpha, beta, s, ierr )
!
!  On machines with limited single-precision exponent range, the routine
!  CHEB_R4 may generate an underflow exception, which however is harmless
!  and can be ignored.
!
    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, iderr )

    write(*,1) ierr,iderr
    1   format(/6x,'ierr in CHEB_R4 = ',i4,' iderr in CHEB_R8 = ',i4/)
!
!  Compute and print the relative errors and their maxima.
!
    eamax = 0.0E+00
    ebmax = 0.0E+00

    write ( *, '(a)' ) ' '
    write ( *, '(a,f5.1)' ) '  SIGMA = ', sigma
    write ( *, '(a)' ) ' '
    write(*,3)
    3   format(3x,'k',13x,'dalpha(k)',25x,'dbeta(k)'/)

    do k = 1, n

      if ( iderr == 0 .or. k - 1 < abs ( iderr ) ) then

        write(*,4) k - 1, dalpha(k), dbeta(k)
    4   format(1x,i3,2d33.25)

        if ( ierr == 0 .or. k - 1 < abs ( ierr ) ) then

          erra = sngl ( abs( real ( alpha(k), kind = 8 ) - dalpha(k) ) / dalpha(k))
          errb = sngl ( abs( real ( beta(k), kind = 8 ) - dbeta(k) ) / dbeta(k))
          write(*,5) erra,errb
    5     format(4x,e12.4,21x,e12.4)

          if ( eamax < erra ) then
            eamax = erra
            kamax = k - 1
          end if

          if ( ebmax < errb ) then
            ebmax = errb
            kbmax = k - 1
          end if

        end if
      end if
    end do

    write(*,6) eamax,kamax,ebmax,kbmax
    6   format(/6x,'eamax  = ',e11.4,' at',i3,9x,'ebmax  = ',e11.4, &
      ' at',i3//)

  end do

  return
end
subroutine mm_02_r4 ( n, modmom, intexp, sigma, fnu )

!*****************************************************************************80
!
!! MM_02_R4 generates ordinary or modified moments for a particular weight function.
!
!  Discussion:
!
!    The routine generates the first  2*n  modified moments (if modmom = .true.)
!    relative to shifted monic Legendre polynomials, using Eq. (3.12) of
!    the companion paper, and the first  2*n  ordinary moments (if modmom
!     = .false.) by Eq. (3.11), of the weight function
!
!      (x**sigma)*ln(1/x)  on (0,1],   sigma > -1,
!
!    for sigma an integer (if intexp = .true.) or a real number (if intexp
!     = .false.). In either case, the input variable  sigma  is of type real.
!
  integer n

  real c
  real fi
  real fiq
  real fk
  real fkm1
  real fnu(2*n)
  logical intexp
  integer iq
  integer isigma
  integer k
  integer kmax
  logical modmom
  integer nd
  real p
  real q
  real s
  real sigma
  real sigp1

  nd = 2 * n
  sigp1 = sigma + 1.0E+00

  if ( modmom ) then

    isigma = int ( sigma )

    if ( intexp .and. isigma + 1 < 2 * n ) then
      kmax = isigma + 1
    else
      kmax = 2 * n
    end if

    c = 1.0E+00

    do k = 1, kmax
      fk = real ( k, kind = 4 )
      p = 1.0E+00
      s = 1.0E+00 / sigp1
      if ( 1 < kmax ) then
        do i = 1, k - 1
          fi = real ( i, kind = 4 )
          p = ( sigp1 - fi ) * p / ( sigp1 + fi )
          s = s + 1.0E+00 / ( sigp1 + fi ) - 1.0E+00 / ( sigp1 - fi )
        end do
      end if
      fnu(k) = c * s * p / sigp1
      c = fk * c / ( 4.0E+00 * fk - 2.0E+00 )
    end do

    if ( .not. intexp .or. 2 * n <= isigma + 1 ) then
      return
    end if

    q = - 0.5E+00

    do iq = 1, isigma
      fiq = real ( iq, kind = 4 )
      q = fiq * fiq * q / ( ( 2.0E+00 * fiq + 1.0E+00 ) * ( 2.0E+00 * fiq + 2.0E+00 ) )
    end do

    fnu(isigma+2) = c * q

    do k = isigma + 3, 2 * n
      fkm1 = real ( k - 1, kind = 4 )
      fnu(k) = - fkm1 * ( fkm1 - sigp1 ) * fnu( k - 1 ) &
        / ( ( 4.0E+00 * fkm1 - 2.0E+00 ) * ( fkm1 + sigp1 ) )
    end do

    return

  else

    do k = 1, 2 * n
      fnu(k) = ( 1.0E+00 / ( sigp1 + real ( k - 1, kind = 4 ) ) )**2
    end do

  end if

  return
end
subroutine mm_02_r8 ( n, modmom, intexp, dsigma, dnu )

!*****************************************************************************80
!
!! MM_02_R8 is a double-precision version of MM_02_R4.
!
  integer n

  real ( kind = 8 ) dc
  real ( kind = 8 ) di
  real ( kind = 8 ) diq
  real ( kind = 8 ) dk
  real ( kind = 8 ) dkm1
  real ( kind = 8 ) dnu(2*n)
  real ( kind = 8 ) dp
  real ( kind = 8 ) dq
  real ( kind = 8 ) ds
  real ( kind = 8 ) dsigma
  real ( kind = 8 ) dsigp1
  integer i
  logical intexp
  integer iq
  integer isigma
  integer k
  integer kmax
  logical modmom
  integer nd

  nd = 2 * n

  dsigp1 = dsigma + 1.0D+00

  if ( modmom) then

    isigma = int ( dsigma )

    if ( intexp .and. isigma + 1 < 2 * n ) then
      kmax = isigma + 1
    else
      kmax = 2 * n
    end if

    dc = 1.0D+00
    do k = 1, kmax
      dk = real ( k, kind = 8 )
      dp = 1.0D+00
      ds = 1.0D+00 / dsigp1
      if ( 1 < kmax ) then
        do i = 1, k - 1
          di = real ( i, kind = 8 )
          dp = ( dsigp1 - di ) * dp / ( dsigp1 + di )
          ds = ds + 1.0D+00 / ( dsigp1 + di ) - 1.0D+00 / ( dsigp1 - di )
        end do
      end if
      dnu(k) = dc * ds * dp / dsigp1
      dc = dk * dc / ( 4.0D+00 * dk - 2.0D+00 )
    end do

    if ( .not. intexp .or. 2 * n <= isigma + 1 ) then
      return
    end if

    dq = - 0.5D+00

    do iq = 1, isigma
      diq = real ( iq, kind = 8 )
      dq = diq * diq * dq / ( ( 2.0D+00 * diq + 1.0D+00 ) * ( 2.0D+00 * diq + 2.0D+00 ) )
    end do

    dnu(isigma+2) = dc * dq

    do k = isigma + 3, 2 * n
      dkm1 = real ( k - 1, kind = 8 )
      dnu(k) = - dkm1 * ( dkm1 - dsigp1 ) * dnu(k-1) &
        / ( ( 4.0D+00 * dkm1 - 2.0D+00 ) * ( dkm1 + dsigp1 ) )
    end do

    return

  else

    do k = 1, 2 * n
      dkm1 = real ( k - 1, kind = 8 )
      dnu(k) = ( 1.0D+00 / ( dsigp1 + dkm1 ) )**2
    end do

  end if

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 applies the Stieltjes procedure and Lanczos algorithm.
!
!  Discussion:
!
!    This test applies both the Stieltjes procedure and the Lanczos algorithm (cf.
!    W.B. Gragg and W.J. Harrod, The numerically stable reconstruction of
!    Jacobi matrices from spectral data'', Numer. Math. 44, 1984, 317-335)
!    to generate the  N  recursion coefficients, N = 40, 80, 160, 320,
!    for the (monic) orthogonal polynomials relative to the discrete inner
!    product supported on  N  equally spaced points on [-1,1] (including
!    the end points) and having equal weights 2/N. The routine prints the
!    absolute errors in the alpha-coefficients and the relative errors in
!    the beta-coefficients, these being computed using known formulae for
!    the coefficients. The maxima of these errors are also printed.
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
  real all(320)
  real als(320)
  real bel(320)
  real bes(320)
  real bexact(320)
  real fncap
  integer in
  integer incap
  integer k
  integer ncap
  real p0(320)
  real p1(320)
  real p2(320)
  real w(320)
  real x(320)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'

  ncap = 20

  do incap = 1, 4

    ncap = 2 * ncap
    fncap = real ( ncap, kind = 4 )
    fncm1 = real ( ncap - 1, kind = 4 )
!
!  Generate the abscissae and weights of the discrete inner product and
!  the exact beta-coefficients.
!
    x(1) = - 1.0E+00
    w(1) = 2.0E+00 / real ( ncap, kind = 4 )
    bexact(1) = 2.0E+00

    do k = 2, ncap
      fkm1 = real ( k - 1, kind = 4 )
      x(k) = - 1.0E+00 + 2.0E+00 * fkm1 / fncm1
      w(k) = w(1)
      bexact(k) = ( 1.0E+00 + 1.0E+00 / fncm1 )**2 &
        * ( 1.0E+00 - ( fkm1 / fncap )**2 ) / ( 4.0E+00 - ( 1.0E+00 / fkm1 )**2 )
    end do
!
!  Compute the desired coefficients, first by the Stieltjes procedure,
!  and then by the Lanczos algorithm. Indicate via the error flag  ierrs
!  whether a critical underflow condition has arisen in Stieltjes's
!  procedure. (There may, in addition, occur harmless underflow, which
!  the routine STI_R4 does not test for.)
!
    call sti_r4 ( ncap, ncap, x, w, als, bes, ierrs, p0, p1, p2 )

    call lancz_r4 ( ncap, ncap, x, w, all, bel, ierrl, p0, p1 )

    write(*,1) ierrs,ierrl
    1   format(/5x,'ierr in STI_R4 = ',i4,11x,'ierr in lancz = ',i3/)
!
!  Compute and print the absolute errors of the alpha-coefficients and
!  the relative errors of the beta-coefficients as well as the maximum
!  respective errors.
!
    erralm = 0.0E+00
    errblm = 0.0E+00
    write(*,2)
    2   format(5x,'k',4x,'erra',8x,'errb',10x,'erra',8x,'errb'/)

    do in = 1, ncap

      erras = abs ( als(in))
      errbs = abs ( (bes(in)-bexact(in))/bexact(in))
      erral = abs ( all(in))
      errbl = abs ( (bel(in)-bexact(in))/bexact(in))
      erralm = max ( erralm, erral )
      errblm = max ( errblm, errbl )

      if ( ierrs == 0 .or. in - 1 < abs ( ierrs ) ) then
        if ( in == 1 ) then
          write(*,3) in - 1,erras,errbs,erral,errbl,ncap
    3     format(1x,i5,2e12.4,2x,2e12.4,'   N  = ',i4)
        else
          write(*,4) in - 1,erras,errbs,erral,errbl
    4     format(1x,i5,2e12.4,2x,2e12.4)
        end if
      else
        if ( in == 1 ) then
          write(*,5) in - 1,erral,errbl,ncap
    5     format(1x,i5,26x,2e12.4,'   N  = ',i4)
        else
          write(*,6) in - 1,erral,errbl
    6     format(1x,i5,26x,2e12.4)
        end if
      end if
    end do

    write(*,7) erralm,errblm
    7   format(/32x,2e12.4//)

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 ...
!
!  Discussion:
!
!    This is a test of the routine MCDIS_R4, which is applied to generate
!    the first 80 recurrence coefficients for the orthogonal polynomials
!    belonging to the weight function
!
!      (1-x**2)**(-1/2) + c   on (-1,1),  c = 1, 10, 100.
!
!    The corresponding inner product is discretized by applying the
!    Gauss-Chebyshev quadrature rule to the first term of the weight
!    function, and the Gauss-Legendre rule to the second term. In
!    addition to the beta-coefficients (all alpha's are zero), the
!    routine prints the variables  ncap  and  kount  to confirm
!    convergence after one iteration.
!
  real a(81)
  real alpha(80)
  real b(81)
  real beta(80)
  real betap(80,3)
  real c
  real e(81)
  real endl(1)
  real endr(1)
  real epsma
  logical finl
  logical finr
  integer ic
  integer idelta
  integer iem
  integer iq
  integer irout
  integer mc
  integer mp
  integer n
  integer ncapm
  external qchle_r4
  real wfer(1)
  real xfer(1)
  real xp(1)
  real yp(1)

  common / common04 / c, a, b, e, epsma

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'

  epsma = epsilon ( epsma )
  iq = 1
  idelta = 2
  irout = 1
  n = 80
  mc = 2
  mp = 0
  ncapm = 81
  eps = 5000.0E+00 * epsma
  iem = 0
  c = 0.1E+00

  do ic = 1, 3

    c = 10.0E+00 * c
!
!  Compute the desired recursion coefficients. On machines with limited
!  exponent range, harmless underflow may occur in the routine  gauss_r4
!  used in  qchle  to generate the Gauss-Legendre quadrature rule.
!
    call mcdis_r4 ( n, ncapm, mc, mp, xp, yp, qchle_r4, eps, iq, idelta, irout, finl, &
      finr, endl, endr, xfer, wfer, alpha, beta, ncap, kount, ierr, ie )

    write(*,2) ncap,kount,ierr,ie,c
    2   format(3x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3, &
      ' ie = ',i3,' for c = ',f5.1)

    iem = max ( iem, abs ( ie ) )

    if ( ie /= 0 .and. abs ( ie ) <=  n ) then

      call mcdis_r4 ( abs(ie)-1, ncapm, mc, mp, xp, yp, qchle_r4, eps, iq, idelta, &
        irout, finl, finr, endl, endr, xfer, wfer, alpha, beta, ncap, kount, &
        ierr, ie )

      write(*,2) ncap,kount,ierr,ie,c
      write ( *, '(a)' ) ' '

    end if
!
!  Assemble the results in an array.
!
    do k = 1, n
      if ( ie == 0 .or. k - 1 < abs ( ie ) ) then
        betap(k,ic) = beta(k)
      else
        betap(k,ic) = 0.0E+00
      end if
    end do

  end do
!
!  Print the results.
!
  write(*,3)
    3 format(//3x,'k',4x,'beta(k), c = 1',6x,'beta(k), c = 10    beta(k), c = 100'/)

  do k = 1, n

    if ( iem == 0 .or. k - 1 < iem ) then
      write(*,4) k-1,betap(k,1), betap(k,2),betap(k,3)
    end if

    4   format(1x,i3,3e18.10)
  end do

  return
end
subroutine qchle_r4 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QCHLE_R4 returns Gauss-Chebyshev or Gauss-Legendre rules.
!
  integer n

  real a(81)
  real b(81)
  real c
  real e(81)
  real epsma
  integer i
  integer ierr
  integer k
  real pi
  real w(n)
  real x(n)

  common / common04 / c, a, b, e, epsma

  pi = 4.0E+00 * atan ( 1.0E+00 )
!
!  Gauss-Chebyshev rule.
!
  if ( i == 1 ) then

    do k = 1, n
      x(k) = cos ( real ( 2 * k - 1, kind = 4 ) * pi / real ( 2 * n, kind = 4 ) )
      w(k) = pi / real ( k, kind = 4 )
    end do
!
!  Gauss-Legendre rule.
!
  else if ( i == 2 ) then

    call recur_r4 ( n, 1, 0.0E+00, 0.0E+00, a, b, ierr )

    call gauss_r4 ( n, a, b, epsma, x, w, ierr, e )

    w(1:n) = c * w(1:n)

  end if

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 ???
!
!  Discussion:
!
!    This test generates the first 40 recursion coefficients for
!    polynomials orthogonal with respect to the Jacobi weight function
!    with parameters  alj = -.8(.2)1., bej = -.8(.2)1.  and an added mass
!    point of strength  y = .5, 1, 2, 4, 8  at the left end point. It
!    also computes the maximum relative errors (absolute errors for alpha-
!    coefficients near zero) of the computed coefficients by comparing
!    them against the exact coefficients known analytically.
!
  real a(41)
  real alpha(40)
  real b(41)
  real beta(40)
  real e(41)
  real endl(1)
  real endr(1)
  real epsma
  logical finl
  logical finr
  integer ia
  integer idelta
  integer iq
  integer irout
  integer iy
  integer mc
  integer mp
  integer n
  integer ncapm
  external qjac_r4
  real wfer(1)
  real xfer(1)
  real xp(1)
  real y
  real yp(1)

  real ( kind = 8 ) dy,dalj,dbej,dalpbe,da(41),db(41),daex(40)
  real ( kind = 8 ) dbex(40),dnum,dd,dkm1,dden

  common / common05 / a, b, e, epsma

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'

  epsma = epsilon ( epsma )
  iq = 1
  idelta = 2
  irout = 1
  n = 40
  mc = 1
  mp = 1
  xp(1) = -1.0E+00
  ncapm = 41
  eps = 5000.0E+00 * epsma
  y = 0.25E+00

  do iy = 1, 5

    y = 2.0E+00 * y
    dy = real ( y, kind = 8 )
    yp(1) = y
    write(*,2) y
    2   format(/1x,'y = ',f6.2/)
    write(*,3)
    3   format(2x,'alj',2x,'bej',5x,'erra',7x,'errb',6x,'alpha', &
      4x,'beta',4x,'ka',2x,'kb',1x,'ierr',1x,'ie',1x,'it'/)

    do ia = 1, 10

      alj = -1.0E+00 + 0.2E+00 * real ( ia, kind = 4 )
      dalj = real ( alj, kind = 8 )

      do ib = 1, 10

        bej = -1.0E+00 + 0.2E+00 * real ( ib, kind = 4 )
        alpbe = alj + bej
        dbej = real ( bej, kind = 8 )
        dalpbe = dalj + dbej
!
!  Generate the Jacobi recurrence coefficients.
!
        call recur_r4 ( ncapm, 6, alj, bej, a, b, ierr )

        call recur_r8 ( ncapm, 6, dalj, dbej, da, db, iderr )
!
!  Compute the desired recursion coefficients.
!
        call mcdis_r4 ( n, ncapm, mc, mp, xp, yp, qjac_r4, eps, iq, idelta, irout, &
          finl, finr, endl, endr, xfer, wfer, alpha, beta, ncap, kount, &
          ierr, ie )
!
!  Compute the exact coefficients by Eqs. (4.19)-(4.21) of the companion
!  paper along with the relative errors (absolute errors for alpha-
!  coefficients close to zero).
!
        daex(1) = ( da(1) - dy ) / ( 1.0D+00 + dy )
        dbex(1) = 1.0D+00 + dy
        erra = abs ( sngl ( real ( alpha(1), kind = 8 ) - daex(1) ) )

        if ( abs ( sngl ( daex(1) ) ) > eps ) then
          erra = erra / abs ( sngl ( daex(1) ) )
        end if

        errb = abs ( sngl ( ( dble ( beta(1) ) - dbex(1) ) / dbex(1) ) )
        erram = erra
        alpham = alpha(1)
        errbm = errb
        betam = beta(1)
        kam = 0
        kbm = 0
        dnum = 1.0D+00 + dy
        dd = 1.0D+00

        do k = 2, n

          dkm1 = real ( k - 1, kind = 8 )
          dden = dnum

          if ( 2 < k ) then
            dd = ( dbej + dkm1 ) * ( dalpbe + dkm1 ) * dd / ( ( dalj + dkm1 &
              - 1.0D+00 ) * ( dkm1 - 1.0D+00 ) )
          end if

          dnum = ( 1.0D+00 + ( dbej + dkm1 + 1.0D+00 ) * ( dalpbe + dkm1 + 1.0D+00 ) * dy * dd / &
            ( dkm1 * ( dalj + dkm1 ) ) ) / ( 1.0D+00 + dy * dd )

          daex(k) = da(k) + 2.0D+00 * dkm1 * ( dkm1 + dalj ) * ( dnum - 1.0D+00 ) &
            / ( ( dalpbe + 2.0D+00 * dkm1 ) * ( dalpbe + 2.0D+00 * dkm1 + 1.0D+00 ) ) &
            + 2.0D+00 * ( dbej + dkm1 + 1.0D+00 ) * ( dalpbe + dkm1 + 1.0D+00 ) * ( ( 1.0D+00 / dnum ) &
            - 1.0D+00 ) / ( ( dalpbe + 2.0D+00 * dkm1 + 1.0D+00 ) * ( dalpbe + 2.0D+00 * dkm1 + 2.0D+00 ) )

          dbex(k) = dnum * db(k) / dden
          erra = abs ( sngl( real ( alpha(k), kind = 8 ) - daex(k) ) )

          if ( eps < abs ( sngl ( daex(k) ) ) ) then
            erra = erra / abs ( sngl ( daex(k) ) )
          end if

          errb = abs ( sngl ( ( real ( beta(k), kind = 8 ) - dbex(k) ) / dbex(k) ) )

          if ( erram < erra ) then
            erram = erra
            alpham = alpha(k)
            kam = k - 1
          end if

          if ( errbm < errb ) then
            errbm = errb
            betam = beta(k)
            kbm = k - 1
          end if

        end do
!
!  Print the results.
!
        write(*,4) alj,bej,erram,errbm,alpham,betam,kam,kbm,ierr, ie,kount
    4   format(1x,2f5.2,2e11.4,2f9.6,4i4,i2)

      end do

      write(*,*)' '

    end do

  end do

  return
end
subroutine qjac_r4 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QJAC_R4 returns a Gauss-Jacobi rule.
!
  integer n

  real a(41)
  real b(41)
  real e(41)
  real epsma
  integer i
  integer ierr
  real w(n)
  real x(n)

  common / common05 / a, b, e, epsma

  call gauss_r4 ( n, a, b, epsma, x, w, ierr, e )

  w(1:n) = w(1:n) / b(1)

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 ???
!
!  Discussion:
!
!    This test generates in single and double precision the first 40
!    recursion coefficients of the orthogonal polynomials belonging to
!    the logistic density function
!
!      exp(-x)/((1+exp(-x))**2)  on (-oo, oo).
!
!    It prints the double-precision beta-coefficients (the alpha's being
!    all zero) along with the absolute and relative errors of the alpha-
!    coefficients or beta-coefficients.
!
  real a(100)
  real alpha(40)
  real b(100)
  real beta(40)
  real e(100)
  real endl(1)
  real endr(1)
  real eps
  logical finl
  logical finld
  logical finr
  logical finrd
  integer idelta
  integer iq
  integer irout
  integer mc
  integer mp
  integer n
  integer ncapm
  external qlag_r4
  external qlag_r8
  real wfer(1)
  real xfer(1)
  real xp(1)
  real yp(1)

  real ( kind = 8 ) da(300),db(300),de(300),depsma,dxp(1)
  real ( kind = 8 ) dyp(1),dendl(1),dendr(1),dxfer(1),dwfer(1)
  real ( kind = 8 ) dalpha(40),dbeta(40)
  real ( kind = 8 ) deps


  common / common06 / a, b, e, epsma
  common / common06_2 / da,db,de,depsma

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'

  iq = 1
  idelta = 1
  irout = 1
  epsma = epsilon ( epsma )
  depsma = epsilon ( depsma )
  n = 40
  mc = 2
  mp = 0
  ncapm = 100
  ncapmm = mc * ncapm + mp
  ncpmd = 300
  ncpmmd = mc * ncpmd + mp
  eps = 5000.0E+00 * epsma
  deps = 1000.0E+00 * depsma
!
!  Compute the desired coefficients. On machines with limited exponent
!  range, some of the weights in the Gauss-Laguerre quadrature rule may
!  underflow.
!
  call mcdis_r4 ( n, ncapm, mc, mp, xp, yp, qlag_r4, eps, iq, idelta, irout, finl, &
    finr, endl, endr, xfer, wfer, alpha, beta, ncap, kount, ierr, ie )

  call mcdis_r8 ( n, ncpmd, mc, mp, dxp, dyp, qlag_r8, deps, iq, idelta, irout, &
    finld, finrd, dendl, dendr, dxfer, dwfer, dalpha, dbeta, ncapd, kountd, &
    ierrd, ied )

  write(*,2) ncap,kount,ierr,ie
    2 format(/1x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3,' ie = ',i5)
  write(*,3) ncapd,kountd,ierrd,ied
    3 format(1x,'ncapd =  ',i3,' kountd =  ',i2,' ierrd =  ',i3,' ied =  ',i5/)
!
! Print the results.
!
  write(*,4)
    4 format(/5x,'k',13x,'dbeta(k)',17x,'erra',8x,'errb'/)

  do k = 1, n

    erra = abs ( alpha(k) )
    errb = abs ( sngl ( ( real ( beta(k), kind = 8 ) - dbeta(k) ) / dbeta(k) ) )

    if ( ied == 0 .or. k - 1 < ied ) then

      if ( ie == 0 .or. k - 1 < ie ) then
        write(*,5) k - 1, dbeta(k), erra, errb
    5   format(1x,i5,d33.25,2e12.4)
      else
        write(*,6) k - 1, dbeta(k)
    6   format(1x,i5,d33.25)
      end if

    end if

  end do

  return
end
subroutine qlag_r4 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QLAG_R4 returns a Gauss-Laguerre rule.
!
  integer n

  real a(100)
  real b(100)
  real e(100)
  real epsma
  integer i
  integer ierr
  integer k
  real w(n)
  real x(n)

  common / common06 / a, b, e, epsma

  call recur_r4 ( n, 7, 0.0E+00, 0.0E+00, a, b, ierr )

  call gauss_r4 ( n, a, b, epsma, x, w, ierr, e )

  do k = 1, n
    w(k) = w(k) / ( ( 1.0E+00 + exp ( - x(k) ) )**2 )
  end do

  if ( i == 1 ) then
    x(1:n) = - x(1:n)
  end if

  return
end
subroutine qlag_r8 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QLAG_R8 returns a Gauss-Laguerre rule.
!
  integer n

  real ( kind = 8 ) a(300)
  real ( kind = 8 ) b(300)
  real ( kind = 8 ) e(300)
  real ( kind = 8 ) epsma
  integer i
  integer ierr
  integer k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  common / common06_2 / a, b, e, epsma

  call recur_r8 ( n, 7, 0.0D+00, 0.0D+00, a, b, ierr )

  call gauss_r8 ( n, a, b, epsma, x, w, ierr, e )

  do k = 1, n
    w(k) = w(k) / ( ( 1.0D+00 + exp ( - x(k) ) )**2 )
  end do

  if ( i == 1 ) then
    x(1:n) = - x(1:n)
  end if

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07
!
!  Discussion:
!
!    This test generates in single and real ( kind = 8 ) the first 40
!    recurrence coefficients of the orthogonal polynomials for the half-
!    range Hermite weight function
!
!      exp(-x**2)  on  (0,oo).
!
!    Printed are the double-precision values of the alpha- and beta-
!    coefficients along with the respective relative errors of the single-
!    precision values.
!
  real alpha(40)
  real beta(40)
  real endl(4)
  real endr(4)
  logical finl
  logical finld
  logical finr
  logical finrd
  integer i
  integer idelta
  integer iq
  integer irout
  integer k
  external qgp_r4
  external qgp_r8
  real wfer(100)
  real xfer(100)
  real xp(1)
  real yp(1)

  real ( kind = 8 ) depsma,dendl(4),dendr(4),dxp(1),dyp(1)
  real ( kind = 8 ) dxfer(250),dwfer(250),dalpha(40),dbeta(40)
  real ( kind = 8 ) di,deps

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'

  finl = .true.
  finr = .false.
  finld = .true.
  finrd = .false.
  epsma = epsilon ( epsma )
  depsma = epsilon ( depsma )

  iq = 2
  idelta = 1
  irout = 1
  n = 40
  mc = 4
  mcd = 4
  mp = 0
  ncapm = 100
  ncpmd = 250
!
!  Set up the partition for the discretization of the inner product.
!
  do i = 1, 4
    endl(i) = 3.0E+00 * real ( i - 1, kind = 4 )
    endr(i) = 3.0E+00 * real ( i, kind = 4 )
    dendl(i) = 3.0D+00 * real ( i - 1, kind = 8 )
    dendr(i) = 3.0D+00 * real ( i, kind = 8 )
  end do

  eps = 50.0E+00 * epsma
  deps = 1000.0E+00 * depsma
!
!  Compute the desired recursion coefficients by the multiple-component
!  discretization procedure. On the third and fourth subinterval of the
!  partition, the quadrature weights produced by QGP_R4 or QGP_R8 may
!  underflow, on the fourth subinterval even on machines with large
!  exponent range.
!
  call mcdis_r4 ( n, ncapm, mc, mp, xp, yp, qgp_r4, eps, iq, idelta, irout, finl, &
    finr, endl, endr, xfer, wfer, alpha, beta, ncap, kount, ierr, ie )

  call mcdis_r8 ( n, ncpmd, mcd, mp, dxp, dyp, qgp_r8, deps, iq, idelta, irout, &
    finld, finrd, dendl, dendr, dxfer, dwfer, dalpha, dbeta, ncapd, kountd, &
    ierrd, ied )

  write(*,1) ncap,kount,ierr,ie
    1 format(/1x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3,' ie = ',i3)
  write(*,2) ncapd,kountd,ierrd,ied
    2 format(1x,'ncapd  = ',i3,' kountd  = ',i2,' ierrd  = ',i12,' ied  = ',i3/)
!
!  Print the results.
!
  write(*,3)
    3 format(/5x,'k',13x,'dalpha(k)',25x,'dbeta(k)',19x,'erra',29x,'errb')

  do k = 1, n
    erra = abs(sngl(( real ( alpha(k), kind = 8 ) -dalpha(k))/dalpha(k)))
    errb = abs(sngl(( real ( beta(k), kind = 8 ) -dbeta(k))/dbeta(k)))
    write(*,5) k - 1,dalpha(k),dbeta(k),erra,errb
    5   format(1x,i5,2d33.25,6x,e12.4,21x,e12.4)
  end do

  return
end
function wf_r4 ( x, i )

!*****************************************************************************80
!
!! WF_R4
!
  integer i
  real wf_r4
  real x

  wf_r4 = exp ( -x**2 )

  return
end
function wf_r8 ( x, i )

!*****************************************************************************80
!
!! WF_R8
!
  integer i
  real ( kind = 8 ) wf_r8
  real ( kind = 8 ) x

  wf_r8 = exp ( - x**2 )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08
!
!  Discussion:
!
!    This test reproduces the results of  test1  in single precision,
!    for n = 40, using the routine MCCHEB_R4 in place of CHEB_R4 and
!    Gauss-Chebyshev quadrature to discretize the modified moments.
!
  real a(79)
  real alpha(40)
  real b(79)
  real be(40)
  real beta(40)
  real endl(1)
  real endr(1)
  real eps
  real epsma
  logical finl
  logical finr
  real fnu(80)
  integer idelta
  integer iq
  integer k
  integer mc
  integer mp
  integer n
  integer ncapm
  real, dimension ( 7 ) :: oom2 = (/ &
    0.1E+00, 0.3E+00, 0.5E+00, 0.7E+00, 0.9E+00, &
    0.99E+00, 0.999E+00 /)
  external qcheb_r4
  real wfer(1)
  real xfer(1)
  real xp(1)
  real yp(1)

  common / common08 / om2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'

  epsma = epsilon ( epsma )
  iq = 1
  idelta = 1
  n = 40
  mc = 1
  mp = 0
  ncapm = 500
  eps = 100.0E+00 * epsma
!
!  Generate the recurrence coefficients for the (Chebyshev) polynomials
!  defining the modified moments.
!
  call recur_r4 ( 2 * n - 1, 3, 0.0E+00, 0.0E+00, a, b, ierr )
!
!  Compute the desired recursion coefficients by the discretized
!  Chebyshev algorithm.
!
  do iom = 1, 7

    om2 = oom2(iom)

    call mccheb_r4 ( n, ncapm, mc, mp, xp, yp, qcheb_r4, eps, iq, idelta, finl, &
      finr, endl, endr, xfer, wfer, a, b, fnu, alpha, beta, ncap, kount, &
      ierr )
!
!  On machines with limited single-precision exponent range, the routine
!  cheb may have generated an underflow exception, which however is
!  harmless and can be ignored.
!
    write(*,2) ncap,kount,ierr
    2   format(/'ncap = ',i3,'  kount = ',i3,'  ierr = ',i3/)
!
!  Print the results.
!
    write ( *, '(a)' ) '     k      beta(k)'
    write ( *, '(a)' ) ' '

    do k = 1, n
      if ( k == 1 ) then
        write(*,4) k - 1,beta(k),om2
    4   format(1x,i5,e18.10,'   om2  = ',f6.3)
      else
        write(*,5) k - 1, beta(k)
    5     format(1x,i5,e18.10)
      end if
    end do

    write ( *, '(a)' ) ' '

  end do

  return
end
subroutine qcheb_r4 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QCHEB_R4 returns a Gauss-Chebyshev rule.
!
  integer n

  integer i
  integer ierr
  integer k
  real pi
  real w(n)
  real x(n)

  common / common08 / om2

  pi = 4.0E+00 * atan ( 1.0E+00 )

  do k = 1, n
    x(k) = cos ( real ( 2 * k - 1, kind = 4 ) * pi / real ( 2 * n, kind = 4 ) )
    w(k) = pi / ( real ( n, kind = 4 ) * sqrt ( 1.0E+00 - om2 * x(k)**2 ) )
  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 ???
!
!  Discussion:
!
!    This test recomputes the results of  test2  for  sigma = .5  by applying
!    the routine  chri  with  iopt = 1, x = 0  to the weight function with
!    parameter  sigma = -.5. Printed are the relative discrepancies in both
!    single and real ( kind = 8 ) between these results and those obtained
!    in  test 2  by the modified Chebyshev algorithm. The test is embedded
!    in the routine  test2, from which all print statements have been
!    removed.
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
  real a(199)
  real alpha(100)
  real alphc(100)
  real b(199)
  real beta(100)
  real betc(100)
  real ( kind = 8 ) da(199)
  real ( kind = 8 ) dalpha(100)
  real ( kind = 8 ) dalphc(100)
  real ( kind = 8 ) db(199)
  real ( kind = 8 ) dbeta(100)
  real ( kind = 8 ) dbetc(100)
  real ( kind = 8 ) dnu(200)
  real ( kind = 8 ) ds(100)
  real ( kind = 8 ) dsigma
  real fnu(200)
  integer ierr
  logical intexp
  integer is
  logical modmom
  integer n
  integer ncd
  real s(100)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'

  modmom = .true.
!
!  Generate the recursion coefficients for the polynomials defining the
!  modified or ordinary moments.
!
  if ( modmom ) then
    n = 100
    call recur_r4 ( 2 * n - 1, 2, 0.0E+00, 0.0E+00, a, b, ierr )
    call recur_r8 ( 2 * n - 1, 2, 0.0D+00, 0.0D+00, da, db, iderr )
  else
    n = 12
    a(1:2*n-1) = 0.0E+00
    b(1:2*n-1) = 0.0E+00
    da(1:2*n-1) = 0.0D+00
    db(1:2*n-1) = 0.0D+00
  end if

  do is = 1, 3

    sigma = -0.5E+00 + 0.5E+00 * real ( is - 1, kind = 4 )
    dsigma = -0.5D+00 + 0.5D+00 * real ( is - 1, kind = 8 )

    if ( is == 2 ) then
      intexp = .true.
    else
      intexp = .false.
    end if
!
!  Compute the modified or ordinary moments using Eqs. (3.12) and
!  (3.11) of the companion paper. On machines with limited exponent
!  range, some of the high-order modified moments may underflow, without
!  this having any deteriorating effect on the accuracy.
!
    call mm_09_r4 ( n, modmom, intexp, sigma, fnu )
    call mm_09_r8 ( n, modmom, intexp, dsigma, dnu )
!
!  Compute the desired recursion coefficients by means of the modified
!  Chebyshev algorithm.
!
    call cheb_r4 ( n, a, b, fnu, alpha, beta, s, ierr )
!
!  On machines with limited single-precision exponent range,
!  CHEB_R4 may generate an underflow exception, which however is harmless
!  and can be ignored.
!
    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, iderr )
!
!  Up to this point the code is identical with the one of test2.
!
    if ( is == 1 ) then

      write(*,1) ierr,iderr
    1     format(/' ierr in CHEB_R4 = ',i4,' iderr in CHEB_R8 = ',i4/)

      if ( ierr /= 0 ) then
        nc = abs ( ierr )
      else
        nc = n
      end if

      if ( iderr /= 0 ) then
        ncd = abs ( iderr )
      else
        ncd = n
      end if
!
!  Compute the desired recursion coefficients by a modification
!  algorithm.
!
      call chri_r4 ( nc - 1, 1, alpha, beta, 0.0, 0.0, 0.0, 0.0, alphc, betc, ierr )

      call chri_r8 ( ncd - 1, 1, dalpha, dbeta, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, dalphc, &
        dbetc, iderr )

    end if

    if ( is == 3 ) then

      write(*,2)
    2     format(/1x,'test of the results for sigma = 1/2'/)
      np = nc

      if ( ncd < nc ) then
        np = ncd
      end if
!
!  Compute and print the relative discrepancies between the results of
!  the modified Chebyshev algorithm and the modification algorithm.
!
      write(*,3)
    3     format(3x,'k',2x,'err alpha',4x,'err beta',12x,'err dalpha', &
    2x,'err dbeta'/)

      do k = 1, np - 1
        errac = abs ( alpha(k) - alphc(k) ) / alpha(k)
        errbc = abs ( beta(k)-betc(k))/beta(k)
        errdac = sngl (abs(dalpha(k)-dalphc(k))/dalpha(k))
        errdbc = sngl (abs(dbeta(k)-dbetc(k))/dbetc(k))
        write(*,4) k - 1,errac,errbc,errdac,errdbc
    4   format(1x,i3,2e12.4,9x,2e12.4)
      end do

      write(*,5)
    5     format(/1x,'end of test'/)
    end if
!
!  The rest of the code is essentially the same as the corresponding
!  piece of code in  test2  with all print statements removed.
!
    eamax = 0.0E+00
    ebmax = 0.0E+00

    do k = 1, n

      erra = sngl ( abs ( real ( alpha(k), kind = 8 ) - dalpha(k) ) / dalpha(k) )
      errb = sngl ( abs ( real ( beta(k), kind = 8 ) - dbeta(k) ) / dbeta(k) )

      if ( eamax < erra ) then
        eamax = erra
        kamax = k - 1
      end if

      if ( ebmax < errb ) then
        ebmax = errb
        kbmax = k - 1
      end if

    end do

  end do

  return
end
subroutine mm_09_r4 ( n, modmom, intexp, sigma, fnu )

!*****************************************************************************80
!
!! MM_09_R4 generates ordinary or modified moments for a particular weight function.
!
!  Discussion:
!
!    The routine generates the first  2*n  modified moments (if modmom = .true.)
!    relative to shifted monic Legendre polynomials, using Eq. (3.12) of
!    the companion paper, and the first  2*n  ordinary moments (if modmom
!     = .false.) by Eq. (3.11), of the weight function
!
!      (x**sigma)*ln(1/x)  on (0,1],   sigma > -1,
!
!    for sigma an integer (if intexp = .true.) or a real number (if intexp
!     = .false.). In either case, the input variable sigma is of type real.
!
  integer n

  real fnu(2*n)
  logical intexp
  integer iq
  integer isigma
  integer k
  integer kmax
  logical modmom
  integer nd
  real sigma
  real sigp1

  nd = 2 * n
  sigp1 = sigma + 1.0E+00

  if ( modmom ) then

    isigma = int ( sigma )

    if ( intexp .and. isigma + 1 < 2 * n ) then
      kmax = isigma + 1
    else
      kmax = 2 * n
    end if

    c = 1.0E+00

    do k = 1, kmax
      fk = real ( k, kind = 4 )
      p = 1.0E+00
      s = 1.0E+00 / sigp1
      if ( 1 < kmax ) then
        do i = 1, k - 1
          fi = real ( i, kind = 4 )
          p = ( sigp1 - fi ) * p / ( sigp1 + fi )
          s = s + 1.0E+00 / ( sigp1 + fi ) - 1.0E+00 / ( sigp1 - fi )
        end do
      end if
      fnu(k) = c * s * p / sigp1
      c = fk * c / ( 4.0E+00 * fk - 2.0E+00 )
    end do

    if ( .not. intexp .or. 2 * n <= isigma + 1 ) then
      return
    end if

    q = - 0.5E+00

    do iq = 1, isigma
      fiq = real ( iq, kind = 4 )
      q = fiq * fiq * q / ( ( 2.0E+00 * fiq + 1.0E+00 ) * ( 2.0E+00 * fiq + 2.0E+00 ) )
    end do

    fnu(isigma+2) = c * q

    do k = isigma + 3, 2 * n
      fkm1 = real ( k - 1, kind = 4 )
      fnu(k) = - fkm1 * ( fkm1 - sigp1 ) * fnu( k - 1 ) &
        / ( ( 4.0E+00 * fkm1 - 2.0E+00 ) * ( fkm1 + sigp1 ) )
    end do

  else

    do k = 1, 2 * n
      fkm1 = real ( k - 1, kind = 4 )
      fnu(k) = ( 1.0E+00 / ( sigp1 + fkm1 ) )**2
    end do

  end if

  return
end
subroutine mm_09_r8 ( n, modmom, intexp, dsigma, dnu )

!*****************************************************************************80
!
!! MM_09_R8 is a double-precision version of the routine MM_09_R4.
!
  integer n

  real ( kind = 8 ) dc
  real ( kind = 8 ) di
  real ( kind = 8 ) diq
  real ( kind = 8 ) dk
  real ( kind = 8 ) dkm1
  real ( kind = 8 ) dnu(2*n)
  real ( kind = 8 ) dp
  real ( kind = 8 ) dq
  real ( kind = 8 ) ds
  real ( kind = 8 ) dsigma
  real ( kind = 8 ) dsigp1
  integer i
  logical intexp
  integer iq
  integer isigma
  integer k
  integer kmax
  logical modmom
  integer nd

  nd = 2 * n
  dsigp1 = dsigma + 1.0D+00

  if ( modmom ) then

    isigma = int ( dsigma )

    if ( intexp .and. isigma + 1 < 2 * n ) then
      kmax = isigma + 1
    else
      kmax = 2 * n
    end if

    dc = 1.0D+00

    do k = 1, kmax

      dk = real ( k, kind = 8 )
      dp = 1.0D+00
      ds = 1.0D+00 / dsigp1

      if ( 1 < kmax ) then
        do i = 1, k - 1
          di = real ( i, kind = 8 )
          dp = ( dsigp1 - di ) * dp / ( dsigp1 + di )
          ds = ds + 1.0D+00 / ( dsigp1 + di ) - 1.0D+00 / ( dsigp1 - di )
        end do
      end if

      dnu(k) = dc * ds * dp / dsigp1
      dc = dk * dc / ( 4.0D+00 * dk - 2.0D+00 )

    end do

    if ( .not. intexp .or. 2 * n <= isigma + 1 ) then
      return
    end if

    dq = - 0.5D+00

    do iq = 1, isigma
      diq = real ( iq, kind = 8 )
      dq = diq * diq * dq / ( ( 2.0D+00 * diq + 1.0D+00 ) * ( 2.0D+00 * diq + 2.0D+00 ) )
    end do

    dnu(isigma+2) = dc * dq

    do k = isigma + 3, 2 * n
      dkm1 = real ( k - 1, kind = 8 )
      dnu(k) = - dkm1 * ( dkm1 - dsigp1 ) * dnu( k - 1 ) &
        / ( ( 4.0D+00 * dkm1 - 2.0D+00 ) * ( dkm1 + dsigp1 ) )
    end do

  else

    do k = 1, 2 * n
      dkm1 = real ( k - 1, kind = 8 )
      dnu(k) = ( 1.0D+00 / ( dsigp1 + dkm1 ) )**2
    end do

  end if

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10
!
!  Discussion:
!
!    This test applies the routines INDP_R4 and INDP_R8 to generate the
!    first 20 recursion coefficients of the induced Legendre polynomials
!    pind(k,m)(.), m = 0,1,2,...,11, that is, of the polynomials orthogonal
!    relative to the weight function
!
!      [p(m)(x)]**2    on [-1,1],
!
!    where  p(m)(.)  is the (monic) Legendre polynomial of degree m.
!    (When m = 0, then  pind(k,0)(.) = p(k)(.).) The routine also prints the
!    absolute and relative errors, respectively, of the alpha- and beta-
!    coefficients.
!
  real a(31)
  real a1(31)
  real alpha(31)
  real b(31)
  real b1(31)
  real beta(31)
  real betap(20,12)
  real ( kind = 8 ) da(31)
  real ( kind = 8 ) da1(31)
  real ( kind = 8 ) dalpha(31)
  real ( kind = 8 ) db(31)
  real ( kind = 8 ) db1(31)
  real ( kind = 8 ) dbeta(31)
  real ( kind = 8 ) de(11)
  real ( kind = 8 ) depsma
  real ( kind = 8 ) dw(11)
  real ( kind = 8 ) dz(11)
  real e(11)
  real epsma
  real erram(12)
  real errbm(12)
  integer im
  integer ip
  integer m
  integer n
  real w(11)
  real z(11)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'

  epsma = epsilon ( epsma )
  depsma = epsilon ( depsma )
  n = 20

  do im = 1, 12

    m = im - 1
!
!  Generate the Legendre recurrence coefficients required in the
!  routines indp and dindp.
!
    call recur_r4 ( n + m, 1, 0.0E+00, 0.0E+00, a, b, ierr )
    call recur_r8 ( n + m, 1, 0.0D+00, 0.0D+00, da, db, ierr )
!
!  Compute the desired recursion coefficients.
!
    call indp_r4 ( n, m, a, b, epsma, alpha, beta, ierr, z, w, e, a1, b1 )

    call indp_r8 ( n, m, da, db, depsma, dalpha, dbeta, ierr, dz, dw, &
      de, da1, db1 )
!
!  Compute and print the respective errors.
!
    erram(im) = 0.0E+00
    errbm(im) = 0.0E+00

    do k = 1, n
      erra = sngl(abs( real ( alpha(k), kind = 8 )-dalpha(k)))
      errb = sngl(abs(( real ( beta(k), kind = 8 )-dbeta(k))/dbeta(k)))
      erram(im) = max ( erram(im), erra )
      errbm(im) = max ( errbm(im), errb )
      betap(k,im) = real ( dbeta(k), kind = 4 )
    end do

  end do

  do ip = 1, 3

    ip4 = 1 + 4 * ( ip - 1 )

    write(*,1) ip4-1,ip4,ip4+1,ip4+2
    1   format(5x,'k',2x,'m = ',i1,'  beta(k)',2x,'m = ',i1,'  beta(k)',2x, &
      'm = ',i2,' beta(k)',2x,'m = ',i2,' beta(k)'/)

    do k = 1, n
      write(*,2) k - 1,betap(k,ip4),betap(k,ip4+1),betap(k,ip4+2), &
        betap(k,ip4+3)
    2     format(1x,i5,4f14.10)
    end do

    write(*,3) erram(ip4),erram(ip4+1),erram(ip4+2),erram(ip4+3)
    3   format(/4x,'erra',e12.4,3e14.4)
    write(*,4) errbm(ip4),errbm(ip4+1),errbm(ip4+2),errbm(ip4+3)
    4   format(4x,'errb',e12.4,3e14.4//)

  end do

  return
end
subroutine indp_r4 ( n, m, a, b, eps, alpha, beta, ierr, z, w, e, a1, b1 )

!*****************************************************************************80
!
!! INDP_R4
!
!  Discussion:
!
!    If  p(m)(.)  denotes the (monic) orthogonal polynomial of degree  m
!    relative to the weight function  w(x), then the corresponding m-th
!    induced orthogonal polynomials  pind(k,m)(.), k = 0,1,2,..., are those
!    orthogonal with respect to the weight function
!
!      (p(m)(x)**2)*w(x).
!
!    For background on induced orthogonal polynomials, including an
!    algorithm for generating their recursion coefficients, see W. Gautschi
!    and S. Li,A set of orthogonal polynomials induced by a given
!    orthogonal polynomial'', Aequationes Math., to appear.
!
!    This routine obtains the first n recurrence coefficients of the m-th
!    induced orthogonal polynomials by an m-fold application of the routine
!    CHRI_R4 with  iopt = 7, the shifts taken being, in succession, the zeros of
!    p(m)(.).
!
  integer m
  integer n

  real a(n+m)
  real a1(n+m)
  real alpha(n+m)
  real b(n+m)
  real b1(n+m)
  real beta(n+m)
  integer imu
  integer mi
  real w(m)
  real x
  real z(m)

  alpha(1:n+m) = a(1:n+m)
  beta(1:n+m) = b(1:n+m)

  if ( m == 0 ) then
    return
  end if

  call gauss_r4 ( m, a, b, eps, z, w, ierr )

  do imu = 1, m
    mi = n + m - imu
    a1(1:mi+1) = alpha(1:mi+1)
    b1(1:mi+1) = beta(1:mi+1)
    x = z(imu)
    call chri_r4 ( mi, 7, a1, b1, x, 0.0, 0.0, 0.0, alpha, beta, ierrc )
  end do

  return
end
subroutine indp_r8 ( n, m, da, db, deps, dalpha, dbeta, ierr, dz, dw, &
  de, da1, db1 )

!*****************************************************************************80
!
!! INDP_R8 is a double-precision version of the routine  indp.
!
  integer m
  integer n

  real ( kind = 8 ) da(n+m)
  real ( kind = 8 ) da1(n+m)
  real ( kind = 8 ) dalpha(n+m)
  real ( kind = 8 ) db(n+m)
  real ( kind = 8 ) db1(n+m)
  real ( kind = 8 ) dbeta(n+m)
  real ( kind = 8 ) deps
  real ( kind = 8 ) dw(m)
  real ( kind = 8 ) dx
  real ( kind = 8 ) dz(m)
  integer imu
  integer mi

  dalpha(1:n+m) = da(1:n+m)
  dbeta(1:n+m) = db(1:n+m)

  if ( m == 0 ) then
    return
  end if

  call gauss_r8 ( m, da, db, deps, dz, dw, ierr )

  do imu = 1, m
    mi = n + m - imu
    da1(1:mi+1) = dalpha(1:mi+1)
    db1(1:mi+1) = dbeta(1:mi+1)
    dx = dz(imu)
    call chri_r8 ( mi, 7, da1, db1, dx, 0.0D+00, 0.0D+00, 0.0D+00, dalpha, &
      dbeta, ierrc )
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11
!
!  Discussion:
!
!    This test is to illustrate the dissimilar performance of the routines
!    CHRI_R4 and GCHRI_R4 in the case of division of the Jacobi weight
!    function  w(t;alj,bej)  with parameters  alj,bej  by either a linear
!    divisor  t-x  or a quadratic divisor  (t-x)**2 + y**2.
!
!    In either case, the parameters selected are  alj = -.8(.4).8,
!    bej = alj(.4).8.  In the former case, x = -1.001, -1.01, -1.04, -1.07
!    and -1.1, whereas in the latter case,  x and  y  are chosen to lie,
!    regularly spaced, on the upper half of an ellipse with foci at +1 and
!    -1 and sum of the semiaxes equal to  rho = 1.05, 1.1625, 1.275, 1.3875
!    and 1.5.
!
!    The routines are run in both single and real ( kind = 8 ) with n = 40, the
!    results of the latter being used to calculate, and print, the maximum
!    absolute and relative error of the single-precision alpha- and beta-
!    coefficients, respectively.
!
!    Also printed are the starting recurrence indexes required in the backward
!    recurrence schemes of GCHRI_R4 and GCHRI_R8 to achieve single and
!    double-precision accuracy.
!
!    This information is contained in the first line of each 3-line block
!    of the output,
!    where in the case of quadratic divisors only average values (averaged
!    over the upper half of the respective ellipse) are shown. The second
!    and third line of each 3-line block display the maximum
!    reconstruction error'', that is, the maximum errors in the alpha's
!    and beta's if the coefficients produced by  gchri,chri  and  dgchri,
!    dchri  are fed back to the routines  chri  and  dchri  with  iopt = 1
!    to recover the original recursion coefficients in single and double
!    precision.
!
  real a(500)
  real alpha(40)
  real alphc(40)
  real alphcr(40)
  real alphr(40)
  real b(500)
  real beta(40)
  real betc(40)
  real betcr(40)
  real betr(40)
  real ( kind = 8 ) da(800)
  real ( kind = 8 ) dal
  real ( kind = 8 ) db(800)
  real ( kind = 8 ) dbe
  real ( kind = 8 ) deps
  real ( kind = 8 ) depsma
  real ( kind = 8 ) dhi
  complex e
  integer ial
  integer ipoly
  integer n
  integer nd
  integer nu0
  integer nu0jac
  integer numax
  integer numaxd
  complex rho(80)
  complex rold(80)
  real, dimension ( 5 ) :: rr = (/ &
    1.05E+00, 1.1625E+00, 1.275E+00, 1.3875E+00, 1.5E+00 /)
  real, dimension ( 5 ) :: xx = (/ &
    1.001E+00, 1.01E+00, 1.04E+00, 1.07E+00, 1.1E+00 /)
  complex z

  real ( kind = 8 ) dx,dy,dalpha(40),dbeta(40),drhor(80),drhoi(80)
  real ( kind = 8 ) droldr(80),droldi(80),dhr
  real ( kind = 8 ) dalphc(40),dbetc(40),dalphr(40),dbetr(40),dalcr(40),dbecr(40)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'

  write(*,1)
    1 format(/)
  epsma = epsilon ( epsma )
  depsma = epsilon ( depsma )
  n = 40
  nd = 2 * n
  numax = 500
  numaxd = 800
  eps = 10.0E+00 * epsma
  deps = 100.0D+00 * depsma
  epsd = sngl ( deps )
  ipoly = 6

  do ial = 1, 5

    al = - 0.8E+00 + 0.4E+00 * real ( ial - 1 )
    dal =  real ( al, kind = 8 )
    ibemax = 6 - ial

    do ibe = 1, ibemax

      be = al + 0.4E+00 * real ( ibe - 1 )
      dbe = real ( be, kind = 8 )
      write(*,2) al,be
    2     format(///1x,'al = ',f6.2,'  be = ',f6.2//)
      hi = 0.0E+00
      dhi = 0.0D+00
!
!  Generate the Jacobi recurrence coefficients to be used in the
!  backward recurrence algorithm of the routines GCHRI_R4 and GCHRI_R8.
!
      call recur_r4 ( numax, ipoly, al, be, a, b, ierr )

      call recur_r8 ( numaxd, ipoly, dal, dbe, da, db, ierrd )

      write(*,3)
    3     format(30x,'gchri',20x,'chri')
      write(*,4)
    4     format(5x,'x',5x,'nu0',2x,'nud0',4x,'erra',8x,'errb',10x, &
    'erra',7x,'errb'/)

      do ix = 1, 5

        x = - xx(ix)
        dx =  real ( x, kind = 8 )
        y = 0.0E+00
        dy = 0.0D+00
        z = cmplx ( x, y )
!
!  Compute the starting index for backward recurrence.
!
        nu0 = nu0jac ( 2 * n - 1, z, eps )
        nu0d = nu0jac ( 2 * n - 1, z, epsd )
!
!  Generate the recurrence coefficients for the Jacobi weight function
!  divided by a linear divisor.
!
        call gchri_r4 ( n, 1, nu0, numax, eps, a, b, x, y, alpha, beta, nu, ierrg, &
          ierrc, rho, rold )
!
!  On machines with limited single-precision exponent range, the routine
!  CHEB_R4 used in GCHRI_R4 may have generated an underflow exception,
!  which however is harmless and can be ignored.
!
        call gchri_r8 ( n, 1, nu0d, numaxd, deps, da, db, dx, dy, dalpha, dbeta, &
          nud, ierrgd, ierrcd, drhor, drhoi, droldr, droldi )

        if ( ierrg /= 0 .or. ierrc /= 0 .or. ierrgd /= 0 .or. ierrcd /= 0 ) then
          write(*,5) ierrg,ierrgd,al,be,x
    5     format(/1x,'ierrg in GCHRI_R4 = ',i4,' ierrg in GCHRI_R8 = ', &
          i4,' for al = ',f6.2,' be = ',f6.2,' x = ',f7.4)
          write(*,6) ierrc,ierrcd,al,be,x
    6     format(1x,'ierrc in GCHRI_R4 = ',i4,' ierrc in GCHRI_R8 = ', &
          i4,' for al = ',f6.2,' be = ',f6.2,' x = ',f7.4/)
          return
        end if
!
!  Generate the recurrence coefficients for the Jacobi weight function
!  divided by a linear divisor, using the routines  chri, dchri.
!
        hr = real ( rho(1) )
        call chri_r4 ( n, 4, a, b, x, y, hr, hi, alphc, betc, ierr )

        dhr = drhor(1)
        call chri_r8 ( n, 4, da, db, dx, dy, dhr, dhi, dalphc, dbetc, ierr )
!
!  Do the reconstruction.
!
        call chri_r4 ( n - 1, 1, alpha, beta, x, y, 0.0E+00, 0.0E+00, alphr, betr, ierr )

        call chri_r8 ( n - 1, 1, dalpha, dbeta, dx, dy, 0.0D+00, 0.0D+00, dalphr, dbetr, ierr )

        call chri_r4 ( n - 1, 1, alphc, betc, x, y, 0.0E+00, 0.0E+00, alphcr, betcr, ierr )

        call chri_r8 ( n - 1, 1, dalphc, dbetc, dx, dy, 0.0D+00, 0.0D+00, dalcr, dbecr, ierr )
!
!  Compute and print the maximum errors.
!
        erragm = 0.0E+00
        errbgm = 0.0E+00
        erracm = 0.0E+00
        errbcm = 0.0E+00
        erram = 0.0E+00
        errbm = 0.0E+00
        errdam = 0.0E+00
        errdbm = 0.0E+00
        eracrm = 0.0E+00
        erbcrm = 0.0E+00
        edacrm = 0.0E+00
        edbcrm = 0.0E+00

        do k = 1, n

          errag = abs (sngl( real ( alpha(k), kind = 8 ) -dalpha(k)))
          errbg = abs (sngl(( real ( beta(k), kind = 8 ) -dbeta(k))/dbeta(k)))
          errac = abs (sngl( real ( alphc(k), kind = 8 ) -dalphc(k)))
          errbc = abs (sngl(( real ( betc(k), kind = 8 ) -dbetc(k))/dbetc(k)))

          if ( k < n ) then

            erra = abs(alphr(k)-a(k))
            errb = abs((betr(k)-b(k))/b(k))
            errda = abs(sngl(dalphr(k)-da(k)))
            errdb = abs(sngl((dbetr(k)-db(k))/db(k)))
            eracr = abs(alphcr(k)-a(k))
            erbcr = abs((betcr(k)-b(k))/b(k))
            edacr = abs(sngl(dalcr(k)-da(k)))
            edbcr = abs(sngl((dbecr(k)-db(k))/db(k)))

            erram = max ( erram, erra )
            errbm = max ( errbm, errb )
            errdam = max ( errdam, errda )
            errdbm = max ( errdbm, errdb )
            eracrm = max ( eracrm, eracr )
            erbcrm = max ( erbcrm, erbcr )
            edacrm = max ( edacrm, edacr )
            edbcrm = max ( edbcrm, edbcr )

          end if

          erragm = max ( erragm, errag )
          errbgm = max ( errbgm, errbg )
          erracm = max ( erracm, errac )
          errbcm = max ( errbcm, errbc )

        end do

        write(*,7) x,nu0,nu0d,erragm,errbgm,erracm,errbcm
    7   format(/1x,f7.4,2i6,2e12.4,2x,2e12.4)

        if ( ix == 1) then
          write(*,8) erram,errbm,errdam,errdbm
    8     format(11x,'reconstr.',2e12.4,2x,2e12.4)
          write(*,9) eracrm,erbcrm,edacrm,edbcrm
    9     format(12x,'errors',2x,2e12.4,2x,2e12.4)
        else
          write(*,11) erram,errbm,errdam,errdbm
   11     format(20x,2e12.4,2x,2e12.4)
          write(*,11) eracrm,erbcrm,edacrm,edbcrm
        end if

      end do

      write(*,1)
      ndiv = 20
      fndiv = real ( ndiv, kind = 4 )
      fndm1 = real ( ndiv - 1, kind = 4 )
      pi = 4.0E+00 * atan ( 1.0E+00 )
      write(*,3)
      write(*,12)
   12     format(4x,'rho',4x,'nu0',2x,'nud0',4x,'erra',8x,'errb',10x, &
       'erra',7x,'errb'/)

      do ir = 1, 5

        r = rr(ir)
        agmv = 0.0E+00
        bgmv = 0.0E+00
        acmv = 0.0E+00
        bcmv = 0.0E+00
        amv = 0.0E+00
        bmv = 0.0E+00
        damv = 0.0E+00
        dbmv = 0.0E+00
        acrmv = 0.0E+00
        bcrmv = 0.0E+00
        dacrmv = 0.0E+00
        dbcrmv = 0.0E+00
        nu0v = 0
        nu0dv = 0

        do ith = 1, ndiv - 1
!
!  Generate the points on the ellipse.
!
          theta = pi * real ( ith ) / fndiv
          e = cmplx ( cos ( theta ), sin ( theta ) )
          z = 0.5E+00 * ( r * e + 1.0E+00 / ( r * e ) )
          x = real ( z )
          y = aimag ( z )
          dx = real ( x, kind = 8 )
          dy = real ( y, kind = 8 )
!
!  Compute the starting index for backward recurrence.
!
          nu0 = nu0jac ( 2 * n - 1, z, eps )
          nu0d = nu0jac ( 2 * n - 1, z, epsd )
!
!  Generate the recurrence coefficients for the Jacobi weight function
!  divided by a quadratic divisor.
!
          call gchri_r4 ( n, 2, nu0, numax, eps, a, b, x, y, alpha, beta, &
            nu, ierrg, ierrc, rho, rold )

          call gchri_r8 ( n, 2, nu0d, numaxd, deps, da, db, dx, dy, dalpha, dbeta, &
            nud, ierrgd, ierrcd, drhor, drhoi, droldr, droldi )

          if ( ierrg /= 0 .or. ierrc /= 0 .or. ierrgd /= 0 .or. ierrcd /= 0 ) then
            write(*,5) ierrg,ierrgd,al,be,x
            write(*,6) ierrc,ierrcd,al,be,x
            cycle
          end if

          nu0v = nu0v + nu0
          nu0dv = nu0dv + nu0d
!
!  Generate the recurrence coefficients for the Jacobi weight function
!  divided by a quadratic divisor, using the routines  chri,dchri.
!
          hr = real ( rho(1), kind = 4 )
          hi = aimag ( rho(1) )

          call chri_r4 ( n, 5, a, b, x, y, hr, hi, alphc, betc, ierr )

          dhr = drhor(1)
          dhi = drhoi(1)

          call chri_r8 ( n, 5, da, db, dx, dy, dhr, dhi, dalphc, dbetc, ierr )
!
!  Do the reconstruction.
!
          call chri_r4 ( n - 1, 2, alpha, beta, x, y, 0.0, 0.0, alphr, betr, ierr )

          call chri_r8 ( n - 1, 2, dalpha, dbeta, dx, dy, 0.0D+00, 0.0D+00, dalphr, &
            dbetr, ierr )

          call chri_r4 ( n - 1, 2, alphc, betc, x, y, 0.0, 0.0, alphcr, betcr, ierr )

          call chri_r8 ( n - 1, 2, dalphc, dbetc, dx, dy, 0.0D+00, 0.0D+00, dalcr, dbecr, ierr )
!
!  Compute and print the maximum average errors.
!
          erragm = 0.0E+00
          errbgm = 0.0E+00
          erracm = 0.0E+00
          errbcm = 0.0E+00
          erram = 0.0E+00
          errbm = 0.0E+00
          errdam = 0.0E+00
          errdbm = 0.0E+00
          eracrm = 0.0E+00
          erbcrm = 0.0E+00
          edacrm = 0.0E+00
          edbcrm = 0.0E+00

          do k = 1, n

            errag = abs ( sngl ( real ( alpha(k), kind = 8 ) - dalpha(k)))
            errbg = abs ( sngl (( real ( beta(k), kind = 8 ) - dbeta(k))/dbeta(k)))
            errac = abs ( sngl ( real ( alphc(k), kind = 8 ) - dalphc(k)))
            errbc = abs ( sngl (( real ( betc(k), kind = 8 ) - dbetc(k))/dbetc(k)))

            if ( k < n ) then

              erra = abs ( alphr(k) - a(k) )
              errb = abs ( ( betr(k) - b(k) ) / b(k) )
              errda = abs ( sngl ( dalphr(k) - da(k) ) )
              errdb = abs ( sngl ( ( dbetr(k) - db(k) ) / db(k) ) )
              eracr = abs ( alphcr(k) - a(k) )
              erbcr = abs ( ( betcr(k) - b(k) ) / b(k) )
              edacr = abs ( sngl ( dalcr(k) - da(k) ) )
              edbcr = abs ( sngl ( ( dbecr(k) - db(k) ) / db(k) ) )

              erra = max ( erra, erra )
              errb = max ( errb, errb )
              errda = max ( errda, errda )
              errdb = max ( errdb, errdb )
              eracr = max ( eracr, eracr )
              erbcr = max ( erbcr, erbcr )
              edacr = max ( edacr, edacr )
              edbcr = max ( edbcr, edbcr )

            end if

            erragm = max ( erragm, errag )
            errbgm = max ( errbgm, errbg )
            erracm = max ( erracm, errac )
            errbcm = max ( errbcm, errbc )

          end do

          agmv   = agmv   + erragm
          bgmv   = bgmv   + errbgm
          acmv   = acmv   + erracm
          bcmv   = bcmv   + errbcm
          amv    = amv    + erram
          bmv    = bmv    + errbm
          damv   = damv   + errdam
          dbmv   = dbmv   + errdbm
          acrmv  = acrmv  + eracrm
          bcrmv  = bcrmv  + erbcrm
          dacrmv = dacrmv + edacrm
          dbcrmv = dbcrmv + edbcrm

        end do

        nu0 = real ( nu0v, kind = 4 ) / fndm1
        nu0d = real ( nu0dv, kind = 4 ) / fndm1
        erragm = agmv / fndm1
        errbgm = bgmv / fndm1
        erracm = acmv / fndm1
        errbcm = bcmv / fndm1
        erram = amv / fndm1
        errbm = bmv / fndm1
        errdam = damv / fndm1
        errdbm = dbmv / fndm1
        eracrm = acrmv / fndm1
        erbcrm = bcrmv / fndm1
        edacrm = dacrmv / fndm1
        edbcrm = dbcrmv / fndm1

        if ( ir == 1 ) then
          write(*,7) r,nu0,nu0d,erragm,errbgm,erracm,errbcm
          write(*,8) erram,errbm,errdam,errdbm
          write(*,9) eracrm,erbcrm,edacrm,edbcrm
        else
          write(*,7) r,nu0,nu0d,erragm,errbgm,erracm,errbcm
          write(*,11) erram,errbm,errdam,errdbm
          write(*,11) eracrm,erbcrm,edacrm,edbcrm
        end if

      end do

    end do

  end do

  return
end
subroutine test99 ( )

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) n
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  n = 63
  allocate ( x(1:n) )
  allocate ( w(1:n) )
  i = 0

  call qlag_r8 ( n, x, w, i, ierr )

  do i = 1, n
    write ( *, * ) i, x(i), w(i)
  end do

  deallocate ( w )
  deallocate ( x )

  return
end
