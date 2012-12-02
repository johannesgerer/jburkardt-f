program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS611_PRB.
!
!  Discussion:
!
!    TOMS611_PRB is a simple test program for SMSNO, SUMSL, and HUMSL.
!
!    In these examples,  n = 4  and  f(x) = (1.0 + 0.5*(x1**t)*a*x1)**0.5,
!    where  x1(i) = d1(i)*x(i) - i,  **t denotes transpose, and  a  is a
!    matrix having fives on the main diagonal and ones everywhere else.
!    the scale vector  d1  is passed to qdrtf, the subroutine that
!    evaluates  f,  as part of urparm.  specifically, the matrix  urp
!    declared below is passed for ufparm, and  d1  is urp(*,1), the first
!    column of urp.  this main program repeatedly minimizes f, starting
!    from  x = 0,  by calling smsno, sumsl, and humsl.  we actually
!    use two different objective functions, since we change  d1  after
!    the first call on sumsl.  all runs but the last use  d = d1.
!
!    f(x) is minimized at  x1 = 0  (a vector of zeros), i.e., at
!    x(i) = i/d1(i).
!
!  Local Parameters:
!
!    qdrtf  - passed for calcf to smsno, sumsl, and humsl.
!
!    qdrtg  - passed for calcg to sumsl.
!
!    qdrtgh - passed for calcgh to humsl.
!
  implicit none

  integer, parameter :: liv = 60
  integer, parameter :: lv = 150

  real ( kind = 8 ) d(4)
  external deflt
  external humsl
  integer i
  integer iv(liv)
  external qdrtf
  external qdrtg
  external qdrtgh
  integer uip(1)
  real ( kind = 8 ) urp(4,3)
  real ( kind = 8 ) v(lv)
  real ( kind = 8 ) x(4)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS611_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS611 library.'
!
!  initialize d, d1, and x.
!
  d(1:4) = 1.0D+00
  urp(1:4,1) = 1.0D+00
  x(1:4) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SMSNO on QDRTF:'
!
!  Before this first call, we set iv(1) to 0 so that all input
!  components of iv and v will be given default values.  before
!  subsequent calls, we set iv(1) to 12 so that the old input values
!  of iv and v are used.
!
!  qdrtf does not make use of ufparm.  in calling smsno, we
!  arbitrarily pass qdrtf for ufparm to satifsy the calling sequence.
!
  iv(1) = 0
  call smsno(4, d, x, qdrtf, iv, liv, lv, v, uip, urp, qdrtf)
!
!  we reinitialize x and minimize  f  again, this time using sumsl.
!  qdrtg, the subroutine passed for calcg, assumes that ufparm is qdrtf.
!
  x(1:4) = 0.0d+0

  write(*,40)
 40   format(/16h sumsl on qdrtf )

  iv(1) = 12
  call sumsl(4, d, x, qdrtf, qdrtg, iv, liv, lv, v, uip,urp,qdrtf)
!
!  now we modify  f  by using a different choice of d1.  we still use
!  d = d1, so the performance of sumsl should stay the same -- only d
!  and the final x and gradient should be affected.
!
  do i = 1, 4
     x(i) = 0.0d+0
     d(i) = 1.0d+02 ** i
     urp(i,1) = d(i)
  end do

  write(*,40)

  iv(1) = 12
  call sumsl(4, d, x, qdrtf, qdrtg, iv, liv, lv, v, uip,urp,qdrtf)
!
!  using the last choice of d and d1, we now use humsl to minimize f.
!  like qdrtg, qdrtgh assumes that ufparm is qdrtf.
!
  x(1:4) = 0.0D+00
  write(*,70)
 70   format(/16h humsl on qdrtf )

  iv(1) = 12
  call humsl(4, d, x, qdrtf, qdrtgh, iv, liv, lv, v, uip,urp,qdrtf)
!
!  we repeat the last run with iv(dtype) = 1 and v(dinit) = 0.0, so
!  that humsl will determine  d  from the diagonal of the hessian.
!  this run also demonstrates the use of subroutine deflt and the
!  passing of nondefault parameters.  (since the iv and v input
!  components still have their default values at his point, it is not
!  really necessary to call deflt.  it is necessary to reset iv(1) to
!  12, however, and deflt does this for us.)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUMSL updating D.'

  x(1:4) = 0.0D+00

  call deflt(2, iv, liv, lv, v)
  iv(16) = 1
  v(38) = 0.0D+00
  call humsl(4, d, x, qdrtf, qdrtgh, iv, liv, lv, v, uip,urp,qdrtf)
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS611_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine qdrtf ( n, x, nf, f, uip, urp, ufp )

!*******************************************************************************
!
!! QDRTF evaluates the objective function f(x).
!
!  Discussion:
!
!    The function f(x) is described in the
!    main program above.  it stores in urp(*,2) and urp(*,3) some
!    information useful in evaluating the gradient and hessian of f, and
!    it stores  nf  in uip(1) to identify the x corresponding to this
!    information.  f(x) has the form  f(x) = phi(q(x)),  where  q  is a
!    quadratic form and  phi(y) = y**0.5.  the gradient of  f  is
!    g(x) = phiprm(q(x))*gq(x),  where  phiprm  is the derivative of phi
!    and  gq  is the gradient of q.  this routine stores phiprm(q(x)) in
!    urp(1,3) and gq(x) in urp(*,2).  the hessian of f is
!    h(x) = phi2prm(q(x))*gq(x)*gq(x)**t + phiprm(q(x))*hq(x),  where
!    phi2prm  is the second derivative of phi, **t denotes transpose,
!    and  hq  is the hessian of q.  this routine stores phi2prm(q(x)) in
!    urp(2,3).  the subroutines qdrtg and qdrtgh given below would work
!    without change on any other choice of phi.  qdrtg would also work
!    with any other differentiable function q.  qdrtgh, on the other
!    hand, assumes that  hq(x)  is the matrix  a  described in the main
!    program above.
!
  implicit none

  integer n

  integer nf, uip(*)
  real ( kind = 8 ) x(n), f, urp(n,3)
  external ufp
  integer i
  real ( kind = 8 ) dn, f2, t, t1

  uip(1) = nf
  dn = n
  t = 0.d+0
  do i = 1, n
    urp(i,2) = urp(i,1) * x(i) - dble ( i )
    t = t + urp(i,2)
  end do

  f2 = 0.d+0

  do i = 1, n
    t1 = dn * urp(i,2) + t
    f2 = f2 + t1 * urp(i,2)
    urp(i,2) = urp(i,1) * t1
  end do

  f2 = 1.d+0  +   0.5d+0 * f2
  f = dsqrt(f2)
  urp(1,3) = 0.5D+00 / f
  urp(2,3) = -0.5D+00 / ( f * f2 )

  return
end
subroutine qdrtg ( n, x, nf, g, uip, urp, qdrtf )

!*******************************************************************************
!
!! QDRTG evaluates the gradient of the objective function f(x).
!
!  Discussion:
!
!    The objective function is described in the main program above.
!    See the comments there and in subroutine qdrtf above.
!
  implicit none

  integer n

  real ( kind = 8 ) f
  real ( kind = 8 ) g(n)
  integer nf
  external qdrtf
  integer uip(*)
  real ( kind = 8 ) urp(n,3)
  real ( kind = 8 ) x(n)

  if ( nf /= uip(1) ) then
    call qdrtf ( n, x, nf, f, uip, urp, qdrtf )
  end if

  g(1:n) = urp(1,3) * urp(1:n,2)

  return
end
subroutine qdrtgh ( n, x, nf, g, h, uip, urp, qdrtf )

!*******************************************************************************
!
!! QTRTGH evaluates the gradient and hessian of the objective function f(x).
!
!  Discussion:
!
!    The objective function is described in the main program above.
!    See the comments there and in subroutine qdrtf above.
!    Note that the H returned is
!    the lower triangle of the hessian, stored row-wise.
!
  implicit none

  integer n

  integer nf, uip(*)
  real ( kind = 8 ) x(n), g(n)
  real ( kind = 8 ) h(n*(n+1)/2)
  real ( kind = 8 ) urp(n,3)
  external qdrtf
  integer i, j, k
  real ( kind = 8 ) dn, f, t1, t2

  if ( nf /= uip(1) ) then
    call qdrtf(n, x, nf, f, uip, urp, qdrtf)
  end if

  k = 0
  dn = n
  do i = 1, n
    g(i) = urp(1,3) * urp(i,2)
    t1 = urp(1,3) * urp(i,1)
    t2 = urp(i,2) * urp(2,3)
    do j = 1, i
      k = k + 1
      h(k) = t2*urp(j,2) + t1*urp(j,1)
    end do
    h(k) = h(k) + dn*urp(i,1)*t1
  end do

  return
end
