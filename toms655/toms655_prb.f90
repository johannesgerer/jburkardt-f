program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS655_PRB.
!
!  Discussion:
!
!    TOMS655_PRB is a test program for TOMS655.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) nt

  write ( *, '(a)' ) ' '
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS655_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS655 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
!
!  Compute 15 points of an example of each rule, with default A, B.
!
  do kind = 1, 9
    nt = 15
    if ( kind == 8 ) then
      alpha = 1.0D+00
      beta = - alpha - 2 * nt - 2
    else
      alpha = 0.0D+00
      beta = 0.0D+00
    end if
    call test10 ( nt, kind, alpha, beta )
  end do
!
!  Compute 15 points of an example of each rule using nondefault A, B.
!
  do kind = 1, 9

    nt = 15

    if ( kind == 1 ) then
      alpha = 0.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kind == 2 ) then
      alpha = 0.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kind == 3 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kind == 4 ) then
      alpha = 1.5D+00
      beta = 0.5D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kind == 5 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
      a = 1.0D+00
      b = 1.0D+00
    else if ( kind == 6 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 0.5D+00
    else if ( kind == 7 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kind == 8 ) then
      alpha = 1.0D+00
      beta = - alpha - 2 * nt - 2
      a = 0.0D+00
      b = 1.0D+00
    else if ( kind == 9 ) then
      alpha = 0.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    end if

    call test11 ( nt, kind, alpha, beta, a, b )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS655_PRB'
  write ( *, '(a)' ) '  Normal end of TOMS655 tests.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CIQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts
  real    ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real    ( kind = 8 ), allocatable :: t(:)
  real    ( kind = 8 ), allocatable :: wts(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test CIQFS.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = 8 ) * pi / real ( 2 * nt, kind = 8 ) )
  end do
!
!  Set the knot multiplicities.
!
  allocate ( mlt(1:nt) )
  mlt(1:nt) = 2
!
!  Set the size of the weights array.
!
  nwts = sum ( mlt(1:nt) )
!
!  Because KEY = 1, NDX will be set up for us.
!
  allocate ( ndx(1:nt) )
!
!  KEY = 1 indicates that the WTS array should hold the weights
!  in the usual order.
!
  key = 1
!
!  Request Legendre weight function.
!
  kind = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  LU controls printing.
!  A positive value requests that we compute and print weights, and
!  conduct a moments check.
!
  lu = 6
!
!  This call returns the WTS array.
!
  allocate ( wts(1:nwts) )

  call ciqfs ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, lu, wts )

  deallocate ( mlt )
  deallocate ( ndx )
  deallocate ( t )
  deallocate ( wts )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CIQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts
  real    ( kind = 8 ), allocatable :: t(:)
  real    ( kind = 8 ), allocatable :: wts(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Test CIQF, CIQFS, CGQF and CGQFS'
  write ( *, '(a)' ) '  with all classical weight functions.'
!
!  Try all weight functions.
!
  do kind = 1, 9
!
!  Number of knots.
!
    nt = 5
!
!  Set parameters ALPHA and BETA.
!
    alpha = 0.5D+00
    if ( kind /= 8 ) then
      beta  = 2.0D+00
    else
      beta = - 16.0D+00
    end if
!
!  Set A and B.
!
    a = - 0.5D+00
    b = 2.0D+00
!
!  Have CGQF compute the knots and weights.
!
    allocate ( t(1:nt) )
    allocate ( wts(1:nt) )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Knots and weights of Gauss quadrature formula'
    write ( *, '(a)' ) '  computed by CGQF.'
    call cgqf ( nt, kind, alpha, beta, a, b, t, wts )
!
!  Now compute the weights for the same knots by CIQF.
!
!  Set the knot multiplicities.
!
    allocate ( mlt(1:nt) )
    mlt(1:nt) = 2
!
!  Set the size of the weights array.
!
    nwts = sum ( mlt(1:nt) )
!
!  We need to deallocate and reallocate WTS because it is now of
!  dimension NWTS rather than NT.
!
    deallocate ( wts )
    allocate ( wts(1:nwts) )
!
!  Because KEY = 1, NDX will be set up for us.
!
    allocate ( ndx(1:nt) )
!
!  KEY = 1 indicates that the WTS array should hold the weights
!  in the usual order.
!
    key = 1

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Weights of Gauss quadrature formula computed from the'
    write ( *, '(a)' ) '  knots by CIQF.'
    call ciqf ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, a, b, lu, wts )

    deallocate ( mlt )
    deallocate ( ndx )
    deallocate ( t )
    deallocate ( wts )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests CEIQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  real    ( kind = 8 ), external :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts
  real    ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real    ( kind = 8 ) qfsum
  real    ( kind = 8 ) qfsx
  real    ( kind = 8 ), allocatable :: t(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Test CEIQFS.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = 8 ) * pi / real ( 2 * nt, kind = 8 ) )
  end do
!
!  Set the knot multiplicities.
!
  allocate ( mlt(1:nt) )
  mlt(1:nt) = 2
!
!  Set KIND to the Legendre weight function.
!
  kind = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  Call CEIQFS to set up the quadrature formula and evaluate it on F.
!
  call ceiqfs ( nt, t, mlt, kind, alpha, beta, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integral of sin(x) on -1, 1 by Fejer type rule'
  write ( *, '(a,i4,a,i4)' ) '  with ', nt, ' points of multiplicity 2.'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = cos ( - 1.0D+00 ) - cos ( 1.0D+00 )
  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  deallocate ( mlt )
  deallocate ( t )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests CEIQF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta
  real    ( kind = 8 ), external :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts
  real    ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real    ( kind = 8 ) qfsum
  real    ( kind = 8 ) qfsx
  real    ( kind = 8 ), allocatable :: t(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test CEIQF.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = 8 ) * pi / real ( 2 * nt, kind = 8 ) )
  end do
!
!  Set the knot multiplicities.
!
  allocate ( mlt(1:nt) )
  mlt(1:nt) = 2
!
!  Set KIND to the Legendre weight function.
!
  kind = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  Set nonstandard interval A, B.
!
  a = -0.5D+00
  b = 2.0D+00
!
!  Shift knots from [-1,1] to [A,B].
!
  do i = 1, nt
    t(i) = ( ( b - a ) * t(i) + ( a + b ) ) / 2.0D+00
  end do
!
!  Call CEIQF to set up the quadrature formula and evaluate it on F.
!
  call ceiqf ( nt, t, mlt, kind, alpha, beta, a, b, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Integral of sin(x) from ', a, ' to ', b
  write ( *, '(a,i4,a)' ) '  by Fejer type rule with ', nt, ' points'
  write ( *, '(a)' ) '  of multiplicity 2.'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = cos ( a ) - cos ( b )
  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  deallocate ( mlt )
  deallocate ( t )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests CLIQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) nt
  real    ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real    ( kind = 8 ), allocatable :: t(:)
  real    ( kind = 8 ), allocatable :: wts(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Test CLIQFS.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = 8 ) * pi / real ( 2 * nt, kind = 8 ) )
  end do
!
!  Request Legendre weight function.
!
  kind = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  LU controls printing.
!  A positive value requests that we compute and print weights, and
!  conduct a moments check.
!
  lu = 6
!
!  This call returns the WTS array.
!
  allocate ( wts(1:nt) )

  call cliqfs ( nt, t, kind, alpha, beta, lu, wts )

  deallocate ( t )
  deallocate ( wts )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests CEIQF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta
  real    ( kind = 8 ), external :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) nt
  real    ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real    ( kind = 8 ) qfsum
  real    ( kind = 8 ) qfsx
  real    ( kind = 8 ), allocatable :: t(:)
  real    ( kind = 8 ), allocatable :: wts(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Test CLIQF and EIQFS.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = 8 ) * pi / real ( 2 * nt, kind = 8 ) )
  end do
!
!  Set KIND to the Legendre weight function.
!
  kind = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  Set nonstandard interval A, B.
!
  a = -0.5D+00
  b = 2.0D+00
!
!  Shift knots from [-1,1] to [A,B].
!
  do i = 1, nt
    t(i) = ( ( b - a ) * t(i) + ( a + b ) ) / 2.0D+00
  end do
!
!  LU controls printout.
!
  lu = 6
!
!  Allocate space for WTS.
!
  allocate ( wts(1:nt) )
!
!  Call CLIQF to set up the quadrature formula.
!
  call cliqf ( nt, t, kind, alpha, beta, a, b, lu, wts )
!
!  Call EIQFS to evaluate the quadrature formula.
!
  call eiqfs ( nt, t, wts, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Integral of sin(x) from ', a, ' to ', b
  write ( *, '(a,i4,a)' ) '  by Fejer type rule with ', nt, ' points'
  write ( *, '(a)' ) '  of multiplicity 1.'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = cos ( a ) - cos ( b )
  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  deallocate ( t )
  deallocate ( wts )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests CEGQF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta
  real    ( kind = 8 ), external :: f
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) nt
  real    ( kind = 8 ) qfsum
  real    ( kind = 8 ) qfsx

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Test CEGQF.'
!
!  Number of knots.
!
  nt = 12
!
!  Request exponential weight function.
!
  kind = 7
!
!  Set ALPHA and BETA.
!
  alpha = 1.0D+00
  beta  = 0.0D+00
!
!  Set interval [A,B].
!
  a = -0.5D+00
  b = 2.0D+00
!
!  Call CEGQF to compute and evaluate the Gauss quadrature formula.
!
  call cegqf ( nt, kind, alpha, beta, a, b, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Integral of x*sin(x) from ', a, ' to ', b
  write ( *, '(a,i4,a)' ) '  by Gauss-exponential rule with ', nt, ' points'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = ( b - a ) * 0.5D+00 * ( cos ( a ) - cos ( b ) ) &
    + sin ( b ) + sin ( a ) - 2.0D+00 * sin ( ( a + b ) / 2.0D+00 )

  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests CEGQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  real    ( kind = 8 ), external :: f
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) nt
  real    ( kind = 8 ) qfsum
  real    ( kind = 8 ) qfsx

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Test CEGQFS.'
!
!  Number of knots.
!
  nt = 12
!
!  Request exponential weight function.
!
  kind = 7
!
!  Set ALPHA and BETA.
!
  alpha = 1.0D+00
  beta  = 0.0D+00
!
!  Call CEGQFS to compute and evaluate the Gauss quadrature formula.
!
  call cegqfs ( nt, kind, alpha, beta, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integral of x*sin(x) from -1 to +1'
  write ( *, '(a,i4,a)' ) '  by Gauss-exponential rule with ', nt, ' points'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = cos ( -1.0D+00 ) - cos ( +1.0D+00 )

  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 calls CGQFS to compute and print generalized Gauss-Hermite rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) io
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) nt
  real    ( kind = 8 ), allocatable :: t(:)
  real    ( kind = 8 ), allocatable :: wts(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  Call CGQFS to compute generalized Hermite rules.'

  nt = 15
  kind = 6
  alpha = 1.0D+00
  beta = 0.0D+00
  io = - 6
  allocate ( t(nt) )
  allocate ( wts(nt) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NT = ', nt
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha

  call cgqfs ( nt, kind, alpha, beta, io, t, wts )

  deallocate ( t )
  deallocate ( wts )

  return
end
subroutine test10 ( nt, kind, alpha, beta )

!*****************************************************************************80
!
!! TEST10 calls CDGQF to compute a quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nt

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real    ( kind = 8 ) t(nt)
  real    ( kind = 8 ) wts(nt)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  Call CDGQF to compute a quadrature formula.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  KIND = ', kind
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA  = ', beta

  call cdgqf ( nt, kind, alpha, beta, t, wts )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Index     Abscissas                 Weights'
  write ( *, '(a)' ) ' '
  do i = 1, nt
    write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) i, t(i), wts(i)
  end do

  return
end
subroutine test11 ( nt, kind, alpha, beta, a, b )

!*****************************************************************************80
!
!! TEST11 calls CGQF to compute a quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nt

  real    ( kind = 8 ) a
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real    ( kind = 8 ) t(nt)
  real    ( kind = 8 ) wts(nt)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Call CGQF to compute a quadrature formula'
  write ( *, '(a)' ) '  with nondefault values of parameters A, B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  KIND = ', kind
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA  = ', beta
  write ( *, '(a,g14.6)' ) '  A =     ', a
  write ( *, '(a,g14.6)' ) '  B =     ', b

  call cgqf ( nt, kind, alpha, beta, a, b, t, wts )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Index     Abscissas                 Weights'
  write ( *, '(a)' ) ' '
  do i = 1, nt
    write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) i, t(i), wts(i)
  end do

  return
end
function f ( x, i )

!*****************************************************************************80
!
!! F returns values of the integrand or its derivatives.
!
!  Discussion:
!
!    This function is an example of an integrand function.
!
!    The package can generate quadrature formulas that use derivative
!    information as well as function values.  Therefore, this routine is
!    set up to provide derivatives of any order as well as the function
!    value.  In an actual application, the highest derivative needed
!    is of order one less than the highest knot multiplicity.
!
!    In other words, in the usual case where knots are not repeated,
!    this routine only needs to return function values, not any derivatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Input, integer ( kind = 4 ) I, the order of the derivative of F to
!    be evaluated.
!
!    Output, real ( kind = 8 ) F, the value of the I-th derivative of F at X.
!
  implicit none

  real    ( kind = 8 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  real    ( kind = 8 ) x

  l = mod ( i, 4 )

  if ( l == 0 ) then
    f = sin ( x )
  else if ( l == 1 ) then
    f = cos ( x )
  else if ( l == 2 ) then
    f = - sin ( x )
  else if ( l == 3 ) then
    f = - cos ( x )
  end if

  return
end
