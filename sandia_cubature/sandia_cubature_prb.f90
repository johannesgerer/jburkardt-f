program main

!*****************************************************************************80
!
!! SANDIA_CUBATURE_PRB tests SANDIA_CUBATURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_CUBATURE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the SANDIA_CUBATURE library.'

  call cn_geg_tests ( )
  call cn_jac_tests ( )
  call cn_leg_tests ( )
  call en_her_tests ( )
  call epn_glg_tests ( )
  call epn_lag_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_CUBATURE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine cn_geg_tests ( )

!*****************************************************************************80
!
!! CN_GEG_TESTS tests the rules for CN with Gegenbauer weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ), dimension ( test_num ) :: alpha_test = (/ &
    - 0.5D+00, 0.0D+00, 0.5D+00, 1.0D+00, 1.5D+00 /)
  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CN_GEG_TESTS'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  CN_GEG, that is, the hypercube [-1,+1]^N, with the'
  write ( *, '(a)' ) '  Gegenbauer weight W(ALPHA;X) = product ( 1 <= I <= N )'
  write ( *, '(a)' ) '    (1-X(I)^2)^ALPHA'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    do test = 1, test_num

      alpha = alpha_test(test)

      expon(1:n) = 0
      call cn_geg_test ( n, alpha, expon )

    end do

    do test = 1, test_num

      alpha = alpha_test(test)

      expon(1:n) = 0
      expon(n) = 1
      call cn_geg_test ( n, alpha, expon )

    end do

    if ( 2 <= n ) then

      do test = 1, test_num

        alpha = alpha_test(test)

        expon(1:n) = 0
        expon(1) = 1
        expon(2) = 1
        call cn_geg_test ( n, alpha, expon )

      end do

    end if

    do test = 1, test_num

      alpha = alpha_test(test)

      expon(1:n) = 0
      expon(1) = 2
      call cn_geg_test ( n, alpha, expon )

    end do

    deallocate ( expon )

  end do

  return
end
subroutine cn_geg_test ( n, alpha, expon )

!*****************************************************************************80
!
!! CN_GEG_TEST tests the rules for CN with Gegenbauer weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real    ( kind = 8 ) delta0
  real    ( kind = 8 ) err
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real    ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) quad
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ), allocatable :: v(:)
  real    ( kind = 8 ) volume_1d
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call cn_geg_monomial_integral ( n, alpha, expon, exact )

  p = 0

  if ( d <= p ) then
    call cn_geg_00_1_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_geg_00_1 ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_GEG_00_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 1

  if ( d <= p ) then
    call cn_geg_01_1_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_geg_01_1 ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_GEG_01_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call cn_geg_02_xiu_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_geg_02_xiu ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_GEG_02_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = 1.0D+00
    delta0 = 0.0D+00
    c1 = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
    volume_1d = sqrt ( pi ) * r8_gamma ( alpha + 1.0D+00 ) &
      / r8_gamma ( alpha + 1.5D+00 )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 3

  if ( d <= p ) then
    call cn_geg_03_xiu_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_geg_03_xiu ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_GEG_03_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine cn_jac_tests ( )

!*****************************************************************************80
!
!! CN_JAC_TESTS tests the rules for CN with Jacobi weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ), dimension ( test_num ) :: alpha_test = (/ &
    0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00 /)
  real    ( kind = 8 ) beta
  real    ( kind = 8 ), dimension ( test_num ) :: beta_test = (/ &
    0.0D+00, 0.0D+00, 2.0D+00, 1.5D+00 /)
  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CN_JAC_TESTS'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  CN_JAC, that is, the hypercube [-1,+1]^N, with the'
  write ( *, '(a)' ) '  Jacobi weight W(ALPHA,BETA;X) = product ( 1 <= I <= N )'
  write ( *, '(a)' ) '    (1-X(I))^ALPHA (1+X(I))^BETA'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    do test = 1, test_num

      alpha = alpha_test(test)
      beta  = beta_test(test)

      expon(1:n) = 0
      call cn_jac_test ( n, alpha, beta, expon )

    end do

    do test = 1, test_num

      alpha = alpha_test(test)
      beta  = beta_test(test)

      expon(1:n) = 0
      expon(n) = 1
      call cn_jac_test ( n, alpha, beta, expon )

    end do

    if ( 2 <= n ) then

      do test = 1, test_num

        alpha = alpha_test(test)
        beta  = beta_test(test)

        expon(1:n) = 0
        expon(1) = 1
        expon(2) = 1
        call cn_jac_test ( n, alpha, beta, expon )

      end do

    end if

    do test = 1, test_num

      alpha = alpha_test(test)
      beta  = beta_test(test)

      expon(1:n) = 0
      expon(1) = 2
      call cn_jac_test ( n, alpha, beta, expon )

    end do

    deallocate ( expon )

  end do

  return
end
subroutine cn_jac_test ( n, alpha, beta, expon )

!*****************************************************************************80
!
!! CN_JAC_TEST tests the rules for CN with Jacobi weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real    ( kind = 8 ) delta0
  real    ( kind = 8 ) err
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real    ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real    ( kind = 8 ) quad
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ), allocatable :: v(:)
  real    ( kind = 8 ) volume_1d
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA =  ', beta
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call cn_jac_monomial_integral ( n, alpha, beta, expon, exact )

  p = 0

  if ( d <= p ) then
    call cn_jac_00_1_size ( n, alpha, beta, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_jac_00_1 ( n, alpha, beta, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_JAC_00_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 1

  if ( d <= p ) then
    call cn_jac_01_1_size ( n, alpha, beta, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_jac_01_1 ( n, alpha, beta, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_JAC_01_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call cn_jac_02_xiu_size ( n, alpha, beta, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_jac_02_xiu ( n, alpha, beta, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_JAC_02_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = ( alpha + beta + 2.0D+00 ) / 2.0D+00
    delta0 = ( alpha - beta ) / 2.0D+00
    c1 = 2.0D+00 * ( alpha + 1.0D+00 ) * ( beta + 1.0D+00 ) &
      / ( alpha + beta + 3.0D+00 ) / ( alpha + beta + 2.0D+00 )
    volume_1d = 2.0D+00 ** ( alpha + beta + 1.0D+00 ) &
      * r8_gamma ( alpha + 1.0D+00 ) * r8_gamma ( beta + 1.0D+00 ) &
      / ( alpha + beta + 1.0D+00 ) / r8_gamma ( alpha + beta + 1.0D+00 )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine cn_leg_tests ( )

!*****************************************************************************80
!
!! CN_LEG_TESTS tests the rules for CN with Legendre weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CN_LEG_TESTS'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  CN_LEG, that is, the hypercube [-1,+1]^N, with the'
  write ( *, '(a)' ) '  Legendre weight W(X) = 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    expon(1:n) = 0
    call cn_leg_test ( n, expon )

    expon(1:n) = 0
    expon(n) = 1
    call cn_leg_test ( n, expon )

    if ( 2 <= n ) then

      expon(1:n) = 0
      expon(1) = 1
      expon(2) = 1
      call cn_leg_test ( n, expon )

    end if

    expon(1:n) = 0
    expon(1) = 2
    call cn_leg_test ( n, expon )

    expon(1:n) = 0
    expon(n) = 3
    call cn_leg_test ( n, expon )

    expon(1:n) = 0
    expon(n) = 4
    call cn_leg_test ( n, expon )

    if ( 2 <= n ) then
      expon(1:n) = 0
      expon(1) = 2
      expon(n) = 3
      call cn_leg_test ( n, expon )
    end if

    deallocate ( expon )

  end do

  return
end
subroutine cn_leg_test ( n, expon )

!*****************************************************************************80
!
!! CN_LEG_TEST tests the rules for CN with Legendre weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real    ( kind = 8 ) delta0
  real    ( kind = 8 ) err
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real    ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real    ( kind = 8 ) quad
  real    ( kind = 8 ), allocatable :: v(:)
  real    ( kind = 8 ) volume_1d
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call cn_leg_monomial_integral ( n, expon, exact )

  p = 1

  if ( d <= p ) then
    call cn_leg_01_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_leg_01_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_01_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call cn_leg_02_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_leg_02_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_02_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = 1.0D+00
    delta0 = 0.0D+00
    c1 = 1.0D+00 / 3.0D+00
    volume_1d = 2.0D+00
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 3

  if ( d <= p ) then

    call cn_leg_03_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_leg_03_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_03_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call cn_leg_03_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_leg_03_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_03_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 5

  if ( d <= p ) then

    if ( 4 <= n .and. n <= 6 ) then
      call cn_leg_05_1_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      option = 1
      call cn_leg_05_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_05_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( 4 <= n .and. n <= 5 ) then
      call cn_leg_05_1_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      option = 2
      call cn_leg_05_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_05_1(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( 2 <= n ) then
      call cn_leg_05_2_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call cn_leg_05_2 ( n, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_05_2:   ', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine en_her_tests ( )

!*****************************************************************************80
!
!! EN_HER_TESTS tests the Stroud EN_HER rules on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EN_HER_TESTS'
  write ( *, '(a)' ) '  Demonstrate the use of Stroud rules for the region'
  write ( *, '(a)' ) '  EN_HER, that is, all of N-dimensional space, with the'
  write ( *, '(a)' ) '  Hermite weight function '
  write ( *, '(a)' ) '    W(X) = product ( 1 <= i <= N ) exp ( - X(I)^2 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 7

    allocate ( expon(1:n) )

    expon(1:n) = 0
    call en_her_test ( n, expon )

    expon(1:n) = 0
    expon(1) = 2
    call en_her_test ( n, expon )

    expon(1:n) = 0
    expon(2) = 4
    call en_her_test ( n, expon )

    expon(1:n) = 0
    i = mod ( 3 - 1, n ) + 1
    expon(i) = 6
    call en_her_test ( n, expon )

    expon(1:n) = 0
    expon(1) = 2
    expon(2) = 4
    call en_her_test ( n, expon )

    expon(1:n) = 0
    i = mod ( 4 - 1, n ) + 1
    expon(i) = 8
    call en_her_test ( n, expon )

    expon(1:n) = 0
    i = mod ( 5 - 1, n ) + 1
    expon(i) = 10
    call en_her_test ( n, expon )

    do i = 1, n
      expon(i) = i - 1
    end do
    call en_her_test ( n, expon )

    expon(1:n) = 2
    call en_her_test ( n, expon )

    deallocate ( expon )

  end do

  return
end
subroutine en_her_test ( n, expon )

!*****************************************************************************80
!
!! EN_HER_TEST tests the Stroud EN_HER rules on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real    ( kind = 8 ) delta0
  real    ( kind = 8 ) err
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real    ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) quad
  real    ( kind = 8 ), allocatable :: v(:)
  real    ( kind = 8 ) volume_1d
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call en_her_monomial_integral ( n, expon, exact )

  p = 1

  if ( d <= p ) then
    call en_her_01_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_her_01_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_01_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call en_her_02_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_her_02_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_02_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = 2.0D+00
    delta0 = 0.0D+00
    c1 = 1.0D+00
    volume_1d = sqrt ( pi )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 3

  if ( d <= p ) then

    call en_her_03_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_her_03_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_03_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call en_her_03_2_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_her_03_2 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_03_2:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call en_her_03_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_her_03_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_03_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 5

  if ( d <= p ) then

    if ( 2 <= n .and. n <= 7 ) then

      option = 1
      call en_her_05_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_05_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_05_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )

      if ( n == 3 .or. n == 5 .or. n == 6 ) then
        option = 2
        call en_her_05_1_size ( n, option, o )
        allocate ( x(n,o) )
        allocate ( w(o) )
        call en_her_05_1 ( n, option, o, x, w )
        allocate ( v(1:o) )
        call monomial_value ( n, o, x, expon, v )
        quad = dot_product ( w(1:o), v(1:o) )
        err = abs ( quad - exact )
        write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_05_1(2):', o, quad, err
        deallocate ( v )
        deallocate ( w )
        deallocate ( x )
      end if

    end if

    call en_her_05_2_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_her_05_2 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_05_2:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    if ( 3 <= n ) then
      call en_her_05_3_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_05_3 ( n, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_05_3:   ', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    call en_her_05_4_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_her_05_4 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_05_4:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call en_her_05_5_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_her_05_5 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_05_5:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    if ( 5 <= n ) then
      call en_her_05_6_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_05_6 ( n, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_05_6:   ', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  p = 7

  if ( d <= p ) then

    if ( n == 3 .or. n == 4 .or. n == 6 .or. n == 7 ) then
      option = 1
      call en_her_07_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_07_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_07_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( n == 3 .or. n == 4 ) then
      option = 2
      call en_her_07_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_07_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_07_1(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( 3 <= n ) then
      call en_her_07_2_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_07_2 ( n, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_07_2:   ', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( 3 <= n .and. n <= 6 ) then
      option = 1
      call en_her_07_3_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_07_3 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_07_3(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )

      option = 2
      call en_her_07_3_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_07_3 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_07_3(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  p = 9

  if ( d <= p ) then

    if ( 3 <= n .and. n <= 6 ) then
      option = 1
      call en_her_09_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_09_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_09_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )

      option = 2
      call en_her_09_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_09_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_09_1(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  p = 11

  if ( d <= p ) then

    if ( 3 <= n .and. n <= 5 ) then
      option = 1
      call en_her_11_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_11_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_11_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )

      option = 2
      call en_her_11_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_her_11_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_HER_11_1(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine epn_glg_tests ( )

!*****************************************************************************80
!
!! EPN_GLG_TESTS tests the rules for EPN with GLG weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ), dimension ( test_num ) :: alpha_test = (/ &
    - 0.5D+00, 0.0D+00, 0.5D+00, 1.0D+00, 2.0D+00 /)
  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EPN_GLG_TESTS'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  EPN_GLG, that is, the positive half space [0,+oo)^N, with the'
  write ( *, '(a)' ) '  Generalized Laguerre weight '
  write ( *, '(a)' ) '  W(ALPHA;X) = product ( 1 <= I <= N ) X(I)^ALPHA exp ( -X(I) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    expon(1:n) = 0
    do test = 1, test_num
      alpha = alpha_test(test)
      call epn_glg_test ( n, expon, alpha )
    end do

    expon(1:n) = 0
    expon(n) = 1
    do test = 1, test_num
      alpha = alpha_test(test)
      call epn_glg_test ( n, expon, alpha )
    end do

    if ( 2 <= n ) then
      expon(1:n) = 0
      expon(1) = 1
      expon(2) = 1
      do test = 1, test_num
        alpha = alpha_test(test)
        call epn_glg_test ( n, expon, alpha )
      end do
    end if

    expon(1:n) = 0
    expon(1) = 2
    do test = 1, test_num
      alpha = alpha_test(test)
      call epn_glg_test ( n, expon, alpha )
    end do

    deallocate ( expon )

  end do

  return
end
subroutine epn_glg_test ( n, expon, alpha )

!*****************************************************************************80
!
!! EPN_GLG_TEST tests the rules for EPN with GLG weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real    ( kind = 8 ) delta0
  real    ( kind = 8 ) err
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real    ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real    ( kind = 8 ) quad
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ), allocatable :: v(:)
  real    ( kind = 8 ) volume_1d
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call epn_glg_monomial_integral ( n, expon, alpha, exact )

  p = 0

  if ( d <= p ) then
    call epn_glg_00_1_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_glg_00_1 ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_GLG_00_1:  ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 1

  if ( d <= p ) then
    call epn_glg_01_1_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_glg_01_1 ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_GLG_01_1:  ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call epn_glg_02_xiu_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_glg_02_xiu ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_GLG_02_XIU:', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = - 1.0D+00
    delta0 = alpha + 1.0D+00
    c1 = - alpha - 1.0D+00
    volume_1d = r8_gamma ( 1.0D+00 + alpha )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine epn_lag_tests ( )

!*****************************************************************************80
!
!! EPN_LAG_TESTS tests the rules for EPN with Laguerre weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EPN_LAG_TESTS'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  EPN_LAG, that is, the positive half space [0,+oo)^N, with'
  write ( *, '(a)' ) '  the Laguerre weight W(X) = product ( 1 <= I <= N ) exp ( -X(I) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    expon(1:n) = 0
    call epn_lag_test ( n, expon )

    expon(1:n) = 0
    expon(n) = 1
    call epn_lag_test ( n, expon )

    if ( 2 <= n ) then
      expon(1:n) = 0
      expon(1) = 1
      expon(2) = 1
      call epn_lag_test ( n, expon )
    end if

    expon(1:n) = 0
    expon(1) = 2
    call epn_lag_test ( n, expon )

    deallocate ( expon )

  end do

  return
end
subroutine epn_lag_test ( n, expon )

!*****************************************************************************80
!
!! EPN_LAG_TEST tests the rules for EPN with Laguerre weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real    ( kind = 8 ) delta0
  real    ( kind = 8 ) err
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real    ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real    ( kind = 8 ) quad
  real    ( kind = 8 ), allocatable :: v(:)
  real    ( kind = 8 ) volume_1d
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call epn_lag_monomial_integral ( n, expon, exact )

  p = 0

  if ( d <= p ) then
    call epn_lag_00_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_lag_00_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_LAG_00_1:  ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 1

  if ( d <= p ) then
    call epn_lag_01_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_lag_01_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_LAG_01_1:  ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call epn_lag_02_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_lag_02_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_LAG_02_XIU:', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = - 1.0D+00
    delta0 = 1.0D+00
    c1 = - 1.0D+00
    volume_1d = 1.0D+00
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
