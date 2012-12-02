program main

!*****************************************************************************80
!
!! MAIN is the main program for TANH_QUAD_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) p

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TANH_QUAD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TANH_QUAD library.'

  call test01 ( )
  call test012 ( )
  call test02 ( )
  call test025 ( )
  call test03 ( )

  do p = 1, 14
    call test032 ( p )
  end do

  do p = 1, 14
    call test04 ( p )
  end do

  do p = 1, 14
    call test05 ( p )
  end do

  do p = 1, 14
    call test06 ( p )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TANH_QUAD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates TANH_M_TO_H and TANH_H_TO_N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real    ( kind = 8 ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  TANH_M_TO_H determines the spacing H from level M'
  write ( *, '(a)' ) '  TANH_H_TO_N determines the quadrature order N from'
  write ( *, '(a)' ) '  the spacing H and tolerance TOL.'

  tol = epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g16.8)' ) '  All tests use TOL = ', tol
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M        H                  N'
  write ( *, '(a)' ) ' '

  do m = -2, 8

    call tanh_m_to_h ( m, h )

    call tanh_h_to_n ( h, tol, n )

    write (  *, '(2x,i8,2x,g16.8,2x,i8)' ) m, h, n
  end do

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 computes a midpoint quadrature rule W, X with N = 0, 1, 3, 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  Determine nested midpoint quadrature rules W, X'
  write ( *, '(a)' ) '  by choosing N = 0, 1, 3, 7.'

  do m = 1, 4

    order = 2**m - 1
    n = ( ( order + 1 ) / 2 ) - 1
    h = 2.0D+00 / real ( order + 1, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,i8,a,i8,a,g16.8)' ) &
      '  M = ', m, '  ORDER = ', order, '  N = ', n, '  H = ', h

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call midpoint_rule ( n, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '         I       Wi                Xi'
    write ( *, '(a)' ) ' '

    do i = - n, n
      write (  *, '(2x,i8,2x,g16.8,2x,g16.8)' ) i, w(i), x(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 computes a tanh quadrature rule W, X with N = 5, 10, 20.
!
!  Discussion:
!
!    These results should match the following values reported in
!    Kahaner, Moler and Nash:
!
!       I    Wi         Xi
!      --  --------  --------
!      -5  0.000471  -0.999737
!      -4  0.002807  -0.998428
!      -3  0.016631  -0.990649
!      -2  0.094844  -0.945434
!      -1  0.439127  -0.713098
!       0  0.893459   0.000000
!       1  0.439127   0.713098
!       2  0.094844   0.945434
!       3  0.016631   0.990649
!       4  0.002807   0.998428
!       5  0.000471   0.999737
!
!    Note that these values do not sum to 2, although they come close!
!    Thus, a fundamental feature of most quadrature rules is ignored here.
!    This rule will not integrate f(x) = 1 exactly.  But it is not a
!    family of rules based on polynomial accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
  implicit none

  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real    ( kind = 8 ) tol
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Determine specific tanh quadrature rules W, X'
  write ( *, '(a)' ) '  by choosing N = 5, 10, 20.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rule for N = 5 appears in the reference'
  write ( *, '(a)' ) '  Kahaner, Moler and Nash.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Also note that the rule for H(N) means that'
  write ( *, '(a)' ) '  rules for doubled N do not nest.'

  tol = epsilon ( tol )

  n = 5

  do j = 1, 3

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Quadrature order N = ', n
    write ( *, '(a,g16.8)' ) '  H = ', h

    call tanh_n_to_h ( n, h )

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call tanh_rule ( n, h, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '         I       Wi                Xi'
    write ( *, '(a)' ) ' '

    do i = - n, n
      write (  *, '(2x,i8,2x,g16.8,2x,g16.8)' ) i, w(i), x(i)
    end do

    deallocate ( w )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Note that, especially for low N, the weights need not sum to 2!:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       H             sum(W)'
  write ( *, '(a)' ) ' '

  n = 5

  do j = 1, 10

    call tanh_n_to_h ( n, h )

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call tanh_rule ( n, h, x, w )

    write ( *, '(2x,i8,2x,g14.6,2x,g16.8)' ) n, h, sum ( w(-n:n) )

    deallocate ( w )
    deallocate ( x )

    n = 2 * n

  end do


  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 computes tanh sinh quadrature rules.
!
!  Discussion:
!
!    We are seeking a family of nested rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
  implicit none

  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  Determine tanh-sinh quadrature rules W, X'
  write ( *, '(a)' ) '  for N = 0, 1, 3, 7, 15, 31, 63.'

  do m = -3, 3

    order = 2**( m + 4 ) - 1
    n = ( ( order + 1 ) / 2 ) - 1
    h = 4.0D+00 / real ( order + 1, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,i8,a,i8,a,g16.8)' ) &
      '  M = ', m, '  ORDER = ', order, '  N = ', n, '  H = ', h

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call tanh_sinh_rule ( n, h, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '         I       Wi                Xi'
    write ( *, '(a)' ) ' '

    do i = - n, n
      write (  *, '(2x,i8,2x,g16.8,2x,g16.8)' ) i, w(i), x(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Note that, especially for low N, the weights need not sum to 2!:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       H              sum(W)'
  write ( *, '(a)' ) ' '

  do m = -3, 10

    order = 2**(m+4) - 1
    n = ( ( order + 1 ) / 2 ) - 1
    h = 1.0D+00 / 2.0D+00**m

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call tanh_sinh_rule ( n, h, x, w )

    write ( *, '(2x,i8,2x,g16.8,2x,g16.8)' ) n, h, sum ( w(-n:n) )

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 computes a quadrature rule W, X based on a tolerance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real    ( kind = 8 ) tol
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Determine a quadrature rule W, X by specifying'
  write ( *, '(a)' ) '  a tolerance.'

  tol = epsilon ( tol )

  write ( *, '(a,g16.8)' ) '  Tolerance TOL = ', tol

  do m = -1, 2

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Level M = ', m
    call tanh_m_to_h ( m, h )
    write ( *, '(a,g16.8)' ) '  Spacing H = ', h
    call tanh_h_to_n ( h, tol, n )
    write ( *, '(a,i8)' ) '  Quadrature order N = ', n

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call tanh_rule ( n, h, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '         I      Wi                Xi'
    write ( *, '(a)' ) ' '

    do i = - n, n
      write (  *, '(2x,i8,2x,g16.8,2x,g16.8)' ) i, w(i), x(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test032 ( p )

!*****************************************************************************80
!
!! TEST032 applies a sequence of midpoint rules to a test integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the problem index.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) error
  real    ( kind = 8 ) exact
  real    ( kind = 8 ), allocatable :: fx(:)
  real    ( kind = 8 ) h
  real    ( kind = 8 ) h_midpoint
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  real    ( kind = 8 ) quad
  real    ( kind = 8 ) tol
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032'
  write ( *, '(a)' ) '  Apply a sequence of midpoint rules to a test integral.'

  tol = epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Problem index P = ', p
  write ( *, '(a,g16.8)' ) '  Tolerance TOL = ', tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M      H_MIDPOINT           N          Exact           Quad            Error'
  write ( *, '(a)' ) ' '

  do m = 1, 8
!
!  Choose N based on the value of H that would be used by the TANH rule.
!
    call tanh_m_to_h ( m, h )
    call tanh_h_to_n ( h, tol, n )

    h_midpoint = 2.0 / real ( 2 * ( n + 1 ), kind = 8 )

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call midpoint_rule ( n, x, w )
!
!  Adjust the rule from [-1,+1] to the actual integration limits.
!
    a = -1.0D+00
    b = +1.0D+00

    call p00_ab ( p, c, d )

    call rule_adjust ( a, b, c, d, 2*n+1, x, w )
!
!  Evaluate the integrand.
!
    allocate ( fx(-n:n) )

    call p00_f ( p, 2*n+1, x(-n:n), fx(-n:n) )
!
!  Form the quadrature estimate.
!
    quad = dot_product ( w(-n:n), fx(-n:n) )

    call p00_e ( p, exact )

    error = abs ( exact - quad )

    write ( *, '(2x,i8,2x,g16.8,2x,i8,2x,g16.8,2x,g16.8,2x,g16.8)' ) &
      m, h_midpoint, n, exact, quad, error

    deallocate ( fx )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test04 ( p )

!*****************************************************************************80
!
!! TEST04 applies a sequence of trapezoid rules to a test integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the problem index.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) error
  real    ( kind = 8 ) exact
  real    ( kind = 8 ), allocatable :: fx(:)
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  real    ( kind = 8 ) quad
  real    ( kind = 8 ) tol
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Apply a sequence of trapezoid rules to a test integral.'

  tol = epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Problem index P = ', p
  write ( *, '(a,g16.8)' ) '  Tolerance TOL = ', tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M      H                    N          Exact           Quad            Error'
  write ( *, '(a)' ) ' '

  do m = 1, 8

    call tanh_m_to_h ( m, h )

    call tanh_h_to_n ( h, tol, n )

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call trap_rule ( n, x, w )
!
!  Adjust the rule from [-1,+1] to the actual integration limits.
!
    a = -1.0D+00
    b = +1.0D+00

    call p00_ab ( p, c, d )

    call rule_adjust ( a, b, c, d, 2*n+1, x, w )
!
!  Evaluate the integrand.
!
    allocate ( fx(-n:n) )

    call p00_f ( p, 2*n+1, x(-n:n), fx(-n:n) )
!
!  Form the quadrature estimate.
!
    quad = dot_product ( w(-n:n), fx(-n:n) )

    call p00_e ( p, exact )

    error = abs ( exact - quad )

    write ( *, '(2x,i8,2x,g16.8,2x,i8,2x,g16.8,2x,g16.8,2x,g16.8)' ) &
      m, h, n, exact, quad, error

    deallocate ( fx )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test05 ( p )

!*****************************************************************************80
!
!! TEST05 applies a sequence of tanh rules to a test integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the problem index.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) error
  real    ( kind = 8 ) exact
  real    ( kind = 8 ), allocatable :: fx(:)
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  real    ( kind = 8 ) quad
  real    ( kind = 8 ) tol
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Apply a sequence of tanh rules to a test integral.'

  tol = epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Problem index P = ', p
  write ( *, '(a,g16.8)' ) '  Tolerance TOL = ', tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M      H                    N        Exact           Quad            Error'
  write ( *, '(a)' ) ' '

  do m = -2, 8

    call tanh_m_to_h ( m, h )

    call tanh_h_to_n ( h, tol, n )

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call tanh_rule ( n, h, x, w )
!
!  Adjust the rule from [-1,+1] to the actual integration limits.
!
    a = -1.0D+00
    b = +1.0D+00

    call p00_ab ( p, c, d )

    call rule_adjust ( a, b, c, d, 2*n+1, x, w )
!
!  Evaluate the integrand.
!
    allocate ( fx(-n:n) )

    call p00_f ( p, 2*n+1, x(-n:n), fx(-n:n) )
!
!  Form the quadrature estimate.
!
    quad = dot_product ( w(-n:n), fx(-n:n) )

    call p00_e ( p, exact )

    error = abs ( exact - quad )

    write ( *, '(2x,i8,2x,g16.8,2x,i8,2x,g16.8,2x,g16.8,2x,g16.8)' ) &
      m, h, n, exact, quad, error

    deallocate ( fx )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test06 ( p )

!*****************************************************************************80
!
!! TEST06 applies a sequence of tanh-sinh rules to a test integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the problem index.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) error
  real    ( kind = 8 ) exact
  real    ( kind = 8 ), allocatable :: fx(:)
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  real    ( kind = 8 ) quad
  real    ( kind = 8 ) tol
  real    ( kind = 8 ), allocatable :: w(:)
  real    ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Apply a sequence of tanh-sinh rules to a test integral.'

  tol = epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Problem index P = ', p
  write ( *, '(a,g16.8)' ) '  Tolerance TOL = ', tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M      H                    N        Exact           Quad            Error'
  write ( *, '(a)' ) ' '

  do m = -2, 8

    call tanh_m_to_h ( m, h )

    call tanh_sinh_h_to_n ( h, tol, n )

    allocate ( w(-n:n) )
    allocate ( x(-n:n) )

    call tanh_sinh_rule ( n, h, x, w )
!
!  Adjust the rule from [-1,+1] to the actual integration limits.
!
    a = -1.0D+00
    b = +1.0D+00

    call p00_ab ( p, c, d )

    call rule_adjust ( a, b, c, d, 2*n+1, x, w )
!
!  Evaluate the integrand.
!
    allocate ( fx(-n:n) )

    call p00_f ( p, 2*n+1, x(-n:n), fx(-n:n) )
!
!  Form the quadrature estimate.
!
    quad = dot_product ( w(-n:n), fx(-n:n) )

    call p00_e ( p, exact )

    error = abs ( exact - quad )

    write ( *, '(2x,i8,2x,g16.8,2x,i8,2x,g16.8,2x,g16.8,2x,g16.8)' ) &
      m, h, n, exact, quad, error

    deallocate ( fx )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine p00_ab ( p, a, b )

!*****************************************************************************80
!
!! P00_AB returns the integration limits for a given problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the problem index.
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  integer ( kind = 4 ) p

  if ( p == 1 ) then
    call p01_ab ( a, b )
  else if ( p == 2 ) then
    call p02_ab ( a, b )
  else if ( p == 3 ) then
    call p03_ab ( a, b )
  else if ( p == 4 ) then
    call p04_ab ( a, b )
  else if ( p == 5 ) then
    call p05_ab ( a, b )
  else if ( p == 6 ) then
    call p06_ab ( a, b )
  else if ( p == 7 ) then
    call p07_ab ( a, b )
  else if ( p == 8 ) then
    call p08_ab ( a, b )
  else if ( p == 9 ) then
    call p09_ab ( a, b )
  else if ( p == 10 ) then
    call p10_ab ( a, b )
  else if ( p == 11 ) then
    call p11_ab ( a, b )
  else if ( p == 12 ) then
    call p12_ab ( a, b )
  else if ( p == 13 ) then
    call p13_ab ( a, b )
  else if ( p == 14 ) then
    call p14_ab ( a, b )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_AB - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of P = ', p
    stop
  end if

  return
end
subroutine p00_e ( p, exact )

!*****************************************************************************80
!
!! P00_E returns the exact value of the integral for a given problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the problem index.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact
  integer ( kind = 4 ) p

  if ( p == 1 ) then
    call p01_e ( exact )
  else if ( p == 2 ) then
    call p02_e ( exact )
  else if ( p == 3 ) then
    call p03_e ( exact )
  else if ( p == 4 ) then
    call p04_e ( exact )
  else if ( p == 5 ) then
    call p05_e ( exact )
  else if ( p == 6 ) then
    call p06_e ( exact )
  else if ( p == 7 ) then
    call p07_e ( exact )
  else if ( p == 8 ) then
    call p08_e ( exact )
  else if ( p == 9 ) then
    call p09_e ( exact )
  else if ( p == 10 ) then
    call p10_e ( exact )
  else if ( p == 11 ) then
    call p11_e ( exact )
  else if ( p == 12 ) then
    call p12_e ( exact )
  else if ( p == 13 ) then
    call p13_e ( exact )
  else if ( p == 14 ) then
    call p14_e ( exact )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_E - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of P = ', p
    stop
  end if

  return
end
subroutine p00_f ( p, n, x, fx )

!*****************************************************************************80
!
!! P00_F evaluates the integrand for a given problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the problem index.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  integer ( kind = 4 ) p
  real    ( kind = 8 ) x(n)

  if ( p == 1 ) then
    call p01_f ( n, x, fx )
  else if ( p == 2 ) then
    call p02_f ( n, x, fx )
  else if ( p == 3 ) then
    call p03_f ( n, x, fx )
  else if ( p == 4 ) then
    call p04_f ( n, x, fx )
  else if ( p == 5 ) then
    call p05_f ( n, x, fx )
  else if ( p == 6 ) then
    call p06_f ( n, x, fx )
  else if ( p == 7 ) then
    call p07_f ( n, x, fx )
  else if ( p == 8 ) then
    call p08_f ( n, x, fx )
  else if ( p == 9 ) then
    call p09_f ( n, x, fx )
  else if ( p == 10 ) then
    call p10_f ( n, x, fx )
  else if ( p == 11 ) then
    call p11_f ( n, x, fx )
  else if ( p == 12 ) then
    call p12_f ( n, x, fx )
  else if ( p == 13 ) then
    call p13_f ( n, x, fx )
  else if ( p == 14 ) then
    call p14_f ( n, x, fx )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of P = ', p
    stop
  end if

  return
end
subroutine p01_ab ( a, b )

!*****************************************************************************80
!
!! P01_AB returns the integration limits for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p01_e ( exact )

!*****************************************************************************80
!
!! P01_E returns the exact value of the integral for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact

  exact = 1.0D+00

  return
end
subroutine p01_f ( n, x, fx )

!*****************************************************************************80
!
!! P01_F evaluates the integrand for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00

  return
end
subroutine p02_ab ( a, b )

!*****************************************************************************80
!
!! P02_AB returns the integration limits for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p02_e ( exact )

!*****************************************************************************80
!
!! P02_E returns the exact value of the integral for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact

  exact = 1.0D+00

  return
end
subroutine p02_f ( n, x, fx )

!*****************************************************************************80
!
!! P02_F evaluates the integrand for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = 2.0D+00 * x(1:n)

  return
end
subroutine p03_ab ( a, b )

!*****************************************************************************80
!
!! P03_AB returns the integration limits for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p03_e ( exact )

!*****************************************************************************80
!
!! P03_E returns the exact value of the integral for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact

  exact = 1.0D+00

  return
end
subroutine p03_f ( n, x, fx )

!*****************************************************************************80
!
!! P03_F evaluates the integrand for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = 3.0D+00 * x(1:n) * x(1:n)

  return
end
subroutine p04_ab ( a, b )

!*****************************************************************************80
!
!! P04_AB returns the integration limits for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = -1.0D+00
  b = +1.0D+00

  return
end
subroutine p04_e ( exact )

!*****************************************************************************80
!
!! P04_E returns the exact value of the integral for problem 04.
!
!  Discussion:
!
!    The 20 digit estimate for the exact value comes from Mathematica.
!
!    N [ Integrate [ Exp[-x*x]*Log[1 + x], { x, -1, 1 } ], 20 ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact

  exact = -0.31671419631358172053D+00

  return
end
subroutine p04_f ( n, x, fx )

!*****************************************************************************80
!
!! P04_F evaluates the integrand for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = exp ( - x(1:n) * x(1:n) ) * log ( 1.0D+00 - x(1:n) )

  return
end
subroutine p05_ab ( a, b )

!*****************************************************************************80
!
!! P05_AB returns the integration limits for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p05_e ( exact )

!*****************************************************************************80
!
!! P05_E returns the exact value of the integral for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact

  exact = 0.25D+00

  return
end
subroutine p05_f ( n, x, fx )

!*****************************************************************************80
!
!! P05_F evaluates the integrand for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = x(1:n) * log ( 1.0D+00 + x(1:n) )

  return
end
subroutine p06_ab ( a, b )

!*****************************************************************************80
!
!! P06_AB  returns the integration limits for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p06_e ( exact )

!*****************************************************************************80
!
!! P06_E returns the exact value of the integral for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = ( pi - 2.0D+00 + 2.0D+00 * log ( 2.0D+00 ) ) / 12.0D+00

  return
end
subroutine p06_f ( n, x, fx )

!*****************************************************************************80
!
!! P06_F evaluates the integrand for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = x(1:n) * x(1:n) * atan ( x(1:n) )

  return
end
subroutine p07_ab ( a, b )

!*****************************************************************************80
!
!! P07_AB  returns the integration limits for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = pi / 2.0D+00

  return
end
subroutine p07_e ( exact )

!*****************************************************************************80
!
!! P07_E returns the exact value of the integral for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = ( exp ( pi / 2.0D+00 ) - 1.0D+00 ) / 2.0D+00

  return
end
subroutine p07_f ( n, x, fx )

!*****************************************************************************80
!
!! P07_F evaluates the integrand for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = exp ( x(1:n) ) * cos ( x(1:n) )

  return
end
subroutine p08_ab ( a, b )

!*****************************************************************************80
!
!! P08_AB  returns the integration limits for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p08_e ( exact )

!*****************************************************************************80
!
!! P08_E returns the exact value of the integral for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = ( 5.0D+00 * pi * pi ) / 96.0D+00

  return
end
subroutine p08_f ( n, x, fx )

!*****************************************************************************80
!
!! P08_F evaluates the integrand for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = atan ( sqrt ( 2.0D+00 + x(1:n) * x(1:n) ) ) &
    / ( 1.0D+00 + x(1:n) * x(1:n) ) &
    / sqrt ( 2.0D+00 + x(1:n) * x(1:n) )

  return
end
subroutine p09_ab ( a, b )

!*****************************************************************************80
!
!! P09_AB  returns the integration limits for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p09_e ( exact )

!*****************************************************************************80
!
!! P09_E returns the exact value of the integral for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact

  exact = - 4.0D+00 / 9.0D+00

  return
end
subroutine p09_f ( n, x, fx )

!*****************************************************************************80
!
!! P09_F evaluates the integrand for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) =  sqrt ( x(1:n) ) * log ( x(1:n) )

  return
end
subroutine p10_ab ( a, b )

!*****************************************************************************80
!
!! P10_AB  returns the integration limits for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p10_e ( exact )

!*****************************************************************************80
!
!! P10_E returns the exact value of the integral for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = pi / 4.0D+00

  return
end
subroutine p10_f ( n, x, fx )

!*****************************************************************************80
!
!! P10_F evaluates the integrand for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = sqrt ( 1.0D+00 - x(1:n) * x(1:n) )

  return
end
subroutine p11_ab ( a, b )

!*****************************************************************************80
!
!! P11_AB  returns the integration limits for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p11_e ( exact )

!*****************************************************************************80
!
!! P11_E returns the exact value of the integral for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = 2.0D+00 * sqrt ( pi ) * gamma ( 0.75D+00 ) / gamma ( 0.25D+00 )

  return
end
subroutine p11_f ( n, x, fx )

!*****************************************************************************80
!
!! P11_F evaluates the integrand for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = sqrt ( x(1:n) ) &
    / sqrt ( 1.0D+00 - x(1:n) * x(1:n) )

  return
end
subroutine p12_ab ( a, b )

!*****************************************************************************80
!
!! P12_AB  returns the integration limits for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p12_e ( exact )

!*****************************************************************************80
!
!! P12_E returns the exact value of the integral for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact

  exact = 2.0D+00

  return
end
subroutine p12_f ( n, x, fx )

!*****************************************************************************80
!
!! P12_F evaluates the integrand for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = ( log ( x(1:n) ) )**2

  return
end
subroutine p13_ab ( a, b )

!*****************************************************************************80
!
!! P13_AB  returns the integration limits for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = pi / 2.0D+00

  return
end
subroutine p13_e ( exact )

!*****************************************************************************80
!
!! P13_E returns the exact value of the integral for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = - pi * log ( 2.0D+00 ) / 2.0D+00

  return
end
subroutine p13_f ( n, x, fx )

!*****************************************************************************80
!
!! P13_F evaluates the integrand for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = log ( cos ( x(1:n) ) )

  return
end
subroutine p14_ab ( a, b )

!*****************************************************************************80
!
!! P14_AB  returns the integration limits for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the integration limits.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = pi / 2.0D+00

  return
end
subroutine p14_e ( exact )

!*****************************************************************************80
!
!! P14_E returns the exact value of the integral for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real    ( kind = 8 ) exact
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = pi * sqrt ( 2.0D+00 ) / 2.0D+00

  return
end
subroutine p14_f ( n, x, fx )

!*****************************************************************************80
!
!! P14_F evaluates the integrand for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Karthik Jeyabalan, Xiaoye Li,
!    A Comparison of Three High-Precision Quadrature Schemes,
!    Experimental Mathematics,
!    Volume 14, Number 3, pages 317-329.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fx(n)
  real    ( kind = 8 ) x(n)

  fx(1:n) = sqrt ( tan ( x(1:n) ) )

  return
end
