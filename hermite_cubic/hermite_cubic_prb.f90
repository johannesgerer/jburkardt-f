program main

!*****************************************************************************80
!
!! MAIN is the main program for HERMITE_CUBIC_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( );
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HERMITE_CUBIC library.'

  call hermite_cubic_test01 ( )
  call hermite_cubic_test02 ( )
  call hermite_cubic_test03 ( )
  call hermite_cubic_test04 ( )
  call hermite_cubic_test05 ( )
  call hermite_cubic_test06 ( )
  call hermite_cubic_test07 ( )
  call hermite_cubic_test08 ( )
  call hermite_cubic_test09 ( )

  call hermite_cubic_test10 ( )
  call hermite_cubic_test11 ( )
  call hermite_cubic_test12 ( )
  call hermite_cubic_test13 ( )
  call hermite_cubic_test14 ( )
  call hermite_cubic_test15 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( );

  stop
end
subroutine hermite_cubic_test01 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST01 tests HERMITE_CUBIC_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 1

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_interval
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST01:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.'
  write ( *, '(a)' ) '  Try out four sets of data:'
  write ( *, '(a)' ) '  (F1,D1,F2,D2) = (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)'
  write ( *, '(a)' ) '  on [0,1] and [1.0,-2.0] (interval reversed)'

  do x_interval = 1, 2

    if ( x_interval == 1 ) then
      x1 = 0.0D+00
      x2 = 1.0D+00
    else
      x1 = 1.0D+00
      x2 = -2.0D+00
    end if

    do i = 1, 4

      f1 = 0.0D+00
      d1 = 0.0D+00
      f2 = 0.0D+00
      d2 = 0.0D+00

      if ( i == 1 ) then
        f1 = 1.0D+00
      else if ( i == 2 ) then
        d1 = 1.0D+00
      else if ( i == 3 ) then
        f2 = 1.0D+00
      else if ( i == 4 ) then
        d2 = 1.0D+00
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    J      X           F           D'
      write ( *, '(a)' ) ' '

      do j = - 3, 12

        x = ( real ( 10 - j, kind = 8 ) * x1   &
            + real (      j, kind = 8 ) * x2 ) &
            / real ( 10,     kind = 8 )

        call hermite_cubic_value ( x1, f1, d1, x2, f2, d2, 1, x, f, d, s, t )

        if ( j == 0 ) then
          write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) '*Data', x1, f1, d1
        end if
        write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,f10.4)' ) j, x(1), f(1), d(1)
        if ( j == 10 ) then
          write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) '*Data', x2, f2, d2
        end if

      end do

    end do

  end do

  return
end
subroutine hermite_cubic_test02 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST02 tests HERMITE_CUBIC_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 1

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dc
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fc
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) j
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) sc
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tc
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_interval
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST02:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.'
  write ( *, '(a)' ) '  Try out data from a cubic function:'
  write ( *, '(a)' ) '  on [0,10] and [-1.0,1.0] and [0.5,0.75]'

  do x_interval = 1, 3

    if ( x_interval == 1 ) then
      x1 = 0.0D+00
      x2 = 10.0D+00
    else if ( x_interval == 2 ) then
      x1 = -1.0D+00
      x2 = +1.0D+00
    else if ( x_interval == 3 ) then
      x1 = 0.5D+00
      x2 = 0.75D+00
    end if

    call cubic_value ( x1, f1, d1, s1, t1 )
    call cubic_value ( x2, f2, d2, s2, t2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    J      X           F           D           S           T'

    do j = - 3, 12

      x(1) = ( ( 10 - j ) * x1   &
             +        j   * x2 ) &
               / 10.0D+00

      call hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, d, s, t )
      call cubic_value ( x, fc, dc, sc, tc )

      write ( *, '(a)' ) ' '
      if ( j == 0 ) then
        write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) '*Data', x1, f1, d1
      end if
      write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 'Exact', x(1), fc, dc, sc, tc
      write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) j, x(1), f(1), d(1), s(1), t(1)
      if ( j == 10 ) then
        write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) '*Data', x2, f2, d2
      end if

    end do

  end do

  return
end
subroutine hermite_cubic_test03 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST03 tests HERMITE_CUBIC_INTEGRATE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) j
  real ( kind = 8 ) q_computed
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  integer ( kind = 4 ) x_interval
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST03:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic'
  write ( *, '(a)' ) '  polynomial from A to B.'

  do x_interval = 1, 3

    if ( x_interval == 1 ) then
      x1 = 0.0D+00
      x2 = 10.0D+00
    else if ( x_interval == 2 ) then
      x1 = -1.0D+00
      x2 = +1.0D+00
    else if ( x_interval == 3 ) then
      x1 = 0.5D+00
      x2 = 0.75D+00
    end if

    call cubic_value ( x1, f1, d1, s1, t1 )
    call cubic_value ( x2, f2, d2, s2, t2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '                                 Exact       Computed'
    write ( *, '(a)' ) '    J      A           B         Integral    Integral'
    write ( *, '(a)' ) ' '

    a = x1 - 1.0D+00

    do j = - 3, 12

      b = ( real ( 10 - j, kind = 8 ) * x1   &
          + real (      j, kind = 8 ) * x2 ) &
          / real ( 10,     kind = 8 )

      call cubic_integrate ( a, b, q_exact )

      call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, q_computed )

      write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) &
        j, a, b, q_exact, q_computed

    end do

  end do

  return
end
subroutine hermite_cubic_test04 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST04 tests HERMITE_CUBIC_SPLINE_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 51
  integer ( kind = 4 ), parameter :: nn = 11

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dn(nn)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fn(nn)
  integer ( kind = 4 ) i
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xn(nn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST04:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.'

  x1 = 0.0D+00
  x2 = 10.0D+00

  call r8vec_even ( nn, x1, x2, xn )

  fn(1:nn) = sin ( xn(1:nn) )
  dn(1:nn) = cos ( xn(1:nn) )

  call r8vec_even ( n, x1, x2, x )

  call hermite_cubic_spline_value ( nn, xn, fn, dn, n, x, f, d, s, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      X       F computed     F exact      Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    u = sin ( x(i) )
    v = abs ( f(i) - u )
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) i, x(i), f(i), u, v
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      X       D computed     D exact      Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    u = cos ( x(i) )
    v = abs ( d(i) - u )
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) i, x(i), d(i), u, v
  end do

  return
end
subroutine hermite_cubic_test05 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST05 tests HERMITE_CUBIC_TO_POWER_CUBIC
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 1

  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) d1
  real ( kind = 8 ) d1r
  real ( kind = 8 ) d2
  real ( kind = 8 ) d2r
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) f1
  real ( kind = 8 ) f1r
  real ( kind = 8 ) f2
  real ( kind = 8 ) f2r
  real ( kind = 8 ) fp
  integer ( kind = 4 ) j
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST05:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_TO_POWER_CUBIC converts the Hermite data'
  write ( *, '(a)' ) '  to the coefficients of the power form of the polynomial'
  write ( *, '(a)' ) '  POWER_CUBIC_TO_HERMITE_CUBIC converts the power form'
  write ( *, '(a)' ) '  to Hermite form'

  x1 = -1.0D+00
  x2 = +1.0D+00

  call cubic_value ( x1, f1, d1, s1, t1 )
  call cubic_value ( x2, f2, d2, s2, t2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Hermite data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) '  X1, F1, D1:', x1, f1, d1
  write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) '  X2, F2, D2:', x2, f2, d2

  call hermite_cubic_to_power_cubic ( x1, f1, d1, x2, f2, d2, c0, c1, c2, c3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Power form:'
  write ( *, '(a,f10.3,a,f10.3,a,f10.3,a,f10.3,a)' ) &
    '    p(x) = ', c0, ' + ', c1, ' * x + ', c2, ' * x^2 + ', c3, ' * x^3'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       F (Hermite)  F (power)'
  write ( *, '(a)' ) ' '

  do j = - 3, 12

    x(1) = ( real ( 10 - j, kind = 8 ) * x1   &
           + real (      j, kind = 8 ) * x2 ) &
           / real ( 10,     kind = 8 )

    call hermite_cubic_value ( x1, f1, d1, x2, f2, d2, 1, x, f, d, s, t )

    fp = c0 + x(1) * ( c1 + x(1) * ( c2 + x(1) * c3 ) )

    write ( *, '(2x,f10.4,2x,f10.4,2x,f10.4)' ) x, f, fp

  end do

  call power_cubic_to_hermite_cubic ( c0, c1, c2, c3, x1, x2, f1r, d1r, &
    f2r, d2r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use POWER_CUBIC_TO_HERMITE_CUBIC to recover the'
  write ( *, '(a)' ) '  original Hermite data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         Original   Recovered'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2x,f10.4,2x,f10.4)' ) '  F1:  ', f1, f1r
  write ( *, '(a,2x,f10.4,2x,f10.4)' ) '  D1:  ', d1, d1r
  write ( *, '(a,2x,f10.4,2x,f10.4)' ) '  F2:  ', f2, f2r
  write ( *, '(a,2x,f10.4,2x,f10.4)' ) '  D2:  ', d2, d2r

  return
end
subroutine hermite_cubic_test06 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST06 tests HERMITE_CUBIC_INTEGRATE using vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) b_hi
  real ( kind = 8 ) b_lo
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) q_computed
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST06:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic'
  write ( *, '(a)' ) '  polynomial from A to B.'
  write ( *, '(a)' ) '  Use A, B vectors for the calculation.'

  x1 = 0.0D+00
  x2 = 10.0D+00

  call cubic_value ( x1, f1, d1, s1, t1 )
  call cubic_value ( x2, f2, d2, s2, t2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                 Exact       Computed'
  write ( *, '(a)' ) '    J      A           B         Integral    Integral'
  write ( *, '(a)' ) ' '

  do i = -3, 12

    a = x1 - 1.0D+00
    b = ( real ( 10 - i, kind = 8 ) * x1 &
        + real (      i, kind = 8 ) * x2 ) &
        / real ( 10,     kind = 8 )

    call cubic_integrate ( a, b, q_exact )

    call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, q_computed )

    write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) &
      i, a, b, q_exact, q_computed

  end do

  return
end
subroutine hermite_cubic_test07 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST07 tests HERMITE_CUBIC_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) q_computed
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  integer ( kind = 4 ) x_interval
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST07:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_INTEGRAL integrates a Hermite cubic'
  write ( *, '(a)' ) '  polynomial over the definition interval [X1,X2].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                            Exact       Computed'
  write ( *, '(a)' ) '     X1          X2         Integral    Integral'
  write ( *, '(a)' ) ' '

  do x_interval = 1, 3

    if ( x_interval == 1 ) then
      x1 = 0.0D+00
      x2 = 10.0D+00
    else if ( x_interval == 2 ) then
      x1 = -1.0D+00
      x2 = +1.0D+00
    else if ( x_interval == 3 ) then
      x1 = 0.5D+00
      x2 = 0.75D+00
    end if

    call cubic_value ( x1, f1, d1, s1, t1 )
    call cubic_value ( x2, f2, d2, s2, t2 )

    call cubic_integrate ( x1, x2, q_exact )

    call hermite_cubic_integral ( x1, f1, d1, x2, f2, d2, q_computed )

    write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) &
      x1, x2, q_exact, q_computed

  end do

  return
end
subroutine hermite_cubic_test08 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST08 tests HERMITE_CUBIC_SPLINE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nn = 11

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) dn(nn)
  real ( kind = 8 ) fn(nn)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) q_computed
  real ( kind = 8 ) q_exact
  integer ( kind = 4 ) test
  real ( kind = 8 ) xn(nn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST08:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite'
  write ( *, '(a)' ) '  cubic spline over the definition interval [X1,XNN].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                            Exact       Computed'
  write ( *, '(a)' ) '     X1          XNN        Integral    Integral'
  write ( *, '(a)' ) ' '

  do test = 1, 3

    if ( test == 1 ) then

      a = 0.0D+00
      b = 1.0D+00

      call r8vec_even ( nn, a, b, xn )
      fn(1:nn) = xn(1:nn) * ( 4.0D+00 * xn(1:nn) - 1.0D+00 ) * ( xn(1:nn) - 1.0D+00 )
      dn(1:nn) = 1.0D+00 + xn(1:nn) * ( - 10.0D+00 + xn(1:nn) * 12.0D+00 )
      q_exact = &
        ( xn(nn) * xn(nn) * ( 0.5D+00 + xn(nn) * ( - ( 5.0D+00 / 3.0D+00 ) + xn(nn) ) ) ) &
      - ( xn(1)  * xn(1)  * ( 0.5D+00 + xn(1)  * ( - ( 5.0D+00 / 3.0D+00 ) + xn(1)  ) ) )
!
!  Use variable spacing.
!
    else if ( test == 2 ) then

      a = 0.0D+00
      b = 1.0D+00

      call r8vec_even ( nn, a, b, xn )
      xn(1:nn) = sqrt ( xn(1:nn) )
      fn(1:nn) = xn(1:nn) * ( 4.0D+00 * xn(1:nn) - 1.0D+00 ) * ( xn(1:nn) - 1.0D+00 )
      dn(1:nn) = 1.0D+00 + xn(1:nn) * ( - 10.0D+00 + xn(1:nn) * 12.0D+00 )
      q_exact = &
        ( xn(nn) * xn(nn) * ( 0.5D+00 + xn(nn) * ( - ( 5.0D+00 / 3.0D+00 ) + xn(nn) ) ) ) &
      - ( xn(1)  * xn(1)  * ( 0.5D+00 + xn(1)  * ( - ( 5.0D+00 / 3.0D+00 ) + xn(1)  ) ) )
!
!  Try a non-cubic.
!
    else if ( test == 3 ) then

      a = 0.0D+00
      b = pi

      call r8vec_even ( nn, a, b, xn )
      fn(1:nn) = sin ( xn(1:nn) )
      dn(1:nn) = cos ( xn(1:nn) )
      q_exact = - cos ( xn(nn) ) + cos ( xn(1) )

    end if

    call hermite_cubic_spline_integral ( nn, xn, fn, dn, q_computed )

    write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) &
      xn(1), xn(nn), q_exact, q_computed

  end do

  return
end
subroutine hermite_cubic_test09 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST09 tests HERMITE_CUBIC_SPLINE_INTEGRATE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nn = 11
  integer ( kind = 4 ), parameter :: n = 25

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) dn(n)
  real ( kind = 8 ) fn(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) q(n)
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) sn(n)
  real ( kind = 8 ) tn(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xn(nn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST09:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_SPLINE_INTEGRATE integrates a Hermite'
  write ( *, '(a)' ) '  cubic spline from A to B.'
!
!  Define the cubic spline.
!
  x1 = 0.0D+00
  x2 = 10.0D+00

  call r8vec_even ( nn, x1, x2, xn )

  do i = 1, nn
    call cubic_value ( xn(i), fn(i), dn(i), sn(i), tn(i) )
  end do

  a(1:n) = 2.5D+00
  call r8vec_even ( n, x1 - 1.0D+00, x2 + 1.0D+00, b )

  call hermite_cubic_spline_integrate ( nn, xn, fn, dn, n, a, b, q )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                 Exact       Computed'
  write ( *, '(a)' ) '    I      A           B         Integral    Integral'
  write ( *, '(a)' ) ' '

  do i = 1, n

    call cubic_integrate ( a(i), b(i), q_exact )

    write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) &
      i, a(i), b(i), q_exact, q(i)

  end do

  return
end
subroutine hermite_cubic_test10 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST10 tests HERMITE_CUBIC_SPLINE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nn = 11

  character ( len = 100 ) comment
  real ( kind = 8 ) dn(nn)
  real ( kind = 8 ) fn(nn)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) q_computed
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  real ( kind = 8 ) xn(nn)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST10:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite'
  write ( *, '(a)' ) '  cubic spline over the definition interval [X1,XNN].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If the subintervals are equally spaced, the derivative'
  write ( *, '(a)' ) '  information has no effect on the result, except for'
  write ( *, '(a)' ) '  the first and last values, DN(1) and DN(NN).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                            Exact       Computed'
  write ( *, '(a)' ) '     X1          XNN        Integral    Integral  Comment'
  write ( *, '(a)' ) ' '

  do test = 1, 5
!
!  Equal spacing.
!
    if ( test == 1 ) then

      call r8vec_even ( nn, 0.0D+00, pi, xn )
      fn(1:nn) = sin ( xn(1:nn) )
      dn(1:nn) = cos ( xn(1:nn) )
      q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
      comment = 'Equal spacing, correct DN'
!
!  Equal spacing, reset DN(2:NN-1) to random numbers.
!
    else if ( test == 2 ) then

      call r8vec_even ( nn, 0.0D+00, pi, xn )
      fn(1:nn) = sin ( xn(1:nn) )
      dn(1:nn) = cos ( xn(1:nn) )
      do i = 2, nn - 1
        dn(i) = 1000.0D+00 * r8_uniform_01 ( seed )
      end do
      q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
      comment = 'Equal spacing, DN(2:N-1) random'
!
!  Equal spacing, now reset all of DN to random numbers.
!
    else if ( test == 3 ) then

      call r8vec_even ( nn, 0.0D+00, pi, xn )
      fn(1:nn) = sin ( xn(1:nn) )
      do i = 1, nn
        dn(i) = 1000.0D+00 * r8_uniform_01 ( seed )
      end do
      q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
      comment = 'Equal spacing, DN(1:N) random'
!
!  Variable spacing, correct data.
!
    else if ( test == 4 ) then

      call r8vec_even ( nn, 0.0D+00, pi**2, xn )
      xn(1:nn) = sqrt ( xn(1:nn) )
      fn(1:nn) = sin ( xn(1:nn) )
      dn(1:nn) = cos ( xn(1:nn) )
      q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
      comment = 'Variable spacing, correct DN'
!
!  Variable spacing, change one entry in DN.
!
    else if ( test == 5 ) then

      call r8vec_even ( nn, 0.0D+00, pi**2, xn )
      xn(1:nn) = sqrt ( xn(1:nn) )
      fn(1:nn) = sin ( xn(1:nn) )
      dn(1:nn) = cos ( xn(1:nn) )
      dn( ( nn + 1 ) / 2 ) = 1000.0D+00 * r8_uniform_01 ( seed )
      q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
      comment = 'Variable spacing, a single internal DN randomized.'

    end if

    call hermite_cubic_spline_integral ( nn, xn, fn, dn, q_computed )

    write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6,2x,a)' ) &
      xn(1), xn(nn), q_exact, q_computed, trim ( comment )

  end do

  return
end
subroutine hermite_cubic_test11 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST11 tests HERMITE_CUBIC_LAGRANGE_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ) d(4,n)
  real ( kind = 8 ) f(4,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) s(4,n)
  real ( kind = 8 ) t(4,n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST11:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_LAGRANGE_VALUE evaluates the four'
  write ( *, '(a)' ) '  Lagrange basis functions associated with F1, D1,'
  write ( *, '(a)' ) '  F2 and D2 such that'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P(X) = F1 * LF1(X) + D1 * LD1(X)'
  write ( *, '(a)' ) '       + F2 * LF2(X) + D2 * LD2(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first, second and third derivatives of these four'
  write ( *, '(a)' ) '  Lagrange basis functions are also computed.'

  x1 = 1.0D+00
  x2 = 2.0D+00
  call r8vec_even ( n, 0.0D+00, 2.5D+00, x );

  call hermite_cubic_lagrange_value ( x1, x2, n, x, f, d, s, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Lagrange basis functions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I        X           LF1         LD1         LF2         LD2'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      j, x(j), f(1:4,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The derivative of the Lagrange basis functions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I        X           LF1         LD1         LF2         LD2'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      j, x(j), d(1:4,j)
  end do

  return
end
subroutine hermite_cubic_test12 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST12 tests HERMITE_CUBIC_LAGRANGE_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) q(4)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST12:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_LAGRANGE_INTEGRAL returns the integrals'
  write ( *, '(a)' ) '  of the four Lagrange basis functions associated'
  write ( *, '(a)' ) '  with F1, D1, F2 and D2 such that'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P(X) = F1 * LF1(X) + D1 * LD1(X)'
  write ( *, '(a)' ) '       + F2 * LF2(X) + D2 * LD2(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Lagrange basis function integrals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        X1          X2          LF1         LD1         LF2         LD2'
  write ( *, '(a)' ) ' '

  x2 = 1.0D+00

  do i = -6, 2
    x1 = real ( i, kind = 8 )
    call hermite_cubic_lagrange_integral ( x1, x2, q )
    write ( *, '(6(2x,f10.4))' ) x1, x2, q(1:4)
  end do

  return
end
subroutine hermite_cubic_test13 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST13 tests HERMITE_CUBIC_LAGRANGE_INTEGRATE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(4)
  real ( kind = 8 ) q(4)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' )  ' '
  write ( *, '(a)' )  'HERMITE_CUBIC_TEST13:'
  write ( *, '(a)' )  '  HERMITE_CUBIC_LAGRANGE_INTEGRATE integrates a Hermite cubic'
  write ( *, '(a)' )  '  Lagrange polynomial from A to B.'
  write ( *, '(a)' )  ' '
  write ( *, '(a)' )  '  Compute each result TWICE:'
  write ( *, '(a)' )  '  First row computed using HERMITE_CUBIC_INTEGRATE.'
  write ( *, '(a)' )  '  Second row computed using HERMITE_CUBIC_LAGRANGE_INTEGRATE.'

  x1 = 0.0D+00
  x2 = 10.0D+00

  write ( *, '(a)' )  ' '
  write ( *, '(a)' )  '        A           B           LF1         LD1         LF2         LD2'
  write ( *, '(a)' )  ' '

  a = x1 - 1.0D+00

  do j = -3, 12

    b = ( real ( 10 - j, kind = 8 ) * x1   &
        + real (      j, kind = 8 ) * x2 ) &
        / real ( 10,     kind = 8 )

    f1 = 1.0D+00
    d1 = 0.0D+00
    f2 = 0.0D+00
    d2 = 0.0D+00
    call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, p(1) )

    f1 = 0.0D+00
    d1 = 1.0D+00
    f2 = 0.0D+00
    d2 = 0.0D+00
    call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, p(2) )

    f1 = 0.0D+00
    d1 = 0.0D+00
    f2 = 1.0D+00
    d2 = 0.0D+00
    call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, p(3) )

    f1 = 0.0D+00
    d1 = 0.0D+00
    f2 = 0.0D+00
    d2 = 1.0D+00
    call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, p(4) )

    call hermite_cubic_lagrange_integrate ( x1, x2, a, b, q )

    write ( *, '(6(2x,f10.4))' ) a, b, p(1:4)
    write ( *, '(24x,4(2x,f10.4))' )       q(1:4)

  end do

  return
end
subroutine hermite_cubic_test14 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST14 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ) dn(n)
  real ( kind = 8 ) fn(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) q
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(2,n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST14:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule'
  write ( *, '(a)' ) '  for Hermite cubic splines.'

  do k = 1, 2

    write ( *, '(a)' ) ' '
    if ( k == 1 ) then
      write ( *, '(a)' ) '  Case 1: Random spacing'

      seed = 123456789

      call r8vec_uniform_01 ( n, seed, r )

      x(1) = r(1)
      do i = 2, n
        x(i) = x(i-1) + r(i)
      end do
    else if ( k == 2 ) then
      write ( *, '(a)' ) '  Case 2: Uniform spacing'
      write ( *, '(a)' ) '  F(2:N-1) have equal weight.'
      write ( *, '(a)' ) '  D(2:N-1) have zero weight.'
      do i = 1, n
        x(i) = real ( 10 + i - 1, kind = 8 ) / 20.0D+00
      end do
    end if

    call hermite_cubic_spline_quad_rule ( n, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   I   J        X         W                Q'
    write ( *, '(a)' ) ' '

    do i = 1, 2

      do j = 1, n

        fn(1:n) = 0.0D+00
        dn(1:n) = 0.0D+00

        if ( i == 1 ) then
          fn(j) = 1.0D+00
        else
          dn(j) = 1.0D+00
        end if

        call hermite_cubic_spline_integral ( n, x, fn, dn, q )

        write ( *, '(2x,i2,2x,i2,2x,f10.4, 2x, g14.6,2x,g14.6)' ) &
          i, j, x(j), w(i,j), q

      end do

    end do

  end do

  return
end
subroutine hermite_cubic_test15 ( )

!*****************************************************************************80
!
!! HERMITE_CUBIC_TEST15 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ) dn(n)
  real ( kind = 8 ) fn(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) q
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) w(2,n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_CUBIC_TEST15:'
  write ( *, '(a)' ) '  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule'
  write ( *, '(a)' ) '  for Hermite cubic splines.'

  seed = 123456789

  call r8vec_uniform_01 ( n, seed, r )

  x(1) = r(1)
  do j = 2, n
    x(j) = x(j-1) + r(j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random spacing'
  write ( *, '(a,i8)' ) '  Number of points N = ', n
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  Interval = [', x(1), ',', x(n), ']'

  call hermite_cubic_spline_quad_rule ( n, x, w )

  do j = 1, n
    call cubic_value ( x(j), fn(j), dn(j), s, t )
  end do

  q = dot_product ( w(1,1:n), fn(1:n) ) + dot_product ( w(2,1:n), dn(1:n) )

  call cubic_integrate ( x(1), x(n), q_exact )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Q         = ', q
  write ( *, '(a,g14.6)' ) '  Q (exact) = ', q_exact

  return
end
function cubic_antiderivative ( x )

!*****************************************************************************80
!
!! CUBIC_ANTIDERIVATIVE evaluates the antiderivative function of a cubic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) CUBIC_ANTIDERIVATIVE, the value.
!
  implicit none

  real ( kind = 8 ) cubic_antiderivative
  real ( kind = 8 ) x

  cubic_antiderivative = x * x * ( 5.0D+00 + x * ( - 7.0D+00 / 3.0D+00 &
    + x * 1.0D+00 / 4.0D+00 ) )

  return
end
subroutine cubic_integrate ( a, b, q )

!*****************************************************************************80
!
!! CUBIC_INTEGRATE integrates the cubic from A to B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the integration interval.
!
!    Output, real ( kind = 8 ) Q, the integral from A to B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cubic_antiderivative
  real ( kind = 8 ) q

  q = cubic_antiderivative ( b ) - cubic_antiderivative ( a )

  return
end
subroutine cubic_value ( x, f, d, s, t )

!*****************************************************************************80
!
!! CUBIC_VALUE evaluates a cubic function.
!
!  Discussion:
!
!    f(x) =   x^3 -  7 x^2 + 10 x
!    d(x) = 3 x^2 - 14 x   + 10
!    s(x) = 6 x   - 14
!    t(x) = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) F, D, S, T, the value and first three
!    derivatives of the cubic function.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  f = 0.0D+00 + x * ( 10.0D+00 + x * (  - 7.0D+00 + x * 1.0D+00 ) )
  d =                 10.0D+00 + x * ( - 14.0D+00 + x * 3.0D+00 )
  s =                                  - 14.0D+00 + x * 6.0D+00
  t =                                                   6.0D+00

  return
end
