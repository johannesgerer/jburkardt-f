program main

!*****************************************************************************80
!
!! MAIN is the main program for DIVDIF_PRB.
!
!  Discussion:
!
!    DIVDIF_PRB tests the DIVDIF library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIVDIF_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DIVDIF library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test085 ( )
  call test09 ( )
  call test095 ( )

  call test10 ( )
  call test105 ( )
  call test11 ( )
  call test115 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIVDIF_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DIF_APPEND, DIF_ANTIDERIV, DIF_DERIV_TABLE, DIF_SHIFT_ZERO, DIF_VAL;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 10

  real ( kind = 8 ) diftab(maxtab)
  real ( kind = 8 ) diftab2(maxtab)
  real ( kind = 8 ) diftab3(maxtab)
  integer ( kind = 4 ) ntab
  integer ( kind = 4 ) ntab2
  integer ( kind = 4 ) ntab3
  real ( kind = 8 ) xtab(maxtab)
  real ( kind = 8 ) xtab2(maxtab)
  real ( kind = 8 ) xtab3(maxtab)
  real ( kind = 8 ) xval
  real ( kind = 8 ) ytab(maxtab)
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a divided difference polynomial:'
  write ( *, '(a)' ) '  DATA_TO_DIF_DISPLAY sets up a difference table'
  write ( *, '(a)' ) '  and displays intermediate calculations;'
  write ( *, '(a)' ) '  DIF_APPEND appends a new data point;'
  write ( *, '(a)' ) '  DIF_ANTIDERIV computes the antiderivative;'
  write ( *, '(a)' ) '  DIF_DERIV_TABLE computes the derivative;'
  write ( *, '(a)' ) '  DIF_SHIFT_ZERO shifts all the abscissas to 0;'
  write ( *, '(a)' ) '  DIF_VAL evaluates at a point.'
  write ( *, '(a)' ) ' '
!
!  Set XTAB, YTAB to X, X^2.
!
  ntab = 4
  call r8vec_indicator ( ntab, xtab )
  ytab(1:ntab) = xtab(1:ntab)**2

  call data_to_dif_display ( ntab, xtab, ytab, diftab )

  call dif_print ( ntab, xtab, diftab, '  The divided difference polynomial:' )
!
!  Append (5,25) to the table.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Append the data (5,25) to the table.'
  write ( *, '(a)' ) ' '

  xval = 5.0D+00
  yval = 25.0D+00

  call dif_append ( ntab, xtab, diftab, xval, yval, ntab, xtab, diftab )

  call dif_print ( ntab, xtab, diftab, &
    '  The augmented divided difference polynomial:' )
!
!  Evaluate the polynomial at 2.5.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Evaluate the table at a point.'
  write ( *, '(a)' ) ' '

  xval = 2.5D+00

  call dif_val ( ntab, xtab, diftab, xval, yval )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  P( ', xval, ' ) = ', yval
!
!  Shift the base to zero.
!
  call dif_shift_zero ( ntab, xtab, diftab )

  call dif_print ( ntab, xtab, diftab, '  The table, rebased at 0:' )
!
!  Compute a divided difference table for the derivative.
!
  call dif_deriv_table ( ntab, xtab, diftab, ntab2, xtab2, diftab2 )

  call dif_print ( ntab2, xtab2, diftab2, '  The derivative:' )

  call dif_val ( ntab2, xtab2, diftab2, xval, yval )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  P''( ', xval, ' ) = ', yval
!
!  Compute the antiderivative.
!
  call dif_antideriv ( ntab, xtab, diftab, ntab3, xtab3, diftab3 )

  call dif_print ( ntab3, xtab3, diftab3, '  The antiderivative:' )

  call dif_val ( ntab3, xtab3, diftab3, xval, yval )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Ant(P)( ', xval, ' ) = ', yval

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DATA_TO_DIF and DIF_VAL.
!
!  Discussion:
!
!    This test demonstrates how a divided difference table can be generated,
!    and then used to approximate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 16

  real ( kind = 8 ) diftab(maxtab)
  real ( kind = 8 ) err
  real ( kind = 8 ) f02
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  integer ( kind = 4 ), parameter :: ntest = 1001
  real ( kind = 8 ) xtab(maxtab)
  real ( kind = 8 ) xtest
  real ( kind = 8 ) ytab(maxtab)
  real ( kind = 8 ) ytest

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  DATA_TO_DIF takes a set of (X,F(X)) data'
  write ( *, '(a)' ) '  and computes a divided difference table.'
  write ( *, '(a)' ) '  DIF_VAL evaluates the corresponding polynomial'
  write ( *, '(a)' ) '  interpolant at an arbitrary point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By increasing the number of data points,'
  write ( *, '(a)' ) '  the approximation should improve for a while.'
  write ( *, '(a)' ) '  However, our function is non-differentiable'
  write ( *, '(a)' ) '  at one point, so the approximation begins to'
  write ( *, '(a)' ) '  misbehave rapidly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our interval is [-1,1].'
  write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We estimate the interpolation error using'
  write ( *, '(a)' ) '  1001 equally spaced sample points.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Order  Interpolation Error'
  write ( *, '(a)' ) ' '

  do ntab = 1, maxtab

    call r8vec_even ( ntab, -1.0D+00, +1.0D+00, xtab )

    do j = 1, ntab
      ytab(j) = f02 ( xtab(j) )
    end do

    call data_to_dif ( ntab, xtab, ytab, diftab )

    do j = 1, ntest
      call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
      call dif_val ( ntab, xtab, diftab, xtest, ytest )
      err = err + ( ytest - f02 ( xtest ) )**2
    end do

    err = 2.0D+00 * sqrt ( err ) / real ( ntab, kind = 8 )

    write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

  end do

  return
end
function f02 ( x )

!*****************************************************************************80
!
!! F02 evaluates the function F(X) = |X| + X/2 - X^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f02
  real ( kind = 8 ) x

  f02 = abs ( x ) + x / 2.0D+00 - x * x

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DIF_BASIS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 5

  real ( kind = 8 ) diftab(ntab,ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: nstep = 9
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  DIF_BASIS computes Lagrange basis polynomials'
  write ( *, '(a)' ) '  in difference form.'
  write ( *, '(a)' ) ' '
!
!  Set the base points.
!
  call r8vec_indicator ( ntab, xtab )

  call r8vec_print ( ntab, xtab, '  The base points:' )
!
!  Get the difference tables for the basis polynomials and print them.
!
  call dif_basis ( ntab, xtab, diftab )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The table of difference vectors defining the basis'
  write ( *, '(a)' ) '  polynomials.  Each column represents a polynomial.'
  write ( *, '(a)' ) ' '
  do i = 1, ntab
    write ( *, '(2x,5g14.6)' ) diftab(i,1:ntab)
  end do
!
!  Evaluate basis polynomial 3 at a set of points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Evaluate basis polynomial #3 at a set of points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X        Y'
  write ( *, '(a)' ) ' '
  xhi = real ( ntab, kind = 8 )
  xlo = 1.0D+00

  do i = 1, nstep

    xval = ( real ( nstep - i,     kind = 8 ) * xlo   &
           + real (         i - 1, kind = 8 ) * xhi ) &
           / real ( nstep     - 1, kind = 8 )

    call dif_val ( ntab, xtab, diftab(1,3), xval, yval )

    write ( *, '(2g14.6)' ) xval, yval

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DIF_ROOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 10

  real ( kind = 8 ) abserr
  real ( kind = 8 ), external :: cubic
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) maxstp
  real ( kind = 8 ) relerr
  real ( kind = 8 ) xroot
  real ( kind = 8 ) xtry1
  real ( kind = 8 ) xtry2
!
!  Seek a root of the function F(X) = (X+3)*(X+1)*(X-1)
!  given starting estimates of 0.5 and 2.0.
!
  abserr = 0.00001D+00
  iprint = 1
  maxstp = 15
  relerr = 0.00001D+00
  xtry1 = 0.5D+00
  xtry2 = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DIF_ROOT seeks a zero of F(x).'
  write ( *, '(a)' ) '  F(X) = (X+3)*(X+1)*(X-1)'
  write ( *, '(a)' ) ' '

  call dif_root ( abserr, cubic, iprint, maxstp, maxtab, &
    relerr, xroot, xtry1, xtry2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated root = ', xroot
  write ( *, '(a,g14.6)' ) '  F(X) = ', cubic ( xroot )

  return
end
function cubic ( x )

!*****************************************************************************80
!
!! CUBIC is the function whose root is sought.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cubic
  real ( kind = 8 ) x

  cubic = ( x + 3.0D+00 ) * ( x + 1.0D+00 ) * ( x - 1.0D+00 )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests DIF_TO_R8POLY and DIF_SHIFT_ZERO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 10

  real ( kind = 8 ) c(maxtab)
  real ( kind = 8 ) diftab1(maxtab)
  real ( kind = 8 ) diftab2(maxtab)
  integer ( kind = 4 ) ntab
  real ( kind = 8 ) xtab1(maxtab)
  real ( kind = 8 ) xtab2(maxtab)
  real ( kind = 8 ) ytab1(maxtab)
  real ( kind = 8 ) ytab2(maxtab)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  DIF_TO_R8POLY converts a difference table to a'
  write ( *, '(a)' ) '  polynomial;'
  write ( *, '(a)' ) '  DIF_SHIFT_ZERO shifts a divided difference '
  write ( *, '(a)' ) '  table to all zero abscissas;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  These are equivalent operations!'
  write ( *, '(a)' ) ' '
!
!  Set XTAB, YTAB to X, F(X)
!
  ntab = 4

  call r8vec_indicator ( ntab, xtab1 )
  ytab1(1:ntab) = xtab1(1:ntab)**3 - 2.0D+00 * xtab1(1:ntab)**2 &
    + 3.0D+00 * xtab1(1:ntab) - 4.0D+00

  call r8vec_indicator ( ntab, xtab2 )
  ytab2(1:ntab) = xtab2(1:ntab)**3 - 2.0D+00 * xtab2(1:ntab)**2 &
    + 3.0D+00 * xtab2(1:ntab) - 4.0D+00
!
!  Compute and display the finite difference table.
!
  call data_to_dif_display ( ntab, xtab1, ytab1, diftab1 )

  call data_to_dif_display ( ntab, xtab2, ytab2, diftab2 )
!
!  Examine corresponding polynomial.
!
  call dif_print ( ntab, xtab1, diftab1, &
    '  The divided difference polynomial:' )
!
!  Shift to zero.
!
  call dif_shift_zero ( ntab, xtab1, diftab1 )

  call r8poly_print ( ntab, diftab1, '  Using DIF_SHIFT_ZERO' )

  call dif_print ( ntab, xtab1, diftab1, &
    '  The divided difference polynomial after DIF_SHIFT_ZERO:' )
!
!  Shift to zero.
!
  call dif_to_r8poly ( ntab, xtab2, diftab2, c )

  call r8poly_print ( ntab, c, '  Using DIF_TO_R8POLY' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests R8POLY_*.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  real ( kind = 8 ) poly_cof(n)
  real ( kind = 8 ) poly_cof2(n+1)
  real ( kind = 8 ) poly_cof3(n-1)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval
  real ( kind = 8 ) yval2
  real ( kind = 8 ) yval3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  R8POLY_ANT_COF computes the coefficients of the'
  write ( *, '(a)' ) '  antiderivative of a polynomial;'
  write ( *, '(a)' ) '  R8POLY_ANT_VAL evaluates the antiderivative of'
  write ( *, '(a)' ) '  a polynomial;'
  write ( *, '(a)' ) '  R8POLY_DER_COF computes the coefficients of the'
  write ( *, '(a)' ) '  derivative of a polynomial;'
  write ( *, '(a)' ) '  R8POLY_DER_VAL evaluates the derivative of'
  write ( *, '(a)' ) '  a polynomial;'
  write ( *, '(a)' ) '  R8POLY_PRINT prints a polynomial;'
  write ( *, '(a)' ) '  R8POLY_VAL evaluates a polynomial.'

  do i = 1, n
    poly_cof(i) = real ( i, kind = 8 )
  end do

  call r8poly_print ( n, poly_cof, '  Our initial polynomial:' )

  call r8poly_ant_cof ( n, poly_cof, poly_cof2 )

  call r8poly_print ( n+1, poly_cof2, '  The antiderivative polynomial:' )

  call r8poly_der_cof ( n, poly_cof, poly_cof3 )

  call r8poly_print ( n-1, poly_cof3, '  The derivative polynomial:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Evaluate the polynomial, antiderivative and'
  write ( *, '(a)' ) '  derivative, using only the original polynomial'
  write ( *, '(a)' ) '  coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X             P(X)         Anti_P(X)     P''(X)'
  write ( *, '(a)' ) ' '

  do i = 0, 2

    xval = real ( i, kind = 8 )

    call r8poly_val_horner ( n, poly_cof, xval, yval )

    call r8poly_ant_val ( n, poly_cof, xval, yval2 )

    call r8poly_der_val ( n, poly_cof, xval, yval3 )

    write ( *, '(4g14.6)' ) xval, yval, yval2, yval3

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests R8POLY_BASIS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 5

  real ( kind = 8 ) polcof(ntab,ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: nstep = 9
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  R8POLY_BASIS computes Lagrange basis polynomials'
  write ( *, '(a)' ) '  in standard form.'
  write ( *, '(a)' ) ' '
!
!  Set the base points.
!
  call r8vec_indicator ( ntab, xtab )
!
!  Get the difference tables for the basis polynomials and print them.
!
  call r8poly_basis ( ntab, xtab, polcof )

  do i = 1, ntab
    write ( *, '(2x,5g14.6)' ) polcof(i,1:ntab)
  end do
!
!  Print basis polynomial 3 in polynomial form.
!
  call r8poly_print ( ntab, polcof(1,3), &
    '  Basis polynomial 3 in standard form:' )
!
!  Evaluate basis polynoimial 3 at a set of points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Evaluate basis polynomial 3 at a set of points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X        Y'
  write ( *, '(a)' ) ' '
  xhi = real ( ntab, kind = 8 )
  xlo = 1.0D+00

  do i = 1, nstep

    xval = ( real ( nstep - i,     kind = 8 ) * xlo   &
           + real (         i - 1, kind = 8 ) * xhi ) &
           / real ( nstep     - 1, kind = 8 )

    call r8poly_val_horner ( ntab, polcof(1,3), xval, yval )

    write ( *, '(2x,2g14.6)' ) xval, yval

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests R8POLY_SHIFT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) i
  real ( kind = 8 ) poly_cof(n)
  real ( kind = 8 ) scale
  real ( kind = 8 ) shift

  scale = 2.0D+00
  shift = +3.0D+00
  poly_cof(1:3) = (/ +6.0D+00, -1.0D+00,  2.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  R8POLY_SHIFT shifts polynomial coefficients.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Polynomial coefficients for argument X'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i5,g14.6)' ) i, poly_cof(i)
  end do

  call r8poly_shift ( scale, shift, n, poly_cof )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  SCALE = ', scale
  write ( *, '(a,g14.6)' ) '  SHIFT = ', shift
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Polynomial coefficients for argument '
  write ( *, '(a)' ) '    Z = SCALE * X + SHIFT'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i5,g14.6)' ) i, poly_cof(i)
  end do

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests LAGRANGE_VAL;
!
!  Discussion:
!
!    This test demonstrates how a divided difference table can be generated,
!    and then used to approximate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 16

  real ( kind = 8 ) err
  real ( kind = 8 ) f02
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  integer ( kind = 4 ), parameter :: ntest = 1001
  real ( kind = 8 ) xtab(maxtab)
  real ( kind = 8 ) xtest
  real ( kind = 8 ) ytab(maxtab)
  real ( kind = 8 ) ytest

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  LAGRANGE_VAL uses naive Lagrange interpolation'
  write ( *, '(a)' ) '  to compute the polynomial interpolant to data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By increasing the number of data points,'
  write ( *, '(a)' ) '  the approximation should improve for a while.'
  write ( *, '(a)' ) '  However, our function is non-differentiable'
  write ( *, '(a)' ) '  at one point, so the approximation begins to'
  write ( *, '(a)' ) '  misbehave rapidly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our interval is [-1,1].'
  write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We estimate the interpolation error using'
  write ( *, '(a)' ) '  1001 equally spaced sample points.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Order  Interpolation Error'
  write ( *, '(a)' ) ' '

  do ntab = 1, maxtab

    call r8vec_even ( ntab, -1.0D+00, +1.0D+00, xtab )
!
!  Now set up some data.
!
    do j = 1, ntab
      ytab(j) = f02 ( xtab(j) )
    end do

    err = 0.0D+00

    do j = 1, ntest
      call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
      call lagrange_val ( ntab, xtab, ytab, xtest, ytest )
      err = err + ( ytest - f02 ( xtest ) )**2
    end do

    err = 2.0D+00 * sqrt ( err ) / real ( ntab, kind = 8 )

    write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests LAGRANGE_RULE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 8

  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  real ( kind = 8 ) wtab(maxtab)
  real ( kind = 8 ) xtab(maxtab)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  LAGRANGE_RULE computes Lagrange interpolation formulas;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lagrange Interpolation Rules on [-1,1]'
  write ( *, '(a)' ) '  using equally spaced abscissas.'

  do ntab = 1, maxtab

    call r8vec_even ( ntab, -1.0D+00, +1.0D+00, xtab )

    call lagrange_rule ( ntab, xtab, wtab )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      Abscissa       Weight'
    write ( *, '(a)' ) ' '
    do j = 1, ntab
      write ( *, '(2x,i3,g14.6,2x,g14.6)' ) j, xtab(j), wtab(j)
    end do

  end do

  return
end
subroutine test095 ( )

!*****************************************************************************80
!
!! TEST095 tests LAGRANGE_RULE and LAGRANGE_SUM.
!
!  Discussion:
!
!    This test demonstrates how a divided difference table can be generated,
!    and then used to approximate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 16

  real ( kind = 8 ) err
  real ( kind = 8 ) f02
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  integer ( kind = 4 ), parameter :: ntest = 1001
  real ( kind = 8 ) wtab(maxtab)
  real ( kind = 8 ) xtab(maxtab)
  real ( kind = 8 ) xtest
  real ( kind = 8 ) ytab(maxtab)
  real ( kind = 8 ) ytest

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST095'
  write ( *, '(a)' ) '  LAGRANGE_RULE sets the weights for a Lagrange rule.'
  write ( *, '(a)' ) '  LAGRANGE_SUM uses the rule to compute the Lagrange'
  write ( *, '(a)' ) '  interpolant to data at a given point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this test, the data abscissas are '
  write ( *, '(a)' ) '  equally spaced.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By increasing the number of data points,'
  write ( *, '(a)' ) '  the approximation should improve for a while.'
  write ( *, '(a)' ) '  However, our function is non-differentiable'
  write ( *, '(a)' ) '  at one point, so the approximation begins to'
  write ( *, '(a)' ) '  misbehave rapidly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our interval is [-1,1].'
  write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We estimate the interpolation error using'
  write ( *, '(a)' ) '  1001 equally spaced sample points.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Order  Interpolation Error'
  write ( *, '(a)' ) ' '

  do ntab = 1, maxtab

    call r8vec_even ( ntab, -1.0D+00, +1.0D+00, xtab )

    call lagrange_rule ( ntab, xtab, wtab )
!
!  Now set up some data.
!
    do j = 1, ntab
      ytab(j) = f02 ( xtab(j) )
    end do

    err = 0.0D+00
    do j = 1, ntest
      call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
      call lagrange_sum ( ntab, xtab, wtab, ytab, xtest, ytest )
      err = err + ( ytest - f02 ( xtest ) )**2
    end do

    err = 2.0D+00 * sqrt ( err ) / real ( ntab, kind = 8 )

    write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests LAGRANGE_RULE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 8

  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  real ( kind = 8 ) wtab(maxtab)
  real ( kind = 8 ) xtab(maxtab)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  LAGRANGE_RULE computes Lagrange interpolation formulas;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lagrange Interpolation Rules on [-1,1]'
  write ( *, '(a)' ) '  using Chebyshev T abscissas.'

  do ntab = 1, maxtab

    call cheby_t_zero ( ntab, xtab )

    call lagrange_rule ( ntab, xtab, wtab )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      Abscissa       Weight'
    write ( *, '(a)' ) ' '
    do j = 1, ntab
      write ( *, '(2x,i3,g14.6,2x,g14.6)' ) j, xtab(j), wtab(j)
    end do

  end do

  return
end
subroutine test105 ( )

!*****************************************************************************80
!
!! TEST105 tests LAGRANGE_RULE and LAGRANGE_SUM.
!
!  Discussion:
!
!    This test demonstrates how a divided difference table can be generated,
!    and then used to approximate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 16

  real ( kind = 8 ) err
  real ( kind = 8 ) f02
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  integer ( kind = 4 ), parameter :: ntest = 1001
  real ( kind = 8 ) wtab(maxtab)
  real ( kind = 8 ) xtab(maxtab)
  real ( kind = 8 ) xtest
  real ( kind = 8 ) ytab(maxtab)
  real ( kind = 8 ) ytest

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) '  LAGRANGE_RULE sets the weights for a Lagrange rule.'
  write ( *, '(a)' ) '  LAGRANGE_SUM uses the rule to compute the Lagrange'
  write ( *, '(a)' ) '  interpolant to data at a given point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this test, the data abscissas are the'
  write ( *, '(a)' ) '  zeroes of the Chebyshev T polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By increasing the number of data points,'
  write ( *, '(a)' ) '  the approximation should improve for a while.'
  write ( *, '(a)' ) '  However, our function is non-differentiable'
  write ( *, '(a)' ) '  at one point, so the approximation begins to'
  write ( *, '(a)' ) '  misbehave rapidly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our interval is [-1,1].'
  write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We estimate the interpolation error using'
  write ( *, '(a)' ) '  1001 equally spaced sample points.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Order  Interpolation Error'
  write ( *, '(a)' ) ' '

  do ntab = 1, maxtab

    call cheby_t_zero ( ntab, xtab )

    call lagrange_rule ( ntab, xtab, wtab )
!
!  Now set up some data.
!
    do j = 1, ntab
      ytab(j) = f02 ( xtab(j) )
    end do

    err = 0.0D+00
    do j = 1, ntest
      call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
      call lagrange_sum ( ntab, xtab, wtab, ytab, xtest, ytest )
      err = err + ( ytest - f02 ( xtest ) )**2
    end do

    err = 2.0D+00 * sqrt ( err ) / real ( ntab, kind = 8 )

    write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests LAGRANGE_RULE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 8

  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  real ( kind = 8 ) wtab(maxtab)
  real ( kind = 8 ) xtab(maxtab)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  LAGRANGE_RULE computes Lagrange interpolation formulas;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lagrange Interpolation Rules on [-1,1]'
  write ( *, '(a)' ) '  using Chebyshev U abscissas.'

  do ntab = 1, maxtab

    call cheby_u_zero ( ntab, xtab )

    call lagrange_rule ( ntab, xtab, wtab )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      Abscissa       Weight'
    write ( *, '(a)' ) ' '
    do j = 1, ntab
      write ( *, '(2x,i3,g14.6,2x,g14.6)' ) j, xtab(j), wtab(j)
    end do

  end do

  return
end
subroutine test115 ( )

!*****************************************************************************80
!
!! TEST115 tests LAGRANGE_RULE and LAGRANGE_SUM.
!
!  Discussion:
!
!    This test demonstrates how a divided difference table can be generated,
!    and then used to approximate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 16

  real ( kind = 8 ) err
  real ( kind = 8 ) f02
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  integer ( kind = 4 ), parameter :: ntest = 1001
  real ( kind = 8 ) wtab(maxtab)
  real ( kind = 8 ) xtab(maxtab)
  real ( kind = 8 ) xtest
  real ( kind = 8 ) ytab(maxtab)
  real ( kind = 8 ) ytest

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115'
  write ( *, '(a)' ) '  LAGRANGE_RULE sets the weights for a Lagrange rule.'
  write ( *, '(a)' ) '  LAGRANGE_SUM uses the rule to compute the Lagrange'
  write ( *, '(a)' ) '  interpolant to data at a given point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this test, the data abscissas are the'
  write ( *, '(a)' ) '  zeroes of the Chebyshev U polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By increasing the number of data points,'
  write ( *, '(a)' ) '  the approximation should improve for a while.'
  write ( *, '(a)' ) '  However, our function is non-differentiable'
  write ( *, '(a)' ) '  at one point, so the approximation begins to'
  write ( *, '(a)' ) '  misbehave rapidly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our interval is [-1,1].'
  write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We estimate the interpolation error using'
  write ( *, '(a)' ) '  1001 equally spaced sample points.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Order  Interpolation Error'
  write ( *, '(a)' ) ' '

  do ntab = 1, maxtab

    call cheby_u_zero ( ntab, xtab )

    call lagrange_rule ( ntab, xtab, wtab )
!
!  Now set up some data.
!
    do j = 1, ntab
      ytab(j) = f02 ( xtab(j) )
    end do

    err = 0.0D+00
    do j = 1, ntest
      call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
      call lagrange_sum ( ntab, xtab, wtab, ytab, xtest, ytest )
      err = err + ( ytest - f02 ( xtest ) )**2
    end do

    err = 2.0D+00 * sqrt ( err ) / real ( ntab, kind = 8 )

    write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests NCC_RULE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: norder = 8

  integer ( kind = 4 ) i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  NCC_RULE computes closed Newton Cotes formulas'
  write ( *, '(a)' ) '  for quadrature (approximate integration).'
  write ( *, '(a)') ' '

  call ncc_rule ( norder, xtab, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Newton-Cotes Closed Quadrature Rule:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Abscissa       Weight'
  write ( *, '(a)' ) ' '

  do i = 1, norder
    write ( *, '(2x,i3,2g14.6)' ) i, xtab(i), weight(i)
  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests NCO_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: norder = 8

  integer ( kind = 4 ) i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  NCO_RULE computes open Newton Cotes formulas'
  write ( *, '(a)' ) '  for quadrature (approximate integration).'
  write ( *, '(a)' ) ' '

  call nco_rule ( norder, xtab, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Newton-Cotes Open Quadrature Rule:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Abscissa       Weight'
  write ( *, '(a)' ) ' '

  do i = 1, norder
    write ( *, '(2x,i3,2g14.6)' ) i, xtab(i), weight(i)
  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests ROOTS_TO_DIF and DIF_TO_R8POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxroots = 4

  real ( kind = 8 ) c(maxroots+1)
  real ( kind = 8 ) diftab(maxroots+1)
  integer ( kind = 4 ) nroots
  integer ( kind = 4 ) ntab
  real ( kind = 8 ) roots(maxroots)
  real ( kind = 8 ) xtab(maxroots+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  ROOTS_TO_DIF computes the divided difference'
  write ( *, '(a)' ) '  polynomial with given roots;'
  write ( *, '(a)' ) '  DIF_TO_R8POLY converts it to a standard form'
  write ( *, '(a)' ) '  polynomial.'
  write ( *, '(a)' ) ' '

  nroots = 1
  roots(1) = 3.0D+00
  call r8vec_print ( nroots, roots, '  The roots:' )

  call roots_to_dif ( nroots, roots, ntab, xtab, diftab )

  call dif_to_r8poly ( ntab, xtab, diftab, c )

  call r8poly_print ( ntab, c, '  The polynomial:' )

  nroots = 2
  roots(1:2) = (/ 3.0D+00, 1.0D+00 /)
  call r8vec_print ( nroots, roots, '  The roots:' )

  call roots_to_dif ( nroots, roots, ntab, xtab, diftab )

  call dif_to_r8poly ( ntab, xtab, diftab, c )

  call r8poly_print ( ntab, c, '  The polynomial:' )

  nroots = 3
  roots(1:3) = (/ 3.0D+00, 1.0D+00, 2.0D+00 /)
  call r8vec_print ( nroots, roots, '  The roots:' )

  call roots_to_dif ( nroots, roots, ntab, xtab, diftab )

  call dif_to_r8poly ( ntab, xtab, diftab, c )

  call r8poly_print ( ntab, c, '  The polynomial:' )

  nroots = 4
  roots(1:4) = (/ 3.0D+00, 1.0D+00, 2.0D+00, 4.0D+00 /)
  call r8vec_print ( nroots, roots, '  The roots:' )

  call roots_to_dif ( nroots, roots, ntab, xtab, diftab )

  call dif_to_r8poly ( ntab, xtab, diftab, c )

  call r8poly_print ( ntab, c, '  The polynomial:' )

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests ROOTS_TO_R8POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxroot = 5

  real ( kind = 8 ) c(maxroot+1)
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nroots
  real ( kind = 8 ) roots(maxroot)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  ROOTS_TO_R8POLY computes polynomial coefficients'
  write ( *, '(a)' ) '  from roots.'
  write ( *, '(a)' ) ' '

  nroots = 1
  roots(1) = 3.0D+00
  call r8vec_print ( nroots, roots, '  The roots:' )

  call roots_to_r8poly ( nroots, roots, nc, c )

  call r8poly_print ( nc, c, '  The polynomial:' )

  nroots = 2
  roots(1:2) = (/ 3.0D+00, 1.0D+00 /)
  call r8vec_print ( nroots, roots, '  The roots:' )

  call roots_to_r8poly ( nroots, roots, nc, c )

  call r8poly_print ( nc, c, '  The polynomial:' )

  nroots = 3
  roots(1:3) = (/ 3.0D+00, 1.0D+00, 2.0D+00 /)
  call r8vec_print ( nroots, roots, '  The roots:' )

  call roots_to_r8poly ( nroots, roots, nc, c )

  call r8poly_print ( nc, c, '  The polynomial:' )

  nroots = 4
  roots(1:4) = (/ 3.0D+00, 1.0D+00, 2.0D+00, 4.0D+00 /)
  call r8vec_print ( nroots, roots, '  The roots:' )

  call roots_to_r8poly ( nroots, roots, nc, c )

  call r8poly_print ( nc, c, '  The polynomial:' )

  return
end
