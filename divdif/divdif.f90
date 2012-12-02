subroutine cheby_t_zero ( n, z )

!*****************************************************************************80
!
!! CHEBY_T_ZERO returns zeroes of the Chebyshev polynomial T(N)(X).
!
!  Discussion:
!
!    The I-th zero of T(N)(X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) Z(N), the zeroes of T(N)(X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) z(n)

  do i = 1, n
    angle = real ( 2 * i - 1, kind = 8 ) * pi &
          / real ( 2 * n,     kind = 8 )
    z(i) = cos ( angle )
  end do

  return
end
subroutine cheby_u_zero ( n, z )

!*****************************************************************************80
!
!! CHEBY_U_ZERO returns zeroes of the Chebyshev polynomial U(N)(X).
!
!  Discussion:
!
!    The I-th zero of U(N)(X) is cos((I-1)*PI/(N-1)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) Z(N), the zeroes of U(N)(X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) z(n)

  do i = 1, n
    angle = real ( i,     kind = 8 ) * pi &
          / real ( n + 1, kind = 8 )
    z(i) = cos ( angle )
  end do

  return
end
subroutine data_to_dif ( ntab, xtab, ytab, diftab )

!*****************************************************************************80
!
!! DATA_TO_DIF sets up a divided difference table from raw data.
!
!  Discussion:
!
!    Space can be saved by using a single array for both the DIFTAB and
!    YTAB dummy parameters.  In that case, the divided difference table will
!    overwrite the Y data without interfering with the computation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of pairs of points
!    (XTAB(I),YTAB(I)) which are to be used as data.  The
!    number of entries to be used in DIFTAB, XTAB and YTAB.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values at which data was taken.
!    These values must be distinct.
!
!    Input, real ( kind = 8 ) YTAB(NTAB), the corresponding Y values.
!
!    Output, real ( kind = 8 ) DIFTAB(NTAB), the divided difference coefficients
!    corresponding to the input (XTAB,YTAB).
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) diftab(ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r8vec_distinct
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) ytab(ntab)

  if ( .not. r8vec_distinct ( ntab, xtab ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_TO_DIF - Fatal error!'
    write ( *, '(a)' ) '  Two entries of XTAB are equal!'
    stop
  end if
!
!  Copy the data values into DIFTAB.
!
  diftab(1:ntab) = ytab(1:ntab)
!
!  Compute the divided differences.
!
  do i = 2, ntab
    do j = ntab, i, -1

      diftab(j) = ( diftab(j) - diftab(j-1) ) / ( xtab(j) - xtab(j+1-i) )

    end do
  end do

  return
end
subroutine data_to_dif_display ( ntab, xtab, ytab, diftab )

!*****************************************************************************80
!
!! DATA_TO_DIF_DISPLAY computes a divided difference table and shows how.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of pairs of points
!    (XTAB(I),YTAB(I)) which are to be used as data.  The
!    number of entries to be used in DIFTAB, XTAB and YTAB.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values at which data was taken.
!
!    Input, real ( kind = 8 ) YTAB(NTAB), the corresponding Y values.
!
!    Output, real ( kind = 8 ) DIFTAB(NTAB), the divided difference
!    coefficients corresponding to the input (XTAB,YTAB).
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) diftab(ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r8vec_distinct
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) ytab(ntab)

  if ( .not. r8vec_distinct ( ntab, xtab ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_TO_DIF_DISPLAY - Fatal error!'
    write ( *, '(a)' ) '  Two entries of XTAB are equal!'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Divided difference table:'
  write ( *, '(a)' ) ' '
  write ( *, '(6x,5g14.6)' ) xtab(1:ntab)
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i3,1x,5g14.6)' ) 0, ytab(1:ntab)
!
!  Copy the data values into DIFTAB.
!
  diftab(1:ntab) = ytab(1:ntab)
!
!  Compute the divided differences.
!
  do i = 2, ntab
    do j = ntab, i, -1

      diftab(j) = ( diftab(j) - diftab(j-1) ) / ( xtab(j) - xtab(j+1-i) )

    end do

    write ( *, '(2x,i3,1x,5g14.6)' ) i-1, diftab(i:ntab)

  end do

  return
end
subroutine data_to_r8poly ( ntab, xtab, ytab, c )

!*****************************************************************************80
!
!! DATA_TO_R8POLY computes the coefficients of a polynomial interpolating data.
!
!  Discussion:
!
!    Space can be saved by using a single array for both the C and
!    YTAB parameters.  In that case, the coefficients will
!    overwrite the Y data without interfering with the computation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of data points.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), YTAB(NTAB), the data values.
!
!    Output, real ( kind = 8 ) C(NTAB), the coefficients of the 
!    polynomial that passes through the data (XTAB,YTAB).  C(1) is the 
!    constant term.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) c(ntab)
  real ( kind = 8 ) diftab(ntab)
  logical r8vec_distinct
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) ytab(ntab)

  if ( .not. r8vec_distinct ( ntab, xtab ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_TO_POLY - Fatal error!'
    write ( *, '(a)' ) '  Two entries of XTAB are equal!'
    stop
  end if

  call data_to_dif ( ntab, xtab, ytab, diftab )

  call dif_to_r8poly ( ntab, xtab, diftab, c )

  return
end
subroutine dif_antideriv ( nd, xd, yd, na, xa, ya )

!*****************************************************************************80
!
!! DIF_ANTIDERIV computes the antiderivative of a divided difference polynomial.
!
!  Discussion:
!
!    The routine uses the divided difference representation 
!    of a polynomial to compute the divided difference representation
!    of the antiderivative of the polynomial.
!
!    The antiderivative of a polynomial P(X) is any polynomial Q(X)
!    with the property that d/dX Q(X) = P(X).
!
!    This routine chooses the antiderivative whose constant term is zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the size of the difference table.
!
!    Input, real ( kind = 8 ) XD(ND), the abscissas of the 
!    difference table.
!
!    Input, real ( kind = 8 ) YD(ND), the difference table.
!
!    Input, integer ( kind = 4 ) NA, the size of the difference table for the
!    antiderivative, which will be ND+1.
!
!    Output, real ( kind = 8 ) XA(NA), the abscissas of the 
!    difference table for the antiderivative.
!
!    Output, real ( kind = 8 ) YA(NA), the difference table
!    for the antiderivative.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) na
  real ( kind = 8 ) x0(nd)
  real ( kind = 8 ) xa(nd+1)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) y0(nd)
  real ( kind = 8 ) ya(nd+1)
  real ( kind = 8 ) yd(nd)
!
!  Copy the input data.
!
  x0(1:nd) = xd(1:nd)
  y0(1:nd) = yd(1:nd)
!
!  Compute an equivalent difference table with all abscissas 0.
!
  call dif_shift_zero ( nd, x0, y0 )
!
!  There is one more abscissas for the antiderivative polynomial.
!
  na = nd + 1
  xa(1:na) = 0.0D+00
!
!  Get the antiderivative of the standard form polynomial.
!
  call r8poly_ant_cof ( nd, y0, ya )

  return
end
subroutine dif_append ( ntab, xtab, diftab, xval, yval, ntab2, xtab2, diftab2 )

!*****************************************************************************80
!
!! DIF_APPEND adds a pair of data values to a divided difference table.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the size of the difference table.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X data values.
!
!    Input, real ( kind = 8 ) DIFTAB(NTAB), the difference table.
!
!    Input, real ( kind = 8 ) XVAL, the X data value to be inserted as XTAB(1).
!
!    Input, real ( kind = 8 ) YVAL, the Y data value to be inserted as YTAB(1).
!
!    Output, integer ( kind = 4 ) NTAB2, the updated size of the 
!    difference table.
!
!    Output, real ( kind = 8 ) XTAB2(NTAB2), the updated abscissas.
!
!    Output, real ( kind = 8 ) DIFTAB2(NTAB2), the updated difference table.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) diftab(ntab)
  real ( kind = 8 ) diftab2(ntab+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ntab2
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) xtab2(ntab+1)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval
!
!  Move the original data up one index.
!
  ntab2 = ntab + 1

  do i = ntab2, 2, -1
    xtab2(i) = xtab(i-1)
  end do

  do i = ntab2, 2, -1
    diftab2(i) = diftab(i-1)
  end do
!
!  Insert the new data.
!
  xtab2(1) = xval
  diftab2(1) = yval
!
!  Recompute the difference table.
!
  do i = 2, ntab2
    diftab2(i) = ( diftab2(i) - diftab2(i-1) ) / ( xtab2(i) - xtab2(1) )
  end do

  return
end
subroutine dif_basis ( ntab, xtab, diftab )

!*****************************************************************************80
!
!! DIF_BASIS computes all Lagrange basis polynomials in divided difference form.
!
!  Discussion:
!
!    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
!    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
!    XTAB(J) for J not equal to I, and 1 when J is equal to I.
!
!    The Lagrange basis polynomials have the property that the interpolating
!    polynomial through a set of NTAB data points (XTAB,YTAB) may be
!    represented as
!
!      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
!
!    Higher order interpolation at selected points may be accomplished
!    using repeated X values, and scaled derivative values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of X data points XTAB, 
!    and the number of basis polynomials to compute.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values upon which the 
!    Lagrange basis polynomials are to be based.
!
!    Output, real ( kind = 8 ) DIFTAB(NTAB,NTAB), the set of divided 
!    difference tables.  Column I of DIFTAB contains the table for 
!    the I-th Lagrange basis polynomial.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) diftab(ntab,ntab)
  integer ( kind = 4 ) i
  real ( kind = 8 ) xtab(ntab)
!
!  Initialize DIFTAB to the identity matrix.
!
  diftab(1:ntab,1:ntab) = 0.0D+00
  do i = 1, ntab
    diftab(i,i) = 1.0D+00
  end do
!
!  Compute each Lagrange basis polynomial.
!
  do i = 1, ntab
    call data_to_dif ( ntab, xtab, diftab(1,i), diftab(1,i) )
  end do

  return
end
subroutine dif_basis_i ( ival, ntab, xtab, diftab )

!*****************************************************************************80
!
!! DIF_BASIS_I: I-th Lagrange basis polynomial in divided difference form.
!
!  Discussion:
!
!    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
!    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
!    XTAB(J) for J not equal to I, and 1 when J is equal to I.
!
!    The Lagrange basis polynomials have the property that the interpolating
!    polynomial through a set of NTAB data points (XTAB,YTAB) may be
!    represented as
!
!      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
!
!    Higher order interpolation at selected points may be accomplished
!    using repeated X values, and scaled derivative values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired Lagrange 
!    basis polynomial.  IVAL should be between 1 and NTAB.
!
!    Input, integer ( kind = 4 ) NTAB, the number of data points XTAB.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values upon which the 
!    Lagrange basis polynomial is to be based.
!
!    Output, real ( kind = 8 ) DIFTAB(NTAB), the divided difference table 
!    for the IVAL-th Lagrange basis polynomial.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) diftab(ntab)
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xtab(ntab)
!
!  Check IVAL.
!
  if ( ival < 1 .or. ntab < ival ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIF_BASIS_I - Fatal error!'
    write ( *, '(a,i6)' ) '  IVAL must be between 1 and ', ntab
    write ( *, '(a,i6)' ) '  but your value is ', ival
    stop
  end if
!
!  Initialize DIFTAB to Delta(I,J).
!
  diftab(1:ntab) = 0.0D+00
  diftab(ival) = 1.0D+00
!
!  Compute the IVAL-th Lagrange basis polynomial.
!
  call data_to_dif ( ntab, xtab, diftab, diftab )

  return
end
subroutine dif_deriv_table ( nd, xd, yd, ndp, xdp, ydp )

!*****************************************************************************80
!
!! DIF_DERIV_TABLE computes the divided difference table for a derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the size of the input table.
!
!    Input, real ( kind = 8 ) XD(ND), the abscissas for the divided
!    difference table.
!
!    Input, real ( kind = 8 ) YD(ND), the divided difference table.
!
!    Output, integer ( kind = 4 ) NDP, the size of the output table, 
!    which is ND-1.
!
!    Input, real ( kind = 8 ) XDP(NDP), the abscissas for the divided
!    difference table for the derivative.
!
!    Output, real ( kind = 8 ) YDP(NDP), the divided difference 
!    table for the derivative.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ndp
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xd_temp(nd)
  real ( kind = 8 ) xdp(nd-1)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yd_temp(nd)
  real ( kind = 8 ) ydp(nd-1)
!
!  Using a temporary copy of the difference table, shift the abscissas to zero.
!
  xd_temp(1:nd) = xd(1:nd)
  yd_temp(1:nd) = yd(1:nd)

  call dif_shift_zero ( nd, xd_temp, yd_temp )
!
!  Now construct the derivative.
!
  ndp = nd - 1

  xdp(1:ndp) = 0.0D+00

  do i = 1, ndp
    ydp(i) = real ( i, kind = 8 ) * yd_temp(i+1)
  end do

  return
end
subroutine dif_print ( nd, xd, yd, title )

!*****************************************************************************80
!
!! DIF_PRINT prints the polynomial represented by a divided difference table.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the order of the difference table.
!
!    Input, real ( kind = 8 ) XD(ND), the X values for the polynomial.
!
!    Input, real ( kind = 8 ) YD(ND), the divided difference table
!    for the polynomial.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) i
  character ( len = * ) title
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '( ''  p(x) =                           '', g14.6 )' ) yd(1)

  do i = 2, nd
    write ( *, '( ''       + ( x - '', g14.6, '') * ( '', g14.6 )' ) &
      xd(i-1), yd(i)
  end do

  write ( *, '(80a1)' ) ( '       )', i = 1, nd - 1 )

  return
end
subroutine dif_root ( abserr, fxname, iprint, maxstp, maxtab, &
  relerr, xroot, xtry1, xtry2 )

!*****************************************************************************80
!
!! DIF_ROOT seeks a zero of F(X) using divided difference techniques.
!
!  Discussion:
!
!    The method uses the idea of considering the related function
!
!      H(X) = 1 / F(X)
!
!    The iteration begins with two estimates for the root supplied by
!    the user.
!
!    From the most recent approximation to the root, X(K), the next
!    approximation X(K+1) is determined by:
!
!      X(K+1) = X(K) + H(X(K-R),...,X(K-1)) / H(X(K-R),...,X(K-1),X(K))
!
!    where K-R = 1 until the maximal order NTAB is reached.
!
!    Generally, the next iterate X(K+1) is the zero of a rational function
!    which passes through the previous data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    FM Larkin,
!    Root Finding by Divided Differences,
!    Numerische Mathematik,
!    Volume 37, pages 93-104, 1981.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ABSERR, a positive absolute error tolerance.  
!    If an estimate X for the root is found with ABS ( F(X) ) <= ABSERR, 
!    the iteration is stopped.
!
!    Input, external FXNAME, the name of the function routine which evaluates
!    F(X).  The form of FXNAME must be similar to the following function which
!    has F(X) = ( X - 1 ) * ( X + 1 ).
!
!    function parab ( x )
!
!      real ( kind = 8 ) parab
!      real ( kind = 8 ) x
!
!      parab = ( x - 1.0D+00 ) * ( x + 1.0D+00 )
!
!      return
!    end
!
!    Input, integer ( kind = 4 ) IPRINT, a switch controlling printed output:
!    0, only print error messages.
!    nonzero, also print a table of the iterative process.
!
!    Input, integer ( kind = 4 ) MAXSTP, the limit on how many iterations 
!    may be tried.
!
!    Input, integer ( kind = 4 ) MAXTAB, the limit on how high an order can be 
!    used in the divided difference table.  MAXTAB must be at least 2, and 
!    probably should not be too large.  Perhaps a value of 5 or 6 is reasonable,
!    20 is too large.
!
!    Input, real ( kind = 8 ) RELERR, a tolerance on the size of the change 
!    in the root estimates.  If a step is taken, and the change in the root
!    estimate is less than RELERR, the iteration will stop.
!
!    Output, real ( kind = 8 ) XROOT, the point which the program has 
!    produced as an approximate root.
!    Either ABS ( F(XROOT) ) <= ABSERR, or the maximum number of steps was
!    reached, or the current estimate of the root could not be significantly
!    improved.
!
!    Input, real ( kind = 8 ) XTRY1, XTRY2, two initial approximations to 
!    the root, supplied by the user, which must be distinct.
!
  implicit none

  integer ( kind = 4 ) maxtab

  real ( kind = 8 ) abserr
  real ( kind = 8 ) diftab(maxtab)
  real ( kind = 8 ) froot
  real ( kind = 8 ) ftemp1
  real ( kind = 8 ) ftemp2
  real ( kind = 8 ), external :: fxname
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) maxstp
  integer ( kind = 4 ) ntab
  real ( kind = 8 ) relerr
  real ( kind = 8 ) xdelt
  real ( kind = 8 ) xold
  real ( kind = 8 ) xroot
  real ( kind = 8 ) xtab(maxtab)
  real ( kind = 8 ) xtry1
  real ( kind = 8 ) xtry2
  real ( kind = 8 ) yval
!
!  Make sure XTRY1 and XTRY2 are not equal.
!
  if ( xtry1 == xtry2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIF_ROOT - Fatal error!'
    write ( *, '(a)' ) '  XTRY1 = XTRY2 on input.'
    stop
  end if
!
!  Make sure MAXTAB is at least 2.
!
  if ( maxtab < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIF_ROOT - Fatal error!'
    write ( *, '(a)' ) '  MAXTAB < 2 on input!'
    stop
  end if

  xtab(1) = xtry1
  xtab(2) = xtry2
  ftemp1 = fxname ( xtry1 )
  ftemp2 = fxname ( xtry2 )

  if ( abs ( ftemp2 ) < abs ( ftemp1 ) ) then
    xtab(1) = xtry2
    xtab(2) = xtry1
    call r8_swap ( ftemp1, ftemp2 )
  end if
!
!  Initialize the number of steps.
!
  istep = 0
!
!  Initialize the number of data points.
!
  ntab = 2

  if ( 0 < iprint ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   Step  NTAB    XROOT        F(XROOT)      XDELT'
    write ( *, '(a)' ) ' '
  end if
!
!  Initialize the divided difference table data.
!
  diftab(1) = 1.0D+00 / ftemp1
  diftab(2) = 1.0D+00 / ftemp2

  call data_to_dif ( ntab, xtab, diftab, diftab )
!
!  Initialize values used in the iteration.
!
  xroot = xtry1
  froot = ftemp1
  xdelt = xtry1 - xtry2
!
!  Does the starting data already satisfy the function norm
!  error tolerance ABSERR, or the interval norm error tolerance
!  RELERR?
!
  do

    if ( 0 < iprint ) then
      write ( *, '(3x,i4,4x,i2, 3g14.6)' ) istep, ntab, xroot, froot, xdelt
    end if

    if ( abs ( froot ) <= abserr ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIF_ROOT - Absolute convergence,'
      write ( *, '(a)' ) '  The function value meets the error tolerance.'
      exit
    end if

    if ( abs ( xdelt ) <= relerr ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIF_ROOT - Relative convergence.'
      write ( *, '(a)' ) '  The stepsize meets the error tolerance.'
      exit
    end if
!
!  Check the number of steps taken.
!
    if ( maxstp <= istep ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIF_ROOT - Nonconvergence!'
      write ( *, '(a)' ) '  The maximum number of steps was taken.'
      exit
    end if
!
!  Generate the next point, XVAL.
!
    xold = xroot
    istep = istep + 1
!
!  Algorithm breakdown: The divisor DIFTAB(NTAB) is zero.
!
    if ( diftab(ntab) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIF_ROOT - Fatal error!'
      write ( *, '(a,i6)' ) '  Algorithm using differences of order ', ntab
      write ( *, '(a)' ) '  A zero-divisor was computed.'
      write ( *, '(a)' ) '  The algorithm has broken down.'
      write ( *, '(a)' ) '  Examine the results.  They may be useful.'
      write ( *, '(a)' ) '  Perhaps a lower value of MAXTAB would help.'
      stop
    end if

    xroot = xtab(ntab) + ( diftab(ntab-1) / diftab(ntab) )
    xdelt = xroot - xold
    froot = fxname ( xroot )

    if ( abs ( froot ) <= abserr ) then
      cycle
    end if

    yval = 1.0D+00 / froot
!
!  If we are now using MAXTAB points, we have to remove an old
!  one before adding the new one.
!
    if ( maxtab <= ntab ) then
      ntab = ntab - 1
    end if

    call dif_append ( ntab, xtab, diftab, xroot, yval, ntab, xtab, diftab )

  end do

  return
end
subroutine dif_shift_x ( nd, xd, yd, xv )

!*****************************************************************************80
!
!! DIF_SHIFT_X replaces one abscissa of a divided difference table.
!
!  Discussion:
!
!    The routine shifts the representation of a divided difference polynomial by
!    dropping the last X value in XD, and adding a new value XV to the
!    beginning of the XD array, suitably modifying the coefficients stored
!    in YD.
!
!    The representation of the polynomial is changed, but the polynomial itself
!    should be identical.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of divided difference 
!    coefficients, and the number of entries in XD.
!
!    Input/output, real ( kind = 8 ) XD(ND), the X values used in 
!    the representation of the divided difference polynomial.
!    After a call to this routine, the last entry of XD has been dropped,
!    the other entries have shifted up one index, and XV has been inserted
!    at the beginning of the array.
!
!    Input/output, real ( kind = 8 ) YD(ND), the divided difference
!    coefficients corresponding to the XD array.  On output, this 
!    array has been adjusted.
!
!    Input, real ( kind = 8 ) XV, a new X value which is to be used in 
!    the representation of the polynomial.  On output, XD(1) equals 
!    XV and the representation of the polynomial has been suitably changed.
!    Note that XV does not have to be distinct from any of the original XD
!    values.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) i
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xv
  real ( kind = 8 ) yd(nd)
!
!  Recompute the divided difference coefficients.
!
  do i = nd - 1, 1, -1
    yd(i) = yd(i) + ( xv - xd(i) ) * yd(i+1)
  end do
!
!  Shift the XD values up one position, and insert XV at the beginning.
!
  xd(2:nd) = xd(1:nd-1)

  xd(1) = xv

  return
end
subroutine dif_shift_zero ( nd, xd, yd )

!*****************************************************************************80
!
!! DIF_SHIFT_ZERO shifts a divided difference table so all abscissas are zero.
!
!  Discussion:
!
!    When the abscissas are changed, the coefficients naturally
!    must also be changed.
!
!    The resulting pair (XD, YD) still represents the
!    same polynomial, but the entries in YD are now the
!    standard polynomial coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the length of the table.
!
!    Input/output, real ( kind = 8 ) XD(ND), the abscissas for the divided
!    difference table.  On output, all entries have been reset to 0, but
!    (XD,YD) can still be regarded as a divided difference table for the input
!    polynomial.
!
!    Input/output, real ( kind = 8 ) YD(ND).  On input, the divided difference 
!    table for the polynomial.  On output, the divided difference table for
!    the polynomial, which has been rebased at 0.  Hence, YD is also simply
!    the coefficient array for the standard representation of the polynomial.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)

  do j = 1, nd

    do i = nd - 1, 1, -1
      yd(i) = yd(i) - xd(i) * yd(i+1)
    end do
!
!  Shift the XD values up one position.
!
    xd(2:nd) = xd(1:nd-1)

    xd(1) = 0.0D+00

  end do

  return
end
subroutine dif_to_r8poly ( n, xd, yd, c )

!*****************************************************************************80
!
!! DIF_TO_R8POLY converts a divided difference polynomial to standard form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of coefficients, and abscissas.
!
!    Input, real ( kind = 8 ) XD(N), the X values used in the divided
!    difference representation of the polynomial.
!
!    Input, real ( kind = 8 ) YD(N) the divided difference table.
!
!    Output, real ( kind = 8 ) C(N), the polyomial coefficients.
!    C(1) is the constant term, and C(N) is the coefficient
!    of X^(N-1).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) xd(n)
  real ( kind = 8 ) yd(n)

  c(1:n) = yd(1:n)
!
!  Recompute the divided difference coefficients.
!
  do j = 1, n - 1
    do i = 1, n - j
      c(n-i) = c(n-i) - xd(n-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
subroutine dif_val ( ntab, xtab, diftab, xv, yv )

!*****************************************************************************80
!
!! DIF_VAL evaluates a divided difference polynomial at a point.
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
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of divided difference
!    coefficients in DIFTAB, and the number of points XTAB.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values upon which the
!    divided difference polynomial is based.
!
!    Input, real ( kind = 8 ) DIFTAB(NTAB), the divided difference 
!    polynomial coefficients.
!
!    Input, real ( kind = 8 ) XV, the value where the polynomial 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) YV, the value of the polynomial at XV.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) diftab(ntab)
  integer ( kind = 4 ) i
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) xv
  real ( kind = 8 ) yv

  yv = diftab(ntab)
  do i = ntab - 1, 1, -1
    yv = diftab(i) + ( xv - xtab(i) ) * yv
  end do

  return
end
subroutine dif_vals ( nd, xd, yd, nv, xv, yv )

!*****************************************************************************80
!
!! DIF_VALS evaluates a divided difference polynomial at a set of points.
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
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the order of the difference table.
!
!    Input, real ( kind = 8 ) XD(ND), the X values of the difference table.
!
!    Input, real ( kind = 8 ) YD(ND), the divided differences.
!
!    Input, integer ( kind = 4 ) NV, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XV(NV), the evaluation points.
!
!    Output, real ( kind = 8 ) YV(NV), the value of the divided difference
!    polynomial at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nv

  integer ( kind = 4 ) i
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xv(nv)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yv(nv)

  yv(1:nv) = yd(nd)
  do i = nd - 1, 1, -1
    yv(1:nv) = yd(i) + ( xv(1:nv) - xd(i) ) * yv(1:nv)
  end do

  return
end
subroutine lagrange_rule ( n, x, w )

!*****************************************************************************80
!
!! LAGRANGE_RULE computes the weights of a Lagrange interpolation rule.
!
!  Discussion:
!
!    Given N abscissas X, an arbitrary function F(X) can be
!    interpolated by a polynomial P(X) of order N (and degree N-1)
!    using weights that depend only on X.
!
!    Standard Lagrange interpolation can be rewritten into this form,
!    which is more economical than evaluating the individual Lagrange
!    basis polynomials.
!
!    If we define
!
!      W(I) = 1 / product ( 1 <= J <= N, J /= I ) ( X(J) - X(I) )
!
!    then
!
!      P(XV) = sum ( 1 <= I <= N ) W(I) * F( X(I) ) / ( XV - X(I) )
!            / sum ( 1 <= I <= N ) W(I)             / ( XV - X(I) )
!
!    except when XV = X(J), for some J, when we set:
!
!      P(X(J)) = F(X(J))
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Paul Berrut, Lloyd Trefethen,
!    Barycentric Lagrange Interpolation,
!    SIAM Review,
!    Volume 46, Number 3, September 2004, pages 501-517.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Input, real ( kind = 8 ) X(N), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) W(N), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 1.0D+00

  do i = 1, n
    do j = 1, i - 1
      w(j) = ( x(i) - x(j) ) * w(j)
    end do
    w(i) = product ( ( x(1:i-1) - x(i) ) )
  end do

  w(1:n) = 1.0D+00 / w(1:n)

  return
end
subroutine lagrange_sum ( n, x, w, y, xv, yv )

!*****************************************************************************80
!
!! LAGRANGE_SUM carries out a Lagrange interpolation rule.
!
!  Discussion:
!
!    It is assumed that LAGRANGE_RULE has already been called to compute
!    the appropriate weights for the given set of abscissas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Paul Berrut, Lloyd Trefethen,
!    Barycentric Lagrange Interpolation,
!    SIAM Review,
!    Volume 46, Number 3, September 2004, pages 501-517.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Input, real ( kind = 8 ) X(N), the abscissas of the rule.
!
!    Input, real ( kind = 8 ) W(N), the weights of the rule.
!
!    Input, real ( kind = 8 ) Y(N), the function values at the abscissas.
!
!    Input, real ( kind = 8 ) XV, a point where an interpolated value is
!    needed.
!
!    Output, real ( kind = 8 ) YV, the interpolated function value.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bot
  integer ( kind = 4 ) i
  real ( kind = 8 ) top
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xv
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yv

  do i = 1, n

    if ( xv == x(i) ) then
      yv = y(i)
      return
    end if

  end do

  top = 0.0D+00
  bot = 0.0D+00
  
  do i = 1, n
    top = top + w(i) * y(i) / ( xv - x(i) )
    bot = bot + w(i)        / ( xv - x(i) )
  end do

  yv = top / bot

  return
end
subroutine lagrange_val ( n, x, y, xv, yv )

!*****************************************************************************80
!
!! LAGRANGE_VAL applies a naive form of Lagrange interpolation.
!
!  Discussion:
!
!    Given N abscissas X, an arbitrary function Y(X) can be
!    interpolated by a polynomial P(X) of order N (and degree N-1)
!    using Lagrange basis polynomials of degree N-1.
!
!    Standard Lagrange interpolation can be rewritten into this form,
!    which is more economical than evaluating the individual Lagrange
!    basis polynomials.
!
!    If we define
!
!      L(I)(XV) = product ( 1 <= J <= N, J /= I )
!        ( XV - X(J) ) / ( X(I) - X(J) )
!
!    then
!
!      P(XV) = sum ( 1 <= I <= N ) Y( X(I) ) * L(I)(XV)
!
!    Applying this form of the interpolation rule directly involves 
!    about N^2 work.  There are more efficient forms of the rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Input, real ( kind = 8 ) Y(N), the function values at the abscissas.
!
!    Input, real ( kind = 8 ) XV, a point where an interpolated value is
!    needed.
!
!    Output, real ( kind = 8 ) YV, the interpolated function value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) poly
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xv
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yv

  yv = 0.0D+00

  do i = 1, n
    poly = 1.0D+00
    do j = 1, n
      if ( j /= i ) then
        poly = poly * ( xv - x(j) ) / ( x(i) - x(j) )
      end if
    end do
    yv = yv + y(i) * poly
  end do

  return
end
subroutine nc_rule ( n, a, b, xtab, weight )

!*****************************************************************************80
!
!! NC_RULE computes the weights of a Newton-Cotes quadrature rule.
!
!  Discussion:
!
!    For the interval [A,B], the Newton-Cotes quadrature rule estimates
!
!      Integral ( A <= X <= B ) F(X) dX
!
!    using N equally spaced abscissas XTAB(I) and a weight vector
!    WEIGHT(I):
!
!      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the CLOSED rule, the abscissas include the points A and B.
!    For the OPEN rule, the abscissas do not include A and B.
!
!    For the common closed and open rules, the abscissas are equally spaced.
!
!    This routine allows the user to specify any set of abscissas;
!    hence, it can compute the standard open and closed rules, and
!    other variations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints of the interval
!    over which the quadrature rule is to be applied.
!
!    Input, real ( kind = 8 ) XTAB(N), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(N), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) poly_cof(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) weight(n)
  real ( kind = 8 ) xtab(n)
  real ( kind = 8 ) yvala
  real ( kind = 8 ) yvalb

  do i = 1, n
!
!  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
!  and zero at the other nodes.
!
    call r8poly_basis_1 ( i, n, xtab, poly_cof )
!
!  Evaluate the antiderivative of the polynomial at the left and
!  right endpoints.
!
    call r8poly_ant_val ( n, poly_cof, a, yvala )

    call r8poly_ant_val ( n, poly_cof, b, yvalb )

    weight(i) = yvalb - yvala

  end do

  return
end
subroutine ncc_rule ( norder, xtab, weight )

!*****************************************************************************80
!
!! NCC_RULE computes the coefficients of a Newton-Cotes closed quadrature rule.
!
!  Discussion:
!
!    For the interval [-1,1], the Newton-Cotes quadrature rule estimates
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    using NORDER equally spaced abscissas XTAB(I) and a weight vector
!    WEIGHT(I):
!
!      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the CLOSED rule, the abscissas include A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) norder

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)
!
!  Compute a closed quadrature rule.
!
  a = -1.0D+00
  b = 1.0D+00

  do i = 1, norder
    xtab(i) = ( real ( norder - i,     kind = 8 ) * a   &
              + real (          i - 1, kind = 8 ) * b ) &
              / real ( norder - 1,     kind = 8 )
  end do

  call nc_rule ( norder, a, b, xtab, weight )

  return
end
subroutine nco_rule ( norder, xtab, weight )

!*****************************************************************************80
!
!! NCO_RULE computes the coefficients of a Newton-Cotes open quadrature rule.
!
!  Discussion:
!
!    For the interval [-1,1], the Newton-Cotes quadrature rule estimates
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    using NORDER equally spaced abscissas XTAB(I) and a weight vector
!    WEIGHT(I):
!
!      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the OPEN rule, the abscissas do not include A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the  rule.
!
  implicit none

  integer ( kind = 4 ) norder

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) weight(norder)
  real ( kind = 8 ) xtab(norder)

  a = -1.0D+00
  b = 1.0D+00

  do i = 1, norder
    xtab(i) = ( real ( norder + 1 - i, kind = 8 ) * a   &
              + real (              i, kind = 8 ) * b ) &
              / real ( norder + 1,     kind = 8 )
  end do

  call nc_rule ( norder, a, b, xtab, weight )

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8poly_ant_cof ( n, poly_cof, poly_cof2 )

!*****************************************************************************80
!
!! R8POLY_ANT_COF integrates a polynomial in standard form.
!
!  Discussion:
!
!    The antiderivative of a polynomial P(X) is any polynomial Q(X)
!    with the property that d/dX Q(X) = P(X).
!
!    This routine chooses the antiderivative whose constant term is zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) POLY_COF(N), the polynomial coefficients.
!    POLY_COF(1) is the constant term, and POLY_COF(N) is the
!    coefficient of X^(N-1).
!
!    Output, real ( kind = 8 ) POLY_COF2(N+1), the coefficients of 
!    the antiderivative polynomial, in standard form.  The constant 
!    term is set to zero.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) poly_cof(n)
  real ( kind = 8 ) poly_cof2(n+1)
  integer ( kind = 4 ) i
!
!  Set the constant term.
!
  poly_cof2(1) = 0.0D+00
!
!  Integrate the polynomial.
!
  do i = 2, n + 1
    poly_cof2(i) = poly_cof(i-1) / real ( i - 1, kind = 8 )
  end do

  return
end
subroutine r8poly_ant_val ( n, c, xv, yv )

!*****************************************************************************80
!
!! R8POLY_ANT_VAL evaluates the antiderivative of a polynomial in standard form.
!
!  Discussion:
!
!    The constant term of the antiderivative is taken to be zero.
!
!    Thus, if 
!      p(x) = c(1)     + c(2) * x   + c(3) * x^2   + ... + c(n) * x^(n-1)
!    then q(x), the antiderivative, is taken to be:
!      q(x) = c(1) * x + c(2) * x/2 + c(3) * x^3/3 + ... + c(n) * x^n/n
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) C(N), the polynomial coefficients.
!    C(1) is the constant term, and C(N) is the coefficient of X^(N-1).
!
!    Input, real ( kind = 8 ) XV, the evaluation point.
!
!    Output, real ( kind = 8 ) YV, the value of the antiderivative of 
!    the polynomial at XV.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) xv
  real ( kind = 8 ) yv

  yv = 0.0D+00

  do i = n, 1, -1
    yv = ( yv + c(i) / real ( i, kind = 8 ) ) * xv
  end do

  return
end
subroutine r8poly_basis ( ntab, xtab, poly_cof )

!*****************************************************************************80
!
!! R8POLY_BASIS computes all Lagrange basis polynomials in standard form.
!
!  Discussion:
!
!    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
!    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
!    XTAB(J) for J not equal to I, and 1 when J is equal to I.
!
!    The Lagrange basis polynomials have the property that the interpolating
!    polynomial through a set of NTAB data points (XTAB,YTAB) may be
!    represented as
!
!      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
!
!    Higher order interpolation at selected points may be accomplished
!    using repeated X values, and scaled derivative values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of data points XTAB.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values upon which the 
!    Lagrange basis polynomial is to be based.
!
!    Output, real ( kind = 8 ) POLY_COF(NTAB,NTAB), the polynomial 
!    coefficients for the I-th Lagrange basis polynomial are stored 
!    in column I.  POLY_COF(1,I) is the constant term, and POLY_COF(1,NTAB) 
!    is the coefficient of X^(NTAB-1).
!
  implicit none

  integer ( kind = 4 ) ntab

  integer ( kind = 4 ) i
  real ( kind = 8 ) poly_cof(ntab,ntab)
  real ( kind = 8 ) xtab(ntab)
!
!  Initialize POLY_COF to the identity matrix.
!
  poly_cof(1:ntab,1:ntab) = 0.0D+00
  do i = 1, ntab
    poly_cof(i,i) = 1.0D+00
  end do
!
!  Compute the divided difference table for the IVAL-th Lagrange basis
!  polynomial.
!
  do i = 1, ntab
    call data_to_dif ( ntab, xtab, poly_cof(1,i), poly_cof(1,i) )
  end do
!
!  Convert the divided difference table coefficients to standard polynomial
!  coefficients.
!
  do i = 1, ntab
    call dif_to_r8poly ( ntab, xtab, poly_cof(1,i), poly_cof(1,i) )
  end do

  return
end
subroutine r8poly_basis_1 ( ival, ntab, xtab, poly_cof )

!*****************************************************************************80
!
!! R8POLY_BASIS_1 computes the I-th Lagrange basis polynomial in standard form.
!
!  Discussion:
!
!    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
!    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
!    XTAB(J) for J not equal to I, and 1 when J is equal to I.
!
!    The Lagrange basis polynomials have the property that the interpolating
!    polynomial through a set of NTAB data points (XTAB,YTAB) may be
!    represented as
!
!      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
!
!    Higher order interpolation at selected points may be accomplished
!    using repeated X values, and scaled derivative values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired Lagrange 
!    basis polynomial.  IVAL should be between 1 and NTAB.
!
!    Input, integer ( kind = 4 ) NTAB, the number of data points XTAB.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values upon which the 
!    Lagrange basis polynomial is to be based.
!
!    Output, real ( kind = 8 ) POLY_COF(NTAB), the polynomial 
!    coefficients for the IVAL-th Lagrange basis polynomial.
!
  implicit none

  integer ( kind = 4 ) ntab

  integer ( kind = 4 ) ival
  real ( kind = 8 ) poly_cof(ntab)
  real ( kind = 8 ) xtab(ntab)
!
!  Check IVAL.
!
  if ( ival < 1 .or. ntab < ival ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_BASE_1 - Fatal error!'
    write ( *, '(a,i6)' ) '  IVAL must be between 1 and ', ntab
    write ( *, '(a,i6)' ) '  but your value is ', ival
    stop
  end if
!
!  Initialize POLY_COF to the IVAL-th column of the identity matrix.
!
  poly_cof(1:ntab) = 0.0D+00
  poly_cof(ival) = 1.0D+00
!
!  Compute the divided difference table for the IVAL-th Lagrange basis
!  polynomial.
!
  call data_to_dif ( ntab, xtab, poly_cof, poly_cof )
!
!  Convert the divided difference table coefficients to standard polynomial
!  coefficients.
!
  call dif_to_r8poly ( ntab, xtab, poly_cof, poly_cof )

  return
end
subroutine r8poly_der_cof ( n, poly_cof, poly_cof2 )

!*****************************************************************************80
!
!! R8POLY_DER_COF computes the coefficients of the derivative of a polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) POLY_COF(N), the coefficients of the 
!    polynomial to be differentiated.  POLY_COF(1) is the constant term, and
!    POLY_COF(N) is the coefficient of X^(N-1).
!
!    Output, real ( kind = 8 ) POLY_COF2(N-1), the coefficients of the
!    derivative of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) poly_cof(n)
  real ( kind = 8 ) poly_cof2(n-1)
  integer ( kind = 4 ) i

  do i = 1, n-1
    poly_cof2(i) = real ( i, kind = 8 ) * poly_cof(i+1)
  end do

  return
end
subroutine r8poly_der_val ( n, poly_cof, xval, yval )

!*****************************************************************************80
!
!! R8POLY_DER_VAL evaluates the derivative of a polynomial in standard form.
!
!  Discussion:
!
!    A polynomial in standard form, with coefficients POLY_COF(*),
!    may be written:
!
!      P(X) = POLY_COF(1)
!           + POLY_COF(2) * X
!           ...
!           + POLY_COF(N) * X^(N-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) POLY_COF(N), the polynomial coefficients.
!    POLY_COF(1) is the constant term, and POLY_COF(N) is the coefficient of
!    X^(N-1).
!
!    Input, real ( kind = 8 ) XVAL, a value where the derivative of the
!    polynomial is to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the derivative of the
!    polynomial at XVAL.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) poly_cof(n)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval

  yval = real ( n - 1, kind = 8 ) * poly_cof(n)
  do i = n - 1, 2, -1
    yval = yval * xval + real ( i - 1, kind = 8 ) * poly_cof(i)
  end do

  return
end
subroutine r8poly_order ( na, a, order )

!*****************************************************************************80
!
!! R8POLY_ORDER returns the order of a polynomial.
!
!  Discussion:
!
!    The order of a polynomial is the degree plus 1.
!
!    The order of a constant polynomial is 1.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) ORDER, the order of A.
!
  implicit none

  integer ( kind = 4 ) na

  real ( kind = 8 ) a(na)
  integer ( kind = 4 ) order

  order = na

  do while ( 1 < order )

    if ( a(order) /= 0.0D+00 ) then
      return
    end if

    order = order - 1

  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) A(N), the polynomial coefficients.
!    A(1) is the constant term and
!    A(N) is the coefficient of X^(N-1).
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  
  call r8poly_order ( n, a, n2 )

  if ( a(n2) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 3 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2-1
  else if ( n2 == 2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) plus_minus, mag
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2 - 1, 1, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 3 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i-1
      else if ( i == 2 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8poly_shift ( scale, shift, n, poly_cof )

!*****************************************************************************80
!
!! R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
!
!  Discussion:
!
!    Assuming P(X) is a polynomial in the argument X, of the form:
!
!      P(X) = C(1)
!           + C(2) * X
!           + ...
!           + C(N) * X^(N-1)
!
!    and that Z is related to X by the formula:
!
!      Z = SCALE * X + SHIFT
!
!    then this routine computes coefficients C for the polynomial Q(Z):
!
!      Q(Z) = C(1)
!           + C(2) * Z
!           + ...
!           + C(N) * Z^(N-1)
!
!    so that:
!
!      Q(Z(X)) = P(X)
!
!  Example:
!
!    P(X) = 2 * X^2 - X + 6
!
!    Z = 2.0D+00 * X + 3.0D+00
!
!    Q(Z) = 0.5 *         Z^2 -  3.5 * Z + 12
!
!    Q(Z(X)) = 0.5 * ( 4.0D+00 * X^2 + 12.0D+00 * X +  9 )
!            - 3.5 * (                  2.0D+00 * X +  3 )
!                                                   + 12
!
!            = 2.0D+00         * X^2 -  1.0D+00 * X +  6
!
!            = P(X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SHIFT, SCALE, the shift and scale applied to X,
!    so that Z = SCALE * X + SHIFT.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial
!
!    Input/output, real ( kind = 8 ) POLY_COF(N).
!    On input, the coefficient array in terms of the X variable.
!    On output, the coefficient array in terms of the Z variable.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) poly_cof(n)
  real ( kind = 8 ) scale
  real ( kind = 8 ) shift

  do i = 1, n
    poly_cof(i+1:n) = poly_cof(i+1:n) / scale
  end do

  do i = 1, n
    do j = n-1, i, -1
      poly_cof(j) = poly_cof(j) - shift * poly_cof(j+1)
    end do
  end do

  return
end
subroutine r8poly_val_horner ( n, poly_cof, xval, yval )

!*****************************************************************************80
!
!! R8POLY_VAL_HORNER evaluates a polynomial in standard form.
!
!  Discussion:
!
!    A polynomial in standard form, with coefficients POLY_COF(*),
!    may be written:
!
!      P(X) = C(1)
!           + C(2) * X
!           ...
!           + C(N) * X^(N-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) POLY_COF(N), the polynomial coefficients.
!    POLY_COF(1) is the constant term, and POLY_COF(N) is the coefficient of
!    X^(N-1).
!
!    Input, real ( kind = 8 ) XVAL, a value where the polynomial is to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the polynomial at XVAL.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) poly_cof(n)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval

  yval = 0.0D+00

  do i = n, 1, -1
    yval = yval * xval + poly_cof(i)
  end do

  return
end
function r8vec_distinct ( n, x )

!*****************************************************************************80
!
!! R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector to be checked.
!
!    Output, logical R8VEC_DISTINCT is .TRUE. if all N elements of X
!    are distinct.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r8vec_distinct
  real ( kind = 8 ) x(n)

  r8vec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1
      if ( x(i) == x(j) ) then
        return
      end if
    end do
  end do

  r8vec_distinct = .true.

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns N values, evenly spaced between ALO and AHI.
!
!  Discussion:
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
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
subroutine r8vec_even_select ( n, xlo, xhi, ival, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
!
!  Discussion:
!
!    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
!
!    Unless N = 1, X(1) = XLO and X(N) = XHI.
!
!    If N = 1, then X(1) = 0.5*(XLO+XHI).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) XLO, XHI, the low and high values.
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any value.
!
!    Output, real ( kind = 8 ) XVAL, the IVAL-th of N evenly spaced values
!    between XLO and XHI.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ival
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xval

  if ( n == 1 ) then

    xval = 0.5D+00 * ( xlo + xhi )

  else

    xval = ( real ( n - ival,     kind = 8 ) * xlo   &
           + real (     ival - 1, kind = 8 ) * xhi ) &
           / real ( n        - 1, kind = 8 )

  end if

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets an R8VEC to the indicator vector A(I)=I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = *  ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, a(i)
  end do

  return
end
subroutine roots_to_dif ( nroots, roots, ntab, xtab, diftab )

!*****************************************************************************80
!
!! ROOTS_TO_DIF sets a divided difference table for a polynomial from its roots.
!
!  Discussion:
!
!    This turns out to be a simple task, because of two facts:
!
!    * The divided difference polynomial of one smaller degree which
!      passes through the values ( ROOT(I), 0 ) is the zero polynomial,
!      and hence has a zero divided difference table.
!
!    * We want a polynomial of one degree higher, but we don't want it
!      to pass through an addditional point.  Instead, we specify that
!      the polynomial is MONIC.  This means that the divided difference
!      table is almost the same as for the zero polynomial, except that
!      there is one more pair of entries, an arbitrary X value, and
!      a Y value of 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROOTS, is the number of roots.
!
!    Input, real ( kind = 8 ) ROOTS(NROOTS), the roots of
!    the polynomial.
!
!    Output, integer ( kind = 4 ) NTAB, is equal to NROOTS+1.
!
!    Output, real ( kind = 8 ) XTAB(NTAB), the abscissas of the divided
!    difference table.
!
!    Output, real ( kind = 8 ) DIFTAB(NTAB), the divided difference
!    table.
!
  implicit none

  integer ( kind = 4 ) nroots

  real ( kind = 8 ) diftab(nroots+1)
  integer ( kind = 4 ) ntab
  real ( kind = 8 ) roots(nroots)
  real ( kind = 8 ) xtab(nroots+1)

  ntab = nroots + 1
!
!  Build the appropriate difference table for the polynomial
!  through ( ROOTS(I), 0 ) of order NTAB-1.
!
  diftab(1:ntab-1) = 0.0D+00
!
!  Append the extra data to make a monic polynomial of order NTAB
!  which is zero at the NTAB-1 roots.
!
  xtab(1:ntab-1) = roots(1:ntab-1)
  xtab(ntab) = 0.0D+00

  diftab(ntab) = 1.0D+00

  return
end
subroutine roots_to_r8poly ( nroots, roots, nc, c )

!*****************************************************************************80
!
!! ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROOTS, the number of roots specified.
!
!    Input, real ( kind = 8 ) ROOTS(NROOTS), the roots.
!
!    Output, integer ( kind = 4 ) NC, the order of the polynomial, which will
!    be NROOTS + 1.
!
!    Output, real ( kind = 8 ) C(NC), the coefficients of the polynomial.
!
  implicit none

  integer ( kind = 4 ) nroots

  real ( kind = 8 ) c(nroots+1)
  integer ( kind = 4 ) nc
  real ( kind = 8 ) roots(nroots)
  real ( kind = 8 ) xtab(nroots+1)

  nc = nroots + 1
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  xtab(1:nroots) = roots(1:nroots)
  xtab(nc) = 0.0

  c(1:nc-1) = 0.0D+00
  c(nc) = 1.0D+00
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  call dif_shift_zero ( nc, xtab, c )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  character ( len = 8  ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5  ) zone

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
