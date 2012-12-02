subroutine dif_deriv ( nd, xd, yd, ndp, xdp, ydp )

!*****************************************************************************80
!
!! DIF_DERIV computes the derivative of a polynomial in divided difference form.
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
!    Input, integer ( kind = 4 ) ND, the length of the XD and YD arrays.
!
!    Input/output, real ( kind = 8 ) XD(ND), the X values that 
!    correspond to the divided difference table.  On output, XD 
!    contains only zeroes.
!
!    Input/output, real ( kind = 8 ) YD(ND), the divided difference table
!    for the polynomial.  On output, YD is also the coefficient array for 
!    the standard representation of the polynomial.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) i
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xv
  real ( kind = 8 ) yd(nd)

  xv = 0.0D+00

  do i = 1, nd
    call dif_shift_x ( nd, xd, yd, xv )
  end do

  return
end
subroutine dif_to_r8poly ( nd, xd, yd, c )

!*****************************************************************************80
!
!! DIF_TO_R8POLY converts a divided difference table to a standard polynomial.
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
!    Input, integer ( kind = 4 ) ND, the number of coefficients, and abscissas.
!
!    Input, real ( kind = 8 ) XD(ND), the X values used in the divided
!    difference representation of the polynomial.
!
!    Input, real ( kind = 8 ) YD(ND) the divided difference table.
!
!    Output, real ( kind = 8 ) C(ND), the polyomial coefficients.
!    C(1) is the constant term, and C(ND) is the coefficient of X^(ND-1).
!
  implicit none

  integer ( kind = 4 ) nd

  real ( kind = 8 ) c(nd)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)

  c(1:nd) = yd(1:nd)
!
!  Recompute the divided difference coefficients.
!
  do j = 1, nd - 1
    do i = 1, nd - j
      c(nd-i) = c(nd-i) - xd(nd-i-j+1) * c(nd-i+1)
    end do
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
!    23 May 2011
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
  do i = 1, nd - 1
    yv(1:nv) = yd(nd-i) + ( xv(1:nv) - xd(nd-i) ) * yv(1:nv)
  end do

  return
end
subroutine hermite_basis_0 ( n, x, i, xv, value )

!*****************************************************************************80
!
!! HERMITE_BASIS_0 evaluates a zero-order Hermite interpolation basis function.
!
!  Discussion:
!
!    Given ND points XD, with values YD and derivative values YPD, the
!    Hermite interpolant can be written as:
!
!      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
!           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
!
!    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
!    and H1(I;X) is the I-th first order Hermite interpolation basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of abscissas.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Input, integer ( kind = 4 ) I, the index of the first-order basis function.
!
!    Input, real ( kind = 8 ) XV, the evaluation point.
!
!    Output, real ( kind = 8 ) VALUE, the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) factor(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) li
  real ( kind = 8 ) lp
  real ( kind = 8 ) lpp
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xv

  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_BASIS_0 - Fatal error!'
    write ( *, '(a)' ) '  I < 1 or N < I.'
    stop
  end if
!
!  L(X) = product ( X - X(1:N) )
!
!  L'(X(I)).
!
  factor(1:n) = x(i) - x(1:n)
  factor(i) = 1.0D+00

  lp = product ( factor(1:n) )
!
!  LI(X) = L(X) / ( X - X(I) ) / L'(X(I))
!
  factor(1:n) = xv - x(1:n)
  factor(i) = 1.0D+00

  li = product ( factor(1:n) ) / lp
!
!  L''(X(I)).
!
  lpp = 0.0D+00
  factor(1:n) = x(i) - x(1:n)
  factor(i) = 1.0D+00

  do j = 1, n
    if ( j /= i ) then
      factor(j) = 1.0D+00
      lpp = lpp + 2.0D+00 * product ( factor )
      factor(j) = x(i) - x(j)
    end if
  end do

  value = ( 1.0D+00 - ( xv - x(i) ) * lpp / lp ) * li * li

  return
end
subroutine hermite_basis_1 ( n, x, i, xv, value )

!*****************************************************************************80
!
!! HERMITE_BASIS_1 evaluates a first-order Hermite interpolation basis function.
!
!  Discussion:
!
!    Given ND points XD, with values YD and derivative values YPD, the
!    Hermite interpolant can be written as:
!
!      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
!           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
!
!    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
!    and H1(I;X) is the I-th first order Hermite interpolation basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of abscissas.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Input, integer ( kind = 4 ) I, the index of the first-order basis function.
!
!    Input, real ( kind = 8 ) XV, the evaluation point.
!
!    Output, real ( kind = 8 ) VALUE, the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bot
  real ( kind = 8 ) factor(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) top
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xv

  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_BASIS_1 - Fatal error!'
    write ( *, '(a)' ) '  I < 1 or N < I.'
    stop
  end if

  factor(1:n) = xv - x(1:n)
  factor(i) = 1.0D+00
  top = product ( factor(1:n) )

  factor(1:n) = x(i) - x(1:n)
  factor(i) = 1.0D+00
  bot = product ( factor(1:n) )

  value = ( xv - x(i) ) * ( top / bot ) * ( top / bot )

  return
end
subroutine hermite_demo ( n, x, y, yp )

!*****************************************************************************80
!
!! HERMITE_DEMO computes and prints Hermite interpolant information for data.
!
!  Discussion:
!
!    Given a set of Hermite data, this routine calls HDATA_TO_DIF to determine
!    and print the divided difference table, and then DIF_TO_R8POLY to 
!    determine and print the coefficients of the polynomial in standard form.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Input, real ( kind = 8 ) Y(N), YP(N), the function and derivative
!    values at the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cd(0:2*n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ndp
  integer ( kind = 4 ) nv
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xd(2*n)
  real ( kind = 8 ) xdp(2*n-1)
  real ( kind = 8 ) xv(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yd(2*n)
  real ( kind = 8 ) ydp(2*n-1)
  real ( kind = 8 ) yp(n)
  real ( kind = 8 ) yv(n)
  real ( kind = 8 ) yvp(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_DEMO'
  write ( *, '(a)' ) '  Compute coefficients CD of the Hermite polynomial'
  write ( *, '(a)' ) '  interpolant to given data (x,y,yp).'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data:'
  write ( *, '(a)' ) '              X           Y           Y'''
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, x(i), y(i), yp(i)
  end do

  call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Difference table for interpolant:'
  write ( *, '(a)' ) '              XD          YD'
  write ( *, '(a)' ) ' '

  nd = 2 * n

  do i = 1, nd
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, xd(i), yd(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Difference table for derivative:'
  write ( *, '(a)' ) '              XDP         YDP'
  write ( *, '(a)' ) ' '

  ndp = 2 * n - 1

  do i = 1, ndp
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, xdp(i), ydp(i)
  end do

  call dif_to_r8poly ( nd, xd, yd, cd )

  call r8poly_print ( nd - 1, cd, '  Hermite interpolating polynomial:' )
!
!  Verify interpolation claim!
!
  nv = n
  xv(1:nv) = x(1:n)

  call hermite_interpolant_value ( nd, xd, yd, xdp, ydp, nv, xv, yv, yvp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data Versus Interpolant:'
  write ( *, '(a)' ) &
    '              X           Y           H           YP          HP'
  write ( *, '(a)' ) ' '
  do i = 1, nv
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, xv(i), y(i), yv(i), yp(i), yvp(i)
  end do

  return
end
subroutine hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT sets up a divided difference table from Hermite data.
!
!  Discussion:
!
!    The polynomial represented by the divided difference table can be
!    evaluated by calling HERMITE_INTERPOLANT_VALUE.
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
!    Input, integer ( kind = 4 ) N, the number of items of data 
!    ( X(I), Y(I), YP(I) ).
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!    These values must be distinct.
!
!    Input, real ( kind = 8 ) Y(N), YP(N), the function and 
!    derivative values.
!
!    Output, real ( kind = 8 ) XD(2*N), YD(2*N), the divided difference table
!    for the interpolant value.
!
!    Output, real ( kind = 8 ) XDP(2*N-1), YDP(2*N-1), the divided difference 
!    table for the interpolant derivative.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ndp
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xd(2*n)
  real ( kind = 8 ) xdp(2*n-1)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yd(2*n)
  real ( kind = 8 ) ydp(2*n-1)
  real ( kind = 8 ) yp(n)
!
!  Copy the data:
!
  nd = 2 * n
  xd(1:nd-1:2) = x(1:n)
  xd(2:nd  :2) = x(1:n)
!
!  Carry out the first step of differencing.
!
  yd(1) = y(1)
  yd(3:nd-1:2) = ( y(2:n) - y(1:n-1) ) / ( x(2:n) - x(1:n-1) )
  yd(2:nd  :2) = yp(1:n)  
!
!  Carry out the remaining steps in the usual way.
!
  do i = 3, nd
    do j = nd, i, -1

      yd(j) = ( yd(j) - yd(j-1) ) / ( xd(j) - xd(j+1-i) )

    end do

  end do
!
!  Compute the difference table for the derivative.
!
  call dif_deriv ( nd, xd, yd, ndp, xdp, ydp )

  return
end
subroutine hermite_interpolant_rule ( n, a, b, x, w )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT_RULE: quadrature rule for a Hermite interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of abscissas.
!
!    Input, real ( kind = 8 ) A, B, the integration limits.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(2*N), the quadrature coefficients, given as
!    pairs for function and derivative values at each abscissa.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) a_value
  real ( kind = 8 ) b
  real ( kind = 8 ) b_value
  real ( kind = 8 ) c(0:2*n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nd
  real ( kind = 8 ) w(2*n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xd(2*n)
  real ( kind = 8 ) xdp(2*n-1)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yd(2*n)
  real ( kind = 8 ) ydp(2*n-1)
  real ( kind = 8 ) yp(n)

  nd = 2 * n

  k = 0

  do i = 1, n

    k = k + 1
    y(1:n) = 0.0D+00
    y(i) = 1.0D+00
    yp(1:n) = 0.0D+00
    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
    call dif_to_r8poly ( nd, xd, yd, c )
    call r8poly_ant_val ( n, c, a, a_value )
    call r8poly_ant_val ( n, c, b, b_value )
    w(k) = b_value - a_value

    k = k + 1
    y(1:n) = 0.0D+00
    yp(1:n) = 0.0D+00
    yp(i) = 1.0D+00
    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
    call dif_to_r8poly ( nd, xd, yd, c )
    call r8poly_ant_val ( n, c, a, a_value )
    call r8poly_ant_val ( n, c, b, b_value )
    w(k) = b_value - a_value

  end do

  return
end
subroutine hermite_interpolant_value ( nd, xd, yd, xdp, ydp, nv, xv, yv, yvp )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT_VALUE evaluates the Hermite interpolant polynomial.
!
!  Discussion:
!
!    In fact, this function will evaluate an arbitrary polynomial that is
!    represented by a difference table.
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
!    Input, real ( kind = 8 ) XD(ND), YD(ND), the difference table for the
!    interpolant value.
!
!    Input, real ( kind = 8 ) XDP(ND-1), YDP(ND-1), the difference table for
!    the interpolant derivative.
!
!    Input, integer ( kind = 4 ) NV, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XV(NV), the evaluation points.
!
!    Output, real ( kind = 8 ) YV(NV), YVP(NV), the value of the interpolant and
!    its derivative at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nv

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ndp
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xdp(nd-1)
  real ( kind = 8 ) xv(nv)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) ydp(nd-1)
  real ( kind = 8 ) yv(nv)
  real ( kind = 8 ) yvp(nv)

  yv(1:nv) = yd(nd)
  do i = nd - 1, 1, -1
    yv(1:nv) = yd(i) + ( xv(1:nv) - xd(i) ) * yv(1:nv)
  end do

  ndp = nd - 1

  yvp(1:nv) = ydp(ndp)
  do i = ndp - 1, 1, -1
    yvp(1:nv) = ydp(i) + ( xv(1:nv) - xdp(i) ) * yvp(1:nv)
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
subroutine r8poly_degree ( na, a, degree )

!*****************************************************************************80
!
!! R8POLY_DEGREE returns the degree of a polynomial.
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) DEGREE, the degree of A.
!
  implicit none

  integer ( kind = 4 ) na

  real ( kind = 8 ) a(0:na)
  integer ( kind = 4 ) degree

  degree = na

  do while ( 0 < degree )

    if ( a(degree) /= 0.0D+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X^N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  call r8poly_degree ( n, a, n2 )

  if ( n2 <= 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( a(n2) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 2 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2 - 1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8vec_chebyshev ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_CHEBYSHEV creates a vector of Chebyshev spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) A(N), a vector of Chebyshev spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n

      theta = real ( n - i, kind = 8 ) * pi &
            / real ( n - 1, kind = 8 )

      c = cos ( theta )

      a(i) = ( ( 1.0D+00 - c ) * a_first  &
             + ( 1.0D+00 + c ) * a_last ) &
             /   2.0D+00

    end do

  end if

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * a_first &
             + real (     i - 1, kind = 8 ) * a_last ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
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
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
