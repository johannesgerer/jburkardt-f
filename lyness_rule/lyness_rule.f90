subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine lyness_order ( rule, order )

!*****************************************************************************80
!
!! LYNESS_ORDER returns the order of a Lyness quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Output, integer ( kind = 4 ) ORDER, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) rule

  if ( rule == 0 ) then
    order = 1
  else if ( rule == 1 ) then
    order = 3
  else if ( rule == 2 ) then
    order = 4
  else if ( rule == 3 ) then
    order = 4
  else if ( rule == 4 ) then
    order = 7
  else if ( rule == 5 ) then
    order = 6
  else if ( rule == 6 ) then
    order = 10
  else if ( rule == 7 ) then
    order = 9
  else if ( rule == 8 ) then
    order = 7
  else if ( rule == 9 ) then
    order = 10
  else if ( rule == 10 ) then
    order = 12
  else if ( rule == 11 ) then
    order = 16
  else if ( rule == 12 ) then
    order = 13
  else if ( rule == 13 ) then
    order = 13
  else if ( rule == 14 ) then
    order = 16
  else if ( rule == 15 ) then
    order = 16
  else if ( rule == 16 ) then
    order = 21
  else if ( rule == 17 ) then
    order = 16
  else if ( rule == 18 ) then
    order = 19
  else if ( rule == 19 ) then
    order = 22
  else if ( rule == 20 ) then
    order = 27
  else if ( rule == 21 ) then
    order = 28
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LYNESS_ORDER - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized rule index.'
    stop
  end if

  return
end
subroutine lyness_precision ( rule, precision )

!*****************************************************************************80
!
!! LYNESS_PRECISION returns the precision of a Lyness quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Output, integer ( kind = 4 ) PRECISION, the precision of the rule.
!
  implicit none

  integer ( kind = 4 ) precision
  integer ( kind = 4 ) rule

  if ( rule == 0 ) then
    precision = 1
  else if ( rule == 1 ) then
    precision = 2
  else if ( rule == 2 ) then
    precision = 2
  else if ( rule == 3 ) then
    precision = 3
  else if ( rule == 4 ) then
    precision = 3
  else if ( rule == 5 ) then
    precision = 4
  else if ( rule == 6 ) then
    precision = 4
  else if ( rule == 7 ) then
    precision = 4
  else if ( rule == 8 ) then
    precision = 5
  else if ( rule == 9 ) then
    precision = 5
  else if ( rule == 10 ) then
    precision = 6
  else if ( rule == 11 ) then
    precision = 6
  else if ( rule == 12 ) then
    precision = 6
  else if ( rule == 13 ) then
    precision = 7
  else if ( rule == 14 ) then
    precision = 7
  else if ( rule == 15 ) then
    precision = 8
  else if ( rule == 16 ) then
    precision = 8
  else if ( rule == 17 ) then
    precision = 8
  else if ( rule == 18 ) then
    precision = 9
  else if ( rule == 19 ) then
    precision = 9
  else if ( rule == 20 ) then
    precision = 11
  else if ( rule == 21 ) then
    precision = 11
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LYNESS_PRECISION - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized rule index.'
    stop
  end if

  return
end
subroutine lyness_rule ( rule, order, w, x )

!*****************************************************************************80
!
!! LYNESS_RULE returns the points and weights of a Lyness quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
!    Output, real ( kind = 8 ) X(2,ORDER), the points.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) k
  integer ( kind = 4 ) o
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) s
  integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
  integer ( kind = 4 ) suborder_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: sub_w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sub_xyz
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(2,order)
!
!  Get the suborder information.
!
  call lyness_suborder_num ( rule, suborder_num )

  allocate ( suborder(suborder_num) )
  allocate ( sub_xyz(3,suborder_num) )
  allocate ( sub_w(suborder_num) )

  call lyness_suborder ( rule, suborder_num, suborder )

  call lyness_subrule ( rule, suborder_num, sub_xyz, sub_w )
!
!  Expand the suborder information to a full order rule.
!
  o = 0

  do s = 1, suborder_num

    if ( suborder(s) == 1 ) then

      o = o + 1
      x(1:2,o) = sub_xyz(1:2,s)
      w(o) = sub_w(s)

    else if ( suborder(s) == 3 ) then

      do k = 1, 3
        o = o + 1
        x(1,o) = sub_xyz ( i4_wrap(k,  1,3), s )
        x(2,o) = sub_xyz ( i4_wrap(k+1,1,3), s )
        w(o) = sub_w(s) / 3.0D+00
      end do

    else if ( suborder(s) == 6 ) then

      do k = 1, 3
        o = o + 1
        x(1,o) = sub_xyz ( i4_wrap(k,  1,3), s )
        x(2,o) = sub_xyz ( i4_wrap(k+1,1,3), s )
        w(o) = sub_w(s) / 6.0D+00
      end do

      do k = 1, 3
        o = o + 1
        x(1,o) = sub_xyz ( i4_wrap(k+1,1,3), s )
        x(2,o) = sub_xyz ( i4_wrap(k,  1,3), s )
        w(o) = sub_w(s) / 6.0D+00
      end do

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LYNESS_RULE - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Illegal SUBORDER(', s, ') = ', suborder(s) 
      stop

    end if

  end do

  deallocate ( suborder )
  deallocate ( sub_xyz )
  deallocate ( sub_w )

  return
end
subroutine lyness_rule_num ( rule_num )

!*****************************************************************************80
!
!! LYNESS_RULE_NUM returns the number of Lyness quadrature rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) RULE_NUM, the number of rules.
!
  implicit none

  integer ( kind = 4 ) rule_num

  rule_num = 21

  return
end
subroutine lyness_suborder ( rule, suborder_num, suborder )

!*****************************************************************************80
!
!! LYNESS_SUBORDER returns the suborders for a Lyness rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Input, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders 
!    of the rule.
!
!    Output, integer ( kind = 4 ) SUBORDER(SUBORDER_NUM), the suborders 
!    of the rule.
!
  implicit none

  integer ( kind = 4 ) suborder_num

  integer ( kind = 4 ) rule
  integer ( kind = 4 ) suborder(suborder_num)

  if ( rule == 0 ) then
    suborder(1:suborder_num) = (/ &
      1 /)
  else if ( rule == 1 ) then
    suborder(1:suborder_num) = (/ &
      3 /)
  else if ( rule == 2 ) then
    suborder(1:suborder_num) = (/ &
      1, 3 /)
  else if ( rule == 3 ) then
    suborder(1:suborder_num) = (/ &
      1, 3 /)
  else if ( rule == 4 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3 /)
  else if ( rule == 5 ) then
    suborder(1:suborder_num) = (/ &
      3, 3  /)
  else if ( rule == 6 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 6  /)
  else if ( rule == 7 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 3 /)
  else if ( rule == 8 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3 /)
  else if ( rule == 9 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3 /)
  else if ( rule == 10 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 6 /)
  else if ( rule == 11 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 6 /)
  else if ( rule == 12 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 6 /)
  else if ( rule == 13 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 6 /)
  else if ( rule == 14 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 6 /)
  else if ( rule == 15 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 6 /)
  else if ( rule == 16 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 3, 3, 3, 6 /)
  else if ( rule == 17 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 6 /)
  else if ( rule == 18 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 6 /)
  else if ( rule == 19 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 3, 6 /)
  else if ( rule == 20 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 3, 3, 3, 6, 6 /)
  else if ( rule == 21 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 3, 6, 6 /)
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LYNESS_SUBORDER - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine lyness_suborder_num ( rule, suborder_num )

!*****************************************************************************80
!
!! LYNESS_SUBORDER_NUM returns the number of suborders for a Lyness rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Output, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders 
!    of the rule.
!
  implicit none

  integer ( kind = 4 ) rule
  integer ( kind = 4 ) suborder_num

  if ( rule == 0 ) then
    suborder_num = 1
  else if ( rule == 1 ) then
    suborder_num = 1
  else if ( rule == 2 ) then
    suborder_num = 2
  else if ( rule == 3 ) then
    suborder_num = 2
  else if ( rule == 4 ) then
    suborder_num = 3
  else if ( rule == 5 ) then
    suborder_num = 2
  else if ( rule == 6 ) then
    suborder_num = 3
  else if ( rule == 7 ) then
    suborder_num = 3
  else if ( rule == 8 ) then
    suborder_num = 3
  else if ( rule == 9 ) then
    suborder_num = 4
  else if ( rule == 10 ) then
    suborder_num = 3
  else if ( rule == 11 ) then
    suborder_num = 5
  else if ( rule == 12 ) then
    suborder_num = 4
  else if ( rule == 13 ) then
    suborder_num = 4
  else if ( rule == 14 ) then
    suborder_num = 5
  else if ( rule == 15 ) then
    suborder_num = 5
  else if ( rule == 16 ) then
    suborder_num = 6
  else if ( rule == 17 ) then
    suborder_num = 5
  else if ( rule == 18 ) then
    suborder_num = 6
  else if ( rule == 19 ) then
    suborder_num = 7
  else if ( rule == 20 ) then
    suborder_num = 7
  else if ( rule == 21 ) then
    suborder_num = 8
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LYNESS_SUBORDER_NUM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine lyness_subrule ( rule, suborder_num, sub_xyz, sub_w )

!*****************************************************************************80
!
!! LYNESS_SUBRULE returns a compressed Lyness rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Input, integer ( kind = 4 ) SUB_ORDER_NUM, the number of suborders 
!    of the rule.
!
!    Output, real ( kind = 8 ) SUBORDER_XYZ(3,SUBORDER_NUM),
!    the barycentric coordinates of the abscissas.
!
!    Output, real ( kind = 8 ) SUBORDER_W(SUBORDER_NUM), the
!    suborder weights.
!
  implicit none

  integer ( kind = 4 ) suborder_num

  integer ( kind = 4 ) rule
  integer ( kind = 4 ) s
  real ( kind = 8 ) sub_w(suborder_num)
  real ( kind = 8 ) sub_xyz(3,suborder_num)

  if ( rule == 0 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      1.0000000000D+00 /)

  else if ( rule == 1 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.0000000000D+00,  0.5000000000D+00, 0.5000000000D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      1.0000000000D+00 /)

  else if ( rule == 2 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
      1.0000000000D+00,  0.0000000000D+00, 0.0000000000D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      0.7500000000D+00, &
      0.2500000000D+00 /)

  else if ( rule == 3 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
      0.6000000000D+00,  0.2000000000D+00, 0.2000000000D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      - 0.5625000000D+00, &
        1.5625000000D+00 /)

  else if ( rule == 4 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.33333333333333333D+00,  0.33333333333333333D+00, 0.33333333333333333D+00, &
      1.00000000000000000D+00,  0.00000000000000000D+00, 0.00000000000000000D+00, &
      0.00000000000000000D+00,  0.50000000000000000D+00, 0.50000000000000000D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      0.45D+00, &
      0.15D+00, &
      0.40D+00 /)

  else if ( rule == 5 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      8.168475729804585D-01, 9.157621350977073D-02, 9.15762135097707569D-02, &
      1.081030181680702D-01, 4.459484909159649D-01, 0.44594849091596489D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      3.298552309659655D-01, &
      6.701447690340345D-01 /)

  else if ( rule == 6 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.33333333333333333D+00,  0.33333333333333333D+00, 0.33333333333333333D+00, &
      1.00000000000000000D+00,  0.00000000000000000D+00, 0.00000000000000000D+00, &
      0.00000000000000000D+00,  0.78867513459481281D+00, 0.21132486540518719D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      0.45D+00, &
    - 0.05D+00, &
      0.60D+00 /)

  else if ( rule == 7 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      1.00000000000000000D+00,  0.00000000000000000D+00, 0.00000000000000000D+00, &
      0.00000000000000000D+00,  0.50000000000000000D+00, 0.50000000000000000D+00, &
      0.62283903060710999D+00,  0.18858048469644506D+00, 0.18858048469644506D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      6.16204060378000851D-02, &
      0.18592649660480146D+00, &     
      0.75245309735739840D+00 /)

  else if ( rule == 8 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.33333333333333333D+00,  0.33333333333333333D+00, 0.33333333333333333D+00, &
      0.79742698535308720D+00,  0.10128650732345633D+00, 0.10128650732345633D+00, &
      5.97158717897698088D-02,  0.47014206410511505D+00, 0.47014206410511505D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      0.22500000000000000D+00, &
      0.37781754163448150D+00, &     
      0.39718245836551852D+00 /)

  else if ( rule == 9 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.33333333333333333D+00,  0.33333333333333333D+00, 0.33333333333333333D+00, &
      1.00000000000000000D+00,  0.00000000000000000D+00, 0.00000000000000000D+00, &
      0.00000000000000000D+00,  0.50000000000000000D+00, 0.50000000000000000D+00, &
      0.71428571428571430D-00,  0.14285714285714285D+00, 0.14285714285714285D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      0.25312500000000000D+00, &
      0.03333333333333333D+00, &
      0.21333333333333333D+00, &     
      0.50020833333333333D+00 /)

  else if ( rule == 10 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      5.014265096581342D-01,  2.492867451709329D-01, 0.24928674517093291D+00, &
      8.738219710169965D-01,  6.308901449150177D-02, 6.30890144915016854D-02, &
      6.365024991213939D-01,  5.314504984483216D-02, 0.31035245103377396D+00 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      3.503588271790222D-01, &
      1.525347191106164D-01, &     
      4.971064537103575D-01 /)

  else if ( rule == 11 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.33333333333333333D+00,  0.33333333333333333D+00, 0.33333333333333333D+00, &
      1.00000000000000000D+00,  0.00000000000000000D+00, 0.00000000000000000D+00, &
      0.00000000000000000D+00,  0.50000000000000000D+00, 0.50000000000000000D+00, &
      0.50000000000000000D+00,  0.25000000000000000D+00, 0.25000000000000000D+00, &
      0.00000000000000000D+00,  0.90824829046386302D+00, 9.17517095361370244D-02 &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      -0.57857142857142863D+00, &   
      -5.95238095238095205D-02, &
       0.16190476190476191D+00, &
       1.2190476190476192D+00, & 
       0.25714285714285712D+00 /)

  else if ( rule == 12 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      3.333333333333333D-01,  3.333333333333333D-01,  0.33333333333333333D+00, &
      5.233837209269747D-02,  4.738308139536513D-01,  0.47383081395365129D+00, &
      6.557646607383649D-01,  1.721176696308175D-01,  0.17211766963081757D+00, &
      0.000000000000000D+00,  8.653073540834571D-01,  0.13469264591654295D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      1.527089667883523D-01, &
      2.944076042366762D-01, &
      3.887052878418766D-01, &
      1.641781411330949D-01 /)

  else if ( rule == 13 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      3.333333333333333D-01,  3.333333333333333D-01,  0.33333333333333333D+00, &
      4.793080678419067D-01,  2.603459660790466D-01,  0.26034596607904670D+00, & 
      8.697397941955675D-01,  6.513010290221623D-02,  6.51301029022163108D-02, &
      6.384441885698096D-01,  4.869031542531756D-02,  0.31286549600487290D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      - 1.495700444677495D-01, &
        5.268457722996828D-01, &
        1.600417068265167D-01, &
        4.626825653415500D-01 /)

  else if ( rule == 14 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      3.333333333333333D-01,  3.333333333333333D-01,  0.33333333333333333D+00, &
      1.000000000000000D+00,  0.000000000000000D+00,  0.00000000000000000D+00, &
      6.901278795524791D-01,  1.549360602237604D-01,  0.15493606022376047D+00, &
      6.169850771237593D-02,  4.691507461438120D-01,  0.46915074614381203D+00, &
      0.000000000000000D+00,  8.392991722729236D-01,  0.16070082772707639D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      1.763126156005252D-01, &
      1.210901532768310D-02, &
      3.499561757697094D-01, &
      3.195119754425220D-01, &
      1.421102178595603D-01 /)
       
  else if ( rule == 15 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      3.333333333333333D-01,  3.333333333333333D-01,  0.33333333333333333D+00, &
      8.141482341455413D-02,  4.592925882927229D-01,  0.45929258829272301D+00, &
      8.989055433659379D-01,  5.054722831703103D-02,  5.05472283170311024D-02, &
      6.588613844964797D-01,  1.705693077517601D-01,  0.17056930775176021D+00, &
      8.394777409957211D-03,  7.284923929554041D-01,  0.26311282963463867D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      1.443156076777862D-01, &
      2.852749028018549D-01, &
      9.737549286959440D-02, &
      3.096521116041552D-01, &
      1.633818850466092D-01 /)

  else if ( rule == 16 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      1.000000000000000D+00,  0.000000000000000D+00,  0.00000000000000000D+00, &
      0.000000000000000D+00,  0.500000000000000D+00,  0.50000000000000000D+00, &
      8.637211648883667D-03,  4.956813941755582D-01,  0.49568139417555818D+00, &
      8.193444849714693D-01,  9.032775751426533D-02,  9.03277575142653888D-02, &
      5.316905005853895D-01,  2.341547497073052D-01,  0.23415474970730532D+00, &
      0.000000000000000D-01,  7.236067977499790D-01,  0.27639320225002095D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      1.207273935292775D-02, &
     -8.491579879151455D-01, &
      1.042367468891334D+00, &
      1.947229791412260D-01, &
      4.511852767201322D-01, &
      1.488095238095238D-01 /)

  else if ( rule == 17 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      3.333333333333333D-01,  3.333333333333333D-01,  0.33333333333333333D+00, &
      4.666912123569507D-02,  4.766654393821525D-01,  0.47666543938215239D+00, &
      9.324563118910393D-01,  3.377184405448033D-02,  3.37718440544803877D-02, &
      4.593042216691921D-01,  2.703478891654040D-01,  0.27034788916540398D+00, &
      5.146433548666149D-02,  7.458294907672514D-01,  0.20270617374608713D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      - 2.834183851113958D-01, &
        2.097208857979572D-01, &
        5.127273801480265D-02, &
        6.564896469913506D-01, &
        3.659351143072855D-01 /)

  else if ( rule == 18 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      3.333333333333333D-01,  3.333333333333333D-01,  0.33333333333333333D+00, &
      2.063496160252593D-02,  4.896825191987370D-01,  0.48968251919873701D+00, &
      1.258208170141290D-01,  4.370895914929355D-01,  0.43708959149293541D+00, &
      6.235929287619356D-01,  1.882035356190322D-01,  0.18820353561903219D+00, &
      9.105409732110941D-01,  4.472951339445297D-02,  4.47295133944529688D-02, &
      3.683841205473626D-02,  7.411985987844980D-01,  0.22196298916076573D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      9.713579628279610D-02, &
      9.400410068141950D-02, &
      2.334826230143263D-01, &
      2.389432167816273D-01, &
      7.673302697609430D-02, &
      2.597012362637364D-01 /)

  else if ( rule == 19 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      3.333333333333333D-01,  3.333333333333333D-01,  0.33333333333333333D+00, &
      1.000000000000000D+00,  0.000000000000000D+00,  0.00000000000000000D+00, &
      0.000000000000000D+00,  0.500000000000000D+00,  0.50000000000000000D+00, &
      1.004413236259677D-01,  4.497793381870162D-01,  0.44977933818701610D+00, &
      9.061051136018193D-01,  4.694744319909033D-02,  4.69474431990903329D-02, &
      6.162561745251021D-01,  1.918719127374489D-01,  0.19187191273744902D+00, &
      3.683841205473626D-02,  7.411985987844980D-01,  0.22196298916076573D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      1.133624844599192D-01, &
      1.062573789846380D-03, &
      4.803411513859279D-02, &
      2.524243006337300D-01, &
      7.819254371487040D-02, &
      2.472227459993048D-01, &
      2.597012362637364D-01 /)

  else if ( rule == 20 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      9.352701037774565D-01,  3.236494811127173D-02,  3.23649481112718157D-02, &
      7.612981754348137D-01,  1.193509122825931D-01,  0.11935091228259319D+00, &
    - 6.922209654151433D-02,  5.346110482707572D-01,  0.53461104827075701D+00, &
      5.933801991374367D-01,  2.033099004312816D-01,  0.20330990043128172D+00, &
      2.020613940682885D-01,  3.989693029658558D-01,  0.39896930296585570D+00, &
      5.017813831049474D-02,  5.932012134282132D-01,  0.35662064826129203D+00, &
      2.102201653616613D-02,  8.074890031597923D-01,  0.17148898030404158D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      4.097919300803106D-02, &
      1.085536215102866D-01, &
      2.781018986881812D-03, &
      1.779689321422668D-01, &
      2.314486047444677D-01, &
      3.140226717732234D-01, &
      1.242459578348437D-01 /)

  else if ( rule == 21 ) then

    sub_xyz(1:3,1:suborder_num) = reshape ( (/ &
      3.333333333333333D-01,  3.333333333333333D-01,  0.33333333333333333D+00, &
      9.480217181434233D-01,  2.598914092828833D-02,  2.59891409282883845D-02, &
      8.114249947041546D-01,  9.428750264792270D-02,  9.42875026479226691D-02, &
      1.072644996557060D-02,  4.946367750172147D-01,  0.49463677501721470D+00, &
      5.853132347709715D-01,  2.073433826145142D-01,  0.20734338261451427D+00, &
      1.221843885990187D-01,  4.389078057004907D-01,  0.43890780570049059D+00, &
      0.000000000000000D+00,  8.588702812826364D-01,  0.14112971871736357D+00, &
      4.484167758913055D-02,  6.779376548825902D-01,  0.27722066752827923D+00  &
    /), (/ 3, suborder_num /) )

    sub_w(1:suborder_num) = (/ &
      8.797730116222190D-02, &
      2.623293466120857D-02, &
      1.142447159818060D-01, &
      5.656634416839376D-02, &
      2.164790926342230D-01, &
      2.079874161166116D-01, &
      4.417430269980344D-02, &
      2.463378925757316D-01 /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LYNESS_SUBRULE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real      ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

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

  character ( len = 8 ) ampm
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
