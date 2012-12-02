subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a9to99.txt'     'a0to00.txt'
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer ( kind = 4 ) change
  integer ( kind = 4 ) digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_INC - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  if ( change == 0 ) then
    file_name = ' '
    return
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
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
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the
!    integer value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
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
subroutine ncc_triangle_degree ( rule, degree )

!*****************************************************************************80
!
!! NCC_TRIANGLE_DEGREE returns the degree of an NCC rule for the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Output, integer ( kind = 4 ) DEGREE, the polynomial degree of exactness of
!    the rule.
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) rule

  if ( 1 <= rule .and. rule <= 9 ) then

    degree = rule - 1

  else

    degree = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCC_TRIANGLE_DEGREE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine ncc_triangle_order_num ( rule, order_num )

!*****************************************************************************80
!
!! NCC_TRIANGLE_ORDER_NUM returns the order of an NCC rule for the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Output, integer ( kind = 4 ) ORDER_NUM, the order (number of points)
!    of the rule.
!
  implicit none

  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
  integer ( kind = 4 ) suborder_num

  call ncc_triangle_suborder_num ( rule, suborder_num )

  allocate ( suborder(1:suborder_num) )

  call ncc_triangle_suborder ( rule, suborder_num, suborder )

  order_num = sum ( suborder(1:suborder_num) )

  deallocate ( suborder )

  return
end
subroutine ncc_triangle_rule ( rule, order_num, xy, w )

!*****************************************************************************80
!
!! NCC_TRIANGLE_RULE returns the points and weights of an NCC rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Input, integer ( kind = 4 ) ORDER_NUM, the order (number of points)
!    of the rule.
!
!    Output, real ( kind = 8 ) XY(2,ORDER_NUM), the points of the rule.
!
!    Output, real ( kind = 8 ) W(ORDER_NUM), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) order_num

  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) k
  integer ( kind = 4 ) o
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) s
  integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
  integer ( kind = 4 ) suborder_num
  real    ( kind = 8 ), allocatable, dimension ( : ) :: suborder_w
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyz
  real    ( kind = 8 ) w(order_num)
  real    ( kind = 8 ) xy(2,order_num)
!
!  Get the suborder information.
!
  call ncc_triangle_suborder_num ( rule, suborder_num )

  allocate ( suborder(suborder_num) )
  allocate ( suborder_xyz(3,suborder_num) )
  allocate ( suborder_w(suborder_num) )

  call ncc_triangle_suborder ( rule, suborder_num, suborder )

  call ncc_triangle_subrule ( rule, suborder_num, suborder_xyz, suborder_w )
!
!  Expand the suborder information to a full order rule.
!
  o = 0

  do s = 1, suborder_num

    if ( suborder(s) == 1 ) then

      o = o + 1
      xy(1:2,o) = suborder_xyz(1:2,s)
      w(o) = suborder_w(s)

    else if ( suborder(s) == 3 ) then

      do k = 1, 3
        o = o + 1
        xy(1,o) = suborder_xyz ( i4_wrap(k,  1,3), s )
        xy(2,o) = suborder_xyz ( i4_wrap(k+1,1,3), s )
        w(o) = suborder_w(s)
      end do

    else if ( suborder(s) == 6 ) then

      do k = 1, 3
        o = o + 1
        xy(1,o) = suborder_xyz ( i4_wrap(k,  1,3), s )
        xy(2,o) = suborder_xyz ( i4_wrap(k+1,1,3), s )
        w(o) = suborder_w(s)
      end do

      do k = 1, 3
        o = o + 1
        xy(1,o) = suborder_xyz ( i4_wrap(k+1,1,3), s )
        xy(2,o) = suborder_xyz ( i4_wrap(k,  1,3), s )
        w(o) = suborder_w(s)
      end do

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NCC_TRIANGLE_RULE - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Illegal SUBORDER(', s, ') = ', suborder(s)
      write ( *, '(a,i8)' ) '  RULE =    ', rule
      write ( *, '(a,i8)' ) '  ORDER_NUM = ', order_num
      stop

    end if

  end do

  deallocate ( suborder )
  deallocate ( suborder_xyz )
  deallocate ( suborder_w )

  return
end
subroutine ncc_triangle_rule_num ( rule_num )

!*****************************************************************************80
!
!! NCC_TRIANGLE_RULE_NUM returns the number of NCC rules available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) RULE_NUM, the number of rules available.
!
  implicit none

  integer ( kind = 4 ) rule_num

  rule_num = 9

  return
end
subroutine ncc_triangle_suborder ( rule, suborder_num, suborder )

!*****************************************************************************80
!
!! NCC_TRIANGLE_SUBORDER returns the suborders for an NCC rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
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

  if ( rule == 1 ) then
    suborder(1:suborder_num) = (/ &
      1 /)
  else if ( rule == 2 ) then
    suborder(1:suborder_num) = (/ &
      3 /)
  else if ( rule == 3 ) then
    suborder(1:suborder_num) = (/ &
      3 /)
  else if ( rule == 4 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 1 /)
  else if ( rule == 5 ) then
    suborder(1:suborder_num) = (/ &
      6, 3, 3 /)
  else if ( rule == 6 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 6, 3, 3 /)
  else if ( rule == 7 ) then
    suborder(1:suborder_num) = (/ &
      6, 6, 3, 3, 6, 1 /)
  else if ( rule == 8 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 6, 3, 6, 6, 3, 3 /)
  else if ( rule == 9 ) then
    suborder(1:suborder_num) = (/ &
      6, 6, 3, 6, 6, 3, 6, 3, 3 /)
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCC_TRIANGLE_SUBORDER - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine ncc_triangle_suborder_num ( rule, suborder_num )

!*****************************************************************************80
!
!! NCC_TRIANGLE_SUBORDER_NUM returns the number of suborders for an NCC rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
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
  integer ( kind = 4 ), dimension(1:9) :: suborder = (/ &
     1, 1, 1, 3, 3, 5, 6, 8, 9 /)

  integer ( kind = 4 ) suborder_num

  if ( 1 <= rule .and. rule <= 9 ) then
    suborder_num = suborder(rule)
  else
    suborder_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCC_TRIANGLE_SUBORDER_NUM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine ncc_triangle_subrule ( rule, suborder_num, suborder_xyz, suborder_w )

!*****************************************************************************80
!
!! NCC_TRIANGLE_SUBRULE returns a compressed NCC rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Input, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders
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

  integer ( kind = 4 ), parameter :: i4_3 = 3
  integer ( kind = 4 ) rule
  real    ( kind = 8 ) suborder_w(suborder_num)
  integer ( kind = 4 ) suborder_w_n(suborder_num)
  integer ( kind = 4 ) suborder_w_d
  real    ( kind = 8 ) suborder_xyz(3,suborder_num)
  integer ( kind = 4 ) suborder_xyz_n(3,suborder_num)
  integer ( kind = 4 ) suborder_xyz_d

  if ( rule == 1 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      1,  1, 1  &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 3

    suborder_w_n(1:suborder_num) = (/ &
      1 /)

    suborder_w_d = 1

  else if ( rule == 2 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      1, 0, 0  &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 1

    suborder_w_n(1:suborder_num) = (/ &
      1 /)

    suborder_w_d = 3

  else if ( rule == 3 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      1, 1, 0  &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 2

    suborder_w_n(1:suborder_num) = (/ &
      1 /)

    suborder_w_d = 3

  else if ( rule == 4 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      3, 0, 0,  &
      2, 1, 0,  &
      1, 1, 1   &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 3

    suborder_w_n(1:suborder_num) = (/ &
      4, 9, 54 /)

    suborder_w_d = 120

  else if ( rule == 5 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      3, 1, 0,  &
      2, 2, 0,  &
      2, 1, 1   &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 4

    suborder_w_n(1:suborder_num) = (/ &
      4, -1, 8 /)

    suborder_w_d = 45

  else if ( rule == 6 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      5, 0, 0,  &
      4, 1, 0,  &
      3, 2, 0,  &
      3, 1, 1,  &
      2, 2, 1  &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 5

    suborder_w_n(1:suborder_num) = (/ &
      11, 25, 25, 200, 25 /)

    suborder_w_d = 1008

  else if ( rule == 7 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      5, 1, 0,  &
      4, 2, 0,  &
      4, 1, 1,  &
      3, 3, 0,  &
      3, 2, 1,  &
      2, 2, 2   &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 6

    suborder_w_n(1:suborder_num) = (/ &
      36, -27, 72, 64, 72, -54 /)

    suborder_w_d = 840

  else if ( rule == 8 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      7, 0, 0,  &
      6, 1, 0,  &
      5, 2, 0,  &
      5, 1, 1,  &
      4, 3, 0,  &
      4, 2, 1,  &
      3, 3, 1,  &
      3, 2, 2   &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 7

    suborder_w_n(1:suborder_num) = (/ &
      1336, 2989, 3577, 32242, 2695, -6860, 44590, 3430 /)

    suborder_w_d = 259200

  else if ( rule == 9 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      7, 1, 0,  &
      6, 2, 0,  &
      6, 1, 1,  &
      5, 3, 0,  &
      5, 2, 1,  &
      4, 4, 0,  &
      4, 3, 1,  &
      4, 2, 2,  &
      3, 3, 2   &
      /), (/ i4_3, suborder_num /) )

    suborder_xyz_d = 8

    suborder_w_n(1:suborder_num) = (/ &
      368, -468, 704, 1136, 832, -1083, 672, -1448, 1472 /)

    suborder_w_d = 14175

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCC_TRIANGLE_SUBRULE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  suborder_xyz(1:3,1:suborder_num) = &
      real ( suborder_xyz_n(1:3,1:suborder_num), kind = 8 ) &
    / real ( suborder_xyz_d,                     kind = 8 )

  suborder_w(1:suborder_num) = &
      real ( suborder_w_n(1:suborder_num), kind = 8 ) &
    / real ( suborder_w_d,                 kind = 8 )

  return
end
subroutine reference_to_physical_t3 ( node_xy, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 3 physical triangle and a point
!    (XSI,ETA) in the reference triangle, the routine computes the value
!    of the corresponding image point (X,Y) in physical space.
!
!    This routine is also appropriate for an order 4 triangle,
!    as long as the fourth node is the centroid of the triangle.
!
!    This routine may also be appropriate for an order 6
!    triangle, if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image triangle are straight and that the "midside" nodes in the
!    physical triangle are literally halfway along the sides of
!    the physical triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) NODE_XY(2,3), the coordinates of the vertices.
!    The vertices are assumed to be the images of (0,0), (1,0) and
!    (0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of objects to transform.
!
!    Input, real ( kind = 8 ) REF(2,N), points in the reference triangle.
!
!    Output, real ( kind = 8 ) PHY(2,N), corresponding points in the
!    physical triangle.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real    ( kind = 8 ) node_xy(2,3)
  real    ( kind = 8 ) phy(2,n)
  real    ( kind = 8 ) ref(2,n)

  do i = 1, 2
    phy(i,1:n) = node_xy(i,1) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) ) &
               + node_xy(i,2) *             ref(1,1:n)                &
               + node_xy(i,3) *                          ref(2,1:n)
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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = '31 May 2001   9:45:54.872 AM'
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
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = * ) string
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

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triangle_area ( node_xy, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA computes the area of a triangle.
!
!  Discussion:
!
!    If the triangle's vertices are given in counterclockwise order,
!    the area will be positive.  If the triangle's vertices are given
!    in clockwise order, the area will be negative!
!
!    If you cannot guarantee counterclockwise order, and you need to
!    have the area positive, then you can simply take the absolute value
!    of the result of this routine.
!
!    An earlier version of this routine always returned the absolute
!    value of the computed area.  I am convinced now that that is
!    a less useful result!  For instance, by returning the signed
!    area of a triangle, it is possible to easily compute the area
!    of a nonconvex polygon as the sum of the (possibly negative)
!    areas of triangles formed by node 1 and successive pairs of vertices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) NODE_XY(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) AREA, the area of the triangle.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) node_xy(2,3)

  area = 0.5D+00 * ( &
      node_xy(1,1) * ( node_xy(2,2) - node_xy(2,3) ) &
    + node_xy(1,2) * ( node_xy(2,3) - node_xy(2,1) ) &
    + node_xy(1,3) * ( node_xy(2,1) - node_xy(2,2) ) )

  return
end
subroutine triangle_points_plot ( file_name, node_xy, node_show, point_num, &
  point_xy, point_show )

!*****************************************************************************80
!
!! TRIANGLE_POINTS_PLOT plots a triangle and some points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, real ( kind = 8 ) NODE_XY(2,3), the coordinates of the nodes
!    of the triangle.
!
!    Input, integer ( kind = 4 ) NODE_SHOW,
!   -1, do not show the triangle, or the nodes.
!    0, show the triangle, do not show the nodes;
!    1, show the triangle and the nodes;
!    2, show the triangle, the nodes and number them.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates of the
!    points.
!
!    Input, integer ( kind = 4 ) POINT_SHOW,
!    0, do not show the points;
!    1, show the points;
!    2, show the points and number them.
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 3
  integer ( kind = 4 ) point_num

  character ( len = 40 ) date_time
  integer ( kind = 4 ) :: circle_size
  integer ( kind = 4 ) delta
  integer ( kind = 4 ) e
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_show
  real    ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_show
  real    ( kind = 8 ) point_xy(2,point_num)
  character ( len = 40 ) string
  real    ( kind = 8 ) x_max
  real    ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  real    ( kind = 8 ) x_scale
  real    ( kind = 8 ) y_max
  real    ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108
  real    ( kind = 8 ) y_scale

  call timestring ( date_time )
!
!  We need to do some figuring here, so that we can determine
!  the range of the data, and hence the height and width
!  of the piece of paper.
!
  x_max = max ( maxval ( node_xy(1,1:node_num) ), &
                maxval ( point_xy(1,1:point_num) ) )
  x_min = min ( minval ( node_xy(1,1:node_num) ), &
                minval ( point_xy(1,1:point_num) ) )
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = max ( maxval ( node_xy(2,1:node_num) ), &
                maxval ( point_xy(2,1:point_num) ) )
  y_min = min ( minval ( node_xy(2,1:node_num) ), &
                minval ( point_xy(2,1:point_num) ) )
  y_scale = y_max - y_min

  y_max = y_max + 0.05D+00 * y_scale
  y_min = y_min - 0.05D+00 * y_scale
  y_scale = y_max - y_min

  if ( x_scale < y_scale ) then

    delta = nint ( real ( x_ps_max - x_ps_min, kind = 8 ) &
      * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

    x_ps_max = x_ps_max - delta
    x_ps_min = x_ps_min + delta

    x_ps_max_clip = x_ps_max_clip - delta
    x_ps_min_clip = x_ps_min_clip + delta

    x_scale = y_scale

  else if ( y_scale < x_scale ) then

    delta = nint ( real ( y_ps_max - y_ps_min, kind = 8 ) &
      * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

    y_ps_max      = y_ps_max - delta
    y_ps_min      = y_ps_min + delta

    y_ps_max_clip = y_ps_max_clip - delta
    y_ps_min_clip = y_ps_min_clip + delta

    y_scale = x_scale

  end if

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_POINTS_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
    return
  end if

  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: triangulation_order3_plot.f90'
  write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( file_unit, '(a)' ) '%%CreationDate: ' // trim ( date_time )
  write ( file_unit, '(a)' ) '%%Pages: 1'
  write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( file_unit, '(a)' ) '%%EndComments'
  write ( file_unit, '(a)' ) '%%BeginProlog'
  write ( file_unit, '(a)' ) '/inch {72 mul} def'
  write ( file_unit, '(a)' ) '%%EndProlog'
  write ( file_unit, '(a)' ) '%%Page: 1 1'
  write ( file_unit, '(a)' ) 'save'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( file_unit, '(a)' ) 'stroke'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the font and its size.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '/Times-Roman findfont'
  write ( file_unit, '(a)' ) '0.50 inch scalefont'
  write ( file_unit, '(a)' ) 'setfont'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Print a title.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  210  702  moveto'
  write ( file_unit, '(a)' ) '%  (Triangulation)  show'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a)' ) 'clip newpath'
!
!  Draw the nodes.
!
  if ( 1 <= node_show ) then

    circle_size = 5

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
    write ( file_unit, '(a)' ) '%'

    do node = 1, 3

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
        circle_size, '0 360 arc closepath fill'

    end do

  end if
!
!  Label the nodes.
!
  if ( 2 <= node_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the nodes:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'
    write ( file_unit, '(a)' ) '%'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
        ' moveto (' // trim ( string ) // ') show'

    end do

  end if
!
!  Draw the points.
!
  if ( point_num <= 200 ) then
    circle_size = 5
  else if ( point_num <= 500 ) then
    circle_size = 4
  else if ( point_num <= 1000 ) then
    circle_size = 3
  else if ( point_num <= 5000 ) then
    circle_size = 2
  else
    circle_size = 1
  end if

  if ( 1 <= point_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw filled dots at the points.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to green.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.150  0.750  0.000 setrgbcolor'
    write ( file_unit, '(a)' ) '%'

    do point = 1, point_num

      x_ps = int ( &
        ( ( x_max - point_xy(1,point)         ) &
        * real ( x_ps_min, kind = 8 )   &
        + (         point_xy(1,point) - x_min ) &
        * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                     - x_min ) )

      y_ps = int ( &
        ( ( y_max - point_xy(2,point)         ) &
        * real ( y_ps_min, kind = 8 )   &
        + (         point_xy(2,point) - y_min ) &
        * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                     - y_min ) )

      write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
        circle_size, '0 360 arc closepath fill'

    end do

  end if
!
!  Label the points.
!
  if ( 2 <= point_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the point:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker green.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.250  0.850  0.000 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'
    write ( file_unit, '(a)' ) '%'

    do point = 1, point_num

      x_ps = int ( &
        ( ( x_max - point_xy(1,point)         ) &
        * real ( x_ps_min, kind = 8 )   &
        + (       + point_xy(1,point) - x_min ) &
        * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                     - x_min ) )

      y_ps = int ( &
        ( ( y_max - point_xy(2,point)         ) &
        * real ( y_ps_min, kind = 8 )   &
        + (         point_xy(2,point) - y_min ) &
        * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                     - y_min ) )

      write ( string, '(i4)' ) point
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
        ' moveto (' // trim ( string ) // ') show'

    end do

  end if
!
!  Draw the triangle.
!
  if ( 0 <= node_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw the triangle.'
    write ( file_unit, '(a)' ) '%'

    write ( file_unit, '(a)' ) 'newpath'

    do i = 1, 4

      node = i4_wrap ( i, 1, 3 )

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      if ( i == 1 ) then
        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'
      else
        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'
      end if

    end do

    write ( file_unit, '(a)' ) 'stroke'

  end if

  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'restore  showpage'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  End of page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%%Trailer'
  write ( file_unit, '(a)' ) '%%EOF'
  close ( unit = file_unit )

  return
end
