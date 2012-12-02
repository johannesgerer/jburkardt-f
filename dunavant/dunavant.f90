subroutine dunavant_degree ( rule, degree )

!*****************************************************************************80
!
!! DUNAVANT_DEGREE returns the degree of a Dunavant rule for the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Dunavant,
!    High Degree Efficient Symmetrical Gaussian Quadrature Rules
!    for the Triangle,
!    International Journal for Numerical Methods in Engineering,
!    Volume 21, 1985, pages 1129-1148.
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
!    Output, integer ( kind = 4 ) DEGREE, the polynomial degree of exactness of
!    the rule.
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) rule

  if ( 1 <= rule .and. rule <= 20 ) then

    degree = rule

  else

    degree = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DUNAVANT_DEGREE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine dunavant_order_num ( rule, order_num )

!*****************************************************************************80
!
!! DUNAVANT_ORDER_NUM returns the order of a Dunavant rule for the triangle.
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
!    David Dunavant,
!    High Degree Efficient Symmetrical Gaussian Quadrature Rules
!    for the Triangle,
!    International Journal for Numerical Methods in Engineering,
!    Volume 21, 1985, pages 1129-1148.
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
!    Output, integer ( kind = 4 ) ORDER_NUM, the order (number of points)
!    of the rule.
!
  implicit none

  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
  integer ( kind = 4 ) suborder_num

  call dunavant_suborder_num ( rule, suborder_num )

  allocate ( suborder(1:suborder_num) )

  call dunavant_suborder ( rule, suborder_num, suborder )

  order_num = sum ( suborder(1:suborder_num) )

  deallocate ( suborder )

  return
end
subroutine dunavant_rule ( rule, order_num, xy, w )

!*****************************************************************************80
!
!! DUNAVANT_RULE returns the points and weights of a Dunavant rule.
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
!    David Dunavant,
!    High Degree Efficient Symmetrical Gaussian Quadrature Rules
!    for the Triangle,
!    International Journal for Numerical Methods in Engineering,
!    Volume 21, 1985, pages 1129-1148.
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
  real ( kind = 8 ), allocatable, dimension ( : ) :: suborder_w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyz
  real ( kind = 8 ) w(order_num)
  real ( kind = 8 ) xy(2,order_num)
!
!  Get the suborder information.
!
  call dunavant_suborder_num ( rule, suborder_num )

  allocate ( suborder(suborder_num) )
  allocate ( suborder_xyz(3,suborder_num) )
  allocate ( suborder_w(suborder_num) )

  call dunavant_suborder ( rule, suborder_num, suborder )

  call dunavant_subrule ( rule, suborder_num, suborder_xyz, suborder_w )
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
      write ( *, '(a)' ) 'DUNAVANT_RULE - Fatal error!'
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
subroutine dunavant_rule_num ( rule_num )

!*****************************************************************************80
!
!! DUNAVANT_RULE_NUM returns the number of Dunavant rules available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Dunavant,
!    High Degree Efficient Symmetrical Gaussian Quadrature Rules
!    for the Triangle,
!    International Journal for Numerical Methods in Engineering,
!    Volume 21, 1985, pages 1129-1148.
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) RULE_NUM, the number of rules available.
!
  implicit none

  integer ( kind = 4 ) rule_num

  rule_num = 20

  return
end
subroutine dunavant_suborder ( rule, suborder_num, suborder )

!*****************************************************************************80
!
!! DUNAVANT_SUBORDER returns the suborders for a Dunavant rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Dunavant,
!    High Degree Efficient Symmetrical Gaussian Quadrature Rules
!    for the Triangle,
!    International Journal for Numerical Methods in Engineering,
!    Volume 21, 1985, pages 1129-1148.
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

  if ( rule == 1 ) then
    suborder(1:suborder_num) = (/ &
      1 /)
  else if ( rule == 2 ) then
    suborder(1:suborder_num) = (/ &
      3 /)
  else if ( rule == 3 ) then
    suborder(1:suborder_num) = (/ &
      1, 3 /)
  else if ( rule == 4 ) then
    suborder(1:suborder_num) = (/ &
      3, 3 /)
  else if ( rule == 5 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3 /)
  else if ( rule == 6 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 6 /)
  else if ( rule == 7 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 6 /)
  else if ( rule == 8 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 6 /)
  else if ( rule == 9 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 6 /)
  else if ( rule == 10 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 6, 6, 6 /)
  else if ( rule == 11 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 3, 3, 3, 6, 6 /)
  else if ( rule == 12 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 3, 3, 3, 6, 6, 6 /)
  else if ( rule == 13 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 3, 3, 6, 6, 6 /)
  else if ( rule == 14 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 3, 3, 3, 3, 6, 6, 6, 6 /)
  else if ( rule == 15 ) then
    suborder(1:suborder_num) = (/ &
      3, 3, 3, 3, 3, 3, 6, 6, 6, 6, &
      6 /)
  else if ( rule == 16 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 3, 3, 3, 6, 6, &
      6, 6, 6 /)
  else if ( rule == 17 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 3, 3, 3, 3, 6, &
      6, 6, 6, 6, 6 /)
  else if ( rule == 18 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
      6, 6, 6, 6, 6, 6, 6 /)
  else if ( rule == 19 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 3, 3, 3, 3, 6, &
      6, 6, 6, 6, 6, 6, 6 /)
  else if ( rule == 20 ) then
    suborder(1:suborder_num) = (/ &
      1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
      3, 6, 6, 6, 6, 6, 6, 6, 6 /)
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DUNAVANT_SUBORDER - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine dunavant_suborder_num ( rule, suborder_num )

!*****************************************************************************80
!
!! DUNAVANT_SUBORDER_NUM returns the number of suborders for a Dunavant rule.
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
!    David Dunavant,
!    High Degree Efficient Symmetrical Gaussian Quadrature Rules
!    for the Triangle,
!    International Journal for Numerical Methods in Engineering,
!    Volume 21, 1985, pages 1129-1148.
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
  integer ( kind = 4 ), dimension(1:20) :: suborder = (/ &
     1,  1 , 2,  2,  3,  3,  4,  5,  6,  6, &
     7,  8, 10, 10, 11, 13, 15, 17, 17, 19 /)

  integer ( kind = 4 ) suborder_num

  if ( 1 <= rule .and. rule <= 20 ) then
    suborder_num = suborder(rule)
  else
    suborder_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DUNAVANT_SUBORDER_NUM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
subroutine dunavant_subrule ( rule, suborder_num, suborder_xyz, suborder_w )

!*****************************************************************************80
!
!! DUNAVANT_SUBRULE returns a compressed Dunavant rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Dunavant,
!    High Degree Efficient Symmetrical Gaussian Quadrature Rules
!    for the Triangle,
!    International Journal for Numerical Methods in Engineering,
!    Volume 21, 1985, pages 1129-1148.
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
!    Output, real ( kind = 8 ) SUBORDER_XYZ(3,SUBORDER_NUM),
!    the barycentric coordinates of the abscissas.
!
!    Output, real ( kind = 8 ) SUBORDER_W(SUBORDER_NUM), the
!    suborder weights.
!
  implicit none

  integer ( kind = 4 ) suborder_num

  integer ( kind = 4 ) rule
  real ( kind = 8 ) suborder_w(suborder_num)
  real ( kind = 8 ) suborder_xyz(3,suborder_num)

  if ( rule == 1 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00,  0.333333333333333D+00, 0.333333333333333D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      1.000000000000000D+00 /)

  else if ( rule == 2 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.666666666666667D+00, 0.166666666666667D+00, 0.166666666666667D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.333333333333333D+00 /)

  else if ( rule == 3 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.600000000000000D+00, 0.200000000000000D+00, 0.200000000000000D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      -0.562500000000000D+00, &
       0.520833333333333D+00 /)

  else if ( rule == 4 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.108103018168070D+00, 0.445948490915965D+00, 0.445948490915965D+00, &
      0.816847572980459D+00, 0.091576213509771D+00, 0.091576213509771D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.223381589678011D+00, &
      0.109951743655322D+00 /)

  else if ( rule == 5 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.059715871789770D+00, 0.470142064105115D+00, 0.470142064105115D+00, &
      0.797426985353087D+00, 0.101286507323456D+00, 0.101286507323456D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.225000000000000D+00, &
      0.132394152788506D+00, &
      0.125939180544827D+00 /)

  else if ( rule == 6 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.501426509658179D+00, 0.249286745170910D+00, 0.249286745170910D+00, &
      0.873821971016996D+00, 0.063089014491502D+00, 0.063089014491502D+00, &
      0.053145049844817D+00, 0.310352451033784D+00, 0.636502499121399D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.116786275726379D+00, &
      0.050844906370207D+00, &
      0.082851075618374D+00 /)

  else if ( rule == 7 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.479308067841920D+00, 0.260345966079040D+00, 0.260345966079040D+00, &
      0.869739794195568D+00, 0.065130102902216D+00, 0.065130102902216D+00, &
      0.048690315425316D+00, 0.312865496004874D+00, 0.638444188569810D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
     -0.149570044467682D+00, &
      0.175615257433208D+00, &
      0.053347235608838D+00, &
      0.077113760890257D+00 /)

  else if ( rule == 8 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.081414823414554D+00, 0.459292588292723D+00, 0.459292588292723D+00, &
      0.658861384496480D+00, 0.170569307751760D+00, 0.170569307751760D+00, &
      0.898905543365938D+00, 0.050547228317031D+00, 0.050547228317031D+00, &
      0.008394777409958D+00, 0.263112829634638D+00, 0.728492392955404D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.144315607677787D+00, &
      0.095091634267285D+00, &
      0.103217370534718D+00, &
      0.032458497623198D+00, &
      0.027230314174435D+00 /)

  else if ( rule == 9 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.020634961602525D+00, 0.489682519198738D+00, 0.489682519198738D+00, &
      0.125820817014127D+00, 0.437089591492937D+00, 0.437089591492937D+00, &
      0.623592928761935D+00, 0.188203535619033D+00, 0.188203535619033D+00, &
      0.910540973211095D+00, 0.044729513394453D+00, 0.044729513394453D+00, &
      0.036838412054736D+00, 0.221962989160766D+00, 0.741198598784498D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.097135796282799D+00, &
      0.031334700227139D+00, &
      0.077827541004774D+00, &
      0.079647738927210D+00, &
      0.025577675658698D+00, &
      0.043283539377289D+00 /)

  else if ( rule == 10 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.028844733232685D+00, 0.485577633383657D+00, 0.485577633383657D+00, &
      0.781036849029926D+00, 0.109481575485037D+00, 0.109481575485037D+00, &
      0.141707219414880D+00, 0.307939838764121D+00, 0.550352941820999D+00, &
      0.025003534762686D+00, 0.246672560639903D+00, 0.728323904597411D+00, &
      0.009540815400299D+00, 0.066803251012200D+00, 0.923655933587500D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.090817990382754D+00, &
      0.036725957756467D+00, &
      0.045321059435528D+00, &
      0.072757916845420D+00, &
      0.028327242531057D+00, &
      0.009421666963733D+00 /)

  else if ( rule == 11 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
     -0.069222096541517D+00, 0.534611048270758D+00, 0.534611048270758D+00, &
      0.202061394068290D+00, 0.398969302965855D+00, 0.398969302965855D+00, &
      0.593380199137435D+00, 0.203309900431282D+00, 0.203309900431282D+00, &
      0.761298175434837D+00, 0.119350912282581D+00, 0.119350912282581D+00, &
      0.935270103777448D+00, 0.032364948111276D+00, 0.032364948111276D+00, &
      0.050178138310495D+00, 0.356620648261293D+00, 0.593201213428213D+00, &
      0.021022016536166D+00, 0.171488980304042D+00, 0.807489003159792D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.000927006328961D+00, &
      0.077149534914813D+00, &
      0.059322977380774D+00, &
      0.036184540503418D+00, &
      0.013659731002678D+00, &
      0.052337111962204D+00, &
      0.020707659639141D+00 /)

  else if ( rule == 12 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.023565220452390D+00, 0.488217389773805D+00, 0.488217389773805D+00, &
      0.120551215411079D+00, 0.439724392294460D+00, 0.439724392294460D+00, &
      0.457579229975768D+00, 0.271210385012116D+00, 0.271210385012116D+00, &
      0.744847708916828D+00, 0.127576145541586D+00, 0.127576145541586D+00, &
      0.957365299093579D+00, 0.021317350453210D+00, 0.021317350453210D+00, &
      0.115343494534698D+00, 0.275713269685514D+00, 0.608943235779788D+00, &
      0.022838332222257D+00, 0.281325580989940D+00, 0.695836086787803D+00, &
      0.025734050548330D+00, 0.116251915907597D+00, 0.858014033544073D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.025731066440455D+00, &
      0.043692544538038D+00, &
      0.062858224217885D+00, &
      0.034796112930709D+00, &
      0.006166261051559D+00, &
      0.040371557766381D+00, &
      0.022356773202303D+00, &
      0.017316231108659D+00 /)

  else if ( rule == 13 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.009903630120591D+00, 0.495048184939705D+00, 0.495048184939705D+00, &
      0.062566729780852D+00, 0.468716635109574D+00, 0.468716635109574D+00, &
      0.170957326397447D+00, 0.414521336801277D+00, 0.414521336801277D+00, &
      0.541200855914337D+00, 0.229399572042831D+00, 0.229399572042831D+00, &
      0.771151009607340D+00, 0.114424495196330D+00, 0.114424495196330D+00, &
      0.950377217273082D+00, 0.024811391363459D+00, 0.024811391363459D+00, &
      0.094853828379579D+00, 0.268794997058761D+00, 0.636351174561660D+00, &
      0.018100773278807D+00, 0.291730066734288D+00, 0.690169159986905D+00, &
      0.022233076674090D+00, 0.126357385491669D+00, 0.851409537834241D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.052520923400802D+00, &
      0.011280145209330D+00, &
      0.031423518362454D+00, &
      0.047072502504194D+00, &
      0.047363586536355D+00, &
      0.031167529045794D+00, &
      0.007975771465074D+00, &
      0.036848402728732D+00, &
      0.017401463303822D+00, &
      0.015521786839045D+00 /)

  else if ( rule == 14 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.022072179275643D+00, 0.488963910362179D+00, 0.488963910362179D+00, &
      0.164710561319092D+00, 0.417644719340454D+00, 0.417644719340454D+00, &
      0.453044943382323D+00, 0.273477528308839D+00, 0.273477528308839D+00, &
      0.645588935174913D+00, 0.177205532412543D+00, 0.177205532412543D+00, &
      0.876400233818255D+00, 0.061799883090873D+00, 0.061799883090873D+00, &
      0.961218077502598D+00, 0.019390961248701D+00, 0.019390961248701D+00, &
      0.057124757403648D+00, 0.172266687821356D+00, 0.770608554774996D+00, &
      0.092916249356972D+00, 0.336861459796345D+00, 0.570222290846683D+00, &
      0.014646950055654D+00, 0.298372882136258D+00, 0.686980167808088D+00, &
      0.001268330932872D+00, 0.118974497696957D+00, 0.879757171370171D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.021883581369429D+00, &
      0.032788353544125D+00, &
      0.051774104507292D+00, &
      0.042162588736993D+00, &
      0.014433699669777D+00, &
      0.004923403602400D+00, &
      0.024665753212564D+00, &
      0.038571510787061D+00, &
      0.014436308113534D+00, &
      0.005010228838501D+00 /)

  else if ( rule == 15 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
     -0.013945833716486D+00, 0.506972916858243D+00, 0.506972916858243D+00, &
      0.137187291433955D+00, 0.431406354283023D+00, 0.431406354283023D+00, &
      0.444612710305711D+00, 0.277693644847144D+00, 0.277693644847144D+00, &
      0.747070217917492D+00, 0.126464891041254D+00, 0.126464891041254D+00, &
      0.858383228050628D+00, 0.070808385974686D+00, 0.070808385974686D+00, &
      0.962069659517853D+00, 0.018965170241073D+00, 0.018965170241073D+00, &
      0.133734161966621D+00, 0.261311371140087D+00, 0.604954466893291D+00, &
      0.036366677396917D+00, 0.388046767090269D+00, 0.575586555512814D+00, &
     -0.010174883126571D+00, 0.285712220049916D+00, 0.724462663076655D+00, &
      0.036843869875878D+00, 0.215599664072284D+00, 0.747556466051838D+00, &
      0.012459809331199D+00, 0.103575616576386D+00, 0.883964574092416D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.001916875642849D+00, &
      0.044249027271145D+00, &
      0.051186548718852D+00, &
      0.023687735870688D+00, &
      0.013289775690021D+00, &
      0.004748916608192D+00, &
      0.038550072599593D+00, &
      0.027215814320624D+00, &
      0.002182077366797D+00, &
      0.021505319847731D+00, &
      0.007673942631049D+00 /)

  else if ( rule == 16 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.005238916103123D+00, 0.497380541948438D+00, 0.497380541948438D+00, &
      0.173061122901295D+00, 0.413469438549352D+00, 0.413469438549352D+00, &
      0.059082801866017D+00, 0.470458599066991D+00, 0.470458599066991D+00, &
      0.518892500060958D+00, 0.240553749969521D+00, 0.240553749969521D+00, &
      0.704068411554854D+00, 0.147965794222573D+00, 0.147965794222573D+00, &
      0.849069624685052D+00, 0.075465187657474D+00, 0.075465187657474D+00, &
      0.966807194753950D+00, 0.016596402623025D+00, 0.016596402623025D+00, &
      0.103575692245252D+00, 0.296555596579887D+00, 0.599868711174861D+00, &
      0.020083411655416D+00, 0.337723063403079D+00, 0.642193524941505D+00, &
     -0.004341002614139D+00, 0.204748281642812D+00, 0.799592720971327D+00, &
      0.041941786468010D+00, 0.189358492130623D+00, 0.768699721401368D+00, &
      0.014317320230681D+00, 0.085283615682657D+00, 0.900399064086661D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.046875697427642D+00, &
      0.006405878578585D+00, &
      0.041710296739387D+00, &
      0.026891484250064D+00, &
      0.042132522761650D+00, &
      0.030000266842773D+00, &
      0.014200098925024D+00, &
      0.003582462351273D+00, &
      0.032773147460627D+00, &
      0.015298306248441D+00, &
      0.002386244192839D+00, &
      0.019084792755899D+00, &
      0.006850054546542D+00 /)

  else if ( rule == 17 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.005658918886452D+00, 0.497170540556774D+00, 0.497170540556774D+00, &
      0.035647354750751D+00, 0.482176322624625D+00, 0.482176322624625D+00, &
      0.099520061958437D+00, 0.450239969020782D+00, 0.450239969020782D+00, &
      0.199467521245206D+00, 0.400266239377397D+00, 0.400266239377397D+00, &
      0.495717464058095D+00, 0.252141267970953D+00, 0.252141267970953D+00, &
      0.675905990683077D+00, 0.162047004658461D+00, 0.162047004658461D+00, &
      0.848248235478508D+00, 0.075875882260746D+00, 0.075875882260746D+00, &
      0.968690546064356D+00, 0.015654726967822D+00, 0.015654726967822D+00, &
      0.010186928826919D+00, 0.334319867363658D+00, 0.655493203809423D+00, &
      0.135440871671036D+00, 0.292221537796944D+00, 0.572337590532020D+00, &
      0.054423924290583D+00, 0.319574885423190D+00, 0.626001190286228D+00, &
      0.012868560833637D+00, 0.190704224192292D+00, 0.796427214974071D+00, &
      0.067165782413524D+00, 0.180483211648746D+00, 0.752351005937729D+00, &
      0.014663182224828D+00, 0.080711313679564D+00, 0.904625504095608D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.033437199290803D+00, &
      0.005093415440507D+00, &
      0.014670864527638D+00, &
      0.024350878353672D+00, &
      0.031107550868969D+00, &
      0.031257111218620D+00, &
      0.024815654339665D+00, &
      0.014056073070557D+00, &
      0.003194676173779D+00, &
      0.008119655318993D+00, &
      0.026805742283163D+00, &
      0.018459993210822D+00, &
      0.008476868534328D+00, &
      0.018292796770025D+00, &
      0.006665632004165D+00 /)

  else if ( rule == 18 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.013310382738157D+00, 0.493344808630921D+00, 0.493344808630921D+00, &
      0.061578811516086D+00, 0.469210594241957D+00, 0.469210594241957D+00, &
      0.127437208225989D+00, 0.436281395887006D+00, 0.436281395887006D+00, &
      0.210307658653168D+00, 0.394846170673416D+00, 0.394846170673416D+00, &
      0.500410862393686D+00, 0.249794568803157D+00, 0.249794568803157D+00, &
      0.677135612512315D+00, 0.161432193743843D+00, 0.161432193743843D+00, &
      0.846803545029257D+00, 0.076598227485371D+00, 0.076598227485371D+00, &
      0.951495121293100D+00, 0.024252439353450D+00, 0.024252439353450D+00, &
      0.913707265566071D+00, 0.043146367216965D+00, 0.043146367216965D+00, &
      0.008430536202420D+00, 0.358911494940944D+00, 0.632657968856636D+00, &
      0.131186551737188D+00, 0.294402476751957D+00, 0.574410971510855D+00, &
      0.050203151565675D+00, 0.325017801641814D+00, 0.624779046792512D+00, &
      0.066329263810916D+00, 0.184737559666046D+00, 0.748933176523037D+00, &
      0.011996194566236D+00, 0.218796800013321D+00, 0.769207005420443D+00, &
      0.014858100590125D+00, 0.101179597136408D+00, 0.883962302273467D+00, &
     -0.035222015287949D+00, 0.020874755282586D+00, 1.014347260005363D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.030809939937647D+00, &
      0.009072436679404D+00, &
      0.018761316939594D+00, &
      0.019441097985477D+00, &
      0.027753948610810D+00, &
      0.032256225351457D+00, &
      0.025074032616922D+00, &
      0.015271927971832D+00, &
      0.006793922022963D+00, &
     -0.002223098729920D+00, &
      0.006331914076406D+00, &
      0.027257538049138D+00, &
      0.017676785649465D+00, &
      0.018379484638070D+00, &
      0.008104732808192D+00, &
      0.007634129070725D+00, &
      0.000046187660794D+00  &
      /)

  else if ( rule == 19 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
      0.020780025853987D+00, 0.489609987073006D+00, 0.489609987073006D+00, &
      0.090926214604215D+00, 0.454536892697893D+00, 0.454536892697893D+00, &
      0.197166638701138D+00, 0.401416680649431D+00, 0.401416680649431D+00, &
      0.488896691193805D+00, 0.255551654403098D+00, 0.255551654403098D+00, &
      0.645844115695741D+00, 0.177077942152130D+00, 0.177077942152130D+00, &
      0.779877893544096D+00, 0.110061053227952D+00, 0.110061053227952D+00, &
      0.888942751496321D+00, 0.055528624251840D+00, 0.055528624251840D+00, &
      0.974756272445543D+00, 0.012621863777229D+00, 0.012621863777229D+00, &
      0.003611417848412D+00, 0.395754787356943D+00, 0.600633794794645D+00, &
      0.134466754530780D+00, 0.307929983880436D+00, 0.557603261588784D+00, &
      0.014446025776115D+00, 0.264566948406520D+00, 0.720987025817365D+00, &
      0.046933578838178D+00, 0.358539352205951D+00, 0.594527068955871D+00, &
      0.002861120350567D+00, 0.157807405968595D+00, 0.839331473680839D+00, &
      0.223861424097916D+00, 0.075050596975911D+00, 0.701087978926173D+00, &
      0.034647074816760D+00, 0.142421601113383D+00, 0.822931324069857D+00, &
      0.010161119296278D+00, 0.065494628082938D+00, 0.924344252620784D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.032906331388919D+00, &
      0.010330731891272D+00, &
      0.022387247263016D+00, &
      0.030266125869468D+00, &
      0.030490967802198D+00, &
      0.024159212741641D+00, &
      0.016050803586801D+00, &
      0.008084580261784D+00, &
      0.002079362027485D+00, &
      0.003884876904981D+00, &
      0.025574160612022D+00, &
      0.008880903573338D+00, &
      0.016124546761731D+00, &
      0.002491941817491D+00, &
      0.018242840118951D+00, &
      0.010258563736199D+00, &
      0.003799928855302D+00  &
      /)

  else if ( rule == 20 ) then

    suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
      0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00, &
     -0.001900928704400D+00, 0.500950464352200D+00, 0.500950464352200D+00, &
      0.023574084130543D+00, 0.488212957934729D+00, 0.488212957934729D+00, &
      0.089726636099435D+00, 0.455136681950283D+00, 0.455136681950283D+00, &
      0.196007481363421D+00, 0.401996259318289D+00, 0.401996259318289D+00, &
      0.488214180481157D+00, 0.255892909759421D+00, 0.255892909759421D+00, &
      0.647023488009788D+00, 0.176488255995106D+00, 0.176488255995106D+00, &
      0.791658289326483D+00, 0.104170855336758D+00, 0.104170855336758D+00, &
      0.893862072318140D+00, 0.053068963840930D+00, 0.053068963840930D+00, &
      0.916762569607942D+00, 0.041618715196029D+00, 0.041618715196029D+00, &
      0.976836157186356D+00, 0.011581921406822D+00, 0.011581921406822D+00, &
      0.048741583664839D+00, 0.344855770229001D+00, 0.606402646106160D+00, &
      0.006314115948605D+00, 0.377843269594854D+00, 0.615842614456541D+00, &
      0.134316520547348D+00, 0.306635479062357D+00, 0.559048000390295D+00, &
      0.013973893962392D+00, 0.249419362774742D+00, 0.736606743262866D+00, &
      0.075549132909764D+00, 0.212775724802802D+00, 0.711675142287434D+00, &
     -0.008368153208227D+00, 0.146965436053239D+00, 0.861402717154987D+00, &
      0.026686063258714D+00, 0.137726978828923D+00, 0.835586957912363D+00, &
      0.010547719294141D+00, 0.059696109149007D+00, 0.929756171556853D+00  &
      /), (/ 3, suborder_num /) )

    suborder_w(1:suborder_num) = (/ &
      0.033057055541624D+00, &
      0.000867019185663D+00, &
      0.011660052716448D+00, &
      0.022876936356421D+00, &
      0.030448982673938D+00, &
      0.030624891725355D+00, &
      0.024368057676800D+00, &
      0.015997432032024D+00, &
      0.007698301815602D+00, &
     -0.000632060497488D+00, &
      0.001751134301193D+00, &
      0.016465839189576D+00, &
      0.004839033540485D+00, &
      0.025804906534650D+00, &
      0.008471091054441D+00, &
      0.018354914106280D+00, &
      0.000704404677908D+00, &
      0.010112684927462D+00, &
      0.003573909385950D+00  &
      /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DUNAVANT_SUBRULE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop

  end if

  return
end
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

  character              c
  integer   ( kind = 4 ) change
  integer   ( kind = 4 ) digit
  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lens

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
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the value.
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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

  real ( kind = 8 ) area
  real ( kind = 8 ) node_xy(2,3)

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
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_show
  real ( kind = 8 ) point_xy(2,point_num)
  character ( len = 40 ) string
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108
  real ( kind = 8 ) y_scale
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
