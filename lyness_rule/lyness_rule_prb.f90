program main

!*****************************************************************************80
!
!! MAIN is the main program for LYNESS_RULE_PRB.
!
!  Discussion:
!
!    LYNESS_RULE_PRB tests the LYNESS_RULE library.
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
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LYNESS_RULE_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LYNESS_RULE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LYNESS_RULE_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests LYNESS_RULE_NUM, LYNESS_DEGREE, LYNESS_ORDER_NUM.
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
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) precision
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  LYNESS_RULE_NUM returns the number of rules;'
  write ( *, '(a)' ) '  LYNESS_DEGREE returns the degree of a rule;'
  write ( *, '(a)' ) '  LYNESS_ORDER_NUM returns the order of a rule.'

  call lyness_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of available rules = ', rule_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule     Order  Precision'
  write ( *, '(a)' ) ' '

  do rule = 0, rule_num
    call lyness_order ( rule, order )
    call lyness_precision ( rule, precision )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) rule, order, precision
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 performs the weight sum test on Lyness rules.
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
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ) w_sum
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  LYNESS_RULE returns the points and weights'
  write ( *, '(a)' ) '  of a Lyness rule for the triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply check that the weights'
  write ( *, '(a)' ) '  sum to 1.'

  call lyness_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of available rules = ', rule_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule    Sum of weights'
  write ( *, '(a)' ) ' '

  do rule = 0, rule_num

    call lyness_order ( rule, order )

    allocate ( x(1:2,1:order) )
    allocate ( w(1:order) )

    call lyness_rule ( rule, order, w, x )

    w_sum = sum ( w(1:order) )

    write ( *, '(2x,i8,2x,g25.16)' ) rule, w_sum

    deallocate ( w )
    deallocate ( x )
    
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 performs the barycentric coordinate sum test on Lyness rules.
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
  implicit none

  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  integer ( kind = 4 ) suborder
  integer ( kind = 4 ) suborder_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: sub_w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sub_xyz
  real ( kind = 8 ) xyz_sum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  LYNESS_RULE returns the points and weights'
  write ( *, '(a)' ) '  of a Lyness rule for the triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply check that, for each'
  write ( *, '(a)' ) '  quadrature point, the barycentric coordinates'
  write ( *, '(a)' ) '  sum to 1.'

  call lyness_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule   Suborder    Sum of coordinates'

  do rule = 0, rule_num

    call lyness_suborder_num ( rule, suborder_num )

    allocate ( sub_w(1:suborder_num) )
    allocate ( sub_xyz(1:3,1:suborder_num) )

    call lyness_subrule ( rule, suborder_num, sub_xyz, sub_w )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8,2x,i8)' ) rule, suborder_num
    do suborder = 1, suborder_num
      xyz_sum = sum ( sub_xyz(1:3,suborder) )
      write ( *, '(20x,2x,g25.16)' ) xyz_sum
    end do

    deallocate ( sub_w )
    deallocate ( sub_xyz )
    
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 prints a rule generated by LYNESS_RULE.
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
  implicit none

  integer ( kind = 4 ) j
  integer ( kind = 4 ) order
  integer ( kind = 4 ) precision
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  LYNESS_RULE returns the points and weights'
  write ( *, '(a)' ) '  of a Lyness rule for the triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply print a rule.'

  rule = 18
  call lyness_order ( rule, order )
  call lyness_precision ( rule, precision )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Rule =      ', rule
  write ( *, '(a,i8)' ) '  Order =     ', order
  write ( *, '(a,i8)' ) '  Precision = ', precision

  allocate ( w(1:order) )
  allocate ( x(1:2,1:order) )

  call lyness_rule ( rule, order, w, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      W               X               Y'
  write ( *, '(a)' ) ' '

  do j = 1, order
    write ( *, '(2x,i8,2x,g24.16,2x,g24.16,2x,g24.16)' ) &
      j, w(j), x(1:2,j)
  end do

  deallocate ( w )
  deallocate ( x )
    
  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 writes a rule created by LYNESS_RULE to a file.
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
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) order
  integer ( kind = 4 ) precision
  real ( kind = 8 ) :: r(2,3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  LYNESS_RULE returns the points and weights'
  write ( *, '(a)' ) '  of a Lyness rule for the triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply print a rule.'

  rule = 18
  call lyness_order ( rule, order )
  call lyness_precision ( rule, precision )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Rule =      ', rule
  write ( *, '(a,i8)' ) '  Order =     ', order
  write ( *, '(a,i8)' ) '  Precision = ', precision

  allocate ( w(1:order) )
  allocate ( x(1:2,1:order) )

  call lyness_rule ( rule, order, w, x )

  call r8mat_write ( 'lyness_18_r.txt', 2, 3, r )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the region file "lyness_18_r.txt".'
  call r8mat_write ( 'lyness_18_w.txt', 1, order, w )
  write ( *, '(a)' ) '  Wrote the weight file "lyness_18_w.txt".'
  call r8mat_write ( 'lyness_18_x.txt', 2, order, x )
  write ( *, '(a)' ) '  Wrote the point file "lyness_18_x.txt".'

  deallocate ( w )
  deallocate ( x )
    
  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests the Lyness rules for exact integration of monomials.
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
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) area
  integer ( kind = 4 ) b
  real ( kind = 8 ) coef
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) :: degree_max = 10
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) order
  real ( kind = 8 ) quad
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  real ( kind = 8 ) value
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  LYNESS_RULE returns the points and weights of'
  write ( *, '(a)' ) '  a Lyness rule for the unit triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This routine uses those rules to estimate the'
  write ( *, '(a)' ) '  integral of monomomials in the unit triangle.'

  call lyness_rule_num ( rule_num )

  area = 0.5D+00

  do degree = 0, degree_max

    do a = 0, degree

      b = degree - a
!
!  Multiplying X^A * Y^B by COEF will give us an integrand
!  whose integral is exactly 1.  This makes the error calculations easy.
!
      coef = real ( a + b + 2, kind = 8 ) * real ( a + b + 1, kind = 8 )
      do i = 1, b
        coef = coef * real ( a + i, kind = 8 ) / real ( i, kind = 8 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i8,a,i8)' ) &
        '  Integrate ', coef , ' * X ^ ', a, ' * Y ^ ', b
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      Rule       QUAD           ERROR'
      write ( *, '(a)' ) ' '

      do rule = 0, rule_num

        call lyness_order ( rule, order )

        allocate ( x(1:2,1:order) )
        allocate ( w(1:order) )

        call lyness_rule ( rule, order, w, x )

        quad = 0.0D+00

        do j = 1, order

          if ( a == 0 .and. b == 0 ) then
            value = coef
          else if ( a == 0 .and. b /= 0 ) then
            value = coef        * x(2,j)**b
          else if ( a /= 0 .and. b == 0 ) then
            value = coef * x(1,j)**a
          else if ( a /= 0 .and. b /= 0 ) then
            value = coef * x(1,j)**a * x(2,j)**b
          end if

          quad = quad + w(j) * value

        end do

        quad = area * quad

        exact = 1.0D+00
        err = abs ( exact - quad )

        write ( *, '(2x,i8,2x,g14.6,2x,g10.2)' ) rule, quad, err

        deallocate ( w )
        deallocate ( x )
     
      end do

    end do

  end do

  return
end
