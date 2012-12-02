program main

!*****************************************************************************80
!
!! MAIN is the main program for KEAST_PRB.
!
!  Discussion:
!
!    KEAST_PRB runs tests on the KEAST library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KEAST_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the KEAST library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KEAST_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests KEAST_RULE_NUM, KEAST_DEGREE, KEAST_ORDER_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  KEAST_RULE_NUM returns the number of rules;'
  write ( *, '(a)' ) '  KEAST_DEGREE returns the degree of a rule;'
  write ( *, '(a)' ) '  KEAST_ORDER_NUM returns the order of a rule.'

  call keast_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of available rules = ', rule_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule      Degree     Order'
  write ( *, '(a)' ) ' '

  do rule = 1, rule_num
    call keast_order_num ( rule, order_num )
    call keast_degree ( rule, degree )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) rule, degree, order_num
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests KEAST_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: wtab
  real ( kind = 8 ) wtab_sum
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyztab

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  KEAST_RULE returns the points and weights'
  write ( *, '(a)' ) '  of a Keast rule for the tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply check that the weights'
  write ( *, '(a)' ) '  sum to 1.'

  call keast_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of available rules = ', rule_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule     Sum of weights'
  write ( *, '(a)' ) ' '

  do rule = 1, rule_num

    call keast_order_num ( rule, order_num )

    allocate ( wtab(1:order_num) )
    allocate ( xyztab(1:3,1:order_num) )

    call keast_rule ( rule, order_num, xyztab, wtab )

    wtab_sum = sum ( wtab(1:order_num) )

    write ( *, '(2x,i8,2x,g25.16)' ) rule, wtab_sum

    deallocate ( wtab )
    deallocate ( xyztab )
    
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests KEAST_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2006
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
  real ( kind = 8 ), allocatable, dimension ( : ) :: suborder_w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyzz
  real ( kind = 8 ) xyzz_sum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  KEAST_RULE returns the points and weights'
  write ( *, '(a)' ) '  of a Keast rule for the tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply check that, for each'
  write ( *, '(a)' ) '  quadrature point, the barycentric coordinates'
  write ( *, '(a)' ) '  sum to 1.'

  call keast_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule   Suborder    Sum of coordinates'
  write ( *, '(a)' ) ' '

  do rule = 1, rule_num

    call keast_suborder_num ( rule, suborder_num )

    allocate ( suborder_w(1:suborder_num) )
    allocate ( suborder_xyzz(1:4,1:suborder_num) )

    call keast_subrule ( rule, suborder_num, suborder_xyzz, suborder_w )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8,2x,i8)' ) rule, suborder_num
    do suborder = 1, suborder_num
      xyzz_sum = sum ( suborder_xyzz(1:4,suborder) )
      write ( *, '(20x,2x,g25.16)' ) xyzz_sum
    end do

    deallocate ( suborder_w )
    deallocate ( suborder_xyzz )
    
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 demonstrates TETRAHEDRON_REFERENCE_TO_PHYSICAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) node
  real ( kind = 8 ), dimension ( 3, node_num ) :: node_xyz = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ 3, node_num /) )
  real ( kind = 8 ), dimension ( 3, node_num ) :: node_xyz2 = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    2.0D+00, 2.0D+00, 3.0D+00, &
    1.0D+00, 3.0D+00, 3.0D+00, &
    1.0D+00, 2.0D+00, 9.0D+00 /), (/ 3, node_num /) )
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume2
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  TETRAHEDRON_REFERENCE_TO_PHYSICAL transforms a rule'
  write ( *, '(a)' ) '  on the unit (reference) tetrahedron to a rule on '
  write ( *, '(a)' ) '  an arbitrary (physical) tetrahedron.'

  rule = 2

  call keast_order_num ( rule, order_num )

  allocate ( xyz(1:3,1:order_num) )
  allocate ( xyz2(1:3,1:order_num) )
  allocate ( w(1:order_num) )

  call keast_rule ( rule, order_num, xyz, w )
!
!  Here is the reference tetrahedron, and its rule.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The reference tetrahedron:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) node, node_xyz(1:3,node)
  end do

  call tetrahedron_volume ( node_xyz, volume )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Rule ', rule, ' for reference tetrahedron'
  write ( *, '(a,g14.6)' ) '  with volume = ', volume
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '                X               Y               Z              W'
  write ( *, '(a)' ) ' '

  do order = 1, order_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      order, xyz(1:3,order), w(order)
  end do
!
!  Transform the rule.
!
  call tetrahedron_reference_to_physical ( node_xyz2, order_num, xyz, xyz2 )
!
!  Here is the physical tetrahedron, and its transformed rule.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The physical tetrahedron:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) node, node_xyz2(1:3,node)
  end do

  call tetrahedron_volume ( node_xyz2, volume2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Rule ', rule, ' for physical tetrahedron'
  write ( *, '(a,g14.6)' ) '  with volume = ', volume2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '                X               Y               Z              W'
  write ( *, '(a)' ) ' '

  do order = 1, order_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      order, xyz2(1:3,order), w(order)
  end do

  deallocate ( w )
  deallocate ( xyz )
  deallocate ( xyz2 )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 demonstrates Keast rules on various monomials in the unit tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) degree
  integer ( kind = 4 ), parameter:: degree_max = 3
  integer ( kind = 4 ) expon(dim_num)
  integer ( kind = 4 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: mono
  logical              more
  integer ( kind = 4 ) order_num
  real ( kind = 8 ) quad
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  integer ( kind = 4 ) t
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Demonstrate the KEAST rules on a sequence of '
  write ( *, '(a)' ) '  monomial integrands X^A Y^B Z^C '
  write ( *, '(a)' ) '  on the unit tetrahedron.'

  call keast_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule     Order     Quad'
  write ( *, '(a)' ) ' '

  do degree = 0, degree_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( degree, dim_num, expon, more, h, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i1,a,i1,a,i1)' ) &
        '  F(X,Y,Z) = X^', expon(1), ' * Y^', expon(2), ' * Z^', expon(3)
      write ( *, '(a)' ) ' '

      do rule = 1, rule_num

        call keast_order_num ( rule, order_num )

        allocate ( mono(1:order_num) )
        allocate ( w(1:order_num) )
        allocate ( xyz(1:3,1:order_num) )

        call keast_rule ( rule, order_num, xyz, w )

        call monomial_value ( dim_num, order_num, xyz, expon, mono )

        quad = dot_product ( w(1:order_num), mono(1:order_num) )

        write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) rule, order_num, quad

        deallocate ( mono )
        deallocate ( w )
        deallocate ( xyz )
    
      end do

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests KEAST_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  KEAST_RULE returns the points and weights'
  write ( *, '(a)' ) '  of a Keast rule for the tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply print a rule.'

  rule = 10
  call keast_degree ( rule, degree )
  call keast_order_num ( rule, order_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Rule =   ', rule
  write ( *, '(a,i8)' ) '  Degree = ', degree
  write ( *, '(a,i8)' ) '  Order =  ', order

  allocate ( w(1:order_num) )
  allocate ( xyz(1:3,1:order_num) )

  call keast_rule ( rule, order_num, xyz, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         I      W               X               Y               Z'
  write ( *, '(a)' ) ' '

  do order = 1, order_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      order, w(order), xyz(1:3,order)
  end do

  deallocate ( w )
  deallocate ( xyz )
    
  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests KEAST_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  character ( len = 12 ) w_file
  integer ( kind = 4 ) w_unit
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz
  character ( len = 12 ) xyz_file
  integer ( kind = 4 ) xyz_unit

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  KEAST_RULE returns the points and weights'
  write ( *, '(a)' ) '  of a Keast rule for the tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we write each rule to a file.'
  write ( *, '(a)' ) ' '

  call keast_rule_num ( rule_num )

  do rule = 1, rule_num

    call keast_degree ( rule, degree )
    call keast_order_num ( rule, order_num )
 
    allocate ( w(1:order_num) )
    allocate ( xyz(1:3,1:order_num) )

    call keast_rule ( rule, order_num, xyz, w )

    call get_unit ( w_unit )

    write ( w_file, '(a5,i1,a6)' ) 'keast', rule - 1, '_w.txt'

    open ( unit = w_unit, file = w_file, status = 'replace' )

    do order = 1, order_num
      write ( w_unit, '(f20.16)' ) w(order)
    end do

    close ( unit = w_unit )

    call get_unit ( xyz_unit )

    write ( xyz_file, '(a5,i1,a6)' ) 'keast', rule - 1, '_x.txt'

    open ( unit = xyz_unit, file = xyz_file, status = 'replace' )

    do order = 1, order_num
      write ( xyz_unit, '(3f20.16)' ) xyz(1:3,order)
    end do

    close ( unit = xyz_unit )

    write ( *, '(a,i2,a)' ) '  Wrote rule ', rule, ' to "' &
      // trim ( w_file ) // '" and "' // trim ( xyz_file ) // '".'

    deallocate ( w )
    deallocate ( xyz )

  end do
    
  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 demonstrates a Keast rule on monomials in a general tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: node_num = 4

  integer ( kind = 4 ) degree
  integer ( kind = 4 ), parameter:: degree_max = 4
  integer ( kind = 4 ) expon(dim_num)
  integer ( kind = 4 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: mono
  logical more
  integer ( kind = 4 ) node
  real ( kind = 8 ), dimension ( 3, node_num ) :: node_xyz2 = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    4.0D+00, 1.0D+00, 2.0D+00, &
    2.0D+00, 4.0D+00, 4.0D+00, &
    3.0D+00, 2.0D+00, 5.0D+00 /), (/ 3, node_num /) )
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_num
  real ( kind = 8 ) quad
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  integer ( kind = 4 ) t
  real ( kind = 8 ) volume
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Transform a Keast rule to a general tetrahedron.'
  write ( *, '(a)' ) '  Apply it to monomial integrands.'
!
!  Get Keast rule #7.
!
  rule = 7

  call keast_order_num ( rule, order_num )

  allocate ( xyz(1:3,1:order_num) )
  allocate ( xyz2(1:3,1:order_num) )
  allocate ( w(1:order_num) )

  call keast_rule ( rule, order_num, xyz, w )
!
!  Transform the rule.
!
  call tetrahedron_reference_to_physical ( node_xyz2, order_num, xyz, xyz2 )
!
!  Here is the physical tetrahedron, and its transformed rule.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The physical tetrahedron:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) node, node_xyz2(1:3,node)
  end do

  call tetrahedron_volume ( node_xyz2, volume )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Rule ', rule, ' for physical tetrahedron'
  write ( *, '(a,g14.6)' ) '  with volume = ', volume
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '                X               Y               Z              W'
  write ( *, '(a)' ) ' '

  do order = 1, order_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      order, xyz2(1:3,order), w(order)
  end do

  allocate ( mono(1:order_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule     Order     Quad'
  write ( *, '(a)' ) ' '

  do degree = 0, degree_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( degree, dim_num, expon, more, h, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i1,a,i1,a,i1)' ) &
        '  F(X,Y,Z) = X^', expon(1), ' * Y^', expon(2), ' * Z^', expon(3)
      write ( *, '(a)' ) ' '

      rule = 7

      call monomial_value ( dim_num, order_num, xyz2, expon, mono )

      quad = dot_product ( w(1:order_num), mono(1:order_num) )

      write ( *, '(2x,i8,2x,i8,2x,g24.16)' ) rule, order_num, quad * volume
    
      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( mono )
  deallocate ( w )
  deallocate ( xyz )
  deallocate ( xyz2 )

  return
end
