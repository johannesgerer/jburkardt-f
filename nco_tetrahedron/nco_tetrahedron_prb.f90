program main

!*****************************************************************************80
!
!! MAIN is the main program for NCO_TETRAHEDRON_PRB.
!
!  Discussion:
!
!    NCO_TETRAHEDRON_PRB runs tests on the NCO_TETRAHEDRON library.
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

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NCO_TETRAHEDRON_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NCO_TETRAHEDRON library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NCO_TETRAHEDRON_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests NCO_TETRAHEDRON_RULE_NUM, *_DEGREE, *_ORDER_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2007
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
  write ( *, '(a)' ) '  NCO_TETRAHEDRON_RULE_NUM returns the number of rules;'
  write ( *, '(a)' ) '  NCO_TETRAHEDRON_DEGREE returns the degree of a rule;'
  write ( *, '(a)' ) '  NCO_TETRAHEDRON_ORDER_NUM returns the order of a rule.'

  call nco_tetrahedron_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of available rules = ', rule_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule      Degree     Order'
  write ( *, '(a)' ) ' '

  do rule = 1, rule_num
    call nco_tetrahedron_order_num ( rule, order_num )
    call nco_tetrahedron_degree ( rule, degree )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) rule, degree, order_num
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests NCO_TETRAHEDRON_RULE.
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
  implicit none

  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: wtab
  real ( kind = 8 ) wtab_sum
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyztab

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  NCO_TETRAHEDRON_RULE returns the points and weights'
  write ( *, '(a)' ) '  of an NCO rule for the tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply check that the weights'
  write ( *, '(a)' ) '  sum to 1.'

  call nco_tetrahedron_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of available rules = ', rule_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule     Sum of weights'
  write ( *, '(a)' ) ' '

  do rule = 1, rule_num

    call nco_tetrahedron_order_num ( rule, order_num )

    allocate ( wtab(1:order_num) )
    allocate ( xyztab(1:3,1:order_num) )

    call nco_tetrahedron_rule ( rule, order_num, xyztab, wtab )

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
!! TEST03 tests NCO_TETRAHEDRON_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2007
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
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyz
  real ( kind = 8 ) xyz_sum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  NCO_TETRAHEDRON_RULE returns the points and weights'
  write ( *, '(a)' ) '  of an NCO rule for the tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we simply check that, for each'
  write ( *, '(a)' ) '  quadrature point, the barycentric coordinates'
  write ( *, '(a)' ) '  sum to 1.'

  call nco_tetrahedron_rule_num ( rule_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule   Suborder    Sum of coordinates'
  write ( *, '(a)' ) ' '

  do rule = 1, rule_num

    call nco_tetrahedron_suborder_num ( rule, suborder_num )

    allocate ( suborder_w(1:suborder_num) )
    allocate ( suborder_xyz(1:4,1:suborder_num) )

    call nco_tetrahedron_subrule ( rule, suborder_num, suborder_xyz, suborder_w )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8,2x,i8)' ) rule, suborder_num
    do suborder = 1, suborder_num
      xyz_sum = sum ( suborder_xyz(1:4,suborder) )
      write ( *, '(20x,2x,g25.16)' ) xyz_sum
    end do

    deallocate ( suborder_w )
    deallocate ( suborder_xyz )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests NCO_TETRAHEDRON_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) coef
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order_num
  real ( kind = 8 ) quad
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  real ( kind = 8 ) value
  real ( kind = 8 ) volume
  real ( kind = 8 ), allocatable, dimension ( : ) :: wtab
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyztab
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  NCO_TETRAHEDRON_RULE returns the points and weights of'
  write ( *, '(a)' ) '  an NCO rule for the unit tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This routine uses those rules to estimate the'
  write ( *, '(a)' ) '  integral of monomomials in the unit tetrahedron.'

  call nco_tetrahedron_rule_num ( rule_num )

  volume = 1.0D+00 / 6.0D+00

  do a = 0, 6

    do b = 0, 6 - a

      do c = 0, 6 - a - b
!
!  Multiplying X**A * Y**B * Z**C by COEF will give us an integrand
!  whose integral is exactly 1.  This makes the error calculations easy.
!
        coef = 1.0D+00

!       do i = 1, a
!         coef = coef * i / i
!       end do
        do i = a + 1, a + b
          coef = coef * real ( i, kind = 8 ) / real ( i - a, kind = 8 )
        end do
        do i = a + b + 1, a + b + c
          coef = coef * real ( i, kind = 8 ) / real ( i - a - b, kind = 8 )
        end do
        do i = a + b + c + 1, a + b + c + 3
          coef = coef * real ( i, kind = 8 )
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6,a,i8,a,i8,a,i8)' ) &
        '  Integrate ', coef , ' * X ** ', a, ' * Y ** ', b, ' * Z ** ', c
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Rule       QUAD           ERROR'
        write ( *, '(a)' ) ' '

        do rule = 1, rule_num

          call nco_tetrahedron_order_num ( rule, order_num )

          allocate ( wtab(1:order_num) )
          allocate ( xyztab(1:3,1:order_num) )

          call nco_tetrahedron_rule ( rule, order_num, xyztab, wtab )

          quad = 0.0D+00

          do i = 1, order_num

            x = xyztab(1,i)
            y = xyztab(2,i)
            z = xyztab(3,i)
!
!  Some tedious calculations to avoid 0**0 complaints.
!
            value = coef

            if ( a /= 0 ) then
              value = value * x**a
            end if

            if ( b /= 0 ) then
              value = value * y**b
            end if

            if ( c /= 0 ) then
              value = value * z**c
            end if

            quad = quad + wtab(i) * value

          end do

          quad = quad * volume

          exact = 1.0D+00
          err = abs ( exact - quad )

          write ( *, '(2x,i8,2x,g14.6,2x,f14.8)' ) rule, quad, err

          deallocate ( wtab )
          deallocate ( xyztab )

        end do

      end do

    end do

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 demonstrates REFERENCE_TO_PHYSICAL_T4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: node_num = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) node
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xyz = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xyz2 = reshape ( (/ &
    4.0D+00, 5.0D+00, 1.0D+00, &
    6.0D+00, 5.0D+00, 1.0D+00, &
    4.0D+00, 8.0D+00, 1.0D+00, &
    4.0D+00, 5.0D+00, 5.0D+00 /), (/ dim_num, node_num /) )
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume2
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  REFERENCE_TO_PHYSICAL_T4 transforms a rule'
  write ( *, '(a)' ) '  on the unit (reference) tetrahedron to a rule on '
  write ( *, '(a)' ) '  an arbitrary (physical) tetrahedron.'

  rule = 3

  call nco_tetrahedron_order_num ( rule, order_num )

  allocate ( xyz(1:dim_num,1:order_num) )
  allocate ( xyz2(1:dim_num,1:order_num) )
  allocate ( w(1:order_num) )

  call nco_tetrahedron_rule ( rule, order_num, xyz, w )
!
!  Here is the reference tetrahedron, and its rule.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The reference tetrahedron:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      node, node_xyz(1:dim_num,node)
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
      order, xyz(1:dim_num,order), w(order)
  end do
!
!  Transform the rule.
!
  call reference_to_physical_t4 ( node_xyz2, order_num, xyz, xyz2 )
!
!  Here is the physical tetrahedron, and its transformed rule.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The physical tetrahedron:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      node, node_xyz2(1:dim_num,node)
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
      order, xyz2(1:dim_num,order), w(order)
  end do

  deallocate ( w )
  deallocate ( xyz )
  deallocate ( xyz2 )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests NCO_TETRAHEDRON_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) o
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) s
  integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
  integer ( kind = 4 ) suborder_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: suborder_w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyz
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  NCO_TETRAHEDRON_RULE returns the points and weights'
  write ( *, '(a)' ) '  of an NCO rule for the tetrahedron.'

  rule = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  In this test, we simply print rule ', rule

  call nco_tetrahedron_suborder_num ( rule, suborder_num )

  allocate ( suborder(1:suborder_num) )

  call nco_tetrahedron_suborder ( rule, suborder_num, suborder )

  allocate ( suborder_w(1:suborder_num) )
  allocate ( suborder_xyz(1:4,1:suborder_num) )

  call nco_tetrahedron_subrule ( rule, suborder_num, suborder_xyz, suborder_w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The compressed rule:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of suborders = ', suborder_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     S   Sub     Weight     Xsi1      Xsi2      Xsi3      Xsi4'
  write ( *, '(a)' ) ' '

  do s = 1, suborder_num
    write ( *, '(2x,i4,2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) &
      s, suborder(s), suborder_w(s), suborder_xyz(1:4,s)
  end do

  call nco_tetrahedron_order_num ( rule, order_num )

  allocate ( xyz(1:dim_num,1:order_num) )
  allocate ( w(1:order_num) )

  call nco_tetrahedron_rule ( rule, order_num, xyz, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The full rule:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Order = ', order_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     O    Weight      X         Y         Z'
  write ( *, '(a)' ) ' '

  do o = 1, order_num
    write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) &
      o, w(o), xyz(1:dim_num,o)
  end do

  deallocate ( suborder )
  deallocate ( suborder_w )
  deallocate ( suborder_xyz )
  deallocate ( w )
  deallocate ( xyz )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests NCO_TETRAHEDRON_RULE.
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

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) o
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) rule_num
  integer ( kind = 4 ) s
  integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
  integer ( kind = 4 ) suborder_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: suborder_w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyz
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  character ( len = 10 ) w_file
  integer ( kind = 4 ) w_unit
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xyz
  character ( len = 10 ) xyz_file
  integer ( kind = 4 ) xyz_unit

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  NCO_TETRAHEDRON_RULE returns the points and weights'
  write ( *, '(a)' ) '  of an NCO rule for the tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we write each rule to a file.'
  write ( *, '(a)' ) ' '

  call nco_tetrahedron_rule_num ( rule_num )

  do rule = 1, rule_num

    call nco_tetrahedron_suborder_num ( rule, suborder_num )

    allocate ( suborder(1:suborder_num) )

    call nco_tetrahedron_suborder ( rule, suborder_num, suborder )

    allocate ( suborder_w(1:suborder_num) )
    allocate ( suborder_xyz(1:4,1:suborder_num) )

    call nco_tetrahedron_subrule ( rule, suborder_num, suborder_xyz, &
      suborder_w )

    call nco_tetrahedron_order_num ( rule, order_num )

    allocate ( xyz(1:dim_num,1:order_num) )
    allocate ( w(1:order_num) )

    call nco_tetrahedron_rule ( rule, order_num, xyz, w )

    call get_unit ( w_unit )

    write ( w_file, '(a3,i1,a6)' ) 'nco', rule - 1, '_w.txt'

    open ( unit = w_unit, file = w_file, status = 'replace' )

    do o = 1, order_num
      write ( w_unit, '(f20.16)' ) w(o)
    end do

    close ( unit = w_unit )

    call get_unit ( xyz_unit )

    write ( xyz_file, '(a3,i1,a6)' ) 'nco', rule - 1, '_x.txt'

    open ( unit = xyz_unit, file = xyz_file, status = 'replace' )

    do o = 1, order_num
      write ( xyz_unit, '(3f20.16)' ) xyz(1:dim_num,o)
    end do

    close ( unit = xyz_unit )

    write ( *, '(a,i1,a)' ) '  Wrote rule ', rule, ' to "' &
      // trim ( w_file ) // '" and "' // trim ( xyz_file ) // '".'
    deallocate ( suborder )
    deallocate ( suborder_w )
    deallocate ( suborder_xyz )
    deallocate ( w )
    deallocate ( xyz )

  end do

  return
end
