program main

!*****************************************************************************80
!
!! MAIN is the main program for GM_RULE_PRB.
!
!  Discussion:
!
!    GM_RULE_PRB calls a set of problems for GM_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GM_RULE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GM_RULE library.'

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
  write ( *, '(a)' ) 'GM_RULE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SIMPLEX_UNIT_TO_GENERAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: vertex_num = dim_num + 1
  integer ( kind = 4 ), parameter :: point_num = 10

  integer ( kind = 4 ) j
  real ( kind = 8 ) phy(dim_num,point_num)
  real ( kind = 8 ) phy_unit(dim_num,dim_num+1)
  real ( kind = 8 ) ref(dim_num,point_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(dim_num,vertex_num) :: t = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, &
    2.0D+00, 5.0D+00 /), (/ dim_num, vertex_num /) )
  real ( kind = 8 ), dimension(dim_num,vertex_num) :: t_unit = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ dim_num, vertex_num /) )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_TO_GENERAL'
  write ( *, '(a)' ) '    maps points in the unit simplex to a general simplex.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we consider a simplex in 2D, a triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vertices of the general triangle are:'
  write ( *, '(a)' ) ' '
  do j = 1, vertex_num
    write ( *, '(2x,f8.4,2x,f8.4)' ) t(1:dim_num,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   (  XSI     ETA )   ( X       Y  )'
  write ( *, '(a)' ) ' '

  call simplex_unit_to_general ( dim_num, dim_num+1, t, t_unit, phy_unit )

  do j = 1, dim_num + 1

    write ( *, '(2x,2f8.4,2x,2f8.4)' ) t_unit(1:dim_num,j), phy_unit(1:dim_num,j)

  end do

  call simplex_unit_sample ( dim_num, point_num, seed, ref )

  call simplex_unit_to_general ( dim_num, point_num, t, ref, phy )

  do j = 1, point_num

    write ( *, '(2x,2f8.4,2x,2f8.4)' ) ref(1:dim_num,j), phy(1:dim_num,j)

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests SIMPLEX_UNIT_TO_GENERAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: vertex_num = dim_num + 1
  integer ( kind = 4 ), parameter :: point_num = 10

  integer ( kind = 4 ) j
  real ( kind = 8 ) phy(dim_num,point_num)
  real ( kind = 8 ) phy_unit(dim_num,dim_num+1)
  real ( kind = 8 ) ref(dim_num,point_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(dim_num,vertex_num) :: t = reshape ( (/ &
    1.0D+00, 1.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 4.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 5.0D+00 /), (/ dim_num, vertex_num /) )
  real ( kind = 8 ), dimension(dim_num,vertex_num) :: t_unit = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, vertex_num /) )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_TO_GENERAL'
  write ( *, '(a)' ) '    maps points in the unit simplex to a general simplex.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we consider a simplex in 3D, a tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vertices of the general tetrahedron are:'
  write ( *, '(a)' ) ' '
  do j = 1, vertex_num
    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4)' ) t(1:dim_num,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   (  XSI     ETA     MU )    ( X       Y       Z )'
  write ( *, '(a)' ) ' '

  call simplex_unit_to_general ( dim_num, dim_num+1, t, t_unit, phy_unit )

  do j = 1, dim_num + 1

    write ( *, '(2x,3f8.4,2x,3f8.4)' ) t_unit(1:dim_num,j), phy_unit(1:dim_num,j)

  end do

  call simplex_unit_sample ( dim_num, point_num, seed, ref )

  call simplex_unit_to_general ( dim_num, point_num, t, ref, phy )

  do j = 1, point_num

    write ( *, '(2x,3f8.4,2x,3f8.4)' ) ref(1:dim_num,j), phy(1:dim_num,j)

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests GM_RULE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), dimension ( test_num ) :: dim_num_test = (/ &
    2, 3, 5, 10 /)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  GM_RULE_SIZE returns POINT_NUM, the number of points'
  write ( *, '(a)' ) '  associated with a Grundmann-Moeller quadrature rule'
  write ( *, '(a)' ) '  for the unit simplex of dimension DIM_NUM'
  write ( *, '(a)' ) '  with rule index RULE'
  write ( *, '(a)' ) '  and degree of exactness DEGREE = 2*RULE+1.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   DIM_NUM      RULE    DEGREE POINT_NUM'

  do test = 1, test_num

    dim_num = dim_num_test(test)

    write ( *, '(a)' ) ' '

    do rule = 0, 5

      call gm_rule_size ( rule, dim_num, point_num )
      degree = 2 * rule + 1

      write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) dim_num, rule, degree, point_num

    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests GM_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  GM_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the DIM_NUM dimensional simplex,'
  write ( *, '(a)' ) '  using a rule of in index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'

  dim_num = 3
  rule = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we use DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  RULE = ', rule
  write ( *, '(a,i8)' ) '  DEGREE = ', 2 * rule + 1

  call gm_rule_size ( rule, dim_num, point_num )

  allocate ( w(1:point_num) )
  allocate ( x(1:dim_num,1:point_num) )

  call gm_rule_set ( rule, dim_num, point_num, w, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     POINT        W             X             Y             Z'
  write ( *, '(a)' ) ' '

  do point = 1, point_num
    write ( *, '(2x,i8,2x,f12.6,2x,f12.6,2x,f12.6,2x,f12.6)' ) &
      point, w(point), x(1:dim_num,point)
  end do

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests GM_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), dimension ( test_num ) :: dim_num_test = (/ &
    2, 3, 5, 10 /)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ) w_sum
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  GM_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the DIM_NUM dimensional simplex,'
  write ( *, '(a)' ) '  using a rule of in index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute various rules, and simply'
  write ( *, '(a)' ) '  report the number of points, and the sum of weights.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   DIM_NUM      RULE    POINT_NUM  WEIGHT SUM'

  do test = 1, test_num

    dim_num = dim_num_test(test)

    write ( *, '(a)' ) ' '

    do rule = 0, 5

      call gm_rule_size ( rule, dim_num, point_num )

      allocate ( w(1:point_num) )
      allocate ( x(1:dim_num,1:point_num) )

      call gm_rule_set ( rule, dim_num, point_num, w, x )

      w_sum = sum ( w(1:point_num) )

      write ( *, '(2x,i8,2x,i8,2x,i8,2x,g24.16)' ) &
        dim_num, rule, point_num, w_sum

      deallocate ( w )
      deallocate ( x )

    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests GM_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  character ( len = 12 ) w_file
  integer ( kind = 4 ) w_unit
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  character ( len = 12 ) x_file
  integer ( kind = 4 ) x_unit

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  GM_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the DIM_NUM dimensional simplex,'
  write ( *, '(a)' ) '  using a rule of in index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we write a rule to a file.'

  dim_num = 3
  rule = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we use DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  RULE = ', rule
  write ( *, '(a,i8)' ) '  DEGREE = ', 2 * rule + 1

  call gm_rule_size ( rule, dim_num, point_num )

  allocate ( w(1:point_num) )
  allocate ( x(1:dim_num,1:point_num) )

  call gm_rule_set ( rule, dim_num, point_num, w, x )

  call get_unit ( w_unit )

  write ( w_file, '(a2,i1,a,i1,a7)' ) 'gm', rule, '_', dim_num, 'd_w.txt'

  open ( unit = w_unit, file = w_file, status = 'replace' )

  do point = 1, point_num
    write ( w_unit, '(f20.16)' ) w(point)
  end do

  close ( unit = w_unit )

  call get_unit ( x_unit )

  write ( x_file, '(a2,i1,a,i1,a7)' ) 'gm', rule, '_', dim_num, 'd_x.txt'

  open ( unit = x_unit, file = x_file, status = 'replace' )

  do point = 1, point_num
    write ( x_unit, '(3f20.16)' ) x(1:dim_num,point)
  end do

  close ( unit = x_unit )

  write ( *, '(a,i2,a)' ) '  Wrote rule ', rule, ' to "' &
    // trim ( w_file ) // '" and "' // trim ( x_file ) // '".'

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests GM_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 5

  integer ( kind = 4 ) degree
  integer ( kind = 4 ), parameter :: degree_max = 4
  integer ( kind = 4 ) expon(dim_num)
  integer ( kind = 4 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: mono
  logical more
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) quad
  real ( kind = 8 ) quad_error
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), parameter :: rule_max = 3
  integer ( kind = 4 ) t
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  GM_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the DIM_NUM dimensional simplex,'
  write ( *, '(a)' ) '  using a rule of in index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, look at all the monomials up to'
  write ( *, '(a)' ) '  some maximum degree, choose a few low order rules'
  write ( *, '(a)' ) '  and determine the quadrature error for each.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we use DIM_NUM = ', dim_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule     Order     Quad_Error'
  write ( *, '(a)' ) ' '

  do degree = 0, degree_max

    more = .false.

    do

      call comp_next ( degree, dim_num, expon, more, h, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i1,a,i1,a,i1,a,i1,a,i1)' ) &
        '  F(X) = X1^', expon(1), ' * X2^', expon(2), ' * X3^', expon(3), &
              ' * X4^', expon(4), ' * X5^', expon(5)

      write ( *, '(a)' ) ' '

      do rule = 0, rule_max

        call gm_rule_size ( rule, dim_num, point_num )

        allocate ( mono(1:point_num) )
        allocate ( w(1:point_num) )
        allocate ( x(1:dim_num,1:point_num) )

        call gm_rule_set ( rule, dim_num, point_num, w, x )

        call simplex_unit_monomial_quadrature ( dim_num, expon, point_num, &
          x, w, quad_error )

        write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) rule, point_num, quad_error

        deallocate ( mono )
        deallocate ( w )
        deallocate ( x )

      end do

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  return
end
