program main

!*****************************************************************************80
!
!! MAIN calls a set of problems for LATTICE_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATTICE_RULE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LATTICE_RULE library.'
 
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
  call test10 ( )

  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATTICE_RULE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
   
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests FIBONACCI_LATTICE_Q;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 78-80, page 145.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  FIBONACCI_LATTICE_Q applies a Fibonacci lattice rule'
  write ( *, '(a)' ) '  to integrate a function over the unit square.'
  write ( *, '(a)' ) '  These Fibonacci rules are only available in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     K     M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    call fibonacci_lattice_q ( k, f_01_2d, quad )

    error = abs ( exact - quad )
    m = fibonacci ( k )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests FIBONACCI_LATTICE_T;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 78-80, page 145.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  FIBONACCI_LATTICE_T applies a '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  symmetric Fibonacci lattice rule'
  write ( *, '(a)' ) '  ---------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  to integrate a function over the unit square.'
  write ( *, '(a)' ) '  These Fibonacci rules are only available in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     K     M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    call fibonacci_lattice_t ( k, f_01_2d, quad )

    error = abs ( exact - quad )
    m = fibonacci ( k )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests FIBONACCI_LATTICE_B;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 78-80, page 145.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  FIBONACCI_LATTICE_B applies an '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  optimal Fibonacci lattice rule'
  write ( *, '(a)' ) '  -------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  to integrate a function over the unit square.'
  write ( *, '(a)' ) '  These Fibonacci rules are only available in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     K     M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    call fibonacci_lattice_b ( k, f_01_2d, quad )

    error = abs ( exact - quad )
    m = fibonacci ( k )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests FIBONACCI_LATTICE_Q1;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  FIBONACCI_LATTICE_Q1 applies a Fibonacci lattice rule'
  write ( *, '(a)' ) '  to integrate a function over the unit square.'
  write ( *, '(a)' ) '  A nonlinear coordinate transformation is applied.'
  write ( *, '(a)' ) '  These Fibonacci rules are only available in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     K     M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    call fibonacci_lattice_q1 ( k, f_01_2d, quad )

    error = abs ( exact - quad )
    m = fibonacci ( k )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests FIBONACCI_LATTICE_Q2;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  FIBONACCI_LATTICE_Q2 applies a Fibonacci lattice rule'
  write ( *, '(a)' ) '  to integrate a function over the unit square.'
  write ( *, '(a)' ) '  A nonlinear coordinate transformation is applied.'
  write ( *, '(a)' ) '  These Fibonacci rules are only available in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     K     M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    call fibonacci_lattice_q2 ( k, f_01_2d, quad )

    error = abs ( exact - quad )
    m = fibonacci ( k )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests FIBONACCI_LATTICE_Q3;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  FIBONACCI_LATTICE_Q3 applies a Fibonacci lattice rule'
  write ( *, '(a)' ) '  to integrate a function over the unit square.'
  write ( *, '(a)' ) '  A nonlinear coordinate transformation is applied.'
  write ( *, '(a)' ) '  These Fibonacci rules are only available in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     K     M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    call fibonacci_lattice_q3 ( k, f_01_2d, quad )

    error = abs ( exact - quad )
    m = fibonacci ( k )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests LATTICE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 18.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) prime
  real ( kind = 8 ) quad
  integer ( kind = 4 ) z(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  LATTICE applies a lattice rule to integrate'
  write ( *, '(a)' ) '  a function over the unit hypercube.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  The lattice rule order M will vary.'

  z(1:dim_num) = (/ 1, 2 /)
  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  call i4vec_print ( dim_num, z, '  The lattice generator vector:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  M  EXACT    ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    m = prime ( 3 * i )

    call lattice ( dim_num, m, z, f_01_2d, quad )

    exact = e_01_2d ( dim_num, a, b )

    error = abs ( exact - quad )

    write ( *, '(2i8,3f10.6)' ) i, m, exact, quad, error

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests LATTICE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 18.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: m = 53
  real ( kind = 8 ) quad
  integer ( kind = 4 ) z(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  LATTICE applies a lattice rule to integrate'
  write ( *, '(a)' ) '  a function over the unit hypercube.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  The lattice rule order M is fixed at ', m
  write ( *, '(a)' ) '  The lattice generator vector Z will vary.'

  z(1) = 1
  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  M   Z(1)  Z(2) EXACT    ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do i = 1, m-1

    z(2) = i

    call lattice ( dim_num, m, z, f_01_2d, quad )

    exact = e_01_2d ( dim_num, a, b )

    error = abs ( exact - quad )

    write ( *, '(i8,2i8,3f10.6)' ) m, z(1), z(2), exact, quad, error

  end do

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests LATTICE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 18.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad
  integer ( kind = 4 ) z(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  LATTICE is a lattice rule for periodic functions.'
  write ( *, '(a)' ) '  However, we apply it to a nonperiodic function'
  write ( *, '(a)' ) '  just to see how it does.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  z(1:dim_num) = (/ 1, 2 /)
  call i4vec_print ( dim_num, z, '  The lattice generator vector:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       K       M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    m = fibonacci ( k )

    call lattice ( dim_num, m, z, f_01_2d, quad )

    error = abs ( exact - quad )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests LATTICE_NP0;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 78-80, pages 32-40, 145-147.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad
  integer ( kind = 4 ) z(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  LATTICE_NP0 applies a lattice rule to a'
  write ( *, '(a)' ) '  nonperiodic function by reflecting the function'
  write ( *, '(a)' ) '  about the midpoint and averaging.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  z(1:dim_num) = (/ 1, 2 /)
  call i4vec_print ( dim_num, z, '  The lattice generator vector:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       K       M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    m = fibonacci ( k )

    call lattice_np0 ( dim_num, m, z, f_01_2d, quad )

    error = abs ( exact - quad )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests LATTICE_NP1;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 78-80, pages 32-40, 145-147.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad
  integer ( kind = 4 ) z(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  LATTICE_NP1 applies a lattice rule to a'
  write ( *, '(a)' ) '  nonperiodic function using a nonlinear transformation,'
  write ( *, '(a)' ) '  to integrate a function over the unit square.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  z(1:dim_num) = (/ 1, 2 /)
  call i4vec_print ( dim_num, z, '  The lattice generator vector:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       K       M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 3, 18

    m = fibonacci ( k )

    call lattice_np1 ( dim_num, m, z, f_01_2d, quad )

    error = abs ( exact - quad )

    write ( *, '(2i8,3f10.6)' ) k, m, exact, quad, error

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests MONTE_CARLO;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ), external :: e_01_2d
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), external :: f_01_2d
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) quad

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  MONTE_CARLO applies a Monte Carlo scheme'
  write ( *, '(a)' ) '  to estimate the integral of a function'
  write ( *, '(a)' ) '  over the unit hypercube.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  exact = e_01_2d ( dim_num, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     M  EXACT     ESTIMATE  ERROR'
  write ( *, '(a)' ) ' '

  do k = 2, 5

    m = 10**k

    call monte_carlo ( dim_num, m, f_01_2d, quad )

    error = abs ( exact - quad )

    write ( *, '(i8,3f10.6)' ) m, exact, quad, error

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests LATTICE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 18.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ) z(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  LATTICE_PRINT prints out the lattice generated'
  write ( *, '(a)' ) '  by a single generator vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num

  z(1:dim_num) = (/ 1, 3 /)

  call i4vec_print ( dim_num, z, '  The generator vector:' )

  call lattice_print ( dim_num, m, z, '  The total lattice:' )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests FIND_Z20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 18.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) z(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  FIND_Z20 finds the optimal lattice generator Z'
  write ( *, '(a)' ) '  with Fourier coefficient smoothness ALPHA = 2,'
  write ( *, '(a)' ) '  and copy exponent 0,'
  write ( *, '(a)' ) '  for a rank 1 "method of good lattice points" rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     M      Z(1)  Z(2)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (M = Fibonacci)'
  write ( *, '(a)' ) ' '

  do i = 3, 10

    m = fibonacci ( i )

    call find_z20 ( dim_num, m, z )

    write ( *, '(i8,4x,2i8)' ) m, z(1:dim_num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (M = 2**K)'
  write ( *, '(a)' ) ' '

  do i = 2, 10

    m = 2**i

    call find_z20 ( dim_num, m, z )

    write ( *, '(i8,4x,2i8)' ) m, z(1:dim_num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (M = 3*2**K)'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    m = 3*2**i

    call find_z20 ( dim_num, m, z )

    write ( *, '(i8,4x,2i8)' ) m, z(1:dim_num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (M = Prime)'
  write ( *, '(a)' ) ' '

  do i = 3, 10

    m = prime ( 10 * i )

    call find_z20 ( dim_num, m, z )

    write ( *, '(i8,4x,2i8)' ) m, z(1:dim_num)

  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests FIBONACCI_LATTICE_Q_NODES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2005.
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ian Sloan, Stephen Joe,
!    Lattice Methods for Multiple Integration,
!    Oxford, 1994, page 18.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) fibonacci
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  k = 12
  m = fibonacci ( k ) 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  FIBONACCI_LATTICE_Q_NODES...'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  The Fibonacci index K =   ', k
  write ( *, '(a,i8)' ) '  The Fibonacci value M =   ', m

  allocate ( x(1:dim_num,1:m) )

  call fibonacci_lattice_q_nodes ( k, x )

  call r8mat_transpose_print ( dim_num, m, x, '  The Fibonacci lattice nodes:' )

  deallocate ( x )

  return
end
