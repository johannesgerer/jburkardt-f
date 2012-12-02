program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_INT_HERMITE_PRB.
!
!  Discussion:
!
!    TEST_INT_HERMITE_PRB demonstrates the use of the TEST_INT_HERMITE
!    integration test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INT_HERMITE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_INT_HERMITE library.'

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
  write ( *, '(a)' ) 'TEST_INT_HERMITE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests P00_PROBLEM_NUM and P00_TITLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  P00_PROBLEM_NUM returns the number of problems.'
  write ( *, '(a)' ) '  P00_TITLE returns the title of a problem.'

  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  P00_PROBLEM_NUM: number of problems is ', problem_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem       Title'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(2x,i8,2x,a)' ) problem, '"' // trim ( title ) // '"'

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests P00_EXACT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) m
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  P00_EXACT returns the "exact" integral.'

  call p00_problem_num ( problem_num )

  m = 4
  call p06_param ( 'S', 'M', m )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem       EXACT'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_exact ( problem, exact )

    write ( *, '(2x,i8,2x,g24.16)' ) problem, exact

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests P00_GAUSS_HERMITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  integer ( kind = 4 ) m
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_log
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  P00_GAUSS_HERMITE applies a Gauss-Hermite rule'
  write ( *, '(a)' ) '  to estimate an integral on (-oo,+oo).'

  call p00_problem_num ( problem_num )

  m = 4
  call p06_param ( 'S', 'M', m )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   Problem     Order      Estimate        Exact           Error'

  do problem = 1, problem_num

    call p00_exact ( problem, exact )

    order = 1

    write ( *, '(a)' ) ' '

    do order_log = 0, 6

      call p00_gauss_hermite ( problem, order, estimate )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        problem, order, estimate, exact, error

      order = order * 2

    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests P00_TURING.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  real ( kind = 8 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_log
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) test
  real ( kind = 8 ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  P00_TURING applies a Turing procedure'
  write ( *, '(a)' ) '  to estimate an integral on (-oo,+oo).'

  call p00_problem_num ( problem_num )

  m = 4
  call p06_param ( 'S', 'M', m )

  do test = 1, 2

    if ( test == 1 ) then
      tol = 1.0D-04
    else if ( test == 2 ) then
      tol = 1.0D-07
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Using a tolerance of TOL = ', tol
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '   Problem      H              N      Estimate        Exact           Error'

    do problem = 1, problem_num

      call p00_exact ( problem, exact )

      h = 1.0D+00

      write ( *, '(a)' ) ' '

      do order_log = 0, 6

        call p00_turing ( problem, h, tol, n, estimate )

        order = 2 * n + 1

        error = abs ( exact - estimate )

        write ( *, '(2x,i8,2x,f10.6,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          problem, h, order, estimate, exact, error

        h = h / 2.0D+00

      end do

    end do

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests P00_GAUSS_HERMITE against the polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  integer ( kind = 4 ) m
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_log
  integer ( kind = 4 ) problem

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  P00_GAUSS_HERMITE applies a Gauss-Hermite rule to'
  write ( *, '(a)' ) '  estimate the integral x^m exp(-x*x) over (-oo,+oo).'

  problem = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         M     Order      Estimate        Exact           Error'
  do m = 0, 6

    call p06_param ( 'S', 'M', m )

    call p00_exact ( problem, exact )

    write ( *, '(a)' ) ' '

    do order = 1, 3 + ( m / 2 )

      call p00_gauss_hermite ( problem, order, estimate )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        m, order, estimate, exact, error

    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests P00_MONTE_CARLO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  integer ( kind = 4 ) m
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_log
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  P00_MONTE_CARLO applies a Monte Carlo rule'
  write ( *, '(a)' ) '  to estimate Hermite-weighted integrals on (-oo,+oo).'

  call p00_problem_num ( problem_num )

  m = 4
  call p06_param ( 'S', 'M', m )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   Problem     Order      Estimate        Exact           Error'

  do problem = 1, problem_num

    call p00_exact ( problem, exact )

    order = 128

    write ( *, '(a)' ) ' '

    do order_log = 0, 6

      call p00_monte_carlo ( problem, order, estimate )

      error = abs ( exact - estimate )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        problem, order, estimate, exact, error

      order = order * 4

    end do

  end do

  return
end
