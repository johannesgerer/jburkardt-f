program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_INT_PRB.
!
!  Discussion:
!
!    TEST_INT_PRB demonstrates the use of the TEST_INT integration
!    test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_INT library.'

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
  write ( *, '(a)' ) 'TEST_INT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 applies a composite midpoint rule to finite interval 1D problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) int_log
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Composite midpoint rule,'
  write ( *, '(a)' ) '  for 1D finite interval problems.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem       Exact'
  write ( *, '(a)' ) '         Ints   Approx       Error'
!
!  Pick a problem.
!
  do prob = 1, prob_num

    call p00_exact ( prob, exact )

    write ( *, '(a)' ) ' '
    write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
!
!  Pick a number of subintervals.
!
    do int_log = 0, 7

      int_num = 2**int_log

      call p00_midpoint ( prob, int_num, result )

      error = abs ( exact - result )

      write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

    end do

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 applies a composite Simpson rule to finite interval 1D problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) int_log
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Composite Simpson rule,'
  write ( *, '(a)' ) '  for 1D finite interval problems.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem       Exact'
  write ( *, '(a)' ) '         Ints   Approx       Error'
!
!  Pick a problem.
!
  do prob = 1, prob_num
!
!  Some problems have singularities that kill the calculation.
!
    call p00_exact ( prob, exact )

    write ( *, '(a)' ) ' '
    write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
!
!  Pick a number of subintervals.
!
    do int_log = 0, 7

      int_num = 2**int_log

      call p00_simpson ( prob, int_num, result )

      error = abs ( exact - result )

      write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

    end do

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 applies a Monte Carlo rule to finite interval 1D problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) int_log
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Monte Carlo rule,'
  write ( *, '(a)' ) '  for 1D finite interval problems.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem       Exact'
  write ( *, '(a)' ) '         Pts    Approx       Error'
!
!  Pick a problem.
!
  do prob = 1, prob_num

    call p00_exact ( prob, exact )

    write ( *, '(a)' ) ' '
    write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
!
!  Pick a number of points.
!
    do int_log = 0, 7

      int_num = 2**int_log

      call p00_montecarlo ( prob, int_num, result )

      error = abs ( exact - result )

      write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 applies a composite Gauss-Legendre rule to finite interval 1D problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) int_log
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Use a composite 4 point Gauss-Legendre rule,'
  write ( *, '(a)' ) '  for 1D finite interval problems.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem       Exact'
  write ( *, '(a)' ) '         Ints   Approx	   Error'
!
!  Pick a problem.
!
  do prob = 1, prob_num

    call p00_exact ( prob, exact )

    write ( *, '(a)' ) ' '
    write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
!
!  Pick a number of subintervals.
!
    do int_log = 0, 7

      int_num = 2**int_log

      call p00_gauss_legendre ( prob, int_num, result )

      error = abs ( exact - result )

      write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

    end do

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 applies a composite trapezoid rule to finite interval 1D problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) int_log
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Composite trapezoid rule,'
  write ( *, '(a)' ) '  for 1D finite interval problems.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem       Exact'
  write ( *, '(a)' ) '         Ints   Approx       Error'
!
!  Pick a problem.
!
  do prob = 1, prob_num

    call p00_exact ( prob, exact )

    write ( *, '(a)' ) ' '
    write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
!
!  Pick a number of subintervals.
!
    do int_log = 0, 7

      int_num = 2**int_log

      call p00_trapezoid ( prob, int_num, result )

      error = abs ( exact - result )

      write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 applies a Halton sequence rule to finite interval 1D problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) int_log
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Halton sequence rule,'
  write ( *, '(a)' ) '  for 1D finite interval problems.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem       Exact'
  write ( *, '(a)' ) '         Pts    Approx       Error'
!
!  Pick a problem.
!
  do prob = 1, prob_num

    call p00_exact ( prob, exact )

    write ( *, '(a)' ) ' '
    write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
!
!  Pick a number of points.
!
    do int_log = 0, 7

      int_num = 2**int_log

      call p00_halton ( prob, int_num, result )

      error = abs ( exact - result )

      write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

    end do

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 applies an evenly spaced point rule to finite interval 1D problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) int_log
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Evenly spaced point sequence rule,'
  write ( *, '(a)' ) '  for 1D finite interval problems.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem       Exact'
  write ( *, '(a)' ) '         Pts    Approx       Error'
!
!  Pick a problem.
!
  do prob = 1, prob_num

    call p00_exact ( prob, exact )

    write ( *, '(a)' ) ' '
    write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
!
!  Pick a number of points.
!
    do int_log = 0, 7

      int_num = 2**int_log

      call p00_even ( prob, int_num, result )

      error = abs ( exact - result )

      write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

    end do

  end do

  return
end
