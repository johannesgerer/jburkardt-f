program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_TRI_INT_PRB.
!
!  Discussion:
!
!    TEST_TRI_INT_PRB demonstrates the use of the TEST_TRI_INT
!    integration test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TRI_INT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_TRI_INT library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TRI_INT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests GET_PROB_NUM and P00_TITLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) problem
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  GET_PROB_NUM reports the number of problems.'
  write ( *, '(a)' ) '  P00_TITLE returns a title for each problem.'

  call get_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of problems available is ', prob_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The problem titles:'
  write ( *, '(a)' ) ' '

  do problem = 1, prob_num

    call p00_title ( problem, title )

    write ( *, '(2x,i8,2x,a)' ) problem, trim ( title )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests P00_MONTE_CARLO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) abs_error
  real ( kind = 8 ), parameter :: exact = 1.0D+00
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log
  integer ( kind = 4 ), parameter :: n_log_max = 15
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  P00_MONTE_CARLO applies a Monte Carlo rule.'

  call get_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem            Exact         Seed'
  write ( *, '(a)' ) '           Pts       Approx        Error'
  write ( *, '(a)' ) ' '
!
!  Pick a problem.
!
  do problem = 1, prob_num

    call p00_title ( problem, title )

    seed = 123456789
!
!  Call RANDOM_INITIALIZE in case we are using the FORTRAN90
!  random number generator inside of P00_MONTE_CARLO!
!
    call random_initialize ( seed )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(i8,6x,2x,f12.6,2x,i12)' ) problem, exact, seed
    write ( *, '(a)' ) ' '
!
!  Pick a number of points.
!
    do n_log = 0, n_log_max

      n = 2**n_log

      call p00_monte_carlo ( problem, n, seed, result )

      abs_error = abs ( exact - result )

      write ( *, '(6x,i8,2x,f12.6,2x,f12.6)' ) n, result, abs_error

    end do

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests P00_VERTEX_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) abs_error
  real ( kind = 8 ), parameter :: exact = 1.0D+00
  integer ( kind = 4 ) level
  integer ( kind = 4 ) :: level_max = 4
  integer ( kind = 4 ) n
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result
  integer ( kind = 4 ) singularity
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  P00_VERTEX_SUB applies a vertex rule with subdivision.'

  call get_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem            Exact'
  write ( *, '(a)' ) '           Pts       Approx        Error'
  write ( *, '(a)' ) ' '
!
!  Pick a problem.
!
  do problem = 1, prob_num

    call p00_title ( problem, title )
    call p00_singularity ( problem, singularity )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(i8,6x,2x,f12.6)' ) problem, exact
    write ( *, '(a)' ) ' '

    if ( singularity == 1 ) then
      write ( *, '(a)' ) '  Skip this problem, it has vertex singularities.'
    else if ( singularity == 2 ) then
      write ( *, '(a)' ) '  Skip this problem, it has edge singularities.'
    else if ( singularity == 3 ) then
      write ( *, '(a)' ) '  Skip this problem, it has internal singularities.'
    else
!
!  Pick a number of points.
!
      n = 0
      result = 0.0D+00

      do level = 0, 4

        call p00_vertex_sub ( problem, level, n, result )

        abs_error = abs ( exact - result )

        write ( *, '(6x,i8,2x,f12.6,2x,f12.6)' ) n, result, abs_error

      end do

    end if

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests P00_WANDZURA05_SUB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) abs_error
  real ( kind = 8 ), parameter :: exact = 1.0D+00
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) result
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_max = 5
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  P00_WANDZURA05_SUB applies a Wandzura rule with subdivision.'

  call get_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem            Exact'
  write ( *, '(a)' ) '           Pts       Approx        Error'
  write ( *, '(a)' ) ' '
!
!  Pick a problem.
!
  do problem = 1, prob_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(i8,6x,2x,f12.6)' ) problem, exact
    write ( *, '(a)' ) ' '
!
!  Pick a number of points.
!
    do test = 0, test_max

      level = 2**test

      call p00_wandzura05_sub ( problem, level, n, result )

      abs_error = abs ( exact - result )

      write ( *, '(6x,i8,2x,f12.6,2x,f12.6)' ) n, result, abs_error

      if ( abs_error < 1.0D-06 ) then
        write ( *, * ) '                            Accuracy acceptable'
        exit
      end if

    end do

  end do

  return
end
