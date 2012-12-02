program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_NINT_PRB.
!
!  Discussion:
!
!    TEST_NINT_PRB demonstrates the TEST_NINT integration test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NINT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_NINT library.'

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
  write ( *, '(a)' ) 'TEST_NINT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 retrieves and prints the name for each problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) name
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  GET_PROBLEM_NUM returns the number of problems.'
  write ( *, '(a)' ) '  P00_NAME(#) returns the name for problem #.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use these two routines to print a directory'
  write ( *, '(a)' ) '  of all the problems.'

  call get_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of problems available is ', problem_num
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_name ( problem, name )
    write ( *, '(2x,i8,2x,a)' ) problem, '"' // trim ( name ) // '".'

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 just prints out the title information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  GET_PROBLEM_NUM returns the number of problems.'
  write ( *, '(a)' ) '  P00_TITLE(#) prints the title for problem #.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use these two routines to print a directory'
  write ( *, '(a)' ) '  of all the problems.'

  call get_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of problems available is ', problem_num
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_title ( problem )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 applies a composite midpoint rule to box regions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) n
  integer ( kind = 4 ) problem_num
  character ( len = 10 ) region
  real ( kind = 8 ) result
  integer ( kind = 4 ) sub_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Use a simple product rule on box regions.'
  write ( *, '(a)' ) '  Use a fixed spatial dimension.'

  call get_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Prob   Dim  Subs       Approx          Exact          Error'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    dim_num = 3
!
!  Set problem data to default values.
!
    call p00_default ( problem, dim_num )
!
!  Get the region type.
!
    call p00_region ( problem, region )

    if ( region(1:3) == 'box' .or. region(1:3) == 'BOX' ) then

      do sub_num = 1, 5, 2

        call p00_box_gl05 ( problem, dim_num, sub_num, result )

        call p00_exact ( problem, dim_num, exact )

        if ( exact == huge ( exact ) ) then

          write ( *, '(2x,i4,2x,i4,2x,i4,2x,g14.6,2x,a14,2x,a14)' ) &
            problem, dim_num, sub_num, result, &
            '--------------', '--------------'

        else

          error = abs ( result - exact )

          write ( *, '(2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
            problem, dim_num, sub_num, result, exact, error

        end if

      end do

      write ( *, '(a)' ) ' '

    end if

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 applies a Monte Carlo rule to box regions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) point_num
  character ( len = 10 ) region
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Use a Monte Carlo rule on box regions.'
  write ( *, '(a)' ) '  Use a fixed spatial dimension.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeatedly multiply the number of points by 16.'

  call get_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Prob   Dim    Points     Approx        Exact          Error'
  write ( *, '(a)' ) ' '

  dim_num = 3

  do problem = 1, problem_num
!
!  Set problem data to default values.
!
    call p00_default ( problem, dim_num )
!
!  Get region type.
!
    call p00_region ( problem, region )

    if ( region(1:3) == 'box' ) then

      do i = 1, 5

        if ( i == 1 ) then
          point_num = 1
        else
          point_num = 16 * point_num
        end if

        seed = 123456789
        call random_initialize ( seed )

        call p00_box_mc ( problem, dim_num, point_num, result )

        call p00_exact ( problem, dim_num, exact )

        if ( exact == huge ( exact ) ) then

          write ( *, '(2x,i4,2x,i4,i10,g14.6,a14,2x,a14)' ) &
            problem, dim_num, point_num, result, &
            '--------------', '--------------'

        else

          error = abs ( result - exact )

          write ( *, '(2x,i4,2x,i4,i10,g14.6,g14.6,2x,g14.6)' ) &
            problem, dim_num, point_num, result, exact, error

        end if

      end do

      write ( *, '(a)' ) ' '

    end if

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 demonstrates how a base point can be moved.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  integer ( kind = 4 ) run
  integer ( kind = 4 ) test
  integer ( kind = 4 ), dimension ( test_num ) :: problem_index = (/ &
    16, 17, 18, 19, 31 /)
  real ( kind = 8 ) z(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Demonstrate problems that use a "base point"'
  write ( *, '(a)' ) '  by moving the base point around.'
  write ( *, '(a)' ) '  Use a Monte Carlo rule on box regions.'
  write ( *, '(a)' ) '  Use a fixed spatial dimension.'

  do test = 1, test_num

    problem = problem_index(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Problem number = ', problem

    call p00_default ( problem, dim_num )

    do run = 1, 3

      call p00_r8vec ( problem, 'R', 'Z', dim_num, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i1)') '  Run number ', run
      write ( *, '(a,2f10.4)' ) '  Basis point Z = ', z(1:dim_num)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) &
        '  Prob   Dim    Points    Approx        Exact           Error'
      write ( *, '(a)' ) ' '

      do i = 1, 3

        if ( i == 1 ) then
          point_num = 10
        else if ( i == 2 ) then
          point_num = 1000
        else if ( i == 3 ) then
          point_num = 100000
        end if

        call p00_box_mc ( problem, dim_num, point_num, result )

        call p00_exact ( problem, dim_num, exact )

        if ( exact == huge ( exact ) ) then

          write ( *, '(2x,i4,2x,i4,i10,g14.6,a14,2x,a14)' ) &
            problem, dim_num, point_num, result, &
            '--------------', '--------------'

        else

          error = abs ( result - exact )

          write ( *, '(2x,i4,2x,i4,i10,g14.6,g14.6,2x,g14.6)' ) &
            problem, dim_num, point_num, result, exact, error

        end if

      end do

      write ( *, '(a)' ) ' '

    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 applies a composite midpoint rule for increasing spatial dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) calls
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) result
  integer ( kind = 4 ) sub_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Use a simple product rule on a box region.'
  write ( *, '(a)' ) '  Use a fixed problem;'
  write ( *, '(a)' ) '  Let the spatial dimension increase.'

  problem = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Prob   Dim  Subs    Approx         Exact        Error        Calls'
  write ( *, '(a)' ) ' '

  do dim_num = 1, 6

    call p00_default ( problem, dim_num )

    do sub_num = 1, 5, 2

      call p00_i4 ( problem, 'S', '#', 0 )

      call p00_box_gl05 ( problem, dim_num, sub_num, result )

      call p00_i4 ( problem, 'G', '#', calls )

      call p00_exact ( problem, dim_num, exact )

      if ( exact == huge ( exact ) ) then

        write ( *, '(2x,i4,2x,i4,2x,i4,g14.6,a14,a10,i12)' ) &
          problem, dim_num, sub_num, result, &
          '--------------', '----------', calls

      else

        error = abs ( result - exact )

        write ( *, '(2x,i4,2x,i4,2x,i4,g14.6,g14.6,f10.6,i12)' ) &
          problem, dim_num, sub_num, result, exact, error, calls

      end if

    end do

    write ( *, '(a)' ) ' '

  end do

  return
end
