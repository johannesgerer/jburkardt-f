program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_OPTIMIZATION_PRB.
!
!  Discussion:
!
!    TEST_OPTIMIZATION_PRB calls the TEST_OPTIMIZATION tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp (  )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_OPTIMIZATION_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_OPTIMIZATION library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_OPTIMIZATION_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp (  )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 simply prints the title of each problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2011
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
  write ( *, '(a)' ) '  For each problem, print the title.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem  Title'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(2x,i7,2x,a)' ) problem, title

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 samples the function at 1,000 points and prints the minimum.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 1000

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) f_min
  integer ( kind = 4 ) know
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) seed
  character ( len = 80 ) title
  real ( kind = 8 ) x(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For each problem, using dimension M = 2'
  write ( *, '(a)' ) '  sample the function at N = 1000 points,'
  write ( *, '(a)' ) '  and print the minimum and maximum.'

  seed = 123456789
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Problem     Minimum  Sample Minimum  Sample Maximum'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    know = 0
    call p00_sol ( problem, m, know, x )
    if ( know /= 0 ) then
      call p00_f ( problem, m, 1, x, f )
      f_min = f(1)
    end if

    call p00_ab ( problem, m, a, b )
    call r8col_uniform ( m, n, a, b, seed, x )
    call p00_f ( problem, m, n, x, f )
    if ( know /= 0 ) then
      write ( *, '(2x,i7,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        problem, f_min, minval ( f(1:n) ), maxval ( f(1:n) )
    else
      write ( *, '(2x,i7,2x,14x,2x,g14.6,2x,g14.6)' ) &
        problem,         minval ( f(1:n) ), maxval ( f(1:n) )
    end if

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tries Compass Search on each problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 1000

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) delta_init
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) know
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) seed
  character ( len = 80 ) title
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)

  delta_tol = 0.000001D+00
  k_max = 20000

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For each problem, using dimension M = 2'
  write ( *, '(a)' ) '  try compass search.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    seed = 123456789

    call p00_ab ( problem, m, a, b )
    call r8col_uniform ( m, 1, a, b, seed, x0 )
    call p00_f ( problem, m, 1, x0, fx )
    delta_init = 0.3D+00 * sqrt ( sum ( x0(1:m)**2 ) ) / real ( m, kind = 8 )
    delta_init = max ( delta_init, 1000.0D+00 * delta_tol )
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a,g14.6)' ) '  Problem ', problem, '  DELTA_INIT = ', delta_init
    write ( *, '(a,2x,g14.6,2x,g14.6,2x,g14.6)' ) '  Initial:', x0, fx
    call p00_compass_search ( problem, m, x0, delta_tol, delta_init, &
      k_max, x, fx, k )
    write ( *, '(a,2x,g14.6,2x,g14.6,2x,g14.6,2x,a,i8)' ) '  Final:  ', x, fx, '  Steps = ', k

    know = 0
    do
      call p00_sol ( problem, m, know, x )
      if ( know == 0 ) then
        exit
      end if
      call p00_f ( problem, m, 1, x, fx )
      write ( *, '(a,2x,g14.6,2x,g14.6,2x,g14.6)' ) '  Exact:  ', x, fx
    end do

  end do

  return
end