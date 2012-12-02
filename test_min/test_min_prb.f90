program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_MIN_PRB.
!
!  Discussion:
!
!    TEST_MIN_PRB calls the TEST_MIN tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp (  )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MIN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_MIN library.'

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
  write ( *, '(a)' ) 'TEST_MIN_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 prints the title of each problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) problem
  character ( len = 50 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For each problem, print the title.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem Title'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(2x,i8,2x,a)' ) problem, trim ( title )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 evaluates the objective function at each starting point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f_sol
  real ( kind = 8 ) f_start
  integer ( kind = 4 ) know
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) problem
  character ( len = 50 ) title
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For each problem, evaluate the function'
  write ( *, '(a)' ) '  at the starting point and the solution.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(a)' ) ' '
 
    call p00_start ( problem, x )

    call p00_f ( problem, x, f_start )

    write ( *, '(4x,a,g16.8)' ) 'F(X_START)=', f_start

    call p00_sol ( problem, know, x )

    if ( 0 < know ) then
      call p00_f ( problem, x, f_sol )
      write ( *, '(4x,a,g16.8)' ) 'F(X_SOL)=  ', f_sol
    end if

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 compares the exact and approximate first derivatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) f1_dif
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) problem
  character ( len = 50 ) title
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For each problem, compare the exact and'
  write ( *, '(a)' ) '  approximate gradients at the starting point.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )

    call p00_start ( problem, x )

    call p00_f1 ( problem, x, f1 )

    call p00_f1_dif ( problem, x, f1_dif )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a)' ) 'X'
    write ( *, '(4x,5g16.8)' ) x
    write ( *, '(2x,a)' ) 'F''(X) (exact)'
    write ( *, '(4x,5g16.8)' ) f1
    write ( *, '(2x,a)' ) 'F''(X) (difference)'
    write ( *, '(4x,5g16.8)' ) f1_dif

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 compares the exact and approximate second derivatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) f2_dif
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) problem
  character ( len = 50 ) title
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For each problem, compare the exact and'
  write ( *, '(a)' ) '  approximate second derivatives at the starting point.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )

    call p00_start ( problem, x )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  X:'
    write ( *, '(4x,5g16.8)' ) x

    call p00_f2 ( problem, x, f2 )

    write ( *, '(a)' ) '  F"(X) (exact):'
    write ( *, '(4x,6g13.5)' ) f2

    call p00_f2_dif ( problem, x, f2_dif )

    write ( *, '(a)' ) '  F"(X) (difference):'
    write ( *, '(4x,6g13.5)' ) f2_dif

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 carries out a simple bisection method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) fd
  real ( kind = 8 ) fe
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: max_step = 10
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) problem
  character ( len = 50 ) title
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) xd
  real ( kind = 8 ) xe

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For each problem, take a few steps of '
  write ( *, '(a)' ) '  the bisection method.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )

    call p00_interval ( problem, xa, xc )
    xb = 0.5D+00 * ( xa + xc )
    call p00_f ( problem, xa, fa )
    call p00_f ( problem, xc, fc )
    call p00_f ( problem, xb, fb )

    i = 0
    write ( *, '(a)' ) ' '
    write ( *, '(i6)' ) i
    write ( *, '(a,3g16.8)' ) '  X:', xa, xb, xc
    write ( *, '(a,3g16.8)' ) '  F:', fa, fb, fc

    do i = 1, max_step

      xd = 0.5D+00 * ( xa + xb )
      call p00_f ( problem, xd, fd )

      xe = 0.5D+00 * ( xb + xc )
      call p00_f ( problem, xe, fe )

      if ( fd <= fb ) then
        xc = xb
        fc = fb
        xb = xd
        fb = fd
      else if ( fe <= fb ) then
        xa = xb
        fa = fb
        xb = xe
        fb = fe
      else
        xa = xd
        fa = fd
        xc = xe
        fc = fe
      end if

      write ( *, '(i6)' ) i
      write ( *, '(a,3g16.8)' ) '  X:', xa, xb, xc
      write ( *, '(a,3g16.8)' ) '  F:', fa, fb, fc

    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 carries out a version of Brent's derivative-free minimizer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fx
  real ( kind = 8 ) p00_fmin
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) problem
  character ( len = 50 ) title
  real ( kind = 8 ), parameter :: tol = 0.000001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For each problem, use Brent''s method.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )

    call p00_interval ( problem, xa, xb )

    call p00_f ( problem, xa, fa )
    call p00_f ( problem, xb, fb )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Initial interval [A,B]:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g16.8,14x,g16.8)' ) '   A,       B:', xa,     xb
    write ( *, '(a,g16.8,14x,g16.8)' ) '  FA,      FB:', fa,     fb

    x = p00_fmin ( xa, xb, problem, tol )

    call p00_f ( problem, xa, fa )
    call p00_f ( problem, xb, fb )
    call p00_f ( problem, x, fx )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Final interval [A,X*,B]:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g16.8,g16.8,g16.8)' ) '   A,  X*,  B:', xa, x,  xb
    write ( *, '(a,g16.8,g16.8,g16.8)' ) '  FA, FX*, FB:', fa, fx, fb

  end do

  return
end
