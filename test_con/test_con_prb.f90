program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_CON_PRB.
!
!  Discussion:
!
!    TEST_CON_PRB demonstrates the TEST_CON continuation test problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) option_num
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CON_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_CON library.'
!
!  Find out how many problems are available.
!
  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  There are ', problem_num, ' test functions.'
!
!  Print the number of options for each problem.
!
  call p00_option_num_test ( problem_num )
!
!  Print the title of each problem.
!
  call p00_title_test ( problem_num )
!
!  Print the size of each problem.
!
  call p00_nvar_test ( problem_num )
!
!  Get the starting point X0 and the norm of F(X0).
!
  call p00_start_test ( problem_num )
!
!  Check the jacobian near the starting point.
!
  call p00_jac_test ( problem_num )
!
!  Check the tangent.
!
  call p00_tan_test ( problem_num )
!
!  Apply Newton's method to a slightly perturbed starting point.
!
  call p00_newton_test ( problem_num )
!
!  Get the starting stepsize.
!
  call p00_stepsize_test ( problem_num )
!
!  Run problem 1 as a target computation, seeking X(3) = 1.0.
!  Try option = 1 and option = 4.
!
  option = 1
  call p01_target_test ( option )

  option = 4
  call p01_target_test ( option )
!
!  Run problem 1 as a limit point computation, seeking TAN(LIM) = 0.
!  Try options 2, 3, 5, 6.
!
  option = 2
  call p01_limit_test ( option )

  option = 3
  call p01_limit_test ( option )

  option = 5
  call p01_limit_test ( option )

  option = 6
  call p01_limit_test ( option )
!
!  Run problem 6 as a limit point search, using each options.
!
  call p06_option_num ( option_num)

  do option = 1, option_num
    call p06_limit_test ( option )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CON_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine p00_jac_test ( problem_num )

!*****************************************************************************80
!
!! P00_JAC_TEST compares the jacobian to a finite difference estimate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) max_adif
  integer ( kind = 4 ) max_adif_i
  integer ( kind = 4 ) max_adif_j
  real ( kind = 8 ) max_rdif
  integer ( kind = 4 ) max_rdif_i
  integer ( kind = 4 ) max_rdif_j
  integer ( kind = 4 ) nvar
  integer ( kind = 4 ) option
  integer ( kind = 4 ) option_num
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: x0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P00_JAC_TEST'
  write ( *, '(a)' ) '  Find the maximum relative difference between the'
  write ( *, '(a)' ) '  jacobian and a finite difference estimate.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem    Option      Diff               I         J'

  seed = 123456789

  do problem = 1, problem_num

    call p00_option_num ( problem, option_num )
    write ( *, '(a)' ) ' '

    do option = 1, option_num

      call p00_nvar ( problem, option, nvar )

      allocate ( x0(1:nvar) )

      call p00_start ( problem, option, nvar, x0 )

      do i = 1, nvar
        x0(i) = x0(i) + r8_uniform_01 ( seed )
      end do

      call p00_jac_check ( problem, option, nvar, x0, max_adif, max_adif_i, &
        max_adif_j, max_rdif, max_rdif_i, max_rdif_j )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,i8,2x,i8)' ) &
        problem, option, max_rdif, max_rdif_i, max_rdif_j

      deallocate ( x0 )

    end do
  end do

  return
end
subroutine p00_newton_test ( problem_num )

!*****************************************************************************80
!
!! P00_NEWTON_TEST applies Newton's method to the perturbed starting point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  real      ( kind = 8 ) dx
  real      ( kind = 8 ), allocatable, dimension ( : ) :: fx0
  real      ( kind = 8 ), allocatable, dimension ( : ) :: fx1
  real      ( kind = 8 ), allocatable, dimension ( : ) :: fx2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) nvar
  integer   ( kind = 4 ) option
  integer   ( kind = 4 ) option_num
  integer   ( kind = 4 ) par_index
  integer   ( kind = 4 ) problem
  integer   ( kind = 4 ) problem_num
  real      ( kind = 8 ) r
  real      ( kind = 8 ) r8_uniform_01
  real      ( kind = 8 ) r8vec_norm_l2
  integer   ( kind = 4 ) seed
  integer   ( kind = 4 ) status
  character ( len = 80 ) title
  real      ( kind = 8 ), allocatable, dimension ( : ) :: x0
  real      ( kind = 8 ), allocatable, dimension ( : ) :: x1
  real      ( kind = 8 ), allocatable, dimension ( : ) :: x2

  do problem = 1, problem_num

    call p00_option_num ( problem, option_num )

    write ( *, '(a)' ) ' '

    do option = 1, option_num

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P00_NEWTON_TEST'
      write ( *, '(a,i8)' ) '  Problem number = ', problem
      write ( *, '(a,i8)' ) '  Using option OPTION = ', option
!
!  Get the title.
!
      call p00_title ( problem, option, title )

      write ( *, '(2x,a)' ) trim ( title )
!
!  Get the number of variables.
!
      call p00_nvar ( problem, option, nvar )

      write ( *, '(a,i8)' ) '  Number of variables is ', nvar

      allocate ( fx0(nvar) )
      allocate ( fx1(nvar) )
      allocate ( fx2(nvar) )
      allocate ( x0(nvar) )
      allocate ( x1(nvar) )
      allocate ( x2(nvar) )
!
!  Get the starting point.
!
      call p00_start ( problem, option, nvar, x0 )
!
!  Perturb the starting point.
!
      do i = 1, nvar
        r = r8_uniform_01 ( seed )
        dx = 0.10D+00 * r * ( x0(i) + sign ( 1.0D+00, x0(i) ) )
        x1(i) = x0(i) + dx
      end do
!
!  Choose a continuation parameter index.
!
      call p00_par_index ( problem, option, nvar, x1, par_index )

      write ( *, '(a,i3,a,g14.6)' ) '  Fixing X(', par_index,  ') = ', x1(par_index)
!
!  Apply Newton's method.
!
      x2(1:nvar) = x1(1:nvar)

      call p00_newton ( problem, option, nvar, x2, par_index, status )

      if ( status == -3 ) then
        write ( *, '(a)' ) '  The convergence test was not satisfied.'
      else if ( status == -2 ) then
        write ( *, '(a)' ) '  The iteration seemed to be diverging, and was halted.'
      else if ( status == -1 ) then
        write ( *, '(a)' ) '  The jacobian was singular, and the iteration was halted.'
      else
        write ( *, '(a,i8,a)' ) '  Convergence was achieved in ', status, ' steps.'
      end if

      if ( nvar <= 10 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '        X0               X1=X0+dX       X2'
        write ( *, '(a)' ) ' '
        do i = 1, nvar
          write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x0(i), x1(i), x2(i)
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      ||X0||           ||X1=X0+dX||    ||X2||'
        write ( *, '(a)' ) ' '
        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          r8vec_norm_l2 ( nvar, x0 ), r8vec_norm_l2 ( nvar, x1 ), r8vec_norm_l2 ( nvar, x2 )

      end if
!
!  Compute the function values.
!
      call p00_fun ( problem, option, nvar, x0, fx0 )
      call p00_fun ( problem, option, nvar, x1, fx1 )
      call p00_fun ( problem, option, nvar, x2, fx2 )

      if ( nvar <= 10 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       F(X0)           F(X1=X0+dX)        F(X2)'
        write ( *, '(a)' ) ' '
        do i = 1, nvar - 1
          write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) fx0(i), fx1(i), fx2(i)
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '     ||F(X0)||       ||F(X1=X0+dX)||      ||F(X2)||'
        write ( *, '(a)' ) ' '
        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          r8vec_norm_l2 ( nvar-1, fx0 ), r8vec_norm_l2 ( nvar-1, fx1 ), &
          r8vec_norm_l2 ( nvar-1, fx2 )

      end if

      deallocate ( fx0 )
      deallocate ( fx1 )
      deallocate ( fx2 )
      deallocate ( x0 )
      deallocate ( x1 )
      deallocate ( x2 )

    end do

  end do

  return
end
subroutine p00_nvar_test ( problem_num )

!*****************************************************************************80
!
!! P00_NVAR_TEST prints the problem size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) nvar
  integer ( kind = 4 ) option
  integer ( kind = 4 ) option_num
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P00_NVAR_TEST'
  write ( *, '(a)' ) '  List the problem size.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem    Option      Size'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num
    call p00_option_num ( problem, option_num )
    write ( *, '(a)' ) ' '
    do option = 1, option_num
      call p00_nvar ( problem, option, nvar )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) problem, option, nvar
    end do
  end do

  return
end
subroutine p00_option_num_test ( problem_num )

!*****************************************************************************80
!
!! P00_OPTION_NUM_TEST lists the number of options.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) option_num
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P00_OPTION_NUM_TEST'
  write ( *, '(a)' ) '  List the number of options for each problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem   Options'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num
    call p00_option_num ( problem, option_num )
    write ( *, '(2x,i8,2x,i8)' ) problem, option_num
  end do

  return
end
subroutine p00_start_test ( problem_num )

!*****************************************************************************80
!
!! P00_START_TEST prints the norm of the starting point and its function value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( : ) :: fx0
  real ( kind = 8 ) fx0_norm
  integer ( kind = 4 ) nvar
  integer ( kind = 4 ) option
  integer ( kind = 4 ) option_num
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  real ( kind = 8 ) r8vec_norm_l2
  real ( kind = 8 ), allocatable, dimension ( : ) :: x0
  real ( kind = 8 ) x0_norm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P00_START_TEST'
  write ( *, '(a)' ) '  Get norms of starting point X0 and F(X0)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem    Option      ||X0||      || F(X0)||'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_option_num ( problem, option_num )
    write ( *, '(a)' ) ' '

    do option = 1, option_num

      call p00_nvar ( problem, option, nvar )

      allocate ( x0(1:nvar) )
      allocate ( fx0(1:nvar-1) )

      call p00_start ( problem, option, nvar, x0 )

      x0_norm = r8vec_norm_l2 ( nvar, x0 )

      call p00_fun ( problem, option, nvar, x0, fx0 )

      fx0_norm = r8vec_norm_l2 ( nvar - 1, fx0 )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) &
        problem, option, x0_norm, fx0_norm

      deallocate ( fx0 )
      deallocate ( x0 )

    end do
  end do

  return
end
subroutine p00_step ( problem, option, nvar, x, par_index, h, hmin, hmax, &
  xt, ht, status )

!*****************************************************************************80
!
!! P00_STEP takes one continuation step.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the starting point.
!
!    Input, integer ( kind = 4 ) PAR_INDEX, the continuation parameter.
!    If the program is free to choose this value, set it to 0.
!
!    Input, real ( kind = 8 ) H, HMIN, HMAX, the suggested step, and
!    the minimum and maximum stepsizes.  H may be negative.
!
!    Output, real ( kind = 8 ) XT(NVAR), the computed point.
!
!    Output, real ( kind = 8 ) HT, the  actual stepsize that was used.
!
!    Output, integer ( kind = 4 ) STATUS, the status of the calculation.
!    nonnegative, successful.
!    negative, the Newton iteration failed repeatedly even when the
!    minimum stepsize was used.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) h
  integer ( kind = 4 ) h_reduction
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  real ( kind = 8 ) ht
  integer ( kind = 4 ) option
  integer ( kind = 4 ) par_index
  integer ( kind = 4 ) problem
  real ( kind = 8 ) r8_sign
  integer ( kind = 4 ) status
  real ( kind = 8 ) tan(nvar)
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xt(nvar)
!
!  Compute the tangent.
!
  call p00_tan ( problem, option, nvar, x, tan )
!
!  Estimate the next point.
!
  h_reduction = 0

  do

    xt(1:nvar) = x(1:nvar) + h * tan(1:nvar)
!
!  Use the Newton method.
!
    call p00_newton ( problem, option, nvar, xt, par_index, status )

    if ( 0 <= status ) then

      if ( h_reduction == 0 ) then

        if ( status <= 1 ) then
          h = h * 4.0D+00
        else if ( status <= 3 ) then
          h = h * 2.0D+00
        end if

        if ( hmax < abs ( h ) ) then
          h = hmax * r8_sign ( h )
        end if

      end if

      exit

    end if

    if ( abs ( h ) <= hmin ) then
      exit
    end if

    if ( hmin < abs ( h ) ) then
      h = h / 2.0D+00
      if ( abs ( h ) < hmin ) then
        h = hmin * r8_sign ( h )
      end if
      h_reduction = h_reduction + 1
    end if

  end do

  return
end
subroutine p00_stepsize_test ( problem_num )

!*****************************************************************************80
!
!! P00_STEPSIZE_TEST prints the stepsizes for each problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option
  integer ( kind = 4 ) option_num
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P00_STEPSIZE_TEST'
  write ( *, '(a)' ) '  Print the stepsizes for each problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   Problem    Option      H               HMIN             HMAX'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_option_num ( problem, option_num )
    write ( *, '(a)' ) ' '

    do option = 1, option_num

      call p00_stepsize ( problem, option, h, hmin, hmax )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        problem, option, h, hmin, hmax

    end do
  end do

  return
end
subroutine p00_tan_test ( problem_num )

!*****************************************************************************80
!
!! P00_TAN_TEST computes and tests the tangent vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROBLEM_NUM, the number of problems.
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: jac
  real ( kind = 8 ), allocatable, dimension ( : ) :: jt
  real ( kind = 8 ) jtd
  real ( kind = 8 ) jtn
  integer ( kind = 4 ) nvar
  integer ( kind = 4 ) option
  integer ( kind = 4 ) option_num
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  real ( kind = 8 ) r8vec_norm_l2
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: tan
  real ( kind = 8 ), allocatable, dimension ( : ) :: x0

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P00_TAN_TEST'
  write ( *, '(a)' ) '  Compute the tangent vector TAN(X) at the starting point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that JAC(X) * TAN(X) = 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that det ( JAC ) > 0'
  write ( *, '(a)' ) '                  ( TAN )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem    Option    ||Jac*Tan||     det(Jac|Tan)'

  do problem = 1, problem_num

    call p00_option_num ( problem, option_num )

    write ( *, '(a)' ) ' '

    do option = 1, option_num

      call p00_nvar ( problem, option, nvar )

      allocate ( x0(1:nvar) )

      call p00_start ( problem, option, nvar, x0 )

      allocate ( jac(1:nvar,1:nvar) )

      call p00_jac ( problem, option, nvar, x0, jac )

      allocate ( tan(1:nvar) )

      call p00_tan ( problem, option, nvar, x0, tan )

      allocate ( jt(1:nvar) )

      jt = matmul ( jac, tan )

      jtn = r8vec_norm_l2 ( nvar, jt )

      jac(nvar,1:nvar) = tan(1:nvar)

      call r8mat_det ( nvar, jac, jtd )

      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) problem, option, jtn, jtd

      deallocate ( jac )
      deallocate ( jt )
      deallocate ( tan )
      deallocate ( x0 )

    end do
  end do

  return
end
subroutine p00_title_test ( problem_num )

!*****************************************************************************80
!
!! P00_TITLE_TEST prints the problem titles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  integer   ( kind = 4 ) option
  integer   ( kind = 4 ) option_num
  integer   ( kind = 4 ) problem
  integer   ( kind = 4 ) problem_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P00_TITLE_TEST'
  write ( *, '(a)' ) '  List the problem title'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem    Option  Title'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num
    call p00_option_num ( problem, option_num )
    write ( *, '(a)' ) ' '
    do option = 1, option_num
      call p00_title ( problem, option, title )
      write ( *, '(2x,i8,2x,i8,2x,a)' ) problem, option, trim ( title )
    end do
  end do

  return
end
subroutine p01_limit_test ( option )

!*****************************************************************************80
!
!! P01_LIMIT_TEST seeks limit points for problem 1.
!
!  Discussion:
!
!    We want to find points X such that TAN(LIM) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ferdinand Freudenstein, Bernhard Roth,
!    Numerical Solutions of Nonlinear Equations,
!    Journal of the Association for Computing Machinery,
!    Volume 10, 1963, Pages 550-556.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!    1, target point search, starting point is (15,-2,0);
!    2, limit points in X(1), starting point is (15,-2,0);
!    3, limit points in X(3), starting point is (15,-2,0);
!    4, target point search, starting point is (4,3,0);
!    5, limit points in X(1), starting point is (4,3,0);
!    6, limit points in X(3), starting point is (4,3,0).
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) h_reduction
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  real ( kind = 8 ) ht
  integer ( kind = 4 ) lim
  integer ( kind = 4 ) lim_num
  integer ( kind = 4 ) nvar
  integer ( kind = 4 ) option
  integer ( kind = 4 ) par_index
  integer ( kind = 4 ) problem
  real ( kind = 8 ) r8_sign
  integer ( kind = 4 ) status
  integer ( kind = 4 ) step
  integer ( kind = 4 ) step_max
  real ( kind = 8 ), allocatable, dimension ( : ) :: tan
  real ( kind = 8 ), allocatable, dimension ( : ) :: tan1
  real ( kind = 8 ), allocatable, dimension ( : ) :: tan2
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ), allocatable, dimension ( : ) :: x1
  real ( kind = 8 ), allocatable, dimension ( : ) :: x2

  problem = 1

  if ( option == 2 ) then
    lim = 1
  else if ( option == 3 ) then
    lim = 3
  else if ( option == 5 ) then
    lim = 1
  else if ( option == 6 ) then
    lim = 3
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_LIMIT_TEST'
    write ( *, '(a,i8)' ) '  Unexpected value of OPTION = ', option
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P01_LIMIT_TEST'
  write ( *, '(a)' ) '  Compute a series of solutions for problem 1.'
  write ( *, '(a)' ) '  We are trying to find limit points X such that'
  write ( *, '(a,i2,a)' ) '  TAN(', lim, ') = 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The option chosen is ', option
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         #    Tan(LIM)           X(1)           X(2)            X(3)'
  write ( *, '(a)' ) ' '
!
!  Get problem size.
!
  call p00_nvar ( problem, option, nvar )

  allocate ( tan(1:nvar) )
  allocate ( tan1(1:nvar) )
  allocate ( tan2(1:nvar) )

  allocate ( x(1:nvar) )
  allocate ( x1(1:nvar) )
  allocate ( x2(1:nvar) )

  lim_num = 0
!
!  Get starting point.
!
  call p00_start ( problem, option, nvar, x2 )
!
!  Get the tangent vector.
!
  call p00_tan ( problem, option, nvar, x2, tan2 )
!
!  Determine the appropriate continuation index.
!
  call r8vec_amax_index ( nvar, tan2, par_index )
!
!  Force F(X) = 0.
!
  step = -1
  write ( *, '(2x,i8,2x,g14.6,3(2x,g14.6))' ) step, tan2(lim), x2(1:nvar)

  call p00_newton ( problem, option, nvar, x2, par_index, status )

  if ( status < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Newton iteration failed on starting point.'
    return
  end if

  call p00_tan ( problem, option, nvar, x2, tan2 )

  step = 0
  write ( *, '(2x,i8,2x,g14.6,3(2x,g14.6))' ) step, tan2(lim), x2(1:nvar)
!
!  Get stepsize.
!
  call p00_stepsize ( problem, option, h, hmin, hmax )
!
!  LOOP:
!
  step_max = 40

  do step = 1, step_max
!
!  Save old data as X1, TAN1.
!
    x1(1:nvar) = x2(1:nvar)
    tan1(1:nvar) = tan2(1:nvar)

    h_reduction = 0
!
!  Use X1 + H * TAN1 as a starting estimate for Newton iteration.
!
    do

      if ( hmax < abs ( h ) ) then
        h = hmax * r8_sign ( h )
      end if

      if ( abs ( h ) < hmin ) then
        h = hmin * r8_sign ( h )
      end if

      x2(1:nvar) = x1(1:nvar) + h * tan1(1:nvar)

      par_index = 0
      call p00_newton ( problem, option, nvar, x2, par_index, status )
!
!  If we didn't get it, can we try again?
!
      if ( status < 0 ) then

        if ( abs ( h ) <= hmin ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Cannot decrease stepsize any more!'
          write ( *, '(a)' ) '  Search must terminate.'
          return
        else
          h = h / 4.0D+00
          h_reduction = h_reduction + 1
        end if
!
!  We computed the point.
!  Should we change the stepsize?
!
      else

        if ( h_reduction == 0 ) then

          if ( status <= 1 ) then
            h = h * 4.0D+00
          else if ( status <= 3 ) then
            h = h * 2.0D+00
          else if ( 12 <= status ) then
            h = h / 4.0D+00
          else if ( 8 <= status ) then
            h = h / 2.0D+00
          end if

        end if

        exit

      end if

    end do
!
!  Compute the tangent vector.
!
    call p00_tan ( problem, option, nvar, x2, tan2 )
!
!  Check for a limit point.
!
    if ( tan1(lim) * tan2(lim) <= 0.0D+00 ) then
      call p00_limit ( problem, option, nvar, x1, tan1, x2, tan2, lim, x, tan, status )
      write ( *, '(2x,a8,2x,g14.6,3(2x,g14.6))' ) '(limit) ',   tan(lim), x(1:nvar)
      lim_num = lim_num + 1
    end if

    write ( *, '(2x,i8,2x,g14.6,3(2x,g14.6))' ) step, tan2(lim), x2(1:nvar)

    if ( step == step_max ) then
      exit
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of limit points found was ', lim_num

  deallocate ( tan )
  deallocate ( tan1 )
  deallocate ( tan2 )
  deallocate ( x )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine p01_target_test ( option )

!*****************************************************************************80
!
!! P01_TARGET_TEST seeks a target solution for problem 1.
!
!  Discussion:
!
!    We want to compute a solution which has X(3) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ferdinand Freudenstein, Bernhard Roth,
!    Numerical Solutions of Nonlinear Equations,
!    Journal of the Association for Computing Machinery,
!    Volume 10, 1963, Pages 550-556.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!    1, target point search, starting point is (15,-2,0);
!    2, limit points in X(1), starting point is (15,-2,0);
!    3, limit points in X(3), starting point is (15,-2,0);
!    4, target point search, starting point is (4,3,0);
!    5, limit points in X(1), starting point is (4,3,0);
!    6, limit points in X(3), starting point is (4,3,0).
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) h
  integer ( kind = 4 ) h_reduction
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  real ( kind = 8 ) ht
  integer ( kind = 4 ) option
  integer ( kind = 4 ) par_index
  integer ( kind = 4 ) problem
  real ( kind = 8 ) r8_sign
  integer ( kind = 4 ) status
  integer ( kind = 4 ) step
  integer ( kind = 4 ) step_max
  real ( kind = 8 ), allocatable, dimension ( : ) :: tan
  integer ( kind = 4 ) tar_index
  real ( kind = 8 ) tar_value
  real ( kind = 8 ), allocatable, dimension ( : ) :: x1
  real ( kind = 8 ), allocatable, dimension ( : ) :: x2
  real ( kind = 8 ), allocatable, dimension ( : ) :: xt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P01_TARGET_TEST'
  write ( *, '(a)' ) '  Compute a series of solutions for problem 1.'
  write ( *, '(a)' ) '  We are trying to find a solution for which X(3) = 1.0'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The option chosen is ', option
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         #       X1              X2              X3'
  write ( *, '(a)' ) ' '

  problem = 1

  tar_index = 3
  tar_value = 1.0D+00
!
!  Get problem size.
!
  call p00_nvar ( problem, option, nvar )

  allocate ( tan(1:nvar) )
  allocate ( x1(1:nvar) )
  allocate ( x2(1:nvar) )
  allocate ( xt(1:nvar) )
!
!  Get starting point.
!
  call p00_start ( problem, option, nvar, x2 )
!
!  Get the tangent vector.
!
  call p00_tan ( problem, option, nvar, x2, tan )
!
!  Determine the appropriate continuation index.
!
  call r8vec_amax_index ( nvar, tan, par_index )
!
!  Force F(X) = 0.
!
  step = -1
  write ( *, '(2x,i8,3(2x,g14.6))' ) step, x2(1:nvar)

  call p00_newton ( problem, option, nvar, x2, par_index, status )

  if ( status < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Newton iteration failed on starting point.'
    return
  end if

  step = 0
  write ( *, '(2x,i8,3(2x,g14.6))' ) step, x2(1:nvar)
!
!  Get stepsize.
!
  call p00_stepsize ( problem, option, h, hmin, hmax )
!
!  LOOP:
!
  step_max = 40

  do step = 1, step_max

    x1(1:nvar) = x2(1:nvar)
!
!  Compute the tangent vector.
!
    call p00_tan ( problem, option, nvar, x1, tan )

    h_reduction = 0
!
!  Use X + H * TAN as a starting estimate for Newton iteration.
!
    do

      if ( hmax < abs ( h ) ) then
        h = hmax * r8_sign ( h )
      end if

      if ( abs ( h ) < hmin ) then
        h = hmin * r8_sign ( h )
      end if

      x2(1:nvar) = x1(1:nvar) + h * tan(1:nvar)

      par_index = 0
      call p00_newton ( problem, option, nvar, x2, par_index, status )
!
!  If we didn't get it, can we try again?
!
      if ( status < 0 ) then

        if ( abs ( h ) <= hmin ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Cannot decrease stepsize any more!'
          write ( *, '(a)' ) '  Did not reach target point.'
          return
        else
          h = h / 4.0D+00
          h_reduction = h_reduction + 1
        end if
!
!  We computed the point.
!  Should we change the stepsize?
!
      else

        if ( h_reduction == 0 ) then

          if ( status <= 1 ) then
            h = h * 4.0D+00
          else if ( status <= 3 ) then
            h = h * 2.0D+00
          else if ( 12 <= status ) then
            h = h / 4.0D+00
          else if ( 8 <= status ) then
            h = h / 2.0D+00
          end if

        end if

        exit

      end if

    end do
!
!  Check for the target point.
!
    if ( ( x1(tar_index) - tar_value ) * ( x2(tar_index) - tar_value ) <= 0.0D+00 ) then
      call p00_target ( problem, option, nvar, x1, x2, tar_index, tar_value, xt, status )
      write ( *, '(2x,a8,3(2x,g14.6))' ) '(target)',   xt(1:nvar)
      write ( *, '(2x,i8,3(2x,g14.6))' ) step, x2(1:nvar)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Reached target point.'
      exit
    else
      write ( *, '(2x,i8,3(2x,g14.6))' ) step, x2(1:nvar)
    end if

    if ( step == step_max ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Did not reach target point.'
      exit
    end if

  end do

  deallocate ( tan )
  deallocate ( x1 )
  deallocate ( x2 )
  deallocate ( xt )

  return
end
subroutine p06_limit_test ( option )

!*****************************************************************************80
!
!! P06_LIMIT_TEST seeks limit points for problem 6.
!
!  Discussion:
!
!    We want to find points X such that TAN(7) = 0.
!
!    The number of limit points that may be expected depends on the option:
!
!    There are five options, which vary in the value they fix the
!    elevator value in function 6:
!
!      Option   Elevator Value    Limit Points
!
!       1        -0.050              1
!       2        -0.008              3
!       3         0.0                2
!       4         0.05               1
!       5         0.1                1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Raman Mehra, William Kessel, James Carroll,
!    Global stability and contral analysis of aircraft at high angles of attack,
!    Technical Report CR-215-248-1, -2, -3,
!    Office of Naval Research, June 1977.
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    Albert Schy, Margery Hannah,
!    Prediction of Jump Phenomena in Roll-coupled Maneuvers of Airplanes,
!    Journal of Aircraft,
!    Volume 14, Number 4, 1977,  pages 375-382.
!
!    John Young, Albert Schy, Katherine Johnson,,
!    Prediction of Jump Phenomena in Aircraft Maneuvers, Including
!    Nonlinear Aerodynamic Effects,
!    Journal of Guidance and Control,
!    Volume 1, Number 1, 1978, pages 26-31.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!    Input, integer ( kind = 4 ) OPTION, the option index.
!    1, X(6) fixed at -0.050, X(7) free, X(8) fixed at 0.0;
!    2, X(6) fixed at -0.008, X(7) free, X(8) fixed at 0.0;
!    3, X(6) fixed at  0.000, X(7) free, X(8) fixed at 0.0;
!    4, X(6) fixed at  0.050, X(7) free, X(8) fixed at 0.0;
!    5, X(6) fixed at  0.100, X(7) free, X(8) fixed at 0.0.
!
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) h_reduction
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  real ( kind = 8 ) ht
  integer ( kind = 4 ) lim
  integer ( kind = 4 ) lim_num
  integer ( kind = 4 ) nvar
  integer ( kind = 4 ) option
  integer ( kind = 4 ) par_index
  integer ( kind = 4 ) problem
  real ( kind = 8 ) r8_sign
  integer ( kind = 4 ) status
  integer ( kind = 4 ) step
  integer ( kind = 4 ) step_max
  real ( kind = 8 ), allocatable, dimension ( : ) :: tan
  real ( kind = 8 ), allocatable, dimension ( : ) :: tan1
  real ( kind = 8 ), allocatable, dimension ( : ) :: tan2
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ), allocatable, dimension ( : ) :: x1
  real ( kind = 8 ), allocatable, dimension ( : ) :: x2

  problem = 6
  lim = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P06_LIMIT_TEST'
  write ( *, '(a)' ) '  Compute a series of solutions for problem 6.'
  write ( *, '(a)' ) '  We are trying to find limit points X such that'
  write ( *, '(a,i2,a)' ) '  TAN(', lim, ') = 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The option chosen is ', option
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   #   Tan(7)     X1       X2       X3       X4       X5       X6       X7       X8'
  write ( *, '(a)' ) '                Roll     Pitch    Yaw      Attack   Sideslip Elevator Aileron  Rudder'
  write ( *, '(a)' ) ' '
!
!  Get problem size.
!
  call p00_nvar ( problem, option, nvar )

  allocate ( tan(1:nvar) )
  allocate ( tan1(1:nvar) )
  allocate ( tan2(1:nvar) )

  allocate ( x(1:nvar) )
  allocate ( x1(1:nvar) )
  allocate ( x2(1:nvar) )

  lim_num = 0
!
!  Get starting point.
!
  call p00_start ( problem, option, nvar, x2 )
!
!  Get the tangent vector.
!
  call p00_tan ( problem, option, nvar, x2, tan2 )
!
!  For correction of initial point, use variable index 7.
!
  par_index = 7
!
!  Force F(X) = 0.
!
  step = -1
  write ( *, '(2x,i2,1x,f8.5,8(1x,f8.5))' ) step, tan2(lim), x2(1:nvar)

  call p00_newton ( problem, option, nvar, x2, par_index, status )

  if ( status < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Newton iteration failed on starting point.'
    return
  end if

  call p00_tan ( problem, option, nvar, x2, tan2 )

  step = 0
  write ( *, '(2x,i2,1x,f8.5,8(1x,f8.5))' ) step, tan2(lim), x2(1:nvar)
!
!  Get stepsize.
!
  call p00_stepsize ( problem, option, h, hmin, hmax )
!
!  LOOP:
!
  step_max = 30

  do step = 1, step_max
!
!  Save old data as X1, TAN1.
!
    x1(1:nvar) = x2(1:nvar)
    tan1(1:nvar) = tan2(1:nvar)

    h_reduction = 0
!
!  Use X1 + H * TAN1 as a starting estimate for Newton iteration.
!
    do

      if ( hmax < abs ( h ) ) then
        h = hmax * r8_sign ( h )
      end if

      if ( abs ( h ) < hmin ) then
        h = hmin * r8_sign ( h )
      end if

      x2(1:nvar) = x1(1:nvar) + h * tan1(1:nvar)

      par_index = 0
      call p00_newton ( problem, option, nvar, x2, par_index, status )
!
!  If we didn't get it, can we try again?
!
      if ( status < 0 ) then

        if ( abs ( h ) <= hmin ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Cannot decrease stepsize any more.'
          write ( *, '(a)' ) '  Cannot complete the computation.'
          return
        else
          h = h / 4.0D+00
          h_reduction = h_reduction + 1
        end if
!
!  We computed the point.
!  Should we change the stepsize?
!
      else

        if ( h_reduction == 0 ) then

          if ( status <= 1 ) then
            h = h * 4.0D+00
          else if ( status <= 3 ) then
            h = h * 2.0D+00
          else if ( 12 <= status ) then
            h = h / 4.0D+00
          else if ( 8 <= status ) then
            h = h / 2.0D+00
          end if

        end if

        exit

      end if

    end do
!
!  Compute the tangent vector.
!
    call p00_tan ( problem, option, nvar, x2, tan2 )
!
!  Check for a limit point.
!
    if ( tan1(lim) * tan2(lim) <= 0.0D+00 ) then
      call p00_limit ( problem, option, nvar, x1, tan1, x2, tan2, lim, x, tan, status )
      write ( *, '(2x,a2,1x,f8.5,8(1x,f8.5))' ) ' L', tan(lim), x(1:nvar)
      lim_num = lim_num + 1
    end if

    write ( *, '(2x,i2,1x,f8.5,8(1x,f8.5))' ) step, tan2(lim), x2(1:nvar)

    if ( step == step_max ) then
      exit
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of limit points found was ', lim_num

  deallocate ( tan )
  deallocate ( tan1 )
  deallocate ( tan2 )
  deallocate ( x )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine p06_test ( option )

!*****************************************************************************80
!
!! P06_TEST seeks solutions for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Raman Mehra, William Kessel, James Carroll,
!    Global stability and contral analysis of aircraft at high angles of attack,
!    Technical Report CR-215-248-1, -2, -3,
!    Office of Naval Research, June 1977.
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    Albert Schy, Margery Hannah,
!    Prediction of Jump Phenomena in Roll-coupled Maneuvers of Airplanes,
!    Journal of Aircraft,
!    Volume 14, Number 4, 1977,  pages 375-382.
!
!    John Young, Albert Schy, Katherine Johnson,,
!    Prediction of Jump Phenomena in Aircraft Maneuvers, Including
!    Nonlinear Aerodynamic Effects,
!    Journal of Guidance and Control,
!    Volume 1, Number 1, 1978, pages 26-31.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!    1, X(6) fixed at -0.050, X(7) free, X(8) fixed at 0.0;
!    2, X(6) fixed at -0.008, X(7) free, X(8) fixed at 0.0;
!    3, X(6) fixed at  0.000, X(7) free, X(8) fixed at 0.0;
!    4, X(6) fixed at  0.050, X(7) free, X(8) fixed at 0.0;
!    5, X(6) fixed at  0.100, X(7) free, X(8) fixed at 0.0.
!
  implicit none

  integer ( kind = 4 ) nvar

  character ( len = 80 ) file_t_name
  character ( len = 80 ) file_x_name
  integer   ( kind = 4 ) file_t_unit
  integer   ( kind = 4 ) file_x_unit
  real      ( kind = 8 ) h
  integer   ( kind = 4 ) h_reduction
  real      ( kind = 8 ) hmax
  real      ( kind = 8 ) hmin
  real      ( kind = 8 ) ht
  integer   ( kind = 8 ) i
  logical, parameter :: make_files = .false.
  integer   ( kind = 4 ) option
  integer   ( kind = 4 ) par_index
  integer   ( kind = 4 ) problem
  real      ( kind = 8 ) r8_sign
  integer   ( kind = 4 ) status
  integer   ( kind = 4 ) step
  integer   ( kind = 4 ) step_max
  real      ( kind = 8 ), allocatable, dimension ( : ) :: tan
  real      ( kind = 8 ), allocatable, dimension ( : ) :: x1
  real      ( kind = 8 ), allocatable, dimension ( : ) :: x2

  if ( make_files ) then
    write ( file_x_name, '(a,i1,a)' ) 'p06_opt0', option, '_x000.txt'
    write ( file_t_name, '(a,i1,a)' ) 'p06_opt0', option, '_t000.txt'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P06_TEST'
  write ( *, '(a)' ) '  Compute a series of solutions for problem 6.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The option chosen is ', option
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   #    X1       X2       X3       X4       X5       X6       X7       X8'
  write ( *, '(a)' ) '      Roll     Pitch    Yaw      Attack   Sideslip Elevator Aileron  Rudder'
  write ( *, '(a)' ) ' '

  problem = 6
!
!  Get problem size.
!
  call p00_nvar ( problem, option, nvar )

  allocate ( tan(1:nvar) )
  allocate ( x1(1:nvar) )
  allocate ( x2(1:nvar) )
!
!  Get starting point.
!
  call p00_start ( problem, option, nvar, x2 )
!
!  Get the tangent vector.
!
  call p00_tan ( problem, option, nvar, x2, tan )
!
!  For correction of initial point, use variable index 7.
!
  par_index = 7
!
!  Force F(X) = 0.
!
  step = -1
  write ( *, '(2x,i2,8(1x,f8.5))' ) step, x2(1:nvar)

  call p00_newton ( problem, option, nvar, x2, par_index, status )

  if ( status < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Newton iteration failed on starting point.'
    return
  end if

  call p00_tan ( problem, option, nvar, x2, tan )

  step = 0
  write ( *, '(2x,i2,8(1x,f8.5))' ) step, x2(1:nvar)

  if ( make_files ) then

    call get_unit ( file_x_unit )
    open ( unit = file_x_unit, file = file_x_name, status = 'replace' )
    do i = 1, nvar
      write ( file_x_unit, '(g16.8)' ) x2(i)
    end do
    close ( unit = file_x_unit )
    call file_name_inc ( file_x_name )

    call get_unit ( file_t_unit )
    open ( unit = file_t_unit, file = file_t_name, status = 'replace' )
    do i = 1, nvar
      write ( file_t_unit, '(g16.8)' ) tan(i)
    end do
    close ( unit = file_t_unit )
    call file_name_inc ( file_t_name )

  end if
!
!  Get stepsize.
!
  call p00_stepsize ( problem, option, h, hmin, hmax )
!
!  LOOP:
!
  step_max = 30

  do step = 1, step_max

    x1(1:nvar) = x2(1:nvar)
!
!  Compute the tangent vector.
!
    call p00_tan ( problem, option, nvar, x1, tan )

    h_reduction = 0
!
!  Use X + H * TAN as a starting estimate for Newton iteration.
!
    do

      if ( hmax < abs ( h ) ) then
        h = hmax * r8_sign ( h )
      end if

      if ( abs ( h ) < hmin ) then
        h = hmin * r8_sign ( h )
      end if

      x2(1:nvar) = x1(1:nvar) + h * tan(1:nvar)

      par_index = 0
      call p00_newton ( problem, option, nvar, x2, par_index, status )
!
!  If we didn't get it, can we try again?
!
      if ( status < 0 ) then

        if ( abs ( h ) <= hmin ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Cannot decrease stepsize any more!'
          write ( *, '(a)' ) '  Did not reach target point.'
          return
        else
          h = h / 4.0D+00
          h_reduction = h_reduction + 1
        end if
!
!  We computed the point.
!  Should we change the stepsize?
!
      else

        if ( h_reduction == 0 ) then

          if ( status <= 1 ) then
            h = h * 4.0D+00
          else if ( status <= 3 ) then
            h = h * 2.0D+00
          else if ( 12 <= status ) then
            h = h / 4.0D+00
          else if ( 8 <= status ) then
            h = h / 2.0D+00
          end if

        end if

        exit

      end if

    end do

    call p00_tan ( problem, option, nvar, x2, tan )
    write ( *, '(2x,i2,8(1x,f8.5))' ) step, x2(1:nvar)

    if ( make_files ) then

      call get_unit ( file_x_unit )
      open ( unit = file_x_unit, file = file_x_name, status = 'replace' )
      do i = 1, nvar
        write ( file_x_unit, '(g16.8)' ) x2(i)
      end do
      close ( unit = file_x_unit )
      call file_name_inc ( file_x_name )

      call get_unit ( file_t_unit )
      open ( unit = file_t_unit, file = file_t_name, status = 'replace' )
      do i = 1, nvar
        write ( file_t_unit, '(g16.8)' ) tan(i)
      end do
      close ( unit = file_t_unit )
      call file_name_inc ( file_t_name )

    end if

    if ( step == step_max ) then
      exit
    end if

  end do

  deallocate ( tan )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
