program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_OPT_PRB.
!
!  Discussion:
!
!    TEST_OPT_PRB calls the TEST_OPT tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp (  )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_OPT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_OPT library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_OPT_PRB'
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
!    08 January 2008
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
  write ( *, '(a)' ) 'Problem	 Title'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(i6,2x,a)' ) problem, title

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
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  real ( kind = 8 ) f_sol
  real ( kind = 8 ) f_start
  integer ( kind = 4 ) know
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) problem
  real ( kind = 8 ) p
  integer ( kind = 4 ) seed
  character ( len = 50 ) title
  real ( kind = 8 ) x(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For each problem, evaluate the function'
  write ( *, '(a)' ) '  at the starting point.'
!
!  Initialize the P25 parameter.
!
  p = 0.0D+00
  call p25_p_set ( p )
!
!  Initialize the P36 problems.
!
  seed = 123456789
  call p36_p_init ( seed )
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem

    call p00_title ( problem, title )

    write ( *, '(2x,a)' ) trim ( title )

    call p00_n ( problem, n )

    n_min = n

    if ( n < 0 ) then
      n_min = abs ( n_min )
      n = max ( n_min, 4 )
    end if

    write ( *, '(2x,a,i8)' ) '  N_MIN = ', n_min
    write ( *, '(2x,a,i8)' ) '  N =     ', n
    write ( *, '(a)' ) ' '

    if ( n_max < n ) then
      write ( *, '(a)' ) '  Problem too large.'
      cycle
    end if
 
    call p00_start ( problem, n, x )

    call p00_f ( problem, n, x, f_start )

    write ( *, '(4x,a,g14.6)' ) '  F(X_START)=', f_start

    know = 0
    do
      call p00_sol ( problem, n, know, x )
      if ( know == 0 ) then
        exit
      end if
      call p00_f ( problem, n, x, f_sol )
      write ( *, '(4x,a,g14.6)' ) '  F(X_SOL)=  ', f_sol
    end do

    if ( problem == 25 ) then

      p = 1.0D+00
      call p25_p_set ( p )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Repeat problem with P = ', p

      call p00_start ( problem, n, x )

      call p00_f ( problem, n, x, f_start )

      write ( *, '(4x,a,g14.6)' ) '  F(X_START)=', f_start

      know = 0
      do
        call p00_sol ( problem, n, know, x )
        if ( know == 0 ) then
          exit
        end if
        call p00_f ( problem, n, x, f_sol )
        write ( *, '(4x,a,g14.6)' ) '  F(X_SOL)=  ', f_sol
      end do
        
    end if

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 compares the exact and approximate gradient vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  real ( kind = 8 ), allocatable, dimension ( : ) :: g
  real ( kind = 8 ), allocatable, dimension ( : ) :: g_dif
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) problem
  real ( kind = 8 ) p
  integer ( kind = 4 ) seed
  character ( len = 50 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
!
!  Initialize the P25 parameter.
!
  p = 0.0D+00
  call p25_p_set ( p )
!
!  Initialize the P36 problems.
!
  seed = 123456789
  call p36_p_init ( seed )

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

    call p00_n ( problem, n )

    n_min = n

    if ( n < 0 ) then
      n_min = abs ( n_min )
      n = max ( n_min, 4 )
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(2x,a,i8)' ) 'N =     ', n

    if ( n_max < n ) then
      write ( *, '(a)' ) '  Problem too big.'
      cycle
    end if

    allocate ( g(1:n) )
    allocate ( g_dif(1:n) )
    allocate ( x(1:n) )

    if ( .true. ) then

      call random_number ( harvest = x(1:n) )

    else

     call p00_start ( problem, n, x )

    end if

    call p00_g ( problem, n, x, g )

    call p00_gdif ( problem, n, x, g_dif )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  X'
    write ( *, '(4x,5g14.6)' ) x(1:n)
    write ( *, '(a)' ) '  G'
    write ( *, '(4x,5g14.6)' ) g(1:n)
    write ( *, '(a)' ) '  G_DIF'
    write ( *, '(4x,5g14.6)' ) g_dif(1:n)

    if ( problem == 25 ) then

      p = 1.0D+00
      call p25_p_set ( p )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Repeat problem with P = ', p

      call p00_start ( problem, n, x )

      call p00_g ( problem, n, x, g )

      call p00_gdif ( problem, n, x, g_dif )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X'
      write ( *, '(4x,5g14.6)' ) x(1:n)
      write ( *, '(a)' ) '  G'
      write ( *, '(4x,5g14.6)' ) g(1:n)
      write ( *, '(a)' ) '  G_DIF'
      write ( *, '(4x,5g14.6)' ) g_dif(1:n)

    end if

    deallocate ( g )
    deallocate ( g_dif )
    deallocate ( x )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 compares the exact and approximate Hessian matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: h
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: h_dif
  integer ( kind = 4 ) i
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) problem
  real ( kind = 8 ) p
  integer ( kind = 4 ) seed
  character ( len = 50 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
!
!  Initialize the P25 parameter.
!
  p = 0.0D+00
  call p25_p_set ( p )
!
!  Initialize the P36 problems.
!
  seed = 123456789
  call p36_p_init ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For each problem, compare the exact and'
  write ( *, '(a)' ) '  approximate Hessians at the starting point.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    call p00_n ( problem, n )

    n_min = n

    if ( n < 0 ) then
      n_min = abs ( n_min )
      n = max ( n_min, 4 )
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(2x,a,i8)' ) '  N =     ', n

    if ( n_max < n ) then
      write ( *, '(a)' ) '  This problem is too large.'
      cycle
    end if

    allocate ( h(1:n,1:n) )
    allocate ( h_dif(1:n,1:n) )
    allocate ( x(1:n) )

    if ( .true. ) then

      call random_number ( harvest = x(1:n) )

    else

      call p00_start ( problem, n, x )

    end if

    call p00_h ( problem, n, x, h )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  X:'
    write ( *, '(4x,5g14.6)' ) x(1:n)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  H:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(4x,6g13.5)' ) h(i,1:n)
    end do

    call p00_hdif ( problem, n, x, h_dif )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  H_DIF:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(4x,6g13.5)' ) h_dif(i,1:n)
    end do

    if ( problem == 25 ) then

      p = 1.0D+00
      call p25_p_set ( p )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Repeat problem with P = ', p

      call p00_start ( problem, n, x )

      call p00_h ( problem, n, x, h )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X:'
      write ( *, '(4x,5g14.6)' ) x(1:n)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  H:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(4x,6g13.5)' ) h(i,1:n)
      end do

      call p00_hdif ( problem, n, x, h_dif )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  HDIF:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(4x,6g13.5)' ) h_dif(i,1:n)
      end do

    end if

    deallocate ( h )
    deallocate ( h_dif )
    deallocate ( x )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 carries out a simple gradient method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  real ( kind = 8 ) g(n_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: max_step = 6
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) problem
  real ( kind = 8 ) p
  logical quit
  integer ( kind = 4 ) reduce
  integer ( kind = 4 ), parameter :: reduce_max = 10
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  character ( len = 50 ) title
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) x2(n_max)
!
!  Initialize the P25 parameter.
!
  p = 0.0D+00
  call p25_p_set ( p )
!
!  Initialize the P36 problems.
!
  seed = 123456789
  call p36_p_init ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For each problem, take a few steps of '
  write ( *, '(a)' ) '  the gradient method.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    call p00_n ( problem, n )

    n_min = n

    if ( n < 0 ) then
      n_min = abs ( n_min )
      n = max ( n_min, 4 )
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(2x,a,i8)' ) 'N =     ', n

    if ( n_max < n ) then
      write ( *, '(a)' )  '  Problem too large.'
      cycle
    end if

    call p00_start ( problem, n, x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a)' ) '  Starting X:'
    write ( *, '(4x,5g14.6)' ) x(1:n)
    call p00_f ( problem, n, x, f )
    write ( *, '(2x,a)' ) '  Starting F:'
    write ( *, '(4x,g14.6)' ) f

    quit = .false.

    do i = 1, max_step

      s = 1.0D+00

      reduce = 0

      call p00_g ( problem, n, x, g )
      write ( *, '(2x,a)' ) '  Gradient:'
      write ( *, '(4x,5g14.6)' ) g(1:n)

      if ( all ( g(1:n) == 0.0D+00 ) ) then
        write ( *, '(2x,a)' ) '  Terminate because of zero gradient.'
        exit
      end if

      do

        x2(1:n) = x(1:n) - s * g(1:n)
 
        call p00_f ( problem, n, x2, f2 )

        if ( f2 < f ) then
!       if ( f2 <= 0.95D+00 * f ) then
          exit
        end if

        reduce = reduce + 1
        write ( *, '(2x,a,g14.6)' ) '  Reject step, F = ', f2

        if ( reduce_max < reduce ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a)' ) '  Repeated step reductions do not help.'
          write ( *, '(2x,a)' ) '  Problem abandoned.'
          quit = .true.
          exit
        end if

        s = s / 4.0D+00

      end do

      if ( quit ) then
        exit
      end if

      x(1:n) = x2(1:n)
      f = f2

      write ( *, '(a)' ) ' '
      write ( *, '(2x,a)' ) '  New X:'
      write ( *, '(4x,5g14.6)' ) x(1:n)
      call p00_f ( problem, n, x, f )
      write ( *, '(2x,a)' ) '  New F:'
      write ( *, '(4x,g14.6)' ) f

    end do

  end do

  return
end
