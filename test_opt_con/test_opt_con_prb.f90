program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_OPT_CON_PRB.
!
!  Discussion:
!
!    TEST_OPT_CON_PRB calls the TEST_OPT_CON tests.
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

  call timestamp (  )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_OPT_CON_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_OPT_CON library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_OPT_CON_PRB'
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
!    18 October 2011
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
  write ( *, '(a)' ) '  Problem    Title'
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
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: f(:)
  real ( kind = 8 ) fs
  integer ( kind = 4 ) i
  integer ( kind = 4 ) know
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: n = 100000
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  integer ( kind = 4 ) seed
  character ( len = 50 ) title
  real ( kind = 8 ), allocatable :: x(:,:)
  real ( kind = 8 ), allocatable :: xs(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For each problem, evaluate the function at many points.'
  write ( *, '(a,i8)' ) '  Number of sample points = ', n
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  allocate ( f(1:n) )

  do problem = 1, problem_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem

    call p00_title ( problem, title )

    write ( *, '(2x,a)' ) trim ( title )

    call p00_m ( problem, m )

    write ( *, '(2x,a,i8)' ) '  M =     ', m

    allocate ( a(1:m) )
    allocate ( b(1:m) )
 
    call p00_ab ( problem, m, a, b )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I      A(i)      B(i)'
    write ( *, '(a)' ) ' '

    do i = 1, m
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, a(i), b(i)
    end do

    seed = 123456789
    allocate ( x(1:m,1:n) )
    call r8col_uniform ( m, n, a, b, seed, x )
    call p00_f ( problem, m, n, x, f )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Max(F) = ', maxval ( f(1:n) )
    write ( *, '(a,g14.6)' ) '  Min(F) = ', minval ( f(1:n) )

    allocate ( xs(1:m) )
    know = 0
    call p00_sol ( problem, m, know, xs )
    if ( know /= 0 ) then
      call p00_f ( problem, m, 1, xs, fs )
      write ( *, '(a,g14.6)' ) '  F(X*)  = ', fs
    else
      write ( *, '(a)' ) '  X* is not given.'
    end if

    deallocate ( a )
    deallocate ( b )
    deallocate ( x )
    deallocate ( xs )

  end do

  deallocate ( f )

  return
end
