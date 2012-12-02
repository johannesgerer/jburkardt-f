program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_INTERP_ND_PRB.
!
!  Discussion:
!
!    TEST_INTERP_ND_PRB calls the TEST_INTERP_ND tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call timestamp (  )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_ND_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_INTERP_ND library.'
  write ( *, '(a)' ) '  The R8LIB library is also needed.'

  call test01 ( )

  n = 10
  do m = 2, 4
    call test02 ( m, n )
  end do

  n = 2
  do m = 2, 4
    call test03 ( m, n )
  end do

  n = 10000
  m = 4
  call test04 ( m, n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_ND_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 prints the title of each test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  P00_TITLE returns the problem title.'

  call p00_prob_num ( prob_num )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i2,a)' ) '  There are a total of ', prob_num, ' problems.'
  write ( *, '(a)' ) ' '

  do prob = 1, prob_num
    call p00_title ( prob, title )
    write ( *, '(2x,i2,2x,a)' ) prob, '"' // trim ( title ) // '"'
  end do

  return
end
subroutine test02 ( m, n )

!*****************************************************************************80
!
!! TEST02 samples each function in M dimensions, at N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  P00_F evaluates the function.'
  write ( *, '(a,i4)' ) '  Here, we use spatial dimension M = ', m
  write ( *, '(a,i6)' ) '  The number of points is N = ', n

  seed = 123456789
  call r8mat_uniform_01 ( m, n, seed, x )

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem  ', prob

    call p00_c ( prob, m, seed, c )
    call r8vec_print ( m, c, '  C parameters:' )

    call p00_w ( prob, m, seed, w )
    call r8vec_print ( m, w, '  W parameters:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      F(X)              X(1)      X(2) ...'
    write ( *, '(a)' ) ' '

    call p00_f ( prob, m, c, w, n, x, f )

    do j = 1, n
      write ( *, '(2x,g14.6,2x,10f10.4)' ) f(j), x(1:m,j)
    end do

  end do

  return
end
subroutine test03 ( m, n )

!*****************************************************************************80
!
!! TEST03 samples each derivative component in M dimensions, at N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) d(m,n)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) id
  integer ( kind = 4 ) j
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  P00_D evaluates derivative components.'
  write ( *, '(a,i4)' ) '  Here, we use spatial dimension M = ', m
  write ( *, '(a,i6)' ) '  The number of points is N = ', n

  seed = 123456789
  call r8mat_uniform_01 ( m, n, seed, x )

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem  ', prob

    call p00_c ( prob, m, seed, c )
    call r8vec_print ( m, c, '  C parameters:' )

    call p00_w ( prob, m, seed, w )
    call r8vec_print ( m, w, '  W parameters:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '                         X(1)      X(2) ...'
    write ( *, '(a)' ) '      F(X)            dFdX(1)   dFdX(2) ...'

    call p00_f ( prob, m, c, w, n, x, f )

    do id = 1, m
      call p00_d ( prob, m, id, c, w, n, x, d(id,1:n) )
    end do

    do j = 1, n
      write ( *, '(a)' ) ' '
      write ( *, '(2x,14x,2x,10f10.4)' )       x(1:m,j)
      write ( *, '(2x,g14.6,2x,10f10.4)' ) f(j), d(1:m,j)
    end do

  end do

  return
end
subroutine test04 ( m, n )

!*****************************************************************************80
!
!! TEST04 estimates integrals in M dimensions, using N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) q1
  real ( kind = 8 ) q2
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  P00_Q returns the integral of F over [0,1]^m.'
  write ( *, '(a,i4)' ) '  Here, we use spatial dimension M = ', m
  write ( *, '(a,i6)' ) '  The number of sample points is N = ', n

  seed = 123456789
  call r8mat_uniform_01 ( m, n, seed, x )

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem  ', prob

    call p00_c ( prob, m, seed, c )
    call r8vec_print ( m, c, '  C parameters:' )

    call p00_w ( prob, m, seed, w )
    call r8vec_print ( m, w, '  W parameters:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      Exact Integral     Q'
    write ( *, '(a)' ) ' '

    call p00_q ( prob, m, c, w, q1 )

    call p00_f ( prob, m, c, w, n, x, f )
    q2 = sum ( f(1:n) ) / real ( n, kind = 8 )

    write ( *, '(2x,g14.6,2x,g14.6)' ) q1, q2

  end do

  return
end
